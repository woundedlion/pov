/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_single.h
 * @brief Single-Teensy POV display driver for the Holosphere.
 *
 * One Teensy controls one LED strip spanning both sides of the ring.
 * The IntervalTimer ISR sweeps columns at a rate derived from RPM and
 * the virtual canvas width.
 */
#pragma once
#include "core/platform/led.h" // PIN_DATA, PIN_CLOCK (FastLED path)
#include "pov_single_map.h"

#ifdef ARDUINO
#include <Arduino.h>
#ifdef USE_DMA_LEDS
#include "dma_led.h"
#else
#include <FastLED.h>
#endif
#include "core/render/canvas.h"
#include "core/math/geometry.h"
#include "core/memory.h"
#include <utility>
#include <new>

/**
 * @brief Manages the display loop for a single-Teensy POV rig.
 * @tparam S Total LED count on the strip (both sides of the ring).
 * @tparam RPM The rotations per minute of the device.
 */
template <int S, int RPM> class POVDisplay {
  static_assert(S > 0 && S % 2 == 0,
                "POVDisplay requires an even, positive LED count S");

public:
  POVDisplay() = delete;
  /** Initializes the transport after Arduino core startup; call once. */
  static void begin() {
#ifdef USE_DMA_LEDS
    ledController.begin();
    ledController.set_correction(hd107s::LINEAR_STRIP_GAIN.r,
                                 hd107s::LINEAR_STRIP_GAIN.g,
                                 hd107s::LINEAR_STRIP_GAIN.b);
    ledController.set_temperature(hd107s::LINEAR_WARM_GAIN.r,
                                  hd107s::LINEAR_WARM_GAIN.g,
                                  hd107s::LINEAR_WARM_GAIN.b);
    ledController.set_brightness(255);
#else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB,
                    DATA_RATE_MHZ(pov::FASTLED_CLOCK_MHZ)>(leds, S);
    restore_correction_baseline();
    FastLED.setBrightness(255);
#endif
  }

  /**
   * @brief Runs a specific Effect for a given duration.
   * @tparam E The Effect class to run.
   * @tparam Args Effect constructor argument types.
   * @param duration The time in seconds to run the effect.
   * @param args Arguments forwarded to the effect constructor.
   * @details Fills the scanline LUTs for E's resolution before the first
   * frame, so the ISR never observes a half-filled table.
   */
  template <typename E, typename... Args>
  static void show(unsigned long duration, Args &&...args) {
    math::GeometryResolution<E>::init();
    configure_arenas_default(); // Reset before init so effects can override
    E *e = new (std::nothrow) E(std::forward<Args>(args)...);
    HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
    e->init();
    run(e, duration);
    delete e;
  }

private:
#if defined(USE_DMA_LEDS)
  /**
   * @brief HD107S SPI clock for the single-board DMA path, in Hz.
   */
  static constexpr uint32_t SPI_CLOCK_HZ =
      DMALEDController<S>::DEFAULT_CLOCK_HZ;

  /**
   * @brief Worst-case duration of one column's LED transfer, in µs.
   * @details Image frame plus the trailing strobe black frame, at
   * SPI_CLOCK_HZ.
   */
  static constexpr unsigned long COLUMN_TRANSFER_US =
      dma::transfer_us(HD107SFrame<S>::COMPOSITE_SIZE, SPI_CLOCK_HZ);
#else
  // FastLED's WS2801Controller waits up to 1000 us before each transmission.
  static constexpr unsigned long FASTLED_SHOW_US =
      pov::fastled_show_us(S, pov::FASTLED_CLOCK_MHZ);
#endif

  /**
   * @brief Non-template core of show(): drives the column ISR for the effect's
   * lifetime.
   * @param e Effect to run; borrowed, not owned.
   * @param duration The time in seconds to run the effect.
   * @details Publishes e to the ISR only while the timer is attached.
   */
  static void run(Effect *e, unsigned long duration) {
    // Unsigned (millis() - start) stays correct across the millis() wrap.
    const unsigned long start = millis();
    HS_CHECK(
        duration <= ~0UL / 1000UL,
        "show duration too long (duration*1000 ms overflows unsigned long)");
    const unsigned long duration_ms = duration * 1000;
    HS_CHECK(effect == nullptr,
             "POVDisplay::run() re-entered while an effect is live");
    effect = e;
    HS_CHECK(effect->height() == S / 2,
             "POVDisplay: effect canvas height must equal S/2");
    // Odd width truncates w/2, misregistering the bottom hemisphere.
    HS_CHECK(effect->width() % 2 == 0,
             "POVDisplay: effect canvas width must be even");
    x = 0;
    IntervalTimer timer;
    static_assert(RPM > 0, "POVDisplay: RPM must be positive");
    const unsigned long cols_per_min =
        static_cast<unsigned long>(RPM) * effect->width();
    HS_CHECK(cols_per_min > 0, "column sweep rate is zero (width is 0)");
    const float interval_us = pov::column_interval_us(cols_per_min);
    HS_CHECK(interval_us >= 1,
             "column interval below 1 µs (RPM/width too high)");
#if defined(USE_DMA_LEDS)
    HS_CHECK(interval_us > COLUMN_TRANSFER_US,
             "LED transfer outlasts the column period (S, RPM and canvas width "
             "would overrun the DMA every column)");
#else
    HS_CHECK(interval_us >
                 FASTLED_SHOW_US * (effect->strobe_columns() ? 2UL : 1UL),
             "FastLED transfer outlasts the column period");
#endif
    HS_CHECK(timer.begin(show_col, interval_us),
             "column IntervalTimer failed to start (no PIT channel)");
#if defined(USE_DMA_LEDS)
    uint32_t last_overrun = ledController.get_overrun_count();
#endif
    while (millis() - start < duration_ms) {
      unsigned long t0 = micros();
      effect->draw_frame();
      unsigned long dt = micros() - t0;
      if (hs::debug) {
        Serial.print("ft ");
        Serial.println(dt);
#if defined(USE_DMA_LEDS)
        const uint32_t overruns = ledController.get_overrun_count();
        if (overruns != last_overrun) {
          Serial.print("overrun ");
          Serial.println(overruns);
          last_overrun = overruns;
        }
#endif
      }
    }
    timer.end();
    effect = nullptr;
  }

  /**
   * @brief Static function called by the IntervalTimer to display one column of
   * the frame.
   */
  static void show_col() {
    const int w = effect->width();
    const bool slow = effect->overrides_get_pixel();
    const Pixel *buf = slow ? nullptr : effect->display_buffer();
#if defined(USE_DMA_LEDS)
    auto &frame = ledController.back_frame();
#endif
    x = pov::run_single_column<S>(
        x, w,
        [&](int column, int row) {
          return slow ? effect->get_pixel(column, row) : buf[row * w + column];
        },
        [&](int led, const Pixel &pixel) {
#if defined(USE_DMA_LEDS)
          frame.pack_pixel(led, pixel);
#else
          leds[led] = static_cast<CRGB>(pixel);
#endif
        },
        [&] {
#if defined(USE_DMA_LEDS)
          (void)ledController.submit_frame(effect->strobe_columns());
#else
          FastLED.show();
          if (effect->strobe_columns())
            FastLED.showColor(CRGB(0, 0, 0));
#endif
        },
        [&] { effect->advance_display(); });
  }

#ifndef USE_DMA_LEDS
  static CRGB
      leds[S]; /**< Array holding the CRGB data for the physical LED strip. */
#endif
  static Effect *effect; /**< Currently running effect; the ISR only reads it.
                              Written only while the column timer is detached. */
  static int x; /**< Current column index being displayed (virtual position).
                     Seeded only while the column timer is detached. */
#if defined(USE_DMA_LEDS)
  static DMALEDController<S>
      ledController; /**< HD107S DMA controller driving the physical strip. */
#endif
};

template <int S, int RPM> int POVDisplay<S, RPM>::x = 0;

template <int S, int RPM> Effect *POVDisplay<S, RPM>::effect = nullptr;

#ifndef USE_DMA_LEDS
template <int S, int RPM> CRGB POVDisplay<S, RPM>::leds[S];
#endif

#if defined(USE_DMA_LEDS)
// DMAMEM survives only on an explicit specialization, so each instantiating
// single-board target invokes HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM) once
// at file scope.
#define HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM)                            \
  template <> DMAMEM DMALEDController<S> POVDisplay<S, RPM>::ledController {   \
    POVDisplay<S, RPM>::SPI_CLOCK_HZ                                           \
  }
#endif

#endif // ARDUINO
