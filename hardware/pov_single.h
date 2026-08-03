/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_single.h
 * @brief Single-Teensy POV display driver for the Holosphere.
 *
 * One Teensy controls one LED strip spanning both sides of the ring.
 * The IntervalTimer ISR sweeps columns at a rate derived from RPM and
 * the virtual canvas width.
 *
 * Include directly from target .ino files:
 *   #include "../../hardware/pov_single.h"
 */
#pragma once
#include "render/led.h" // PIN constants, NoColorCorrection, NoTempCorrection, USE_DMA_LEDS
#include "pov_single_map.h" // pure strip index math (host-tested)

// Arduino-only: depends on IntervalTimer, FastLED/DMA, and the Teensy runtime.
// Pure strip index math is in pov_single_map.h (host-tested).
#ifdef ARDUINO
#include <Arduino.h>
#ifdef USE_DMA_LEDS
#include "dma_led.h"
#else
#include <FastLED.h>
#endif
#include "render/canvas.h"
#include "math/geometry.h"
#include "engine/memory.h"
#include <new> // std::nothrow — fail-fast OOM check on the effect allocation

/**
 * @brief Manages the display loop for a single-Teensy POV rig.
 * @tparam S The number of physical segments/pixels.
 * @tparam RPM The rotations per minute of the device.
 *
 * Used by Holosphere (96x20 with one Teensy owning the full strip).
 */
template <int S, int RPM> class POVDisplay {
  static_assert(S > 0 && S % 2 == 0,
                "POVDisplay requires an even, positive segment count S");

public:
  /**
   * @brief Constructs the driver, initializing the LED strip and hardware-
   * specific optimizations (correction, temperature, brightness).
   * @details Seeds FastLED's legacy LCG to 1337; modern effects draw from the
   * separate hs::random() mt19937(1337) reproduced by the simulator.
   */
  POVDisplay() {
    randomSeed(
        1337); // FastLED LCG only; modern effects use hs::random() (platform.h)
#ifdef USE_DMA_LEDS
    ledController.begin();
    ledController.setCorrection(255, 176, 240); // TypicalLEDStrip
    ledController.setTemperature(255, 147, 41); // Candle
    ledController.setBrightness(255);
#else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(leds,
                                                                        S);
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
    FastLED.setBrightness(255);
#endif
  }

  /**
   * @brief Runs a specific Effect for a given duration.
   * @tparam E The Effect class to run.
   * @param duration The time in seconds to run the effect.
   * @details Eagerly fills the scanline LUTs for E's resolution before the
   * first frame so the ISR never observes a half-filled table, then configures
   * the arenas, constructs, runs, and deletes the effect.
   */
  template <typename E> void show(unsigned long duration) {
    GeometryResolution<E>::init();
    configure_arenas_default(); // Reset before init so effects can override
    E *e = new (std::nothrow) E();
    HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
    e->init();
    run(e, duration);
    delete e;
  }

private:
#if defined(USE_DMA_LEDS)
  /**
   * @brief Worst-case duration of one column's LED transfer, in µs.
   * @details Image frame plus the trailing black frame strobe_columns() appends,
   * eight clocks per byte at the clock the default-constructed ledController
   * runs at. Rounded up, so run()'s column-period check never under-counts.
   */
  static constexpr unsigned long COLUMN_TRANSFER_US =
      static_cast<unsigned long>(
          (static_cast<uint64_t>(HD107SFrame<S>::COMPOSITE_SIZE) * 8u *
               1000000u +
           DMALEDController<S>::DEFAULT_CLOCK_HZ - 1) /
          DMALEDController<S>::DEFAULT_CLOCK_HZ);
#endif

  /**
   * @brief Non-template core of show(): drives the column ISR for the effect's
   * lifetime.
   * @param e Effect to run; borrowed, not owned. The caller (show) retains
   * ownership and deletes it.
   * @param duration The time in seconds to run the effect.
   * @details Publishes e to the ISR-visible effect only while the timer ISR is
   * attached, then unpublishes it; does not delete e. Traps a driver/effect
   * resolution mismatch and a failed timer start before either becomes a dark
   * strip or an unguarded OOB read in the column ISR.
   */
  void run(Effect *e, unsigned long duration) {
    // Unsigned start + (millis() - start) stays correct across the millis()
    // wraparound; a signed start would mis-compare on overflow.
    const unsigned long start = millis();
    // duration * 1000 must fit in unsigned long; trap the overflow for symmetry
    // with the segmented driver's invariant checks (unreachable in practice —
    // past ~49.7 days on a 32-bit millis() clock).
    HS_CHECK(
        duration <= ~0UL / 1000UL,
        "show duration too long (duration*1000 ms overflows unsigned long)");
    const unsigned long duration_ms = duration * 1000;
    effect = e;
    // show_col() indexes buf[y * width + x] for y in [0, S/2), in-bounds only
    // when the effect's canvas height equals the strip's half-height.
    HS_CHECK(effect->height() == S / 2,
             "POVDisplay: effect canvas height must equal S/2");
    // strip_opposite_col's (x + w/2) % w antipode and the x == w/2 swap cadence
    // both truncate w/2 on odd width, misregistering the bottom hemisphere.
    HS_CHECK(effect->width() % 2 == 0,
             "POVDisplay: effect canvas width must be even");
    x = 0;
    IntervalTimer timer;
    // One column sweep period (µs), rounded to nearest; a pathological RPM/width
    // could round it to 0 µs, an undefined IntervalTimer period.
    static_assert(RPM > 0, "POVDisplay: RPM must be positive");
    const unsigned long cols_per_min =
        static_cast<unsigned long>(RPM) * effect->width();
    HS_CHECK(cols_per_min > 0, "column sweep rate is zero (width is 0)");
    const unsigned long interval_us =
        (60000000UL + cols_per_min / 2) / cols_per_min;
    HS_CHECK(interval_us >= 1,
             "column interval rounded to 0 µs (RPM/width too high)");
#if defined(USE_DMA_LEDS)
    // show_col() discards submitFrame()'s overrun return, and unlike the
    // segmented driver it has no retry latch and no dark fallback — a dropped
    // column would freeze the strip on the last accepted frame. Sound only
    // while one composite transfer fits inside a column period.
    HS_CHECK(interval_us > COLUMN_TRANSFER_US,
             "LED transfer outlasts the column period (S, RPM and canvas width "
             "would overrun the DMA every column)");
#endif
    HS_CHECK(timer.begin(show_col, interval_us),
             "column IntervalTimer failed to start (no PIT channel)");
#if defined(USE_DMA_LEDS)
    uint32_t last_overrun = ledController.getOverrunCount();
#endif
    while (millis() - start < duration_ms) {
      unsigned long t0 = micros();
      effect->draw_frame();
      unsigned long dt = micros() - t0;
      if (hs::debug) {
        Serial.print("ft ");
        Serial.println(dt);
#if defined(USE_DMA_LEDS)
        const uint32_t overruns = ledController.getOverrunCount();
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

private:
  /**
   * @brief Static function called by the IntervalTimer to display one column of
   * the frame.
   */
  static FASTRUN void show_col() {
    // Top half reads column x; bottom half reads the opposite column (x+W/2).
    const int w = effect->width();
    const int x_top = x;
    const int x_bot = pov::strip_opposite_col(x, w);

    // ISR fast path: index the display buffer directly, dropping the per-column
    // virtual get_pixel() dispatches. The one effect that overrides get_pixel
    // (RingTwist) routes to the slow path via this probe.
    const bool slow = effect->overrides_get_pixel();
    const Pixel *buf = slow ? nullptr : effect->display_buffer();

#if defined(USE_DMA_LEDS)
    auto &frame = ledController.backFrame();
    for (int y = 0; y < S / 2; ++y) {
      // Top half is wired reversed, bottom half straight.
      frame.packPixel(pov::strip_top_led(y, S),
                      slow ? effect->get_pixel(x_top, y) : buf[y * w + x_top]);
      frame.packPixel(pov::strip_bottom_led(y, S),
                      slow ? effect->get_pixel(x_bot, y) : buf[y * w + x_bot]);
    }
    // Overrun result discarded: run() rejects any configuration whose column
    // period does not clear COLUMN_TRANSFER_US, so no overrun watchdog here.
    (void)ledController.submitFrame(effect->strobe_columns());
#else
    for (int y = 0; y < S / 2; ++y) {
      leds[pov::strip_top_led(y, S)] = static_cast<CRGB>(
          slow ? effect->get_pixel(x_top, y) : buf[y * w + x_top]);
      leds[pov::strip_bottom_led(y, S)] = static_cast<CRGB>(
          slow ? effect->get_pixel(x_bot, y) : buf[y * w + x_bot]);
    }
    FastLED.show();
    if (effect->strobe_columns()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }
#endif

    x = (x + 1) % w;
    // Advance the display buffer at each half-revolution boundary.
    if (x == 0 || x == w / 2) {
      effect->advance_display();
    }
  }

#ifndef USE_DMA_LEDS
  static CRGB
      leds[S]; /**< Array holding the CRGB data for the physical LED strip. */
#endif
  static Effect *effect; /**< Currently running effect; the ISR only reads it.
                              Written by run() solely while the column timer is
                              detached (before begin(), after end()), never
                              mid-show. */
  static int x; /**< Current column index being displayed (virtual position). */
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
// ledController has no out-of-line definition: DMAMEM survives only on an
// explicit specialization (see DMALEDController in dma_led_controller.h). Each
// instantiating single-board target invokes
// HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM) once at file scope.
#define HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM)                            \
  template <> DMAMEM DMALEDController<S> POVDisplay<S, RPM>::ledController {}
#endif

#endif // ARDUINO
