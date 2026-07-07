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
#include "render/led.h"   // PIN constants, NoColorCorrection, NoTempCorrection, USE_DMA_LEDS
#include "pov_single_map.h"  // pure strip index math (host-tested)

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

/**
 * @brief Manages the display loop for a single-Teensy POV rig.
 * @tparam S The number of physical segments/pixels.
 * @tparam RPM The rotations per minute of the device.
 *
 * Used by Holosphere (96x20 / 288x144 with one Teensy owning the full strip).
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
    randomSeed(1337); // FastLED LCG only; modern effects use hs::random() (platform.h)
  #ifdef USE_DMA_LEDS
    ledController_.begin();
    ledController_.setCorrection(255, 176, 240);  // TypicalLEDStrip
    ledController_.setTemperature(255, 147, 41);   // Candle
    ledController_.setBrightness(255);
  #else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(leds_,
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
   * first frame so the ISR never observes a half-filled table, then constructs,
   * configures, runs, and deletes the effect.
   */
  template <typename E> void show(unsigned long duration) {
    GeometryResolution<E>::init();
    E *e = new E();
    configure_arenas_default(); // Reset before init so effects can override
    e->init();
    run(e, duration);
    delete e;
  }

private:
  /**
   * @brief Non-template core of show(): drives the column ISR for the effect's
   * lifetime.
   * @param e Effect to run; borrowed, not owned. The caller (show) retains
   * ownership and deletes it.
   * @param duration The time in seconds to run the effect.
   * @details Publishes e to the ISR-visible effect_ only while the timer ISR is
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
    HS_CHECK(duration <= ~0UL / 1000UL,
             "show duration too long (duration*1000 ms overflows unsigned long)");
    const unsigned long duration_ms = duration * 1000;
    effect_ = e;
    // show_col() indexes buf[y * width + x] for y in [0, S/2), in-bounds only
    // when the effect's canvas height equals the strip's half-height.
    HS_CHECK(effect_->height() == S / 2,
             "POVDisplay: effect canvas height must equal S/2");
    x_ = 0;
    IntervalTimer timer;
    // One column sweep period (µs), rounded to nearest; a pathological RPM/width
    // could round it to 0 µs, an undefined IntervalTimer period.
    static_assert(RPM > 0, "POVDisplay: RPM must be positive");
    const unsigned long cols_per_min =
        static_cast<unsigned long>(RPM) * effect_->width();
    HS_CHECK(cols_per_min > 0, "column sweep rate is zero (width is 0)");
    const unsigned long interval_us =
        (60000000UL + cols_per_min / 2) / cols_per_min;
    HS_CHECK(interval_us >= 1,
        "column interval rounded to 0 µs (RPM/width too high)");
    HS_CHECK(timer.begin(show_col, interval_us),
        "column IntervalTimer failed to start (no PIT channel)");
    while (millis() - start < duration_ms) {
      unsigned long t0 = micros();
      effect_->draw_frame();
      unsigned long dt = micros() - t0;
      if (hs::debug) {
        Serial.print("ft ");
        Serial.println(dt);
      }
    }
    timer.end();
    effect_ = nullptr;
  }

private:
  /**
   * @brief Static function called by the IntervalTimer to display one column of
   * the frame.
   */
  static FASTRUN void show_col() {
    // Top half reads column x_; bottom half reads the opposite column (x_+W/2).
    const int w = effect_->width();
    const int x_top = x_;
    const int x_bot = pov::strip_opposite_col(x_, w);

    // ISR fast path: index the display buffer directly, dropping the per-column
    // virtual get_pixel() dispatches. The one effect that overrides get_pixel
    // (RingTwist) routes to the slow path via this probe.
    const bool slow = effect_->overrides_get_pixel();
    const Pixel *buf = slow ? nullptr : effect_->display_buffer();

#if defined(USE_DMA_LEDS)
    auto& frame = ledController_.backFrame();
    for (int y = 0; y < S / 2; ++y) {
      // Top half is wired reversed, bottom half straight.
      frame.packPixel(pov::strip_top_led(y, S),
          slow ? effect_->get_pixel(x_top, y) : buf[y * w + x_top]);
      frame.packPixel(pov::strip_bottom_led(y, S),
          slow ? effect_->get_pixel(x_bot, y) : buf[y * w + x_bot]);
    }
    // Overrun result discarded: at the ~1.3 ms single-board column period a DMA
    // overrun cannot occur, so no overrun watchdog here.
    (void)ledController_.submitFrame(effect_->strobe_columns());
#else
    for (int y = 0; y < S / 2; ++y) {
      leds_[pov::strip_top_led(y, S)] = static_cast<CRGB>(
          slow ? effect_->get_pixel(x_top, y) : buf[y * w + x_top]);
      leds_[pov::strip_bottom_led(y, S)] = static_cast<CRGB>(
          slow ? effect_->get_pixel(x_bot, y) : buf[y * w + x_bot]);
    }
    FastLED.show();
    if (effect_->strobe_columns()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }
#endif

    x_ = (x_ + 1) % w;
    // Advance the display buffer at each half-revolution boundary.
    if (x_ == 0 || x_ == w / 2) {
      effect_->advance_display();
    }
  }

#ifndef USE_DMA_LEDS
  static CRGB
      leds_[S]; /**< Array holding the CRGB data for the physical LED strip. */
#endif
  static Effect
      *effect_; /**< Pointer to the currently running effect instance. */
  static int
      x_; /**< Current column index being displayed (virtual position). */
#if defined(USE_DMA_LEDS)
  static DMALEDController<S>
      ledController_; /**< HD107S DMA controller driving the physical strip. */
#endif
};

template <int S, int RPM> int POVDisplay<S, RPM>::x_ = 0;

template <int S, int RPM> Effect *POVDisplay<S, RPM>::effect_ = nullptr;

#ifndef USE_DMA_LEDS
template <int S, int RPM> CRGB POVDisplay<S, RPM>::leds_[S];
#endif

#if defined(USE_DMA_LEDS)
// ledController_ is intentionally NOT defined out-of-line here. Its HD107SFrame
// buffers are the eDMA TX source and belong in cached, DMA-reachable OCRAM (where
// HD107SFrame's arm_dcache_flush() write-back is meaningful — in DTCM it is a
// dead no-op). DMAMEM (a section attribute) is silently dropped by GCC on a
// vague-linkage template static member, so a generic definition would land in
// DTCM regardless. Each instantiating single-board target must instead invoke
// HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM) once at file scope; it emits the
// required explicit specialization, whose ordinary strong linkage keeps the
// DMAMEM section attribute.
// POVSegmented carries the same contract (see targets/Phantasm/Phantasm.ino).
#define HS_DEFINE_POV_SINGLE_LED_CONTROLLER(S, RPM)                            \
  template <>                                                                  \
  DMAMEM DMALEDController<S> POVDisplay<S, RPM>::ledController_ {}
#endif

#endif // ARDUINO
