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
#include "led.h"   // PIN constants, NoColorCorrection, NoTempCorrection, USE_DMA_LEDS
#include "pov_single_map.h"  // pure strip index math (host-tested)

// Like POVSegmented, this driver is Arduino-only: it depends on IntervalTimer,
// FastLED/DMA, and the Teensy runtime, and is instantiated solely by the
// Holosphere .ino target. The pure strip index math it relies on lives in
// pov_single_map.h (above, host-tested); the class itself is compiled out off
// device so there is no half-built non-Arduino path to misuse.
#ifdef ARDUINO
#include <Arduino.h>
  #ifdef USE_DMA_LEDS
    #include "dma_led.h"
  #else
    #include <FastLED.h>
  #endif
#include "canvas.h"
#include "geometry.h"
#include "memory.h"

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

  /**
   * @brief HD107S SPI clock for the single-board Holosphere DMA path, in Hz.
   *
   * The Holosphere runs one Teensy over the full strip at a wide column period
   * (~1.3 ms for 96 columns at 480 RPM), so the TeensySPIDMA 12 MHz default
   * clears the per-column transfer with ample headroom; the 24 MHz clock used
   * by the segmented driver is not needed here. Stated explicitly rather than
   * left to the constructor default so the per-board clock choice is on the
   * record at the call site.
   */
  static constexpr uint32_t SPI_CLOCK_HZ = 12000000;

public:
  /**
   * @brief Constructs the driver, initializing the LED strip and hardware-
   * specific optimizations (correction, temperature, brightness).
   * @details Seeds FastLED's legacy LCG to 1337; modern effects draw from the
   * separate hs::random() mt19937(1337) reproduced by the simulator.
   */
  POVDisplay() {
    // Seeds FastLED's LCG only (legacy effects); modern effects draw from
    // hs::random() mt19937(1337). See the determinism contract in platform.h.
    randomSeed(1337);
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
    // One column sweep period (µs), rounded; a pathological RPM/width could
    // floor it to 0 µs, an undefined IntervalTimer period.
    const unsigned long interval_us = static_cast<unsigned long>(
        1000000.0f / (RPM / 60.0f) / effect_->width() + 0.5f);
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
    // submitFrame()'s overrun result is intentionally discarded here: at the
    // single-board column period (~1.3 ms) a DMA overrun is realistically
    // impossible, so this path runs without the overrun watchdog the segmented
    // driver keeps (it tracks getOverrunCount()).
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
// DTCM regardless. Each instantiating target therefore defines it as an explicit
// specialization (ordinary strong linkage, so DMAMEM sticks) — see Phantasm.ino.
#endif

#endif // ARDUINO
