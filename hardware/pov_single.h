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
  // The whole driver tiles the strip as two equal arms of S/2 rows: rows = S/2
  // truncates, show_col() loops y in [0, S/2), and the resolution HS_CHECK only
  // validates effect height == S/2. An odd S would silently drop the middle LED
  // (S=41 -> writes 40) and mis-tile the strip with no runtime trap. Guard the
  // precondition at compile time, in keeping with the fail-fast doctrine.
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
    // Seeds FastLED's LCG only (bare random8/16/random(), used by legacy
    // effects); modern effects draw from hs::random() mt19937(1337). See the
    // determinism contract in platform.h.
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
    FastLED.setBrightness(255);  // explicit (FastLED defaults to 255), so both
                                 // branches set correction+temperature+brightness
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
    // Unsigned start + unsigned elapsed (millis() - start) is the overflow-safe
    // timing idiom: the modular subtraction stays correct across the ~49.7-day
    // millis() wraparound. A signed start would mis-compare on overflow.
    const unsigned long start = millis();
    const unsigned long duration_ms = duration * 1000;
    effect_ = e;
    // ISR seam (project doctrine: trap at the cold bind site, not the hot ISR):
    // show_col() walks y in [0, S/2) and indexes the display buffer as
    // buf[y * width + x], in-bounds only when the effect's canvas height equals
    // the strip's half-height. Trap a driver/effect resolution mismatch here,
    // before it becomes an unguarded OOB read inside the column ISR.
    HS_CHECK(effect_->height() == S / 2,
             "POVDisplay: effect canvas height must equal S/2");
    x_ = 0;
    IntervalTimer timer;
    // One column sweep period (µs), rounded. A pathological RPM/width could floor
    // this to 0 µs, which gives IntervalTimer::begin an undefined period — trap
    // the degenerate timing here, beside the no-PIT-channel guard.
    const unsigned long interval_us = static_cast<unsigned long>(
        1000000.0f / (RPM / 60.0f) / effect_->width() + 0.5f);
    HS_CHECK(interval_us >= 1,
        "column interval rounded to 0 µs (RPM/width too high)");
    // sweep the width once per rotation; begin() returns false if no PIT
    // channel is free — an unstarted timer means the column ISR never fires
    // and the strip stays dark, so trap rather than render to a dead timer.
    HS_CHECK(timer.begin(show_col, interval_us),
        "column IntervalTimer failed to start (no PIT channel)");
    while (millis() - start < duration_ms) {
      unsigned long t0 = micros();
      effect_->draw_frame();
      unsigned long dt = micros() - t0;
      if (hs::debug) {
        // dt is micros() elapsed; the unit-neutral "ft " label matches the
        // segmented driver.
        Serial.print("ft ");
        Serial.println(dt);
      }
    }
    timer.end();
    effect_ = nullptr; // ISR detached above — unpublish; caller deletes e
  }

private:
  /**
   * @brief Static function called by the IntervalTimer to display one column of
   * the frame.
   */
  static FASTRUN void show_col() {
    // Top half reads column x_; bottom half reads the opposite column (x_+W/2).
    // The physical-LED ↔ canvas mapping (pov::strip_*) is host-unit-tested in
    // tests/test_pov_single.h.
    const int w = effect_->width();
    const int x_top = x_;
    const int x_bot = pov::strip_opposite_col(x_, w);

    // ISR fast path: fetch the display buffer base once and index it directly,
    // dropping the S per-column virtual get_pixel() dispatches. Sound because
    // (1) prev_ is stable for this whole column — advance_display() runs below
    // after the loop — and (2) buf[y * w + x] == get_pixel(x, y) for any effect
    // that does not override get_pixel. Only the legacy scroller (RingTwist)
    // does; one virtual probe per column routes it to the correct slow path.
    const bool slow = effect_->overrides_get_pixel();
    const Pixel *buf = slow ? nullptr : effect_->display_buffer();

#if defined(USE_DMA_LEDS)
    // Direct Pixel16 → HD107S wire packing (no intermediate CRGB array).
    auto& frame = ledController_.backFrame();
    for (int y = 0; y < S / 2; ++y) {
      // Map to physical strip: top half is reversed, bottom half is straight.
      frame.packPixel(pov::strip_top_led(y, S),
          slow ? effect_->get_pixel(x_top, y) : buf[y * w + x_top]);
      frame.packPixel(pov::strip_bottom_led(y, S),
          slow ? effect_->get_pixel(x_bot, y) : buf[y * w + x_bot]);
    }
    // Steady-state column path intentionally drops the accept/overrun result:
    // a dropped image column self-heals next tick (see submitFrame's doc).
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
    // When the POV sweep completes one full virtual revolution (x_ = 0 or x_ =
    // width/2), advance the display buffer.
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
// DMAMEM (OCRAM): the controller's HD107SFrame buffers are the actual eDMA TX
// source, so they must live in DMA-reachable, cached OCRAM — which is exactly
// what HD107SFrame's arm_dcache_flush() assumes. Default placement is DTCM,
// where that flush is a dead no-op; here it does the required write-back.
template <int S, int RPM>
DMAMEM DMALEDController<S>
    POVDisplay<S, RPM>::ledController_{POVDisplay<S, RPM>::SPI_CLOCK_HZ};
#endif

#endif // ARDUINO
