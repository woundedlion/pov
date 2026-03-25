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

#ifdef ARDUINO
#include <Arduino.h>
  #ifdef USE_DMA_LEDS
    #include "dma_led.h"
  #else
    #include <FastLED.h>
  #endif
#endif
#include "canvas.h"
#include "memory.h"

/**
 * @brief Manages the display loop for a single-Teensy POV rig.
 * @tparam S The number of physical segments/pixels.
 * @tparam RPM The rotations per minute of the device.
 *
 * Used by Holosphere (96x20 / 288x144 with one Teensy owning the full strip).
 */
template <int S, int RPM> class POVDisplay {
public:
  /**
   * @brief Initializes the LED strip and sets up hardware specific
   * optimizations.
   */
  POVDisplay() {
#ifdef ARDUINO
    randomSeed(1337);
  #ifdef USE_DMA_LEDS
    ledController_.begin();
    ledController_.setCorrection(255, 176, 240);  // TypicalLEDStrip
    ledController_.setTemperature(255, 147, 41);   // Candle
    ledController_.setBrightness(255);
    // Enable slew rate limiting for signal integrity
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 & ~IOMUXC_PAD_SRE; // Pin 11
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 & ~IOMUXC_PAD_SRE; // Pin 13
  #else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(leds_,
                                                                        S);
    // enable slew rate limiting
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 & ~IOMUXC_PAD_SRE; // Pin 11
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 & ~IOMUXC_PAD_SRE; // Pin 13
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  #endif
#else
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
#endif
  }

  /**
   * @brief Runs a specific Effect for a given duration.
   * @tparam E The Effect class to run.
   * @param duration The time in seconds to run the effect.
   */
  template <typename E> void show(unsigned long duration) {
    effect_ = new E();
    effect_->init();  // Static dispatch — FLASHMEM honored via explicit instantiation
    run(effect_, duration);
    delete effect_;
    effect_ = nullptr;
  }

private:
  /**
   * @brief Non-template core of show(). Takes ownership of e.
   */
  void run(Effect *e, unsigned long duration) {
    configure_arenas_default();

    long start = millis();
    effect_ = e;
    x_ = 0;
#ifdef ARDUINO
    IntervalTimer timer;
    // sweep the width once per rotation
    timer.begin(show_col,
        static_cast<unsigned long>(1000000.0f / (RPM / 60.0f) / effect_->width() + 0.5f));
    while (millis() - start < duration * 1000) {
      unsigned long t0 = micros();
      effect_->draw_frame();
      effect_->advance_display();
//      delay(125);
      unsigned long dt = micros() - t0;
      Serial.print("frame ms: ");
      Serial.println(dt);
    }
    timer.end();
#else
    while (millis() - start < duration * 1000) {
      unsigned long t0 = micros();
      effect_->draw_frame();
      unsigned long dt = micros() - t0;
      hs::log("ft %lu", dt);
    }
#endif
  }

private:
  /**
   * @brief Static function called by the IntervalTimer to display one column of
   * the frame.
   */
  static FASTRUN void show_col() {
    return;

#if defined(ARDUINO) && defined(USE_DMA_LEDS)
    // Direct Pixel16 → HD107S wire packing (no intermediate CRGB array).
    auto& frame = ledController_.backFrame();
    for (int y = 0; y < S / 2; ++y) {
      frame.packPixel(S / 2 - y - 1, effect_->get_pixel(x_, y));
      frame.packPixel(S / 2 + y, effect_->get_pixel(
          (x_ + (effect_->width() / 2)) % effect_->width(), y));
    }
    ledController_.submitFrame(effect_->show_bg());
#else
    for (int y = 0; y < S / 2; ++y) {
      // Map to physical strip: top half is inverted, bottom half is straight.
      leds_[S / 2 - y - 1] = effect_->get_pixel(x_, y);
      leds_[S / 2 + y] = effect_->get_pixel(
          (x_ + (effect_->width() / 2)) % effect_->width(), y);
    }
    FastLED.show();
    if (effect_->show_bg()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }
#endif

    x_ = (x_ + 1) % effect_->width();
    // When the POV sweep completes one full virtual revolution (x_ = 0 or x_ =
    // width/2), advance the display buffer.
    if (x_ == 0 || x_ == effect_->width() / 2) {
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
#if defined(ARDUINO) && defined(USE_DMA_LEDS)
  static DMALEDController<S> ledController_;
#endif
};

template <int S, int RPM> int POVDisplay<S, RPM>::x_ = 0;

template <int S, int RPM> Effect *POVDisplay<S, RPM>::effect_ = nullptr;

#ifndef USE_DMA_LEDS
template <int S, int RPM> CRGB POVDisplay<S, RPM>::leds_[S];
#endif

#if defined(ARDUINO) && defined(USE_DMA_LEDS)
template <int S, int RPM>
DMALEDController<S> POVDisplay<S, RPM>::ledController_;
#endif
