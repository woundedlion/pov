/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wvolatile"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include <Arduino.h>
#include <FastLED.h>
#include "constants.h"
#include "canvas.h"
#pragma GCC diagnostic pop

/**
 * @brief Analog pin used for seeding the random number generator.
 */
static constexpr int PIN_RANDOM = 15;
/**
 * @brief Data pin for the LED strip
 */
static constexpr int PIN_DATA = 11;
/**
 * @brief Clock pin for the LED strip.
 */
static constexpr int PIN_CLOCK = 13;


/**
 * @brief RAII guard to temporarily disable FastLED's color correction.
 */
struct NoColorCorrection {
  /**
   * @brief Disables standard correction in the constructor.
   */
  NoColorCorrection()
  {
    FastLED.setCorrection(UncorrectedColor);
    FastLED.setTemperature(UncorrectedTemperature);
  }

  /**
   * @brief Restores standard correction in the destructor.
   */
  ~NoColorCorrection() {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  }
};

/**
 * @brief RAII guard to temporarily disable FastLED's temperature correction.
 */
struct NoTempCorrection {
  /**
   * @brief Disables temperature correction in the constructor.
   */
  NoTempCorrection()
  {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(UncorrectedTemperature);
  }

  /**
   * @brief Restores temperature correction in the destructor.
   */
  ~NoTempCorrection() {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  }
};

/**
 * @brief Manages the display loop for the POV hardware.
 * @tparam S The number of physical segments/pixels.
 * @tparam RPM The rotations per minute of the device.
 */
template <int S, int RPM>
class POVDisplay
{
public:

  /**
   * @brief Initializes the LED strip and sets up hardware specific optimizations.
   */
  POVDisplay()
  {
    randomSeed(analogRead(PIN_RANDOM));
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(leds_, S);

    // enable slew rate limiting
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 =
      IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 & ~IOMUXC_PAD_SRE; // Pin 11
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 =
      IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 & ~IOMUXC_PAD_SRE; // Pin 13

    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  }

  /**
   * @brief Runs a specific Effect for a given duration.
   * @tparam E The Effect class to run.
   * @param duration The time in seconds to run the effect.
   */
  template <typename E>
  void show(unsigned long duration)
  {
    long start = millis();
    effect_ = new E();
    x_ = 0;
    IntervalTimer timer;
    // The timer interval is calculated to sweep the width exactly once per rotation.
    timer.begin(show_col, 1000000 / (RPM / 60) / effect_->width());
    while (millis() - start < duration * 1000) {
      effect_->draw_frame();
    }
    timer.end();
    delete effect_;
    effect_ = nullptr;
  }

private:

  /**
   * @brief Static function called by the IntervalTimer to display one column of the frame.
   */
  static inline void show_col() {

    for (int y = 0; y < S / 2; ++y) {
      // Map to physical strip: top half is inverted, bottom half is straight.
      leds_[S / 2 - y - 1] = effect_->get_pixel(x_, y);
      leds_[S / 2 + y] =
        effect_->get_pixel((x_ + (effect_->width() / 2)) % effect_->width(), y);
    }

    FastLED.show();
    if (effect_->show_bg()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }

    x_ = (x_ + 1) % effect_->width();
    // When the POV sweep completes one full virtual revolution (x_ = 0 or x_ = width/2), 
    // advance the display buffer.
    if (x_ == 0 || x_ == effect_->width() / 2) {
      effect_->advance_display();
    }
  }

  static CRGB leds_[S]; /**< Array holding the CRGB data for the physical LED strip. */
  static Effect* effect_; /**< Pointer to the currently running effect instance. */
  static int x_; /**< Current column index being displayed (virtual position). */
};

template<int S, int RPM>
int POVDisplay<S, RPM>::x_ = 0;

template<int S, int RPM>
Effect* POVDisplay<S, RPM>::effect_ = nullptr;

template<int S, int RPM>
CRGB POVDisplay<S, RPM>::leds_[S];
