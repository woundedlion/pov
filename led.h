/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <Arduino.h>
#include <FastLED.h>
#pragma GCC diagnostic pop

/**
 * @brief Rotations Per Minute of the POV display.
 */
static constexpr unsigned int RPM = 480;
/**
 * @brief The number of physical LED pixels on the strip.
 */
static constexpr int NUM_PIXELS = 40;
/**
 * @brief Half the number of pixels (physical height).
 */
static constexpr int H = NUM_PIXELS / 2;
/**
 * @brief Virtual height used for mapping (H + buffer).
 */
static constexpr int H_VIRT = H + 3;
/**
 * @brief Maximum horizontal resolution (width) for effects.
 */
static constexpr int MAX_W = 96;
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
 * @brief Macro to calculate the 1D array index from 2D coordinates (x, y).
 * @param x The horizontal coordinate (column).
 * @param y The vertical coordinate (row).
 * @return The 1D index.
 */
inline constexpr int XY(int x, int y) { return x * H + y; }
class Canvas;

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for drawing a frame.
 */
class Effect {
  friend class Canvas;

public:
  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect.
   */
  Effect(int W) :
    persist_pixels(true),
    width_(W)
  {
    bufs_[0] = buffer_a;
    memset(bufs_[0], 0, sizeof(CRGB) * W * H);
    bufs_[1] = buffer_b;
    memset(bufs_[1], 0, sizeof(CRGB) * W * H);
  }

  virtual ~Effect() {
  };

  /**
   * @brief Abstract method to be implemented by derived classes to generate a frame of graphics data.
   */
  virtual void draw_frame() = 0;
  /**
   * @brief Indicates whether the display should show black between POV sweeps.
   * @return True if the background should be shown (often black).
   */
  virtual bool show_bg() const = 0;

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The CRGB color reference.
   */
  virtual const CRGB& get_pixel(int x, int y) const {
    return bufs_[prev_][XY(x, y)];
  }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  inline int width() const { return width_; }
  /**
   * @brief Checks if the display buffer and drawing buffer are synced.
   * @return True if the current frame is ready to be shown.
   */
  inline bool buffer_free() const { return prev_ == next_; }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   */
  inline void advance_display() { prev_ = next_; }
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content to the new buffer.
   */
  inline void advance_buffer() {
    noInterrupts();
    cur_ = cur_ ? 0 : 1;
    interrupts();
    if (persist_pixels) {
      memcpy(bufs_[cur_], bufs_[prev_], sizeof(CRGB) * width_ * H);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   */
  inline void queue_frame() {
    noInterrupts();
    next_ = cur_;
    interrupts();
  }

protected:
  /**
   * @brief Flag indicating if the previous frame's pixels should be copied to the new buffer (for trails/decay).
   */
  bool persist_pixels;

private:
  volatile int prev_ = 0, cur_ = 0, next_ = 0; /**< Pointers/indices for double buffering logic. */
  int width_; /**< The width of the effect. */
  inline static CRGB buffer_a[MAX_W * H]; /**< Static storage for buffer A. */
  inline static CRGB buffer_b[MAX_W * H]; /**< Static storage for buffer B. */
  CRGB* bufs_[2]; /**< Pointers to the two buffer storage locations. */
};

/**
 * @brief Context class providing a safe, scoped interface to the current drawing buffer.
 * @details Ensures the buffer advances correctly when the object is constructed and queued when destroyed.
 */
class Canvas {
public:
  /**
   * @brief Constructs the Canvas, advancing the effect buffer and optionally clearing it.
   * @param effect The effect instance owning the buffer.
   */
  Canvas(Effect& effect) : effect_(effect) {
    while (!effect_.buffer_free()) {}
    //    start_time = millis();
    effect_.advance_buffer();
    if (!effect_.persist_pixels) {
      clear_buffer();
    }
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() {
    //    Serial.printf("draw_frame_duration: %d\n", (millis() - start_time));
    effect_.queue_frame();
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Reference to the CRGB pixel.
   */
  inline CRGB& operator()(int x, int y) {
    return effect_.bufs_[effect_.cur_][XY(x, y)];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the CRGB pixel.
   */
  inline CRGB& operator()(int xy) {
    return effect_.bufs_[effect_.cur_][xy];
  }

  /**
   * @brief Gets the width of the canvas.
   * @return The width.
   */
  const int width() { return effect_.width(); }

  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    memset(effect_.bufs_[effect_.cur_], 0, sizeof(CRGB) * effect_.width_ * H);
  }

private:
  Effect& effect_; /**< Reference to the owning Effect instance. */
  size_t start_time; /**< Tracks frame drawing duration (debug). */
};

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
      FastLED.showColor(CRGB::Black);
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
