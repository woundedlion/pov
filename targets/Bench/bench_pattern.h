/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file bench_pattern.h
 * @brief Stationary bench test pattern: one colour across the whole canvas,
 *        cycling through the primaries and white.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Whole-canvas colour cycle that holds on red, green, blue and white.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Uniform in x, so a rotor at rest shows a steady colour rather than a
 * smear of columns, and uniform in y, so every segment of every arm shows the
 * same colour at the same instant. The phase is the frame counter and each
 * board draws exactly one frame per display window from a shared epoch, so a
 * board whose sync has slipped reads directly as a mismatched arm. The holds
 * on full primaries expose a dead channel or a swapped LED colour order, which
 * a continuous sweep hides.
 */
template <int W, int H> class BenchPattern : public Effect {
public:
  /**
   * @brief Display windows per second at the rotor speed this pattern's
   *        timings assume; the sketch checks it against its own RPM.
   */
  static constexpr int FRAMES_PER_SECOND = 16;
  /** @brief Frames each key colour is held. */
  static constexpr int HOLD_FRAMES = 2 * FRAMES_PER_SECOND;
  /** @brief Frames spent ramping from one key colour to the next. */
  static constexpr int RAMP_FRAMES = 4 * FRAMES_PER_SECOND;

  /** @brief One colour the cycle holds on, in 8-bit sRGB. */
  struct Key {
    uint8_t r, g, b;
  };
  /** @brief The held colours, in cycle order. */
  static constexpr Key KEYS[] = {
      {255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 255}};
  static constexpr int KEY_COUNT = static_cast<int>(std::size(KEYS));
  static constexpr int STEP_FRAMES = HOLD_FRAMES + RAMP_FRAMES;
  /** @brief Frames in one full pass over KEYS. */
  static constexpr uint32_t CYCLE_FRAMES =
      static_cast<uint32_t>(KEY_COUNT) * STEP_FRAMES;

  /**
   * @brief Fraction of full LED drive the pattern runs at, in percent.
   * @details Scaled in linear space, so it is also the fraction of the strip
   * current a white hold draws. Override with -D HS_BENCH_BRIGHTNESS_PERCENT.
   */
  static constexpr uint32_t BRIGHTNESS_PERCENT =
#ifdef HS_BENCH_BRIGHTNESS_PERCENT
      HS_BENCH_BRIGHTNESS_PERCENT;
#else
      25;
#endif
  static_assert(BRIGHTNESS_PERCENT > 0 && BRIGHTNESS_PERCENT <= 100,
                "HS_BENCH_BRIGHTNESS_PERCENT must be in 1..100");

  /**
   * @brief Constructs the pattern at the W x H canvas resolution.
   */
  HS_COLD_MEMBER BenchPattern() : Effect(W, H) {}

  /**
   * @brief Fills this board's segment band with the current cycle colour.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    const Pixel colour = colour_at(frame++);
    const ClipRegion &band = canvas.clip();
    Pixel *const buffer = canvas.data();
    for (int y = band.y_start; y < band.y_end; ++y)
      std::fill_n(buffer + y * W + band.x_start, band.x_end - band.x_start,
                  colour);
  }

  /**
   * @brief The colour the cycle shows on one frame.
   * @param f Frames since construction.
   * @return The linear-space pixel every LED on every board shows that frame.
   */
  static Pixel colour_at(uint32_t f) {
    const uint32_t t = f % CYCLE_FRAMES;
    const int step = static_cast<int>(t / STEP_FRAMES);
    const int held = static_cast<int>(t % STEP_FRAMES);
    const Key from = KEYS[step];
    if (held < HOLD_FRAMES)
      return dimmed(from.r, from.g, from.b);
    const Key to = KEYS[(step + 1) % KEY_COUNT];
    const int n = held - HOLD_FRAMES;
    return dimmed(ramp(from.r, to.r, n), ramp(from.g, to.g, n),
                  ramp(from.b, to.b, n));
  }

private:
  static uint8_t ramp(uint8_t from, uint8_t to, int n) {
    return static_cast<uint8_t>(
        from +
        (static_cast<int>(to) - static_cast<int>(from)) * n / RAMP_FRAMES);
  }

  static Pixel dimmed(uint8_t r, uint8_t g, uint8_t b) {
    const Pixel full{CRGB(r, g, b)};
    return Pixel(scale(full.r), scale(full.g), scale(full.b));
  }

  static uint16_t scale(uint16_t v) {
    return static_cast<uint16_t>(static_cast<uint32_t>(v) * BRIGHTNESS_PERCENT /
                                 100u);
  }

  uint32_t frame = 0;
};
