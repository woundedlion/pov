/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "effects_engine.h"

/**
 * @brief Intermediate base class for per-pixel "shader" effects.
 * @details Provides the full draw_frame loop with 4x SSAA, timeline,
 * dual orientation (view + global), noise generator, and palette.
 * Leaf classes only implement transform() and sample().
 *
 * @tparam W Width of the LED display.
 * @tparam H Height of the LED display.
 */
template <int W, int H> class ShaderEffect : public Effect {
public:
  /**
   * @brief Constructs a ShaderEffect.
   * @param palette The generative palette for color output.
   */
  ShaderEffect(const GenerativePalette &palette)
      : Effect(W, H), myPalette(palette) {
    persist_pixels = false;
    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));
    timeline.add(0, Animation::RandomWalk<W>(global_orientation, UP, noise));
  }

  bool show_bg() const override { return true; }

  /**
   * @brief Map a world-space unit vector to a complex sample coordinate.
   * @param v World-space unit vector (already oriented by global_orientation).
   * @return Complex coordinate in the sampling domain.
   */
  virtual Complex transform(const Vector &v) const = 0;

  /**
   * @brief Evaluate the pattern at a sample point.
   * @param z Complex coordinate from transform().
   * @param warp_x Noise-based warp offset (x).
   * @param warp_y Noise-based warp offset (y).
   * @param t Scaled time value.
   * @return Pattern intensity, typically in [-1, 1].
   */
  virtual float sample(const Complex &z, float warp_x, float warp_y,
                       float t) const = 0;

  /**
   * @brief Compute noise-based warp offsets for a sample point.
   * @details Default implementation produces zero warp. Override to add
   * noise-driven displacement.
   */
  virtual void calc_warp(const Complex &z, float t, float &warp_x,
                         float &warp_y) const {
    warp_x = 0.0f;
    warp_y = 0.0f;
  }

  /**
   * @brief Return the effective time for this frame.
   * @details Default returns raw timeline frame count. Override to apply
   * effect-specific time scaling (e.g. multiply by a speed parameter).
   */
  virtual float get_time() const { return static_cast<float>(timeline.t); }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float t = get_time();

    constexpr float eps = 0.5f;
    constexpr float offsets_x[4] = {eps, -eps, eps, -eps};
    constexpr float offsets_y[4] = {eps, eps, -eps, -eps};

    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);

    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {

        // Calc noise warp once per pixel at the center
        Vector true_center_v = pixel_to_vector<W, H>(x, y);
        Complex center_z = transform(true_center_v);

        float warp_x, warp_y;
        calc_warp(center_z, t, warp_x, warp_y);

        float total_pattern = 0.0f;

        // 4x SSAA loop
        for (int i = 0; i < 4; ++i) {
          float px = static_cast<float>(x) + offsets_x[i];
          float py = static_cast<float>(y) + offsets_y[i];

          // Direct trig to bypass Spherical struct overhead
          float theta = (px * 2.0f * PI_F) / w_float;
          float phi = (py * PI_F) / h_virt_minus_1;
          float sin_phi = sinf(phi);
          Vector true_v(sin_phi * cosf(theta), cosf(phi),
                        sin_phi * sinf(theta));

          Complex z = transform(true_v);
          total_pattern += sample(z, warp_x, warp_y, t);
        }

        float avg_pattern = total_pattern * 0.25f;
        float normalized_pattern = (avg_pattern + 1.0f) * 0.5f;
        canvas(x, y) = myPalette.get(normalized_pattern).color;
      }
    }
  }

protected:
  Timeline<W> timeline;
  Orientation<W> orientation;
  Orientation<W> global_orientation;
  FastNoiseLite noise;
  GenerativePalette myPalette;
};
