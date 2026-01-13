#pragma once
/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <vector>
#include <cmath>
#include "../effects_engine.h"

template <int W>
class PetalFlow : public Effect {
public:
  PetalFlow() :
    Effect(W),
    palette(
      { 0.029f, 0.029f, 0.029f },
      { 0.500f, 0.500f, 0.500f },
      { 0.461f, 0.461f, 0.461f },
      { 0.539f, 0.701f, 0.809f }
    ),
    filters(
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>()
    )
  {
    persist_pixels = false;

    timeline
      .add(0, Rotation<W>(orientation, Y_AXIS, PI_F / 4.0f, 160, ease_mid, true))
      .add(0, Mutation(twist_factor, sin_wave(2.0f, 2.5f, 1.0f, 0.0f), 160, ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    float time = (millis() / 1000.0f) * (speed * 0.015f);
    float current_spacing = spacing;
    float loop_val = time / current_spacing;
    int loop_count = static_cast<int>(std::floor(loop_val));
    float loop_t = (loop_val - loop_count) * current_spacing;

    draw_petals(canvas, loop_count, loop_t);
  }

private:
  void draw_petals(Canvas& canvas, int loop_count, float loop_t) {
    const float log_min = -3.75f;
    const float log_max = 3.75f;
    const float current_spacing = spacing;

    int min_k = static_cast<int>(std::floor(log_min / current_spacing)) - 1;
    int max_k = static_cast<int>(std::ceil(log_max / current_spacing)) + 1;
    float progress = loop_t / current_spacing;

    auto get_shift = [](float angle_normal) -> float {
      return 0.6f * std::abs(sinf(3.0f * PI_F * angle_normal));
      };

    const int num_samples = W;
    const float step = 2.0f * PI_F / num_samples;

    for (int k = min_k; k <= max_k; ++k) {
      float log_r = k * current_spacing;
      float effective_log_r = log_r + loop_t;

      float dist = std::abs(effective_log_r);
      float opacity = 1.0f;
      if (dist > 2.5f) {
        opacity = std::max(0.0f, 1.0f - (dist - 2.5f) / 1.0f);
      }
      if (opacity <= 0.01f) continue;

      float twist_angle = (k + progress) * twist_factor;
      Points points;

      // Generate Ring directly in Complex Plane
      for (int i = 0; i < num_samples; ++i) {
        float t = static_cast<float>(i) / num_samples;
        float theta = i * step;

        // Apply Petal Wiggle to the Radius (rho)
        float rho = effective_log_r + get_shift(t);
        float final_theta = theta + twist_angle;

        // Convert to Complex Number z = e^(rho + i*theta)
        float R = expf(rho);
        Complex z(
          R * cosf(final_theta),
          R * sinf(final_theta)
        );

        points.push_back(orientation.orient(inv_stereo(z)));
      }

      // Rasterize & Color
      int color_index = (k - loop_count) + 10000;
      float hue = wrap(color_index * 0.13f, 1.0f);

      Plot<W>::rasterize(filters, canvas, points, [&](const Vector&, float) {
        Color4 res = palette.get(hue);
        res.alpha *= opacity;
        return res;
        }, true);
    }
  }

  float alpha = 0.2f;
  float spacing = 0.3f;
  float twist_factor = 2.15f;
  float speed = 8.0f;

  ProceduralPalette palette;
  Orientation orientation;
  Pipeline<W,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;
  Timeline timeline;
};