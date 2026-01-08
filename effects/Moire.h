/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W>
class Moire : public Effect {
public:
  Moire() :
    Effect(W),
    base_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
    base_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
    int_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
    int_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
    filters(
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>()
    )
  {
    persist_pixels = false;

    timeline
      .add(0, PeriodicTimer(80, [this](auto&) { color_wipe(); }))
      .add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300, ease_mid, true))
      .add(0, Transition(rotation, 2 * PI_F, 160, ease_mid, false, true)
        .then([this]() { rotation = 0.0f; }))
      .add(0, Mutation(amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f), 160, ease_mid, true)
      );
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    dots.clear();

    draw_layer([this](const Vector& p) { return inv_transform(p); }, base_palette);
    draw_layer([this](const Vector& p) { return transform(p); }, int_palette);

    plot_dots<W>(dots, filters, canvas, 0, alpha);
  }

private:

  void color_wipe() {
    base_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);
    int_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);

    timeline.add(0, ColorWipe(base_palette, base_next_palette, 80, ease_mid));
    timeline.add(0, ColorWipe(int_palette, int_next_palette, 80, ease_mid));
  }

  Vector transform(Vector p) {
    Quaternion q_z = make_rotation(Z_AXIS, rotation);
    Quaternion q_x = make_rotation(X_AXIS, rotation);
    p = rotate(p, q_z);
    p = rotate(p, q_x);
    return p;
  }

  Vector inv_transform(Vector p) {
    Quaternion q_nx = make_rotation(-X_AXIS, rotation);
    Quaternion q_nz = make_rotation(-Z_AXIS, rotation);
    p = rotate(p, q_nx);
    p = rotate(p, q_nz);
    return p;
  }

  void draw_layer(TransformFn auto trans_fn, const GenerativePalette& pal) {
    int count = static_cast<int>(std::ceil(density));
    for (int i = 0; i <= count; ++i) {
      float t = static_cast<float>(i) / count;
      float r = t * 2.0f;

      Points points;
      sample_fn<W>(points, Quaternion(1, 0, 0, 0), Z_AXIS, r,
        sin_wave(-amp, amp, 4.0f, 0.0f));

      Points transformed_points;
      for (const auto& p : points) {
        transformed_points.push_back(trans_fn(p));
      }

      rasterize<W>(dots, transformed_points,
        [&](const Vector&, float) { return pal.get(t); }, true);
    }
  }

  float alpha = 0.2f;
  float density = 10.0f;
  float rotation = 0.0f;
  float amp = 0.0f;

  GenerativePalette base_palette;
  GenerativePalette base_next_palette;
  GenerativePalette int_palette;
  GenerativePalette int_next_palette;

  Orientation orientation;
  Timeline timeline;
  Dots dots;

  Pipeline<W,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;
};