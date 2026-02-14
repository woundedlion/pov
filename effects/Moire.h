/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W, int H>
class Moire : public Effect {
public:
  Moire() :
    Effect(W, H),
    base_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
    base_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
    int_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
    int_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
    filters(
      Filter::World::Orient<W>(orientation),
      Filter::Screen::AntiAlias<W, H>()
    )
  {
    persist_pixels = false;

    timeline
      .add(0, Animation::PeriodicTimer(80, [this](auto&) { color_wipe(); }))
      .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300, ease_mid, true))
      .add(0, Animation::Transition(rotation, 2 * PI_F, 160, ease_mid, false, true)
        .then([this]() { rotation = 0.0f; }))
      .add(0, Animation::Mutation(amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f), 160, ease_mid, true)
      );
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    draw_layer(canvas, [this](const Vector& p) { return inv_transform(p); }, base_palette);
    draw_layer(canvas, [this](const Vector& p) { return transform(p); }, int_palette);
  }

private:

  void color_wipe() {
    base_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);
    int_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);

    timeline.add(0, Animation::ColorWipe(base_palette, base_next_palette, 80, ease_mid));
    timeline.add(0, Animation::ColorWipe(int_palette, int_next_palette, 80, ease_mid));
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

    void draw_layer(Canvas& canvas, TransformFn auto trans_fn, const GenerativePalette& pal) {
    int count = static_cast<int>(std::ceil(density));
    for (int i = 0; i <= count; ++i) {
      float t = static_cast<float>(i) / count;
      float r = t * 2.0f;

      Basis basis = make_basis(Quaternion(), Z_AXIS);
      auto fragment_shader = [&](const Vector&, Fragment& f) { 
          f.color = pal.get(t);
      };
      Plot::DistortedRing::draw<W, H>(filters, canvas, basis, r,
        sin_wave(-amp, amp, 4.0f, 0.0f), 
        fragment_shader);
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

  Orientation<W> orientation;
  Timeline<W> timeline;

  Pipeline<W, H,
    Filter::World::Orient<W>,
    Filter::Screen::AntiAlias<W, H>
  > filters;
};