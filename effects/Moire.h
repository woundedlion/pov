/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include <functional>
#include <memory>

#include <map>
#include "core/effects_engine.h"

template <int W, int H> class Moire : public Effect {
public:
  FLASHMEM Moire()
      : Effect(W, H),
        base_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::BELL),
        base_next_palette(GradientShape::CIRCULAR,
                          HarmonyType::SPLIT_COMPLEMENTARY,
                          BrightnessProfile::BELL),
        int_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                    BrightnessProfile::CUP),
        int_next_palette(GradientShape::CIRCULAR,
                         HarmonyType::SPLIT_COMPLEMENTARY,
                         BrightnessProfile::CUP),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  bool show_bg() const override { return false; }

#ifdef __EMSCRIPTEN__
  void init() override {
#else
  FLASHMEM void init() {
#endif
    // Moire relies on DistortedRing parsing loops.
    // It allocates ~9.2KB for the Ring points and ~8.1KB for rasterization
    // caches. Give it 32KB of Scratch A to safely handle iterations, and 16KB
    // of Scratch B for AA filters.
    configure_arenas(GLOBAL_ARENA_SIZE - 48 * 1024, 32 * 1024, 16 * 1024);

    params.density = W <= 96 ? 10.0f : 45.0f;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Density", &params.density, 1.0f, 100.0f);
    registerParam("Amp", &params.amp, 0.0f, 1.0f);

    timeline
        .add(0, Animation::PeriodicTimer(80, [this](auto &) { color_wipe(); }))
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300,
                                       ease_mid, true))
        .add(0, Animation::Transition(rotation, 2 * PI_F, 160, ease_mid, false,
                                      true)
                    .then([this]() { rotation = 0.0f; }))
        .add(0,
             Animation::Mutation(params.amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f),
                                 160, ease_mid, true));
  }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Calculate rotations (matching inv_transform/transform logic conceptually)
    Quaternion q_base =
        make_rotation(-X_AXIS, rotation) * make_rotation(-Z_AXIS, rotation);
    Quaternion q_int =
        make_rotation(X_AXIS, rotation) * make_rotation(Z_AXIS, rotation);

    draw_layer(canvas, q_base, base_palette);
    draw_layer(canvas, q_int, int_palette);
  }

private:
  void color_wipe() {
    base_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    int_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);

    timeline.add(
        0, Animation::ColorWipe(base_palette, base_next_palette, 80, ease_mid));
    timeline.add(
        0, Animation::ColorWipe(int_palette, int_next_palette, 80, ease_mid));
  }

  void draw_layer(Canvas &canvas, Quaternion layer_rotation,
                  const GenerativePalette &pal) {
    int count = static_cast<int>(std::ceil(params.density));
    for (int i = 0; i <= count; ++i) {
      float t = static_cast<float>(i) / count;
      float r = t * 2.0f;

      Basis basis = make_basis(layer_rotation, Z_AXIS);
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        Color4 c = pal.get(f.v0);
        c.alpha *= params.alpha;
        f.color = c;
      };

      Plot::DistortedRing::draw<W, H>(
          filters, canvas, basis, r,
          sin_wave(-params.amp, params.amp, 4.0f, 0.0f), fragment_shader,
          rotation);
    }
  }

  struct Params {
    float alpha = 0.2f;
    float density = 10.0f;
    float amp = 0.0f;
  } params;

  float rotation = 0.0f;

  GenerativePalette base_palette;
  GenerativePalette base_next_palette;
  GenerativePalette int_palette;
  GenerativePalette int_next_palette;

  Orientation<W> orientation;
  Timeline<W> timeline;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Moire)
