/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/effects_engine.h"

// Moire interference effect: two counter-rotating stacks of concentric
// DistortedRing layers beat against each other to produce a shimmering moire
// pattern, tinted by two wipe-animated generative palettes.
template <int W, int H> class Moire : public Effect {
public:
  // Builds the two palette pairs with deliberately contrasting brightness
  // profiles (BELL base vs CUP interference) so the layers stay visually
  // distinct, and wires the world-orient + screen anti-alias filter pipeline.
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

  // Rings are drawn over a transparent background.
  bool show_bg() const override { return false; }

  // One-time setup: size the scratch arenas, bake the palette LUTs, register
  // the user params, and arm the rotation / wipe / amplitude animations.
  void init() override {
    // Each frame draws many concentric DistortedRing layers, each of which
    // binds W+2 ring points and rasterization buffers into Scratch A. Give it
    // 32KB of Scratch A for those, and 16KB of Scratch B for the AA filter.
    configure_arenas(GLOBAL_ARENA_SIZE - 48 * 1024, 32 * 1024, 16 * 1024);

    // Allocate the two palette LUTs (rebaked each frame in draw_frame()).
    base_baked.bake(persistent_arena, base_palette);
    int_baked.bake(persistent_arena, int_palette);

    params.density = W <= 96 ? 10.0f : 45.0f;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Density", &params.density, 1.0f, 100.0f);
    registerParam("Amp", &params.amp, 0.0f, 1.0f);
    // Amp is driven by the Mutation below; flag it so the GUI auto-pauses the
    // animation when the user grabs the slider.
    markAnimated("Amp");

    timeline
        .add(0, Animation::PeriodicTimer(80, [this](auto &) { color_wipe(); }))
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300,
                                       ease_mid, true))
        .add(0, Animation::Transition(rotation, 2 * PI_F, 160, ease_mid, false,
                                      true)
                    .then([this]() { rotation = 0.0f; }))
        .add(0,
             Animation::Mutation(params.amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f),
                                 160, ease_mid, true, &anims_paused_));
  }

  // Per-frame: advance animations, refresh the palette LUTs, then draw the two
  // counter-rotating ring layers so they interfere.
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Refresh the LUTs from the (continuously wipe-animated) source palettes.
    // Once per frame, vs. a full OKLCH interpolation per ring fragment.
    base_baked.rebake(base_palette);
    int_baked.rebake(int_palette);

    // Two layers counter-rotate (opposite-signed axes) so their rings beat
    // against each other, producing the moire interference.
    Quaternion q_base =
        make_rotation(-X_AXIS, rotation) * make_rotation(-Z_AXIS, rotation);
    Quaternion q_int =
        make_rotation(X_AXIS, rotation) * make_rotation(Z_AXIS, rotation);

    draw_layer(canvas, q_base, base_baked);
    draw_layer(canvas, q_int, int_baked);
  }

private:
  // Picks fresh target palettes and schedules a wipe to them; called
  // periodically to keep the colors drifting.
  void color_wipe() {
    // The two counter-rotating layers are deliberately contrasting: the base
    // layer starts BELL, the interference layer CUP (see the constructor). Keep
    // that inter-layer contrast through the wipe by giving the new palettes
    // opposite brightness gradients — otherwise both layers wipe to the same
    // profile and the moire beat washes out to incidental hue noise.
    base_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    int_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::DESCENDING);

    timeline.add(
        0, Animation::ColorWipe(base_palette, base_next_palette, 80, ease_mid));
    timeline.add(
        0, Animation::ColorWipe(int_palette, int_next_palette, 80, ease_mid));
  }

  // Draws one stack of `density` concentric DistortedRings under the given
  // layer rotation, radius stepping from 0 to 2, shaded by `pal`. The rings are
  // wobbled by a sine wave whose amplitude is the animated `amp` param.
  void draw_layer(Canvas &canvas, Quaternion layer_rotation,
                  const BakedPalette &pal) {
    int count = static_cast<int>(std::ceil(params.density));
    // Loop-invariant: layer_rotation and Z_AXIS do not change across rings.
    Basis basis = make_basis(layer_rotation, Z_AXIS);
    auto fragment_shader = [&](const Vector &, Fragment &f) {
      Color4 c = pal.get(f.v0);
      c.alpha *= params.alpha;
      f.color = c;
    };
    for (int i = 0; i <= count; ++i) {
      float r = static_cast<float>(i) / count * 2.0f;

      Plot::DistortedRing::draw<W, H>(
          filters, canvas, basis, r,
          sin_wave(-params.amp, params.amp, 4.0f, 0.0f), fragment_shader,
          rotation);
    }
  }

  // User-tunable params: per-ring opacity, ring count, and wobble amplitude.
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

  // 256-entry LUTs rebaked from the two source palettes once per frame (cheap:
  // 2 x 256 OKLCH evaluations) so the per-fragment shader does a LUT lookup
  // instead of a full sRGB->OKLCH interpolation on every ring fragment.
  BakedPalette base_baked;
  BakedPalette int_baked;

  Orientation<> orientation;
  Timeline timeline;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Moire)
