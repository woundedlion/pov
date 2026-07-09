/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/engine/engine.h"

/**
 * @brief Moire interference effect: two counter-rotating stacks of concentric
 *        DistortedRing layers beating against each other.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The two stacks produce a shimmering moire pattern, tinted by two
 *          wipe-animated generative palettes.
 */
template <int W, int H> class Moire : public Effect {
public:
  /**
   * @brief Constructs the effect, building the two palette pairs and the
   *        filter pipeline.
   * @details Builds the initial palette pairs with contrasting brightness
   *          profiles (BELL base vs CUP interference) so the layers start
   *          visually distinct; color_wipe() later swaps both to a different
   *          contrasting pair (ASCENDING vs DESCENDING). Wires the world-orient
   *          + screen anti-alias filter pipeline.
   */
  FLASHMEM Moire()
      : Effect(W, H, {.strobe = true}),
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

  /**
   * @brief One-time setup of arenas, palette LUTs, params, and animations.
   * @details Sizes the scratch arenas, bakes the palette LUTs, registers the
   *          user params, and arms the rotation / wipe / amplitude animations.
   */
  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 48 * 1024, 32 * 1024, 16 * 1024);

    base_baked.bake(persistent_arena, base_palette);
    int_baked.bake(persistent_arena, int_palette);

    // Set before register_param, which snapshots *ptr as the reset-to-default:
    // the resolution-dependent default must be in place first.
    params.density = W <= 96 ? 10.0f : 45.0f;

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Density", &params.density, 1.0f, 100.0f);
    register_animated_param("Amp", &params.amp, 0.0f, 1.0f);

    timeline
        .add(0,
             Animation::PeriodicTimer(80, [this](Canvas &) { color_wipe(); },
                                      true))
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300,
                                       ease_linear, true))
        // Ramp `rotation` 0 -> 2pi then snap back to 0. The snap is seamless
        // only because every consumer of `rotation` is exactly 2pi-periodic
        // (make_rotation about fixed axes; draw_layer's wobble is sin_wave at an
        // integer frequency).
        .add(0, Animation::Transition(rotation, 2 * PI_F, 160, ease_linear, false,
                                      true)
                    .then([this]() { rotation = 0.0f; }))
        .add(0,
             Animation::Mutation(params.amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f),
                                 160, ease_linear, true, &anims_paused_));
  }

  /**
   * @brief Renders one frame of the moire effect.
   * @details Advances animations, refreshes the palette LUTs, then draws the
   *          two counter-rotating ring layers so they interfere.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    base_baked.rebake(base_palette);
    int_baked.rebake(int_palette);

    // Counter-rotate the two layers (opposite-signed axes) so their rings beat
    // against each other, producing the moire interference.
    Quaternion q_base =
        make_rotation(-X_AXIS, rotation) * make_rotation(-Z_AXIS, rotation);
    Quaternion q_int =
        make_rotation(X_AXIS, rotation) * make_rotation(Z_AXIS, rotation);

    draw_layer(canvas, q_base, base_baked);
    draw_layer(canvas, q_int, int_baked);
  }

private:
  /**
   * @brief Picks fresh target palettes and schedules a wipe to them.
   * @details Called periodically to keep the colors drifting.
   */
  void color_wipe() {
    // Opposite brightness gradients (ASCENDING vs DESCENDING) keep the two
    // layers visually distinct through the wipe.
    base_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    int_next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::DESCENDING);

    timeline.add(0, Animation::ColorWipe(base_palette, base_next_palette,
                                         WIPE_FRAMES, ease_linear));
    timeline.add(0, Animation::ColorWipe(int_palette, int_next_palette,
                                         WIPE_FRAMES, ease_linear));
  }

  static constexpr int WIPE_FRAMES = 80; /**< Duration of a palette cross-fade ColorWipe, in frames. */

  /**
   * @brief Rings finer than one display row apart collapse on the raster, so
   * the per-layer ring count is capped at one ring per row to skip overdraw at
   * high Density.
   */
  static constexpr int MAX_RINGS = H > 1 ? H : 1;

  /**
   * @brief Draws one stack of concentric DistortedRings for a layer.
   * @param canvas Render target for this layer's rings.
   * @param layer_rotation Quaternion rotation applied to this stack.
   * @param pal Baked palette LUT used to shade the rings.
   * @details Draws `density` rings with radius stepping from 0 to 2, each
   *          wobbled by a sine wave whose amplitude is the animated `amp` param.
   */
  void draw_layer(Canvas &canvas, Quaternion layer_rotation,
                  const BakedPalette &pal) {
    int count = static_cast<int>(std::ceil(params.density));
    if (count > MAX_RINGS)
      count = MAX_RINGS;
    Basis basis = make_basis(layer_rotation, Z_AXIS);
    auto fragment_shader = [&](const Vector &, Fragment &f) {
      Color4 c = pal.get(f.v0);
      c.alpha *= params.alpha;
      f.color = c;
    };
    // i runs 1..count: i == count closes the pattern at r == 2.0; i == 0 would
    // be a degenerate r == 0 ring.
    for (int i = 1; i <= count; ++i) {
      float r = static_cast<float>(i) / count * 2.0f;

      Plot::DistortedRing::draw<W, H>(
          filters, canvas, basis, r,
          sin_wave(-params.amp, params.amp, 4.0f, 0.0f), fragment_shader,
          rotation);
    }
  }

  /**
   * @brief User-tunable parameters for the effect.
   */
  struct Params {
    float alpha = 0.2f;   /**< Per-ring opacity in [0, 1]. */
    float density = 10.0f; /**< Ring count per layer in [1, 100]. Overridden per-resolution in init(), so this default never reaches a live frame. */
    float amp = 0.0f;     /**< Wobble amplitude in [0, 1]. */
  } params;

  float rotation = 0.0f; /**< Current layer rotation angle in radians. */

  GenerativePalette base_palette;      /**< Live source palette for the base layer. */
  GenerativePalette base_next_palette; /**< Wipe target palette for the base layer. */
  GenerativePalette int_palette;       /**< Live source palette for the interference layer. */
  GenerativePalette int_next_palette;  /**< Wipe target palette for the interference layer. */

  /**
   * @brief 256-entry LUT rebaked from base_palette once per frame.
   * @details Cheap (256 OKLCH evaluations) so the per-fragment shader does a
   *          LUT lookup instead of a full sRGB->OKLCH interpolation on every
   *          ring fragment.
   */
  BakedPalette base_baked;
  BakedPalette int_baked; /**< 256-entry LUT rebaked from int_palette once per frame. */

  Orientation<> orientation; /**< World-space orientation shared with the Orient filter. */
  Timeline timeline;         /**< Drives the rotation, wipe, and amplitude animations. */

  /**
   * @brief World-orient + screen anti-alias filter pipeline applied to rings.
   */
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Moire)
