/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Stack of evenly spaced plotted rings displaced by a 3D noise field.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings share one axis and are spaced evenly in colatitude across the
 * whole sphere. Each ring vertex is displaced along the stack axis by a
 * two-octave OpenSimplex field (independent spatial scale per octave) sampled
 * at the vertex's world-space position, so the displacement is coherent across
 * rings and drifts as the field animates. The displacement direction is
 * uniform across the whole sphere (not mirrored per hemisphere). Fragments are
 * shaded from a circular analogous palette that spins across the stack, with
 * hue rotated proportionally to the local displacement magnitude. Orientation
 * random-walks over time.
 */
template <int W, int H> class NoiseRings : public Effect {
public:
  /**
   * @brief Builds the effect with its palette and ring-stack axis.
   */
  HS_COLD_MEMBER NoiseRings()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        ring_palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                     BrightnessProfile::FLAT),
        normal(X_AXIS) {}

  /**
   * @brief Registers params, seeds the noise field, and builds the timeline.
   */
  void init() override {
    ring_baked.bake(persistent_arena, ring_palette);

    field_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    field_noise.SetSeed(hs::rand_int(0, 65536));
    field_noise.SetFrequency(1.0f);

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Rings", &params.num_rings, 1.0f, 72.0f);
    register_param("Amplitude", &params.amplitude, 0.0f, 0.8f);
    register_param("Scale 1", &params.scale1, 0.5f, 4.0f);
    register_param("Scale 2", &params.scale2, 0.5f, 8.0f);
    register_param("Hue Rotate", &params.hue_scale, 0.0f, 3.0f);
    register_param("Flow Speed", &params.flow_speed, 0.0f, 0.15f);

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->draw_fn(canvas, opacity);
                        },
                        -1, 24, ease_linear, 0, ease_linear));

    timeline.add(0, Animation::RandomWalk<W>(orientation, normal, walk_noise));
  }

  /**
   * @brief Advances the noise field and renders the timeline for one frame.
   */
  void draw_frame() override {
    field_time += params.flow_speed;
    color_spin = wrap_t(color_spin + COLOR_SPIN_RATE);
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /**
   * @brief Draws all rings for this frame.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details Per ring, the two-octave displacement is baked once per azimuth
   * column into shift_lut (sampled in the antipode work frame so the field
   * tracks the plotted vertices), then read by both the geometry shift and the
   * per-fragment hue rotation. get_antipode mirrors the polar frame past the
   * equator (radius > 1); hemi_sign negates the far-side shift so a noise crest
   * displaces every ring the same world direction along the stack axis. The
   * fragment hue is rotated by the local |displacement| scaled by Hue Rotate.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      auto res = get_antipode(basis, radius);
      Basis work = res.first;
      float theta_eq = res.second * (PI_F / 2.0f);
      float cos_eq = cosf(theta_eq);
      float sin_eq = sinf(theta_eq);
      float hemi_sign = radius > 1.0f ? -1.0f : 1.0f;

      for (int x = 0; x <= W; ++x) {
        float angle = 2.0f * PI_F * x / W;
        Vector p = (work.v * cos_eq) +
                   ((work.u * cosf(angle)) + (work.w * sinf(angle))) * sin_eq;
        float n1 = field_noise.GetNoise(p.x * params.scale1, p.y * params.scale1,
                                        p.z * params.scale1 + field_time);
        float n2 = field_noise.GetNoise(p.x * params.scale2 + OCTAVE2_OFFSET,
                                        p.y * params.scale2,
                                        p.z * params.scale2 + field_time);
        float n = (n1 + OCTAVE2_WEIGHT * n2) / (1.0f + OCTAVE2_WEIGHT);
        shift_lut[x] = params.amplitude * hemi_sign * n;
      }

      auto sample_lut = [this](float t) {
        float x = wrap_t(t) * W;
        int j = static_cast<int>(x);
        return shift_lut[j] + (x - j) * (shift_lut[j + 1] - shift_lut[j]);
      };

      Color4 ring_color =
          ring_baked.get(wrap_t((i + 0.5f) / n_rings + color_spin));
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = hue_rotate(ring_color,
                             std::abs(sample_lut(f.v0)) * params.hue_scale);
        f.color.alpha *= opacity * params.alpha;
      };

      Plot::DistortedRing::draw<W, H>(filters, canvas, basis, radius, sample_lut,
                                      fragment_shader);
    }
  }

private:
  FastNoiseLite walk_noise;
  FastNoiseLite field_noise;
  Timeline timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  GenerativePalette ring_palette;
  BakedPalette ring_baked;
  Vector normal;
  Orientation<> orientation;

  static constexpr float COLOR_SPIN_RATE = 0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float OCTAVE2_WEIGHT = 0.5f;     /**< Detail-octave weight relative to the base octave. */
  static constexpr float OCTAVE2_OFFSET = 50.0f;    /**< Spatial offset decorrelating octave 2 from octave 1 at equal scales. */

  float field_time = 0.0f; /**< Noise-field time offset, advanced by Flow Speed each frame. */
  float color_spin = 0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  float shift_lut[W + 1] = {}; /**< Per-ring displacement at azimuth resolution; entry W repeats entry 0 for seamless lerp. */

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;       /**< Overall ring opacity multiplier in [0, 1]. */
    float num_rings = 8.0f;   /**< Number of evenly spaced rings (truncated to int when drawn). */
    float amplitude = 0.25f;  /**< Peak polar displacement (radians) scaling the noise field. */
    float scale1 = 1.5f;      /**< Spatial frequency of octave 1 over the unit sphere. */
    float scale2 = 3.0f;      /**< Spatial frequency of octave 2 over the unit sphere. */
    float hue_scale = 1.0f;   /**< Hue rotation (turns) per radian of displacement magnitude. */
    float flow_speed = 0.03f; /**< Noise-field time advance per frame. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(NoiseRings)
