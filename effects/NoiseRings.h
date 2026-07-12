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
 * whole sphere. Each ring vertex is displaced along the stack axis by the
 * product of two OpenSimplex octaves (independent spatial scale per octave)
 * sampled at the vertex's world-space position: octave 1 envelopes octave 2,
 * so perturbations bunch where the envelope is strong and vanish where it
 * crosses zero. Displacement is coherent across rings, drifts as the field
 * animates, and its direction is
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
   * @details Per ring, the two-octave displacement product is baked once per
   * azimuth column into shift_lut (sampled in the antipode work frame so the field
   * tracks the plotted vertices), then read by both the geometry shift and the
   * per-fragment hue rotation. The LUT resolution is adaptive: enough samples for
   * the finer octave along the ring's actual circumference, so small rings and low
   * scales bake far fewer noise columns. Under a partial clip (segmented drivers),
   * rings whose displaced band cannot touch the clip are skipped whole — they
   * could contribute no fragment, so the skip is output-identical.
   * get_antipode mirrors the polar frame past the equator (radius > 1); hemi_sign
   * negates the far-side shift so a noise crest displaces every ring the same
   * world direction along the stack axis. The fragment hue is rotated by the
   * local |displacement| scaled by Hue Rotate.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing stroke width + AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float max_scale = std::max(params.scale1, params.scale2);

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      auto res = get_antipode(basis, radius);
      Basis work = res.first;
      float theta_eq = res.second * (PI_F / 2.0f);
      if (try_cull &&
          !ring_may_touch_clip(work.v, theta_eq, params.amplitude + pad))
        continue;
      float cos_eq = cosf(theta_eq);
      float sin_eq = sinf(theta_eq);
      float hemi_sign = radius > 1.0f ? -1.0f : 1.0f;

      // Nyquist-safe column count: LUT_SAMPLES_PER_UNIT samples per noise-space
      // unit along the ring's circumference in the finer octave.
      int lut_n = hs::clamp(
          static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * PI_F *
                                 max_scale * sin_eq)),
          LUT_MIN_SAMPLES, W);
      for (int x = 0; x < lut_n; ++x) {
        float angle = 2.0f * PI_F * x / lut_n;
        Vector p = (work.v * cos_eq) +
                   ((work.u * cosf(angle)) + (work.w * sinf(angle))) * sin_eq;
        float n1 = field_noise.GetNoise(p.x * params.scale1, p.y * params.scale1,
                                        p.z * params.scale1 + field_time);
        float n2 = field_noise.GetNoise(p.x * params.scale2 + OCTAVE2_OFFSET,
                                        p.y * params.scale2,
                                        p.z * params.scale2 + field_time);
        shift_lut[x] = params.amplitude * hemi_sign * n1 * n2;
      }
      shift_lut[lut_n] = shift_lut[0];

      auto sample_lut = [this, lut_n](float t) {
        float x = wrap_t(t) * lut_n;
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
  /**
   * @brief Conservative test: can the displaced ring reach the current clip?
   * @param axis Ring axis in the antipode work frame (unit vector).
   * @param theta_eq Ring colatitude about @p axis (radians).
   * @param delta Displacement bound plus stroke/AA pad (radians).
   * @return False only when no fragment of the ring can land inside the clip's
   * render region; true is always safe.
   * @details The displaced ring lies in the band theta_eq ± delta around the
   * axis. Rows: the band's polar range about the display's Y axis is
   * [beta - t2, beta + t2] (beta = axis colatitude). Columns: a band that
   * reaches either display pole spans all longitudes; otherwise its longitude
   * half-width about the axis longitude is asin(sin t2 / sin beta) (spherical
   * cap extent), compared against the clip's column wedge with the clip margin
   * plus one pixel of slack.
   */
  bool ring_may_touch_clip(const Vector &axis, float theta_eq,
                           float delta) const {
    const ClipRegion &cr = clip();
    float t2 = std::min(theta_eq + delta, PI_F);
    float beta = acosf(hs::clamp(axis.y, -1.0f, 1.0f));

    float phi_lo = std::max(beta - t2, 0.0f);
    float phi_hi = std::min(beta + t2, PI_F);
    if (phi_to_y<H>(phi_hi) < cr.render_y_start() ||
        phi_to_y<H>(phi_lo) >= cr.render_y_end())
      return false;

    if (cr.x_start == 0 && cr.x_end == cr.w)
      return true;
    if (beta <= t2 || PI_F - beta <= t2)
      return true;
    float dlam = asinf(hs::clamp(sinf(t2) / sinf(beta), 0.0f, 1.0f));
    float lam_v = atan2f(axis.z, axis.x);
    float width_px = static_cast<float>(cr.x_end - cr.x_start);
    float half_w =
        (width_px * 0.5f + cr.margin + 1.0f) * (2.0f * PI_F) / cr.w;
    float lam_c = (cr.x_start + width_px * 0.5f) * (2.0f * PI_F) / cr.w;
    float d = std::fabs(wrap_t((lam_v - lam_c) / (2.0f * PI_F) + 0.5f) - 0.5f) *
              (2.0f * PI_F);
    return d <= dlam + half_w;
  }

  FastNoiseLite walk_noise;
  FastNoiseLite field_noise;
  Timeline timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  GenerativePalette ring_palette;
  BakedPalette ring_baked;
  Vector normal;
  Orientation<> orientation;

  static constexpr float COLOR_SPIN_RATE = 0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float OCTAVE2_OFFSET = 50.0f;    /**< Spatial offset decorrelating octave 2 from octave 1 at equal scales. */
  static constexpr float LUT_SAMPLES_PER_UNIT = 4.0f; /**< Bake columns per noise-space unit of ring circumference. */
  static constexpr int LUT_MIN_SAMPLES = 16;          /**< Bake-column floor for tiny/low-scale rings. */

  float field_time = 0.0f; /**< Noise-field time offset, advanced by Flow Speed each frame. */
  float color_spin = 0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  float shift_lut[W + 1] = {}; /**< Per-ring displacement at azimuth resolution; entry W repeats entry 0 for seamless lerp. */

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;       /**< Overall ring opacity multiplier in [0, 1]. */
    float num_rings = 55.0f;  /**< Number of evenly spaced rings (truncated to int when drawn). */
    float amplitude = 0.25f;  /**< Peak polar displacement (radians) scaling the noise field. */
    float scale1 = 1.5f;      /**< Spatial frequency of the envelope octave; its zero regions leave rings undisturbed. */
    float scale2 = 3.0f;      /**< Spatial frequency of the detail octave, scaled by the envelope octave. */
    float hue_scale = 2.0f;   /**< Hue rotation (turns) per radian of displacement magnitude. */
    float flow_speed = 0.03f; /**< Noise-field time advance per frame. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(NoiseRings)
