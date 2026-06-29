/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Full-sphere liquid effect over a breathing generative palette.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Domain-warped OpenSimplex noise feeds a cross-coupled sinusoidal
 * pattern, stereographically projected and colored through a breathing
 * generative palette. Presets cycle on a random timer.
 * @note `Flyby` is the sibling stereographic effect. Both build on the same
 *       core primitives — `stereo()`, `stereo_noise_warp()`, `pole_attenuation()`,
 *       `pole_normalize_pattern()` — but the per-effect `project()`, the
 *       `sample()` pattern, and `Params::lerp` are intentionally different (here:
 *       dual orientation + a glitch lens, a cross-coupled pattern, and a
 *       staggered per-field lerp). They are independent effects, NOT hand-synced
 *       copies of one shader: a change here is not expected to propagate to
 *       Flyby. The common per-pixel code (the `STEREO_PATTERN_ARG_LIMIT` clamp
 *       and the pole-attenuated normalize) lives in core.
 */
template <int W, int H> class Liquid2D : public Effect {
public:
  /**
   * @brief Constructs the effect and disables inter-frame pixel persistence.
   * @details The shader writes every pixel each frame, so no persistence is
   * needed.
   */
  FLASHMEM Liquid2D() : Effect(W, H) { persist_pixels = false; }

  /**
   * @brief Initializes sliders, timeline drivers, palette, and preset cycling.
   * @details Registers sliders, seeds the timeline drivers/orientations, bakes
   * the palette, and arms the preset-cycling timer; then loads the first preset.
   */
  void init() override {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Time Speed", &params.time_speed, 0.1f, 5.0f);
    registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);
    registerParam("Cycle Speed", &params.cycle_speed, 0.0f, 1.0f);
    // Flag every preset-driven param so "Pause Animation" lets the user take a
    // slider over.
    for (const char *n :
         {"Warp Scale", "Warp Strength", "Pattern Freq", "Time Speed",
          "Complexity", "Pole Fade", "Cycle Speed"})
      markAnimated(n);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));
    timeline.add(0, Animation::RandomWalk<W>(global_orientation, UP, noise));
    timeline.add(0, Animation::Driver(accumulated_time, &params.time_speed, 1.0f,
                                      false));
    // wrap=false: cycle_phase is wrapped by hand to 2pi in draw_frame; the
    // Driver's [0,1) wrap is the wrong period for fast_sinf.
    timeline.add(0, Animation::Driver(cycle_phase, &params.cycle_speed, 1.0f,
                                      false));

    palette.bake(persistent_arena,
                 GenerativePalette{
                     GradientShape::STRAIGHT, HarmonyType::COMPLEMENTARY,
                     BrightnessProfile::CUP, SaturationProfile::VIBRANT, 75});
    static_palette.bind(&palette, &breathe_mod);

    // Cycle presets every 3-5 seconds via a 2 second lerp.
    timeline.add(0, Animation::RandomTimer(
                        90, 150,
                        [this](Canvas &) {
                          if (animationsPaused())
                            return;
                          presets.next();
                          timeline.add(0, Animation::Lerp(params,
                                                          presets.prev_get(),
                                                          presets.get(), 60,
                                                          ease_in_out_sin,
                                                          &anims_paused_));
                        },
                        true));

    params = presets.get();
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Advances the timeline and shades every pixel for one frame.
   * @details Advances the timeline, maintains the wrapped trig phases, then
   * shades every pixel: projects to the sphere, warps in stereographic space,
   * samples the pattern, attenuates near the poles, and maps through the
   * palette.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // accumulated_time is the unbounded noise-time axis: OpenSimplex2 is
    // aperiodic so it cannot be wrapped without a visible jump.
    float t = accumulated_time;
    // The pattern's fast_sinf/fast_cosf time terms ride phases wrapped to 2pi.
    // dt is the live Time Speed slider directly, not (t - prev_time): the Driver
    // adds that speed each step, so it IS this frame's advance of t. Differencing
    // the unbounded accumulator would lose the increment to ULP after multi-day
    // uptime.
    constexpr float kTwoPi = 2.0f * PI_F;
    float dt = params.time_speed;
    if (!std::isfinite(dt))
      dt = 0.0f;
    sin_phase = fmodf(sin_phase + dt, kTwoPi);
    cos_phase = fmodf(cos_phase + 0.8f * dt, kTwoPi);
    // cycle_phase feeds BreatheModifier's fast_sinf, so wrap to 2pi by hand (the
    // Driver's [0,1) wrap is the wrong domain for a radians consumer).
    cycle_phase = fmodf(cycle_phase, kTwoPi);

    auto shader = [&](const Vector &v) -> Color4 {
      Complex z = project(v);
      float r_sq = z.re * z.re + z.im * z.im;
      Complex w = warp(z, r_sq, t).coords; // only the warped coords feed sample
      float pattern = sample(w, sin_phase, cos_phase);
      float value = pole_normalize_pattern(pattern, r_sq, params.pole_fade);
      return static_palette.get(value);
    };

    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

private:
  /**
   * @brief Projects a sphere direction to the stereographic plane.
   * @param v Unit direction vector on the sphere.
   * @return Complex stereographic coordinate after glitch lens and dual
   * orientation are applied.
   */
  Complex project(const Vector &v) const {
    Vector rv = global_orientation.unorient(v);
    Vector sv = apply_glitch_lens(rv);
    return stereo(orientation.orient(sv));
  }

  /**
   * @brief Warps a stereographic coordinate using noise, faded near the pole.
   * @param z Stereographic coordinate to warp.
   * @param r_sq Squared radius |z|^2, used for pole fading.
   * @param t Noise-time axis (seconds); halved internally for the noise sample.
   * @return Warp result whose coords are the displaced stereographic position.
   */
  StereoWarpResult warp(const Complex &z, float r_sq, float t) const {
    return stereo_noise_warp(z, r_sq, noise, params.warp_scale,
                             params.warp_strength, params.pole_fade, t * 0.5f);
  }

  /**
   * @brief Evaluates the cross-coupled sinusoidal pattern at a warped point.
   * @param w Warped stereographic coordinate to sample.
   * @param sin_phase Wrapped +t time term in [0, 2pi) radians.
   * @param cos_phase Wrapped 0.8*t time term in [0, 2pi) radians.
   * @return Pattern value in [-1, 1], modulated by params.complexity.
   */
  float sample(const Complex &w, float sin_phase, float cos_phase) const {
    // Near the pole |w| -> STEREO_INF, so w*pattern_freq can reach ~2e5 where
    // fast_sinf range reduction bands; clamp the (pole-attenuated) argument.
    float pu = hs::clamp(w.re * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    float pv = hs::clamp(w.im * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    return fast_sinf(pu + params.complexity * fast_sinf(pv + sin_phase)) *
           fast_cosf(pv + params.complexity * fast_cosf(pu - cos_phase));
  }

  /**
   * @brief Applies a trig-free glitch lens to a sphere direction.
   * @param v Unit direction vector on the sphere.
   * @return Transformed direction after mirror, squish, and triple-theta steps;
   * returns the up vector when near the lens axis (R^2 < 1e-6).
   */
  static Vector apply_glitch_lens(Vector v) {
    if (v.y < 0.0f) {
      v.y = -v.y;
      v.z = -v.z;
    }

    float x2 = v.x * v.x;
    float z2 = v.z * v.z;
    float R2 = x2 + z2;

    // Pole guard: on the rotation axis (x≈z≈0) the lens map divides by R², so
    // return the pole direction directly.
    constexpr float MIN_AXIS_RADIUS2 = 1e-6f;
    if (R2 < MIN_AXIS_RADIUS2)
      return Vector(0.0f, 1.0f, 0.0f);

    float inv_R2 = 1.0f / R2;
    float y2 = 2.0f * v.y;

    return Vector(y2 * v.x * (4.0f * x2 * inv_R2 - 3.0f), y2 * v.y - 1.0f,
                  y2 * v.z * (3.0f - 4.0f * z2 * inv_R2));
  }

  Timeline timeline; /**< Drives orientations, time/cycle drivers, and presets. */
  Orientation<> orientation;        /**< Inner per-pixel sphere orientation. */
  Orientation<> global_orientation; /**< Outer whole-sphere orientation. */
  FastNoiseLite noise;              /**< OpenSimplex2 source for warp and walks. */

  BakedPalette palette; /**< 16-bit LUT baked from the generative palette. */
  BreatheModifier breathe_mod{&cycle_phase, 0.15f}; /**< Cycle-phase breathing modulator. */
  StaticPalette<BakedPalette, Coords<BreatheModifier>> static_palette; /**< Palette plus breathe modulation used for shading. */

  /**
   * @brief Tunable per-frame parameters exposed as sliders and lerped between
   * presets.
   * @details The default member initializers double as preset 0.
   */
  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 5.0f;
    float time_speed = 0.1f;
    float complexity = 0.5f;
    float pole_fade = 1.4f;
    float cycle_speed = 0.05f;

    /**
     * @brief Staggered interpolation of this Params from a to b.
     * @param a Source parameter set (value at t = 0).
     * @param b Target parameter set (value at t = 1).
     * @param t Normalized progress in [0, 1].
     * @details Only the fields that actually change are animated, each in its
     * own equal time slice, so they transition one after another rather than
     * all at once.
     */
    void lerp(const Params &a, const Params &b, float t) {
      constexpr int N = 7;
      // Trips if the field set changes, so the dst/src/tgt arrays and N below
      // can't silently fall out of sync with Params.
      static_assert(sizeof(Params) == N * sizeof(float),
                    "Liquid2D::Params field set changed — update lerp's "
                    "dst/src/tgt arrays and N to match");
      float *dst[N] = {&warp_scale, &warp_strength, &pattern_freq, &time_speed,
                       &complexity,  &pole_fade,     &cycle_speed};
      const float src[N] = {a.warp_scale, a.warp_strength, a.pattern_freq,
                            a.time_speed, a.complexity,    a.pole_fade,
                            a.cycle_speed};
      const float tgt[N] = {b.warp_scale, b.warp_strength, b.pattern_freq,
                            b.time_speed, b.complexity,    b.pole_fade,
                            b.cycle_speed};
      int active = 0;
      for (int i = 0; i < N; ++i)
        if (src[i] != tgt[i])
          ++active;
      if (active == 0)
        return;
      float slice = 1.0f / active;
      int slot = 0;
      for (int i = 0; i < N; ++i) {
        if (src[i] == tgt[i]) {
          *dst[i] = tgt[i];
          continue;
        }
        float tl = hs::clamp((t - slot * slice) / slice, 0.0f, 1.0f);
        *dst[i] = hs::lerp(src[i], tgt[i], tl);
        ++slot;
      }
    }
  };
  Params params; /**< Live per-frame parameters, lerped between presets. */
  float accumulated_time = 0.0f; /**< Unbounded noise-time axis (see draw_frame). */
  float cycle_phase = 0.0f;      /**< Wrapped to [0, 2pi) each frame for breathe. */
  float sin_phase = 0.0f;   /**< Wrapped to [0, 2pi): pattern's +t term. */
  float cos_phase = 0.0f;   /**< Wrapped to [0, 2pi): pattern's 0.8*t term. */

  /**
   * @brief Preset bank cycled by the timeline timer.
   * @details All 7 Params fields are listed explicitly per preset. Omitting the
   * trailing cycle_speed would silently fall back to its default member
   * initializer rather than the preset's intent.
   */
  Presets<Params, 2> presets = {{{
      {{1.5f, 0.5f, 5.0f, 0.1f, 0.5f, 1.4f, 0.05f}},
      {{1.5f, 0.5f, 1.2f, 0.05f, 3.0f, 1.4f, 0.05f}},
  }}};
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Liquid2D)
