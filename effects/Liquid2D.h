/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

// Unit-test accessor for the noise-time and trig-phase wrap invariants.
namespace hs_test {
namespace effects_tests {
struct Liquid2DWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Full-sphere liquid effect over a breathing generative palette.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Domain-warped OpenSimplex noise feeds a cross-coupled sinusoidal
 * pattern, stereographically projected and colored through a breathing
 * generative palette. Presets cycle on a random timer.
 * @note `Flyby` is the sibling stereographic effect, sharing the core primitives
 *       but with its own `project()` (dual orientation + glitch lens), a
 *       cross-coupled `sample()`, and a staggered `Params::lerp`; the two are
 *       independent, so changes need not propagate. Common per-pixel code lives
 *       in core.
 */
template <int W, int H> class Liquid2D : public Effect {
public:
  /**
   * @brief Constructs the effect and disables inter-frame pixel persistence.
   * @details The shader writes every pixel each frame, so no persistence is
   * needed.
   */
  HS_COLD_MEMBER Liquid2D() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Initializes sliders, timeline drivers, palette, and preset cycling.
   * @details Registers sliders, seeds the timeline drivers/orientations, bakes
   * the palette, and arms the preset-cycling timer; then loads the first preset.
   */
  void init() override {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    register_animated_param("Warp Scale", &params.warp_scale, WARP_SCALE_MIN, WARP_SCALE_MAX);
    register_animated_param("Warp Strength", &params.warp_strength, WARP_STRENGTH_MIN, WARP_STRENGTH_MAX);
    register_animated_param("Pattern Freq", &params.pattern_freq, PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_animated_param("Time Speed", &params.time_speed, TIME_SPEED_MIN, TIME_SPEED_MAX);
    register_animated_param("Complexity", &params.complexity, COMPLEXITY_MIN, COMPLEXITY_MAX);
    register_animated_param("Pole Fade", &params.pole_fade, POLE_FADE_MIN, POLE_FADE_MAX);
    register_animated_param("Cycle Speed", &params.cycle_speed, CYCLE_SPEED_MIN, CYCLE_SPEED_MAX);

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
                          if (animations_paused())
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

  /**
   * @brief Advances the timeline and shades every pixel for one frame.
   * @details Advances the timeline, maintains the wrapped trig phases, then
   * shades every pixel: projects to the sphere, warps in stereographic space,
   * samples the pattern, attenuates near the poles, and maps through the
   * palette.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(lq_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(lq_timeline_step);
      timeline.step(canvas);
    }

    // Wrap the noise-time accumulator so the float ULP never swallows the
    // increment and freezes the warp; OpenSimplex2 is aperiodic so the wrap pops
    // the field once per period (far apart at TIME_PERIOD).
    accumulated_time = fmodf(accumulated_time, TIME_PERIOD);
    float t = accumulated_time;
    // dt is the live Time Speed slider, not (t - prev): the Driver adds that speed
    // each step so it IS this frame's advance, and differencing across the wrap
    // seam would spike.
    constexpr float TWO_PI_F = 2.0f * PI_F;
    float dt = params.time_speed;
    sin_phase = fmodf(sin_phase + dt, TWO_PI_F);
    cos_phase = fmodf(cos_phase + 0.8f * dt, TWO_PI_F);
    // cycle_phase feeds BreatheModifier's fast_sinf, so wrap to 2pi by hand (the
    // Driver's [0,1) wrap is the wrong domain for a radians consumer).
    cycle_phase = fmodf(cycle_phase, TWO_PI_F);

    auto shader = [&](const Vector &v) -> Color4 {
      Complex z = project(v);
      float r_sq = z.re * z.re + z.im * z.im;
      Complex w = warp(z, r_sq, t).coords; // only the warped coords feed sample
      float pattern = sample(w, sin_phase, cos_phase);
      float value = pole_normalize_pattern(pattern, r_sq, params.pole_fade);
      return static_palette.get(value);
    };

    {
      HS_PROFILE(lq_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
  }

private:
  // Test seam: reaches the noise-time and trig-phase wrap invariants.
  friend struct ::hs_test::effects_tests::Liquid2DWhiteBox;

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
  static constexpr float TIME_PERIOD = 65536.0f;
  Params params; /**< Live per-frame parameters, lerped between presets. */
  float accumulated_time = 0.0f; /**< Noise-time axis, wrapped to TIME_PERIOD (see draw_frame). */
  float cycle_phase = 0.0f;      /**< Wrapped to [0, 2pi) each frame for breathe. */
  float sin_phase = 0.0f;   /**< Wrapped to [0, 2pi): pattern's +t term. */
  float cos_phase = 0.0f;   /**< Wrapped to [0, 2pi): pattern's 0.8*t term. */

  static constexpr float WARP_SCALE_MIN = 0.1f, WARP_SCALE_MAX = 10.0f;
  static constexpr float WARP_STRENGTH_MIN = 0.0f, WARP_STRENGTH_MAX = 3.0f;
  static constexpr float PATTERN_FREQ_MIN = 1.0f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float TIME_SPEED_MIN = 0.05f, TIME_SPEED_MAX = 5.0f;
  static constexpr float COMPLEXITY_MIN = 0.5f, COMPLEXITY_MAX = 3.0f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float CYCLE_SPEED_MIN = 0.0f, CYCLE_SPEED_MAX = 1.0f;

  /** @brief True iff every preset-driven field of @p p lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Params &p) {
    return p.warp_scale >= WARP_SCALE_MIN && p.warp_scale <= WARP_SCALE_MAX &&
           p.warp_strength >= WARP_STRENGTH_MIN &&
           p.warp_strength <= WARP_STRENGTH_MAX &&
           p.pattern_freq >= PATTERN_FREQ_MIN &&
           p.pattern_freq <= PATTERN_FREQ_MAX &&
           p.time_speed >= TIME_SPEED_MIN && p.time_speed <= TIME_SPEED_MAX &&
           p.complexity >= COMPLEXITY_MIN && p.complexity <= COMPLEXITY_MAX &&
           p.pole_fade >= POLE_FADE_MIN && p.pole_fade <= POLE_FADE_MAX &&
           p.cycle_speed >= CYCLE_SPEED_MIN && p.cycle_speed <= CYCLE_SPEED_MAX;
  }

  /**
   * @brief Preset bank cycled by the timeline timer.
   * @details All 7 Params fields are listed explicitly; a trailing omission would
   * silently fall back to the default member initializer, not the preset's value.
   */
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS = {{
      {{1.5f, 0.5f, 5.0f, 0.1f, 0.5f, 1.4f, 0.05f}},
      {{1.5f, 0.5f, 1.2f, 0.05f, 3.0f, 1.4f, 0.05f}},
  }};
  static_assert(preset_in_ranges(PRESETS[0].params) &&
                    preset_in_ranges(PRESETS[1].params),
                "a Liquid2D preset drives a param outside its registered slider "
                "range; widen the range to accommodate the preset (the range "
                "exposes the presets, it does not clamp them)");

  Presets<Params, 2> presets{PRESETS};
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Liquid2D)
