/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ShaderBall.h
 * @brief Stereographic shader spanning liquid domain-warp and grid fly-through
 *        looks from one continuous parameter space.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"

// Unit-test accessor for the wrap invariants, lens, and pattern formula.
namespace hs_test {
namespace effects_tests {
struct ShaderBallWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Stereographic shader effect over a two-slot bank of palette cyclers.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Projects a noise-warped sinusoidal pattern through stereographic
 * space. Every look axis — camera spin vs. random-walk wander, glitch-lens
 * blend, pattern cross-coupling vs. direct phase feed, palette slot, breathe
 * depth, hue shift, value fade — is a continuous preset-lerped float, so the
 * preset timeline morphs between any two looks without a discrete pop.
 */
template <int W, int H> class ShaderBall : public Effect {
public:
  // Gamut boundary bracket grid this effect buys from the persistent arena
  // (131,072 B; only the palette bakes share the partition). The flash
  // master's own resolution, so the copy is verbatim and buys only the RAM read
  // latency.
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  /**
   * @brief Constructs the effect at W x H and disables pixel persistence.
   */
  HS_COLD_MEMBER ShaderBall() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Initializes animated params, walks, palette bank, and preset loop.
   * @details Loads the first preset, registers animated params, seeds the
   * random walks and cycle driver, bakes the palette bank and gamut LUT, and
   * enters the preset choreography.
   */
  void init() override {
    warp_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    // Dedicated warp generator, never handed to a walk: RandomWalk's ctor
    // reseeds and re-frequencies whatever generator it is given, which would
    // make the warp field nondeterministic.
    warp_noise.SetFrequency(WARP_NOISE_FREQUENCY);

    // register_param range-checks *ptr and captures it as the GUI default, so
    // params has to hold the first preset before any registration below.
    blend.params = presets.get();

    // Every preset-driven param is flagged animated so "Pause Animation" lets
    // the user take a slider over.
    Params &params = blend.params;
    register_animated_param("Warp Scale", &params.warp_scale, WARP_SCALE_MIN,
                            WARP_SCALE_MAX);
    register_animated_param("Warp Strength", &params.warp_strength,
                            WARP_STRENGTH_MIN, WARP_STRENGTH_MAX);
    register_animated_param("Warp Time", &params.warp_time_scale, WARP_TIME_MIN,
                            WARP_TIME_MAX);
    register_animated_param("Pattern Freq", &params.pattern_freq,
                            PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("Complexity", &params.complexity, COMPLEXITY_MIN,
                            COMPLEXITY_MAX);
    register_animated_param("Phase Direct", &params.phase_direct,
                            PHASE_DIRECT_MIN, PHASE_DIRECT_MAX);
    register_animated_param("Drift", &params.phase2_rate, PHASE2_RATE_MIN,
                            PHASE2_RATE_MAX);
    register_animated_param("Pole Fade", &params.pole_fade, POLE_FADE_MIN,
                            POLE_FADE_MAX);
    register_animated_param("Spin Rate", &params.spin_rate, SPIN_RATE_MIN,
                            SPIN_RATE_MAX);
    register_animated_param("Wander", &params.wander, WANDER_MIN, WANDER_MAX);
    register_animated_param("Lens Mix", &params.lens_mix, LENS_MIX_MIN,
                            LENS_MIX_MAX);
    register_animated_param("Palette", &params.palette_pos, PALETTE_POS_MIN,
                            PALETTE_POS_MAX);
    register_animated_param("Breathe Depth", &params.breathe_depth, BREATHE_MIN,
                            BREATHE_MAX);
    register_animated_param("Cycle Speed", &params.cycle_speed, CYCLE_SPEED_MIN,
                            CYCLE_SPEED_MAX);
    register_animated_param("Hue Shift", &params.hue_shift, HUE_SHIFT_MIN,
                            HUE_SHIFT_MAX);
    register_animated_param("Value Fade", &params.value_fade, VALUE_FADE_MIN,
                            VALUE_FADE_MAX);

    // The walks run permanently on shadow orientations; `wander` mixes their
    // contribution in per frame (update_camera), so a transition into wander
    // picks up a live, well-mixed trajectory rather than a cold start.
    timeline.add(0, Animation::RandomWalk<W>(inner_walk, UP, inner_walk_noise));
    timeline.add(0, Animation::RandomWalk<W>(outer_walk, UP, outer_walk_noise));
    // wrap=false: cycle_phase is wrapped by hand to 2pi in draw_frame; the
    // Driver's [0,1) wrap is the wrong period for fast_sinf.
    timeline.add(0, Animation::Driver(cycle_phase, &blend.params.cycle_speed,
                                      1.0f, false));

    const auto liquid_recipes = EffectPaletteRecipes::shader_ball_liquid_set();
    const auto flyby_recipes = EffectPaletteRecipes::shader_ball_flyby_set();
    for (int i = 0; i < PALETTE_SET_SIZE; ++i) {
      liquid_palettes[i] = GenerativePalette{liquid_recipes[i]};
      flyby_palettes[i] = GenerativePalette{flyby_recipes[i]};
      liquid_entries[i] = liquid_palettes[i];
      flyby_entries[i] = flyby_palettes[i];
    }
    cyclers[0].init(persistent_arena, liquid_entries.data(), PALETTE_SET_SIZE,
                    PALETTE_DWELL_FRAMES, PALETTE_FADE_FRAMES, ease_in_out_sin,
                    &anims_paused);
    cyclers[1].init(persistent_arena, flyby_entries.data(), PALETTE_SET_SIZE,
                    PALETTE_DWELL_FRAMES, PALETTE_FADE_FRAMES, ease_in_out_sin,
                    &anims_paused);
    // The shader's hue_rotate clips per pixel, so the RAM copy pays for itself
    // over reading the flash master.
    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);

    enter_preset();
  }

  /**
   * @brief Advances animation phases and shades every pixel for one frame.
   * @details Steps the timeline, maintains the wrapped accumulators, rebuilds
   * the per-frame camera quaternions, then shades every pixel: outer orient,
   * lens blend, inner orient, stereographic warp, pattern sample, pole
   * normalize, palette-bank lookup with breathe, value fade, hue shift.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sb_timeline_step);
      timeline.step(canvas);
    }

    const Params &P = blend.params;
    // Wrap the noise-time accumulator so the float ULP never swallows the
    // increment and freezes the warp; OpenSimplex2 is aperiodic so the wrap pops
    // the field once per period (far apart at STEREO_NOISE_TIME_PERIOD).
    noise_time = fmodf(noise_time + P.speed, STEREO_NOISE_TIME_PERIOD);
    // Wrap the trig phases to 2pi so fast_sinf/fast_cosf keep precise range
    // reduction.
    sin_phase = fmodf(sin_phase + P.speed, TWO_PI_F);
    phase2 = fmodf(phase2 + P.speed * P.phase2_rate, TWO_PI_F);
    spin_phase = fmodf(spin_phase + P.spin_rate, TWO_PI_F);
    // cycle_phase feeds fast_sinf below, so wrap to 2pi by hand (the Driver's
    // [0,1) wrap is the wrong domain for a radians consumer).
    cycle_phase = fmodf(cycle_phase, TWO_PI_F);

    update_camera(P.wander);

    cyclers[0].step();
    cyclers[1].step();

    // Frame-constant hoists; the per-pixel skip branches below all test these,
    // so the branch predictor makes the skips nearly free.
    const float warp_time = noise_time * P.warp_time_scale;
    const float direct1 = P.phase_direct * sin_phase;
    const float direct2 = P.phase_direct * phase2;
    const float breathe_offset = fast_sinf(cycle_phase) * P.breathe_depth;
    const int pal_lo =
        std::min(static_cast<int>(P.palette_pos), PALETTE_COUNT - 1);
    const float pal_frac = P.palette_pos - static_cast<float>(pal_lo);

    auto shader = [&](const Vector &v) -> Color4 {
      Vector rv = rotate(v, cam_outer_conj);
      Vector sv = rv;
      if (P.lens_mix != 0.0f) {
        Vector g = apply_glitch_lens(rv);
        sv = P.lens_mix == 1.0f ? g : nlerp_dir(rv, g, P.lens_mix);
      }
      Complex z = stereo(rotate(sv, cam_inner_conj));
      float r_sq = z.re * z.re + z.im * z.im;
      auto [w, displacement] =
          stereo_noise_warp(z, r_sq, warp_noise, P.warp_scale, P.warp_strength,
                            P.pole_fade, warp_time);
      float pattern =
          sample_pattern(stereo_pattern_args(w, P.pattern_freq), P.complexity,
                         direct1, direct2, sin_phase, phase2);
      float value = pole_normalize_pattern(pattern, r_sq, P.pole_fade);
      // The Wrap lookup below folds an exact 1.0 onto the palette's other end.
      value = std::min(value, ONE_BELOW_UNIT);
      float u = wrap_t(value + breathe_offset);
      Color4 c = pal_frac == 0.0f
                     ? cyclers[pal_lo].palette().get(u)
                     : cyclers[pal_lo].palette().get(u).lerp(
                           cyclers[pal_lo + 1].palette().get(u), pal_frac);
      c.alpha *= (1.0f - value * P.value_fade);
      if (P.hue_shift != 0.0f)
        c = hue_rotate(c, -displacement * P.hue_shift);
      return c;
    };

    {
      HS_PROFILE(sb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
  }

private:
  // Test seam: reaches the accumulators, lens, and pattern formula.
  friend struct ::hs_test::effects_tests::ShaderBallWhiteBox;

  /**
   * @brief Evaluates the generalized sinusoidal pattern at a bounded point.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param complexity Cross-coupling depth c.
   * @param direct1 Pre-multiplied direct phase feed s * sin_phase.
   * @param direct2 Pre-multiplied direct phase feed s * phase2.
   * @param sin_phase Wrapped +t phase in [0, 2pi).
   * @param phase2 Wrapped drift phase in [0, 2pi).
   * @return sin(p.re + c*sin(p.im + phase1) + s*phase1)
   *         * cos(p.im + c*cos(p.re - phase2) - s*phase2) in [-1, 1].
   * @details Exact superset of both classic forms: c != 0, s = 0 is the
   * cross-coupled liquid pattern; c = 0, s = 1 is the separable grid pattern.
   */
  static float sample_pattern(const Complex &p, float complexity, float direct1,
                              float direct2, float sin_phase, float phase2) {
    float re = p.re + direct1;
    float im = p.im - direct2;
    if (complexity != 0.0f) {
      re += complexity * fast_sinf(p.im + sin_phase);
      im += complexity * fast_cosf(p.re - phase2);
    }
    return fast_sinf(re) * fast_cosf(im);
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

  /**
   * @brief Normalized linear blend between two unit directions.
   * @param a Blend source (t = 0), unit length.
   * @param b Blend target (t = 1), unit length.
   * @param t Blend factor in (0, 1).
   * @return Unit direction along the chord from a to b.
   * @details Guards on the squared magnitude of the blend vector (never an
   * acos-derived angle): near-antipodal endpoints collapse it to zero, in
   * which case the target direction is returned.
   */
  static Vector nlerp_dir(const Vector &a, const Vector &b, float t) {
    Vector m = a + (b - a) * t;
    float m_sq = dot(m, m);
    if (m_sq < LENS_BLEND_MIN_SQ)
      return b;
    return m / sqrtf(m_sq);
  }

  /**
   * @brief Rebuilds the per-frame camera quaternions from walks, spin, wander.
   * @param wander Random-walk motion gain in [0, 1]; 0 freezes the walk camera,
   * 1 applies its full per-frame motion.
   * @details Integrates a fraction of each walk's incremental rotation. This
   * keeps partial wander continuous through every full turn and while the
   * wander amount changes.
   */
  // Per-frame, not per-pixel: the inlined slerps are too big for ITCM.
  HS_COLD_MEMBER void update_camera(float wander) {
    const Quaternion inner = inner_walk.get();
    const Quaternion inner_delta = inner * inner_walk_prev.conjugate();
    inner_wander =
        (slerp(Quaternion(), inner_delta.normalized(), wander) * inner_wander)
            .normalized();
    inner_walk_prev = inner;

    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_prev.conjugate();
    outer_wander =
        (slerp(Quaternion(), outer_delta.normalized(), wander) * outer_wander)
            .normalized();
    outer_walk_prev = outer;

    cam_inner_conj =
        (make_rotation(Y_AXIS, spin_phase) * base_orientation * inner_wander)
            .conjugate();
    cam_outer_conj = outer_wander.conjugate();
  }

  /**
   * @brief Runs the choreography step for the preset just entered.
   * @details Consumes the entered preset's Choreo entry: with no dwell the
   * blend into the next preset chains immediately (a repeating RandomTimer(0,
   * 0) would spawn a fresh Lerp every frame); otherwise a one-shot RandomTimer
   * holds for dwell_min..dwell_max frames first.
   */
  HS_COLD_MEMBER void enter_preset() {
    const Choreo &ch = CHOREO[presets.current_index()];
    if (ch.dwell_max == 0) {
      begin_blend();
    } else {
      timeline.add_pausable(
          0,
          Animation::RandomTimer(ch.dwell_min, ch.dwell_max,
                                 [this](Canvas &) { begin_blend(); }),
          &anims_paused);
    }
  }

  /**
   * @brief Schedules the blend from the current preset into the next.
   * @details Advances the preset selector and lerps params into it over the
   * outgoing preset's blend_frames, staggered or parallel per its Choreo
   * entry; the blend's completion re-enters the state machine.
   */
  HS_COLD_MEMBER void begin_blend() {
    const Choreo &ch = CHOREO[presets.current_index()];
    presets.next();
    blend.staggered = ch.staggered;
    blend_from.params = presets.prev_get();
    blend_to.params = presets.get();
    timeline.add_pausable(0,
                          Animation::Lerp(blend, blend_from, blend_to,
                                          ch.blend_frames, ease_in_out_sin)
                              .then([this]() { enter_preset(); }),
                          &anims_paused);
  }

  /**
   * @brief Inner per-pixel walk orientation.
   * @details Declared before `timeline` so it outlives the RandomWalk that
   * points here, which ~Timeline clears on teardown.
   */
  Orientation<> inner_walk;
  /**
   * @brief Outer whole-sphere walk orientation.
   * @details Declared before `timeline` so it outlives the RandomWalk that
   * points here, which ~Timeline clears on teardown.
   */
  Orientation<> outer_walk;
  Timeline timeline; /**< Drives walks, cycle driver, and preset blends. */
  /** @brief OpenSimplex2 source for the warp field; deterministic (default
   *  seed, fixed frequency) and never handed to a RandomWalk. */
  FastNoiseLite warp_noise;
  FastNoiseLite inner_walk_noise; /**< OpenSimplex2 source, `inner_walk`. */
  /** @brief OpenSimplex2 source for the `outer_walk`; a generator shared with
   *  the inner walk would drive both along the same trajectory. */
  FastNoiseLite outer_walk_noise;

  /** @brief Fixed camera pre-rotation under the Y spin. */
  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  Quaternion inner_walk_prev; /**< Previous inner walk orientation. */
  Quaternion outer_walk_prev; /**< Previous outer walk orientation. */
  Quaternion inner_wander;    /**< Integrated inner walk camera. */
  Quaternion outer_wander;    /**< Integrated outer walk camera. */
  Quaternion cam_inner_conj;  /**< Per-frame inverse of the inner camera. */
  Quaternion cam_outer_conj;  /**< Per-frame inverse of the outer camera. */

  float noise_time = 0.0f;  /**< Noise-time axis, wrapped to
                               STEREO_NOISE_TIME_PERIOD. */
  float sin_phase = 0.0f;   /**< Wrapped to [0, 2pi): the pattern's +t term. */
  float phase2 = 0.0f;      /**< Wrapped to [0, 2pi): pattern's drift term. */
  float spin_phase = 0.0f;  /**< Wrapped to [0, 2pi): Y-spin angle. */
  float cycle_phase = 0.0f; /**< Wrapped to [0, 2pi) each frame for breathe. */

  static constexpr int PALETTE_COUNT = 2;
  static constexpr int PALETTE_SET_SIZE = 3;
  /** @brief Frames each cycler dwells on a palette between fades. */
  static constexpr int PALETTE_DWELL_FRAMES = 900;
  /** @brief Frames each palette fade spans. */
  static constexpr int PALETTE_FADE_FRAMES = 300;
  std::array<GenerativePalette, PALETTE_SET_SIZE>
      liquid_palettes; /**< Cycled sources for slot 0. */
  std::array<GenerativePalette, PALETTE_SET_SIZE>
      flyby_palettes; /**< Cycled sources for slot 1. */
  std::array<PaletteCycler::Entry, PALETTE_SET_SIZE> liquid_entries;
  std::array<PaletteCycler::Entry, PALETTE_SET_SIZE> flyby_entries;
  std::array<PaletteCycler, PALETTE_COUNT>
      cyclers; /**< Per-slot palette cycles: liquid set, flyby set. */

  /** @brief Base spatial frequency of the warp generator; the `Warp Scale`
   *  slider multiplies it. */
  static constexpr float WARP_NOISE_FREQUENCY = 0.01f;
  /** @brief Largest float below 1.0; keeps palette coordinates off the Wrap
   *  fold. */
  static constexpr float ONE_BELOW_UNIT = 0x1.fffffep-1f;
  /** @brief Squared-magnitude floor for the lens nlerp blend vector. */
  static constexpr float LENS_BLEND_MIN_SQ = 1e-12f;

  /**
   * @brief Tunable warp/pattern/camera/color state, one snapshot per preset.
   */
  struct Params {
    float warp_scale;
    float warp_strength;
    float warp_time_scale;
    float pattern_freq;
    float speed;
    float complexity;
    float phase_direct;
    float phase2_rate;
    float pole_fade;
    float spin_rate;
    float wander;
    float lens_mix;
    float palette_pos;
    float breathe_depth;
    float cycle_speed;
    float hue_shift;
    float value_fade;

    static constexpr int N = 17;

    /**
     * @brief Interpolates every field in parallel from a to b.
     * @param a Source params (t = 0).
     * @param b Destination params (t = 1).
     * @param t Blend factor in [0, 1].
     */
    HS_COLD_MEMBER void lerp(const Params &a, const Params &b, float t) {
      for (auto f : FIELDS)
        this->*f = hs::lerp(a.*f, b.*f, t);
    }

    /**
     * @brief Staggered interpolation of this Params from a to b.
     * @param a Source parameter set (value at t = 0).
     * @param b Target parameter set (value at t = 1).
     * @param t Normalized progress in [0, 1].
     * @details Only the fields that actually change are animated, each in its
     * own equal time slice, so they transition one after another rather than
     * all at once.
     */
    HS_COLD_MEMBER void lerp_staggered(const Params &a, const Params &b,
                                       float t) {
      int active = 0;
      for (auto f : FIELDS)
        if (a.*f != b.*f)
          ++active;
      if (active == 0)
        return;
      int slot = 0;
      for (auto f : FIELDS) {
        if (a.*f == b.*f) {
          this->*f = b.*f;
          continue;
        }
        // t * active - slot, not (t - slot/active)/(1/active): the product is
        // exact at t = 1, so every slot's slice reaches exactly 1 and gate
        // params land bit-exactly on their endpoints.
        float tl = hs::clamp(t * active - slot, 0.0f, 1.0f);
        this->*f = hs::lerp(a.*f, b.*f, tl);
        ++slot;
      }
    }

    /** @brief Field table both lerps walk, so neither can drop a member. */
    static constexpr std::array<float Params::*, N> FIELDS = {
        &Params::warp_scale,   &Params::warp_strength, &Params::warp_time_scale,
        &Params::pattern_freq, &Params::speed,         &Params::complexity,
        &Params::phase_direct, &Params::phase2_rate,   &Params::pole_fade,
        &Params::spin_rate,    &Params::wander,        &Params::lens_mix,
        &Params::palette_pos,  &Params::breathe_depth, &Params::cycle_speed,
        &Params::hue_shift,    &Params::value_fade};
  };
  // Trips if the field set changes, so FIELDS and the bare-float preset rows
  // can't silently fall out of sync with Params.
  static_assert(sizeof(Params) == Params::N * sizeof(float),
                "ShaderBall::Params field set changed — update FIELDS and the "
                "preset float lists to match");

  /**
   * @brief Live params plus the blend-style flag consumed by Animation::Lerp.
   */
  struct Blend {
    Params params;
    bool staggered = false;

    /**
     * @brief Forwards to the staggered or parallel Params lerp.
     * @param a Source blend state (t = 0).
     * @param b Destination blend state (t = 1).
     * @param t Blend factor in [0, 1].
     */
    HS_COLD_MEMBER void lerp(const Blend &a, const Blend &b, float t) {
      if (staggered)
        params.lerp_staggered(a.params, b.params, t);
      else
        params.lerp(a.params, b.params, t);
    }
  };

  /**
   * @brief Choreography entry consumed when entering the matching preset:
   * sequential dwell (hold), then a blend into the next preset.
   */
  struct Choreo {
    uint16_t dwell_min;    /**< Minimum hold before blending, frames. */
    uint16_t dwell_max;    /**< Maximum hold; 0 chains the blend directly. */
    uint16_t blend_frames; /**< Blend duration into the next preset. */
    bool staggered;        /**< Staggered vs. parallel field interpolation. */
  };

  static constexpr float WARP_SCALE_MIN = 0.1f, WARP_SCALE_MAX = 100.0f;
  static constexpr float WARP_STRENGTH_MIN = 0.0f, WARP_STRENGTH_MAX = 30.0f;
  static constexpr float WARP_TIME_MIN = 0.05f, WARP_TIME_MAX = 1.0f;
  static constexpr float PATTERN_FREQ_MIN = 1.0f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float COMPLEXITY_MIN = 0.0f, COMPLEXITY_MAX = 3.0f;
  static constexpr float PHASE_DIRECT_MIN = 0.0f, PHASE_DIRECT_MAX = 1.0f;
  static constexpr float PHASE2_RATE_MIN = 0.0f, PHASE2_RATE_MAX = 2.0f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float SPIN_RATE_MIN = 0.0f, SPIN_RATE_MAX = 0.05f;
  static constexpr float WANDER_MIN = 0.0f, WANDER_MAX = 1.0f;
  static constexpr float LENS_MIX_MIN = 0.0f, LENS_MIX_MAX = 1.0f;
  static constexpr float PALETTE_POS_MIN = 0.0f;
  static constexpr float PALETTE_POS_MAX =
      static_cast<float>(PALETTE_COUNT - 1);
  static constexpr float BREATHE_MIN = 0.0f, BREATHE_MAX = 0.3f;
  static constexpr float CYCLE_SPEED_MIN = 0.0f, CYCLE_SPEED_MAX = 1.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 1.0f;
  static constexpr float VALUE_FADE_MIN = 0.0f, VALUE_FADE_MAX = 1.0f;

  /** @brief Y-spin rate reproducing the classic 300-frame grid orbit. */
  static constexpr float ORBIT_SPIN_RATE = TWO_PI_F / 300.0f;
  /** @brief High-contrast, liquid-biased palette blend used by every preset. */
  static constexpr float PRESET_PALETTE_POS = 0.02f;

  /** @brief True iff every preset-driven field of @p p lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Params &p) {
    return p.warp_scale >= WARP_SCALE_MIN && p.warp_scale <= WARP_SCALE_MAX &&
           p.warp_strength >= WARP_STRENGTH_MIN &&
           p.warp_strength <= WARP_STRENGTH_MAX &&
           p.warp_time_scale >= WARP_TIME_MIN &&
           p.warp_time_scale <= WARP_TIME_MAX &&
           p.pattern_freq >= PATTERN_FREQ_MIN &&
           p.pattern_freq <= PATTERN_FREQ_MAX && p.speed >= SPEED_MIN &&
           p.speed <= SPEED_MAX && p.complexity >= COMPLEXITY_MIN &&
           p.complexity <= COMPLEXITY_MAX &&
           p.phase_direct >= PHASE_DIRECT_MIN &&
           p.phase_direct <= PHASE_DIRECT_MAX &&
           p.phase2_rate >= PHASE2_RATE_MIN &&
           p.phase2_rate <= PHASE2_RATE_MAX && p.pole_fade >= POLE_FADE_MIN &&
           p.pole_fade <= POLE_FADE_MAX && p.spin_rate >= SPIN_RATE_MIN &&
           p.spin_rate <= SPIN_RATE_MAX && p.wander >= WANDER_MIN &&
           p.wander <= WANDER_MAX && p.lens_mix >= LENS_MIX_MIN &&
           p.lens_mix <= LENS_MIX_MAX && p.palette_pos >= PALETTE_POS_MIN &&
           p.palette_pos <= PALETTE_POS_MAX && p.breathe_depth >= BREATHE_MIN &&
           p.breathe_depth <= BREATHE_MAX && p.cycle_speed >= CYCLE_SPEED_MIN &&
           p.cycle_speed <= CYCLE_SPEED_MAX && p.hue_shift >= HUE_SHIFT_MIN &&
           p.hue_shift <= HUE_SHIFT_MAX && p.value_fade >= VALUE_FADE_MIN &&
           p.value_fade <= VALUE_FADE_MAX;
  }

  // Field order per row: warp_scale, warp_strength, warp_time_scale,
  // pattern_freq, speed, complexity, phase_direct, phase2_rate, pole_fade,
  // spin_rate, wander, lens_mix, palette_pos, breathe_depth, cycle_speed,
  // hue_shift, value_fade. Gate params (complexity, phase_direct, wander,
  // lens_mix, hue_shift, value_fade) must sit exactly on 0 or 1 whenever the
  // matching per-pixel skip is intended: Animation::Lerp lands bit-exactly
  // only on those endpoints, and a near-0 value un-latches the skip for good.
  static constexpr std::array<PresetEntry<Params>, 12> PRESETS = {{
      // Wandering liquid: mild, deep, then fine-grained cross-coupling.
      {{3.0f, 0.5f, 0.5f, 5.0f, 0.1f, 0.5f, 0.0f, 0.8f, 1.4f, 0.0f, 1.0f, 1.0f,
        PRESET_PALETTE_POS, 0.15f, 0.05f, 0.0f, 0.0f}},
      {{3.0f, 0.5f, 0.5f, 1.2f, 0.05f, 3.0f, 0.0f, 0.8f, 1.4f, 0.0f, 1.0f, 1.0f,
        PRESET_PALETTE_POS, 0.15f, 0.05f, 0.0f, 0.0f}},
      {{3.0f, 1.479f, 0.5f, 14.528f, 0.1f, 0.5f, 0.0f, 0.8f, 1.0f, 0.0f, 1.0f,
        1.0f, PRESET_PALETTE_POS, 0.15f, 0.05f, 0.0f, 0.0f}},
      {{3.0f, 0.0f, 0.5f, 15.763f, 0.1f, 2.950552f, 0.0f, 0.8f, 1.0f, 0.0f,
        1.0f, 0.0f, PRESET_PALETTE_POS, 0.15f, 0.05f, 0.0f, 0.0f}},
      {{0.1f, 13.47f, 0.5f, 3.28f, 0.1f, 2.463f, 0.0f, 0.8f, 1.209f, 0.03725f,
        0.252f, 0.066f, 0.022f, 0.19710001f, 0.02f, 0.011f, 0.0f}},
      // Spinning grid fly-throughs.
      {{47.752f, 11.55f, 0.3f, 2.7f, 0.586f, 0.0f, 1.0f, 0.7f, 1.55f,
        ORBIT_SPIN_RATE, 0.0f, 0.0f, PRESET_PALETTE_POS, 0.0f, 0.0f, 0.097f,
        1.0f}},
      {{0.1f, 0.87f, 0.3f, 14.262f, 0.586f, 0.0f, 1.0f, 0.7f, 3.527f,
        ORBIT_SPIN_RATE, 0.0f, 0.0f, PRESET_PALETTE_POS, 0.0f, 0.0f, 0.097f,
        1.0f}},
      {{1.5f, 0.5f, 0.3f, 8.0f, 0.30f, 0.0f, 1.0f, 0.7f, 2.0f, ORBIT_SPIN_RATE,
        0.0f, 0.0f, PRESET_PALETTE_POS, 0.0f, 0.0f, 0.15f, 1.0f}},
      {{47.752f, 2.55f, 0.3f, 7.878f, 0.562f, 0.0f, 1.0f, 0.7f, 2.843f,
        ORBIT_SPIN_RATE, 0.0f, 0.0f, PRESET_PALETTE_POS, 0.0f, 0.0f, 0.0f,
        1.0f}},
      {{100.0f, 8.67f, 0.3f, 1.0f, 0.586f, 0.0f, 1.0f, 0.7f, 3.432f,
        ORBIT_SPIN_RATE, 0.0f, 0.0f, PRESET_PALETTE_POS, 0.0f, 0.0f, 0.636f,
        1.0f}},
      // Grid look, liquid palette.
      {{50.749298f, 30.0f, 0.4699f, 1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f,
        1.5482996f, 0.020879198f, 0.0030917525f, 0.0f, PRESET_PALETTE_POS,
        0.25410002f, 0.00015458837f, 0.201f, 0.847f}},
      {{38.761299f, 30.0f, 0.4699f, 1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f,
        1.5482996f, 0.020879198f, 0.0030917525f, 0.0f, PRESET_PALETTE_POS,
        0.25410002f, 0.00015458837f, 0.201f, 0.847f}},
  }};
  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "a ShaderBall preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset (the "
                "range exposes the presets, it does not clamp them)");

  /** @brief Per-preset choreography, consumed on entry; the cross-family
   *  boundaries (rows 3 and 11) blend parallel. */
  static constexpr std::array<Choreo, 12> CHOREO = {{
      {30, 90, 60, true},
      {30, 90, 60, true},
      {30, 90, 60, true},
      {30, 90, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
  }};
  static_assert(CHOREO.size() == PRESETS.size(),
                "CHOREO must carry one entry per preset");

  Presets<Params, 12> presets{PRESETS};

  Blend blend{PRESETS[0].params}; /**< Live params; init() reloads from
                                     presets. */
  Blend blend_from; /**< Blend start snapshot owned across the Lerp. */
  Blend blend_to;   /**< Blend target snapshot owned across the Lerp. */

  // init() buys the gamut bracket grid and the palette-bank LUTs from the
  // persistent arena. Effect keeps the default arena split, so the total must
  // fit the device persistent partition.
  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      PALETTE_COUNT * PaletteCycler::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "ShaderBall persistent footprint exceeds the default "
                "partition; coarsen the gamut grid or carve arenas");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShaderBall)
