/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file DisplacementField.h
 * @brief Stack of soft-stroked rings displaced by a noise and falling-ball
 *        displacement-field stack.
 */

#include "core/platform/effect_identity.h"
#include <numeric>
#include <bit>

#include "core/animation/orientation.h"
#include "effects/common/palette_recipes.h"
#include "core/engine/engine.h"

namespace hs_test {
namespace effects_tests {
struct DisplacementFieldWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Stack of evenly spaced soft-stroked rings displaced by a
 * displacement-field stack that alternates a 3D noise field and falling ball
 * bumps.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each ring vertex is displaced along the stack axis by the dominant
 * ball field plus the noise field. A noise phase alternates with a ball phase
 * in which cap-shaped bumps fall from world +Y to -Y on random meridians.
 */
template <int W, int H> class DisplacementField : public Effect {
  friend struct ::hs_test::effects_tests::DisplacementFieldWhiteBox;

public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::DisplacementField::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /** Stable preset identities, in automatic sequence order. */
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"noise", "balls"};
  /** @return Safe full-cycle frame budget at shipping defaults, including ball drain. */
  static consteval size_t cycle_frames() {
    constexpr auto DEFAULTS = PRESETS[BALLS_PRESET];
    return 2 * NOISE_FADE_FRAMES + NOISE_HOLD_FRAMES + BALL_PHASE_FRAMES +
           static_cast<int>(BALL_RATE_FPS / DEFAULTS.ball_speed_min) +
           PHASE_HANDOFF_FRAMES;
  }

  /** @return Number of authored displacement presets. */
  static consteval size_t authored_preset_count() { return PRESET_IDS.size(); }

#if HS_ENABLE_EFFECT_CONTROL_API
  /** @brief Selects one authored preset for profiling.
   * @param index Authored preset index. */
  void profile_select_preset(size_t index) {
    HS_CHECK(index < authored_preset_count(),
             "profile preset index out of range");
    HS_CHECK(selectPreset(index), "profile preset selection failed");
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(authored_preset_count()));
  }
#endif

  /**
   * @brief Builds the effect with its palette.
   */
  HS_COLD_MEMBER DisplacementField()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        balls(timeline), noise_field(timeline), palette(make_palette()) {
    // The authored downstream RNG stream starts one draw past the palette.
    static_cast<void>(hs::rand_int(0, 256));
  }

  /**
   * @brief Allocates the bake LUTs, registers params, seeds the noise field,
   * and builds the timeline.
   */
  void init() override {
    configure_presets(PRESET_IDS.size());
    hue_table.knots = persistent_arena.allocate_n<Pixel>(HueTable::KNOTS);
    ball_colat = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_reach = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_scale = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_params =
        persistent_arena.allocate_n<const Animation::BumpParams *>(MAX_BALLS);
    ball_local = persistent_arena.allocate_n<int>(MAX_BALLS);
    ball_cv = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_rho = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_azimuth = persistent_arena.allocate_n<float>(MAX_BALLS);
    rings.init_storage(persistent_arena);
    candidates = persistent_arena.allocate_n<CandidateTable>(1);
    octave1 = persistent_arena.allocate_n<float>(W + 1);
    octave2 = persistent_arena.allocate_n<float>(W + 1);
    alternate_num = persistent_arena.allocate_n<float>(W + 1);
    alternate_den = persistent_arena.allocate_n<float>(W + 1);
    knot_visible = persistent_arena.allocate_n<uint8_t>(W + 1);
    knot_pos = persistent_arena.allocate_n<math::Vector>(W + 1);
    knot_near = persistent_arena.allocate_n<uint8_t>(W + 1);
    drawn_colors = persistent_arena.allocate_n<Color4>(RING_SLOTS);
    drawn_color_ready = persistent_arena.make_n<uint8_t>(RING_SLOTS);
    chunk_cos = persistent_arena.allocate_n<float>(BAKE_CHUNKS);
    chunk_sin = persistent_arena.allocate_n<float>(BAKE_CHUNKS);
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const float a = (2.0f * c + 1.0f) * (math::PI_F / BAKE_CHUNKS);
      chunk_cos[c] = cosf(a);
      chunk_sin[c] = sinf(a);
    }
    balls.init_storage(persistent_arena);
    noise_field.init_storage(persistent_arena);

    noise_field.template_params.noise.SetSeed(hs::rand_int(0, 65536));

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_int_param("Rings", &params.num_rings, 1, RING_SLOTS);
    register_param("Thickness", &params.thickness, 0.4f * THICKNESS_PX,
                   6.0f * THICKNESS_PX);
    register_param("Ball Amp", &params.ball_amp, 0.0f, 0.8f);
    register_param("Noise Amp", &params.noise_amp, 0.0f, 0.8f);
    register_param("Scale 1", &params.scale1, 0.5f, 4.0f);
    register_param("Scale 2", &params.scale2, 0.5f, 8.0f);
    register_param("Hue Rotate", &params.hue_scale, 0.0f, 3.0f);
    register_param("Flow Speed", &params.flow_speed, 0.0f, 0.15f);
    register_param("Ball Min", &params.ball_min, 0.05f, 1.0f);
    register_param("Ball Max", &params.ball_max, 0.05f, 1.0f);
    register_param("Ball Rate", &params.ball_rate, 0.5f, 50.0f);
    register_param("Speed Min", &params.ball_speed_min, 0.1f, 3.0f);
    register_param("Speed Max", &params.ball_speed_max, 0.1f, 3.0f);

    // Pinned first, before any finite event can precede it in the buffer.
    noise_field.spawn_pinned(0);

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->draw_rings(canvas, opacity);
                        },
                        -1, {.fade_in = {24}}));

    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, STACK_AXIS, walk_noise));

    timeline.add(
        0, Animation::PeriodicTimer(PALETTE_CYCLE_FRAMES,
                                    [this](Canvas &) { this->roll_palette(); },
                                    {.repeat = true}));

    enter_noise();
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, advances the
   * NOISE -> BALLS phase machine, then renders one frame.
   * @details Preset pause holds phase departures; entry fades and ambient motion continue.
   */
  void draw_frame() override {
    color_spin = math::wrap_t(color_spin + COLOR_SPIN_RATE);

    if (!seed_balls_pending) {
      switch (preset_index) {
      case BALLS_PRESET:
        // Slider-paced spawner: the cooldown re-reads Ball Rate on every spawn.
        if (ball_phase_left > 0 || anims_paused) {
          if (!anims_paused)
            --ball_phase_left;
          spawn_cooldown -= 1.0f;
          while (spawn_cooldown <= 0.0f) {
            spawn_ball();
            spawn_cooldown +=
                hs::rand_f(0.5f, 1.5f) * BALL_RATE_FPS / params.ball_rate;
          }
        } else if (!anims_paused && balls.active_count() == 0) {
          HS_CHECK(advance_preset(), "DisplacementField: noise advance failed");
        }
        break;
      case NOISE_PRESET:
        if (anims_paused)
          break;
        if (master_gain >= 1.0f && noise_hold > 0 && --noise_hold == 0)
          schedule_gain(0.0f, true);
        if (noise_hold == 0 && master_gain <= 0.0f)
          HS_CHECK(advance_preset(), "DisplacementField: balls advance failed");
        break;
      }
    }

    balls.template_params.amplitude =
        params.ball_amp * BALL_DRAPE_PER_AMPLITUDE;
    noise_field.template_params.amplitude =
        preset_index == NOISE_PRESET ? params.noise_amp * master_gain : 0.0f;
    noise_field.template_params.scale1 = params.scale1;
    noise_field.template_params.scale2 = params.scale2;
    noise_field.template_params.speed = params.flow_speed;
    {
      HS_PROFILE(df_prepare_fields);
      balls.prepare_frame();
      noise_field.prepare_frame();
    }

    Canvas canvas(*this);
    {
      HS_PROFILE(df_timeline_step);
      timeline.step(canvas);
    }
    if (seed_balls_pending) {
      seed_balls_pending = false;
      for (int i = 0; i < PREVIEW_BALLS; ++i)
        spawn_ball((i + 1.0f) / (PREVIEW_BALLS + 1.0f), &canvas);
    }
  }

#ifdef HS_PROFILE_ENABLE
  /** @brief Logs automatic preset transitions during profiling.
   * @param change Committed preset transition. */
  HS_COLD_MEMBER void preset_changed(const PresetChange &change) override {
    hs::log("Preset: %u/%u", static_cast<unsigned>(change.to + 1),
            static_cast<unsigned>(authored_preset_count()));
  }
#endif

private:
  static int noise_lut_samples(float feature_scale, float sin_theta,
                               bool octave_path) {
    const int samples =
        hs::clamp(static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * math::TWO_PI_F *
                                         feature_scale * sin_theta)),
                  LUT_MIN_SAMPLES, W);
    return octave_path ? (samples + OCTAVE_GRID - 1) / OCTAVE_GRID * OCTAVE_GRID
                       : samples;
  }

  /** @brief Evaluates the active ball fields using cached ring geometry. */
  HS_O3_FN float ball_field(const math::Vector &p, const int *ks, int n,
                            float theta) const {
    DominantFieldAccumulator accumulator;
    for (int j = 0; j < n; ++j) {
      const int k = ks[j];
      float f = bump_field_with_y(p, *ball_params[k], theta - ball_colat[k]);
      accumulator.add(f);
    }
    return accumulator.value();
  }

  struct HueMode {
    bool use_table = false;  ///< Sample the hue table instead of rotating.
    bool precompute = false; ///< Populate table cells before sampling.
    float domain = 0.0f;     ///< Hue-turn span covered by the table.
    bool cyclic = false;     ///< Wrap table lookups at the domain boundary.
  };

  /**
   * @brief Chooses how one ring's bake evaluates its hue rotation.
   * @param lut_n Bake columns for the ring.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param hue_extent Signed hue turns the displacement range covers.
   * @return Table sampling mode and domain.
   * @details Culled chunks bake no columns, so the precompute decision
   * amortizes the table over the visible chunks' columns rather than lut_n.
   */
  __attribute__((always_inline)) HueMode select_hue_mode(int lut_n,
                                                         uint32_t visible,
                                                         float hue_extent) {
    HueMode mode;
    mode.use_table = params.hue_scale != 0.0f && lut_n > 2 * HUE_TABLE_SIZE;
#if HS_ENABLE_TEST_ORACLES
    mode.use_table = mode.use_table && !force_exact_hue;
#endif
    if (!mode.use_table)
      return mode;
    int visible_samples = 0;
    Plot::for_each_chunk_span<BAKE_CHUNKS>(
        lut_n, visible,
        [&](int begin, int end)
            __attribute__((always_inline)) { visible_samples += end - begin; },
        [](int, int) __attribute__((always_inline)) {});
    mode.precompute = visible_samples > 2 * HUE_TABLE_SIZE;
    mode.cyclic = std::fabs(hue_extent) > 1.0f;
    mode.domain = mode.cyclic ? std::copysign(1.0f, hue_extent) : hue_extent;
#if HS_ENABLE_TEST_HOOKS
    ++hue_table_uses;
#endif
    return mode;
  }

  /** @brief Alternate sweep-window shifts and shading state. */
  struct SweepCache {
    static constexpr int CAPACITY = 3840;       ///< Packed alternate shifts.
    static constexpr int WORDS = (W + 31) / 32; ///< Mask words per ring.
    float *values;            ///< Second-window shifts in ring and knot order.
    uint32_t *masks;          ///< Per-ring alternate masks.
    uint16_t *starts;         ///< First alternate per mask word.
    HueMode *modes;           ///< Per-ring hue sampling mode.
    Color4 *colors;           ///< Second-window colors.
    float *alphas;            ///< Second-window fragment alpha.
    uint8_t *boundary;        ///< Current ring's conservative boundary set.
    int ring = 0;             ///< Current ring index.
    int count = 0;            ///< Retained alternates.
    bool table_ready = false; ///< Candidate table includes both windows.
  };
  static_assert(SweepCache::CAPACITY >= W,
                "DisplacementField: one ring's knots must fit the alternate "
                "buffer");

  /**
   * @brief A frame's endpoints for a draw that evaluates the displacement
   * fields at each knot's sweep time.
   */
  struct Sweep {
    /// Rotation vector carrying a point from where the stack held it at the
    /// previous endpoint to where it holds it at the current one.
    math::Vector turn;
    float noise_time[2]; ///< Noise time at the previous and current endpoint.
    float noise_amplitude[2];    ///< Noise amplitude at both endpoints.
    float opacity[2];            ///< Sprite opacity at both endpoints.
    SweepCache *cache = nullptr; ///< Single-pass alternate storage.
    float middle = 0.5f;         ///< Middle sweep fraction of the drawn window.
  };

  /**
   * @brief A ball's path over the frame, in the frame the stack holds at the
   * current endpoint.
   */
  struct BallPath {
    math::Vector from; ///< Center at the previous endpoint.
    math::Vector to;   ///< Center at the current endpoint.
    float axis_from;   ///< `from` along the stack axis.
    float axis_to;     ///< `to` along the stack axis.
    float reach_from;  ///< Footprint radius at the previous endpoint.
    float reach_to;    ///< Footprint radius at the current endpoint.
    float cos_reach;   ///< Cosine of the larger footprint radius.
    float cos_cover;   ///< Cosine of the cap covering the whole path.
  };

  /** @brief Sweep windows a swept draw splits the half-turn into. */
  static constexpr int SWEEP_WINDOWS = 2;
  static_assert(W % (2 * SWEEP_WINDOWS) == 0,
                "DisplacementField: sweep windows must tile the half-turn");

  /**
   * @brief Frame fraction at which a swept draw evaluates a ring point.
   * @param p The point where the stack holds it at the current endpoint.
   * @param middle Middle sweep fraction of the drawn window.
   * @return The sweep fraction of the column showing @p p, taken within half
   *         a frame of @p middle and clamped to the frame.
   * @details The window's columns read only points near them, so points
   * shown across an arm boundary from the window take the boundary's time.
   */
  HS_HOT_INLINE static float knot_fraction(const math::Vector &p,
                                           float middle) {
    float u = math::fast_atan2(p.z, p.x) * (1.0f / math::PI_F);
    u -= floorf(u);
    if (u - middle >= 0.5f)
      u -= 1.0f;
    else if (middle - u > 0.5f)
      u += 1.0f;
    return hs::clamp(u, 0.0f, 1.0f);
  }

  /** @brief Returns a knot's sweep fraction, or the reference column's time.
   * @param p Knot position at the current endpoint.
   * @param sweep Frame endpoints and window midpoint. */
  HS_HOT_INLINE float sample_fraction(const math::Vector &p,
                                      const Sweep &sweep) const {
#if HS_ENABLE_TEST_ORACLES
    if (force_column_reference)
      return sweep.middle;
#endif
    return knot_fraction(p, sweep.middle);
  }

  /**
   * @brief One ball's displacement of a ring point at a frame fraction.
   * @param p The point where the stack holds it at the current endpoint.
   * @param bp The ball's params.
   * @param path The ball's path over the frame.
   * @param theta Ring colatitude.
   * @param t Frame fraction.
   * @return The ball's field at @p p, with its center and footprint
   *         interpolated to @p t.
   */
  HS_HOT_INLINE static float ball_at(const math::Vector &p,
                                     const Animation::BumpParams &bp,
                                     const BallPath &path, float theta,
                                     float t) {
    const float reach = path.reach_from + (path.reach_to - path.reach_from) * t;
    if (reach <= 1e-3f || bp.amplitude <= 0.001f)
      return 0.0f;
    const math::Vector c = path.from + (path.to - path.from) * t;
    const float inverse = 1.0f / sqrtf(math::dot(c, c));
    const float cos_d = math::dot(p, c) * inverse;
    if (cos_d <= path.cos_reach)
      return 0.0f;
    const float d = math::fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
    if (!(d < reach))
      return 0.0f;
    const float colat = math::fast_acos(hs::clamp(
        (path.axis_from + (path.axis_to - path.axis_from) * t) * inverse, -1.0f,
        1.0f));
    return bump_field_profile(bp, reach, d, hs::clamp(theta - colat, -d, d));
  }

  /** @brief ball_field() with each ball at a frame fraction. */
  HS_O3_FN float ball_field_at(const math::Vector &p, const int *ks, int n,
                               float theta, float t) const {
    DominantFieldAccumulator accumulator;
    for (int j = 0; j < n; ++j) {
      const int k = ks[j];
      accumulator.add(ball_at(p, *ball_params[k], ball_paths[k], theta, t));
    }
    return accumulator.value();
  }

  /**
   * @brief Bakes every ring over the active ball pool, then rasterizes the
   * stack in one fused scan.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite fade multiplied into each fragment's alpha.
   * @details Under a partial clip, rings that cannot touch the clip and
   * invisible azimuth chunks skip the bake. Each column sees the stack
   * orientation at its sweep time; while a displacement field is present,
   * the fields, ring colors and opacity move with it (draw_swept()).
   */
  HS_HOT_FLASH_MEMBER void draw_rings(Canvas &canvas, float opacity) {
    HS_PROFILE(df_draw_rings);
    const int n_rings = params.num_rings;
    HS_CHECK(n_rings <= RING_SLOTS,
             "DisplacementField: Rings slider exceeds the baked-ring pool");
    math::Basis basis = math::make_basis(orientation.get(), STACK_AXIS);
    // The walk steps after this sprite draws, so the previous endpoint is the
    // orientation the last draw used.
    const math::Quaternion &current = orientation.get();
    const math::Quaternion previous =
        orientation_drawn ? drawn_orientation : current;
    drawn_orientation = current;
    ScratchScope column_guard(scratch_arena_a);
    const Temporal::ColumnTime column_time =
        Scan::prepare_column_basis<W, H>(scratch_arena_a, previous, current);
    Scan::ColumnTimeScope column_scope(canvas, column_time);

    const int n_balls = preset_index == BALLS_PRESET && !seed_balls_pending
                            ? balls.active_count()
                            : 0;
    Sweep sweep;
    const bool swept =
        prepare_sweep(previous, current, opacity, n_balls, sweep);
    orientation_drawn = true;
    if (swept) {
      draw_swept(canvas, basis, column_time, n_balls, sweep);
      return;
    }
    drawn_colors_valid = false;
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      ball_params[b] = &sp;
      const float cv = math::dot(basis.v, sp.center);
      ball_colat[b] = math::fast_acos(hs::clamp(cv, -1.0f, 1.0f));
      ball_reach[b] = sp.field_bound();
      ball_scale[b] = 2.0f / sp.radius;
      const float cu = math::dot(basis.u, sp.center);
      const float cw = math::dot(basis.w, sp.center);
      ball_cv[b] = cv;
      ball_rho[b] = sqrtf(cu * cu + cw * cw);
      const float a0 = atan2f(cw, cu);
      ball_azimuth[b] = a0 < 0.0f ? a0 + 2.0f * math::PI_F : a0;
    }
    draw_pass<false>(canvas, basis, column_time, opacity, n_balls, nullptr);
  }

  /**
   * @brief Reads the frame's endpoints and retains the current ones for the
   * next draw.
   * @param previous Stack orientation at the previous endpoint.
   * @param current Stack orientation at the current endpoint.
   * @param opacity Sprite opacity this frame.
   * @param n_balls Active balls drawn this frame.
   * @param sweep Receives the endpoints.
   * @return Whether the frame draws swept: temporal sampling is on, a
   *         previous draw supplies the first endpoint and a displacement field
   *         is present.
   */
  HS_FLASH_MEMBER bool prepare_sweep(const math::Quaternion &previous,
                                     const math::Quaternion &current,
                                     float opacity, int n_balls, Sweep &sweep) {
    const bool noise_live = noise_field.active_count() == 1;
    const float time = noise_live ? noise_field.active_params(0).time : 0.0f;
    const float amplitude =
        noise_live ? noise_field.active_params(0).amplitude : 0.0f;
    const bool primed = orientation_drawn;
    sweep.noise_time[0] = primed ? drawn_noise_time : time;
    sweep.noise_time[1] = time;
    sweep.noise_amplitude[0] = primed ? drawn_noise_amplitude : amplitude;
    sweep.noise_amplitude[1] = amplitude;
    sweep.opacity[0] = primed ? drawn_opacity : opacity;
    sweep.opacity[1] = opacity;
    drawn_noise_time = time;
    drawn_noise_amplitude = amplitude;
    drawn_opacity = opacity;
    if (!Temporal::temporal_sampling_enabled() || !primed)
      return false;
    const bool noise = noise_live && fmaxf(std::fabs(sweep.noise_amplitude[0]),
                                           std::fabs(amplitude)) > 0.001f;
    if (!noise && n_balls == 0)
      return false;
    math::Quaternion step = current * previous.conjugate();
    if (step.r < 0.0f)
      step = -step;
    const float half_sine = sqrtf(math::dot(step.v, step.v));
    const float angle = 2.0f * asinf(fminf(half_sine, 1.0f));
    if (half_sine > 1e-6f)
      sweep.turn = step.v * (angle / half_sine);
    else
      sweep.turn = step.v * 2.0f;
    return true;
  }

  /**
   * @brief Bakes column-time knots once and scans the stack in two windows.
   * @param canvas Render target.
   * @param basis Stack frame at the current endpoint.
   * @param column_time The frame's column basis.
   * @param n_balls Active balls drawn this frame.
   * @param sweep The frame's endpoints.
   * @details Boundary knots retain an alternate for the other arm pass.
   * Colors and opacity take each window's middle time. A ring whose boundary
   * knots would overfill the alternate buffer first scans the rings before
   * it, which keeps every pixel's blend order.
   */
  HS_HOT_FLASH_MEMBER void draw_swept(Canvas &canvas, const math::Basis &basis,
                                      const Temporal::ColumnTime &column_time,
                                      int n_balls, Sweep &sweep) {
#if HS_ENABLE_TEST_HOOKS
    primary_knots = duplicate_knots = alternate_noise_samples = primary_bakes =
        stack_flushes = 0;
#endif
    ScratchScope sweep_guard(scratch_arena_a);
    ball_paths = scratch_arena_a.allocate_n<BallPath>(MAX_BALLS);
    knot_times = scratch_arena_a.allocate_n<float>(W + 1);
    ring_colors = scratch_arena_a.allocate_n<Color4>(RING_SLOTS);
    ring_color_ready = scratch_arena_a.allocate_n<uint8_t>(RING_SLOTS);
    for (int i = 0; i < RING_SLOTS; ++i)
      ring_color_ready[i] = 0;
    if (drawn_ring_count != params.num_rings)
      drawn_colors_valid = false;
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      ball_params[b] = &sp;
      BallPath &path = ball_paths[b];
      path.from = Temporal::rotate_by(sp.previous_center, sweep.turn);
      path.to = sp.center;
      path.axis_from = math::dot(basis.v, path.from);
      path.axis_to = math::dot(basis.v, path.to);
      path.reach_from = sp.radius * sp.previous_envelope;
      path.reach_to = sp.radius * sp.envelope;
      const float reach = fmaxf(path.reach_from, path.reach_to);
      if (reach <= 1e-3f) {
        // No field at either endpoint, as before the fall's first step.
        ball_colat[b] = 0.0f;
        ball_reach[b] = -1.0f;
        continue;
      }
      path.cos_reach = cosf(fminf(reach, math::PI_F));
      // Every center on the path lies within half its chord, padded for the
      // arc, of the chord's midpoint.
      const math::Vector chord = path.to - path.from;
      const float cover = reach + 0.5f * sqrtf(math::dot(chord, chord)) + 1e-4f;
      path.cos_cover = cosf(fminf(cover, math::PI_F));
      const math::Vector mid = (path.from + path.to).normalized();
      const float cv = math::dot(basis.v, mid);
      ball_colat[b] = math::fast_acos(hs::clamp(cv, -1.0f, 1.0f));
      ball_reach[b] = cover;
      ball_scale[b] = 2.0f / sp.radius;
      const float cu = math::dot(basis.u, mid);
      const float cw = math::dot(basis.w, mid);
      ball_cv[b] = cv;
      ball_rho[b] = sqrtf(cu * cu + cw * cw);
      const float a0 = atan2f(cw, cu);
      ball_azimuth[b] = a0 < 0.0f ? a0 + 2.0f * math::PI_F : a0;
    }
    ScratchScope cache_guard(scratch_arena_b);
    SweepCache cache{
        scratch_arena_b.allocate_n<float>(SweepCache::CAPACITY),
        scratch_arena_a.make_n<uint32_t>(RING_SLOTS * SweepCache::WORDS),
        scratch_arena_a.make_n<uint16_t>(RING_SLOTS * SweepCache::WORDS),
        scratch_arena_a.make_n<HueMode>(RING_SLOTS),
        scratch_arena_a.allocate_n<Color4>(RING_SLOTS),
        scratch_arena_a.allocate_n<float>(RING_SLOTS),
        scratch_arena_a.allocate_n<uint8_t>(W + 1)};
    sweep.cache = &cache;
    active_cache = &cache;
#if HS_ENABLE_TEST_ORACLES
    if (force_window_bakes || force_column_reference)
      sweep.cache = active_cache = nullptr;
#endif
    if (sweep.cache)
      draw_pass<true>(canvas, basis, column_time, sweep.opacity[1], n_balls,
                      &sweep);
#if HS_ENABLE_TEST_ORACLES
    if (!sweep.cache) {
      const int windows = force_column_reference ? W / 2 : SWEEP_WINDOWS;
      for (int k = 0; k < windows; ++k) {
        sweep.middle = force_column_reference ? static_cast<float>(k) / windows
                                              : (k + 0.5f) / windows;
        RenderXWindow window(canvas, k * (W / 2) / windows,
                             (k + 1) * (W / 2) / windows);
        draw_pass<true>(canvas, basis, column_time, sweep.opacity[1], n_balls,
                        &sweep);
      }
    }
#endif
    active_cache = nullptr;
    sweep.cache = nullptr;
    for (int i = 0; i < RING_SLOTS; ++i) {
      if (ring_color_ready[i])
        drawn_colors[i] = ring_colors[i];
      drawn_color_ready[i] = ring_color_ready[i];
    }
    drawn_colors_valid = true;
    drawn_ring_count = params.num_rings;
    ball_paths = nullptr;
    knot_times = nullptr;
    ring_colors = nullptr;
    ring_color_ready = nullptr;
  }

  /**
   * @brief A ring's color at the current endpoint, from the palette once per
   * frame.
   * @param i Ring index.
   * @param n_rings Stack size.
   * @return The ring's color.
   */
  HS_FLASH_MEMBER Color4 current_ring_color(int i, int n_rings) {
    if (!ring_color_ready[i]) {
      ring_colors[i] =
          palette.get(math::wrap_t((i + 0.5f) / n_rings + color_spin));
      ring_color_ready[i] = 1;
    }
    return ring_colors[i];
  }

  /**
   * @brief Bakes and scans the rings of one pass over the stack.
   * @tparam SWEPT Whether knots are baked at their sweep times.
   * @param canvas Render target.
   * @param basis Stack frame at the current endpoint.
   * @param column_time The frame's column basis.
   * @param opacity Sprite opacity; a swept pass takes the window's own.
   * @param n_balls Active balls drawn this frame, with their caches filled.
   * @param sweep The frame's endpoints; null unless SWEPT.
   */
  template <bool SWEPT>
  HS_HOT_INLINE void draw_pass(Canvas &canvas, const math::Basis &basis,
                               const Temporal::ColumnTime &column_time,
                               float opacity, int n_balls, const Sweep *sweep) {
    const int n_rings = params.num_rings;
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the rasterizer's soft stroke cross-section past
    // the clip's margin and the stack's rotation over the frame.
    const float pad = 3.0f * math::PI_F / H + column_time.sweep;
    float noise_bound;
    if constexpr (SWEPT)
      noise_bound = fmaxf(std::fabs(sweep->noise_amplitude[0]),
                          std::fabs(sweep->noise_amplitude[1]));
    else
      noise_bound = noise_field.field_bound();

    const float noise_feature =
        noise_bound > 0.0f ? params.scale1 + params.scale2 : 0.0f;

    rings.begin_frame(SWEPT ? active_cache : nullptr);

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      float theta = radius * (math::PI_F / 2.0f);

      int n_local = 0;
      float band = 0.0f;
      float ball_feature = 0.0f;
      for (int b = 0; b < n_balls; ++b) {
        if (std::fabs(theta - ball_colat[b]) < ball_reach[b] + BALL_TOUCH_EPS) {
          ball_local[n_local++] = b;
          band = fmaxf(band, ball_reach[b]);
          ball_feature = fmaxf(ball_feature, ball_scale[b]);
        }
      }

      if (try_cull && !Plot::cap_may_touch_clip<H>(clip(), basis.v,
                                                   theta + band + noise_bound +
                                                       params.thickness + pad))
        continue;

      Color4 ring_color;
      float ring_opacity = opacity;
      if constexpr (SWEPT) {
        const Color4 now = current_ring_color(i, n_rings);
        const Color4 before =
            drawn_colors_valid && drawn_color_ready[i] ? drawn_colors[i] : now;
        const float middle = sweep->cache ? 0.25f : sweep->middle;
        ring_color = before.lerp(now, middle);
        if (sweep->cache) {
          sweep->cache->colors[i] = before.lerp(now, 0.75f);
          sweep->cache->alphas[i] =
              sweep->cache->colors[i].alpha *
              (sweep->opacity[0] +
               (sweep->opacity[1] - sweep->opacity[0]) * 0.75f) *
              params.alpha;
          sweep->cache->ring = i;
          sweep->cache->starts[i * SweepCache::WORDS] =
              static_cast<uint16_t>(sweep->cache->count);
        }
        ring_opacity = sweep->opacity[0] +
                       (sweep->opacity[1] - sweep->opacity[0]) * middle;
      } else {
        ring_color =
            palette.get(math::wrap_t((i + 0.5f) / n_rings + color_spin));
      }
      HueRotateBase hue_base = make_hue_rotate_base(ring_color);

      const typename RingPool::Pending rows = rings.next();
      float *slut = rows.shift_row;
      Pixel *hlut = rows.hue_row;

      int lut_n;
      if (band + noise_bound <= 0.0f) {
        // Flat rings still take a LUT: the fused scan needs every ring to
        // keep per-pixel blend order.
        Pixel flat = hue_rotate(hue_base, 0.0f).color;
        lut_n = LUT_MIN_SAMPLES;
        for (int x = 0; x <= lut_n; ++x) {
          slut[x] = 0.0f;
          hlut[x] = flat;
        }
      } else {
        float feature_scale = fmaxf(noise_feature, ball_feature);
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);

        float octave_amplitude;
        if constexpr (SWEPT)
          octave_amplitude = noise_bound;
        else
          octave_amplitude = noise_field.active_count() == 1
                                 ? noise_field.active_params(0).amplitude
                                 : 0.0f;
        const Animation::NoiseProductParams *octaves =
            noise_field.active_count() == 1 &&
                    std::fabs(octave_amplitude) > 0.001f
                ? &noise_field.active_params(0)
                : nullptr;
        lut_n = noise_lut_samples(feature_scale, sin_t, octaves != nullptr);

        uint32_t visible = CHUNK_MASK;
        if constexpr (SWEPT) {
          const float band_r = band + noise_bound + params.thickness + pad;
          const ColumnWindow &window = canvas.render_x_window();
          {
            HS_PROFILE(df_chunk_cull);
            visible = Plot::chunk_mask_where<BAKE_CHUNKS>(
                [&](const math::Vector &mid, float reach,
                    float sin_reach) HS_HOT_INLINE {
                  return (sweep->cache || Plot::cap_may_touch_window<W>(
                                              window, mid, reach, sin_reach)) &&
                         (!try_cull || Plot::cap_may_touch_clip<H>(
                                           clip(), mid, reach, sin_reach));
                },
                basis, theta, cos_t, sin_t, band_r, params.thickness, chunk_cos,
                chunk_sin);
          }
          if (!visible)
            continue;
        } else if (try_cull) {
          const float band_r = band + noise_bound + params.thickness + pad;
          {
            HS_PROFILE(df_chunk_cull);
            visible = Plot::visible_chunk_mask<H, BAKE_CHUNKS>(
                clip(), basis, theta, cos_t, sin_t, band_r, params.thickness,
                chunk_cos, chunk_sin);
          }
          if (!visible)
            continue;
        }

        HS_PROFILE(df_lut_bake);
#if HS_ENABLE_TEST_HOOKS
        ++primary_bakes;
#endif

        const float dphi = 2.0f * math::PI_F / lut_n;
        const float cos_d = cosf(dphi);
        const float sin_d = sinf(dphi);
        if constexpr (SWEPT) {
          if (sweep->cache) {
            SweepCache &cache = *sweep->cache;
            // A ring retains at most one alternate per knot.
            if (cache.count + lut_n > SweepCache::CAPACITY) {
              scan_swept_stack(canvas, cache);
              rings.begin_frame(&cache);
              cache.starts[i * SweepCache::WORDS] = 0;
              const typename RingPool::Pending fresh = rings.next();
              slut = fresh.shift_row;
              hlut = fresh.hue_row;
            }
            mark_boundary(basis, {cos_t, sin_t, cos_d, sin_d}, theta, lut_n,
                          band + noise_bound, column_time.sweep, cache);
          }
        }

        if (octaves) {
          bake_noise_octaves<SWEPT>(*octaves, basis,
                                    {cos_t, sin_t, cos_d, sin_d}, theta, lut_n,
                                    visible, n_local, slut, sweep);
        } else {
          bake_ball_spans<SWEPT>(basis, {cos_t, sin_t, cos_d, sin_d}, theta,
                                 lut_n, visible, n_local, slut, sweep);
        }

        // Culled chunks keep hlut stale; the padded `visible` mask keeps
        // rasterized pixels off their columns. Their shifts are zeroed:
        // DistortedRing scans every knot for its shift bounds.
        float max_shift = 0.0f;
        Plot::for_each_chunk_span<BAKE_CHUNKS>(
            lut_n, visible,
            [&](int begin, int end) __attribute__((always_inline)) {
              for (int x = begin; x < end; ++x)
                max_shift = fmaxf(max_shift, std::fabs(slut[x]));
            },
            [&](int begin, int end) __attribute__((always_inline)) {
              for (int x = begin; x < end; ++x)
                slut[x] = 0.0f;
            });

        const HueMode mode = select_hue_mode(
            lut_n, visible, (band + noise_bound) * params.hue_scale);
        if constexpr (SWEPT)
          if (sweep->cache)
            sweep->cache->modes[i] = mode;
        if (mode.use_table)
          hue_table.bind(hue_base, mode.domain, mode.cyclic);
        if (mode.precompute) {
          HS_PROFILE(df_hue_table_prep);
          hue_table.prepare(hue_table.cells(max_shift * params.hue_scale));
        }
        Pixel zero_hue;
        if (params.hue_scale == 0.0f)
          zero_hue = hue_rotate(hue_base, 0.0f).color;

        auto hue_for_shift = [&](float shift) HS_HOT_INLINE {
          if (params.hue_scale == 0.0f)
            return zero_hue;
          const float amount = std::fabs(shift) * params.hue_scale;
          if (mode.precompute)
            return hue_table.sample(amount);
          if (mode.use_table)
            return hue_table.sample_lazy(amount);
          return hue_rotate(hue_base, amount).color;
        };

        Plot::for_each_chunk_span<BAKE_CHUNKS>(
            lut_n, visible,
            [&](int begin, int end) __attribute__((always_inline)) {
              for (int x = begin; x < end; ++x)
                hlut[x] = hue_for_shift(slut[x]);
            },
            [](int, int) __attribute__((always_inline)) {});
        slut[lut_n] = slut[0];
        hlut[lut_n] = hlut[0];
      }

      if constexpr (SWEPT) {
        if (sweep->cache) {
          SweepCache &cache = *sweep->cache;
          int index = cache.starts[i * SweepCache::WORDS];
          for (int word = 0; word < SweepCache::WORDS; ++word) {
            cache.starts[i * SweepCache::WORDS + word] =
                static_cast<uint16_t>(index);
            index += std::popcount(cache.masks[i * SweepCache::WORDS + word]);
          }
        }
      }
      rings.commit(i, ring_color.alpha * ring_opacity * params.alpha, lut_n,
                   basis, radius, params.thickness, slut, lut_n, 0.0f, nullptr);
    }

    if (rings.size() == 0)
      return;

    if constexpr (SWEPT) {
      if (sweep->cache) {
        scan_swept_stack(canvas, *sweep->cache);
        return;
      }
    }
    RingShader ring_shader{&rings};
    HS_PROFILE(df_fused_scan);
    Scan::DistortedRingStack::draw<W, H>(
        filters, canvas, n_rings, rings.shapes(), rings.slot_map(),
        rings.size(), *candidates, ring_shader);
    rings.release();
  }

  /**
   * @brief Scans the committed rings in both sweep windows and empties the
   * pool and the alternate buffer.
   * @param canvas Render target.
   * @param cache The frame's alternate shifts, holding the committed rings'.
   */
  HS_HOT_FLASH_MEMBER void scan_swept_stack(Canvas &canvas, SweepCache &cache) {
    if (rings.size() > 0) {
      RingShader ring_shader{&rings};
      HS_PROFILE(df_fused_scan);
      for (int k = 0; k < SWEEP_WINDOWS; ++k) {
        if (k)
          apply_second_window(cache);
        RenderXWindow window(canvas, k * W / 4, (k + 1) * W / 4);
        Scan::DistortedRingStack::draw<W, H>(
            filters, canvas, params.num_rings, rings.shapes(), rings.slot_map(),
            rings.size(), *candidates, ring_shader);
      }
    }
    rings.release();
#if HS_ENABLE_TEST_HOOKS
    duplicate_knots += cache.count;
    ++stack_flushes;
#endif
    cache.count = 0;
    cache.table_ready = false;
  }

  /**
   * @brief Sets knot_visible[x] to whether knot x lies in a visible chunk.
   * @param lut_n Knot count.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   */
  __attribute__((always_inline)) void mark_visible_knots(int lut_n,
                                                         uint32_t visible) {
    int x = 0;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const int x_end = Plot::chunk_end<BAKE_CHUNKS>(c, lut_n);
      const uint8_t v = static_cast<uint8_t>((visible >> c) & 1u);
      for (; x < x_end; ++x)
        knot_visible[x] = v;
    }
  }

  /** @brief Ring and knot-step trigonometry. */
  struct KnotStep;
  /** @brief Recurrence over a ring's knot positions. */
  struct RingKnotWalk;

  /** @brief Marks knots whose supported geometry can cross an arm boundary.
   * @param basis Ring frame.
   * @param step Ring and knot-step trigonometry.
   * @param theta Ring colatitude.
   * @param n Knot count.
   * @param displacement Bound on polar displacement.
   * @param motion Frame rotation bound in radians.
   * @param cache Receives the boundary mask.
   * @details The supported azimuth interval covers the spline stencil and
   * bounded segment search. Its full displaced polar band is tested against
   * the arm plane widened by frame rotation and one column. */
  HS_HOT_FLASH_MEMBER void mark_boundary(const math::Basis &basis,
                                         KnotStep step, float theta, int n,
                                         float displacement, float motion,
                                         SweepCache &cache) {
    const float low = fmaxf(0.0f, theta - displacement - params.thickness);
    const float high =
        fminf(math::PI_F, theta + displacement + params.thickness);
    const float min_sine = fmaxf(0.001f, fminf(sinf(low), sinf(high)));
    const float cell = math::TWO_PI_F / n;
    const float limit =
        std::min(SDF::DistortedRing::MAX_SEARCH_CELLS, n / 2 + 1) + 1.0f;
    const float search =
        fminf(limit, ceilf(params.thickness / (min_sine * cell)) + 1.0f);
    constexpr int SUPPORT = 2 * std::max(OCTAVE1_STRIDE, OCTAVE2_STRIDE) - 1;
    const float delta = fminf(math::PI_F, (search + SUPPORT) * cell);
    const float epsilon =
        sinf(fminf(math::PI_F / 2, motion + math::TWO_PI_F / W)) + 1e-4f;
    const float axis = basis.v.z;
    const bool pole =
        (low == 0.0f || high == math::PI_F) && std::fabs(axis) <= epsilon;
    auto threshold = [&](float edge, bool minimum) HS_HOT_INLINE {
      auto at = [&](float angle) HS_HOT_INLINE {
        return (edge - axis * cosf(angle)) / fmaxf(1e-6f, sinf(angle));
      };
      float value =
          minimum ? fminf(at(low), at(high)) : fmaxf(at(low), at(high));
      const float cosine = axis / edge;
      if (std::fabs(cosine) < 1.0f) {
        const float angle = acosf(cosine);
        if (angle > low && angle < high)
          value = minimum ? fminf(value, at(angle)) : fmaxf(value, at(angle));
      }
      return value;
    };
    const float target_low = threshold(-epsilon, true);
    const float target_high = threshold(epsilon, false);
    const float radius = sqrtf(basis.u.z * basis.u.z + basis.w.z * basis.w.z);
    const float cosine = cosf(delta), sine = sinf(delta);
    RingKnotWalk walk{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
    for (int k = 0; k < n; ++k, walk.advance()) {
      const float center = basis.u.z * walk.cos_a + basis.w.z * walk.sin_a;
      const float tangent = basis.w.z * walk.cos_a - basis.u.z * walk.sin_a;
      float lo = center * cosine - std::fabs(tangent * sine);
      float hi = center * cosine + std::fabs(tangent * sine);
      if (center >= radius * cosine)
        hi = radius;
      if (center <= -radius * cosine)
        lo = -radius;
      cache.boundary[k] =
          pole || (lo <= target_high + 1e-4f && hi >= target_low - 1e-4f);
    }
  }

  /** @brief Retains a knot's opposite-window value and installs its early value.
   * @param k Knot index.
   * @param primary Value in the knot's own sweep window.
   * @param alternate Value in the other window.
   * @param time Knot's own sweep fraction.
   * @param row Early-window row.
   * @param cache Packed late-window values. */
  HS_HOT_INLINE void retain_alternate(int k, float primary, float alternate,
                                      float time, float *row,
                                      SweepCache &cache) {
    row[k] = primary;
    if (primary == alternate)
      return;
    cache.masks[cache.ring * SweepCache::WORDS + k / 32] |= uint32_t{1}
                                                            << (k % 32);
    row[k] = time < 0.5f ? primary : alternate;
    cache.values[cache.count++] = time < 0.5f ? alternate : primary;
  }

  /** @brief Installs late-window shifts and shades its knot rows.
   * @param cache Frame's packed alternate shifts and colors. */
  HS_HOT_FLASH_MEMBER void apply_second_window(const SweepCache &cache) {
    for (int i = 0; i < params.num_rings; ++i) {
      const int slot = rings.slot_of(i);
      if (slot < 0)
        continue;
      const int n = static_cast<int>(rings.lut_columns(slot));
      float *row = rings.shift_row(slot);
      Pixel *hues = rings.mutable_hue_row(slot);
      const uint32_t *mask = cache.masks + i * SweepCache::WORDS;
      int index = cache.starts[i * SweepCache::WORDS];
      const HueMode &mode = cache.modes[i];
      const HueRotateBase base = make_hue_rotate_base(cache.colors[i]);
      if (mode.use_table)
        hue_table.bind(base, mode.domain, mode.cyclic);
      for (int k = 0; k < n; ++k) {
        if (mask[k / 32] & (uint32_t{1} << (k % 32)))
          row[k] = cache.values[index++];
        const float amount = std::fabs(row[k]) * params.hue_scale;
        hues[k] = mode.use_table ? hue_table.sample_lazy(amount)
                                 : hue_rotate(base, amount).color;
      }
      row[n] = row[0];
      hues[n] = hues[0];
      rings.set_alpha(slot, cache.alphas[i]);
    }
  }

  /** @brief A ring's colatitude and knot-step trig, passed in VFP registers. */
  struct KnotStep {
    float cos_t; /**< Cosine of the ring colatitude. */
    float sin_t; /**< Sine of the ring colatitude. */
    float cos_d; /**< Cosine of one knot cell's azimuth step. */
    float sin_d; /**< Sine of one knot cell's azimuth step. */
  };

  /**
   * @brief Steps a ring's knot positions around its azimuth by the
   * angle-addition recurrence, starting at knot 0.
   */
  struct RingKnotWalk {
    const math::Basis &basis; /**< Ring frame; basis.v is the stack axis. */
    float cos_t;              /**< Cosine of the ring colatitude. */
    float sin_t;              /**< Sine of the ring colatitude. */
    float cos_d;              /**< Cosine of one knot cell's azimuth step. */
    float sin_d;              /**< Sine of one knot cell's azimuth step. */
    float cos_a = 1.0f;       /**< Cosine of the current knot's azimuth. */
    float sin_a = 0.0f;       /**< Sine of the current knot's azimuth. */

    /** @brief Unit position of the current knot. */
    __attribute__((always_inline)) math::Vector position() const {
      return ((basis.v * cos_t) +
              ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t)
          .normalized();
    }

    /** @brief Moves to the next knot. */
    __attribute__((always_inline)) void advance() {
      const float next_cos = cos_a * cos_d - sin_a * sin_d;
      sin_a = sin_a * cos_d + cos_a * sin_d;
      cos_a = next_cos;
    }
  };

  /**
   * @brief Bakes one ring's centerline shifts with each noise octave sampled on
   *        its own knot grid.
   * @param np The noise field's single active entity.
   * @param basis Ring frame; basis.v is the stack axis.
   * @param step Trig of theta and of one knot cell's azimuth step.
   * @param theta Ring colatitude.
   * @param lut_n Knot count, a multiple of OCTAVE_GRID.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param n_local Balls that can reach the ring (ball_local).
   * @param slut Receives the shift of every knot in a visible chunk.
   * @details The two octaves use OCTAVE1_STRIDE and OCTAVE2_STRIDE knot spacing
   * and fill between samples with a Catmull-Rom spline over four grid knots, so a
   * grid knot is evaluated when any visible knot lies within two strides of it.
   * Knot positions come from the same azimuth recurrence as the exact bake.
   * A swept bake samples each knot where the stack held it at the knot's
   * time, at the noise time and amplitude of that time.
   * @tparam SWEPT Whether knots take their sweep times.
   * @param sweep The frame's endpoints; read only when SWEPT.
   */
  template <bool SWEPT = false>
  HS_HOT_FLASH_MEMBER void
  bake_noise_octaves(const Animation::NoiseProductParams &np,
                     const math::Basis &basis, KnotStep step, float theta,
                     int lut_n, uint32_t visible, int n_local, float *slut,
                     const Sweep *sweep = nullptr) {
    constexpr int D1 = OCTAVE1_STRIDE;
    constexpr int D2 = OCTAVE2_STRIDE;
    mark_visible_knots(lut_n, visible);
    auto wrap = [lut_n](int k) HS_HOT_INLINE {
      return k < 0 ? k + lut_n : k >= lut_n ? k - lut_n : k;
    };
    // knot_near[x]: a visible knot lies within the widest spline reach of x.
    constexpr int REACH = 2 * (D1 > D2 ? D1 : D2) - 1;
    int in_window = 0;
    for (int o = -REACH; o <= REACH; ++o)
      in_window += knot_visible[wrap(o)];
    for (int x = 0; x < lut_n; ++x) {
      knot_near[x] = in_window > 0;
      in_window +=
          knot_visible[wrap(x + REACH + 1)] - knot_visible[wrap(x - REACH)];
    }

    HS_PROFILE(df_octave_noise);
    RingKnotWalk knots{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
    for (int x = 0; x < lut_n; ++x, knots.advance()) {
      const bool vis = knot_visible[x] != 0;
#if HS_ENABLE_TEST_HOOKS
      primary_knots += vis;
#endif
      const bool near = knot_near[x] != 0;
      const bool g1 = near && x % D1 == 0;
      const bool g2 = D2 == 1 ? vis : near && x % D2 == 0;
      if (vis || g1 || g2) {
        math::Vector p = knots.position();
        if constexpr (SWEPT) {
          const float t = sample_fraction(p, *sweep);
          knot_times[x] = t;
          const math::Vector q =
              Temporal::rotate_by(p, sweep->turn * (t - 1.0f));
          const float time = sweep->noise_time[0] +
                             (sweep->noise_time[1] - sweep->noise_time[0]) * t;
          if (g1)
            octave1[x] = np.noise.GetNoise(q.x * np.scale1, q.y * np.scale1,
                                           q.z * np.scale1 + time);
          if (g2)
            octave2[x] = np.noise.GetNoise(
                q.x * np.scale2 + Animation::NoiseProductParams::OCTAVE2_OFFSET,
                q.y * np.scale2, q.z * np.scale2 + time);
          if (sweep->cache && (g1 || g2)) {
            knot_pos[x].x = g1 ? octave1[x] : 0.0f;
            knot_pos[x].y = g2 ? octave2[x] : 0.0f;
            if (sweep->cache->boundary[x] && (t < 0.25f || t >= 0.75f)) {
#if HS_ENABLE_TEST_HOOKS
              alternate_noise_samples +=
                  static_cast<int>(g1) + static_cast<int>(g2);
#endif
              const float other = t < 0.5f ? 1.0f : 0.0f;
              const math::Vector aq =
                  Temporal::rotate_by(p, sweep->turn * (other - 1.0f));
              const float at =
                  sweep->noise_time[0] +
                  (sweep->noise_time[1] - sweep->noise_time[0]) * other;
              if (g1)
                knot_pos[x].x = np.noise.GetNoise(
                    aq.x * np.scale1, aq.y * np.scale1, aq.z * np.scale1 + at);
              if (g2)
                knot_pos[x].y = np.noise.GetNoise(
                    aq.x * np.scale2 +
                        Animation::NoiseProductParams::OCTAVE2_OFFSET,
                    aq.y * np.scale2, aq.z * np.scale2 + at);
            }
          }
          if (vis)
            slut[x] = ball_field_at(p, ball_local, n_local, theta, t);
        } else {
          if (g1)
            octave1[x] = np.noise.GetNoise(p.x * np.scale1, p.y * np.scale1,
                                           p.z * np.scale1 + np.time);
          if (g2)
            octave2[x] = np.noise.GetNoise(
                p.x * np.scale2 + Animation::NoiseProductParams::OCTAVE2_OFFSET,
                p.y * np.scale2, p.z * np.scale2 + np.time);
          if (vis)
            slut[x] = ball_field(p, ball_local, n_local, theta);
        }
      }
    }

    bool late = false;
    auto sample = [&](const float *oct, int d, int k) HS_HOT_INLINE {
      auto value = [&](int index) HS_HOT_INLINE {
        if constexpr (SWEPT) {
          if (sweep->cache && ((late && knot_times[index] < 0.25f) ||
                               (!late && knot_times[index] >= 0.75f)))
            return oct == octave1 ? knot_pos[index].x : knot_pos[index].y;
        }
        return oct[index];
      };
      const int r = k % d;
      if (r == 0)
        return value(k);
      const int b = k - r;
      const float t = static_cast<float>(r) / d;
      const float p0 = value(wrap(b - d)), p1 = value(b),
                  p2 = value(wrap(b + d)), p3 = value(wrap(b + 2 * d));
      return p1 + 0.5f * t *
                      (p2 - p0 +
                       t * (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3 +
                            t * (3.0f * (p1 - p2) + p3 - p0)));
    };
    RingKnotWalk final_walk{basis, step.cos_t, step.sin_t, step.cos_d,
                            step.sin_d};
    for (int x = 0; x < lut_n; ++x, final_walk.advance()) {
      if (!knot_visible[x])
        continue;
      float amplitude = np.amplitude;
      if constexpr (SWEPT)
        amplitude = sweep->noise_amplitude[0] +
                    (sweep->noise_amplitude[1] - sweep->noise_amplitude[0]) *
                        knot_times[x];
      if constexpr (SWEPT)
        late = knot_times[x] >= 0.5f;
      const float primary =
          slut[x] + amplitude * sample(octave1, D1, x) * sample(octave2, D2, x);
      if constexpr (SWEPT) {
        if (sweep->cache && sweep->cache->boundary[x]) {
          late = !late;
          const float t = knot_times[x];
          const float other = t < 0.25f ? 1.0f : t >= 0.75f ? 0.0f : t;
          const float amp =
              sweep->noise_amplitude[0] +
              (sweep->noise_amplitude[1] - sweep->noise_amplitude[0]) * other;
          const float ball =
              n_local ? ball_field_at(final_walk.position(), ball_local,
                                      n_local, theta, other)
                      : 0.0f;
          const float alternate =
              ball + amp * sample(octave1, D1, x) * sample(octave2, D2, x);
          retain_alternate(x, primary, alternate, t, slut, *sweep->cache);
          continue;
        }
      }
      slut[x] = primary;
    }
  }

  /**
   * @brief Bakes one ring's centerline shifts from the balls, each evaluated
   *        only across the knots its cap can cover.
   * @param basis Ring frame; basis.v is the stack axis.
   * @param step Trig of theta and of one knot cell's azimuth step.
   * @param theta Ring colatitude.
   * @param lut_n Knot count.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param n_local Balls that can reach the ring (ball_local).
   * @param slut Receives the shift of every knot in a visible chunk.
   * @details On the ring, a ball's cap test dot(p, center) > cos_radius reads
   * cos_t * cv + sin_t * rho * cos(a - a0) > cos_radius, which holds on one
   * azimuth arc around a0. Knots off every padded arc get exactly the zero a
   * rejected ball adds, and each knot folds its balls in ball_local order, so
   * the shifts match evaluating every ball at every knot. A swept bake moves
   * each ball along its path to the knot's time, inside the arc of the cap
   * covering the path.
   * @tparam SWEPT Whether knots take their sweep times.
   * @param sweep The frame's endpoints; read only when SWEPT.
   */
  template <bool SWEPT = false>
  HS_HOT_FLASH_MEMBER void
  bake_ball_spans(const math::Basis &basis, KnotStep step, float theta,
                  int lut_n, uint32_t visible, int n_local, float *slut,
                  const Sweep *sweep = nullptr) {
    // Covers the recurrence's drift from the analytic knot azimuth and the
    // float slop in the cap dot product.
    constexpr float COS_MARGIN = 2e-4f;
    float *num = octave1;
    float *den = octave2;
    mark_visible_knots(lut_n, visible);
    RingKnotWalk knots{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
    for (int x = 0; x < lut_n; ++x, knots.advance()) {
      if (knot_visible[x]) {
#if HS_ENABLE_TEST_HOOKS
        ++primary_knots;
#endif
        knot_pos[x] = knots.position();
        if constexpr (SWEPT)
          knot_times[x] = sample_fraction(knot_pos[x], *sweep);
        num[x] = 0.0f;
        den[x] = 0.0f;
        if constexpr (SWEPT) {
          if (sweep->cache)
            alternate_num[x] = alternate_den[x] = 0.0f;
        }
      }
    }

    const float cells_per_radian = lut_n / (2.0f * math::PI_F);
    for (int j = 0; j < n_local; ++j) {
      const int k = ball_local[j];
      const Animation::BumpParams &bp = *ball_params[k];
      const float a_term = step.cos_t * ball_cv[k];
      const float b_term = step.sin_t * ball_rho[k];
      int first = 0;
      int count = lut_n;
      float cos_radius = bp.cos_radius;
      if constexpr (SWEPT)
        cos_radius = ball_paths[k].cos_cover;
      if (b_term > COS_MARGIN) {
        const float q = (cos_radius - COS_MARGIN - a_term) / b_term;
        if (q >= 1.0f)
          continue;
        if (q > -1.0f) {
          const float half = acosf(q) * cells_per_radian;
          const float center = ball_azimuth[k] * cells_per_radian;
          first = static_cast<int>(floorf(center - half)) - 1;
          count = static_cast<int>(ceilf(center + half)) + 1 - first + 1;
          if (count > lut_n)
            count = lut_n;
          first = first % lut_n;
          if (first < 0)
            first += lut_n;
        }
      }
      const float y = theta - ball_colat[k];
      for (int i = 0, xi = first; i < count; ++i) {
        if (knot_visible[xi]) {
          float f;
          if constexpr (SWEPT)
            f = ball_at(knot_pos[xi], bp, ball_paths[k], theta, knot_times[xi]);
          else
            f = bump_field_with_y(knot_pos[xi], bp, y);
          DominantFieldAccumulator::accumulate(num[xi], den[xi], f);
          if constexpr (SWEPT) {
            const float t = knot_times[xi];
            if (sweep->cache && sweep->cache->boundary[xi] &&
                (t < 0.25f || t >= 0.75f)) {
              const float other = t < 0.5f ? 1.0f : 0.0f;
              const float alternate =
                  ball_at(knot_pos[xi], bp, ball_paths[k], theta, other);
              DominantFieldAccumulator::accumulate(
                  alternate_num[xi], alternate_den[xi], alternate);
            }
          }
        }
        if (++xi == lut_n)
          xi = 0;
      }
    }

    for (int x = 0; x < lut_n; ++x)
      if (knot_visible[x]) {
        slut[x] = DominantFieldAccumulator::resolve(num[x], den[x]) +
                  noise_field.field(knot_pos[x]);
        if constexpr (SWEPT) {
          const float t = knot_times[x];
          if (sweep->cache && sweep->cache->boundary[x] &&
              (t < 0.25f || t >= 0.75f)) {
            const float alternate = DominantFieldAccumulator::resolve(
                                        alternate_num[x], alternate_den[x]) +
                                    noise_field.field(knot_pos[x]);
            retain_alternate(x, slut[x], alternate, t, slut, *sweep->cache);
          }
        }
      }
  }

  /**
   * @brief Builds a fresh random palette for the next wipe.
   * @details Draws a base hue from the shared RNG; hues may repeat.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::displacement_field(
        EffectPaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Spawns one falling ball with a random meridian, footprint, and
   * speed drawn from the Speed Min/Max sliders; dropped safely if the ball pool
   * or the timeline is full.
   * @param progress Initial fall progress; zero starts at the north pole.
   * @param canvas Active frame canvas for nonzero initial progress.
   * @details Only the first pool-full drop of each ball phase is logged.
   */
  HS_COLD_MEMBER void spawn_ball(float progress = 0.0f,
                                 Canvas *canvas = nullptr) {
    balls.template_params.radius =
        hs::rand_f(std::min(params.ball_min, params.ball_max),
                   std::max(params.ball_min, params.ball_max));
    float speed =
        hs::rand_f(std::min(params.ball_speed_min, params.ball_speed_max),
                   std::max(params.ball_speed_min, params.ball_speed_max));
    int fall_frames = std::max(2, static_cast<int>(BALL_RATE_FPS / speed));
    if (balls.active_count() >= MAX_BALLS && !logged_pool_full) {
      logged_pool_full = true;
      hs::log("DisplacementField: ball pool full, dropping spawn");
    }
    auto *fall = balls.spawn(0, orientation, STACK_AXIS,
                             hs::rand_f(0.0f, 2.0f * math::PI_F), fall_frames);
    if (fall && progress > 0.0f) {
      HS_CHECK(canvas != nullptr, "DisplacementField: preview needs a canvas");
      for (int frame = 0; frame < static_cast<int>(progress * fall_frames);
           ++frame)
        fall->step(*canvas);
      fall->settle();
    }
  }

  /** @brief Schedules a cancellable noise gain transition.
   * @param target Final noise gain.
   * @param pausable Whether preset pause holds the transition. */
  void schedule_gain(float target, bool pausable) {
    timeline.add_get(0,
                     Animation::Transition(master_gain, target,
                                           NOISE_FADE_FRAMES,
                                           math::ease_in_out_sin),
                     Timeline::Pin::UNPINNED,
                     pausable ? &anims_paused : nullptr, &phase_events);
  }

  /** @brief Applies an authored row and starts its displacement mode.
   * @param change Accepted preset change request.
   * @return True after installing the target defaults. */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    params = PRESETS[change.to];
    if (change.origin == PresetChangeOrigin::AUTOMATIC) {
      if (change.to == NOISE_PRESET)
        enter_noise();
      else
        enter_balls();
    } else {
      timeline.cancel_owner(&phase_events);
      timeline.cancel_owner(&balls);
      orientation_drawn = false;
      master_gain = change.to == NOISE_PRESET ? 1.0f : 0.0f;
      noise_hold = NOISE_HOLD_FRAMES;
      ball_phase_left = BALL_PHASE_FRAMES;
      spawn_cooldown = 0.0f;
      logged_pool_full = false;
      seed_balls_pending = change.to == BALLS_PRESET;
    }
    return true;
  }

  /** @brief Cancels a departure and restarts the tuned preset dwell. */
  HS_COLD_MEMBER void parameter_written() override {
    if (preset_index == NOISE_PRESET && noise_hold == 0) {
      timeline.cancel_owner(&phase_events);
      master_gain = 1.0f;
    }
    noise_hold = NOISE_HOLD_FRAMES;
    ball_phase_left = BALL_PHASE_FRAMES;
  }

  /**
   * @brief Enters the noise phase and fades it in from zero.
   */
  HS_COLD_MEMBER void enter_noise() {
    master_gain = 0.0f;
    noise_hold = NOISE_HOLD_FRAMES;
    schedule_gain(1.0f, false);
  }

  /**
   * @brief Enters the ball phase, opening a fresh spawning window.
   */
  HS_COLD_MEMBER void enter_balls() {
    ball_phase_left = BALL_PHASE_FRAMES;
    spawn_cooldown = 0;
    logged_pool_full = false;
  }

  /**
   * @brief Rolls the palette toward a freshly generated one via a ColorWipe.
   */
  void roll_palette() {
    palette_start = palette.snapshot();
    palette_target = make_palette().snapshot();
    timeline.add(0,
                 Animation::ColorWipe(palette, palette_start, palette_target,
                                      PALETTE_WIPE_FRAMES, math::ease_linear));
  }

  FastNoiseLite walk_noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  // Each in-flight ball is one timeline event; a saturated pool drops spawns.
  static constexpr int MAX_BALLS =
      56; /**< Concurrent falling-ball pool slots. */
  static constexpr int RESERVED_EVENTS =
      6; /**< Non-ball timeline events the effect schedules. */
  static_assert(MAX_BALLS + RESERVED_EVENTS <= Timeline::MAX_EVENTS,
                "DisplacementField: a full ball pool plus the effect's own "
                "events exceeds the shared timeline budget");
  static constexpr int BALL_PHASE_FRAMES =
      900; /**< Ball-phase spawning window. */
  static constexpr float BALL_RATE_FPS =
      60.0f; /**< Frames per Ball Rate / Speed slider unit. */
  static constexpr float BALL_DRAPE_PER_AMPLITUDE =
      4.0f; /**< Drape gain per Ball Amp unit. */
  // Phase departures run before their timeline completion is observed.
  static constexpr int PHASE_HANDOFF_FRAMES = 3;
  static constexpr int NOISE_FADE_FRAMES =
      150; /**< Noise amplitude ramp on each phase handoff. */
  static constexpr int NOISE_HOLD_FRAMES =
      600; /**< Full-noise dwell before fading out into the next ball phase. */
  /** @brief Un-oriented ring stack and ball displacement axis. */
  static constexpr math::Vector STACK_AXIS = math::X_AXIS;

  BallDropTransformer<MAX_BALLS>
      balls; /**< Falling-ball displacement fields. */
  NoiseProductTransformer<1>
      noise_field; /**< Two-octave noise displacement field. */

  GenerativePalette
      palette; /**< Active palette (mutated by an in-flight ColorWipe). */
  GenerativePalette::Snapshot palette_start;  /**< Current wipe's start. */
  GenerativePalette::Snapshot palette_target; /**< Current wipe's target. */
  math::Orientation<> orientation;
  math::Quaternion drawn_orientation; ///< Orientation of the last ring draw.
  /// Whether a ring draw since the last history reset supplies the previous
  /// endpoint.
  bool orientation_drawn = false;
  float drawn_noise_time = 0.0f;      ///< Noise time of the last ring draw.
  float drawn_noise_amplitude = 0.0f; ///< Noise amplitude of the last draw.
  float drawn_opacity = 0.0f;         ///< Sprite opacity of the last draw.
  Color4 *drawn_colors =
      nullptr; ///< RING_SLOTS ring colors of the last swept draw.
  uint8_t *drawn_color_ready =
      nullptr; ///< RING_SLOTS flags: the last swept draw colored the ring.
  bool drawn_colors_valid = false; ///< The last draw was swept.
  int drawn_ring_count = 0;        ///< Stack size of the last swept draw.
  // Frame scratch of a swept draw, null outside it.
  BallPath *ball_paths = nullptr; ///< MAX_BALLS active-ball paths.
  float *knot_times = nullptr;    ///< W + 1 frame fractions of baked knots.
  Color4 *ring_colors = nullptr;  ///< RING_SLOTS current ring colors.
  uint8_t *ring_color_ready =
      nullptr; ///< RING_SLOTS flags: ring_colors holds the ring's color.

  /** Noise row in PRESETS and PRESET_IDS. */
  static constexpr size_t NOISE_PRESET = 0;
  /** Balls row in PRESETS and PRESET_IDS. */
  static constexpr size_t BALLS_PRESET = 1;
  /** Visible falls seeded by a manual Balls selection. */
  static constexpr int PREVIEW_BALLS = 8;
  bool phase_events = false; ///< Lifetime token for phase gain animations.
  bool seed_balls_pending =
      false; ///< Selected balls seed after old events retire.
  int ball_phase_left =
      BALL_PHASE_FRAMES; /**< Frames left in this ball phase's spawning window. */
  float spawn_cooldown = 0.0f; /**< Frames until the next ball spawn. */
  bool logged_pool_full =
      false; /**< Pool-full log latch; cleared on each ball phase. */
  float master_gain =
      0.0f; /**< Noise fade envelope in [0, 1]; gates the noise field, animated by Transitions. */
  int noise_hold =
      0; /**< Frames until the noise phase begins fading back out. */

  static constexpr int PALETTE_CYCLE_FRAMES =
      180; /**< Palette rollover period. */
  static constexpr int PALETTE_WIPE_FRAMES =
      168; /**< Wipe duration, shorter than the rollover period. */
  // The wipe is armed mid-step and first steps on the next frame, so it spans
  // PALETTE_WIPE_FRAMES + 1 frames.
  static_assert(PALETTE_CYCLE_FRAMES > PALETTE_WIPE_FRAMES + 1,
                "DisplacementField: a rollover firing mid-wipe would clobber "
                "the snapshots the live ColorWipe still references");
  static constexpr float COLOR_SPIN_RATE =
      0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float LUT_SAMPLES_PER_UNIT =
      8.0f; /**< Bake columns per feature-space unit of ring circumference. */
  static constexpr int LUT_MIN_SAMPLES =
      16; /**< Bake-column floor for tiny/low-scale rings. */
  static constexpr int OCTAVE1_STRIDE =
      2; /**< Knots between envelope-octave noise samples. */
  static constexpr int OCTAVE2_STRIDE =
      2; /**< Knots between detail-octave noise samples. */
  static constexpr int OCTAVE_GRID = std::lcm(
      OCTAVE1_STRIDE, OCTAVE2_STRIDE); /**< Knot-count multiple holding both
                                          octave grids. */
  static_assert(W % OCTAVE_GRID == 0 && LUT_MIN_SAMPLES % OCTAVE_GRID == 0,
                "rounding lut_n up to the octave grid must stay within W");
  static constexpr float THICKNESS_PX = math::RADIANS_PER_COLUMN<
      W>; /**< One pixel of azimuth in ring-space; the Thickness
                            slider range is authored in multiples of it. */
  static constexpr int HUE_TABLE_SIZE =
      64; /**< Hue-turn interpolation cells per ring. */
  static constexpr int BAKE_CHUNKS =
      16; /**< Azimuth chunks clip-tested during the bake. */
  static constexpr uint32_t CHUNK_MASK =
      (1u << BAKE_CHUNKS) - 1; /**< All-chunks-visible bake mask. */
  static_assert(LUT_MIN_SAMPLES >= BAKE_CHUNKS,
                "the one-chunk visibility pad needs at least one LUT column "
                "per chunk");

  float color_spin =
      0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  static constexpr int RING_SLOTS = std::min(
      72, H); /**< Baked-ring pool capacity and Rings slider maximum. */

  /**
   * @brief The frame's drawn rings: placement-built DistortedRings with their
   * shift and hue rows, shading values and the ring-to-slot map.
   * @details A frame runs begin_frame(), then next() and commit() per drawn
   * ring, then release(). next() changes no state, so a ring abandoned after
   * its bake costs nothing; only commit() constructs a ring and counts it.
   */
  class RingPool {
    union RingSlot {
      char empty;
      SDF::DistortedRing ring;
      RingSlot() : empty{} {}
      ~RingSlot() {}
    };
    static_assert(sizeof(RingSlot) == sizeof(SDF::DistortedRing));

  public:
    static constexpr int SLOTS = RING_SLOTS; /**< Slot capacity. */
    static constexpr int ROW = W + 1;    /**< Entries per shift or hue row. */
    static constexpr int8_t CULLED = -1; /**< slot_of() for an undrawn ring. */
    static_assert(SLOTS <= INT8_MAX,
                  "the ring-to-slot map is int8_t with CULLED = -1; a larger "
                  "pool wraps slot indices negative");
    using ShapeStorage = std::array<RingSlot, SLOTS>; /**< Ring storage. */
    /** @brief Persistent bytes init_storage() takes, alignment included. */
    static constexpr size_t BYTES =
        SLOTS * ROW * (sizeof(float) + sizeof(Pixel)) +
        SLOTS * (2 * sizeof(float) + sizeof(int8_t)) + sizeof(ShapeStorage) +
        3 * alignof(float) + alignof(Pixel) + alignof(int8_t) +
        alignof(ShapeStorage);

    /** @brief The rows of the next free slot. */
    struct Pending {
      float *shift_row; /**< ROW centerline shifts; entry lut_n repeats entry
                           0 to close the polyline. */
      Pixel *hue_row;   /**< ROW hue-rotated ring colors. */
    };

    /** @brief Constructed rings indexed by slot. */
    struct ShapeView {
      ShapeStorage *storage; ///< Pool slots; not owned.
      const RingPool *pool;  ///< Slot map and alternate storage.
      /** @brief Bounds both sweep-window values of a knot.
       * @param slot Constructed ring slot.
       * @param knot Knot index.
       * @param lo Receives the smaller shift.
       * @param hi Receives the larger shift. */
      void knot_range(int slot, int knot, float &lo, float &hi) const {
        lo = hi = (*storage)[slot].ring.knots[knot];
        if (!pool->cache)
          return;
        const int ring = pool->ring_by_slot[slot];
        const auto &cache = *pool->cache;
        const uint32_t *mask = cache.masks + ring * SweepCache::WORDS;
        const uint32_t bit = uint32_t{1} << (knot % 32);
        if (!(mask[knot / 32] & bit))
          return;
        int index = cache.starts[ring * SweepCache::WORDS + knot / 32];
        index += std::popcount(mask[knot / 32] & (bit - 1));
        lo = fminf(lo, cache.values[index]);
        hi = fmaxf(hi, cache.values[index]);
      }
      /** @return Whether the frame's shared candidate table is filled. */
      bool table_ready() const {
        return pool->cache && pool->cache->table_ready;
      }
      /** @brief Retains the candidate table across sweep windows. */
      void table_built() const {
        if (pool->cache)
          pool->cache->table_ready = true;
      }
      /**
       * @brief Ring constructed in a slot.
       * @param index Slot index in [0, SLOTS).
       * @return The slot's ring.
       */
      SDF::DistortedRing &operator[](size_t index) const {
        return (*storage)[index].ring;
      }
    };

    /** @brief Allocates the slots from @p arena; hue rows start zeroed. */
    HS_COLD_MEMBER void init_storage(Arena &arena) {
      shift_rows = arena.allocate_n<float>(SLOTS * ROW);
      // Zeroed, not raw: a culled azimuth chunk leaves its columns unbaked.
      hue_rows = arena.make_n<Pixel>(SLOTS * ROW);
      alphas = arena.allocate_n<float>(SLOTS);
      columns = arena.allocate_n<float>(SLOTS);
      slot_by_ring = arena.allocate_n<int8_t>(SLOTS);
      storage = arena.make<ShapeStorage>();
    }

    /** @brief Starts a frame with every ring culled; the pool must be empty. */
    void begin_frame(SweepCache *frame_cache = nullptr) {
      cache = frame_cache;
      HS_CHECK(size() == 0,
               "DisplacementField: ring pool begins a frame with live rings");
      for (int i = 0; i < SLOTS; ++i)
        slot_by_ring[i] = CULLED;
    }

    /** @brief The next free slot's rows; repeated calls return the same rows. */
    __attribute__((always_inline)) Pending next() const {
      HS_CHECK(count < SLOTS, "DisplacementField: ring pool is full");
      return {shift_rows + count * ROW, hue_rows + count * ROW};
    }

    /**
     * @brief Constructs a ring in the slot next() returned and maps @p ring
     * to it; the pool must not be full.
     * @param ring Ring index in the stack.
     * @param frag_alpha Fragment alpha the shader applies to the ring.
     * @param lut_n Bake columns in the slot's rows.
     * @param ring_args DistortedRing constructor arguments.
     */
    template <typename... RingArgs>
    HS_HOT_FLASH_MEMBER void commit(int ring, float frag_alpha, int lut_n,
                                    RingArgs &&...ring_args) {
      assert(count < SLOTS && "DisplacementField: ring pool is full");
      const int s = count;
      ::new (static_cast<void *>(&(*storage)[s].ring))
          SDF::DistortedRing(std::forward<RingArgs>(ring_args)...);
      columns[s] = static_cast<float>(lut_n);
      alphas[s] = frag_alpha;
      slot_by_ring[ring] = static_cast<int8_t>(s);
      ring_by_slot[s] = ring;
      if (cache) {
        auto &shape = (*storage)[s].ring;
        float low = 0.0f, high = 0.0f;
        for (int k = 0; k < lut_n; ++k) {
          float lo, hi;
          shapes().knot_range(s, k, lo, hi);
          low = fminf(low, lo);
          high = fmaxf(high, hi);
        }
        shape.max_distortion = fmaxf(std::fabs(low), std::fabs(high));
        shape.max_thickness = shape.thickness + shape.max_distortion;
        SDF::band_edge_cosines(shape.target_angle + low - shape.thickness,
                               shape.target_angle + high + shape.thickness,
                               shape.cos_max_limit, shape.cos_min_limit);
      }
      count = s + 1;
    }

    /**
     * @brief Destroys every committed ring and empties the pool.
     * @details The device ScalarFn member is not trivially destructible.
     */
    void release() {
      for (int s = 0; s < count; ++s) {
        (*storage)[s].ring.~DistortedRing();
#if HS_ENABLE_TEST_HOOKS
        ++destroyed;
#endif
      }
      count = 0;
    }

    /** @brief Committed rings this frame. */
    int size() const { return count; }
    /** @brief Slot of @p ring, or CULLED. */
    int8_t slot_of(int ring) const { return slot_by_ring[ring]; }
    /** @brief The ring-to-slot map, SLOTS entries. */
    const int8_t *slot_map() const { return slot_by_ring; }
    /** @brief The committed rings. */
    ShapeView shapes() const { return {storage, this}; }
    /** @brief Hue row of slot @p s. */
    __attribute__((always_inline)) const Pixel *hue_row(int s) const {
      return hue_rows + s * ROW;
    }
    /** @brief Bake column count of slot @p s. */
    __attribute__((always_inline)) float lut_columns(int s) const {
      return columns[s];
    }
    /** @brief Fragment alpha of slot @p s. */
    __attribute__((always_inline)) float frag_alpha(int s) const {
      return alphas[s];
    }

  public:
    SweepCache *cache = nullptr; ///< Current frame's alternate shifts.
    int ring_by_slot[SLOTS]{};   ///< Authored ring index of each slot.
    /** @brief Replaces one slot's fragment alpha. */
    void set_alpha(int slot, float alpha) { alphas[slot] = alpha; }
    /** @return Writable shift row of a committed slot. */
    float *shift_row(int slot) { return shift_rows + slot * ROW; }
    /** @return Writable hue row of a committed slot. */
    Pixel *mutable_hue_row(int slot) { return hue_rows + slot * ROW; }

  private:
    float *shift_rows = nullptr;     /**< SLOTS shift rows. */
    Pixel *hue_rows = nullptr;       /**< SLOTS hue rows. */
    float *alphas = nullptr;         /**< Per-slot fragment alpha. */
    float *columns = nullptr;        /**< Per-slot bake column count. */
    int8_t *slot_by_ring = nullptr;  /**< Ring index to slot, or CULLED. */
    ShapeStorage *storage = nullptr; /**< Placement storage for the rings. */
    int count = 0;                   /**< Committed rings. */
#if HS_ENABLE_TEST_HOOKS
  public:
    int destroyed = 0; /**< Rings release() has destroyed. */
#endif
  };
  SweepCache *active_cache = nullptr; ///< Active swept draw; not retained.
  RingPool rings;                     /**< Rings drawn this frame. */

  /**
   * @brief Shades a stack fragment from its slot's hue row.
   * @details v2 is the stroke coverage the scan applies again on plot, so the
   * ring edge ramps as coverage squared. The stack hands v0 over in [0, 1).
   */
  struct RingShader {
    const RingPool *pool; ///< The frame's drawn rings.

    /**
     * @brief Shades one fragment.
     * @param s Slot of the ring lighting the fragment.
     * @param f Fragment to shade; v0 is the ring azimuth in [0, 1).
     */
    __attribute__((always_inline)) void operator()(int s, const math::Vector &,
                                                   Fragment &f) const {
      const Pixel *hue = pool->hue_row(s);
      float x = f.v0 * pool->lut_columns(s);
      int j = static_cast<int>(x);
      f.color = Color4(
          hue[j].lerp16(hue[j + 1], frac_to_q16(math::quintic_kernel(x - j))),
          pool->frag_alpha(s) * f.v2);
    }
  };
  using CandidateTable = Scan::DistortedRingStack::CandidateTable<W, H>;
  CandidateTable *candidates =
      nullptr; /**< Fused scan's per-frame ring candidate map. */
  float *octave1 =
      nullptr; /**< W + 1 envelope-octave samples on its knot grid. */
  float *octave2 =
      nullptr; /**< W + 1 detail-octave samples on its knot grid. */
  float *alternate_num = nullptr; ///< Alternate-window ball numerator per knot.
  float *alternate_den =
      nullptr; ///< Alternate-window ball denominator per knot.
  uint8_t *knot_visible =
      nullptr; /**< W + 1 flags: knot lies in a visible bake chunk. */
  math::Vector *knot_pos =
      nullptr; /**< W + 1 knot positions of the ring being ball-baked. */
  float *ball_cv =
      nullptr; /**< MAX_BALLS ball-center components along the stack axis. */
  float *ball_rho =
      nullptr; /**< MAX_BALLS ball-center distances from the stack axis. */
  float *ball_azimuth =
      nullptr; /**< MAX_BALLS ball-center azimuths in the ring frame. */
  uint8_t *knot_near =
      nullptr; /**< W + 1 flags: a visible knot lies within spline reach. */
  float *chunk_cos =
      nullptr; /**< cos of each bake chunk's mid-azimuth; baked once at init. */
  float *chunk_sin =
      nullptr; /**< sin of each bake chunk's mid-azimuth; baked once at init. */
  /**
   * @brief One ring's hue-rotation table: HUE_TABLE_SIZE cells over a bound
   * hue-turn domain, baked up front or knot by knot on demand.
   * @details bind() starts each ring and clears knot validity, so no sample
   * reads a knot baked for an earlier ring.
   */
  struct HueTable {
    static constexpr int KNOTS = HUE_TABLE_SIZE + 1; /**< Knots per table. */
    /** @brief Persistent bytes of the knot storage, alignment included. */
    static constexpr size_t BYTES = KNOTS * sizeof(Pixel) + alignof(Pixel);
    Pixel *knots = nullptr; /**< KNOTS hue-rotated colors of the bound ring. */
    uint64_t valid[(HUE_TABLE_SIZE + 64) / 64] =
        {}; /**< Knots sample_lazy() has baked since the last bind(). */
    const HueRotateBase *base = nullptr; /**< Bound ring color. */
    float domain = 0.0f; /**< Hue-turn interval the knots span. */
    bool cyclic = false; /**< Amounts past the domain wrap. */

    /**
     * @brief Binds the table to one ring and clears knot validity.
     * @param ring_base Ring color's precomputed OKLab base; must outlive
     * the ring's samples.
     * @param hue_domain Hue-turn interval the knots span.
     * @param wraps Whether amounts past the domain wrap.
     */
    __attribute__((always_inline)) void bind(const HueRotateBase &ring_base,
                                             float hue_domain, bool wraps) {
      base = &ring_base;
      domain = hue_domain;
      cyclic = wraps;
      for (uint64_t &word : valid)
        word = 0;
    }

    /**
     * @brief Counts the knots a sample at or below @p max_amount can read.
     * @param max_amount Largest hue offset (turns) any sample asks for.
     * @return Knots 0..count - 1 are the only ones sample() reads.
     * @details A wrapping lookup that can reach a full turn, or a domain that
     * is not positive, reads anywhere.
     */
    int cells(float max_amount) const {
      if (!(domain > 0.0f) || (cyclic && max_amount >= domain))
        return KNOTS;
      const float x =
          hs::clamp(max_amount / domain, 0.0f, 1.0f) * HUE_TABLE_SIZE;
      return std::min(static_cast<int>(x) + 2, KNOTS);
    }

    /**
     * @brief Bakes the first @p count knots for the bound ring.
     * @param count Knots to bake, at most KNOTS.
     */
    HS_HOT_FLASH_MEMBER void prepare(int count = KNOTS) {
      for (int i = 0; i < count; ++i)
        knots[i] =
            hue_rotate(*base, domain * (static_cast<float>(i) / HUE_TABLE_SIZE))
                .color;
    }

    /** @brief Samples knots already baked by prepare(). */
    HS_O3_FN Pixel sample(float amount) const {
      return sample_with(amount, [](int) {});
    }

    /** @brief Samples the table, baking each knot it reads on first use. */
    HS_O3_FN Pixel sample_lazy(float amount) {
      return sample_with(amount, [&](int index) HS_HOT_INLINE {
        const uint64_t bit = uint64_t{1} << (index & 63);
        uint64_t &word = valid[index >> 6];
        if (!(word & bit)) {
          knots[index] = hue_rotate(*base, domain * (static_cast<float>(index) /
                                                     HUE_TABLE_SIZE))
                             .color;
          word |= bit;
        }
      });
    }

  private:
    /**
     * @brief Interpolates the table at a hue offset.
     * @param amount Hue offset in turns.
     * @param ensure Called with every knot index read, before the read.
     */
    template <typename Ensure>
    Pixel sample_with(float amount, Ensure ensure) const {
      assert(base && "DisplacementField::HueTable sampled before bind()");
      float t = amount / domain;
      t = cyclic ? math::wrap_t(t) : hs::clamp(t, 0.0f, 1.0f);
      float x = t * HUE_TABLE_SIZE;
      if (x >= HUE_TABLE_SIZE) {
        ensure(HUE_TABLE_SIZE);
        return knots[HUE_TABLE_SIZE];
      }
      int i = static_cast<int>(x);
      ensure(i);
      ensure(i + 1);
      return knots[i].lerp16(knots[i + 1], frac_to_q16(x - i));
    }
  };
  HueTable hue_table; /**< Hue table of the ring being baked. */
  float *ball_colat =
      nullptr; /**< MAX_BALLS active-ball center colatitudes about the stack axis (radians), rebuilt per frame. */
  float *ball_reach =
      nullptr; /**< MAX_BALLS active-ball support extents (radians): both the reach prefilter bound and the per-ring band bound. */
  float *ball_scale =
      nullptr; /**< MAX_BALLS active-ball LUT feature scales (2/radius). */
  const Animation::BumpParams **ball_params =
      nullptr; /**< Active-ball params validated and cached once per frame. */
  int *ball_local =
      nullptr; /**< MAX_BALLS scratch: active indices of the balls that can reach the current ring. */
#if HS_ENABLE_TEST_ORACLES
  bool force_exact_hue = false;
  bool force_column_reference = false; ///< Bake each column at its exact time.
  bool force_window_bakes = false;     ///< Two-window field-bake oracle.
#endif
#if HS_ENABLE_TEST_HOOKS
  int primary_knots =
      0; ///< Primary visible field knots in the last swept frame.
  int duplicate_knots = 0;         ///< Distinct retained alternate shifts.
  int alternate_noise_samples = 0; ///< Additional scalar octave evaluations.
  int primary_bakes = 0; ///< Non-flat rings baked in the last swept frame.
  int stack_flushes = 0; ///< Stack scans in the last swept frame.
  int hue_table_uses = 0;
#endif

  /** Ring-to-ball prefilter pad absorbing fast_acos and tangent-recurrence
   *  rounding (radians); a ball excluded despite the pad fields exactly 0
   *  everywhere on the ring. */
  static constexpr float BALL_TOUCH_EPS = 1e-3f;

public:
  /** @brief Slider-backed parameters. */
  struct Params {
    float alpha = 0.3f; /**< Overall ring opacity multiplier in [0, 1]. */
    int num_rings =
        std::min(48, RING_SLOTS); /**< Number of evenly spaced rings. */
    float thickness = 0.03f;      /**< Stroke half-width (radians). */
    float ball_amp =
        0.1f; /**< Ball drape strength; scaled by BALL_DRAPE_PER_AMPLITUDE into the drape gain. */
    float noise_amp =
        0.2f; /**< Peak polar displacement (radians) of the noise phase. */
    float scale1 =
        1.5f; /**< Spatial frequency of the envelope octave; its zero regions leave rings undisturbed. */
    float scale2 =
        3.0f; /**< Spatial frequency of the detail octave, scaled by the envelope octave. */
    float hue_scale =
        2.0f; /**< Hue rotation (turns) per radian of displacement magnitude. */
    float flow_speed = 0.03f; /**< Noise-field time advance per frame. */
    float ball_min = 0.15f;   /**< Smallest ball footprint (radians). */
    float ball_max = 0.3f;    /**< Largest ball footprint (radians). */
    float ball_rate =
        20.0f; /**< Ball spawns per BALL_RATE_FPS frames (jittered ±50%). */
    float ball_speed_min =
        0.45f; /**< Slowest fall (pole-to-pole traversals per BALL_RATE_FPS frames). */
    float ball_speed_max =
        0.85f; /**< Fastest fall (pole-to-pole traversals per BALL_RATE_FPS frames). */
  };

  /** Authored parameter defaults for Noise and Balls. */
  static constexpr std::array<Params, 2> PRESETS{Params{}, Params{}};

private:
  Params params = PRESETS[0];

  static_assert(Params{}.thickness >= 0.4f * THICKNESS_PX &&
                    Params{}.thickness <= 6.0f * THICKNESS_PX,
                "DisplacementField Thickness default falls outside its "
                "W-scaled slider range at this build resolution");

  /** @brief Persistent bytes of the knot bake scratch, alignment included. */
  static constexpr size_t KNOT_SCRATCH_BYTES =
      (W + 1) *
          (4 * sizeof(float) + 2 * sizeof(uint8_t) + sizeof(math::Vector)) +
      4 * alignof(float) + 2 * alignof(uint8_t) + alignof(math::Vector);
  /** @brief Persistent bytes of the per-frame ball cache. */
  static constexpr size_t BALL_CACHE_BYTES =
      MAX_BALLS * (6 * sizeof(float) + sizeof(int) +
                   sizeof(const Animation::BumpParams *)) +
      6 * alignof(float) + alignof(int) +
      alignof(const Animation::BumpParams *);
  /** @brief Persistent bytes of the last swept draw's ring colors. */
  static constexpr size_t COLOR_HISTORY_BYTES =
      RING_SLOTS * (sizeof(Color4) + sizeof(uint8_t)) + alignof(Color4) +
      alignof(uint8_t);
  /** @brief Persistent bytes of the chunk mid-azimuth trig tables. */
  static constexpr size_t CHUNK_TRIG_BYTES =
      2 * BAKE_CHUNKS * sizeof(float) + 2 * alignof(float);
  /** @brief Persistent bytes of the fused scan's candidate table. */
  static constexpr size_t CANDIDATE_BYTES =
      sizeof(CandidateTable) + alignof(CandidateTable);
  /** @brief Persistent bytes of the two displacement-field pools. */
  static constexpr size_t FIELD_BYTES =
      MAX_BALLS * (sizeof(typename decltype(balls)::Entity) + sizeof(int)) +
      alignof(typename decltype(balls)::Entity) + alignof(int) +
      sizeof(typename decltype(noise_field)::Entity) + sizeof(int) +
      alignof(typename decltype(noise_field)::Entity) + alignof(int);
  /** @brief Every persistent allocation init() makes. */
  static constexpr size_t FOOTPRINT_BYTES =
      RingPool::BYTES + HueTable::BYTES + KNOT_SCRATCH_BYTES +
      BALL_CACHE_BYTES + COLOR_HISTORY_BYTES + CHUNK_TRIG_BYTES +
      CANDIDATE_BYTES + FIELD_BYTES;
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "DisplacementField persistent footprint exceeds the default "
                "partition; retune RING_SLOTS/MAX_BALLS or carve arenas");
};
