/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/render/sdf/ring_knots.h"

/**
 * @file DisplacementField.h
 * @brief Stack of soft-stroked rings displaced by a noise and falling-ball
 *        displacement-field stack.
 */

#include "core/color/hue_rotation_row.h"
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
  HS_FLASH_MEMBER DisplacementField()
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
    alternate_a = persistent_arena.allocate_n<float>(W + 1);
    alternate_b = persistent_arena.allocate_n<float>(W + 1);
    knot_visible = persistent_arena.allocate_n<uint8_t>(W + 1);
    knot_pos = persistent_arena.allocate_n<math::Vector>(W + 1);
    knot_near = persistent_arena.allocate_n<uint8_t>(W + 1);
    drawn_colors = persistent_arena.allocate_n<Color4>(RING_SLOTS);
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

    register_described_params<&DisplacementField::parameter_fields>(params);

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
  HS_FLASH_MEMBER void preset_changed(const PresetChange &change) override {
    hs::log("Preset: %u/%u", static_cast<unsigned>(change.to + 1),
            static_cast<unsigned>(authored_preset_count()));
  }
#endif

private:
  HS_NOINLINE_NOCLONE static int
  noise_lut_samples(float feature_scale, float sin_theta, bool octave_path) {
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
   * @param visible Visible-chunk mask, one bit per bake chunk.
   * @param hue_extent Hue turns the displacement range covers.
   * @return Table sampling mode and domain.
   * @details Culled chunks bake no columns, so the precompute decision
   * amortizes the table over the visible chunks' columns rather than lut_n.
   */
  HS_HOT_INLINE HueMode select_hue_mode(int lut_n, uint32_t visible,
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
            HS_HOT_INLINE { visible_samples += end - begin; },
        [](int, int) HS_HOT_INLINE {});
    mode.precompute = visible_samples > 2 * HUE_TABLE_SIZE;
    mode.cyclic = hue_extent > 1.0f;
    mode.domain = mode.cyclic ? 1.0f : hue_extent;
#if HS_ENABLE_TEST_HOOKS
    ++hue_table_uses;
#endif
    return mode;
  }

  static constexpr int ALTERNATE_CAPACITY = 3840; ///< Packed alternate shifts.
  using Alternates = Scan::KnotAlternates<W>;
  /** @brief Second-window hue and opacity state. */
  struct SecondWindowShading {
    HueMode *modes; ///< Per-ring hue sampling mode.
    Color4 *colors; ///< Second-window colors.
    float *alphas;  ///< Second-window fragment alpha.
  };
  Alternates alternates; ///< Current frame's alternate knots and counters.
  SecondWindowShading shading; ///< Current frame's late shading rows.

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
    Alternates *cache = nullptr; ///< Single-pass alternate storage.
    float middle = 0.5f;         ///< Middle sweep fraction of the drawn window.

    /** @brief Whether alternate samples use the cache. */
    HS_HOT_INLINE bool uses_cache() const {
#if HS_ENABLE_TEST_ORACLES
      return cache != nullptr;
#else
      return true;
#endif
    }
  };

  static constexpr float SWEEP_QUARTER = Scan::SweptRingStack::QUARTER;
  static constexpr int SWEEP_WINDOWS = Scan::SweptRingStack::WINDOWS;
  using ScanSweep = Scan::SweptRingStack;
  /** @brief Returns a knot's sweep fraction, or the reference column's time.
   * @param p Knot position at the current endpoint.
   * @param sweep Frame endpoints and window midpoint. */
  HS_HOT_INLINE float sample_fraction(const math::Vector &p,
                                      const Sweep &sweep) const {
#if HS_ENABLE_TEST_ORACLES
    if (force_column_reference)
      return sweep.middle;
#endif
    return ScanSweep::sweep_fraction(p, sweep.middle);
  }

  /** @brief ball_field() with each ball at a frame fraction. */
  HS_O3_FN float ball_field_at(const math::Vector &p, const int *ks, int n,
                               float theta, float t) const {
    DominantFieldAccumulator accumulator;
    for (int j = 0; j < n; ++j) {
      const int k = ks[j];
      accumulator.add(Animation::bump_field_at(p, *ball_params[k],
                                               ball_paths[k], theta, t));
    }
    return accumulator.value();
  }

  /**
   * @brief Caches one ball's bounds in the ring frame.
   * @param b Ball index.
   * @param basis Ring frame.
   * @param centre Ball bound centre.
   * @param reach Angular bound radius.
   * @param radius Ball radius.
   */
  HS_HOT_INLINE void cache_ball_frame(int b, const math::Basis &basis,
                                      const math::Vector &centre, float reach,
                                      float radius) {
    const float cv = math::dot(basis.v, centre);
    ball_colat[b] = math::fast_acos(hs::clamp(cv, -1.0f, 1.0f));
    ball_reach[b] = reach;
    ball_scale[b] = 2.0f / radius;
    const float cu = math::dot(basis.u, centre);
    const float cw = math::dot(basis.w, centre);
    ball_cv[b] = cv;
    ball_rho[b] = sqrtf(cu * cu + cw * cw);
    const float a0 = atan2f(cw, cu);
    ball_azimuth[b] = a0 < 0.0f ? a0 + 2.0f * math::PI_F : a0;
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
    const math::Quaternion &current = orientation.get();
    const math::Quaternion previous =
        draw_history.start(DrawEndpoint{.pose = current}).pose;
    ScratchScope column_guard(scratch_arena_a);
    const Scan::ColumnBasis column_basis =
        Scan::prepare_column_basis<W, H>(scratch_arena_a, previous, current);
    Scan::ColumnBasisScope column_scope(canvas, column_basis);

    const int n_balls = preset_index == BALLS_PRESET && !seed_balls_pending
                            ? balls.active_count()
                            : 0;
    Sweep sweep;
    const bool swept =
        prepare_sweep(previous, current, opacity, n_balls, sweep);
    if (swept) {
      draw_swept(canvas, basis, column_basis, n_balls, sweep);
      return;
    }
    drawn_colors_valid = false;
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      ball_params[b] = &sp;
      cache_ball_frame(b, basis, sp.center, sp.field_bound(), sp.radius);
    }
    draw_pass<false>(canvas, basis, column_basis, opacity, n_balls, nullptr);
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
    const bool primed = draw_history.initialized;
    const DrawEndpoint end{current, time, amplitude, opacity};
    const DrawEndpoint start = draw_history.start(end);
    sweep.noise_time[0] = start.noise_time;
    sweep.noise_time[1] = time;
    sweep.noise_amplitude[0] = start.noise_amplitude;
    sweep.noise_amplitude[1] = amplitude;
    sweep.opacity[0] = start.opacity;
    sweep.opacity[1] = opacity;
    draw_history.capture(end);
    if (Spacetime::flat_surface() || !primed)
      return false;
    const bool noise = noise_live && fmaxf(std::fabs(sweep.noise_amplitude[0]),
                                           std::fabs(amplitude)) >
                                         Animation::FIELD_AMPLITUDE_CUT;
    if (!noise && n_balls == 0)
      return false;
    sweep.turn = Spacetime::rotation_vector(previous, current);
    return true;
  }

  /**
   * @brief Bakes column-time knots once and scans the stack in two windows.
   * @param canvas Render target.
   * @param basis Stack frame at the current endpoint.
   * @param column_basis The frame's column basis.
   * @param n_balls Active balls drawn this frame.
   * @param sweep The frame's endpoints.
   * @details Boundary knots retain an alternate for the other arm pass.
   * Colors and opacity take each window's middle time. A ring whose boundary
   * knots would overfill the alternate buffer first scans the rings before
   * it, which keeps every pixel's blend order.
   */
  HS_HOT_FLASH_MEMBER void draw_swept(Canvas &canvas, const math::Basis &basis,
                                      const Scan::ColumnBasis &column_basis,
                                      int n_balls, Sweep &sweep) {
#if HS_ENABLE_TEST_HOOKS
    primary_knots = alternate_noise_samples = primary_bakes = 0;
    alternates.duplicate_knots = alternates.stack_flushes = 0;
#endif
    ScratchScope sweep_guard(scratch_arena_a);
    ball_paths = scratch_arena_a.allocate_n<Animation::BumpPath>(MAX_BALLS);
    knot_times = scratch_arena_a.allocate_n<float>(W + 1);
    ring_colors = scratch_arena_a.allocate_n<Color4>(RING_SLOTS);
    for (int i = 0; i < params.num_rings; ++i)
      ring_colors[i] = ring_color(i, params.num_rings);
    if (drawn_ring_count != params.num_rings)
      drawn_colors_valid = false;
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      ball_params[b] = &sp;
      Animation::BumpCover cover;
      if (!Animation::make_bump_path(sp, sweep.turn, basis.v, ball_paths[b],
                                     cover)) {
        ball_colat[b] = 0.0f;
        ball_reach[b] = -1.0f;
        continue;
      }
      cache_ball_frame(b, basis, cover.centre, cover.reach, sp.radius);
    }
    ScratchScope cache_guard(scratch_arena_b);
    Alternates &cache = alternates;
    cache.allocate(scratch_arena_b, scratch_arena_a, RING_SLOTS,
                   ALTERNATE_CAPACITY, [&](Arena &arena) {
                     shading = {arena.make_n<HueMode>(RING_SLOTS),
                                arena.allocate_n<Color4>(RING_SLOTS),
                                arena.allocate_n<float>(RING_SLOTS)};
                   });
    sweep.cache = &cache;
    active_cache = &cache;
#if HS_ENABLE_TEST_ORACLES
    if (force_window_bakes || force_column_reference)
      sweep.cache = active_cache = nullptr;
#endif
    if (sweep.uses_cache())
      draw_pass<true>(canvas, basis, column_basis, sweep.opacity[1], n_balls,
                      &sweep);
#if HS_ENABLE_TEST_ORACLES
    if (!sweep.uses_cache()) {
      const int windows = force_column_reference ? W / 2 : SWEEP_WINDOWS;
      for (int k = 0; k < windows; ++k) {
        sweep.middle = force_column_reference ? static_cast<float>(k) / windows
                                              : (k + 0.5f) / windows;
        TimeWindowScope window(canvas, k * (W / 2) / windows,
                               (k + 1) * (W / 2) / windows);
        draw_pass<true>(canvas, basis, column_basis, sweep.opacity[1], n_balls,
                        &sweep);
      }
    }
#endif
    active_cache = nullptr;
    sweep.cache = nullptr;
    for (int i = 0; i < params.num_rings; ++i)
      drawn_colors[i] = ring_colors[i];
    drawn_colors_valid = true;
    drawn_ring_count = params.num_rings;
    ball_paths = nullptr;
    knot_times = nullptr;
    ring_colors = nullptr;
  }

  /**
   * @brief A ring's palette color at the current endpoint.
   * @param i Ring index.
   * @param n_rings Stack size.
   * @return The ring's color.
   */
  HS_FLASH_MEMBER Color4 ring_color(int i, int n_rings) const {
    return palette.get(math::wrap_t((i + 0.5f) / n_rings + color_spin));
  }

  /**
   * @brief Bakes and scans the rings of one pass over the stack.
   * @tparam SWEPT Whether knots are baked at their sweep times.
   * @param canvas Render target.
   * @param basis Stack frame at the current endpoint.
   * @param column_basis The frame's column basis.
   * @param opacity Sprite opacity; a swept pass takes the window's own.
   * @param n_balls Active balls drawn this frame, with their caches filled.
   * @param sweep The frame's endpoints; null unless SWEPT.
   */
  template <bool SWEPT>
  HS_HOT_INLINE void draw_pass(Canvas &canvas, const math::Basis &basis,
                               const Scan::ColumnBasis &column_basis,
                               float opacity, int n_balls, const Sweep *sweep) {
    const int n_rings = params.num_rings;
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the rasterizer's soft stroke cross-section past
    // the clip's margin and the stack's rotation over the frame.
    const float pad = 3.0f * math::PI_F / H + column_basis.sweep_arc;
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

      Color4 color;
      float ring_opacity = opacity;
      if constexpr (SWEPT) {
        const Color4 now = ring_colors[i];
        const Color4 before = drawn_colors_valid ? drawn_colors[i] : now;
        const float middle =
            sweep->uses_cache() ? SWEEP_QUARTER : sweep->middle;
        color = before.lerp(now, middle);
        if (sweep->uses_cache()) {
          shading.colors[i] = before.lerp(now, 1.0f - SWEEP_QUARTER);
          shading.alphas[i] =
              shading.colors[i].alpha *
              (sweep->opacity[0] + (sweep->opacity[1] - sweep->opacity[0]) *
                                       (1.0f - SWEEP_QUARTER)) *
              params.alpha;
          sweep->cache->begin_ring(i);
        }
        ring_opacity = sweep->opacity[0] +
                       (sweep->opacity[1] - sweep->opacity[0]) * middle;
      } else {
        color = ring_color(i, n_rings);
      }
      HueRotateBase hue_base = make_hue_rotate_base(color);

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
                    std::fabs(octave_amplitude) > Animation::FIELD_AMPLITUDE_CUT
                ? &noise_field.active_params(0)
                : nullptr;
        lut_n = noise_lut_samples(feature_scale, sin_t, octaves != nullptr);

        uint32_t visible = CHUNK_MASK;
        if constexpr (SWEPT) {
          const float band_r = band + noise_bound + params.thickness + pad;
          const Spacetime::TimeWindow &window = canvas.time_window();
          {
            HS_PROFILE(df_chunk_cull);
            visible = Plot::chunk_mask_where<BAKE_CHUNKS>(
                [&](const math::Vector &mid, float reach, float sin_reach)
                    HS_HOT_INLINE {
                      return (sweep->uses_cache() ||
                              Plot::cap_may_touch_window<W>(window, mid, reach,
                                                            sin_reach)) &&
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
          if (sweep->uses_cache()) {
            Alternates &cache = *sweep->cache;
            // A ring retains at most one alternate per knot.
            if (!cache.fits(lut_n)) {
              scan_swept_stack(canvas, cache);
              rings.begin_frame(&cache);
              cache.begin_ring(i);
              const typename RingPool::Pending fresh = rings.next();
              slut = fresh.shift_row;
              hlut = fresh.hue_row;
            }
            ScanSweep::mark_arm_boundary<W>(
                basis, {cos_t, sin_t, cos_d, sin_d}, theta, lut_n,
                band + noise_bound, params.thickness, column_basis.sweep_arc,
                2 * std::max(OCTAVE1_STRIDE, OCTAVE2_STRIDE) - 1,
                cache.boundary);
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
            [&](int begin, int end) HS_HOT_INLINE {
              for (int x = begin; x < end; ++x)
                max_shift = fmaxf(max_shift, std::fabs(slut[x]));
            },
            [&](int begin, int end) HS_HOT_INLINE {
              for (int x = begin; x < end; ++x)
                slut[x] = 0.0f;
            });

        const HueMode mode = select_hue_mode(
            lut_n, visible, (band + noise_bound) * params.hue_scale);
        if constexpr (SWEPT)
          if (sweep->uses_cache())
            shading.modes[i] = mode;
        if (mode.use_table)
          hue_table.bind(hue_base, mode.domain, mode.cyclic);
        if (mode.precompute) {
          HS_PROFILE(df_hue_table_prep);
          hue_table.prepare(
              hue_table.knots_through(max_shift * params.hue_scale));
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
            [&](int begin, int end) HS_HOT_INLINE {
              for (int x = begin; x < end; ++x)
                hlut[x] = hue_for_shift(slut[x]);
            },
            [](int, int) HS_HOT_INLINE {});
        slut[lut_n] = slut[0];
        hlut[lut_n] = hlut[0];
      }

      if constexpr (SWEPT) {
        if (sweep->uses_cache()) {
          Alternates &cache = *sweep->cache;
          cache.finish_ring(i);
        }
      }
      rings.commit(i, color.alpha * ring_opacity * params.alpha, lut_n, basis,
                   radius, params.thickness, slut, lut_n, 0.0f, nullptr);
    }

    if (rings.size() == 0)
      return;

    if constexpr (SWEPT) {
      if (sweep->uses_cache()) {
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
  HS_HOT_FLASH_MEMBER void scan_swept_stack(Canvas &canvas, Alternates &cache) {
    RingShader ring_shader{&rings};
    struct Reshade {
      DisplacementField &effect;
      Pixel *hues = nullptr;
      HueMode mode{};
      HueRotateBase base{};
      Pixel zero_hue{};
      bool zero_ready = false;
      int n = 0;
      int ring = 0;
      int slot = 0;

      HS_HOT_INLINE void begin(int i, int s) {
        ring = i;
        slot = s;
        hues = effect.rings.mutable_hue_row(slot);
        n = static_cast<int>(effect.rings.lut_columns(slot));
        mode = effect.shading.modes[i];
        base = make_hue_rotate_base(effect.shading.colors[i]);
        if (mode.use_table)
          effect.hue_table.bind(base, mode.domain, mode.cyclic);
        zero_ready = false;
      }
      HS_HOT_INLINE void operator()(int k, float shift) {
        const float amount = std::fabs(shift) * effect.params.hue_scale;
        if (amount == 0.0f) {
          if (!zero_ready) {
            zero_hue = mode.use_table ? effect.hue_table.sample_lazy(amount)
                                      : hue_rotate(base, amount).color;
            zero_ready = true;
          }
          hues[k] = zero_hue;
        } else {
          hues[k] = mode.use_table ? effect.hue_table.sample_lazy(amount)
                                   : hue_rotate(base, amount).color;
        }
      }
      HS_HOT_INLINE void finish() {
        hues[n] = hues[0];
        effect.rings.set_alpha(slot, effect.shading.alphas[ring]);
      }
    };
    HS_PROFILE(df_fused_scan);
    ScanSweep::scan_two_windows<W, H>(filters, canvas, params.num_rings, rings,
                                      cache, *candidates, ring_shader,
                                      Reshade{*this});
  }

  /**
   * @brief Sets knot_visible[x] to whether knot x lies in a visible chunk.
   * @param lut_n Knot count.
   * @param visible Visible-chunk mask, one bit per bake chunk.
   */
  HS_HOT_INLINE void mark_visible_knots(int lut_n, uint32_t visible) {
    int x = 0;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const int x_end = Plot::chunk_end<BAKE_CHUNKS>(c, lut_n);
      const uint8_t v = static_cast<uint8_t>((visible >> c) & 1u);
      for (; x < x_end; ++x)
        knot_visible[x] = v;
    }
  }

  /**
   * @brief Bakes one ring's centerline shifts with each noise octave sampled on
   *        its own knot grid.
   * @param np The noise field's single active entity.
   * @param basis Ring frame; basis.v is the stack axis.
   * @param step Trig of theta and of one knot cell's azimuth step.
   * @param theta Ring colatitude.
   * @param lut_n Knot count, a multiple of OCTAVE_GRID.
   * @param visible Visible-chunk mask, one bit per bake chunk.
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
                     const math::Basis &basis, SDF::KnotStep step, float theta,
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
    SDF::KnotWalk knots{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
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
              Spacetime::rotate_by(p, sweep->turn * (t - 1.0f));
          const float time = sweep->noise_time[0] +
                             (sweep->noise_time[1] - sweep->noise_time[0]) * t;
          if (g1)
            octave1[x] = np.octave1(q, time);
          if (g2)
            octave2[x] = np.octave2(q, time);
          if (sweep->uses_cache() && (g1 || g2)) {
            alternate_a[x] = g1 ? octave1[x] : 0.0f;
            alternate_b[x] = g2 ? octave2[x] : 0.0f;
            if (sweep->cache->boundary[x] && ScanSweep::retains_endpoint(t)) {
#if HS_ENABLE_TEST_HOOKS
              alternate_noise_samples +=
                  static_cast<int>(g1) + static_cast<int>(g2);
#endif
              const float other = ScanSweep::other_endpoint(t);
              const math::Vector aq =
                  Spacetime::rotate_by(p, sweep->turn * (other - 1.0f));
              const float at =
                  sweep->noise_time[0] +
                  (sweep->noise_time[1] - sweep->noise_time[0]) * other;
              if (g1)
                alternate_a[x] = np.octave1(aq, at);
              if (g2)
                alternate_b[x] = np.octave2(aq, at);
            }
          }
          if (vis)
            slut[x] = ball_field_at(p, ball_local, n_local, theta, t);
        } else {
          if (g1)
            octave1[x] = np.octave1(p, np.time);
          if (g2)
            octave2[x] = np.octave2(p, np.time);
          if (vis)
            slut[x] = ball_field(p, ball_local, n_local, theta);
        }
      }
    }

    bool late = false;
    auto sample = [&](const float *oct, const float *alternate, int d,
                      int k) HS_HOT_INLINE {
      auto value = [&](int index) HS_HOT_INLINE {
        if constexpr (SWEPT) {
          if (sweep->uses_cache() &&
              ((late && knot_times[index] < SWEEP_QUARTER) ||
               (!late && knot_times[index] >= 1.0f - SWEEP_QUARTER)))
            return alternate[index];
        }
        return oct[index];
      };
      return SDF::catmull_rom_strided(value, d, k, lut_n);
    };
    SDF::KnotWalk final_walk{basis, step.cos_t, step.sin_t, step.cos_d,
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
      const float primary = slut[x] + amplitude *
                                          sample(octave1, alternate_a, D1, x) *
                                          sample(octave2, alternate_b, D2, x);
      if constexpr (SWEPT) {
        if (sweep->uses_cache() && sweep->cache->boundary[x]) {
          late = !late;
          const float t = knot_times[x];
          const float other =
              ScanSweep::retains_endpoint(t) ? ScanSweep::other_endpoint(t) : t;
          const float amp =
              sweep->noise_amplitude[0] +
              (sweep->noise_amplitude[1] - sweep->noise_amplitude[0]) * other;
          const float ball =
              n_local ? ball_field_at(final_walk.position(), ball_local,
                                      n_local, theta, other)
                      : 0.0f;
          const float alternate =
              ball + amp * sample(octave1, alternate_a, D1, x) *
                         sample(octave2, alternate_b, D2, x);
          sweep->cache->retain(x, primary, alternate, t, slut);
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
   * @param visible Visible-chunk mask, one bit per bake chunk.
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
  bake_ball_spans(const math::Basis &basis, SDF::KnotStep step, float theta,
                  int lut_n, uint32_t visible, int n_local, float *slut,
                  const Sweep *sweep = nullptr) {
    // Covers the recurrence's drift from the analytic knot azimuth and the
    // float slop in the cap dot product.
    constexpr float COS_MARGIN = 2e-4f;
    float *num = octave1;
    float *den = octave2;
    mark_visible_knots(lut_n, visible);
    SDF::KnotWalk knots{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
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
          if (sweep->uses_cache())
            alternate_a[x] = alternate_b[x] = 0.0f;
        }
      }
    }

    for (int j = 0; j < n_local; ++j) {
      const int k = ball_local[j];
      const Animation::BumpParams &bp = *ball_params[k];
      float cos_radius = bp.cos_radius;
      if constexpr (SWEPT)
        cos_radius = ball_paths[k].cos_cover;
      const auto arc =
          Plot::cap_ring_arc(step.cos_t, step.sin_t, ball_cv[k], ball_rho[k],
                             ball_azimuth[k], cos_radius, lut_n, COS_MARGIN);
      if (arc.empty)
        continue;
      const int first = arc.first;
      const int count = arc.count;
      const float y = theta - ball_colat[k];
      for (int i = 0, xi = first; i < count; ++i) {
        if (knot_visible[xi]) {
          float f;
          if constexpr (SWEPT)
            f = Animation::bump_field_at(knot_pos[xi], bp, ball_paths[k], theta,
                                         knot_times[xi]);
          else
            f = bump_field_with_y(knot_pos[xi], bp, y);
          DominantFieldAccumulator::accumulate(num[xi], den[xi], f);
          if constexpr (SWEPT) {
            const float t = knot_times[xi];
            if (sweep->uses_cache() && sweep->cache->boundary[xi] &&
                ScanSweep::retains_endpoint(t)) {
              const float other = ScanSweep::other_endpoint(t);
              const float alternate = Animation::bump_field_at(
                  knot_pos[xi], bp, ball_paths[k], theta, other);
              DominantFieldAccumulator::accumulate(alternate_a[xi],
                                                   alternate_b[xi], alternate);
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
          if (sweep->uses_cache() && sweep->cache->boundary[x] &&
              ScanSweep::retains_endpoint(t)) {
            const float alternate = DominantFieldAccumulator::resolve(
                                        alternate_a[x], alternate_b[x]) +
                                    noise_field.field(knot_pos[x]);
            sweep->cache->retain(x, slut[x], alternate, t, slut);
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
        PaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Spawns one falling ball with a random meridian, footprint, and
   * speed drawn from the Speed Min/Max sliders; dropped safely if the ball pool
   * or the timeline is full.
   * @param progress Initial fall progress; zero starts at the north pole.
   * @param canvas Active frame canvas for nonzero initial progress.
   * @details Only the first pool-full drop of each ball phase is logged.
   */
  HS_FLASH_MEMBER void spawn_ball(float progress = 0.0f,
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
  HS_FLASH_MEMBER bool apply_preset(const PresetChange &change) override {
    params = PRESETS[change.to];
    if (change.origin == PresetChangeOrigin::AUTOMATIC) {
      if (change.to == NOISE_PRESET)
        enter_noise();
      else
        enter_balls();
    } else {
      timeline.cancel_owner(&phase_events);
      timeline.cancel_owner(&balls);
      draw_history.reset();
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
  HS_FLASH_MEMBER void parameter_written() override {
    if (preset_index == NOISE_PRESET && noise_hold == 0) {
      timeline.cancel_owner(&phase_events);
      schedule_gain(1.0f, false);
    }
    noise_hold = NOISE_HOLD_FRAMES;
    ball_phase_left = BALL_PHASE_FRAMES;
  }

  /**
   * @brief Enters the noise phase and fades it in from zero.
   */
  HS_FLASH_MEMBER void enter_noise() {
    master_gain = 0.0f;
    noise_hold = NOISE_HOLD_FRAMES;
    schedule_gain(1.0f, false);
  }

  /**
   * @brief Enters the ball phase, opening a fresh spawning window.
   */
  HS_FLASH_MEMBER void enter_balls() {
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
  /** @brief Ring-stack pose and shading controls at one sprite draw. */
  struct DrawEndpoint {
    math::Quaternion pose;        ///< Stack orientation used by the draw.
    float noise_time = 0.0f;      ///< Noise field phase.
    float noise_amplitude = 0.0f; ///< Noise displacement amplitude.
    float opacity = 0.0f;         ///< Sprite opacity multiplier.
  };
  Spacetime::Retained<DrawEndpoint> draw_history; ///< Last ring-stack draw.
  Color4 *drawn_colors =
      nullptr; ///< RING_SLOTS ring colors of the last swept draw.
  bool drawn_colors_valid = false; ///< The last draw was swept.
  int drawn_ring_count = 0;        ///< Stack size of the last swept draw.
  // Frame scratch of a swept draw, null outside it.
  Animation::BumpPath *ball_paths = nullptr; ///< MAX_BALLS active-ball paths.
  float *knot_times = nullptr;   ///< W + 1 frame fractions of baked knots.
  Color4 *ring_colors = nullptr; ///< RING_SLOTS current ring colors.

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

  Alternates *active_cache = nullptr; ///< Active swept draw; not retained.
  using RingPool = Scan::RingStackPool<RING_SLOTS, W + 1>;
  using RingShader = Scan::HueRowShader<RingPool>;
  RingPool rings; /**< Rings drawn this frame. */

  using CandidateTable = Scan::DistortedRingStack::CandidateTable<W, H>;
  CandidateTable *candidates =
      nullptr; /**< Fused scan's per-frame ring candidate map. */
  float *octave1 =
      nullptr; /**< W + 1 envelope-octave samples on its knot grid. */
  float *octave2 =
      nullptr; /**< W + 1 detail-octave samples on its knot grid. */
  float *alternate_a =
      nullptr; ///< Alternate-window octave 1 or ball numerator per knot.
  float *alternate_b =
      nullptr; ///< Alternate-window octave 2 or ball denominator per knot.
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
  using HueTable = HueRotationRow<HUE_TABLE_SIZE>;
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
  int alternate_noise_samples = 0; ///< Additional scalar octave evaluations.
  int primary_bakes = 0; ///< Non-flat rings baked in the last swept frame.
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

  /** @brief GUI parameter descriptions, in registration order.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto parameter_fields() {
    return std::tuple{
        Control::Field<Params, float>{.id = "alpha",
                                      .member = &Params::alpha,
                                      .name = "Alpha",
                                      .spec = {.min = 0.0f, .max = 1.0f}},
        Control::Field<Params, int>{.id = "num_rings",
                                    .member = &Params::num_rings,
                                    .name = "Rings",
                                    .spec = {.min = 1, .max = RING_SLOTS}},
        Control::Field<Params, float>{
            .id = "thickness",
            .member = &Params::thickness,
            .name = "Thickness",
            .spec = {.min = 0.4f * THICKNESS_PX, .max = 6.0f * THICKNESS_PX}},
        Control::Field<Params, float>{.id = "ball_amp",
                                      .member = &Params::ball_amp,
                                      .name = "Ball Amp",
                                      .spec = {.min = 0.0f, .max = 0.8f}},
        Control::Field<Params, float>{.id = "noise_amp",
                                      .member = &Params::noise_amp,
                                      .name = "Noise Amp",
                                      .spec = {.min = 0.0f, .max = 0.8f}},
        Control::Field<Params, float>{.id = "scale1",
                                      .member = &Params::scale1,
                                      .name = "Scale 1",
                                      .spec = {.min = 0.5f, .max = 4.0f}},
        Control::Field<Params, float>{.id = "scale2",
                                      .member = &Params::scale2,
                                      .name = "Scale 2",
                                      .spec = {.min = 0.5f, .max = 8.0f}},
        Control::Field<Params, float>{.id = "hue_scale",
                                      .member = &Params::hue_scale,
                                      .name = "Hue Rotate",
                                      .spec = {.min = 0.0f, .max = 3.0f}},
        Control::Field<Params, float>{.id = "flow_speed",
                                      .member = &Params::flow_speed,
                                      .name = "Flow Speed",
                                      .spec = {.min = 0.0f, .max = 0.15f}},
        Control::Field<Params, float>{.id = "ball_min",
                                      .member = &Params::ball_min,
                                      .name = "Ball Min",
                                      .spec = {.min = 0.05f, .max = 1.0f}},
        Control::Field<Params, float>{.id = "ball_max",
                                      .member = &Params::ball_max,
                                      .name = "Ball Max",
                                      .spec = {.min = 0.05f, .max = 1.0f}},
        Control::Field<Params, float>{.id = "ball_rate",
                                      .member = &Params::ball_rate,
                                      .name = "Ball Rate",
                                      .spec = {.min = 0.5f, .max = 50.0f}},
        Control::Field<Params, float>{.id = "ball_speed_min",
                                      .member = &Params::ball_speed_min,
                                      .name = "Speed Min",
                                      .spec = {.min = 0.1f, .max = 3.0f}},
        Control::Field<Params, float>{.id = "ball_speed_max",
                                      .member = &Params::ball_speed_max,
                                      .name = "Speed Max",
                                      .spec = {.min = 0.1f, .max = 3.0f}}};
  }

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
      RING_SLOTS * sizeof(Color4) + alignof(Color4);
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

  /** @brief Peak swept staging with the column basis live in scratch A. */
  static constexpr size_t SCRATCH_A_PEAK_BYTES = [] {
    ArenaFootprint peak;
    peak.template add<Scan::ColumnBasisTable<W, H>>(1);
    peak.template add<Animation::BumpPath>(MAX_BALLS);
    peak.template add<float>(W + 1);
    peak.template add<Color4>(RING_SLOTS);
    ArenaFootprint values;
    Alternates::add_footprint(values, peak, RING_SLOTS, ALTERNATE_CAPACITY,
                              [](ArenaFootprint &table) {
                                table.add<HueMode>(RING_SLOTS);
                                table.add<Color4>(RING_SLOTS);
                                table.add<float>(RING_SLOTS);
                              });
    return peak.bytes;
  }();
  /** @brief Alternate shifts with nested candidate-table arrays in scratch B.
   */
  static constexpr size_t SCRATCH_B_PEAK_BYTES = [] {
    ArenaFootprint peak;
    ArenaFootprint table;
    Alternates::add_footprint(peak, table, RING_SLOTS, ALTERNATE_CAPACITY);
    peak.template add<float>(CandidateTable::CHUNKS);
    peak.template add<float>(CandidateTable::CHUNKS);
    return peak.bytes;
  }();
  static_assert(SCRATCH_A_PEAK_BYTES <= DEFAULT_SCRATCH_A_SIZE,
                "DisplacementField swept staging exceeds scratch_a budget");
  static_assert(SCRATCH_B_PEAK_BYTES <= DEFAULT_SCRATCH_B_SIZE,
                "DisplacementField alternate shifts and candidates exceed "
                "scratch_b budget");
};
