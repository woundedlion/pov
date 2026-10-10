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

#include <numeric>

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
  /// Stable persisted effect ID; seeds the effect RNG stream.
  static constexpr const char *EFFECT_ID = "DisplacementField";

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
    knot_visible = persistent_arena.allocate_n<uint8_t>(W + 1);
    knot_pos = persistent_arena.allocate_n<math::Vector>(W + 1);
    knot_near = persistent_arena.allocate_n<uint8_t>(W + 1);
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

    timeline.add(0, Animation::PeriodicTimer(
                        PALETTE_CYCLE_FRAMES,
                        [this](Canvas &) { this->roll_palette(); }, true));

    enter_noise();
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, advances the
   * NOISE -> BALLS phase machine, then renders one frame.
   * @details The phase machine holds while animations are paused.
   */
  void draw_frame() override {
    color_spin = math::wrap_t(color_spin + COLOR_SPIN_RATE);

    balls.template_params.amplitude =
        params.ball_amp * BALL_DRAPE_PER_AMPLITUDE;
    noise_field.template_params.amplitude =
        phase == Phase::NOISE ? params.noise_amp * master_gain : 0.0f;
    noise_field.template_params.scale1 = params.scale1;
    noise_field.template_params.scale2 = params.scale2;
    noise_field.template_params.speed = params.flow_speed;
    {
      HS_PROFILE(df_prepare_fields);
      balls.prepare_frame();
      noise_field.prepare_frame();
    }

    // Ball events freeze while paused and never retire, so the phase machine
    // freezes with them.
    if (!anims_paused) {
      switch (phase) {
      case Phase::BALLS:
        // Slider-paced spawner: the cooldown re-reads Ball Rate on every spawn.
        if (ball_phase_left > 0) {
          --ball_phase_left;
          spawn_cooldown -= 1.0f;
          while (spawn_cooldown <= 0.0f) {
            spawn_ball();
            spawn_cooldown +=
                hs::rand_f(0.5f, 1.5f) * BALL_RATE_FPS / params.ball_rate;
          }
        } else if (balls.active_count() == 0) {
          enter_noise();
        }
        break;
      case Phase::NOISE:
        if (master_gain >= 1.0f && noise_hold > 0 && --noise_hold == 0)
          timeline.add(0, Animation::Transition(master_gain, 0.0f,
                                                NOISE_FADE_FRAMES,
                                                math::ease_in_out_sin));
        if (noise_hold == 0 && master_gain <= 0.0f)
          enter_balls();
        break;
      }
    }

    Canvas canvas(*this);
    {
      HS_PROFILE(df_timeline_step);
      timeline.step(canvas);
    }
  }

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

  /**
   * @brief Chooses how one ring's bake evaluates its hue rotation.
   * @param lut_n Bake columns for the ring.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param hue_extent Signed hue turns the ring's displacement range covers.
   * @param use_hue_table Out: sample a HUE_TABLE_SIZE-cell table instead of
   * rotating per column.
   * @param precompute_hue_table Out: fill the table up front rather than
   * populating cells on demand.
   * @param hue_domain Out: hue-turn span the table's cells cover.
   * @param cyclic_hue_table Out: the domain wraps, so lookups wrap with it.
   * @details Culled chunks bake no columns, so the precompute decision
   * amortizes the table over the visible chunks' columns rather than lut_n.
   */
  __attribute__((always_inline)) void
  select_hue_mode(int lut_n, uint32_t visible, float hue_extent,
                  bool &use_hue_table, bool &precompute_hue_table,
                  float &hue_domain, bool &cyclic_hue_table) {
    use_hue_table = params.hue_scale != 0.0f && lut_n > 2 * HUE_TABLE_SIZE;
#if HS_ENABLE_TEST_ORACLES
    use_hue_table = use_hue_table && !force_exact_hue;
#endif
    precompute_hue_table = false;
    hue_domain = 0.0f;
    cyclic_hue_table = false;
    if (!use_hue_table)
      return;
    int visible_samples = 0;
    Plot::for_each_chunk_span<BAKE_CHUNKS>(
        lut_n, visible,
        [&](int begin, int end)
            __attribute__((always_inline)) { visible_samples += end - begin; },
        [](int, int) __attribute__((always_inline)) {});
    precompute_hue_table = visible_samples > 2 * HUE_TABLE_SIZE;
    cyclic_hue_table = std::fabs(hue_extent) > 1.0f;
    hue_domain =
        cyclic_hue_table ? std::copysign(1.0f, hue_extent) : hue_extent;
#if HS_ENABLE_TEST_HOOKS
    ++hue_table_uses;
#endif
  }

  /**
   * @brief Bakes every ring over the active ball pool, then rasterizes the
   * stack in one fused scan.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite fade multiplied into each fragment's alpha.
   * @details Under a partial clip, rings that cannot touch the clip and
   * invisible azimuth chunks skip the bake.
   */
  HS_O3_FN void draw_rings(Canvas &canvas, float opacity) {
    HS_PROFILE(df_draw_rings);
    const int n_rings = params.num_rings;
    HS_CHECK(n_rings <= RING_SLOTS,
             "DisplacementField: Rings slider exceeds the baked-ring pool");
    math::Basis basis = math::make_basis(orientation.get(), STACK_AXIS);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the rasterizer's soft stroke cross-section past
    // the clip's margin.
    const float pad = 3.0f * math::PI_F / H;
    const float noise_bound = noise_field.field_bound();

    const int n_balls = balls.active_count();
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

    const float noise_feature =
        noise_bound > 0.0f ? params.scale1 + params.scale2 : 0.0f;

    rings.begin_frame();

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

      Color4 ring_color =
          palette.get(math::wrap_t((i + 0.5f) / n_rings + color_spin));
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

        const Animation::NoiseProductParams *octaves =
            noise_field.active_count() == 1 &&
                    std::fabs(noise_field.active_params(0).amplitude) > 0.001f
                ? &noise_field.active_params(0)
                : nullptr;
        lut_n = noise_lut_samples(feature_scale, sin_t, octaves != nullptr);

        uint32_t visible = CHUNK_MASK;
        if (try_cull) {
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

        const float dphi = 2.0f * math::PI_F / lut_n;
        const float cos_d = cosf(dphi);
        const float sin_d = sinf(dphi);

        if (octaves) {
          bake_noise_octaves(*octaves, basis, theta, cos_t, sin_t, cos_d, sin_d,
                             lut_n, visible, n_local, slut);
        } else {
          bake_ball_spans(basis, theta, cos_t, sin_t, cos_d, sin_d, lut_n,
                          visible, n_local, slut);
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

        bool use_hue_table;
        bool precompute_hue_table;
        float hue_domain;
        bool cyclic_hue_table;
        select_hue_mode(lut_n, visible, (band + noise_bound) * params.hue_scale,
                        use_hue_table, precompute_hue_table, hue_domain,
                        cyclic_hue_table);
        if (use_hue_table)
          hue_table.bind(hue_base, hue_domain, cyclic_hue_table);
        if (precompute_hue_table) {
          HS_PROFILE(df_hue_table_prep);
          hue_table.prepare(hue_table.cells(max_shift * params.hue_scale));
        }
        Pixel zero_hue;
        if (params.hue_scale == 0.0f)
          zero_hue = hue_rotate(hue_base, 0.0f).color;

        auto hue_for_shift = [&](float shift) HS_O3_FN {
          if (params.hue_scale == 0.0f)
            return zero_hue;
          const float amount = std::fabs(shift) * params.hue_scale;
          if (precompute_hue_table)
            return hue_table.sample(amount);
          if (use_hue_table)
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

      rings.commit(i, ring_color.alpha * opacity * params.alpha, lut_n, basis,
                   radius, params.thickness, slut, lut_n, 0.0f, nullptr);
    }

    if (rings.size() == 0)
      return;

    // v2 is the stroke coverage the scan applies again on plot, so the ring
    // edge ramps as coverage squared. The stack hands v0 over in [0, 1).
    auto ring_shader = [this](int s, const math::Vector &, Fragment &f) {
      const Pixel *hue = rings.hue_row(s);
      float x = f.v0 * rings.lut_columns(s);
      int j = static_cast<int>(x);
      f.color = Color4(
          hue[j].lerp16(hue[j + 1], frac_to_q16(math::quintic_kernel(x - j))),
          rings.frag_alpha(s) * f.v2);
    };
    HS_PROFILE(df_fused_scan);
    Scan::DistortedRingStack::draw<W, H>(
        filters, canvas, n_rings, rings.shapes(), rings.slot_map(),
        rings.size(), *candidates, ring_shader);
    rings.release();
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
      return (basis.v * cos_t) +
             ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t;
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
   * @param theta Ring colatitude.
   * @param cos_t Cosine of theta.
   * @param sin_t Sine of theta.
   * @param cos_d Cosine of one knot cell's azimuth step.
   * @param sin_d Sine of one knot cell's azimuth step.
   * @param lut_n Knot count, a multiple of OCTAVE_GRID.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param n_local Balls that can reach the ring (ball_local).
   * @param slut Receives the shift of every knot in a visible chunk.
   * @details The two octaves use OCTAVE1_STRIDE and OCTAVE2_STRIDE knot spacing
   * and fill between samples with a Catmull-Rom spline over four grid knots, so a
   * grid knot is evaluated when any visible knot lies within two strides of it.
   * Knot positions come from the same azimuth recurrence as the exact bake.
   */
  HS_HOT_FLASH_MEMBER void
  bake_noise_octaves(const Animation::NoiseProductParams &np,
                     const math::Basis &basis, float theta, float cos_t,
                     float sin_t, float cos_d, float sin_d, int lut_n,
                     uint32_t visible, int n_local, float *slut) {
    constexpr int D1 = OCTAVE1_STRIDE;
    constexpr int D2 = OCTAVE2_STRIDE;
    mark_visible_knots(lut_n, visible);
    auto wrap = [lut_n](int k) {
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
    RingKnotWalk knots{basis, cos_t, sin_t, cos_d, sin_d};
    for (int x = 0; x < lut_n; ++x, knots.advance()) {
      const bool vis = knot_visible[x] != 0;
      const bool near = knot_near[x] != 0;
      const bool g1 = near && x % D1 == 0;
      const bool g2 = D2 == 1 ? vis : near && x % D2 == 0;
      if (vis || g1 || g2) {
        math::Vector p = knots.position();
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

    auto sample = [&](const float *oct, int d, int k) {
      const int r = k % d;
      if (r == 0)
        return oct[k];
      const int b = k - r;
      const float t = static_cast<float>(r) / d;
      const float p0 = oct[wrap(b - d)], p1 = oct[b], p2 = oct[wrap(b + d)],
                  p3 = oct[wrap(b + 2 * d)];
      return p1 + 0.5f * t *
                      (p2 - p0 +
                       t * (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3 +
                            t * (3.0f * (p1 - p2) + p3 - p0)));
    };
    for (int x = 0; x < lut_n; ++x)
      if (knot_visible[x])
        slut[x] = slut[x] + np.amplitude * sample(octave1, D1, x) *
                                sample(octave2, D2, x);
  }

  /**
   * @brief Bakes one ring's centerline shifts from the balls, each evaluated
   *        only across the knots its cap can cover.
   * @param basis Ring frame; basis.v is the stack axis.
   * @param theta Ring colatitude.
   * @param cos_t Cosine of theta.
   * @param sin_t Sine of theta.
   * @param cos_d Cosine of one knot cell's azimuth step.
   * @param sin_d Sine of one knot cell's azimuth step.
   * @param lut_n Knot count.
   * @param visible Chunk mask from Plot::visible_chunk_mask.
   * @param n_local Balls that can reach the ring (ball_local).
   * @param slut Receives the shift of every knot in a visible chunk.
   * @details On the ring, a ball's cap test dot(p, center) > cos_radius reads
   * cos_t * cv + sin_t * rho * cos(a - a0) > cos_radius, which holds on one
   * azimuth arc around a0. Knots off every padded arc get exactly the zero a
   * rejected ball adds, and each knot folds its balls in ball_local order, so
   * the shifts match evaluating every ball at every knot.
   */
  HS_HOT_FLASH_MEMBER void
  bake_ball_spans(const math::Basis &basis, float theta, float cos_t,
                  float sin_t, float cos_d, float sin_d, int lut_n,
                  uint32_t visible, int n_local, float *slut) {
    // Covers the recurrence's drift from the analytic knot azimuth and the
    // float slop in the cap dot product.
    constexpr float COS_MARGIN = 2e-4f;
    float *num = octave1;
    float *den = octave2;
    mark_visible_knots(lut_n, visible);
    RingKnotWalk knots{basis, cos_t, sin_t, cos_d, sin_d};
    for (int x = 0; x < lut_n; ++x, knots.advance()) {
      if (knot_visible[x]) {
        knot_pos[x] = knots.position();
        num[x] = 0.0f;
        den[x] = 0.0f;
      }
    }

    const float cells_per_radian = lut_n / (2.0f * math::PI_F);
    for (int j = 0; j < n_local; ++j) {
      const int k = ball_local[j];
      const Animation::BumpParams &bp = *ball_params[k];
      const float a_term = cos_t * ball_cv[k];
      const float b_term = sin_t * ball_rho[k];
      int first = 0;
      int count = lut_n;
      if (b_term > COS_MARGIN) {
        const float q = (bp.cos_radius - COS_MARGIN - a_term) / b_term;
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
          const float f = bump_field_with_y(knot_pos[xi], bp, y);
          DominantFieldAccumulator::accumulate(num[xi], den[xi], f);
        }
        if (++xi == lut_n)
          xi = 0;
      }
    }

    for (int x = 0; x < lut_n; ++x)
      if (knot_visible[x])
        slut[x] = DominantFieldAccumulator::resolve(num[x], den[x]) +
                  noise_field.field(knot_pos[x]);
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
   * @details Only the first pool-full drop of each ball phase is logged.
   */
  HS_COLD_MEMBER void spawn_ball() {
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
    balls.spawn_pausable(&anims_paused, 0, orientation, STACK_AXIS,
                         hs::rand_f(0.0f, 2.0f * math::PI_F), fall_frames);
  }

  /**
   * @brief Enters the noise phase and fades it in from zero.
   */
  HS_COLD_MEMBER void enter_noise() {
    phase = Phase::NOISE;
    master_gain = 0.0f;
    noise_hold = NOISE_HOLD_FRAMES;
    timeline.add(0, Animation::Transition(master_gain, 1.0f, NOISE_FADE_FRAMES,
                                          math::ease_in_out_sin));
  }

  /**
   * @brief Enters the ball phase, opening a fresh spawning window.
   */
  HS_COLD_MEMBER void enter_balls() {
    phase = Phase::BALLS;
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

  /** @brief Displacement-phase state: the noise field or falling balls. */
  enum class Phase { BALLS, NOISE };

  Phase phase = Phase::
      NOISE; /**< Current displacement phase; the effect opens on noise. */
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
    void begin_frame() {
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
    __attribute__((always_inline)) void
    commit(int ring, float frag_alpha, int lut_n, RingArgs &&...ring_args) {
      assert(count < SLOTS && "DisplacementField: ring pool is full");
      const int s = count;
      ::new (static_cast<void *>(&(*storage)[s].ring))
          SDF::DistortedRing(std::forward<RingArgs>(ring_args)...);
      columns[s] = static_cast<float>(lut_n);
      alphas[s] = frag_alpha;
      slot_by_ring[ring] = static_cast<int8_t>(s);
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
    ShapeView shapes() const { return {storage}; }
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
  RingPool rings; /**< Rings drawn this frame. */
  using CandidateTable = Scan::DistortedRingStack::CandidateTable<W, H>;
  CandidateTable *candidates =
      nullptr; /**< Fused scan's per-frame ring candidate map. */
  float *octave1 =
      nullptr; /**< W + 1 envelope-octave samples on its knot grid. */
  float *octave2 =
      nullptr; /**< W + 1 detail-octave samples on its knot grid. */
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
      return sample_with(amount, [&](int index) {
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
#endif
#if HS_ENABLE_TEST_HOOKS
  int hue_table_uses = 0;
#endif

  /** Ring-to-ball prefilter pad absorbing fast_acos and tangent-recurrence
   *  rounding (radians); a ball excluded despite the pad fields exactly 0
   *  everywhere on the ring. */
  static constexpr float BALL_TOUCH_EPS = 1e-3f;

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
  } params;

  static_assert(Params{}.thickness >= 0.4f * THICKNESS_PX &&
                    Params{}.thickness <= 6.0f * THICKNESS_PX,
                "DisplacementField Thickness default falls outside its "
                "W-scaled slider range at this build resolution");

  /** @brief Persistent bytes of the knot bake scratch, alignment included. */
  static constexpr size_t KNOT_SCRATCH_BYTES =
      (W + 1) *
          (2 * sizeof(float) + 2 * sizeof(uint8_t) + sizeof(math::Vector)) +
      2 * alignof(float) + 2 * alignof(uint8_t) + alignof(math::Vector);
  /** @brief Persistent bytes of the per-frame ball cache. */
  static constexpr size_t BALL_CACHE_BYTES =
      MAX_BALLS * (6 * sizeof(float) + sizeof(int) +
                   sizeof(const Animation::BumpParams *)) +
      6 * alignof(float) + alignof(int) +
      alignof(const Animation::BumpParams *);
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
      BALL_CACHE_BYTES + CHUNK_TRIG_BYTES + CANDIDATE_BYTES + FIELD_BYTES;
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "DisplacementField persistent footprint exceeds the default "
                "partition; retune RING_SLOTS/MAX_BALLS or carve arenas");
};
