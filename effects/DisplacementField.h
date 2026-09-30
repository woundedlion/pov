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
#include "core/color/effect_palette_recipes.h"
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
 * @details Rings share one axis and are spaced evenly in colatitude across the
 * whole sphere. Each ring vertex is displaced along
 * the stack axis by the summed displacement fields sampled at the vertex's
 * world-space position. The noise phase opens the effect and fades in from
 * zero before dwelling at full strength, then fades out into a ball phase.
 * Cap-shaped ball bumps spawn at the stack pole on random meridians and fall to
 * the opposite pole, bowing the rings away from each ball's center.
 * Fragments are shaded from a circular analogous palette that spins across the
 * stack, with hue rotated proportionally to the local displacement magnitude; a
 * ColorWipe slowly fades the palette to a freshly generated one every ~11
 * seconds. Orientation random-walks over time.
 */
template <int W, int H> class DisplacementField : public Effect {
  friend struct ::hs_test::effects_tests::DisplacementFieldWhiteBox;

public:
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
    hue_table = persistent_arena.allocate_n<Pixel>(HUE_TABLE_SIZE + 1);
    ball_colat = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_reach = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_scale = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_params =
        persistent_arena.allocate_n<const Animation::BumpParams *>(MAX_BALLS);
    ball_local = persistent_arena.allocate_n<int>(MAX_BALLS);
    ball_cv = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_rho = persistent_arena.allocate_n<float>(MAX_BALLS);
    ball_azimuth = persistent_arena.allocate_n<float>(MAX_BALLS);
    shift_pool = persistent_arena.allocate_n<float>(RING_SLOTS * (W + 1));
    // Zeroed, not raw: a culled azimuth chunk leaves its columns unbaked.
    hue_pool = persistent_arena.make_n<Pixel>(RING_SLOTS * (W + 1));
    slot_frag_alpha = persistent_arena.allocate_n<float>(RING_SLOTS);
    slot_lut_nf = persistent_arena.allocate_n<float>(RING_SLOTS);
    slot_by_ring = persistent_arena.allocate_n<int8_t>(RING_SLOTS);
    shape_storage = persistent_arena.make<ShapeStorage>();
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
   * NOISE -> BALLS phase machine under the master-gain fade, then renders one
   * frame. The phase machine holds while animations are paused; the rings still
   * render.
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

    // Ball events freeze with anims_paused and never retire, so the phase
    // machine freezes with them: spawning against frozen balls saturates the
    // pool and the zero-active exit never arrives.
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
   * @brief Marks the azimuth chunks of one ring's bake that can reach the clip.
   * @param basis Ring frame; basis.v is the stack axis.
   * @param theta Ring colatitude.
   * @param cos_t Cosine of theta.
   * @param sin_t Sine of theta.
   * @param band_r World-angle radius of the ring's displaced band.
   * @return Bit c set for chunk c of BAKE_CHUNKS, 0 when the whole ring misses
   * the clip and CHUNK_MASK when the pad spans the ring.
   * @details Chunk c owns bake columns [ceil(c * lut_n / BAKE_CHUNKS),
   * ceil((c + 1) * lut_n / BAKE_CHUNKS)); the boundaries round up, so a chunk's
   * span sits up to one column later than the even split whose midpoint this
   * test samples. A clear bit skips that span's field and hue bake and leaves
   * its columns stale. The raw per-chunk clip test is therefore widened by
   * pad_chunks neighbors on both sides, since the rasterizer's soft stroke
   * reaches params.thickness of azimuth away from a knot — that many chunks at
   * the band's smallest circumference. pad_chunks is at least 1, so the same
   * widening also absorbs the one-column rounding overhang.
   */
  __attribute__((always_inline)) uint32_t
  visible_chunk_mask(const math::Basis &basis, float theta, float cos_t,
                     float sin_t, float band_r) const {
    HS_PROFILE(df_chunk_cull);
    const float chunk_reach = (math::PI_F / BAKE_CHUNKS) * sin_t + band_r;
    const float sin_reach = sinf(fminf(chunk_reach, math::PI_F));
    uint32_t raw = 0u;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      math::Vector mid =
          (basis.v * cos_t) +
          ((basis.u * chunk_cos[c]) + (basis.w * chunk_sin[c])) * sin_t;
      if (Plot::cap_may_touch_clip<H>(clip(), mid, chunk_reach, sin_reach))
        raw |= 1u << c;
    }
    if (!raw)
      return 0u;
    const float th_lo = theta - band_r;
    const float th_hi = theta + band_r;
    int pad_chunks = BAKE_CHUNKS;
    if (th_lo > 0.0f && th_hi < math::PI_F) {
      float sin_lo = fminf(sinf(th_lo), sinf(th_hi));
      // A band hugging a pole drives sin_lo to zero; clamp before the cast.
      const float pad_f =
          ceilf(params.thickness * BAKE_CHUNKS / (2.0f * math::PI_F * sin_lo));
      pad_chunks = 1 + static_cast<int>(hs::clamp(
                           pad_f, 0.0f, static_cast<float>(BAKE_CHUNKS)));
    }
    if (2 * pad_chunks >= BAKE_CHUNKS)
      return CHUNK_MASK;
    uint32_t visible = raw;
    for (int k = 1; k <= pad_chunks; ++k)
      visible |= (raw << k) | (raw >> (BAKE_CHUNKS - k)) | (raw >> k) |
                 (raw << (BAKE_CHUNKS - k));
    return visible & CHUNK_MASK;
  }

  /**
   * @brief Chooses how one ring's bake evaluates its hue rotation.
   * @param lut_n Bake columns for the ring.
   * @param visible Chunk mask from visible_chunk_mask.
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
    int x_begin = 0;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const int x_end = ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
      if (visible & (1u << c))
        visible_samples += x_end - x_begin;
      x_begin = x_end;
    }
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
   * @details Per ring, the displacement stack is baked per azimuth column into
   * a pooled slot: the centerline shift knots together with the hue-rotated
   * ring color. The bake evaluates only the balls whose support can reach the
   * ring's colatitude (centers and rings share the stack axis), each only
   * across the azimuth arc its cap covers; the noise octaves are sampled on
   * every other knot and spline-filled between. A ring nothing can displace
   * takes a constant LUT. The LUT resolution is adaptive: enough samples for
   * the finest active feature along the ring's actual circumference. Under a
   * partial clip, rings whose displaced band cannot touch the clip are skipped
   * whole, and invisible azimuth chunks skip the field/hue bake. The hue table
   * is filled only up to the largest shift the ring actually reaches. The
   * baked rings then rasterize as soft SDF strokes with a quintic
   * cross-section falloff against the exact distance to each knot polyline, in
   * a single fused scan (Scan::DistortedRingStack) that hoists the shared-axis
   * pixel frame out of the per-ring distance evaluation.
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

    const ShapeView shapes{shape_storage};
    int n_slots = 0;
    for (int i = 0; i < n_rings; ++i)
      slot_by_ring[i] = -1;

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

      float *slut = shift_pool + n_slots * (W + 1);
      Pixel *hlut = hue_pool + n_slots * (W + 1);

      int lut_n;
      if (band + noise_bound <= 0.0f) {
        // Flat rings take the zero-knot LUT path even under -Os: the fused
        // candidate loop needs every ring in the shared pass to keep per-pixel
        // blend order.
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

        lut_n = hs::clamp(
            static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * math::PI_F *
                                   feature_scale * sin_t)),
            LUT_MIN_SAMPLES, W);
        const Animation::NoiseProductParams *octaves =
            noise_field.active_count() == 1 &&
                    std::fabs(noise_field.active_params(0).amplitude) > 0.001f
                ? &noise_field.active_params(0)
                : nullptr;
        if (octaves)
          lut_n = (lut_n + OCTAVE_GRID - 1) / OCTAVE_GRID * OCTAVE_GRID;

        uint32_t visible = CHUNK_MASK;
        if (try_cull) {
          const float band_r = band + noise_bound + params.thickness + pad;
          visible = visible_chunk_mask(basis, theta, cos_t, sin_t, band_r);
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

        // Culled chunks keep hlut stale; the pad_chunks widening of `visible`
        // above keeps a rasterized pixel from ever sampling their columns.
        // Their shifts are zeroed, since DistortedRing's constructor scans
        // every knot and stale cells would perturb its shift bounds.
        float max_shift = 0.0f;
        {
          int x = 0;
          for (int c = 0; c < BAKE_CHUNKS; ++c) {
            const int x_end = ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
            if (visible & (1u << c)) {
              for (; x < x_end; ++x)
                max_shift = fmaxf(max_shift, std::fabs(slut[x]));
            } else {
              for (; x < x_end; ++x)
                slut[x] = 0.0f;
            }
          }
        }

        bool use_hue_table;
        bool precompute_hue_table;
        float hue_domain;
        bool cyclic_hue_table;
        uint64_t hue_table_valid[(HUE_TABLE_SIZE + 64) / 64] = {};
        select_hue_mode(lut_n, visible, (band + noise_bound) * params.hue_scale,
                        use_hue_table, precompute_hue_table, hue_domain,
                        cyclic_hue_table);
        if (precompute_hue_table) {
          HS_PROFILE(df_hue_table_prep);
          prepare_hue_table(hue_base, hue_domain,
                            hue_table_cells(max_shift * params.hue_scale,
                                            hue_domain, cyclic_hue_table));
        }
        Pixel zero_hue;
        if (params.hue_scale == 0.0f)
          zero_hue = hue_rotate(hue_base, 0.0f).color;

        auto hue_for_shift = [&](float shift) HS_O3_FN {
          if (params.hue_scale == 0.0f)
            return zero_hue;
          const float amount = std::fabs(shift) * params.hue_scale;
          if (precompute_hue_table)
            return sample_hue_table(amount, hue_domain, cyclic_hue_table);
          if (use_hue_table)
            return sample_hue_table_cached(amount, hue_domain, cyclic_hue_table,
                                           hue_base, hue_table_valid);
          return hue_rotate(hue_base, amount).color;
        };

        {
          int x = 0;
          for (int c = 0; c < BAKE_CHUNKS; ++c) {
            const int x_end = ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
            if (visible & (1u << c)) {
              for (; x < x_end; ++x)
                hlut[x] = hue_for_shift(slut[x]);
            } else {
              x = x_end;
            }
          }
        }
        slut[lut_n] = slut[0];
        hlut[lut_n] = hlut[0];
      }

      ::new (static_cast<void *>(&(*shape_storage)[n_slots].ring))
          SDF::DistortedRing(basis, radius, params.thickness, slut, lut_n, 0.0f,
                             nullptr);
      slot_lut_nf[n_slots] = static_cast<float>(lut_n);
      slot_frag_alpha[n_slots] = ring_color.alpha * opacity * params.alpha;
      slot_by_ring[i] = static_cast<int8_t>(n_slots);
      ++n_slots;
    }

    if (n_slots == 0)
      return;

    // v2 is the stroke coverage the scan applies again on plot, so the ring
    // edge ramps as coverage squared. The stack hands v0 over in [0, 1).
    auto ring_shader = [this](int s, const math::Vector &, Fragment &f) {
      const Pixel *hue = hue_pool + s * (W + 1);
      float x = f.v0 * slot_lut_nf[s];
      int j = static_cast<int>(x);
      f.color = Color4(
          hue[j].lerp16(hue[j + 1], frac_to_q16(math::quintic_kernel(x - j))),
          slot_frag_alpha[s] * f.v2);
    };
    HS_PROFILE(df_fused_scan);
    Scan::DistortedRingStack::draw<W, H>(filters, canvas, n_rings, shapes,
                                         slot_by_ring, n_slots, *candidates,
                                         ring_shader);

    // ScalarFn's inplace_function member is not trivially destructible;
    // placement-built shapes must be destroyed before the storage is reused.
    for (int s = 0; s < n_slots; ++s)
      shapes[s].~DistortedRing();
  }

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
   * @param visible Chunk mask from visible_chunk_mask.
   * @param n_local Balls that can reach the ring (ball_local).
   * @param slut Receives the shift of every knot in a visible chunk.
   * @details Octave k is evaluated every OCTAVE_STRIDE[k] knots and filled in
   * between by a Catmull-Rom spline over the four surrounding grid knots, so a
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
    int x = 0;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const int x_end = ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
      const uint8_t v = static_cast<uint8_t>((visible >> c) & 1u);
      for (; x < x_end; ++x)
        knot_visible[x] = v;
    }
    auto wrap = [lut_n](int k) {
      return k < 0 ? k + lut_n : k >= lut_n ? k - lut_n : k;
    };
    // knot_near[x]: a visible knot lies within the widest spline reach of x.
    constexpr int REACH = 2 * (D1 > D2 ? D1 : D2) - 1;
    int in_window = 0;
    for (int o = -REACH; o <= REACH; ++o)
      in_window += knot_visible[wrap(o)];
    for (x = 0; x < lut_n; ++x) {
      knot_near[x] = in_window > 0;
      in_window +=
          knot_visible[wrap(x + REACH + 1)] - knot_visible[wrap(x - REACH)];
    }

    HS_PROFILE(df_octave_noise);
    float cos_a = 1.0f;
    float sin_a = 0.0f;
    for (x = 0; x < lut_n; ++x) {
      const bool vis = knot_visible[x] != 0;
      const bool near = knot_near[x] != 0;
      const bool g1 = near && x % D1 == 0;
      const bool g2 = D2 == 1 ? vis : near && x % D2 == 0;
      if (vis || g1 || g2) {
        math::Vector p =
            (basis.v * cos_t) + ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t;
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
      float next_cos = cos_a * cos_d - sin_a * sin_d;
      sin_a = sin_a * cos_d + cos_a * sin_d;
      cos_a = next_cos;
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
    for (x = 0; x < lut_n; ++x)
      if (knot_visible[x])
        slut[x] = slut[x] + np.amplitude * sample(octave1, D1, x) *
                                sample(octave2, D2, x);
  }

  /**
   * @brief Counts the hue-table knots a ring's samples can read.
   * @param max_amount Largest hue offset (turns) any sample asks for.
   * @param domain Hue-turn interval covered by the table.
   * @param cyclic Whether amounts past the domain wrap.
   * @return Knots 0..count - 1 are the only ones sample_hue_table() reads.
   * @details Samples at amount <= max_amount land at or below its cell, so the
   * table past that cell's upper knot is never read. A wrapping lookup that
   * can reach a full turn, or a domain that is not positive, reads anywhere.
   */
  static int hue_table_cells(float max_amount, float domain, bool cyclic) {
    if (!(domain > 0.0f) || (cyclic && max_amount >= domain))
      return HUE_TABLE_SIZE + 1;
    const float x = hs::clamp(max_amount / domain, 0.0f, 1.0f) * HUE_TABLE_SIZE;
    return std::min(static_cast<int>(x) + 2, HUE_TABLE_SIZE + 1);
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
   * @param visible Chunk mask from visible_chunk_mask.
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
    int x = 0;
    for (int c = 0; c < BAKE_CHUNKS; ++c) {
      const int x_end = ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
      const uint8_t v = static_cast<uint8_t>((visible >> c) & 1u);
      for (; x < x_end; ++x)
        knot_visible[x] = v;
    }
    float cos_a = 1.0f;
    float sin_a = 0.0f;
    for (x = 0; x < lut_n; ++x) {
      if (knot_visible[x]) {
        knot_pos[x] =
            (basis.v * cos_t) + ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t;
        num[x] = 0.0f;
        den[x] = 0.0f;
      }
      float next_cos = cos_a * cos_d - sin_a * sin_d;
      sin_a = sin_a * cos_d + cos_a * sin_d;
      cos_a = next_cos;
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
          num[xi] += f * f * f;
          den[xi] += f * f;
        }
        if (++xi == lut_n)
          xi = 0;
      }
    }

    for (x = 0; x < lut_n; ++x)
      if (knot_visible[x])
        slut[x] = (den[x] > FIELD_DOMINANT_DEN_EPS ? num[x] / den[x] : 0.0f) +
                  noise_field.field(knot_pos[x]);
  }

  /**
   * @brief Fills the first @p count knots of the hue table.
   * @param base Ring color's precomputed OKLab base.
   * @param domain Hue-turn interval the full table covers.
   * @param count Knots to fill, at most HUE_TABLE_SIZE + 1.
   */
  HS_O3_FN __attribute__((noinline)) void
  prepare_hue_table(const HueRotateBase &base, float domain,
                    int count = HUE_TABLE_SIZE + 1) {
    for (int i = 0; i < count; ++i)
      hue_table[i] =
          hue_rotate(base, domain * (static_cast<float>(i) / HUE_TABLE_SIZE))
              .color;
  }

  /**
   * @brief Interpolates the hue table at a hue offset.
   * @param amount Hue offset in turns.
   * @param domain Hue-turn interval covered by the table.
   * @param cyclic Whether amounts past the domain wrap instead of clamping.
   * @param ensure Called with every knot index read, before the read.
   * @return The interpolated hue-rotated ring color.
   */
  template <typename Ensure>
  Pixel sample_hue_table_with(float amount, float domain, bool cyclic,
                              Ensure ensure) const {
    float t = amount / domain;
    t = cyclic ? math::wrap_t(t) : hs::clamp(t, 0.0f, 1.0f);
    float x = t * HUE_TABLE_SIZE;
    if (x >= HUE_TABLE_SIZE) {
      ensure(HUE_TABLE_SIZE);
      return hue_table[HUE_TABLE_SIZE];
    }
    int i = static_cast<int>(x);
    ensure(i);
    ensure(i + 1);
    return hue_table[i].lerp16(hue_table[i + 1], frac_to_q16(x - i));
  }

  /** @brief Samples a hue table already fully baked by prepare_hue_table(). */
  HS_O3_FN Pixel sample_hue_table(float amount, float domain,
                                  bool cyclic) const {
    return sample_hue_table_with(amount, domain, cyclic, [](int) {});
  }

  /** @brief Samples the hue table, baking only the knots it reads and marking
   * them in the `valid` bitset. */
  HS_O3_FN Pixel sample_hue_table_cached(float amount, float domain,
                                         bool cyclic, const HueRotateBase &base,
                                         uint64_t *valid) {
    return sample_hue_table_with(amount, domain, cyclic, [&](int index) {
      const uint64_t bit = uint64_t{1} << (index & 63);
      uint64_t &word = valid[index >> 6];
      if (!(word & bit)) {
        hue_table[index] =
            hue_rotate(base,
                       domain * (static_cast<float>(index) / HUE_TABLE_SIZE))
                .color;
        word |= bit;
      }
    });
  }

  /**
   * @brief Builds a fresh random palette for the next wipe.
   * @details Each construction reseeds, so every cycle fades toward a distinct
   * palette.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::displacement_field(
        EffectPaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Spawns one falling ball with a random meridian, footprint, and
   * speed drawn from the Speed Min/Max sliders; dropped safely if the ball pool
   * or the timeline is full.
   * @details A full pool is logged before the spawn rather than by testing
   * spawn()'s result: consuming that return value costs ~960 B of ITCM, which
   * the phantasm budget cannot spare. A timeline-full drop is already logged by
   * the timeline itself. A saturated pool drops most of a phase's spawns, and
   * hs::log blocks on the serial write, so only the first drop of each ball
   * phase is logged.
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

  // Each in-flight ball is one timeline event, so at high Ball Rate x slow
  // Speed the spawner saturates the pool and drops spawns safely instead of
  // starving the effect's own events.
  static constexpr int MAX_BALLS =
      56; /**< Concurrent falling-ball pool slots. */
  static constexpr int RESERVED_EVENTS =
      6; /**< Non-ball timeline events: pinned noise field, ring Sprite, orientation RandomWalk, palette PeriodicTimer, one master-gain Transition, one palette ColorWipe. */
  static_assert(MAX_BALLS + RESERVED_EVENTS <= Timeline::MAX_EVENTS,
                "DisplacementField: a full ball pool plus the effect's own "
                "events exceeds the shared timeline budget");
  static constexpr int BALL_PHASE_FRAMES =
      900; /**< Ball-phase spawning window (~56 s at 16 fps); balls keep coming the whole window. */
  static constexpr float BALL_RATE_FPS =
      60.0f; /**< Frames per Ball Rate / Speed slider unit; one unit spans ~3.75 s at the 16 fps device cadence. */
  static constexpr float BALL_DRAPE_PER_AMPLITUDE =
      4.0f; /**< Drape gain per Ball Amp unit: the 0.1 default gives gain 0.4. */
  static constexpr int NOISE_FADE_FRAMES =
      150; /**< Noise amplitude ramp on each phase handoff. */
  static constexpr int NOISE_HOLD_FRAMES =
      600; /**< Full-noise dwell before fading out into the next ball phase. */
  /** @brief Un-oriented axis the ring stack and every ball fall share. */
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
      180; /**< Palette rollover period (~11 s at the 16 fps cadence). */
  static constexpr int PALETTE_WIPE_FRAMES =
      168; /**< Wipe duration; slightly under the cycle so a wipe is never still in flight when the next rollover fires. */
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
  static_assert(RING_SLOTS <= INT8_MAX,
                "slot_by_ring is int8_t with -1 as the culled sentinel; a "
                "larger pool wraps slot indices negative");
  float *shift_pool =
      nullptr; /**< RING_SLOTS x (W + 1) pooled shift LUTs, one slot per drawn ring; entry lut_n repeats entry 0 to close the polyline. */
  Pixel *hue_pool =
      nullptr; /**< RING_SLOTS x (W + 1) pooled hue-rotated ring colors, aligned with shift_pool. */
  float *slot_frag_alpha =
      nullptr; /**< Per-slot fragment alpha (ring alpha x sprite fade x Alpha slider). */
  float *slot_lut_nf = nullptr; /**< Per-slot bake column count. */
  int8_t *slot_by_ring =
      nullptr; /**< Ring index -> slot, -1 for culled rings; rebuilt per frame. */
  union RingSlot {
    char empty;
    SDF::DistortedRing ring;
    RingSlot() : empty{} {}
    ~RingSlot() {}
  };
  static_assert(sizeof(RingSlot) == sizeof(SDF::DistortedRing));
  using ShapeStorage = std::array<RingSlot, RING_SLOTS>;
  struct ShapeView {
    ShapeStorage *storage;
    SDF::DistortedRing &operator[](size_t index) const {
      return (*storage)[index].ring;
    }
  };
  ShapeStorage *shape_storage = nullptr;
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
  Pixel *hue_table =
      nullptr; /**< HUE_TABLE_SIZE + 1 dynamic or cyclic hue samples for the current ring. */
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

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
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

  // init() allocates the per-slot bake pools, the ball prefilter scratch, the
  // hue table, the chunk-azimuth table, the ring shapes, the scan's candidate
  // table, and both transformer pools from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      RING_SLOTS * (W + 1) * (sizeof(float) + sizeof(Pixel)) +
      RING_SLOTS *
          (2 * sizeof(float) + sizeof(int8_t) + sizeof(SDF::DistortedRing)) +
      sizeof(CandidateTable) +
      (W + 1) * (2 * sizeof(float) + 2 + sizeof(math::Vector)) +
      MAX_BALLS * (6 * sizeof(float) + sizeof(int) +
                   sizeof(const Animation::BumpParams *)) +
      (HUE_TABLE_SIZE + 1) * sizeof(Pixel) + 2 * BAKE_CHUNKS * sizeof(float) +
      MAX_BALLS * (sizeof(typename decltype(balls)::Entity) + sizeof(int)) +
      (sizeof(typename decltype(noise_field)::Entity) + sizeof(int));
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "DisplacementField persistent footprint exceeds the default "
                "partition; retune RING_SLOTS/MAX_BALLS or carve arenas");
};
