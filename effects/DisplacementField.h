/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

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
 * ColorWipe slowly fades the palette to a freshly generated one every few
 * seconds. Orientation random-walks over time.
 */
template <int W, int H> class DisplacementField : public Effect {
  friend struct ::hs_test::effects_tests::DisplacementFieldWhiteBox;

public:
  /**
   * @brief Builds the effect with its palette and ring-stack axis.
   */
  HS_COLD_MEMBER DisplacementField()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        balls(timeline), noise_field(timeline),
        palette(make_palette()), next_palette(make_palette()), normal(X_AXIS) {}

  /**
   * @brief Allocates the bake LUTs, registers params, seeds the noise field,
   * and builds the timeline.
   */
  void init() override {
    hue_table = persistent_arena.allocate_n<Pixel>(HUE_TABLE_SIZE + 1);
    solid_colat = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_reach = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_shift = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_scale = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_local = persistent_arena.allocate_n<int>(SOLID_MAX);
    shift_pool = persistent_arena.allocate_n<float>(RING_SLOTS * (W + 1));
    hue_pool = persistent_arena.allocate_n<Pixel>(RING_SLOTS * (W + 1));
    slot_frag_alpha = persistent_arena.allocate_n<float>(RING_SLOTS);
    slot_lut_n = persistent_arena.allocate_n<int>(RING_SLOTS);
    slot_by_ring = persistent_arena.allocate_n<int8_t>(RING_SLOTS);
    shapes_raw = persistent_arena.allocate(
        RING_SLOTS * sizeof(SDF::DistortedRing), alignof(SDF::DistortedRing));
    balls.init_storage(persistent_arena);
    noise_field.init_storage(persistent_arena);

    noise_field.template_params.noise.SetSeed(hs::rand_int(0, 65536));

    // One pixel of azimuth in ring-space.
    const float px = 2.0f * PI_F / W;
    params.thickness = 0.04f;

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Rings", &params.num_rings, 1.0f, 72.0f);
    register_param("Thickness", &params.thickness, 0.5f * px, 6.0f * px);
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
                          this->draw_fn(canvas, opacity);
                        },
                        -1, 24, ease_linear, 0, ease_linear));

    timeline.add(0, Animation::RandomWalk<W>(orientation, normal, walk_noise));

    timeline.add(0, Animation::PeriodicTimer(
                        PALETTE_CYCLE_FRAMES,
                        [this](Canvas &) { this->roll_palette(); }, true));

    enter_noise();
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, advances the
   * NOISE -> BALLS phase machine under the master-gain fade, advances the
   * palette wipe, then renders one frame.
   */
  void draw_frame() override {
    color_spin = wrap_t(color_spin + COLOR_SPIN_RATE);
    step_palette_wipe();

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

    switch (phase) {
    case Phase::BALLS:
      // Slider-paced spawner: the cooldown re-reads Ball Rate on every spawn.
      if (ball_phase_left > 0) {
        --ball_phase_left;
        if (--spawn_cooldown <= 0) {
          spawn_ball();
          spawn_cooldown = static_cast<int>(
              hs::rand_f(0.5f, 1.5f) * BALL_RATE_FPS / params.ball_rate);
        }
      } else if (balls.active_count() == 0) {
        enter_noise();
      }
      break;
    case Phase::NOISE:
      if (master_gain >= 1.0f && noise_hold > 0 && --noise_hold == 0)
        timeline.add(0, Animation::Transition(master_gain, 0.0f,
                                              NOISE_FADE_FRAMES,
                                              ease_in_out_sin));
      if (noise_hold == 0 && master_gain <= 0.0f)
        enter_balls();
      break;
    }

    // IIFE so the HS_PROFILE scope measures the buffer_free() spin-wait inside
    // the Canvas constructor without also covering the timeline step.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(df_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(df_timeline_step);
      timeline.step(canvas);
    }
  }

  /**
   * @brief Draws all rings for this frame, over whichever solid-body pool is
   * active.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details The shared per-ring machinery runs over the active ball pool.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    draw_rings(canvas, opacity);
  }

private:
  /** @brief Evaluates the active ball fields using cached ring geometry. */
  float ball_field(const Vector &p, const int *ks, int n, float theta) const {
    float num = 0.0f;
    float den = 0.0f;
    for (int j = 0; j < n; ++j) {
      const int k = ks[j];
      float f = bump_field_with_y(p, balls.active_params(k),
                                  theta - solid_colat[k]);
      num += f * f * f;
      den += f * f;
    }
    return den > 1e-9f ? num / den : 0.0f;
  }

  /**
   * @brief Bakes every ring over one solid-body pool, then rasterizes the
   * stack in one fused scan.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite fade multiplied into each fragment's alpha.
   * @details Per ring, the displacement stack is baked per azimuth column into
   * a pooled slot: the centerline shift knots together with the hue-rotated
   * ring color. The bake evaluates only the bodies whose support can reach the
   * ring's colatitude (centers and rings share the stack axis); a ring nothing
   * can displace takes a constant LUT. The LUT resolution is adaptive: enough
   * samples for the finest active feature along the ring's actual
   * circumference. Under a partial clip, rings whose displaced band cannot
   * touch the clip are skipped whole, and invisible azimuth chunks skip the
   * field/hue bake. The baked rings then rasterize as soft SDF strokes with a
   * quintic cross-section falloff against the exact distance to each knot
   * polyline, in a single fused scan (Scan::DistortedRingStack) that hoists
   * the shared-axis pixel frame out of the per-ring distance evaluation.
   * Debug visuals fall back to per-ring rasterizes so the bounding-box tint
   * keeps per-shape scan bounds.
   */
  void draw_rings(Canvas &canvas, float opacity) {
    HS_PROFILE(df_draw_rings);
    int n_rings = static_cast<int>(params.num_rings);
    assert(n_rings <= RING_SLOTS);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float noise_bound = noise_field.field_bound();

    const int n_balls = balls.active_count();
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      solid_colat[b] =
          fast_acos(hs::clamp(dot(basis.v, sp.center), -1.0f, 1.0f));
      solid_reach[b] = sp.field_bound();
      solid_shift[b] = sp.field_bound();
      solid_scale[b] = 2.0f / sp.radius;
    }

    const float noise_feature =
        noise_bound > 0.0f ? params.scale1 + params.scale2 : 0.0f;

    auto *shapes =
        std::launder(reinterpret_cast<SDF::DistortedRing *>(shapes_raw));
    int n_slots = 0;
    for (int i = 0; i < n_rings; ++i)
      slot_by_ring[i] = -1;

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      float theta = radius * (PI_F / 2.0f);

      int n_local = 0;
      float band = 0.0f;
      float solid_feature = 0.0f;
      for (int b = 0; b < n_balls; ++b) {
        if (std::fabs(theta - solid_colat[b]) <
            solid_reach[b] + BALL_TOUCH_EPS) {
          solid_local[n_local++] = b;
          band = std::max(band, solid_shift[b]);
          solid_feature = std::max(solid_feature, solid_scale[b]);
        }
      }

      if (try_cull &&
          !Plot::cap_may_touch_clip<H>(clip(), basis.v,
                                       theta + band + noise_bound +
                                           params.thickness + pad))
        continue;

      Color4 ring_color =
          palette.get(wrap_t((i + 0.5f) / n_rings + color_spin));
      HueRotateBase hue_base = make_hue_rotate_base(ring_color);

      float *slut = shift_pool + n_slots * (W + 1);
      Pixel *hlut = hue_pool + n_slots * (W + 1);

      float ring_bound = 0.0f;
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
        float feature_scale = std::max(noise_feature, solid_feature);
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);

        lut_n = hs::clamp(
            static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * PI_F *
                                   feature_scale * sin_t)),
            LUT_MIN_SAMPLES, W);

        uint32_t visible = CHUNK_MASK;
        if (try_cull) {
          HS_PROFILE(df_chunk_cull);
          const float band_r = band + noise_bound + params.thickness + pad;
          const float chunk_reach = (PI_F / BAKE_CHUNKS) * sin_t + band_r;
          uint32_t raw = 0u;
          for (int c = 0; c < BAKE_CHUNKS; ++c) {
            float a = (2.0f * c + 1.0f) * (PI_F / BAKE_CHUNKS);
            Vector mid = (basis.v * cos_t) +
                         ((basis.u * cosf(a)) + (basis.w * sinf(a))) * sin_t;
            if (Plot::cap_may_touch_clip<H>(clip(), mid, chunk_reach))
              raw |= 1u << c;
          }
          if (!raw)
            continue;
          const float th_lo = theta - band_r;
          const float th_hi = theta + band_r;
          int pad_chunks = BAKE_CHUNKS;
          if (th_lo > 0.0f && th_hi < PI_F) {
            float sin_lo = std::min(sinf(th_lo), sinf(th_hi));
            pad_chunks = 1 + static_cast<int>(ceilf(
                             params.thickness * BAKE_CHUNKS /
                             (2.0f * PI_F * sin_lo)));
          }
          if (2 * pad_chunks >= BAKE_CHUNKS) {
            visible = CHUNK_MASK;
          } else {
            visible = raw;
            for (int k = 1; k <= pad_chunks; ++k)
              visible |= (raw << k) | (raw >> (BAKE_CHUNKS - k)) |
                         (raw >> k) | (raw << (BAKE_CHUNKS - k));
            visible &= CHUNK_MASK;
          }
        }

        bool use_hue_table =
            params.hue_scale != 0.0f && lut_n > 2 * HUE_TABLE_SIZE;
#ifdef HS_TEST_BUILD
        use_hue_table = use_hue_table && !force_exact_hue;
#endif
        bool precompute_hue_table = false;
        float hue_domain = 0.0f;
        bool cyclic_hue_table = false;
        uint64_t hue_table_valid[2];
        if (use_hue_table) {
          int visible_samples = 0;
          int x_begin = 0;
          for (int c = 0; c < BAKE_CHUNKS; ++c) {
            const int x_end =
                ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
            if (visible & (1u << c))
              visible_samples += x_end - x_begin;
            x_begin = x_end;
          }
          precompute_hue_table = visible_samples > 2 * HUE_TABLE_SIZE;
          const float hue_extent =
              (band + noise_bound) * params.hue_scale;
          cyclic_hue_table = std::fabs(hue_extent) > 1.0f;
          hue_domain = cyclic_hue_table
                           ? std::copysign(1.0f, hue_extent)
                           : hue_extent;
#ifdef HS_TEST_BUILD
          ++hue_table_uses;
#endif
          if (precompute_hue_table) {
            HS_PROFILE(df_hue_table_prep);
            prepare_hue_table(hue_base, hue_domain);
          } else {
            hue_table_valid[0] = 0;
            hue_table_valid[1] = 0;
          }
        }
        Pixel zero_hue;
        if (params.hue_scale == 0.0f)
          zero_hue = hue_rotate(hue_base, 0.0f).color;

        auto hue_for_shift = [&](float shift) {
          if (params.hue_scale == 0.0f)
            return zero_hue;
          const float amount = std::fabs(shift) * params.hue_scale;
          if (precompute_hue_table)
            return sample_hue_table(amount, hue_domain, cyclic_hue_table);
          if (use_hue_table)
            return sample_hue_table_cached(amount, hue_domain,
                                           cyclic_hue_table, hue_base,
                                           hue_table_valid);
          return hue_rotate(hue_base, amount).color;
        };

        HS_PROFILE(df_lut_bake);

        const float dphi = 2.0f * PI_F / lut_n;
        const float cos_d = cosf(dphi);
        const float sin_d = sinf(dphi);
        float cos_a = 1.0f;
        float sin_a = 0.0f;

        int x = 0;
        for (int c = 0; c < BAKE_CHUNKS; ++c) {
          const int x_end =
              ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
          if (visible & (1u << c)) {
            for (; x < x_end; ++x) {
              Vector p = (basis.v * cos_t) +
                         ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t;
              float s = ball_field(p, solid_local, n_local, theta) +
                        noise_field.field(p);
              ring_bound = std::max(ring_bound, std::fabs(s));
              slut[x] = s;
              hlut[x] = hue_for_shift(s);
              float next_cos = cos_a * cos_d - sin_a * sin_d;
              sin_a = sin_a * cos_d + cos_a * sin_d;
              cos_a = next_cos;
            }
          } else {
            for (; x < x_end; ++x) {
              float next_cos = cos_a * cos_d - sin_a * sin_d;
              sin_a = sin_a * cos_d + cos_a * sin_d;
              cos_a = next_cos;
            }
          }
        }
        slut[lut_n] = slut[0];
        hlut[lut_n] = hlut[0];
      }

      ::new (static_cast<void *>(shapes + n_slots)) SDF::DistortedRing(
          basis, radius, params.thickness, slut, lut_n, ring_bound, 0.0f);
      slot_lut_n[n_slots] = lut_n;
      slot_frag_alpha[n_slots] = ring_color.alpha * opacity * params.alpha;
      slot_by_ring[i] = static_cast<int8_t>(n_slots);
      ++n_slots;
    }

    if (n_slots == 0)
      return;

    auto ring_shader = [this](int s, const Vector &, Fragment &f) {
      const Pixel *hue = hue_pool + s * (W + 1);
      float x = wrap_t(f.v0) * slot_lut_n[s];
      int j = static_cast<int>(x);
      f.color = Color4(
          hue[j].lerp16(hue[j + 1], frac_to_q16(quintic_kernel(x - j))),
          slot_frag_alpha[s] * f.v2);
    };
    if (canvas.debug()) {
      // Per-ring rasterizes so the bounding-box tint has per-shape scan
      // bounds; ascending slot order keeps the fused pass's blend order.
      for (int s = 0; s < n_slots; ++s) {
        shapes[s].suppress_pole_fill = true;
        Scan::rasterize<W, H>(
            filters, canvas, shapes[s],
            [&, s](const Vector &p, Fragment &f) { ring_shader(s, p, f); });
      }
    } else {
      HS_PROFILE(df_fused_scan);
      Scan::DistortedRingStack::draw<W, H>(filters, canvas, n_rings, shapes,
                                           slot_by_ring, n_slots,
                                           ring_shader);
    }

    // ScalarFn's inplace_function member is not trivially destructible;
    // placement-built shapes must be destroyed before the storage is reused.
    for (int s = 0; s < n_slots; ++s)
      shapes[s].~DistortedRing();
  }

  __attribute__((noinline)) void prepare_hue_table(const HueRotateBase &base,
                                                   float domain) {
    for (int i = 0; i <= HUE_TABLE_SIZE; ++i)
      hue_table[i] = hue_rotate(
          base, domain * (static_cast<float>(i) / HUE_TABLE_SIZE)).color;
  }

  Pixel sample_hue_table(float amount, float domain, bool cyclic) const {
    float t = amount / domain;
    t = cyclic ? wrap_t(t) : hs::clamp(t, 0.0f, 1.0f);
    float x = t * HUE_TABLE_SIZE;
    if (x >= HUE_TABLE_SIZE)
      return hue_table[HUE_TABLE_SIZE];
    int i = static_cast<int>(x);
    return hue_table[i].lerp16(hue_table[i + 1], frac_to_q16(x - i));
  }

  Pixel sample_hue_table_cached(float amount, float domain, bool cyclic,
                                const HueRotateBase &base,
                                uint64_t *valid) {
    auto ensure = [&](int index) {
      const uint64_t bit = uint64_t{1} << (index & 63);
      uint64_t &word = valid[index >> 6];
      if (!(word & bit)) {
        hue_table[index] = hue_rotate(
            base, domain * (static_cast<float>(index) / HUE_TABLE_SIZE)).color;
        word |= bit;
      }
    };

    float t = amount / domain;
    t = cyclic ? wrap_t(t) : hs::clamp(t, 0.0f, 1.0f);
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

  /**
   * @brief Builds a fresh random palette for the next wipe.
   * @details Each construction reseeds, so every cycle fades toward a distinct
   * palette.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                             BrightnessProfile::FLAT);
  }

  /**
   * @brief Spawns one falling ball with a random meridian, footprint, and
   * speed drawn from the Speed Min/Max sliders; dropped safely if the pool is
   * full.
   */
  void spawn_ball() {
    balls.template_params.radius =
        hs::rand_f(std::min(params.ball_min, params.ball_max),
                   std::max(params.ball_min, params.ball_max));
    float speed =
        hs::rand_f(std::min(params.ball_speed_min, params.ball_speed_max),
                   std::max(params.ball_speed_min, params.ball_speed_max));
    int fall_frames =
        std::max(2, static_cast<int>(BALL_RATE_FPS / speed));
    balls.spawn(0, orientation, normal, hs::rand_f(0.0f, 2.0f * PI_F),
                fall_frames);
  }

  /**
   * @brief Enters the noise phase and fades it in from zero.
   */
  HS_COLD_MEMBER void enter_noise() {
    phase = Phase::NOISE;
    master_gain = 0.0f;
    noise_hold = NOISE_HOLD_FRAMES;
    timeline.add(0, Animation::Transition(master_gain, 1.0f, NOISE_FADE_FRAMES,
                                          ease_in_out_sin));
  }

  /**
   * @brief Enters the ball phase, opening a fresh spawning window.
   */
  HS_COLD_MEMBER void enter_balls() {
    phase = Phase::BALLS;
    ball_phase_left = BALL_PHASE_FRAMES;
    spawn_cooldown = 0;
  }

  /**
   * @brief Rolls the palette toward a freshly generated one via a ColorWipe.
   * @details Skips while a previous wipe is still in flight so a second wipe
   * cannot clobber the target the live one references.
   */
  void roll_palette() {
    if (wipe_frames_remaining > 0)
      return;
    next_palette = make_palette();
    timeline.add(0, Animation::ColorWipe(palette, next_palette,
                                         PALETTE_WIPE_FRAMES, ease_linear));
    wipe_frames_remaining = PALETTE_WIPE_FRAMES;
    wipe_pending = true;
  }

  void step_palette_wipe() {
    if (wipe_pending) {
      wipe_pending = false;
    } else if (wipe_frames_remaining > 0) {
      --wipe_frames_remaining;
    }
  }

  FastNoiseLite walk_noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  // Near the timeline's 64-event budget: each in-flight ball is one event, so
  // at high Ball Rate x slow Speed the spawner saturates here and drops spawns
  // safely instead of overflowing the shared event buffer.
  static constexpr int MAX_BALLS = 56;  /**< Concurrent falling-ball pool slots. */
  static constexpr int SOLID_MAX = MAX_BALLS; /**< Shared prefilter scratch size. */
  static constexpr int BALL_PHASE_FRAMES = 900; /**< Ball-phase spawning window (~15 s); balls keep coming the whole window. */
  static constexpr float BALL_RATE_FPS = 60.0f;    /**< Frame cadence assumed by the Ball Rate and Speed sliders' per-second units. */
  static constexpr float BALL_DRAPE_PER_AMPLITUDE = 4.0f; /**< Drape gain per Ball Amp unit: the 0.25 default = gain 1, the 0.8 max saturates toward full clearance. */
  static constexpr int NOISE_FADE_FRAMES = 150; /**< Noise amplitude ramp on each phase handoff. */
  static constexpr int NOISE_HOLD_FRAMES = 600; /**< Full-noise dwell before fading out into the next solid phase. */

  BallDropTransformer<MAX_BALLS> balls;   /**< Falling-ball displacement fields. */
  NoiseProductTransformer<1> noise_field; /**< Two-octave noise displacement field. */

  GenerativePalette palette;      /**< Active palette (mutated by an in-flight ColorWipe). */
  GenerativePalette next_palette; /**< Target palette the current wipe fades toward. */
  Vector normal;
  Orientation<> orientation;

  int wipe_frames_remaining = 0; /**< Frames left in the in-flight palette wipe. */
  bool wipe_pending = false;     /**< Wipe armed this frame; it first steps next frame. */

  /** @brief Displacement-phase state: the noise field or falling balls. */
  enum class Phase { BALLS, NOISE };

  Phase phase = Phase::NOISE; /**< Current displacement phase; the effect opens on noise. */
  int ball_phase_left = BALL_PHASE_FRAMES; /**< Frames left in this ball phase's spawning window. */
  int spawn_cooldown = 0;  /**< Frames until the next ball spawn. */
  float master_gain = 0.0f; /**< Noise fade envelope in [0, 1]; gates the noise field, animated by Transitions. */
  int noise_hold = 0;      /**< Frames until the noise phase begins fading back out. */

  static constexpr int PALETTE_CYCLE_FRAMES = 180; /**< Palette rollover period (~3 s at the ~60 fps cadence). */
  static constexpr int PALETTE_WIPE_FRAMES = 168;  /**< Wipe duration; slightly under the cycle so a wipe is never still in flight when the next rollover fires. */
  static constexpr float COLOR_SPIN_RATE = 0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float LUT_SAMPLES_PER_UNIT = 8.0f; /**< Bake columns per feature-space unit of ring circumference. */
  static constexpr int LUT_MIN_SAMPLES = 16;          /**< Bake-column floor for tiny/low-scale rings. */
  static constexpr int HUE_TABLE_SIZE = 64;            /**< Hue-turn interpolation cells per ring. */
  static constexpr int BAKE_CHUNKS = 16; /**< Azimuth chunks clip-tested during the bake. */
  static constexpr uint32_t CHUNK_MASK = (1u << BAKE_CHUNKS) - 1; /**< All-chunks-visible bake mask. */
  static_assert(LUT_MIN_SAMPLES >= BAKE_CHUNKS,
                "the one-chunk visibility pad needs at least one LUT column "
                "per chunk");

  float color_spin = 0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  static constexpr int RING_SLOTS = 72; /**< Baked-ring pool capacity; matches the Rings slider max. */
  float *shift_pool = nullptr;  /**< RING_SLOTS x (W + 1) pooled shift LUTs, one slot per drawn ring; entry lut_n repeats entry 0 to close the polyline. */
  Pixel *hue_pool = nullptr;    /**< RING_SLOTS x (W + 1) pooled hue-rotated ring colors, aligned with shift_pool. */
  float *slot_frag_alpha = nullptr; /**< Per-slot fragment alpha (ring alpha x sprite fade x Alpha slider). */
  int *slot_lut_n = nullptr;        /**< Per-slot bake column count. */
  int8_t *slot_by_ring = nullptr;   /**< Ring index -> slot, -1 for culled rings; rebuilt per frame. */
  void *shapes_raw = nullptr;       /**< Raw storage for RING_SLOTS placement-built SDF::DistortedRing shapes. */
  Pixel *hue_table = nullptr; /**< HUE_TABLE_SIZE + 1 dynamic or cyclic hue samples for the current ring. */
  float *solid_colat = nullptr; /**< SOLID_MAX active-body center colatitudes about the stack axis (radians), rebuilt per frame. */
  float *solid_reach = nullptr; /**< SOLID_MAX active-body support extents (radians): the reach prefilter bound. */
  float *solid_shift = nullptr; /**< SOLID_MAX active-body shift bounds (radians): the per-ring band bound. */
  float *solid_scale = nullptr; /**< SOLID_MAX active-body LUT feature scales (2/radius). */
  int *solid_local = nullptr;   /**< SOLID_MAX scratch: active indices of the bodies that can reach the current ring. */
#ifdef HS_TEST_BUILD
  bool force_exact_hue = false;
  int hue_table_uses = 0;
#endif

  /** Ring-to-body prefilter pad absorbing fast_acos and tangent-recurrence
   *  rounding (radians); a body excluded despite the pad fields exactly 0
   *  everywhere on the ring. */
  static constexpr float BALL_TOUCH_EPS = 1e-3f;

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;       /**< Overall ring opacity multiplier in [0, 1]. */
    float num_rings = 48.0f;  /**< Number of evenly spaced rings (truncated to int when drawn). */
    float thickness = 0.035f;  /**< Stroke half-width (radians), set by init(). */
    float ball_amp = 0.1f;    /**< Ball drape strength; scaled by BALL_DRAPE_PER_AMPLITUDE into the drape gain. */
    float noise_amp = 0.2f;   /**< Peak polar displacement (radians) of the noise phase. */
    float scale1 = 1.5f;      /**< Spatial frequency of the envelope octave; its zero regions leave rings undisturbed. */
    float scale2 = 3.0f;      /**< Spatial frequency of the detail octave, scaled by the envelope octave. */
    float hue_scale = 2.0f;   /**< Hue rotation (turns) per radian of displacement magnitude. */
    float flow_speed = 0.03f; /**< Noise-field time advance per frame. */
    float ball_min = 0.15f;   /**< Smallest ball footprint (radians). */
    float ball_max = 0.3f;    /**< Largest ball footprint (radians). */
    float ball_rate = 20.0f;  /**< Ball spawns per second (jittered ±50%). */
    float ball_speed_min = 0.45f; /**< Slowest fall (pole-to-pole traversals per second). */
    float ball_speed_max = 0.85f; /**< Fastest fall (pole-to-pole traversals per second). */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(DisplacementField)
