/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

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
    shift_lut = persistent_arena.allocate_n<float>(W + 1);
    hue_lut = persistent_arena.allocate_n<Pixel>(W + 1);
    solid_colat = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_reach = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_shift = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_scale = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_local = persistent_arena.allocate_n<int>(SOLID_MAX);
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
    balls.prepare_frame();
    noise_field.prepare_frame();

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

    Canvas canvas(*this);
    timeline.step(canvas);
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
   * @brief Bakes and rasterizes every ring over one solid-body pool.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite fade multiplied into each fragment's alpha.
   * @details Per ring, the displacement stack is baked per azimuth column into
   * shift_lut (the ring's centerline shift) together with the hue-rotated ring
   * color in hue_lut. The bake evaluates only the bodies whose support can reach
   * the ring's colatitude (centers and rings share the stack axis); a ring
   * nothing can displace skips the field/hue bake for a constant LUT. The reach
   * prefilter uses each body's support_bound() (cap extent) while the per-ring
   * band widening uses its field_bound() (shift bound) — for balls the two
   * coincide. The LUT resolution is adaptive: enough samples for the finest
   * active feature along the ring's actual circumference. Rings rasterize as
   * soft SDF strokes (Scan::DistortedRing) with a quintic cross-section falloff
   * against the exact distance to the baked knot polyline, so steeply displaced
   * or sharply curved segments keep full thickness; under a partial clip, rings
   * whose displaced band cannot touch the clip are skipped whole.
   */
  void draw_rings(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float noise_bound = noise_field.field_bound();

    // Body centers and rings share the stack axis, so a center's colatitude
    // difference from a ring's is the min ring-to-center distance: body b can
    // reach ring theta only when |theta - colat_b| < support_b.
    const int n_balls = balls.active_count();
    for (int b = 0; b < n_balls; ++b) {
      const auto &sp = balls.active_params(b);
      solid_colat[b] =
          fast_acos(hs::clamp(dot(basis.v, sp.center), -1.0f, 1.0f));
      solid_reach[b] = sp.field_bound();
      solid_shift[b] = sp.field_bound();
      solid_scale[b] = 2.0f / sp.radius;
    }

    // The octave product carries modulation sidebands out to scale1 + scale2.
    const float noise_feature =
        noise_bound > 0.0f ? params.scale1 + params.scale2 : 0.0f;

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      float theta = radius * (PI_F / 2.0f);

      // Bodies that can reach this ring; the dominant blend never exceeds its
      // largest input, so the max touching field_bound bounds the shift.
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

      float ring_bound = 0.0f;
      int lut_n;
      if (band + noise_bound <= 0.0f) {
        Pixel flat = hue_rotate(hue_base, 0.0f).color;
#ifdef __OPTIMIZE_SIZE__
        // The explicit SDF wins under -Os; -O3 folds the zero-knot path.
        float frag_alpha = ring_color.alpha * opacity * params.alpha;
        auto flat_shader = [flat, frag_alpha](const Vector &, Fragment &f) {
          f.color = Color4(flat, frag_alpha * f.v2);
        };
        Scan::DistortedRing::draw_flat<W, H, false>(
            filters, canvas, basis, radius, params.thickness, flat_shader);
        continue;
#else
        lut_n = LUT_MIN_SAMPLES;
        for (int x = 0; x <= lut_n; ++x) {
          shift_lut[x] = 0.0f;
          hue_lut[x] = flat;
        }
#endif
      } else {
        float feature_scale = std::max(noise_feature, solid_feature);
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);

        // Nyquist-safe column count: LUT_SAMPLES_PER_UNIT samples per
        // noise-space unit along the ring's circumference.
        lut_n = hs::clamp(
            static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * PI_F *
                                   feature_scale * sin_t)),
            LUT_MIN_SAMPLES, W);

        // Azimuth chunk cull under a partial clip: a chunk's cap (its ring
        // arc widened by the displaced band) contains every pixel whose
        // ring-frame azimuth falls in the chunk, so a chunk whose cap misses
        // the clip is never sampled and its columns skip the field/hue bake.
        // Visibility is padded by the polyline search's reach: a visible
        // pixel reads knots up to `thickness` of arc to each side, and the
        // band's pixels can sit closer to a pole than the ring itself, where
        // that arc spans more chunks; +1 chunk covers the hue lerp.
        uint32_t visible = CHUNK_MASK;
        if (try_cull) {
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

        // Angle-addition recurrence advancing the tangent (cos, sin) by one
        // column step, dropping a libm cosf/sinf per bake column. It advances
        // through skipped columns too, so baked values are independent of the
        // visible set.
        const float dphi = 2.0f * PI_F / lut_n;
        const float cos_d = cosf(dphi);
        const float sin_d = sinf(dphi);
        float cos_a = 1.0f;
        float sin_a = 0.0f;

        // ring_bound: true max |shift| over this ring's LUT (the polyline
        // never exceeds its knots); the global field bound would widen every
        // ring's scan band to the sum of all active bodies.
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
              shift_lut[x] = s;
              hue_lut[x] =
                  hue_rotate(hue_base, std::fabs(s) * params.hue_scale).color;
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
        shift_lut[lut_n] = shift_lut[0];
        hue_lut[lut_n] = hue_lut[0];
      }

      float frag_alpha = ring_color.alpha * opacity * params.alpha;
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        float x = wrap_t(f.v0) * lut_n;
        int j = static_cast<int>(x);
        f.color = Color4(
            hue_lut[j].lerp16(hue_lut[j + 1], frac_to_q16(quintic_kernel(x - j))),
            frag_alpha * f.v2);
      };

      Scan::DistortedRing::draw<W, H>(filters, canvas, basis, radius,
                                      params.thickness, shift_lut, lut_n,
                                      ring_bound, fragment_shader);
    }
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
  void enter_noise() {
    phase = Phase::NOISE;
    master_gain = 0.0f;
    noise_hold = NOISE_HOLD_FRAMES;
    timeline.add(0, Animation::Transition(master_gain, 1.0f, NOISE_FADE_FRAMES,
                                          ease_in_out_sin));
  }

  /**
   * @brief Enters the ball phase, opening a fresh spawning window.
   */
  void enter_balls() {
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
  static constexpr int BAKE_CHUNKS = 16; /**< Azimuth chunks clip-tested during the bake. */
  static constexpr uint32_t CHUNK_MASK = (1u << BAKE_CHUNKS) - 1; /**< All-chunks-visible bake mask. */
  static_assert(LUT_MIN_SAMPLES >= BAKE_CHUNKS,
                "the one-chunk visibility pad needs at least one LUT column "
                "per chunk");

  float color_spin = 0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  float *shift_lut = nullptr; /**< W + 1 arena-baked displacement knots, one per bake column; entry lut_n repeats entry 0 to close the polyline. */
  Pixel *hue_lut = nullptr;   /**< W + 1 arena-baked hue-rotated ring colors, aligned with shift_lut. */
  float *solid_colat = nullptr; /**< SOLID_MAX active-body center colatitudes about the stack axis (radians), rebuilt per frame. */
  float *solid_reach = nullptr; /**< SOLID_MAX active-body support extents (radians): the reach prefilter bound. */
  float *solid_shift = nullptr; /**< SOLID_MAX active-body shift bounds (radians): the per-ring band bound. */
  float *solid_scale = nullptr; /**< SOLID_MAX active-body LUT feature scales (2/radius). */
  int *solid_local = nullptr;   /**< SOLID_MAX scratch: active indices of the bodies that can reach the current ring. */

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
    float thickness = 0.04f;  /**< Stroke half-width (radians), set by init(). */
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
