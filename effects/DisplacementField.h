/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

namespace hs_test {
namespace poi_tests {
struct PoiWhiteBox;
} // namespace poi_tests
} // namespace hs_test

/**
 * @brief Stack of evenly spaced soft-stroked rings displaced by a
 * displacement-field stack that alternates falling ball bumps, a 3D noise
 * field, and a troupe of dancing poi domes.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings share one axis and are spaced evenly in colatitude across the
 * whole sphere. Each ring vertex is displaced along
 * the stack axis by the summed displacement fields sampled at the vertex's
 * world-space position. Three phases cycle BALLS -> NOISE -> POI -> NOISE ->
 * BALLS: cap-shaped bumps spawn at the stack pole on random meridians and fall
 * to the opposite pole, bowing the rings away from each ball's center (rings
 * part around each ball); the noise field is the product of two OpenSimplex
 * octaves (octave 1 envelopes octave 2) and fades in from zero between the two
 * solid-body phases, dwelling at full strength then fading back out; the poi
 * phase is 12 solid domes anchored at the icosahedron vertices in antipodal
 * pairs, tracing epicyclic tangent-plane curves (spirals, figure-8s, flowers)
 * under the sheet, each pushing the rings one way along the stack axis so the
 * sheet slides over the dancers like silk (a comet wake, never a parting).
 * Fragments are shaded from a circular analogous palette that spins across the
 * stack, with hue rotated proportionally to the local displacement magnitude; a
 * ColorWipe slowly fades the palette to a freshly generated one every few
 * seconds. Orientation random-walks over time.
 */
template <int W, int H> class DisplacementField : public Effect {
  friend struct ::hs_test::poi_tests::PoiWhiteBox;

public:
  /**
   * @brief Builds the effect with its palette and ring-stack axis.
   */
  HS_COLD_MEMBER DisplacementField()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        balls(timeline), noise_field(timeline), pois(timeline),
        palette(make_palette()), next_palette(make_palette()), normal(X_AXIS) {}

  /**
   * @brief Allocates the bake LUTs, registers params, seeds the noise field,
   * and builds the timeline.
   */
  void init() override {
    ring_baked.bake(persistent_arena, palette);
    shift_lut = persistent_arena.allocate_n<float>(W + 1);
    hue_lut = persistent_arena.allocate_n<Pixel>(W + 1);
    solid_colat = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_reach = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_shift = persistent_arena.allocate_n<float>(SOLID_MAX);
    solid_local = persistent_arena.allocate_n<int>(SOLID_MAX);
    balls.init_storage(persistent_arena);
    noise_field.init_storage(persistent_arena);
    pois.init_storage(persistent_arena);
    poi_choreo.set_speed_source(&params.dance_speed);

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
    register_param("Poi Amp", &params.poi_amp, 0.0f, 0.8f);
    register_param("Poi Size", &params.poi_size, 0.15f, 0.5f);
    register_param("Dance Speed", &params.dance_speed, 0.25f, 2.0f);

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
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, advances the
   * BALLS -> NOISE -> POI phase machine under the master-gain fade, steps the
   * poi choreography, advances the palette wipe, then renders one frame.
   * @details master_gain fades the noise phase in and out; the balls and pois
   * run at full drape (each carries its own lifecycle envelope). The instant the
   * last solid body is gone the noise fades straight in from zero — no stillness
   * beat. A next_is_poi flag, toggled at each NOISE entry, routes the NOISE exit
   * to the poi phase or back to the ball phase in turn.
   */
  void draw_frame() override {
    color_spin = wrap_t(color_spin + COLOR_SPIN_RATE);
    step_wipe_rebake(wipe_pending, wipe_frames_remaining, ring_baked, palette);

    if (phase == Phase::POI)
      poi_choreo.step();

    balls.template_params.amplitude =
        params.ball_amp * BALL_DRAPE_PER_AMPLITUDE;
    noise_field.template_params.amplitude =
        phase == Phase::NOISE ? params.noise_amp * master_gain : 0.0f;
    noise_field.template_params.scale1 = params.scale1;
    noise_field.template_params.scale2 = params.scale2;
    noise_field.template_params.speed = params.flow_speed;
    pois.template_params.amplitude = params.poi_amp * POI_DRAPE_PER_AMPLITUDE;
    balls.prepare_frame();
    noise_field.prepare_frame();
    pois.prepare_frame();

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
    case Phase::POI:
      // Dancers were all spawned at phase entry; wait for the last to deflate.
      if (pois.active_count() == 0)
        enter_noise();
      break;
    case Phase::NOISE:
      if (master_gain >= 1.0f && noise_hold > 0 && --noise_hold == 0)
        timeline.add(0, Animation::Transition(master_gain, 0.0f,
                                              NOISE_FADE_FRAMES,
                                              ease_in_out_sin));
      if (noise_hold == 0 && master_gain <= 0.0f) {
        if (next_is_poi)
          enter_poi();
        else
          enter_balls();
      }
      break;
    }

    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /**
   * @brief Draws all rings for this frame, over whichever solid-body pool is
   * active (balls or pois; the phases are temporally disjoint).
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details The two solid pools never overlap in time, so a single active pool
   * is chosen per frame and the shared per-ring machinery runs over it; the poi
   * feature scale follows the same 2/R rule as the balls.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    if (pois.active_count() > 0) {
      draw_rings(canvas, opacity, pois,
                 2.0f / (params.poi_size * (1.0f - POI_SIZE_JITTER)));
    } else {
      float ball_feature =
          balls.active_count() > 0
              ? 2.0f / std::min(params.ball_min, params.ball_max)
              : 0.0f;
      draw_rings(canvas, opacity, balls, ball_feature);
    }
  }

private:
  /**
   * @brief Bakes and rasterizes every ring over one solid-body pool.
   * @tparam Pool The active FieldTransformer pool (balls or pois).
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite fade multiplied into each fragment's alpha.
   * @param pool The active solid-body pool sampled by field_dominant.
   * @param solid_feature The pool's 2/R LUT feature scale (used when a body
   * touches the ring).
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
   * and slope-compensated width, so steeply displaced segments keep full
   * thickness; under a partial clip, rings whose displaced band cannot touch the
   * clip are skipped whole.
   */
  template <typename Pool>
  void draw_rings(Canvas &canvas, float opacity, const Pool &pool,
                  float solid_feature) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float noise_bound = noise_field.field_bound();

    // Body centers and rings share the stack axis, so a center's colatitude
    // difference from a ring's is the min ring-to-center distance: body b can
    // reach ring theta only when |theta - colat_b| < support_b.
    const int n_solid = pool.active_count();
    for (int b = 0; b < n_solid; ++b) {
      const auto &sp = pool.active_params(b);
      solid_colat[b] =
          fast_acos(hs::clamp(dot(basis.v, sp.center), -1.0f, 1.0f));
      solid_reach[b] = sp.support_bound();
      solid_shift[b] = sp.field_bound();
    }

    // The octave product carries modulation sidebands out to scale1 + scale2.
    const float noise_feature = params.scale1 + params.scale2;

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      float theta = radius * (PI_F / 2.0f);

      // Bodies that can reach this ring; the dominant blend never exceeds its
      // largest input, so the max touching field_bound bounds the shift.
      int n_local = 0;
      float band = 0.0f;
      for (int b = 0; b < n_solid; ++b) {
        if (std::fabs(theta - solid_colat[b]) <
            solid_reach[b] + BALL_TOUCH_EPS) {
          solid_local[n_local++] = b;
          band = std::max(band, solid_shift[b]);
        }
      }

      if (try_cull &&
          !Plot::cap_may_touch_clip<H>(clip(), basis.v,
                                       theta + band + noise_bound +
                                           params.thickness + pad))
        continue;

      Color4 ring_color =
          ring_baked.get(wrap_t((i + 0.5f) / n_rings + color_spin));
      HueRotateBase hue_base = make_hue_rotate_base(ring_color);

      float ring_bound = 0.0f;
      int lut_n;
      if (band + noise_bound <= 0.0f) {
        // Nothing can displace this ring: constant zero-shift LUT.
        lut_n = LUT_MIN_SAMPLES;
        Pixel flat = hue_rotate(hue_base, 0.0f).color;
        for (int x = 0; x <= lut_n; ++x) {
          shift_lut[x] = 0.0f;
          hue_lut[x] = flat;
        }
      } else {
        float feature_scale = noise_feature;
        if (n_local > 0)
          feature_scale = std::max(feature_scale, solid_feature);
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
        // Visibility is padded one chunk per side for boundary interpolation.
        uint32_t visible = CHUNK_MASK;
        if (try_cull) {
          const float chunk_reach = (PI_F / BAKE_CHUNKS) * sin_t + band +
                                    noise_bound + params.thickness + pad;
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
          visible = (raw | (raw << 1) | (raw >> (BAKE_CHUNKS - 1)) |
                     (raw >> 1) | (raw << (BAKE_CHUNKS - 1))) &
                    CHUNK_MASK;
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

        // ring_bound: true max |shift| over this ring's LUT (lerp never
        // exceeds it); the global field bound would widen every ring's scan
        // band to the sum of all active bodies.
        int x = 0;
        for (int c = 0; c < BAKE_CHUNKS; ++c) {
          const int x_end =
              ((c + 1) * lut_n + BAKE_CHUNKS - 1) / BAKE_CHUNKS;
          if (visible & (1u << c)) {
            for (; x < x_end; ++x) {
              Vector p = (basis.v * cos_t) +
                         ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t;
              float s = pool.field_dominant(p, solid_local, n_local) +
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

      // Quintic-remapped fractions keep the reconstruction smooth across LUT
      // knots and never exceed the knot values, so ring_bound stays a true max.
      // The slope out-param feeds the rasterizer's slope compensation: the
      // per-cell secant, not the quintic's instantaneous derivative, which dips
      // to zero at every knot and would pulse the compensated stroke width.
      auto sample_lut = [this, lut_n](float t, float &slope) {
        float x = wrap_t(t) * lut_n;
        int j = static_cast<int>(x);
        float d = shift_lut[j + 1] - shift_lut[j];
        slope = d * lut_n;
        return shift_lut[j] + quintic_kernel(x - j) * d;
      };

      float frag_alpha = ring_color.alpha * opacity * params.alpha;
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        float x = wrap_t(f.v0) * lut_n;
        int j = static_cast<int>(x);
        float norm_dist = hs::clamp(f.v1 / f.size, 0.0f, 1.0f);
        float falloff = quintic_kernel(1.0f - norm_dist);
        f.color = Color4(
            hue_lut[j].lerp16(hue_lut[j + 1], frac_to_q16(quintic_kernel(x - j))),
            frag_alpha * falloff);
      };

      Scan::DistortedRing::draw<W, H>(filters, canvas, basis, radius,
                                      params.thickness, sample_lut, ring_bound,
                                      fragment_shader);
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
   * @brief Enters the noise phase: the noise fades straight in from zero and the
   * next-solid-body selector toggles, so BALLS and POI alternate across the two
   * NOISE bridges.
   */
  void enter_noise() {
    phase = Phase::NOISE;
    next_is_poi = !next_is_poi;
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
   * @brief Enters the poi phase: draws a fresh anchor rotation and push sign,
   * reshuffles the move table, and spawns all 12 dancers with a per-pair
   * cascading entrance and footprint jitter.
   */
  void enter_poi() {
    phase = Phase::POI;
    Quaternion q_dance = random_rotation();
    Vector primaries[Animation::PoiChoreography::NUM_PAIRS];
    for (int p = 0; p < Animation::PoiChoreography::NUM_PAIRS; ++p)
      primaries[p] = Solids::Icosahedron::vertices[POI_PRIMARY[p]];
    poi_choreo.begin(q_dance, primaries);
    pois.template_params.direction = hs::rand_f() < 0.5f ? -1.0f : 1.0f;

    float jitter[Animation::PoiChoreography::NUM_PAIRS];
    for (int k = 0; k < N_POIS; ++k) {
      int p = k % Animation::PoiChoreography::NUM_PAIRS;
      if (k < Animation::PoiChoreography::NUM_PAIRS)
        jitter[p] =
            hs::rand_f(1.0f - POI_SIZE_JITTER, 1.0f + POI_SIZE_JITTER);
      pois.template_params.radius = params.poi_size * jitter[p];
      float canon[3];
      poi_choreo.compute_canon(p, params.dance_speed, canon);
      pois.spawn(p * PAIR_STAGGER, poi_choreo, orientation, normal, k, canon,
                 POI_PHASE_FRAMES);
    }
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

  FastNoiseLite walk_noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  // Near the timeline's 64-event budget: each in-flight ball is one event, so
  // at high Ball Rate x slow Speed the spawner saturates here and drops spawns
  // safely instead of overflowing the shared event buffer.
  static constexpr int MAX_BALLS = 56;  /**< Concurrent falling-ball pool slots. */
  static constexpr int N_POIS = 12;     /**< Poi dancer pool slots (6 antipodal pairs). */
  static constexpr int SOLID_MAX = MAX_BALLS > N_POIS ? MAX_BALLS : N_POIS; /**< Shared prefilter scratch size; the phase machine keeps at most one solid pool active. */
  static constexpr int BALL_PHASE_FRAMES = 900; /**< Ball-phase spawning window (~15 s); balls keep coming the whole window. */
  static constexpr int POI_PHASE_FRAMES = 1500; /**< Poi-phase dance duration (~25 s at the 60 fps cadence). */
  static constexpr float BALL_RATE_FPS = 60.0f;    /**< Frame cadence assumed by the Ball Rate and Speed sliders' per-second units. */
  static constexpr float BALL_DRAPE_PER_AMPLITUDE = 4.0f; /**< Drape gain per Ball Amp unit: the 0.25 default = gain 1, the 0.8 max saturates toward full clearance. */
  static constexpr float POI_DRAPE_PER_AMPLITUDE = 4.0f;  /**< Drape gain per Poi Amp unit, mirroring the ball drape. */
  static constexpr float POI_SIZE_JITTER = 0.15f; /**< Per-pair footprint jitter fraction of Poi Size. */
  static constexpr int PAIR_STAGGER = 8; /**< Spawn in_frames step per poi pair (cascading entrance). */
  static constexpr int NOISE_FADE_FRAMES = 150; /**< Noise amplitude ramp on each phase handoff. */
  static constexpr int NOISE_HOLD_FRAMES = 600; /**< Full-noise dwell before fading out into the next solid phase. */

  /** Icosahedron vertex per antipodal pair; the partner dancer negates the
   *  primary's center, so only these six are needed. */
  static constexpr int POI_PRIMARY[N_POIS / 2] = {0, 1, 4, 5, 8, 9};

  BallDropTransformer<MAX_BALLS> balls;   /**< Falling-ball displacement fields. */
  NoiseProductTransformer<1> noise_field; /**< Two-octave noise displacement field. */
  PoiTransformer<N_POIS> pois;            /**< Dancing-poi displacement fields. */
  Animation::PoiChoreography poi_choreo;  /**< Shared poi dance state, stepped per frame. */

  GenerativePalette palette;      /**< Active palette (mutated by an in-flight ColorWipe). */
  GenerativePalette next_palette; /**< Target palette the current wipe fades toward. */
  BakedPalette ring_baked;        /**< LUT the shader samples; rebaked while a wipe runs. */
  Vector normal;
  Orientation<> orientation;

  int wipe_frames_remaining = 0; /**< Frames left to rebake ring_baked for the in-flight wipe. */
  bool wipe_pending = false;     /**< Wipe armed this frame; it first steps next frame. */

  /**
   * @brief Displacement-phase state: falling balls, the noise field (fading in,
   * dwelling, fading out under master_gain), or the dancing poi troupe.
   */
  enum class Phase { BALLS, NOISE, POI };

  Phase phase = Phase::BALLS; /**< Current displacement phase. */
  bool next_is_poi = false;   /**< Whether the next solid phase after NOISE is POI. */
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
  float *shift_lut = nullptr; /**< W + 1 arena-baked displacements, one per bake column; entry lut_n repeats entry 0 for seamless lerp. */
  Pixel *hue_lut = nullptr;   /**< W + 1 arena-baked hue-rotated ring colors, aligned with shift_lut. */
  float *solid_colat = nullptr; /**< SOLID_MAX active-body center colatitudes about the stack axis (radians), rebuilt per frame. */
  float *solid_reach = nullptr; /**< SOLID_MAX active-body support extents (radians): the reach prefilter bound. */
  float *solid_shift = nullptr; /**< SOLID_MAX active-body shift bounds (radians): the per-ring band bound. */
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
    float poi_amp = 0.25f;    /**< Poi drape strength; scaled by POI_DRAPE_PER_AMPLITUDE into the drape gain. */
    float poi_size = 0.3f;    /**< Poi footprint radius (radians), ±15% pair jitter at spawn. */
    float dance_speed = 1.0f; /**< Poi base angular speed scale (live). */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(DisplacementField)
