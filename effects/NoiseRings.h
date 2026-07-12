/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Stack of evenly spaced soft-stroked rings displaced by a
 * displacement-field stack: falling ball bumps alternating with a 3D noise
 * field.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings share one axis and are spaced evenly in colatitude across the
 * whole sphere. Each ring vertex is displaced along
 * the stack axis by the summed displacement fields sampled at the vertex's
 * world-space position. Ball and noise phases alternate: cap-shaped bumps
 * spawn at the stack pole on random meridians and fall to the opposite pole at
 * varying speeds, bowing the rings away from each ball's center. Once the last
 * ball has fallen, the noise field fades in from zero amplitude — the product
 * of two OpenSimplex octaves (independent spatial scale per octave), where
 * octave 1 envelopes octave 2, so perturbations bunch where the envelope is
 * strong and vanish where it crosses zero — dwells at full strength, then
 * fades back out into the next ball phase. Noise displacement is coherent
 * across rings, drifts as the field animates, and its direction is uniform
 * across the whole sphere (not
 * mirrored per hemisphere). Fragments are shaded from a circular analogous
 * palette that spins across the stack, with hue rotated proportionally to the
 * local displacement magnitude; a ColorWipe slowly fades the palette to a
 * freshly generated one every few seconds. Orientation random-walks over time.
 */
template <int W, int H> class NoiseRings : public Effect {
public:
  /**
   * @brief Builds the effect with its palette and ring-stack axis.
   */
  HS_COLD_MEMBER NoiseRings()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        balls(timeline), noise_field(timeline), palette(make_palette()),
        next_palette(make_palette()), normal(X_AXIS) {}

  /**
   * @brief Registers params, seeds the noise field, and builds the timeline.
   */
  void init() override {
    ring_baked.bake(persistent_arena, palette);

    noise_field.template_params.noise.SetSeed(hs::rand_int(0, 65536));

    // One pixel of azimuth in ring-space.
    const float px = 2.0f * PI_F / W;
    params.thickness = 1.5f * px;

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
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, alternates the
   * ball and noise phases under the master-gain fade, advances the palette
   * wipe, then renders the timeline for one frame.
   * @details master_gain fades the noise phase in and out; the balls run at
   * full drape (each ball carries its own pole-edge envelope). The instant the
   * last ball is gone the noise fades straight in from zero — no stillness beat.
   */
  void draw_frame() override {
    color_spin = wrap_t(color_spin + COLOR_SPIN_RATE);
    step_wipe_rebake(wipe_pending, wipe_frames_remaining, ring_baked, palette);

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
        // Last ball gone: fade the noise straight in. master_gain gates only the
        // noise (never the balls), so it starts from 0 with no stillness beat.
        phase = Phase::NOISE;
        master_gain = 0.0f;
        noise_hold = NOISE_HOLD_FRAMES;
        timeline.add(0, Animation::Transition(master_gain, 1.0f,
                                              NOISE_FADE_FRAMES,
                                              ease_in_out_sin));
      }
      break;
    case Phase::NOISE:
      if (master_gain >= 1.0f && noise_hold > 0 && --noise_hold == 0)
        timeline.add(0, Animation::Transition(master_gain, 0.0f,
                                              NOISE_FADE_FRAMES,
                                              ease_in_out_sin));
      if (noise_hold == 0 && master_gain <= 0.0f) {
        phase = Phase::BALLS;
        ball_phase_left = BALL_PHASE_FRAMES;
        spawn_cooldown = 0;
      }
      break;
    }

    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /**
   * @brief Draws all rings for this frame.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details Per ring, the summed displacement stack is baked per azimuth
   * column into shift_lut (the ring's centerline shift) together with the
   * hue-rotated ring color in hue_lut, so fragments lerp a precomputed color
   * instead of building a hue rotation each. The LUT resolution is adaptive:
   * enough samples for the finest active feature (noise octave or ball
   * footprint) along the ring's actual circumference. Rings rasterize as soft
   * SDF strokes (Scan::DistortedRing) with a quintic cross-section falloff;
   * the scan confines work to each ring's own latitude band, and under a
   * partial clip (segmented drivers) rings whose displaced band cannot touch
   * the clip are skipped whole.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing the AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float bound = balls.field_bound() + noise_field.field_bound();
    const float reach = bound + params.thickness + pad;
    float max_scale = std::max(params.scale1, params.scale2);
    // 2/R: the clearance arcs vary over half a footprint.
    if (balls.active_count() > 0)
      max_scale = std::max(
          max_scale, 2.0f / std::min(params.ball_min, params.ball_max));

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      float theta = radius * (PI_F / 2.0f);
      if (try_cull &&
          !Plot::cap_may_touch_clip<H>(clip(), basis.v, theta + reach))
        continue;
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);

      // Nyquist-safe column count: LUT_SAMPLES_PER_UNIT samples per noise-space
      // unit along the ring's circumference in the finer octave.
      int lut_n = hs::clamp(
          static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * PI_F *
                                 max_scale * sin_t)),
          LUT_MIN_SAMPLES, W);

      Color4 ring_color =
          ring_baked.get(wrap_t((i + 0.5f) / n_rings + color_spin));

      // ring_bound: true max |shift| over this ring's LUT (lerp never exceeds
      // it); the global field bound would widen every ring's scan band to the
      // sum of all active bumps.
      float ring_bound = 0.0f;
      for (int x = 0; x < lut_n; ++x) {
        float angle = 2.0f * PI_F * x / lut_n;
        Vector p = (basis.v * cos_t) +
                   ((basis.u * cosf(angle)) + (basis.w * sinf(angle))) * sin_t;
        float s = balls.field_dominant(p) + noise_field.field(p);
        ring_bound = std::max(ring_bound, std::fabs(s));
        shift_lut[x] = s;
        hue_lut[x] =
            hue_rotate(ring_color, std::fabs(s) * params.hue_scale).color;
      }
      shift_lut[lut_n] = shift_lut[0];
      hue_lut[lut_n] = hue_lut[0];

      auto sample_lut = [this, lut_n](float t) {
        float x = wrap_t(t) * lut_n;
        int j = static_cast<int>(x);
        return shift_lut[j] + (x - j) * (shift_lut[j + 1] - shift_lut[j]);
      };

      float frag_alpha = ring_color.alpha * opacity * params.alpha;
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        float x = wrap_t(f.v0) * lut_n;
        int j = static_cast<int>(x);
        float norm_dist = hs::clamp(f.v1 / f.size, 0.0f, 1.0f);
        float falloff = quintic_kernel(1.0f - norm_dist);
        f.color = Color4(hue_lut[j].lerp16(hue_lut[j + 1], frac_to_q16(x - j)),
                         frag_alpha * falloff);
      };

      Scan::DistortedRing::draw<W, H>(filters, canvas, basis, radius,
                                      params.thickness, sample_lut, ring_bound,
                                      fragment_shader);
    }
  }

private:
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
  static constexpr int BALL_PHASE_FRAMES = 900; /**< Ball-phase spawning window (~15 s); balls keep coming the whole window. */
  static constexpr float BALL_RATE_FPS = 60.0f;    /**< Frame cadence assumed by the Ball Rate and Speed sliders' per-second units. */
  static constexpr float BALL_DRAPE_PER_AMPLITUDE = 4.0f; /**< Drape gain per Amplitude unit: the 0.25 default = gain 1, the 0.8 max saturates toward full clearance. */
  static constexpr int NOISE_FADE_FRAMES = 150; /**< Noise amplitude ramp on each phase handoff. */
  static constexpr int NOISE_HOLD_FRAMES = 600; /**< Full-noise dwell before fading out into the next ball phase. */

  BallDropTransformer<MAX_BALLS> balls;   /**< Falling-ball displacement fields. */
  NoiseProductTransformer<1> noise_field; /**< Two-octave noise displacement field. */

  GenerativePalette palette;      /**< Active palette (mutated by an in-flight ColorWipe). */
  GenerativePalette next_palette; /**< Target palette the current wipe fades toward. */
  BakedPalette ring_baked;        /**< LUT the shader samples; rebaked while a wipe runs. */
  Vector normal;
  Orientation<> orientation;

  int wipe_frames_remaining = 0; /**< Frames left to rebake ring_baked for the in-flight wipe. */
  bool wipe_pending = false;     /**< Wipe armed this frame; it first steps next frame. */

  /**
   * @brief Displacement-phase state: falling balls, or the noise field (fading
   * in, dwelling, and fading out under master_gain).
   */
  enum class Phase { BALLS, NOISE };

  Phase phase = Phase::BALLS; /**< Current displacement phase. */
  int ball_phase_left = BALL_PHASE_FRAMES; /**< Frames left in this ball phase's spawning window. */
  int spawn_cooldown = 0;  /**< Frames until the next ball spawn. */
  float master_gain = 0.0f; /**< Noise fade envelope in [0, 1]; gates the noise field, animated by Transitions. */
  int noise_hold = 0;      /**< Frames until the noise phase begins fading back out. */

  static constexpr int PALETTE_CYCLE_FRAMES = 180; /**< Palette rollover period (~3 s at the ~60 fps cadence). */
  static constexpr int PALETTE_WIPE_FRAMES = 168;  /**< Wipe duration; slightly under the cycle so a wipe is never still in flight when the next rollover fires. */
  static constexpr float COLOR_SPIN_RATE = 0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float LUT_SAMPLES_PER_UNIT = 4.0f; /**< Bake columns per feature-space unit of ring circumference. */
  static constexpr int LUT_MIN_SAMPLES = 16;          /**< Bake-column floor for tiny/low-scale rings. */

  float color_spin = 0.0f; /**< Palette offset across the stack (turns, [0,1)). */
  float shift_lut[W + 1] = {}; /**< Per-ring displacement per bake column; entry lut_n repeats entry 0 for seamless lerp. */
  Pixel hue_lut[W + 1] = {};   /**< Hue-rotated ring color per bake column, aligned with shift_lut. */

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;       /**< Overall ring opacity multiplier in [0, 1]. */
    float num_rings = 55.0f;  /**< Number of evenly spaced rings (truncated to int when drawn). */
    float thickness = 0.03f;  /**< Stroke half-width (radians); init() reseeds it to 1.5 px scaled by resolution. */
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
REGISTER_EFFECT(NoiseRings)
