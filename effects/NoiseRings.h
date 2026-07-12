/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Stack of evenly spaced plotted rings displaced by a displacement-field
 * stack: falling ball bumps, then a 3D noise field.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings share one axis and are spaced evenly in colatitude across the
 * whole sphere. Each ring vertex is displaced along the stack axis by the
 * summed displacement fields sampled at the vertex's world-space position. The
 * effect opens with a ball-drop phase: cap-shaped bumps spawn at the stack pole
 * on random meridians and fall to the opposite pole at varying speeds, reading
 * as balls pushing the rings up from beneath. Once the last ball has fallen,
 * the noise field fades in from zero amplitude: the product of two OpenSimplex
 * octaves (independent spatial scale per octave), where octave 1 envelopes
 * octave 2, so perturbations bunch where the envelope is strong and vanish
 * where it crosses zero. Displacement is coherent across rings, drifts as the
 * field animates, and its direction is uniform across the whole sphere (not
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

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Rings", &params.num_rings, 1.0f, 72.0f);
    register_param("Amplitude", &params.amplitude, 0.0f, 0.8f);
    register_param("Scale 1", &params.scale1, 0.5f, 4.0f);
    register_param("Scale 2", &params.scale2, 0.5f, 8.0f);
    register_param("Hue Rotate", &params.hue_scale, 0.0f, 3.0f);
    register_param("Flow Speed", &params.flow_speed, 0.0f, 0.15f);

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

    timeline.add(0, Animation::RandomTimer(
                        BALL_SPAWN_MIN_FRAMES, BALL_SPAWN_MAX_FRAMES,
                        [this](Canvas &) { this->spawn_ball(); }, true));
  }

  /**
   * @brief Refreshes the displacement stack from the sliders, hands off from
   * the ball phase to the noise phase, advances the palette wipe, then renders
   * the timeline for one frame.
   */
  void draw_frame() override {
    color_spin = wrap_t(color_spin + COLOR_SPIN_RATE);
    step_wipe_rebake(wipe_pending, wipe_frames_remaining, ring_baked, palette);

    balls.template_params.amplitude = params.amplitude;
    noise_field.template_params.amplitude = params.amplitude * noise_gain;
    noise_field.template_params.scale1 = params.scale1;
    noise_field.template_params.scale2 = params.scale2;
    noise_field.template_params.speed = params.flow_speed;
    balls.prepare_frame();
    noise_field.prepare_frame();

    if (!noise_armed && balls_spawned >= BALL_COUNT &&
        balls.active_count() == 0) {
      noise_armed = true;
      timeline.add(0, Animation::Transition(noise_gain, 1.0f,
                                            NOISE_FADE_FRAMES, ease_in_out_sin));
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
   * column into shift_lut (sampled in the antipode work frame so the fields
   * track the plotted vertices) together with the hue-rotated ring color in
   * hue_lut, so fragments lerp a precomputed color instead of building a hue
   * rotation each. The LUT resolution is adaptive: enough samples for the
   * finest active feature (noise octave or ball footprint) along the ring's
   * actual circumference. Under a partial clip
   * (segmented drivers), rings whose displaced band cannot touch the clip are
   * skipped whole, and surviving rings are cut into CULL_CHUNKS azimuth chunks
   * whose bounding caps gate both the bake and the vertex sampling — invisible
   * chunks bake nothing and emit nothing, and visible runs draw as open arcs
   * whose vertices are bit-identical to the closed ring's, so clipped output
   * matches full-frame output exactly. get_antipode mirrors the polar frame
   * past the equator (radius > 1); hemi_sign negates the far-side shift so a
   * noise crest displaces every ring the same world direction along the stack
   * axis.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);
    Basis basis = make_basis(orientation.get(), normal);
    const bool try_cull = !clip().is_full();
    // World-angle pad absorbing stroke width + AA splat past the clip's margin.
    const float pad = 3.0f * PI_F / H;
    const float reach = balls.field_bound() + noise_field.field_bound() + pad;
    float max_scale = std::max(params.scale1, params.scale2);
    if (balls.active_count() > 0)
      max_scale = std::max(max_scale, 1.0f / BALL_RADIUS_MIN);

    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      auto res = get_antipode(basis, radius);
      Basis work = res.first;
      float theta_eq = res.second * (PI_F / 2.0f);
      if (try_cull && !may_touch_clip(work.v, theta_eq + reach))
        continue;
      float cos_eq = cosf(theta_eq);
      float sin_eq = sinf(theta_eq);
      float hemi_sign = radius > 1.0f ? -1.0f : 1.0f;

      // Nyquist-safe column count: LUT_SAMPLES_PER_UNIT samples per noise-space
      // unit along the ring's circumference in the finer octave.
      int lut_n = hs::clamp(
          static_cast<int>(ceilf(LUT_SAMPLES_PER_UNIT * 2.0f * PI_F *
                                 max_scale * sin_eq)),
          LUT_MIN_SAMPLES, W);

      Color4 ring_color =
          ring_baked.get(wrap_t((i + 0.5f) / n_rings + color_spin));

      auto bake_col = [&](int x) {
        float angle = 2.0f * PI_F * x / lut_n;
        Vector p = (work.v * cos_eq) +
                   ((work.u * cosf(angle)) + (work.w * sinf(angle))) * sin_eq;
        float s = hemi_sign * (balls.field(p) + noise_field.field(p));
        shift_lut[x] = s;
        hue_lut[x] =
            hue_rotate(ring_color, std::fabs(s) * params.hue_scale).color;
      };

      auto sample_lut = [this, lut_n](float t) {
        float x = wrap_t(t) * lut_n;
        int j = static_cast<int>(x);
        return shift_lut[j] + (x - j) * (shift_lut[j + 1] - shift_lut[j]);
      };

      float frag_alpha = ring_color.alpha * opacity * params.alpha;
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        float x = wrap_t(f.v0) * lut_n;
        int j = static_cast<int>(x);
        f.color = Color4(hue_lut[j].lerp16(hue_lut[j + 1], frac_to_q16(x - j)),
                         frag_alpha);
      };

      uint32_t vis = CHUNK_FULL_MASK;
      if (try_cull) {
        vis = 0;
        float chunk_reach = (PI_F / CULL_CHUNKS) * sin_eq + reach;
        for (int c = 0; c < CULL_CHUNKS; ++c) {
          float a = (c + 0.5f) * (2.0f * PI_F / CULL_CHUNKS);
          Vector mid = (work.v * cos_eq) +
                       ((work.u * cosf(a)) + (work.w * sinf(a))) * sin_eq;
          if (may_touch_clip(mid, chunk_reach))
            vis |= 1u << c;
        }
        if (vis == 0)
          continue;
      }

      if (vis == CHUNK_FULL_MASK) {
        for (int x = 0; x < lut_n; ++x)
          bake_col(x);
        shift_lut[lut_n] = shift_lut[0];
        hue_lut[lut_n] = hue_lut[0];
        Plot::DistortedRing::draw<W, H>(filters, canvas, basis, radius,
                                        sample_lut, fragment_shader);
        continue;
      }

      // Walk visible chunk runs in ascending order — the same segment order
      // the closed draw rasterizes, which alpha blending's order dependence
      // requires wherever the ring's own taps overlap. A run touching chunk 23
      // simply ends at vertex W, so no arc crosses the seam.
      constexpr int V = W / CULL_CHUNKS;
      int c = 0;
      while (c < CULL_CHUNKS) {
        while (c < CULL_CHUNKS && !(vis & (1u << c)))
          ++c;
        if (c == CULL_CHUNKS)
          break;
        int c0 = c;
        while (c < CULL_CHUNKS && (vis & (1u << c)))
          ++c;
        int i0 = c0 * V;
        int i1 = c * V;
        int lo = static_cast<int>(floorf(static_cast<float>(i0) * lut_n / W));
        int hi =
            static_cast<int>(ceilf(static_cast<float>(i1) * lut_n / W)) + 1;
        for (int m = lo; m <= hi; ++m)
          bake_col(m % lut_n);
        shift_lut[lut_n] = shift_lut[0];
        hue_lut[lut_n] = hue_lut[0];
        Plot::DistortedRing::draw_arc<W, H>(filters, canvas, basis, radius,
                                            sample_lut, fragment_shader, i0,
                                            i1);
      }
    }
  }

private:
  /**
   * @brief Conservative test: can a spherical cap reach the current clip?
   * @param dir Cap center direction (unit vector). A displaced ring passes its
   * work-frame axis with half_angle = colatitude + displacement bound (the cap
   * of that radius contains the ring's band); a ring chunk passes its midpoint
   * with half_angle = chunk half-arc + displacement bound.
   * @param half_angle Cap angular radius including stroke/AA pad (radians).
   * @return False only when no fragment inside the cap can land in the clip's
   * render region; true is always safe.
   * @details Rows: the cap's polar range about the display's Y axis is
   * [beta - t2, beta + t2] (beta = center colatitude). Columns: a cap that
   * reaches either display pole spans all longitudes; otherwise its longitude
   * half-width about the center longitude is asin(sin t2 / sin beta), compared
   * against the clip's column wedge with the clip margin plus one pixel of
   * slack.
   */
  bool may_touch_clip(const Vector &dir, float half_angle) const {
    const ClipRegion &cr = clip();
    float t2 = std::min(half_angle, PI_F);
    float beta = acosf(hs::clamp(dir.y, -1.0f, 1.0f));

    float phi_lo = std::max(beta - t2, 0.0f);
    float phi_hi = std::min(beta + t2, PI_F);
    if (phi_to_y<H>(phi_hi) < cr.render_y_start() ||
        phi_to_y<H>(phi_lo) >= cr.render_y_end())
      return false;

    if (cr.x_start == 0 && cr.x_end == cr.w)
      return true;
    if (beta <= t2 || PI_F - beta <= t2)
      return true;
    float dlam = asinf(hs::clamp(sinf(t2) / sinf(beta), 0.0f, 1.0f));
    float lam_v = atan2f(dir.z, dir.x);
    float width_px = static_cast<float>(cr.x_end - cr.x_start);
    float half_w =
        (width_px * 0.5f + cr.margin + 1.0f) * (2.0f * PI_F) / cr.w;
    float lam_c = (cr.x_start + width_px * 0.5f) * (2.0f * PI_F) / cr.w;
    float d = std::fabs(wrap_t((lam_v - lam_c) / (2.0f * PI_F) + 0.5f) - 0.5f) *
              (2.0f * PI_F);
    return d <= dlam + half_w;
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
   * @brief Spawns one falling ball with a random meridian, footprint, and fall
   * time, until the phase's ball budget is spent.
   */
  void spawn_ball() {
    if (balls_spawned >= BALL_COUNT)
      return;
    balls.template_params.radius =
        hs::rand_f(BALL_RADIUS_MIN, BALL_RADIUS_MAX);
    if (balls.spawn(0, orientation, normal, hs::rand_f(0.0f, 2.0f * PI_F),
                    hs::rand_int(BALL_FALL_MIN_FRAMES, BALL_FALL_MAX_FRAMES)))
      ++balls_spawned;
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
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  static constexpr int MAX_BALLS = 8;   /**< Concurrent falling-ball pool slots. */
  static constexpr int BALL_COUNT = 12; /**< Total balls the intro phase spawns. */
  static constexpr int BALL_SPAWN_MIN_FRAMES = 20; /**< Shortest gap between spawns. */
  static constexpr int BALL_SPAWN_MAX_FRAMES = 55; /**< Longest gap between spawns. */
  static constexpr int BALL_FALL_MIN_FRAMES = 130; /**< Fastest pole-to-pole fall. */
  static constexpr int BALL_FALL_MAX_FRAMES = 280; /**< Slowest pole-to-pole fall. */
  static constexpr float BALL_RADIUS_MIN = 0.35f;  /**< Smallest bump footprint (radians). */
  static constexpr float BALL_RADIUS_MAX = 0.7f;   /**< Largest bump footprint (radians). */
  static constexpr int NOISE_FADE_FRAMES = 150; /**< Noise amplitude ramp after the last ball. */

  BallDropTransformer<MAX_BALLS> balls;   /**< Falling-ball displacement fields. */
  NoiseProductTransformer<1> noise_field; /**< Two-octave noise displacement field. */

  GenerativePalette palette;      /**< Active palette (mutated by an in-flight ColorWipe). */
  GenerativePalette next_palette; /**< Target palette the current wipe fades toward. */
  BakedPalette ring_baked;        /**< LUT the shader samples; rebaked while a wipe runs. */
  Vector normal;
  Orientation<> orientation;

  int wipe_frames_remaining = 0; /**< Frames left to rebake ring_baked for the in-flight wipe. */
  bool wipe_pending = false;     /**< Wipe armed this frame; it first steps next frame. */

  int balls_spawned = 0;   /**< Balls spawned so far, capped at BALL_COUNT. */
  bool noise_armed = false; /**< Noise fade-in scheduled (fires once, after the last ball). */
  float noise_gain = 0.0f; /**< Noise amplitude ramp in [0, 1], animated by a Transition. */

  static constexpr int PALETTE_CYCLE_FRAMES = 180; /**< Palette rollover period (~3 s at the ~60 fps cadence). */
  static constexpr int PALETTE_WIPE_FRAMES = 168;  /**< Wipe duration; slightly under the cycle so a wipe is never still in flight when the next rollover fires. */
  static constexpr float COLOR_SPIN_RATE = 0.0015f; /**< Palette spin across the stack, in turns per frame. */
  static constexpr float LUT_SAMPLES_PER_UNIT = 4.0f; /**< Bake columns per feature-space unit of ring circumference. */
  static constexpr int LUT_MIN_SAMPLES = 16;          /**< Bake-column floor for tiny/low-scale rings. */
  static constexpr int CULL_CHUNKS = 24; /**< Azimuth chunks per ring for clip culling; must divide W. */
  static constexpr uint32_t CHUNK_FULL_MASK = (1u << CULL_CHUNKS) - 1;
  static_assert(W % CULL_CHUNKS == 0,
                "chunk-to-vertex mapping requires W divisible by CULL_CHUNKS");

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
    float amplitude = 0.25f;  /**< Peak polar displacement (radians) scaling every displacement field. */
    float scale1 = 1.5f;      /**< Spatial frequency of the envelope octave; its zero regions leave rings undisturbed. */
    float scale2 = 3.0f;      /**< Spatial frequency of the detail octave, scaled by the envelope octave. */
    float hue_scale = 2.0f;   /**< Hue rotation (turns) per radian of displacement magnitude. */
    float flow_speed = 0.03f; /**< Noise-field time advance per frame. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(NoiseRings)
