/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"
#include <array>

/**
 * @brief Animated Hopf fibration effect.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @details A grid of S2 base points is lifted to S3 fibers, folded/twisted/
 * tumbled in 4D, stereographically projected to R3, and drawn as fading
 * palette-colored trail polylines.
 */
template <int W, int H> class HopfFibration : public Effect {
public:
  static constexpr int TRAIL_LEN = 40; /**< Number of points retained per fiber trail. */

  /**
   * @brief Constructs the effect and configures the trail pipeline.
   * @details Disables pixel persistence and seeds the trail pipeline with a
   * screen-space anti-alias filter sized to W x H.
   */
  FLASHMEM HopfFibration()
      : Effect(W, H), trail_pipeline(Filter::Screen::AntiAlias<W, H>()) {
    persist_pixels = false;
  }

  /**
   * @brief Initializes params, buffers, fibers, and timeline animations.
   * @details Registers tuning params, bakes the trail palette, allocates the
   * fiber and trail arrays from the persistent arena, seeds fibers, and wires
   * the ambient spin plus the phase-driver animations onto the timeline.
   */
  void init() override {
    registerParam("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    registerParam("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    registerParam("Folding", &params.folding, 0.0f, 2.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Bake palette for trail coloring
    baked_sunset.bake(persistent_arena, Palettes::richSunset);

    fibers = static_cast<Vector *>(persistent_arena.allocate(
        ACTUAL_FIBERS * sizeof(Vector), alignof(Vector)));

    // Allocate VectorTrails for each fiber
    trails = static_cast<Animation::VectorTrail<TRAIL_LEN> *>(
        persistent_arena.allocate(ACTUAL_FIBERS *
                                      sizeof(Animation::VectorTrail<TRAIL_LEN>),
                                  alignof(Animation::VectorTrail<TRAIL_LEN>)));
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      new (&trails[i]) Animation::VectorTrail<TRAIL_LEN>();
    }

    init_fibers();
    // Ambient Y-axis spin is intentionally left ungated (keeps running while
    // paused, per the setAnimationsPaused contract). The flow/tumble Drivers
    // ARE gated so "Pause Animation" actually freezes the fiber motion.
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_mid, true));
    // Bound Drivers: each pulls its tuning slider (× per-unit rate) every step,
    // so the flow/tumble speeds stay live without a per-frame set_speed re-sync.
    timeline.add(0, Animation::Driver(flow_offset, &params.flow_speed, FLOW_RATE,
                                      true, &anims_paused_));
    timeline.add(0, Animation::Driver(tumble_angle_x, &params.tumble_speed,
                                      TUMBLE_X_RATE, true, &anims_paused_));
    timeline.add(0, Animation::Driver(tumble_angle_y, &params.tumble_speed,
                                      TUMBLE_Y_RATE, true, &anims_paused_));
  }

  /**
   * @brief Reports whether the effect wants the background drawn.
   * @return Always false; trails are drawn over a cleared frame.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances one frame: steps animations, records and renders trails.
   * @details Steps the timeline, refreshes cached tumble/phase values, records
   * each fiber's projected world-space point into its trail, then renders all
   * trails.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    advance_tumble();

    // Project fibers and record trail positions (oriented world space)
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      Vector v = hopf_project(i);
      // Store in world space (unoriented) — orient at render time
      trails[i].record(v);
    }

    // Render trails as polylines (already oriented, skip Orient filter)
    render_trails(canvas);
  }

  /**
   * @brief User-tunable parameters controlling fiber motion and shading.
   */
  struct Params {
    float flow_speed = 10.0f;  /**< Flow phase advance rate, in tuning units. */
    float folding = 0.2f;      /**< Fold modulation depth, dimensionless. */
    float tumble_speed = 2.0f; /**< 4D tumble rate, in tuning units. */
    float twist = 4.0f;        /**< Twist applied to azimuth per fold, dimensionless. */
    float alpha = 1.0f;        /**< Global trail opacity multiplier in [0, 1]. */
  } params;

private:
  static constexpr int RINGS = 15;
  static constexpr int PER_RING = 14;
  static constexpr size_t ACTUAL_FIBERS = RINGS * PER_RING;

  // Stereographic-projection guard: the denominator is (1 + eps) - q3, so the
  // projection pole (q3 == 1) maps to a finite 1/eps instead of dividing by
  // zero. normalized_or() then absorbs the degenerate direction.
  static constexpr float STEREO_POLE_EPSILON = 0.001f;

  // Phases accumulate as wrapped fractions of their natural period ("turns")
  // and are scaled back to radians at the use site, so the arguments handed to
  // fast_sinf/cosf stay bounded. An unbounded accumulator grows float ULPs
  // without limit and visibly steps over multi-day installation runs; Flyby and
  // DreamBalls wrap their trig phases for exactly this reason.
  //
  // flow_offset and tumble_angle_y feed only full-angle terms, so their period
  // is 2pi. tumble_angle_x additionally feeds the half-angle fold_base term
  // (period 4pi); wrapping it over 4pi keeps BOTH the full-angle (cx/sx) and the
  // half-angle (fold_base) terms continuous across the wrap.
  static constexpr float FLOW_PERIOD     = 2 * PI_F;
  static constexpr float TUMBLE_X_PERIOD = 4 * PI_F;
  static constexpr float TUMBLE_Y_PERIOD = 2 * PI_F;

  // Per-unit driver rates, in turns of the corresponding period (the original
  // radians-per-unit-speed rate divided by that period). Multiplied by the
  // tuning param via the bound Drivers set up in init().
  static constexpr float FLOW_RATE     = (0.02f * 0.2f) / FLOW_PERIOD;
  static constexpr float TUMBLE_X_RATE = 0.003f / TUMBLE_X_PERIOD;
  static constexpr float TUMBLE_Y_RATE = 0.005f / TUMBLE_Y_PERIOD;

  // Wrapped phase accumulators in [0, 1) turns (driven in init()).
  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  Vector *fibers = nullptr;
  Animation::VectorTrail<TRAIL_LEN> *trails = nullptr;

  // Cached per-frame values: tumble rotation (cx/sx/cy/sy), the fold_base phase
  // offset, and the radian-scaled flow/tumble-y phases reused across all fibers.
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;
  float flow_rad = 0.0f, ty_rad = 0.0f;

  Orientation<> orientation;
  Timeline timeline;
  BakedPalette baked_sunset;

  // trail_pipeline applies AA only; trail points are oriented by hand in
  // render_trails before rasterizing, so no Orient filter is needed.
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> trail_pipeline;

  /**
   * @brief Seeds the base fibers as a spherical lattice of unit vectors.
   * @details Lays out a RINGS x PER_RING grid of unit vectors, evenly spaced in
   * polar bands and azimuth around each band.
   */
  void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float polar = PI_F * (i + 0.5f) / RINGS;
      for (int j = 0; j < PER_RING; ++j) {
        float azimuth = 2 * PI_F * j / PER_RING;
        fibers[idx++] = Vector(Spherical(azimuth, polar));
      }
    }
  }

  /**
   * @brief Refreshes the per-frame cache from the wrapped phase accumulators.
   * @details Recomputes the tumble rotation sines/cosines, the fold_base offset,
   * and the radian-scaled flow and tumble-y phases shared across all fibers.
   */
  void advance_tumble() {
    // Scale the wrapped [0,1)-turn phases back to radians at use; the arguments
    // are now bounded (ax < 4pi, ty_rad < 2pi) instead of growing without limit.
    const float ax = tumble_angle_x * TUMBLE_X_PERIOD;
    ty_rad = tumble_angle_y * TUMBLE_Y_PERIOD;
    flow_rad = flow_offset * FLOW_PERIOD;
    cx = fast_cosf(ax);
    sx = fast_sinf(ax);
    cy = fast_cosf(ty_rad);
    sy = fast_sinf(ty_rad);
    // fold_base needs sin(ax/2): ax = tumble_angle_x*4pi, so ax*0.5 spans [0,2pi)
    // continuously across the 4pi wrap.
    fold_base = fast_sinf(ax * 0.5f) * 0.5f;
  }

  /**
   * @brief Projects a base fiber into oriented world space.
   * @param i Fiber index in [0, ACTUAL_FIBERS).
   * @return World-space point after folding, twist, S3 lift, tumble, and
   * stereographic projection; falls back to (1, 0, 0) at a degenerate pole.
   * @details Pipeline: folding -> twist -> S3 -> tumble -> stereographic R3.
   */
  Vector hopf_project(size_t i) const {
    const Vector &base = fibers[i];
    Spherical sph(base);
    float theta = sph.phi; // polar angle (co-latitude)
    float phi = sph.theta; // azimuthal angle

    // Folding
    float eta = theta / 2.0f;
    eta += fast_sinf(phi * 2.0f + ty_rad + fold_base) * 0.1f *
           params.tumble_speed * params.folding;

    // Twist
    phi += eta * params.twist;

    // S3 point
    float phase = i * (PI_F / ACTUAL_FIBERS);
    float beta = flow_rad + phase;
    float cos_eta = fast_cosf(eta), sin_eta = fast_sinf(eta);
    float q0 = cos_eta * fast_cosf(phi + beta);
    float q1 = cos_eta * fast_sinf(phi + beta);
    float q2 = sin_eta * fast_cosf(beta);
    float q3 = sin_eta * fast_sinf(beta);

    // Tumble (R_xw, R_yz)
    float q0_r = q0 * cx - q3 * sx;
    q3 = q0 * sx + q3 * cx;
    q0 = q0_r;
    float q1_r = q1 * cy - q2 * sy;
    q2 = q1 * sy + q2 * cy;
    q1 = q1_r;

    // Stereographic S3 → R3. At a fiber pole (q0=q1=q2=0) the projected
    // direction is undefined; fall back to a stable axis instead of trapping.
    float factor = 1.0f / ((1.0f + STEREO_POLE_EPSILON) - q3);
    return normalized_or(Vector(q0 * factor, q1 * factor, q2 * factor),
                         Vector(1, 0, 0));
  }

  /**
   * @brief Renders all fiber trails as anti-aliased polylines.
   * @param canvas Target canvas to rasterize the trail polylines onto.
   * @details Orients each trail's stored world-space points into the current
   * view and shades them with the sunset palette so the tail fades from newest
   * (opaque) to oldest (transparent).
   */
  void render_trails(Canvas &canvas) {
    // Every fragment's alpha is scaled by params.alpha (see the shader below),
    // so alpha == 0 paints nothing. Skip rasterizing all 210 fibers in that
    // case — output-identical, just cheaper. Trails are still recorded in
    // draw_frame(), so motion continues and resumes correctly when alpha > 0.
    if (params.alpha <= 0.0f)
      return;
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      const auto &trail = trails[i];
      size_t len = trail.length();
      if (len < 2)
        continue;

      ScratchScope _sc(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, len);

      for (size_t j = 0; j < len; ++j) {
        Fragment f;
        // Orient from world space to current view at render time
        f.pos = orientation.orient(trail.get(j));
        f.v0 = (len > 1) ? static_cast<float>(j) / (len - 1) : 0.0f;
        f.age = 0;
        points.push_back(f);
      }

      // Shade each trail point by its normalized age: sunset palette lookup
      // plus an alpha ramp that fades the tail from opaque (newest) to clear.
      auto shader = [this](const Vector &, Fragment &f) {
        float t = f.v0; // 0 = oldest, 1 = newest
        Color4 c = baked_sunset.get(1.0f - t);
        c.alpha *= t * params.alpha; // fade out tail
        f.color = c;
      };

      // Rasterize polyline through trail_pipeline (AA only, no Orient)
      Plot::rasterize<W, H>(trail_pipeline, canvas, points, shader);
    }
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(HopfFibration)
