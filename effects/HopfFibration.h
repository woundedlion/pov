/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"
#include <array>

// Unit-test accessor reaching the private per-frame cache and hopf_project() to
// pin the S3-lift + stereographic projection math directly.
namespace hs_test {
namespace effects_tests {
struct HopfWhiteBox;
} // namespace effects_tests
} // namespace hs_test

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
  HS_COLD_MEMBER HopfFibration()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(trail_pipeline)::any_crosses_segments}),
        trail_pipeline(Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Initializes params, buffers, fibers, and timeline animations.
   * @details Registers tuning params, bakes the trail palette, allocates the
   * fiber, phase, and trail arrays from the persistent arena, seeds fibers, and
   * wires the ambient spin plus the phase-driver animations onto the timeline.
   */
  void init() override {
    register_param("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    register_param("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    register_param("Folding", &params.folding, 0.0f, 2.0f);
    register_param("Twist", &params.twist, -5.0f, 5.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    baked_sunset.bake(persistent_arena, Palettes::RICH_SUNSET);

    // Whole footprint lives in the persistent arena with no per-frame scratch, so
    // keep the default split (no configure_arenas()); FOOTPRINT_BYTES asserts it
    // fits the device default partition.
    fibers = static_cast<Spherical *>(persistent_arena.allocate(
        ACTUAL_FIBERS * sizeof(Spherical), alignof(Spherical)));

    trails = static_cast<Animation::VectorTrail<TRAIL_LEN> *>(
        persistent_arena.allocate(ACTUAL_FIBERS *
                                      sizeof(Animation::VectorTrail<TRAIL_LEN>),
                                  alignof(Animation::VectorTrail<TRAIL_LEN>)));
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      new (&trails[i]) Animation::VectorTrail<TRAIL_LEN>();
    }

    init_fibers();
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_linear, true));
    timeline.add(0, Animation::Driver(flow_offset, &params.flow_speed, FLOW_RATE,
                                      true));
    timeline.add(0, Animation::Driver(tumble_angle_x, &params.tumble_speed,
                                      TUMBLE_X_RATE, true));
    timeline.add(0, Animation::Driver(tumble_angle_y, &params.tumble_speed,
                                      TUMBLE_Y_RATE, true));
  }

  /**
   * @brief Advances one frame: steps animations, records and renders trails.
   * @details Steps the timeline, refreshes cached tumble/phase values, records
   * each fiber's projected world-space point into its trail, then renders all
   * trails.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(hf_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(hf_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(hf_advance_tumble);
      advance_tumble();
    }

    {
      HS_PROFILE(hf_project_record);
      for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
        Vector v = hopf_project(i);
        // Store unoriented — oriented at render time.
        trails[i].record(v);
      }
    }

    render_trails(canvas);
  }

private:
  friend struct ::hs_test::effects_tests::HopfWhiteBox;

  static constexpr int RINGS = 15;
  static constexpr int PER_RING = 14;
  static constexpr size_t ACTUAL_FIBERS = RINGS * PER_RING;
  static constexpr float PHASE_STEP = PI_F / ACTUAL_FIBERS;

  // Persistent allocations: palette LUT + one Spherical and one trail per fiber.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::LUT_SIZE * sizeof(Color4) +
      ACTUAL_FIBERS * sizeof(Spherical) +
      ACTUAL_FIBERS * sizeof(Animation::VectorTrail<TRAIL_LEN>);
  // Effect keeps the default arena split, so the footprint must fit the device
  // persistent partition. Guards a RINGS/PER_RING/TRAIL_LEN retune.
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HopfFibration persistent footprint exceeds the default "
                "partition; retune RINGS/PER_RING/TRAIL_LEN or carve arenas");

  // render_trails stages one fiber's points (up to TRAIL_LEN) in scratch_a at a
  // time, so the trail length is bounded by the default scratch budget, not by
  // FOOTPRINT_BYTES; a retune past it should fail to compile, not overflow.
  static_assert(TRAIL_LEN * sizeof(Fragment) <= DEFAULT_SCRATCH_A_SIZE,
                "HopfFibration trail staging exceeds the default scratch_a "
                "budget; retune TRAIL_LEN or carve a larger scratch arena");

  // The trail gate holds every trail's per-edge visibility bytes at once.
  static_assert(ACTUAL_FIBERS * (TRAIL_LEN - 1) <= DEFAULT_SCRATCH_B_SIZE,
                "HopfFibration trail-gate bytes exceed the default scratch_b "
                "budget; retune RINGS/PER_RING/TRAIL_LEN");

  // Keeps the projection denominator (1 + eps) - q3 positive at the pole
  // (q3 == 1); only needed to keep the pre-normalize direction NaN-free, since
  // normalized_or() substitutes its fallback axis there.
  static constexpr float STEREO_POLE_EPSILON = 0.001f;

  // One 8-bit LSB at full-scale color: a shader alpha below this cannot move a
  // pixel, so such trail segments are skipped.
  static constexpr float MIN_VISIBLE_ALPHA = 1.0f / 255.0f;

  // Phases accumulate as wrapped fractions of their period ("turns") and scale
  // back to radians at use, keeping the trig arguments bounded. tumble_angle_x
  // also feeds the half-angle fold_base term, so it wraps over 4pi to keep both
  // the full-angle and half-angle terms continuous.
  static constexpr float FLOW_PERIOD     = 2 * PI_F;
  static constexpr float TUMBLE_X_PERIOD = 4 * PI_F;
  static constexpr float TUMBLE_Y_PERIOD = 2 * PI_F;

  // Per-unit driver rates, in turns (radians-per-unit-speed / period).
  static constexpr float FLOW_RATE     = (0.02f * 0.2f) / FLOW_PERIOD;
  static constexpr float TUMBLE_X_RATE = 0.003f / TUMBLE_X_PERIOD;
  static constexpr float TUMBLE_Y_RATE = 0.005f / TUMBLE_Y_PERIOD;

  // Wrapped phase accumulators in [0, 1) turns (driven in init()).
  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  Spherical *fibers = nullptr;
  Animation::VectorTrail<TRAIL_LEN> *trails = nullptr;

  // Per-frame values reused across all fibers (refreshed in advance_tumble).
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;
  float flow_rad = 0.0f, ty_rad = 0.0f;

  Orientation<> orientation;
  Timeline timeline;
  BakedPalette baked_sunset;

  // AA only; trail points are oriented by hand in render_trails, so no Orient
  // filter is needed.
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> trail_pipeline;

  /**
   * @brief Seeds the base fibers as a spherical lattice.
   * @details Lays out a RINGS x PER_RING grid of spherical coordinates, evenly
   * spaced in polar bands and azimuth around each band. Stored as Spherical;
   * hopf_project derives each fiber's S3 phase inline.
   */
  void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float polar = PI_F * (i + 0.5f) / RINGS;
      for (int j = 0; j < PER_RING; ++j) {
        float azimuth = 2 * PI_F * j / PER_RING;
        std::construct_at(&fibers[idx], azimuth, polar);
        ++idx;
      }
    }
  }

  /**
   * @brief Refreshes the per-frame cache from the wrapped phase accumulators.
   * @details Recomputes the tumble rotation sines/cosines, the fold_base offset,
   * and the radian-scaled flow and tumble-y phases shared across all fibers.
   */
  void advance_tumble() {
    // Scale wrapped [0,1)-turn phases back to radians (ax < 4pi, ty_rad < 2pi).
    const float ax = tumble_angle_x * TUMBLE_X_PERIOD;
    ty_rad = tumble_angle_y * TUMBLE_Y_PERIOD;
    flow_rad = flow_offset * FLOW_PERIOD;
    cx = fast_cosf(ax);
    sx = fast_sinf(ax);
    cy = fast_cosf(ty_rad);
    sy = fast_sinf(ty_rad);
    // fold_base needs sin(ax/2); ax*0.5 spans [0,2pi) continuously across the
    // 4pi wrap.
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
    const Spherical &sph = fibers[i];
    float polar = sph.phi;
    float azimuth = sph.theta;

    // Folding: amplitude gated by the Folding slider so it persists when tumble
    // is frozen, though its phase still tracks ty_rad.
    float eta = polar / 2.0f;
    eta += fast_sinf(azimuth * 2.0f + ty_rad + fold_base) * 0.2f * params.folding;

    // Twist
    azimuth += eta * params.twist;

    // S3 point
    float beta = flow_rad + static_cast<float>(i) * PHASE_STEP;
    float cos_eta = fast_cosf(eta), sin_eta = fast_sinf(eta);
    float q0 = cos_eta * fast_cosf(azimuth + beta);
    float q1 = cos_eta * fast_sinf(azimuth + beta);
    float q2 = sin_eta * fast_cosf(beta);
    float q3 = sin_eta * fast_sinf(beta);

    // Tumble (R_xw, R_yz)
    float q0_r = q0 * cx - q3 * sx;
    q3 = q0 * sx + q3 * cx;
    q0 = q0_r;
    float q1_r = q1 * cy - q2 * sy;
    q2 = q1 * sy + q2 * cy;
    q1 = q1_r;

    // Stereographic S3 -> R3; at a fiber pole the direction is undefined, so
    // fall back to a stable axis.
    float factor = 1.0f / ((1.0f + STEREO_POLE_EPSILON) - q3);
    return normalized_or(Vector(q0 * factor, q1 * factor, q2 * factor),
                         Vector(1, 0, 0));
  }

  /**
   * @brief Gates every trail's edges against the segment clip in one pass.
   * @param canvas Target canvas (supplies the active clip band).
   * @param bits Output: per-edge visibility bytes (gate_trail_edges),
   *        TRAIL_LEN - 1 per trail, allocated from scratch_arena_b (caller
   *        owns the scope).
   * @param skip Output, one entry per fiber: true when no edge is
   *        clip-visible, so the trail's stage + rasterize can be skipped.
   */
  HS_O3_FN void gate_trails(Canvas &canvas, uint8_t *&bits, bool *skip) {
    HS_PROFILE(hf_trail_gate);
    constexpr size_t STRIDE = TRAIL_LEN - 1;
    const ClipRegion &cr = canvas.clip();
    const auto xc = cr.x_clip();
    const int band_len = xc.wrap ? xc.re - xc.rs + W : xc.re - xc.rs;
    bits = static_cast<uint8_t *>(
        scratch_arena_b.allocate(ACTUAL_FIBERS * STRIDE, alignof(uint8_t)));

    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      skip[i] = false;
      const auto &trail = trails[i];
      const size_t len = trail.length();
      if (len < 2)
        continue;
      ScratchScope sc_guard(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, len);
      for (size_t j = 0; j < len; ++j) {
        Fragment f;
        f.pos = orientation.orient(trail.get(j));
        points.push_back(f);
      }
      skip[i] = !Plot::gate_trail_edges<W, H>(trail_pipeline, cr, xc, band_len,
                                              points, &bits[i * STRIDE]);
    }
  }

  /**
   * @brief Renders all fiber trails as anti-aliased polylines.
   * @param canvas Target canvas to rasterize the trail polylines onto.
   * @details Orients each trail's stored world-space points into the current
   * view and shades them with the sunset palette so the tail fades from newest
   * (opaque) to oldest (transparent). Under an active segment clip the trail
   * gate runs first: fully-invisible trails are skipped whole and the per-edge
   * bits feed rasterize as its cull.
   */
  HS_O3_FN void render_trails(Canvas &canvas) {
    // A sub-LSB alpha paints nothing; skip rasterizing. Trails are still
    // recorded in draw_frame(), so motion resumes when alpha rises.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;
    HS_PROFILE(hf_render_trails);

    constexpr size_t STRIDE = TRAIL_LEN - 1;
    const bool clip_active = !canvas.clip().is_full();
    ScratchScope gate_guard(scratch_arena_b);
    uint8_t *bits = nullptr;
    bool skip[ACTUAL_FIBERS] = {};
    if (clip_active)
      gate_trails(canvas, bits, skip);

    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      const auto &trail = trails[i];
      size_t len = trail.length();
      if (len < 2)
        continue;
      if (clip_active && skip[i])
        continue;

      // The tail fades to transparent; drop leading points whose outgoing
      // segment peaks (at its newer endpoint) below visibility. v0 keeps the
      // full-trail parameterization, so kept segments render unchanged.
      size_t first = 0;
      while (first + 1 < len &&
             static_cast<float>(first + 1) / (len - 1) * params.alpha <
                 MIN_VISIBLE_ALPHA)
        ++first;

      ScratchScope sc_guard(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, len - first);

      for (size_t j = first; j < len; ++j) {
        Fragment f;
        f.pos = orientation.orient(trail.get(j));
        f.v0 = static_cast<float>(j) / (len - 1);
        f.age = 0;
        points.push_back(f);
      }

      auto shader = [this](const Vector &, Fragment &f) {
        float t = f.v0; // 0 = oldest, 1 = newest
        Color4 c = baked_sunset.get(1.0f - t);
        c.alpha *= t * params.alpha; // fade out tail
        f.color = c;
      };

      {
        HS_PROFILE(hf_trail_raster);
        Plot::rasterize<W, H>(trail_pipeline, canvas, points, shader, false,
                              nullptr, false,
                              clip_active ? &bits[i * STRIDE + first]
                                          : nullptr);
      }
    }
  }

  /**
   * @brief User-tunable parameters controlling fiber motion and shading.
   */
  struct Params {
    float flow_speed = 10.0f;  /**< Flow phase advance rate, in tuning units. */
    float tumble_speed = 2.0f; /**< 4D tumble rate, in tuning units. */
    float folding = 0.2f;      /**< Fold modulation depth, dimensionless. */
    float twist = 4.0f;        /**< Twist applied to azimuth per fold, dimensionless. */
    float alpha = 1.0f;        /**< Global trail opacity multiplier in [0, 1]. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(HopfFibration)
