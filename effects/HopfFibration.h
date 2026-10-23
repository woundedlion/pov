/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file HopfFibration.h
 * @brief Animated Hopf fibration: S2 base points lifted to S3 fibers, tumbled
 *        in 4D and stereographically projected.
 */

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include <span>

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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::HopfFibration::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  static constexpr int TRAIL_LEN =
      40; /**< Number of points retained per fiber trail. */

  /**
   * @brief Constructs the effect and configures the trail pipeline.
   */
  HS_COLD_MEMBER HopfFibration()
      : Effect(W, H,
               pipeline_config<decltype(trail_pipeline)>({.strobe = true})),
        trail_pipeline(Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Initializes params, buffers, fibers, and timeline animations.
   */
  HS_COLD_MEMBER void init() override {
    register_param("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    register_param("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    register_param("Folding", &params.folding, 0.0f, 2.0f);
    register_param("Twist", &params.twist, -5.0f, 5.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    baked_sunset.bake(persistent_arena, Palettes::RICH_SUNSET);
    timed_pipeline.init_storage(persistent_arena);

    fibers = persistent_arena.allocate_n<math::Spherical>(ACTUAL_FIBERS);

    trails = persistent_arena.make_n<Animation::VectorTrail<TRAIL_LEN>>(
        ACTUAL_FIBERS);

    retired = persistent_arena.allocate_n<math::Vector>(ACTUAL_FIBERS);
    init_fibers();
    timeline.add(
        0, Animation::Rotation<W>(orientation, math::Y_AXIS, 2 * math::PI_F,
                                  600, math::ease_linear, {.repeat = true}));
    timeline.add(0,
                 Animation::Driver(flow_offset, &params.flow_speed, FLOW_RATE));
    timeline.add(0, Animation::Driver(tumble_angle_x, &params.tumble_speed,
                                      TUMBLE_X_RATE));
    timeline.add(0, Animation::Driver(tumble_angle_y, &params.tumble_speed,
                                      TUMBLE_Y_RATE));
  }

  /**
   * @brief Advances one frame: steps animations, records and renders trails.
   */
  HS_HOT_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(hf_timeline_step);
      timeline.step(canvas);
      if (!history_ready) {
        if (!Spacetime::flat_surface())
          orientation.collapse();
        history_ready = true;
      }
    }
    {
      HS_PROFILE(hf_advance_tumble);
      advance_tumble();
    }

    {
      HS_PROFILE(hf_project_record);
      trail_was_full = trails[0].length() == TRAIL_LEN;
      for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
        if (trail_was_full)
          retired[i] = trails[i].get(0);
        math::Vector v = hopf_project(i);
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
  static constexpr float PHASE_STEP = math::PI_F / ACTUAL_FIBERS;

  // Persistent allocations: frame poses, palette LUT, one Spherical and one
  // trail per fiber.
  static constexpr size_t FOOTPRINT_BYTES =
      Filter::World::Orient::storage_bytes(W) +
      BakedPalette::required_arena_bytes() +
      ACTUAL_FIBERS * sizeof(math::Spherical) + alignof(math::Spherical) +
      ACTUAL_FIBERS * sizeof(Animation::VectorTrail<TRAIL_LEN>) +
      alignof(Animation::VectorTrail<TRAIL_LEN>) +
      ACTUAL_FIBERS * sizeof(math::Vector) + alignof(math::Vector);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HopfFibration persistent footprint exceeds the default "
                "partition; retune RINGS/PER_RING/TRAIL_LEN or carve arenas");

  // scratch_a stages one fiber's fragments and rasterize's sub-step cache.
  static_assert(TRAIL_LEN * sizeof(Fragment) +
                        Plot::rasterize_scratch_a_bytes<W>(0, TRAIL_LEN) <=
                    DEFAULT_SCRATCH_A_SIZE,
                "HopfFibration trail staging exceeds the default scratch_a "
                "budget; retune TRAIL_LEN or carve a larger scratch arena");

  // Keeps the projection denominator (1 + eps) - q3 positive at q3 == 1.
  static constexpr float STEREO_POLE_EPSILON = 0.001f;

  // Phases accumulate as wrapped fractions of their period ("turns") and scale
  // back to radians at use. tumble_angle_x also feeds the half-angle fold_base
  // term, so it wraps over 4pi to keep both terms continuous.
  static constexpr float FLOW_PERIOD = 2 * math::PI_F;
  static constexpr float TUMBLE_X_PERIOD = 4 * math::PI_F;
  static constexpr float TUMBLE_Y_PERIOD = 2 * math::PI_F;

  // Per-unit driver rates, in turns (radians-per-unit-speed / period).
  static constexpr float FLOW_RATE = (0.02f * 0.2f) / FLOW_PERIOD;
  static constexpr float TUMBLE_X_RATE = 0.003f / TUMBLE_X_PERIOD;
  static constexpr float TUMBLE_Y_RATE = 0.005f / TUMBLE_Y_PERIOD;

  // Wrapped phase accumulators in [0, 1) turns.
  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  math::Spherical *fibers = nullptr;
  Animation::VectorTrail<TRAIL_LEN> *trails = nullptr;

  // Per-frame values reused across all fibers (refreshed in advance_tumble).
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;
  float flow_rad = 0.0f, ty_rad = 0.0f;

  /**
   * @brief User-tunable parameters controlling fiber motion and shading.
   */
  struct Params {
    float flow_speed = 10.0f;  /**< Flow phase advance rate, in tuning units. */
    float tumble_speed = 2.0f; /**< 4D tumble rate, in tuning units. */
    float folding = 0.2f;      /**< Fold modulation depth, dimensionless. */
    float twist =
        4.0f; /**< Twist applied to azimuth per fold, dimensionless. */
    float alpha = 1.0f; /**< Global trail opacity multiplier in [0, 1]. */
  } params;

  // Borrowed by timeline-resident animations, with the phase accumulators and
  // params; declared before the Timeline to outlive it.
  math::Orientation<> orientation;
  /** @brief Timed trails: orientation at column times, then anti-aliasing. */
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      timed_pipeline{Filter::World::Orient(orientation),
                     Filter::Screen::AntiAlias<W, H>()};
  math::Vector *retired = nullptr; /**< Outgoing zero-opacity tail nodes. */
  bool trail_was_full = false; /**< Whether this frame evicted a tail node. */
  bool history_ready =
      false; /**< Whether initial orientation history is primed. */
  Timeline timeline;
  BakedPaletteStorage baked_sunset;

  // AA only; render_trails orients the points itself.
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> trail_pipeline;

  /**
   * @brief Seeds the base fibers as a spherical lattice.
   * @details A RINGS x PER_RING grid, evenly spaced in polar bands and in
   * azimuth around each band.
   */
  HS_COLD_MEMBER void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float polar = math::PI_F * (i + 0.5f) / RINGS;
      for (int j = 0; j < PER_RING; ++j) {
        float azimuth = 2 * math::PI_F * j / PER_RING;
        std::construct_at(&fibers[idx], azimuth, polar);
        ++idx;
      }
    }
  }

  /**
   * @brief Refreshes the per-frame cache from the wrapped phase accumulators.
   */
  HS_HOT_FLASH_MEMBER void advance_tumble() {
    // Scale wrapped [0,1)-turn phases back to radians (ax < 4pi, ty_rad < 2pi).
    const float ax = tumble_angle_x * TUMBLE_X_PERIOD;
    ty_rad = tumble_angle_y * TUMBLE_Y_PERIOD;
    flow_rad = flow_offset * FLOW_PERIOD;
    cx = math::fast_cosf(ax);
    sx = math::fast_sinf(ax);
    cy = math::fast_cosf(ty_rad);
    sy = math::fast_sinf(ty_rad);
    // fold_base needs sin(ax/2); ax*0.5 spans [0,2pi) continuously across the
    // 4pi wrap.
    fold_base = math::fast_sinf(ax * 0.5f) * 0.5f;
  }

  /**
   * @brief Projects a base fiber to an un-oriented direction.
   * @param i Fiber index in [0, ACTUAL_FIBERS).
   * @return Direction before render_trails orientation, after folding, twist, S3 lift, tumble, and
   * stereographic projection; falls back to (1, 0, 0) at a degenerate pole.
   * @details Pipeline: folding -> twist -> S3 -> tumble -> stereographic R3.
   */
  math::Vector hopf_project(size_t i) const {
    const math::Spherical &sph = fibers[i];
    float polar = sph.phi;
    float azimuth = sph.theta;

    // Folding: amplitude gated by the Folding slider so it persists when tumble
    // is frozen, though its phase still tracks ty_rad.
    float eta = polar / 2.0f;
    eta += math::fast_sinf(azimuth * 2.0f + ty_rad + fold_base) * 0.2f *
           params.folding;

    // Twist
    azimuth += eta * params.twist;

    // S3 point
    float beta = flow_rad + static_cast<float>(i) * PHASE_STEP;
    float cos_eta = math::fast_cosf(eta), sin_eta = math::fast_sinf(eta);
    float q0 = cos_eta * math::fast_cosf(azimuth + beta);
    float q1 = cos_eta * math::fast_sinf(azimuth + beta);
    float q2 = sin_eta * math::fast_cosf(beta);
    float q3 = sin_eta * math::fast_sinf(beta);

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
    return math::normalized_or(
        math::Vector(q0 * factor, q1 * factor, q2 * factor),
        math::Vector(1, 0, 0));
  }

  /**
   * @brief First trail index whose outgoing segment clears the alpha floor.
   * @param len Recorded trail length, at least 2.
   * @param alpha Master alpha, already known to clear MIN_VISIBLE_ALPHA.
   * @return Index in [0, len - 2]; the newest segment always survives, so the
   * staged polyline is never empty.
   * @details The tail fades to transparent, so leading points whose outgoing
   * segment peaks (at its newer endpoint) below MIN_ENCODABLE_ALPHA are
   * dropped. Per-channel rounding can still make a sub-floor peak encode nonzero.
   */
  static size_t trail_trim_start(size_t len, float alpha) {
    size_t first = 0;
    while (first + 1 < len &&
           static_cast<float>(first + 1) / (len - 1) * alpha <
               MIN_ENCODABLE_ALPHA)
      ++first;
    return first;
  }

  /**
   * @brief Renders all fiber trails as anti-aliased polylines.
   * @param canvas Target canvas to rasterize the trail polylines onto.
   * @details The tail fades from newest (opaque) to oldest (transparent).
   * Under an active segment clip, a trail with no visible edge is skipped and
   * the per-edge gate bits feed rasterize as its cull.
   */
  HS_HOT_FLASH_MEMBER void render_trails(Canvas &canvas) {
    // Alpha below one slider LSB: skip rasterizing; trails still record.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;
    HS_PROFILE(hf_render_trails);
    if (!Spacetime::flat_surface()) {
      render_timed_trails(canvas);
      return;
    }

    const ClipRegion &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const ClipRegion::XClip xc = cr.x_clip();

    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      const auto &trail = trails[i];
      size_t len = trail.length();
      if (len < 2)
        continue;

      // v0 keeps the full-trail parameterization, so kept segments render
      // unchanged.
      const size_t first = trail_trim_start(len, params.alpha);

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

      uint8_t bits[TRAIL_LEN];
      std::span<const uint8_t> edge_flags;
      if (clip_active) {
        HS_PROFILE(hf_trail_gate);
        if (!Plot::gate_trail_edges<W, H>(trail_pipeline, cr, xc, points, bits,
                                          0.0f))
          continue;
        edge_flags = {bits, points.size() - 1};
      }

      auto shader = [&](const math::Vector &, Fragment &f) HS_HOT_INLINE {
        const float t = fmaxf(0.0f, f.v0);
        Color4 c = baked_sunset.get(1.0f - t);
        c.alpha *= t * params.alpha; // fade out tail
        f.color = c;
      };
      {
        HS_PROFILE(hf_trail_raster);
        Plot::rasterize<W, H>(
            trail_pipeline, canvas, points, shader,
            {.projection = Plot::RasterProjection::geodesic(edge_flags)});
      }
    }
  }

  /** @brief Trail nodes one column walk draws. */
  static constexpr size_t COLUMN_CHUNK_VERTICES = 384;
  /** @brief Most lattice points of a steady trail edge the point route plots
   *  in place of the column walk. */
  static constexpr int STEADY_POINTS = 4;

  static_assert(COLUMN_CHUNK_VERTICES *
                            (sizeof(math::Vector) + sizeof(float) +
                             sizeof(uint8_t) + sizeof(Plot::Mesh::Edge)) +
                        (TRAIL_LEN + 1) *
                            (3 * sizeof(math::Vector) + sizeof(int)) +
                        8 * alignof(math::Vector) +
                        Plot::Mesh::column_scratch_a_bytes() <=
                    DEFAULT_SCRATCH_A_SIZE,
                "HopfFibration column chunk exceeds scratch_a");
  static_assert(Plot::Mesh::column_scratch_b_bytes(COLUMN_CHUNK_VERTICES) <=
                    DEFAULT_SCRATCH_B_SIZE,
                "HopfFibration node paths exceed scratch_b");

  /** @brief A chunk of trail nodes, as Plot::Mesh reads them. */
  struct NodeMesh {
    std::span<const math::Vector> vertices; /**< Unoriented trail nodes. */
  };

  /**
   * @brief Trail color at a progress and frame fraction.
   * @param progress Node progress, 0 at the full trail's tail and 1 at its
   * head.
   * @param time Frame fraction.
   * @param len Recorded trail length.
   * @return The faded palette color.
   * @details Over the frame every node slides one slot toward the tail, or
   * stretches while the trail fills, and the newest segment grows from the
   * previous head to the new one.
   */
  HS_HOT_INLINE Color4 trail_color(float progress, float time,
                                   size_t len) const {
    float t = progress;
    if (time != 1.0f) {
      const float previous =
          trail_was_full
              ? fminf(1.0f, t + 1.0f / (len - 1))
              : (len > 2 ? fminf(1.0f, t * (len - 1) / (len - 2)) : t);
      t = fmaxf(0.0f, previous + (t - previous) * time);
    }
    t = fmaxf(0.0f, t);
    Color4 c = baked_sunset.get(1.0f - t);
    c.alpha *= t * params.alpha;
    const float birth = fmaxf(0.0f, progress * static_cast<float>(len - 1) -
                                        static_cast<float>(len - 2));
    if (birth > time)
      c.alpha = 0.0f;
    return c;
  }

  /**
   * @brief Renders all fiber trails at their columns' times.
   * @param canvas Target canvas.
   * @details Trail nodes hold still in the source; the world stage turns them
   * by the orientation at each column's time. Each trail is a chain from its
   * head to its retired tail, whose fully faded last node is not plotted. The
   * walk's swept bounds cull each trail whole, then edge by edge; edges that
   * keep a short lattice over the frame are their lattice points, single
   * samples the point route plots, and the column walk draws the rest in
   * chunks.
   */
  HS_HOT_FLASH_MEMBER void render_timed_trails(Canvas &canvas) {
    const auto cull = Plot::SweptEdgeCull<W>::template of<H>(canvas);
    const float spin = timed_pipeline.sweep_arc();
    const Spacetime::PoseTable &poses = *timed_pipeline.pose_table(canvas);
    const math::RotationMatrix start_pose = poses.at(0);
    const math::RotationMatrix end_pose = poses.at(poses.times());
    math::RotationMatrix end_to_start;
    end_to_start.r0 = {math::dot(start_pose.r0, end_pose.r0),
                       math::dot(start_pose.r0, end_pose.r1),
                       math::dot(start_pose.r0, end_pose.r2)};
    end_to_start.r1 = {math::dot(start_pose.r1, end_pose.r0),
                       math::dot(start_pose.r1, end_pose.r1),
                       math::dot(start_pose.r1, end_pose.r2)};
    end_to_start.r2 = {math::dot(start_pose.r2, end_pose.r0),
                       math::dot(start_pose.r2, end_pose.r1),
                       math::dot(start_pose.r2, end_pose.r2)};
    Plot::PointRoute<W, H> route;
    route.prepare(canvas, spin + 2.0f * math::PI_F / W);

    ScratchScope chunk_scope(scratch_arena_a);
    auto *vertices =
        scratch_arena_a.allocate_n<math::Vector>(COLUMN_CHUNK_VERTICES);
    auto *progress = scratch_arena_a.allocate_n<float>(COLUMN_CHUNK_VERTICES);
    auto *lengths = scratch_arena_a.allocate_n<uint8_t>(COLUMN_CHUNK_VERTICES);
    ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_a, COLUMN_CHUNK_VERTICES);
    constexpr int TRAIL_POINTS = TRAIL_LEN + 1;
    auto *sources = scratch_arena_a.allocate_n<math::Vector>(TRAIL_POINTS);
    auto *starts = scratch_arena_a.allocate_n<math::Vector>(TRAIL_POINTS);
    auto *world = scratch_arena_a.allocate_n<math::Vector>(TRAIL_POINTS);
    auto *slots = scratch_arena_a.allocate_n<int>(TRAIL_POINTS);
    size_t count = 0;

    auto shader = [&](const math::Vector &, Fragment &f,
                      float time) HS_HOT_FLASH_MEMBER {
      const Plot::Mesh::Edge &edge = edges[static_cast<size_t>(f.v2)];
      f.color = trail_color(progress[edge.u] +
                                (progress[edge.v] - progress[edge.u]) * f.v0,
                            time, lengths[edge.u]);
    };
    Plot::GeometryCallback geometry;
    geometry.source_positions = false;
    geometry.chained = true;
    auto flush = [&]() HS_HOT_INLINE {
      if (count == 0)
        return;
      Plot::Mesh::draw_column_edges<W, H>(timed_pipeline, canvas,
                                          NodeMesh{{vertices, count}}, edges,
                                          shader, geometry, 0.0f);
      count = 0;
      edges.clear();
    };

    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      const auto &trail = trails[i];
      const size_t len = trail.length();
      if (len < 2)
        continue;
      const size_t first = trail_trim_start(len, params.alpha);
      const bool keep_tail = trail_was_full && first == 0;
      const int points = static_cast<int>(len - first) + (keep_tail ? 1 : 0);
      for (int j = 0; j < points; ++j) {
        sources[j] = keep_tail && j == 0
                         ? retired[i]
                         : trail.get(first + static_cast<size_t>(j) -
                                     (keep_tail ? 1 : 0));
        world[j] = Plot::EdgeWorldsheet::normalize(end_pose.apply(sources[j]));
      }
      float reach_cosine = 1.0f;
      const math::Vector &centre = world[points / 2];
      for (int j = 0; j < points; ++j)
        reach_cosine = fminf(reach_cosine, math::dot(centre, world[j]));
      if (cull.rejects(&centre, 1,
                       spin + sqrtf(fmaxf(0.0f, 2.0f - 2.0f * reach_cosine))))
        continue;
      auto node_progress = [&](int j) HS_HOT_INLINE {
        return static_cast<float>(static_cast<int>(first) + j -
                                  (keep_tail ? 1 : 0)) /
               static_cast<float>(len - 1);
      };
      for (int j = 0; j < points; ++j) {
        starts[j] =
            Plot::EdgeWorldsheet::normalize(start_pose.apply(sources[j]));
        slots[j] = -1;
      }
      if (count + static_cast<size_t>(points) > COLUMN_CHUNK_VERTICES)
        flush();
      for (int j = points - 1; j > 0; --j) {
        const math::Vector pair[2] = {world[j], world[j - 1]};
        if (cull.rejects(pair, 2, spin))
          continue;
        const float head = node_progress(j), tail = node_progress(j - 1);
        if (Plot::Mesh::plot_steady_edge<W, H>(
                route, trail_pipeline, canvas, {starts[j], starts[j - 1]},
                {world[j], world[j - 1]}, STEADY_POINTS,
                [&](float s, float time) HS_HOT_INLINE {
                  return trail_color(head + (tail - head) * s, time, len);
                },
                &end_to_start))
          continue;
        for (const int k : {j, j - 1})
          if (slots[k] < 0) {
            slots[k] = static_cast<int>(count++);
            vertices[slots[k]] = sources[k];
            progress[slots[k]] = node_progress(k);
            lengths[slots[k]] = static_cast<uint8_t>(len);
          }
        edges.push_back({static_cast<uint16_t>(slots[j]),
                         static_cast<uint16_t>(slots[j - 1])});
      }
    }
    flush();
  }
};
