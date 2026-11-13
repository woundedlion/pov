/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
/**
 * @file mesh_motion.h
 * @brief Capture-time motion planning for a mesh draw whose vertices move
 *        along a model path under optional ripples.
 */
#include <cstddef>
#include <cstdint>
#include "core/animation/animation.h"
#include "core/animation/ripple_motion.h"
#include "core/animation/transformer.h"
#include "core/render/scan/mesh.h"

namespace Scan {

/** @brief A ripple pool with no ripples. */
struct NoRipples {
  /** @return 0. */
  int active_count() const { return 0; }
  /**
   * @brief Never called; a pool with no active ripples has no parameters.
   * @return Default parameters.
   */
  const Animation::RippleParams &active_params(int) const { return NONE; }

private:
  inline static const Animation::RippleParams NONE{};
};

/** @brief A draw's camera at the current endpoint. */
struct EndpointCamera {
  math::Quaternion q; ///< The orientation's current rotation.
  bool turns = true;  ///< False for a path already in camera space.

  /**
   * @brief Captures an orientation's current rotation.
   * @tparam CAPACITY History capacity of the orientation.
   * @param orientation The draw's orientation.
   */
  template <int CAPACITY>
  HS_HOT_INLINE explicit EndpointCamera(
      const math::Orientation<CAPACITY> &orientation)
      : q(orientation.get()) {}

  /** @return A camera for a path already in camera space. */
  HS_HOT_INLINE static EndpointCamera none() { return EndpointCamera(); }

  /**
   * @param v A path vector.
   * @return @p v in camera space, as math::Orientation::orient() turns it.
   */
  HS_HOT_INLINE math::Vector operator()(const math::Vector &v) const {
    return turns ? math::rotate(v, q) : v;
  }

private:
  HS_HOT_INLINE EndpointCamera() : turns(false) {}
};

/**
 * @brief Where a draw's model-space vertices stand over the frame.
 * @details Endpoint vertices are unit length.
 */
struct ModelPath {
  const math::Vector *current = nullptr;  ///< At the current endpoint.
  const math::Vector *previous = nullptr; ///< At the previous endpoint.
  /** Shading whose leg moves the model vertex by vertex, or null. */
  const Animation::OpLeg::Shading *leg = nullptr;
  /// At the current endpoint, packed, in place of current.
  const math::Snorm3 *packed = nullptr;

  /**
   * @brief The path of a leg frame's drawn mesh.
   * @param current The mesh's vertices at the current endpoint.
   * @param sh The leg's shading for the frame.
   * @return Along the leg when it moves vertex by vertex, else from its
   *         previous vertices.
   */
  static ModelPath of(const math::Vector *current,
                      const Animation::OpLeg::Shading &sh) {
    return {current, sh.previous_vertices, sh.leg ? &sh : nullptr};
  }

  /** @return Whether the model moves over the frame. */
  bool moves() const { return leg || previous; }

  /** @brief One fraction of the frame, with the leg's parameter there. */
  struct Moment {
    float u = 0.0f; ///< Frame fraction; 0 and 1 are the endpoints.
    float k = 0.0f; ///< The leg's sweep parameter at `u`; 0 without a leg.
  };

  /**
   * @brief A fraction of the frame, resolved for at().
   * @param u Frame fraction; past either endpoint the model continues
   *        along its path, holding at a leg's own ends.
   * @return The moment.
   */
  Moment moment(float u) const {
    return {u, leg ? leg->parameter_at(u) : 0.0f};
  }

  /**
   * @brief A model vertex at a moment of the frame.
   * @param i Vertex index.
   * @param m The moment, as moment() gives it.
   * @return The vertex on its leg's path, normalized-lerped from
   * `previous`, or `current` when the model holds.
   */
  HS_HOT_INLINE math::Vector at(size_t i, const Moment &m) const {
    if (leg)
      return leg->vertex_at_parameter(i, m.k);
    const math::Vector now = packed ? packed[i].decode() : current[i];
    if (previous)
      return (previous[i] + (now - previous[i]) * m.u).normalized();
    return now;
  }

  /**
   * @brief A model vertex at its capture fraction.
   * @param i Vertex index.
   * @param u Frame fraction.
   * @return Vertex on the model path.
   */
  HS_HOT_INLINE math::Vector at(size_t i, float u) const {
    return at(i, moment(u));
  }
};

/**
 * @brief Bound on how far a pool's active ripples move any point over the
 *        frame.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @param ripples The pool.
 * @return Radians.
 */
template <typename Pool> HS_FLASH_MEMBER float ripple_arc(const Pool &ripples) {
  float arc = 0.0f;
  for (int k = 0; k < ripples.active_count(); ++k)
    arc += ripples.active_params(k).step_arc();
  return arc;
}

/**
 * @brief Maximum displacement from the undeformed model over capture times.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @param ripples The pool deforming the model.
 * @param motion motion_bound() of the draw.
 * @return Angular displacement in radians.
 */
template <typename Pool>
HS_HOT_FLASH_MEMBER float visibility_reach(const Pool &ripples, float motion) {
  float reach = 2.0f * fmaxf(motion - ripple_arc(ripples), 0.0f);
  for (int k = 0; k < ripples.active_count(); ++k) {
    const auto &p = ripples.active_params(k);
    reach += fmaxf(fabsf(p.amplitude), fabsf(p.previous_amplitude)) +
             fabsf(p.amplitude - p.previous_amplitude);
  }
  return reach;
}

/**
 * @brief Bound on how far a mesh draw's vertices move over the frame.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @param path The draw's model path.
 * @param count Number of vertices.
 * @param ripples The pool deforming the model.
 * @return Radians; the ripples' bound plus the model's largest motion.
 */
template <typename Pool>
HS_FLASH_MEMBER float motion_bound(const ModelPath &path, size_t count,
                                   const Pool &ripples) {
  float min_cos = 1.0f;
  if (path.moves()) {
    const auto start = path.moment(0.0f);
    const auto end = path.moment(1.0f);
    for (size_t i = 0; i < count; ++i)
      min_cos = fminf(min_cos, math::dot(path.at(i, start), path.at(i, end)));
  }
  return ripple_arc(ripples) + acosf(hs::clamp(min_cos, -1.0f, 1.0f));
}

/// Arc, radians, of one unit of a vertex's recorded model motion.
inline constexpr float MOTION_UNIT = math::PI_F / 65535.0f;

/**
 * @brief Flags the vertices the model or an active ripple's last step can
 *        move.
 * @tparam W Canvas width in pixels.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @param path The draw's model path.
 * @param count Number of vertices.
 * @param leg_arc The model's largest motion (radians).
 * @param moving Receives 1 for each vertex the model moves more than
 *        HELD_VERTEX_ARC, or that stands within @p leg_arc of a moving
 *        ripple's swept band; 0 for the others.
 * @param motion Receives each vertex's model motion over the frame, in
 *        MOTION_UNIT arcs rounded up.
 * @param arena Arena holding the bands for the call.
 * @param visible Per-vertex board visibility bits, or null.
 * @param ripples The pool deforming the model.
 * @details A vertex outside every band is left in place by each ripple in
 * turn.
 */
template <int W, typename Pool>
HS_FLASH_MEMBER void flag_moving(const ModelPath &path, size_t count,
                                 float leg_arc, uint8_t *moving,
                                 uint16_t *motion, Arena &arena,
                                 const uint8_t *visible, const Pool &ripples) {
  ScratchScope scope(arena);
  int band_count = 0;
  const RippleMotion::RippleBand *bands =
      RippleMotion::ripple_bands(ripples, leg_arc, arena, band_count);
  const float held_cos = cosf(HELD_VERTEX_ARC<W>);
  const auto start = path.moment(0.0f);
  const auto end = path.moment(1.0f);
  for (size_t i = 0; i < count; ++i) {
    if (!mesh_visible(visible, i)) {
      moving[i] = 0;
      motion[i] = 0;
      continue;
    }
    const math::Vector model = path.at(i, end);
    bool moves = false;
    motion[i] = 0;
    if (path.moves()) {
      const float c = math::dot(path.at(i, start), model);
      motion[i] = static_cast<uint16_t>(
          fminf(ceilf(cosine_arc(c) * (1.0f / MOTION_UNIT)), 65535.0f));
      moves = c < held_cos;
    }
    for (int b = 0; b < band_count && !moves; ++b)
      moves = bands[b].holds(model);
    moving[i] = moves;
  }
}

#if HS_ENABLE_TEST_HOOKS
/** @brief Host work counts for planned mesh motion. */
struct MeshMotionCounters {
  size_t passes = 0;         ///< capture_position() evaluations.
  size_t candidates = 0;     ///< Ripple bands visited.
  size_t interpolations = 0; ///< Ripple samples interpolated.
};
/** @brief Instrumentation for the current mesh motion measurement. */
inline MeshMotionCounters mesh_motion_counters;
/// capture_position() takes the uncached ripple loop.
inline bool mesh_motion_reference = false;
#endif

/** @brief Ripple capture counters for MeshMotionCounters. */
struct MotionRippleCounters {
  /** @brief Counts a visited ripple band. */
  static void candidate() {
#if HS_ENABLE_TEST_HOOKS
    ++mesh_motion_counters.candidates;
#endif
  }
  /** @brief Counts an interpolated ripple sample. */
  static void interpolated() {
#if HS_ENABLE_TEST_HOOKS
    ++mesh_motion_counters.interpolations;
#endif
  }
};

/** @brief One active ripple's capture state in a MotionPlan. */
using MotionRipple = RippleMotion::MotionRipple<MotionRippleCounters>;

/** @brief Vertex deformation state of one mesh draw. */
struct MotionPlan {
  bool active = false;            ///< Whether vertices deform.
  uint8_t *moving = nullptr;      ///< Nonzero per moving vertex.
  uint16_t *start = nullptr;      ///< Moving vertices' columns at the
                                  ///< frame's start (vertex_column()).
  uint16_t *columns = nullptr;    ///< Per-face seeing columns.
  MotionRipple *during = nullptr; ///< One per active ripple.
};

/**
 * @brief Evaluates one vertex on its model and ripple paths in camera space.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @tparam Camera Callable mapping a model-space vector to camera space.
 * @param plan Frame's ripple cache and motion selection.
 * @param path Model-space vertex path.
 * @param ripples The pool the plan was made with.
 * @param camera The draw's camera at the current endpoint.
 * @param i Vertex index.
 * @param u Frame fraction.
 * @return Vertex at the requested frame fraction.
 */
template <typename Pool, typename Camera>
HS_HOT_INLINE math::Vector
capture_position(const MotionPlan &plan, const ModelPath &path,
                 const Pool &ripples, const Camera &camera, size_t i, float u) {
#if HS_ENABLE_TEST_HOOKS
  ++mesh_motion_counters.passes;
  if (mesh_motion_reference) {
    math::Vector vertex = path.at(i, u);
    for (int k = 0; k < ripples.active_count(); ++k) {
      const auto &p = ripples.active_params(k);
      auto &during = plan.during[k].params;
      during.phase = p.previous_phase + (p.phase - p.previous_phase) * u;
      during.amplitude =
          p.previous_amplitude + (p.amplitude - p.previous_amplitude) * u;
      ++mesh_motion_counters.candidates;
      ++mesh_motion_counters.interpolations;
      vertex = ripple_params_transform(vertex, during);
    }
    return camera(vertex);
  }
#endif
  math::Vector vertex = path.at(i, u);
  for (int k = 0; k < ripples.active_count(); ++k)
    vertex = plan.during[k].transform(vertex, u);
  return camera(vertex);
}

/**
 * @brief Scans a mesh with its planned model and ripple motion.
 * @tparam W Display width.
 * @tparam H Display height.
 * @tparam Path Model-space vertex path type.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @tparam ScanFn Callable () scanning the mesh into the canvas's open window.
 * @param canvas Canvas with column-time state when the plan is active.
 * @param mesh Endpoint mesh, overwritten by the captured vertices.
 * @param plan Frame's ripple cache and motion selection.
 * @param path Model-space vertex path; may be null when the plan is inactive.
 * @param ripples The pool the plan was made with.
 * @param camera The draw's camera at the current endpoint.
 * @param scratch Arena for face selection and vertex placement tables.
 * @param reach Bound on vertex motion in radians.
 * @param scan Scans the mesh with its current vertices and selected faces.
 */
template <int W, int H, typename Path, typename Pool, typename ScanFn>
HS_HOT_INLINE void deform_planned(Canvas &canvas, MeshState &mesh,
                                  const MotionPlan &plan, const Path *path,
                                  const Pool &ripples,
                                  const EndpointCamera &camera, Arena &scratch,
                                  float reach, ScanFn &&scan) {
  deform_mesh<W, H>(
      canvas, mesh, plan.active, plan.columns, plan.moving, plan.start, scratch,
      reach,
      [&](size_t i, float u) HS_HOT_INLINE {
        return capture_position(plan, *path, ripples, camera, i, u);
      },
      static_cast<ScanFn &&>(scan));
}

/**
 * @brief Plans the vertex deformation of one mesh draw.
 * @tparam W Canvas width in pixels.
 * @tparam Pool Ripple pool exposing active_count() and active_params().
 * @tparam Camera Callable mapping a model-space vector to camera space.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The mesh, placed at the current endpoint.
 * @param path The draw's model path.
 * @param bound motion_bound() of the draw.
 * @param ripples The pool deforming the model.
 * @param camera The draw's camera at the current endpoint.
 * @param arena Arena backing the plan's tables for the caller's scope.
 * @param start_arena Arena backing the plan's frame-start columns for the
 *        caller's scope.
 * @return Inactive with temporal sampling off or when every visible vertex
 *         holds.
 * @details The faces the clip can see with a moving vertex are picked; the
 * largest vertex motion activates deformation. A picked face whose vertices
 * all move at most HELD_VERTEX_ARC draws once. A ripple crossing disjoint
 * endpoint supports retains the full frame motion bound. A placed vertex strays from
 * the current endpoint by at most its motion plus the bend of its path
 * between the endpoints, which a column of arc covers. Each candidate vertex
 * is evaluated once, at the frame's start; the scan captures it from its
 * column there.
 */
template <int W, typename Pool, typename Camera>
HS_FLASH_MEMBER MotionPlan plan_motion(const Canvas &canvas,
                                       const MeshState &mesh,
                                       const ModelPath &path, float bound,
                                       const Pool &ripples,
                                       const Camera &camera, Arena &arena,
                                       Arena &start_arena) {
  MotionPlan plan;
  if (!mesh_motion_active<W>(bound))
    return plan;
  const float leg_arc = bound - ripple_arc(ripples);
  const size_t count = mesh.vertices.size();
  plan.moving = arena.allocate_n<uint8_t>(count);
  // Each vertex's model motion stands in its start column's place until
  // the vertex's column replaces it.
  plan.start = start_arena.allocate_n<uint16_t>(count);
  flag_moving<W>(path, count, leg_arc, plan.moving, plan.start, arena,
                 canvas.column_basis()
                     ? canvas.column_basis()->mesh.visible_vertices
                     : nullptr,
                 ripples);
  plan.columns = arena.allocate_n<uint16_t>(2 * mesh.num_faces());
  plan.during = arena.allocate_n<MotionRipple>(ripples.active_count() + 1);
  ScratchScope radii_scope(arena);
  float *radii = arena.allocate_n<float>(mesh.num_faces());
  if (pick_moving_faces<W>(canvas, mesh, plan.moving, bound, plan.columns,
                           radii) == 0)
    return plan;
  HS_CHECK(canvas.column_basis() && canvas.column_basis()->azimuth,
           "mesh motion planning needs the frame's column-time azimuth state");
  bool interior_crossing = false;
  for (int k = 0; k < ripples.active_count(); ++k) {
    const auto &p = ripples.active_params(k);
    interior_crossing |= fabsf(p.phase - p.previous_phase) >
                         4.0f * p.half_width() + 2.0f * leg_arc;
    plan.during[k].prepare(ripples.active_params(k));
    plan.during[k].widen(ripples.active_params(k));
  }
  // The stepping ripples' bands, turned with the camera to the scanned
  // mesh.
  int band_count = 0;
  RippleMotion::RippleBand *bands =
      RippleMotion::ripple_bands(ripples, leg_arc, arena, band_count);
  for (int b = 0; b < band_count; ++b)
    bands[b].center = camera(bands[b].center);
  const float held_cos = cosf(HELD_VERTEX_ARC<W>);
  float min_cos = 1.0f;
  for (size_t i = 0; i < count; ++i) {
    if (!plan.moving[i])
      continue;
    const float model_arc = plan.start[i] * MOTION_UNIT;
    const math::Vector start =
        capture_position(plan, path, ripples, camera, i, 0.0f);
    plan.start[i] = vertex_column(start, *canvas.column_basis(), 0.0f, W);
    const math::Vector &end = mesh.vertices[i];
    const float c = math::dot(start, end);
    min_cos = fminf(min_cos, c);
    if (!(c < held_cos) && !interior_crossing) {
      plan.moving[i] = 0;
      continue;
    }
    // Its model motion over the frame, with the steps of the ripples whose
    // bands it stands in.
    float reach = model_arc;
    for (int b = 0; b < band_count; ++b)
      if (bands[b].holds(end))
        reach += bands[b].step;
    plan.moving[i] = vertex_flags(reach);
  }
  const float arc =
      interior_crossing ? bound : acosf(hs::clamp(min_cos, -1.0f, 1.0f));
  if (drop_held_faces(mesh, plan.moving, plan.columns) == 0)
    return plan;
  plan.active = mesh_motion_active<W>(arc);
  if (plan.active)
    moving_face_columns<W>(canvas, mesh, arc + math::TWO_PI_F / W, plan.columns,
                           radii);
  return plan;
}

} // namespace Scan
