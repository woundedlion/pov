/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/render/scan/mesh.h"

/** @file mesh_capture.cpp
 *  @brief Per-face column-time capture of a deforming mesh scan.
 */
namespace Scan {

HS_HOT_FLASH_MEMBER FaceCaptureInterval
face_capture_interval(const float *angles, size_t count, float reach,
                      float nearest, float sweep) {
  if (count == 0)
    return {-1.0f, 0.0f};
  const bool held =
      reach > 0.0f && reach + sweep >= nearest * (math::PI_F / 4.0f);
  float winding = 0.0f;
  float lo = angles[0], hi = lo;
  for (size_t k = 0; k < count; ++k) {
    const float delta = angles[(k + 1) % count] - angles[k];
    winding += delta - 2.0f * roundf(delta * 0.5f);
    const float angle =
        angles[k] + 2.0f * roundf((angles[0] - angles[k]) * 0.5f);
    lo = fminf(lo, angle);
    hi = fmaxf(hi, angle);
  }
  if (fabsf(winding) > 1.0f)
    return {-2.0f, 2.0f};
  if (hi - lo >= 1.0f || held)
    return {-1.0f, hi - lo};
  const float shift = floorf(lo);
  return {lo - shift, hi - shift};
}

} // namespace Scan

namespace Scan::detail {

namespace {

/** @brief How a deformation scan draws a face. */
enum FaceKind : uint8_t {
  HELD_FACE,   ///< Holds over the frame; drawn in the shared pass.
  SHARED_FACE, ///< Clear of both arm boundaries; drawn in the shared pass.
  SEAM_FACE,   ///< Near one arm boundary; drawn in each arm's shared pass.
  CROSS_FACE,  ///< Reaches past a seam face's reach; one copy per arm.
  POLAR_FACE,  ///< Reaches a pole; drawn in each half-arm window's pass.
  CULLED_FACE, ///< Outside the board's conservative visibility zone.
  SWEPT_FACE,  ///< Reaches a pole and moves too far for one placement;
               ///< swept one column pair at a time.
};

/// Square pixels of step a polar face's single placement may show before
/// it sweeps: the pixels the placement strays times the pixels the face
/// reaches from the polar axis.
constexpr float SWEEP_STEP_AREA = 20.0f;

/**
 * @brief How a face reaching a pole draws.
 * @param d Deformation scan, its pool at the current endpoint.
 * @param face Face index.
 * @param span Its vertices' azimuth span in half turns; 2 for a face that
 *        winds about the pole.
 * @param motion_cache Per-vertex frame motion, negative until evaluated.
 * @return SWEPT_FACE when one placement per half arm would show a step of
 *         more than SWEEP_STEP_AREA, else POLAR_FACE.
 * @details A placement stands within half a frame of a column that sees the
 * face, and a half turn of azimuth is a frame of column time, so it strays
 * by the face's fastest motion over that part of a frame. The step shows
 * along the face from the pole out to its farthest vertex. A vertex's
 * motion over the frame is bounded by its chord from the frame's start to
 * its end plus twice its path's bulge at the frame's middle.
 */
HS_HOT_FLASH_MEMBER FaceKind polar_kind(const SweptMesh &d, size_t face,
                                        float span, float *motion_cache) {
  const MeshState &mesh = d.mesh;
  const size_t count = mesh.get_face_counts_data()[face];
  const uint16_t *indices =
      mesh.get_faces_data() + mesh.get_face_offsets_data()[face];
  float motion = 0.0f, axial = 0.0f;
  for (size_t k = 0; k < count; ++k) {
    const size_t i = indices[k];
    const math::Vector end = mesh.vertices[i];
    axial = fmaxf(axial, end.x * end.x + end.z * end.z);
    if (!d.moving[i])
      continue;
    float vertex_motion = motion_cache[i];
    if (vertex_motion < 0.0f) {
      const math::Vector start =
          d.callbacks.evaluate(d.callbacks.vertex_context, i, 0.0f);
      const math::Vector middle =
          d.callbacks.evaluate(d.callbacks.vertex_context, i, 0.5f);
#if HS_ENABLE_TEST_HOOKS
      mesh_capture_counters.evaluations += 2;
#endif
      const math::Vector bulge = middle - (start + end) * 0.5f;
      vertex_motion = sqrtf(math::distance_squared(start, end)) +
                      2.0f * sqrtf(math::dot(bulge, bulge));
      motion_cache[i] = vertex_motion;
    }
    motion = fmaxf(motion, vertex_motion);
  }
  // Radians to pixels: a column's arc, which is also a row's.
  const float pixel = d.canvas.width() / math::TWO_PI_F;
  const float stray = motion * fminf(span, 0.5f) * pixel;
  return stray * sqrtf(axial) * pixel > SWEEP_STEP_AREA ? SWEPT_FACE
                                                        : POLAR_FACE;
}

/// Vertex state: its endpoint stands in the first half of its arm's sweep.
constexpr uint8_t VERTEX_EARLY = 1;
/// Vertex state: its endpoint stands on the second arm, past half a turn.
constexpr uint8_t VERTEX_ARM1 = 2;
/// Vertex state: the pool holds its capture on its own arm.
constexpr uint8_t POOL_OWN = 4;
/// Vertex state: the pool holds its capture across the arm boundary nearer
/// its endpoint.
constexpr uint8_t POOL_ACROSS = 8;

/// Vertex state: classification has stored the vertex's endpoint tables.
constexpr uint8_t VERTEX_SEEN = 16;
/// Vertex state: its reach can carry it a quarter of a right angle about the
/// polar axis, so it holds at its endpoint.
constexpr uint8_t VERTEX_HELD = 32;
/// Vertex state: the pool holds its capture nearest a quarter of the arm.
constexpr uint8_t POOL_FIRST_HALF = 64;
/// Vertex state: the pool holds its capture nearest three quarters of the
/// arm.
constexpr uint8_t POOL_SECOND_HALF = 128;

/// Half turns per unit of a stored vertex azimuth.
constexpr float TURNS_SCALE = 1.0f / 32768.0f;

/** @brief One deformation scan's per-vertex and per-face tables. */
struct FaceTables {
  uint8_t *kinds;  ///< FaceKind per face.
  uint8_t *state;  ///< Vertex state bits per vertex.
  uint16_t *turns; ///< Endpoint azimuth per classified vertex, in half turns
                   ///< from -1, scaled by 1 / TURNS_SCALE.
  uint16_t *ends;  ///< Endpoint column per classified vertex
                   ///< (vertex_column() at u = 1).
};

/**
 * @brief A classified vertex's endpoint azimuth.
 * @param tables Scan tables.
 * @param i Vertex index.
 * @return Half turns in [-1, 1).
 */
HS_HOT_INLINE float vertex_turns(const FaceTables &tables, size_t i) {
  return tables.turns[i] * TURNS_SCALE - 1.0f;
}

/**
 * @brief Whether the pool holds the requested unwrapped capture branch.
 * @param d Deformation scan.
 * @param tables Per-vertex capture labels.
 * @param i Vertex index.
 * @param centre Centre of the requested branch.
 * @return True when an existing label has the same unwrapped start.
 */
HS_O3_FN HS_HOT_INLINE bool capture_pool_matches(const SweptMesh &d,
                                                 const FaceTables &tables,
                                                 size_t i, float centre) {
  const uint8_t STATE = tables.state[i];
  if (!(STATE & (POOL_OWN | POOL_ACROSS | POOL_FIRST_HALF | POOL_SECOND_HALF)))
    return false;
  const float START = d.start[i] * (1.0f / 65536.0f);
  const float FIRST = unwrap_capture(START, centre);
  const bool EARLY = STATE & VERTEX_EARLY;
  return ((STATE & POOL_OWN) &&
          FIRST == unwrap_capture(START, EARLY ? 0.25f : 0.75f)) ||
         ((STATE & POOL_ACROSS) &&
          FIRST == unwrap_capture(START, EARLY ? 1.25f : -0.25f)) ||
         ((STATE & POOL_FIRST_HALF) && FIRST == unwrap_capture(START, 0.25f)) ||
         ((STATE & POOL_SECOND_HALF) && FIRST == unwrap_capture(START, 0.75f));
}

/**
 * @brief Captures a vertex on a temporal branch.
 * @param d Deformation scan.
 * @param tables Scan tables, holding the vertex's endpoint column.
 * @param centre Arm fraction the branch lies within half a frame of.
 * @param i Vertex index.
 * @return Captured vertex.
 */
HS_HOT_FLASH_MEMBER math::Vector capture_branch(const SweptMesh &d,
                                                const FaceTables &tables,
                                                float centre, size_t i) {
  return capture_vertex(d.start[i], tables.ends[i], centre, d.canvas.width(),
                        [&](float u) HS_HOT_INLINE {
                          return d.callbacks.evaluate(
                              d.callbacks.vertex_context, i, u);
                        });
}

/**
 * @brief Places a moving vertex's capture on its own arm or across the arm
 *        boundary nearer its endpoint, unless the pool already holds it.
 * @param d Deformation scan.
 * @param tables Scan tables, holding the vertex's state.
 * @param i Vertex index.
 * @param across Whether the capture continues across the boundary.
 */
HS_HOT_FLASH_MEMBER void place_vertex(const SweptMesh &d,
                                      const FaceTables &tables, size_t i,
                                      bool across) {
  uint8_t &state = tables.state[i];
  const uint8_t pool = across ? POOL_ACROSS : POOL_OWN;
  if (!d.moving[i] || (state & (pool | VERTEX_HELD)))
    return;
  const bool early = state & VERTEX_EARLY;
  const float centre =
      across ? (early ? 1.25f : -0.25f) : (early ? 0.25f : 0.75f);
  if (capture_pool_matches(d, tables, i, centre)) {
    state |= pool;
    return;
  }
  d.mesh.vertices[i] = capture_branch(d, tables, centre, i);
  state = (state & (VERTEX_EARLY | VERTEX_ARM1)) | pool;
}

/**
 * @brief Places a moving vertex's capture nearest the middle of one half of
 *        the arm, unless the pool already holds it or the vertex holds.
 * @param d Deformation scan.
 * @param tables Scan tables, holding the vertex's state.
 * @param i Vertex index.
 * @param half 0 for the arm's first half, 1 for its second.
 * @details The capture is the vertex's own-arm one when its endpoint stands
 * in that half.
 */
HS_HOT_FLASH_MEMBER void
place_half(const SweptMesh &d, const FaceTables &tables, size_t i, int half) {
  uint8_t &state = tables.state[i];
  const uint8_t pool = half ? POOL_SECOND_HALF : POOL_FIRST_HALF;
  if (!d.moving[i] || (state & (pool | VERTEX_HELD)))
    return;
  const bool own = bool(state & VERTEX_EARLY) == (half == 0);
  if (capture_pool_matches(d, tables, i, half ? 0.75f : 0.25f)) {
    state |= pool | (own ? POOL_OWN : uint8_t(0));
    return;
  }
  d.mesh.vertices[i] = capture_branch(d, tables, half ? 0.75f : 0.25f, i);
  state = (state & (VERTEX_EARLY | VERTEX_ARM1)) | pool |
          (own ? POOL_OWN : uint8_t(0));
}

/**
 * @brief A moving vertex's capture at an arm time, from the pool's arm-side
 *        captures where one holds that time.
 * @param d Deformation scan.
 * @param tables Scan tables.
 * @param i Vertex index.
 * @param time Arm fraction of the capture's branch, past 1 or before 0 for
 *        a time across an arm boundary.
 * @return The captured vertex, or the pool's vertex for one that does not
 *         move.
 */
HS_HOT_FLASH_MEMBER math::Vector
vertex_at(const SweptMesh &d, const FaceTables &tables, size_t i, float time) {
  if (!d.moving[i] || (tables.state[i] & VERTEX_HELD))
    return d.mesh.vertices[i];
  const float turns = vertex_turns(tables, i);
  const float shift = roundf(time - (turns - floorf(turns)));
  const bool early = tables.state[i] & VERTEX_EARLY;
  if (shift == 0.0f || shift == (early ? 1.0f : -1.0f)) {
    place_vertex(d, tables, i, shift != 0.0f);
    return d.mesh.vertices[i];
  }
  return capture_branch(d, tables, time, i);
}

/**
 * @brief Classifies a moving face by its vertices' arm fractions.
 * @param d Deformation scan, its pool at the current endpoint.
 * @param face Face index.
 * @param tables Scan tables; receive the face's vertex azimuths and states.
 * @param sweep The column basis's sweep over the frame.
 * @param motion_cache Per-vertex frame motion for polar classification.
 * @return The face's kind.
 * @details A face reaches a pole when it winds around it, when its fastest
 * vertex can reach the pole from one of its vertices (face_capture_interval()),
 * or when it spans half a turn of azimuth with the arc its scan
 * can reach past its vertices, d.pad at its lowest vertex's distance from
 * the polar axis: every column pair's meridian can cross it. A face is near
 * an arm boundary when that arc meets one. A seam face's vertices lie within
 * a quarter frame of the boundary.
 */
HS_HOT_FLASH_MEMBER FaceKind classify_face(const SweptMesh &d, size_t face,
                                           const FaceTables &tables,
                                           float sweep, float *motion_cache) {
  const MeshState &mesh = d.mesh;
  const size_t count = mesh.get_face_counts_data()[face];
  const uint16_t *indices =
      mesh.get_faces_data() + mesh.get_face_offsets_data()[face];
  ScratchScope scope(d.scratch);
  float *angles = d.scratch.allocate_n<float>(count);
  float nearest = 1.0f, reach = 0.0f;
  for (size_t k = 0; k < count; ++k) {
    const size_t i = indices[k];
    const math::Vector &v = mesh.vertices[i];
    // Squared distance from the polar axis.
    const float axial = v.x * v.x + v.z * v.z;
    nearest = fminf(nearest, axial);
    const float motion = vertex_reach(d.moving[i]);
    reach = fmaxf(reach, motion);
    uint8_t &state = tables.state[i];
    if (!(state & VERTEX_SEEN)) {
      const float angle = math::half_turns(v);
      const float turns = angle < 1.0f ? angle : angle - 2.0f;
      tables.turns[i] = static_cast<uint16_t>(
          hs::clamp((turns + 1.0f) / TURNS_SCALE, 0.0f, 65535.0f));
      // Only a moving vertex is captured from its column.
      bool held = false;
      if (motion > 0.0f) {
        tables.ends[i] =
            vertex_column(v, *d.canvas.column_basis(), 1.0f, d.canvas.width());
        held = motion + sweep >= sqrtf(axial) * (math::PI_F / 4.0f);
      }
      state = VERTEX_SEEN | (turns < 0.0f ? VERTEX_ARM1 : uint8_t(0)) |
              (turns - floorf(turns) < 0.5f ? VERTEX_EARLY : uint8_t(0)) |
              (held ? VERTEX_HELD : uint8_t(0));
    }
    const float turns = vertex_turns(tables, i);
    angles[k] = turns < 0.0f ? turns + 2.0f : turns;
  }
  const float horizontal = sqrtf(nearest);
  const FaceCaptureInterval interval =
      face_capture_interval(angles, count, reach, horizontal, sweep);
  if (interval.lo < 0.0f)
    return polar_kind(d, face, interval.hi, motion_cache);
  if (d.pad >= horizontal)
    return CROSS_FACE;
  // asin(x) / pi, as a half turn less acos(x) / pi.
  const float pad = 0.5f - math::fast_acos(d.pad / horizontal) / math::PI_F;
  if (interval.hi - interval.lo + 2.0f * pad >= 1.0f)
    return polar_kind(d, face, interval.hi - interval.lo + 2.0f * pad,
                      motion_cache);
  float lo = interval.lo, hi = interval.hi;
  if (hi + pad < 1.0f) {
    if (lo >= pad)
      return SHARED_FACE;
    lo += 1.0f;
    hi += 1.0f;
  }
  return lo - pad >= 0.5f && hi + pad < 1.5f ? SEAM_FACE : CROSS_FACE;
}

/**
 * @brief Scans one copy of a face alone, on one arm.
 * @param d Deformation scan.
 * @param copy The copy.
 * @param arm The arm, 0 or 1.
 */
HS_HOT_FLASH_MEMBER void scan_copy(const SweptMesh &d,
                                   const Scan::FaceArmCopy &copy, int arm) {
  Scan::ColumnBasis face_time = *d.canvas.column_basis();
  face_time.mesh.face_copy = &copy;
  ColumnBasisScope scope(d.canvas, face_time);
  d.selected[copy.face] = Scan::detail::FacePass::ARM_COPY;
  {
    TimeWindow window(d.canvas, 0, d.canvas.width() / 2, arm);
    d.callbacks.draw(d.callbacks.scan_context);
  }
  d.selected[copy.face] = Scan::detail::FacePass::SKIPPED;
#if HS_ENABLE_TEST_HOOKS
  ++mesh_capture_counters.scans;
  ++mesh_capture_counters.face_copies;
#endif
}

/**
 * @brief Scans a face reaching past the half-arm reach of a seam face, once
 *        per arm it reaches, each copy at that arm's times.
 * @param d Deformation scan.
 * @param tables Scan tables.
 * @param face Face index.
 * @param arm The arm, 0 or 1.
 * @details The face's vertices on the arm stand at their own times. Each
 * run of vertices on the other arm continues the times of the edges that
 * cross into it: its first half from the vertex before the run, its second
 * from the vertex after, so the two meet at an edge on the other arm, which
 * the copy's window leaves out. A one-vertex run whose two continuations
 * part draws twice, once for each edge. A face with no vertex on the arm
 * continues its times across the boundary nearer its middle.
 */
HS_HOT_FLASH_MEMBER void draw_arm_copy(const SweptMesh &d,
                                       const FaceTables &tables, size_t face,
                                       int arm) {
  const MeshState &mesh = d.mesh;
  const size_t count = mesh.get_face_counts_data()[face];
  const uint16_t *indices =
      mesh.get_faces_data() + mesh.get_face_offsets_data()[face];
  ScratchScope scope(d.scratch);
  // A vertex draws at most twice, once for each of its edges.
  math::Vector *vertices = d.scratch.allocate_n<math::Vector>(2 * count);
  uint16_t *order = d.scratch.allocate_n<uint16_t>(2 * count);
  float *times = d.scratch.allocate_n<float>(2 * count);
  const auto on_arm = [&](size_t k) HS_HOT_INLINE {
    return bool(tables.state[indices[k]] & VERTEX_ARM1) == (arm == 1);
  };
  const auto turns = [&](size_t k) HS_HOT_INLINE {
    return vertex_turns(tables, indices[k]);
  };
  const auto own = [&](size_t k) HS_HOT_INLINE {
    const float t = turns(k);
    return t - floorf(t);
  };
  // Time at vertex k continued from vertex j across their edge.
  const auto continued = [&](size_t j, float time, size_t k) HS_HOT_INLINE {
    const float delta = turns(k) - turns(j);
    return time + delta - 2.0f * roundf(delta * 0.5f);
  };
  size_t first = count;
  for (size_t k = 0; k < count && first == count; ++k)
    if (on_arm(k))
      first = k;
  size_t n = 0;
  const auto emit = [&](size_t k, float time) HS_HOT_INLINE {
    order[n] = static_cast<uint16_t>(k);
    times[n++] = time;
  };
  if (first == count) {
    float middle = 0.0f;
    for (size_t k = 0; k < count; ++k)
      middle += own(k);
    const float shift = middle < 0.5f * count ? 1.0f : -1.0f;
    for (size_t k = 0; k < count; ++k)
      emit(k, own(k) + shift);
  }
  for (size_t step = 0; first < count && step < count;) {
    const size_t k = (first + step) % count;
    if (on_arm(k)) {
      emit(k, own(k));
      ++step;
      continue;
    }
    size_t length = 0;
    while (!on_arm((k + length) % count))
      ++length;
    const auto run = [&](size_t r) HS_HOT_INLINE { return (k + r) % count; };
    const size_t before = (k + count - 1) % count;
    const size_t after = run(length);
    const size_t forward_count = (length + 1) / 2;
    const size_t backward_count = length - forward_count;
    float time = own(before);
    size_t previous = before;
    for (size_t r = 0; r < forward_count; ++r) {
      time = continued(previous, time, run(r));
      emit(run(r), time);
      previous = run(r);
    }
    const size_t start = n;
    n += backward_count;
    float back = own(after);
    size_t next = after;
    for (size_t r = length; r-- > length - backward_count;) {
      back = continued(next, back, run(r));
      order[start + r - (length - backward_count)] =
          static_cast<uint16_t>(run(r));
      times[start + r - (length - backward_count)] = back;
      next = run(r);
    }
    if (length == 1) {
      back = continued(after, own(after), k);
      if (fabsf(back - time) > 0.5f) {
        emit(k, back);
#if HS_ENABLE_TEST_HOOKS
        ++mesh_capture_counters.split_vertices;
#endif
      }
    }
    step += length;
  }
  uint16_t *identity = d.scratch.allocate_n<uint16_t>(n);
  {
    HS_PROFILE(mesh_copy_capture);
    for (size_t m = 0; m < n; ++m) {
      order[m] = indices[order[m]];
      vertices[m] = vertex_at(d, tables, order[m], times[m]);
      identity[m] = static_cast<uint16_t>(m);
    }
  }
  const Scan::FaceArmCopy copy{face, vertices, identity, order, n};
  scan_copy(d, copy, arm);
}

/**
 * @brief Draws a face reaching a pole one column pair at a time.
 * @param d Deformation scan, its pool at the current endpoint.
 * @param face Face index.
 * @details Each moving vertex's path is the quadratic through its closed
 * form at the frame's start and middle and its current endpoint
 * (Spacetime::Worldline::quadratic()).
 */
HS_HOT_FLASH_MEMBER void sweep_face(const SweptMesh &d, size_t face) {
  const MeshState &mesh = d.mesh;
  const size_t count = mesh.get_face_counts_data()[face];
  const uint16_t *indices =
      mesh.get_faces_data() + mesh.get_face_offsets_data()[face];
  ScratchScope scope(d.scratch);
  Spacetime::Worldline *paths =
      d.scratch.allocate_n<Spacetime::Worldline>(count);
  math::Vector middle(0.0f, 0.0f, 0.0f);
  {
    HS_PROFILE(mesh_sweep_paths);
    for (size_t k = 0; k < count; ++k) {
      const size_t i = indices[k];
      const math::Vector &end = mesh.vertices[i];
      if (d.moving[i]) {
        paths[k] = Spacetime::Worldline::quadratic(
            d.callbacks.evaluate(d.callbacks.vertex_context, i, 0.0f),
            d.callbacks.evaluate(d.callbacks.vertex_context, i, 0.5f), end);
#if HS_ENABLE_TEST_HOOKS
        mesh_capture_counters.evaluations += 2;
#endif
      } else {
        paths[k] = Spacetime::Worldline::chord(end, end);
      }
      middle = middle + paths[k].middle;
    }
  }
  const Scan::PoleFaceSweep sweep{face, paths, count, middle.y > 0.0f};
  Scan::ColumnBasis face_time = *d.canvas.column_basis();
  face_time.mesh.face_sweep = &sweep;
  ColumnBasisScope time_scope(d.canvas, face_time);
  d.callbacks.draw(d.callbacks.scan_context);
#if HS_ENABLE_TEST_HOOKS
  ++mesh_capture_counters.scans;
#endif
}

/**
 * @brief Places the captures of the moving faces one pass draws.
 * @param d Deformation scan.
 * @param tables Scan tables.
 * @param kind The pass's faces: SHARED_FACE, SEAM_FACE for an arm's, or
 *        POLAR_FACE for a half-arm window's.
 * @param arm The seam pass's arm, 0 or 1, or the polar pass's half of the
 *        arm.
 * @return Whether any face draws in the pass. A seam face draws in its arm's
 *         window, and a polar face in its half's, when its columns there meet
 *         the clip.
 * @details A polar face's vertices stand at their times nearest the middle
 * of the half (place_half()).
 */
HS_HOT_FLASH_MEMBER bool place_pass(const SweptMesh &d,
                                    const FaceTables &tables, FaceKind kind,
                                    int arm) {
  HS_PROFILE(mesh_capture);
  const MeshState &mesh = d.mesh;
  const size_t faces = mesh.num_faces();
  const int width = d.canvas.width(), half = width / 2, quarter = half / 2;
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  bool any = false;
  for (size_t f = 0; f < faces; ++f) {
    if (tables.kinds[f] != kind)
      continue;
    if (kind != SHARED_FACE) {
      const bool draw =
          kind == SEAM_FACE
              ? detail::face_columns_meet(d.columns + 2 * f, 0, half, d.clip,
                                          width, arm)
              : detail::face_columns_meet(d.columns + 2 * f, arm * quarter,
                                          (arm + 1) * quarter, d.clip, width);
      d.selected[f] = draw ? Scan::detail::FacePass::ARM_COPY
                           : Scan::detail::FacePass::SKIPPED;
      if (!draw)
        continue;
#if HS_ENABLE_TEST_HOOKS
      mesh_capture_counters.seam_faces += kind == SEAM_FACE;
#endif
    }
    any = true;
    for (size_t k = 0; k < counts[f]; ++k) {
      const size_t i = indices[offsets[f] + k];
      if (kind == POLAR_FACE)
        place_half(d, tables, i, arm);
      else
        place_vertex(d, tables, i,
                     kind == SEAM_FACE &&
                         bool(tables.state[i] & VERTEX_ARM1) != (arm == 1));
    }
  }
  return any;
}

} // namespace

HS_HOT_FLASH_MEMBER void scan_swept_faces(const SweptMesh &d) {
  MeshState &mesh = d.mesh;
  const size_t faces = mesh.num_faces(), vertices = mesh.vertices.size();
  const int width = d.canvas.width(), half = width / 2;
  Mesh::check_face_index_domain(mesh.get_faces_data(), mesh.get_faces_size(),
                                vertices);
  ScratchScope guard(d.scratch);
  const FaceTables tables{d.scratch.allocate_n<uint8_t>(faces),
                          d.scratch.allocate_n<uint8_t>(vertices),
                          d.scratch.allocate_n<uint16_t>(vertices),
                          d.scratch.allocate_n<uint16_t>(vertices)};
  const float sweep = d.canvas.column_basis()->sweep_arc;
  std::fill_n(tables.state, vertices, uint8_t(0));
#if HS_ENABLE_TEST_HOOKS
  ++mesh_capture_counters.deformations;
  for (size_t i = 0; i < vertices; ++i)
    mesh_capture_counters.vertices += d.moving[i] != 0;
#endif
  bool shared = false, polar = false;
  {
    HS_PROFILE(mesh_classify);
    ScratchScope classify_scope(d.scratch);
    float *motion_cache = d.scratch.allocate_n<float>(vertices);
    std::fill_n(motion_cache, vertices, -1.0f);
    for (size_t f = 0; f < faces; ++f) {
      if (!mesh_visible(d.canvas.column_basis()->mesh.visible_faces, f)) {
        tables.kinds[f] = CULLED_FACE;
        d.selected[f] = Scan::detail::FacePass::SKIPPED;
        continue;
      }
      const uint8_t kind =
          d.columns[2 * f] == HELD_FACE_COLUMN
              ? HELD_FACE
              : classify_face(d, f, tables, sweep, motion_cache);
      tables.kinds[f] = kind;
      shared |= kind == HELD_FACE || kind == SHARED_FACE;
      polar |= kind == POLAR_FACE;
      d.selected[f] = kind == HELD_FACE || kind == SHARED_FACE
                          ? Scan::detail::FacePass::SHARED
                          : Scan::detail::FacePass::SKIPPED;
#if HS_ENABLE_TEST_HOOKS
      mesh_capture_counters.shared_faces += kind == SHARED_FACE;
      mesh_capture_counters.cross_faces += kind == CROSS_FACE;
      mesh_capture_counters.polar_faces +=
          kind == POLAR_FACE || kind == SWEPT_FACE;
      mesh_capture_counters.swept_faces += kind == SWEPT_FACE;
#endif
    }
  }
  // Faces too fast for one placement sweep from the endpoint pool, before
  // any capture replaces it.
  for (size_t f = 0; f < faces; ++f)
    if (tables.kinds[f] == SWEPT_FACE)
      sweep_face(d, f);
  // Faces reaching a pole draw through the two halves of each arm.
  for (int part = 0; polar && part < 2; ++part) {
    if (place_pass(d, tables, POLAR_FACE, part)) {
      TimeWindow window(d.canvas, part * (half / 2), (part + 1) * (half / 2));
      d.callbacks.draw(d.callbacks.scan_context);
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.scans;
#endif
    }
  }
  for (size_t f = 0; polar && f < faces; ++f)
    if (tables.kinds[f] == POLAR_FACE)
      d.selected[f] = Scan::detail::FacePass::SKIPPED;
  if (shared) {
    place_pass(d, tables, SHARED_FACE, 0);
    d.callbacks.draw(d.callbacks.scan_context);
#if HS_ENABLE_TEST_HOOKS
    ++mesh_capture_counters.scans;
#endif
  }
  for (const int arm : {0, 1}) {
    if (place_pass(d, tables, SEAM_FACE, arm)) {
      TimeWindow window(d.canvas, 0, half, arm);
      d.callbacks.draw(d.callbacks.scan_context);
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.scans;
#endif
    }
    for (size_t f = 0; f < faces; ++f) {
      if (tables.kinds[f] == SEAM_FACE)
        d.selected[f] = Scan::detail::FacePass::SKIPPED;
      if (tables.kinds[f] == CROSS_FACE &&
          detail::face_columns_meet(d.columns + 2 * f, 0, half, d.clip, width,
                                    arm))
        draw_arm_copy(d, tables, f, arm);
    }
  }
}

} // namespace Scan::detail
