/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file mesh_capture.h
 * @brief Face visibility and column-time capture of deforming mesh vertices. */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include "core/memory.h"
#include "core/mesh/mesh_state.h"
#include "core/render/canvas.h"
#include "core/render/sdf/face_classes.h"
#include "core/render/scan/face_raster.h"
#include "core/render/scan/swept_face.h"
#include "core/render/spacetime/crossing.h"

namespace Scan {

/**
 * @brief Traps unless every flat face index addresses a vertex in the pool.
 * @param faces Flat per-face vertex index array.
 * @param num_indices Entries in @p faces.
 * @param num_verts Vertices the mesh carries.
 * @details SDF::PreparedFace indexes the vertex pool through an unchecked std::span.
 */
HS_O3_FN HS_NOINLINE_NOCLONE inline void
check_face_index_domain(const uint16_t *faces, size_t num_indices,
                        size_t num_verts) {
  for (size_t k = 0; k < num_indices; ++k)
    HS_CHECK(static_cast<size_t>(faces[k]) < num_verts,
             "mesh face index exceeds the vertex pool");
}

/**
 * @brief Whether an indexed element belongs to an optional visibility bitset.
 * @param bits Visibility bits, or null to admit every element.
 * @param i Element index.
 * @return True when the element is admitted.
 */
HS_HOT_INLINE inline bool mesh_visible(const uint8_t *bits, size_t i) {
  return !bits || (bits[i / 8] & (1u << (i % 8)));
}

/// Seeing-column entry of a face that holds over the frame.
inline constexpr uint16_t HELD_FACE_COLUMN = 0xFFFF;

/**
 * @brief Canvas columns that can see a cap over one frame.
 * @tparam W Canvas width in pixels.
 * @param axis Cap axis.
 * @param radius Cap radius (radians), or UNBOUNDED_FACE_CAP.
 * @param reach Largest arc the cap or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the columns as detail::face_lit_columns() gives them.
 */
template <int W>
HS_HOT_INLINE void cap_lit_columns(const math::Vector &axis, float radius,
                                   float reach, uint16_t *columns) {
  int lo, hi;
  if (radius == UNBOUNDED_FACE_CAP ||
      !lit_span<W>(axis, radius, reach, lo, hi) || hi - lo >= W - 1)
    return;
  columns[0] = static_cast<uint16_t>(math::wrap(lo, W));
  columns[1] = static_cast<uint16_t>(hi - lo);
}

namespace detail {
/** @brief Checks offset count and audits the solid mesh scan's face spans.
 * @param mesh Mesh to scan. */
HS_HOT_INLINE inline void check_face_offsets(const MeshState &mesh) {
  HS_CHECK(mesh.get_face_offsets_size() == mesh.num_faces(),
           "solid mesh scan requires one face offset per face");
#if HS_ENABLE_STRUCTURAL_AUDITS
  const size_t indices = mesh.get_faces_size();
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  for (size_t f = 0; f < mesh.num_faces(); ++f)
    HS_AUDIT_CHECK(static_cast<size_t>(offsets[f]) + counts[f] <= indices,
                   "mesh face span exceeds face index array");
#endif
}

/**
 * @brief Cap about a face's vertex mean that holds every vertex.
 * @param vertices Mesh vertex pool.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param axis Receives the cap axis, when the vertex mean gives a direction.
 * @return The cap radius (radians), or UNBOUNDED_FACE_CAP for a cap past a
 *         hemisphere or a vertex mean too short to give a direction.
 */
HS_HOT_INLINE inline float face_cap(const math::Vector *vertices,
                                    const uint16_t *indices, size_t count,
                                    math::Vector &axis) {
  const float cosine = face_cap_cosine(vertices, indices, count, axis);
  return cosine > 0.0f ? cosine_arc(cosine) : UNBOUNDED_FACE_CAP;
}
} // namespace detail

/**
 * @brief Selects the model-space faces and vertices a board can see this frame.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam Camera Model-to-endpoint rotation callable.
 * @param clip Board render clip.
 * @param mesh Undeformed model mesh.
 * @param camera Model-to-endpoint rotation.
 * @param reach Deformation displacement and capture-motion bound in radians.
 * @param arena Frame arena holding the selection bits.
 * @param time Column basis receiving the selection bits.
 */
template <int W, int H, typename Camera>
HS_HOT_FLASH_MEMBER void
select_mesh_visibility(const ClipRegion &clip, const MeshState &mesh,
                       const Camera &camera, float reach, Arena &arena,
                       ColumnBasis &time) {
  detail::check_face_offsets(mesh);
  check_face_index_domain(mesh.get_faces_data(), mesh.get_faces_size(),
                          mesh.vertices.size());
  HS_PROFILE(mesh_visibility);
#if HS_ENABLE_TEST_HOOKS
  if (face_audit_reference)
    return;
  visibility_counts[0] += mesh.num_faces();
  visibility_counts[2] += mesh.vertices.size();
  visibility_audit_reach = reach;
#endif
  const size_t faces = mesh.num_faces(), vertices = mesh.vertices.size();
  auto *face_bits = arena.allocate_n<uint8_t>((faces + 7) / 8);
  auto *vertex_bits = arena.allocate_n<uint8_t>((vertices + 7) / 8);
  std::fill_n(face_bits, (faces + 7) / 8, uint8_t(0));
  std::fill_n(vertex_bits, (vertices + 7) / 8, uint8_t(0));
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const float pad =
      reach + time.sweep_arc + SDF::BOUNDS_MARGIN_WIDE + SDF::MOVED_AXIS_SLACK;
  const FaceClipPlanes<W, H> planes(clip, pad);
  float max_length_sq = 0.0f;
  for (const auto &v : mesh.vertices)
    max_length_sq = fmaxf(max_length_sq, math::dot(v, v));
  const float inverse_max_length =
      max_length_sq > 0.0f ? 1.0f / sqrtf(max_length_sq) : 0.0f;
  for (size_t f = 0; f < faces; ++f) {
    const size_t count = mesh.get_face_counts_data()[f];
    const uint16_t *indices =
        mesh.get_faces_data() + mesh.get_face_offsets_data()[f];
    math::Vector axis;
    if (pad < 0.5f * math::PI_F &&
        face_cap_axis(mesh.vertices.data(), indices, count, axis)) {
      float cosine = 1.0f;
      for (size_t k = 0; k < count; ++k)
        cosine = fminf(cosine, math::dot(axis, mesh.vertices[indices[k]]) *
                                   inverse_max_length);
      if (cosine > 0.0f && planes.rejects(camera(axis), cosine))
        continue;
    }
#if HS_ENABLE_TEST_HOOKS
    ++visibility_counts[1];
#endif
    face_bits[f / 8] |= 1u << (f % 8);
    for (size_t k = 0; k < count; ++k)
      vertex_bits[indices[k] / 8] |= 1u << (indices[k] % 8);
  }
#if HS_ENABLE_TEST_HOOKS
  for (size_t i = 0; i < vertices; ++i)
    visibility_counts[3] += mesh_visible(vertex_bits, i);
#endif
  time.mesh.visibility_face_count = faces;
  time.mesh.visibility_vertex_count = vertices;
  time.mesh.visible_faces = face_bits;
  time.mesh.visible_vertices = vertex_bits;
}

namespace detail {
/**
 * @brief Canvas columns that can see a face over one frame.
 * @tparam W Canvas width in pixels.
 * @param vertices Mesh vertex pool, at the current endpoint.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param reach Largest arc the face or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the first seeing column, in [0, W), and how many
 *        more follow it, wrapping past the last column.
 * @param radius Receives the face's cap radius for face_cap_columns(), when
 *        not null.
 * @details Judged by detail::face_cap() widened by @p reach; a cap past a hemisphere
 * or holding a pole is seen by every column.
 */
template <int W>
HS_FLASH_MEMBER void face_lit_columns(const math::Vector *vertices,
                                      const uint16_t *indices, size_t count,
                                      float reach, uint16_t *columns,
                                      float *radius = nullptr) {
  columns[0] = 0;
  columns[1] = W - 1;
  math::Vector axis;
  const float cap = detail::face_cap(vertices, indices, count, axis);
  if (radius)
    *radius = cap;
  cap_lit_columns<W>(axis, cap, reach, columns);
}
} // namespace detail

/**
 * @brief detail::face_lit_columns() of a face whose cap radius it recorded.
 * @tparam W Canvas width in pixels.
 * @param vertices Mesh vertex pool, as when the radius was recorded.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param radius The recorded cap radius.
 * @param reach Largest arc the face or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the columns as detail::face_lit_columns() gives them.
 */
template <int W>
HS_FLASH_MEMBER void
face_cap_columns(const math::Vector *vertices, const uint16_t *indices,
                 size_t count, float radius, float reach, uint16_t *columns) {
  columns[0] = 0;
  columns[1] = W - 1;
  math::Vector axis;
  if (face_cap_axis(vertices, indices, count, axis))
    cap_lit_columns<W>(axis, radius, reach, columns);
}

namespace detail {
/**
 * @brief Whether a face's seeing columns meet a capture window's columns
 *        inside the render clip.
 * @param face The face's seeing columns, as detail::face_lit_columns() gives them.
 * @param begin The window's first half-turn column.
 * @param end One past the window's last half-turn column.
 * @param clip The render clip's column runs.
 * @param width Canvas width in pixels.
 * @param arm_only The one arm, 0 or 1, to test, or -1 for both.
 * @return True when a column of the window, on a tested arm, sees the face
 *         and the clip admits it.
 */
HS_FLASH_MEMBER inline bool face_columns_meet(const uint16_t *face, int begin,
                                              int end, const ClipPieces &clip,
                                              int width, int arm_only = -1) {
  const int half = width / 2;
  const int start = face[0];
  const int stop = start + face[1] + 1;
  // The face's columns past the last one wrap to the first.
  for (const int shift : {0, width})
    for (const int arm : {0, half}) {
      if (arm_only >= 0 && arm != arm_only * half)
        continue;
      const int lo = std::max(start - shift, begin + arm);
      const int hi = std::min(stop - shift, end + arm);
      for (int i = 0; i < clip.count; ++i)
        if (std::max(lo, clip.lo[i]) < std::min(hi, clip.hi[i]))
          return true;
    }
  return false;
}
/**
 * @brief detail::face_columns_meet() under a render clip.
 * @tparam W Canvas width in pixels.
 * @param face The face's seeing columns, as detail::face_lit_columns() gives them.
 * @param begin The window's first half-turn column.
 * @param end One past the window's last half-turn column.
 * @param xc Render clip.
 * @return True when a column of the window, on either arm, sees the face and
 *         the clip admits it.
 */
template <int W>
HS_FLASH_MEMBER bool face_columns_meet(const uint16_t *face, int begin, int end,
                                       ClipRegion::XClip xc) {
  return detail::face_columns_meet(face, begin, end, clip_pieces<W>(xc), W);
}
} // namespace detail

/// Largest held vertex motion over a frame, in pixels of arc.
inline constexpr float HELD_VERTEX_PIXELS = 0.5f;

/// Arc a vertex moves over a frame, radians, up to which a deformation scan
/// holds it at the current endpoint.
template <int W>
inline constexpr float HELD_VERTEX_ARC =
    HELD_VERTEX_PIXELS * math::TWO_PI_F / W;

/// Moving-vertex flag: the vertex's path is captured.
inline constexpr uint8_t VERTEX_MOVES = 1;
/// Moving-vertex flags: shift of the vertex's reach level (vertex_flags()).
inline constexpr int VERTEX_REACH_SHIFT = 2;
/// Arc, radians, the lowest reach level stands for.
inline constexpr float VERTEX_REACH_UNIT = 1.0f / 1024.0f;
/// The arc each reach level stands for, 2^(1/3) apart, the top one
/// unbounded.
inline constexpr auto VERTEX_REACH_LEVELS = [] {
  std::array<float, 32> levels{};
  // 2^(1/3) to float precision.
  constexpr double STEP = 1.2599210498948732;
  double arc = VERTEX_REACH_UNIT;
  for (float &level : levels) {
    level = static_cast<float>(arc);
    arc *= STEP;
  }
  // The top level stands for any reach past the others.
  levels.back() = FLT_MAX;
  return levels;
}();

/**
 * @brief Flags of a moving vertex.
 * @param reach Bound on its distance from its endpoint over the frame
 *        (radians).
 * @return VERTEX_MOVES with the lowest reach level at least @p reach.
 */
HS_HOT_INLINE inline uint8_t vertex_flags(float reach) {
  // Binary search of the ascending levels; the top one bounds any reach.
  int level = 0;
  for (int step = 16; step > 0; step >>= 1)
    if (VERTEX_REACH_LEVELS[level + step - 1] < reach)
      level += step;
  return static_cast<uint8_t>(VERTEX_MOVES | (level << VERTEX_REACH_SHIFT));
}

/**
 * @brief A vertex's reach.
 * @param flags Its flags, as vertex_flags() gives them, or 0.
 * @return Radians: the bound its level stands for, 0 for a held vertex.
 */
HS_HOT_INLINE inline float vertex_reach(uint8_t flags) {
  return flags & VERTEX_MOVES
             ? VERTEX_REACH_LEVELS[(flags >> VERTEX_REACH_SHIFT) & 31]
             : 0.0f;
}

/**
 * @brief Bound on the arc between two positions of a vertex.
 * @param a,b The positions, any nonzero length.
 * @return Radians, at least the angle between them.
 */
HS_HOT_INLINE inline float vertex_arc(const math::Vector &a,
                                      const math::Vector &b) {
  const float norms = math::dot(a, a) * math::dot(b, b);
  return cosine_arc(norms > 0.0f ? math::dot(a, b) / sqrtf(norms) : 1.0f);
}

/// Vertex motion over one frame, in pixels of arc, past which the frame is a
/// jump and draws at its current endpoint.
inline constexpr float MAX_VERTEX_MOTION_PIXELS = 8.0f;

/**
 * @brief Largest arc any vertex moves between two endpoints.
 * @param previous Vertices at the previous endpoint.
 * @param current The same vertices at the current endpoint.
 * @param count Number of vertices.
 * @return The largest angle between a vertex's two positions (radians).
 */
HS_FLASH_MEMBER inline float max_vertex_arc(const math::Vector *previous,
                                            const math::Vector *current,
                                            size_t count) {
  float min_cos = 1.0f;
  for (size_t i = 0; i < count; ++i) {
    const float norms =
        math::dot(previous[i], previous[i]) * math::dot(current[i], current[i]);
    if (norms > 0.0f)
      min_cos =
          fminf(min_cos, math::dot(previous[i], current[i]) / sqrtf(norms));
  }
  return acosf(hs::clamp(min_cos, -1.0f, 1.0f));
}

/**
 * @brief Whether a frame's mesh motion uses capture-time deformation.
 * @tparam W Display width.
 * @param arc Largest vertex step over the frame, in radians.
 * @return True with temporal sampling enabled for a non-held, non-jump step.
 */
template <int W> HS_FLASH_MEMBER bool mesh_motion_active(float arc) {
  const float pixels = arc * (W / math::TWO_PI_F);
  return !Spacetime::flat_surface() && pixels > HELD_VERTEX_PIXELS &&
         pixels <= MAX_VERTEX_MOTION_PIXELS;
}

/** @brief Angular bounds margin for a mesh scan.
 * @tparam W Canvas width in pixels.
 * @param sweep Column-basis sweep arc in radians.
 * @return Base face margin widened by the sweep. */
template <int W> HS_HOT_INLINE float scan_bounds_margin(float sweep) {
  return std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W) + sweep;
}

/**
 * @brief Arc a mesh scan pads a face's rows by past its vertices.
 * @tparam W Canvas width in pixels.
 * @param sweep Largest arc a column's basis moves over the frame (radians).
 * @return Scan::Mesh's bounds margin with a column of slack.
 */
template <int W> HS_HOT_INLINE float scan_row_pad(float sweep) {
  return scan_bounds_margin<W>(sweep) + math::TWO_PI_F / W;
}

/** @brief Vertex heights past which a moving face reaches no clip row. */
struct ClipRows {
  float south = -2.0f; ///< Reached by no face whose vertices all stand below.
  float north = 2.0f;  ///< Reached by no face whose vertices all stand above.

  /**
   * @brief The heights of a render clip's rows.
   * @tparam W Canvas width in pixels.
   * @param clip Render clip.
   * @param reach Bound on how far a placed vertex strays from the current
   *        endpoint over the frame (radians).
   * @param sweep Largest arc a column's basis moves over the frame (radians).
   * @return The heights; each lies past the equator by at least the padded
   *         reach, so a face beyond it keeps every placement in its
   *         hemisphere, where its edges bulge away from the clip.
   */
  template <int W>
  HS_FLASH_MEMBER static ClipRows of(const ClipRegion &clip, float reach,
                                     float sweep) {
    const math::LatitudeGeometry geometry(clip.h);
    const float pad = reach + scan_row_pad<W>(sweep);
    ClipRows rows;
    if (clip.render_y_end() < clip.h) {
      const float phi =
          std::max(geometry.row_to_phi(static_cast<float>(clip.render_y_end())),
                   0.5f * math::PI_F) +
          pad;
      if (phi < math::PI_F)
        rows.south = cosf(phi);
    }
    if (clip.render_y_start() > 0) {
      const float phi = std::min(geometry.row_to_phi(static_cast<float>(
                                     clip.render_y_start() - 1)),
                                 0.5f * math::PI_F) -
                        pad;
      if (phi > 0.0f)
        rows.north = cosf(phi);
    }
    return rows;
  }

  /**
   * @brief Whether a face can reach a clip row.
   * @param vertices Mesh vertex pool, at the current endpoint.
   * @param indices The face's vertex indices.
   * @param count Number of face vertices.
   * @return False when every vertex stands below `south` or above `north`.
   */
  HS_HOT_INLINE bool meet(const math::Vector *vertices, const uint16_t *indices,
                          size_t count) const {
    float low = 2.0f, high = -2.0f;
    for (size_t k = 0; k < count; ++k) {
      const float y = vertices[indices[k]].y;
      low = fminf(low, y);
      high = fmaxf(high, y);
    }
    return high >= south && low <= north;
  }
};

/**
 * @brief Recomputes the seeing columns of a deformation scan's moving faces.
 * @tparam W Canvas width in pixels.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The scanned mesh, at the current endpoint.
 * @param reach Bound on how far a placed vertex strays from the current
 *        endpoint over the frame (radians).
 * @param columns Per face two entries, HELD_FACE_COLUMN first for a face
 *        that holds; the others receive the columns that can see the face,
 *        as detail::face_lit_columns() gives them, or HELD_FACE_COLUMN and 0 when no
 *        placement of the face reaches a row the render clip admits.
 * @param radii Per-face cap radii pick_moving_faces() recorded from the same
 *        vertices, or null to measure the caps.
 * @details Call after pick_moving_faces(), optionally after drop_held_faces(),
 * when tightening the motion reach. Mesh vertices and
 * topology must match the selection. Tables remain caller-owned and are borrowed
 * only for this call; arena-backed tables must outlive it.
 */
template <int W>
HS_FLASH_MEMBER void
moving_face_columns(const Canvas &canvas, const MeshState &mesh, float reach,
                    uint16_t *columns, const float *radii = nullptr) {
  detail::check_face_offsets(mesh);
  const float sweep =
      canvas.column_basis() ? canvas.column_basis()->sweep_arc : 0.0f;
  if (canvas.column_basis())
    check_mesh_pass(canvas.column_basis()->mesh, mesh.num_faces(),
                    mesh.vertices.size());
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  const ClipRows rows = ClipRows::of<W>(canvas.clip(), reach, sweep);
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    uint16_t *face = columns + 2 * f;
    if (face[0] == HELD_FACE_COLUMN)
      continue;
    if (!rows.meet(mesh.vertices.data(), indices + offsets[f], counts[f])) {
      face[0] = HELD_FACE_COLUMN;
      face[1] = 0;
      continue;
    }
    if (radii)
      face_cap_columns<W>(mesh.vertices.data(), indices + offsets[f], counts[f],
                          radii[f], sweep + reach, face);
    else
      detail::face_lit_columns<W>(mesh.vertices.data(), indices + offsets[f],
                                  counts[f], sweep + reach, face);
  }
}

/**
 * @brief Picks the faces whose vertices move during a scan.
 * @tparam W Canvas width in pixels.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The scanned mesh, at the current endpoint.
 * @param moving_vertices Nonzero per moving vertex on entry; receives 1 on
 *        every vertex of a picked face, held ones included, and 0 elsewhere.
 * @param reach Bound on how far a placed vertex strays from the current
 *        endpoint over the frame (radians).
 * @param columns Receives per face two entries: the columns that can see a
 *        picked face, as detail::face_lit_columns() gives them, or HELD_FACE_COLUMN
 *        and 0 for the others.
 * @param radii When not null, receives the cap radius of each face holding a
 *        moving vertex that can reach a clip row, for moving_face_columns().
 * @return The number of picked faces.
 * @details A face is picked when it holds a moving vertex and a column and a
 * row the render clip admits can see it; the others draw once at the current
 * endpoint. Call before drop_held_faces(); optionally
 * tighten the reach with moving_face_columns() before deform_mesh(). The caller
 * owns the moving-vertex flags, column pairs and optional radii; the function
 * borrows them only for this call. Keep flags and columns alive through
 * deform_mesh(), and radii through moving_face_columns(), including their arena
 * scopes. The mesh and installed column-time state must describe the same draw.
 */
template <int W>
HS_FLASH_MEMBER size_t pick_moving_faces(const Canvas &canvas,
                                         const MeshState &mesh,
                                         uint8_t *moving_vertices, float reach,
                                         uint16_t *columns,
                                         float *radii = nullptr) {
  detail::check_face_offsets(mesh);
  check_face_index_domain(mesh.get_faces_data(), mesh.get_faces_size(),
                          mesh.vertices.size());
  if (canvas.column_basis())
    check_mesh_pass(canvas.column_basis()->mesh, mesh.num_faces(),
                    mesh.vertices.size());
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  const size_t faces = mesh.num_faces();
  const float sweep =
      canvas.column_basis() ? canvas.column_basis()->sweep_arc : 0.0f;
  const ClipRows rows = ClipRows::of<W>(canvas.clip(), reach, sweep);
  for (size_t f = 0; f < faces; ++f) {
    uint16_t *face = columns + 2 * f;
    uint8_t moves = 0;
    for (size_t k = 0; k < counts[f]; ++k)
      moves |= moving_vertices[indices[offsets[f] + k]];
    face[0] = HELD_FACE_COLUMN;
    face[1] = 0;
    if (!moves ||
        (canvas.column_basis() &&
         !mesh_visible(canvas.column_basis()->mesh.visible_faces, f)) ||
        !rows.meet(mesh.vertices.data(), indices + offsets[f], counts[f]))
      continue;
    detail::face_lit_columns<W>(mesh.vertices.data(), indices + offsets[f],
                                counts[f], sweep + reach, face,
                                radii ? radii + f : nullptr);
  }
  std::fill_n(moving_vertices, mesh.vertices.size(), uint8_t(0));
  size_t picked = 0;
  for (size_t f = 0; f < faces; ++f) {
    uint16_t *face = columns + 2 * f;
    if (face[0] == HELD_FACE_COLUMN)
      continue;
    if (!detail::face_columns_meet<W>(face, 0, W / 2, canvas.render_x_clip())) {
      face[0] = HELD_FACE_COLUMN;
      face[1] = 0;
      continue;
    }
    ++picked;
    for (size_t k = 0; k < counts[f]; ++k)
      moving_vertices[indices[offsets[f] + k]] = 1;
  }
  return picked;
}

/**
 * @brief Returns to holding the picked faces none of whose vertices moves.
 * @param mesh The scanned mesh.
 * @param moving Nonzero per vertex that moves.
 * @param columns The faces' seeing columns, as pick_moving_faces() gives them;
 *        a dropped face's become HELD_FACE_COLUMN and 0.
 * @return The number of faces still picked.
 */
HS_FLASH_MEMBER inline size_t drop_held_faces(const MeshState &mesh,
                                              const uint8_t *moving,
                                              uint16_t *columns) {
  detail::check_face_offsets(mesh);
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  size_t picked = 0;
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    if (columns[2 * f] == HELD_FACE_COLUMN)
      continue;
    uint8_t moves = 0;
    for (size_t k = 0; k < counts[f]; ++k)
      moves |= moving[indices[offsets[f] + k]];
    if (moves) {
      ++picked;
    } else {
      columns[2 * f] = HELD_FACE_COLUMN;
      columns[2 * f + 1] = 0;
    }
  }
  return picked;
}

/**
 * @brief Unwraps a capture fraction around a window's centre.
 * @param u Periodic arm fraction.
 * @param centre Window centre in arm fractions.
 * @return Capture fraction within half a frame of @p centre; past 1 or before
 *         0 for a vertex the next or previous arm's columns see.
 */
HS_HOT_INLINE inline float unwrap_capture(float u, float centre) {
  return u + roundf(centre - u);
}

/**
 * @brief Projects a vertex into the frame's moving column basis.
 * @param p Vertex in endpoint camera space.
 * @param time Frame's column basis.
 * @param u Fraction selecting the local column basis.
 * @param width Display width.
 * @return Periodic arm fraction of the projected vertex, to math::fast_atan2's
 *         accuracy.
 */
HS_HOT_INLINE inline float capture_fraction(const math::Vector &p,
                                            const Scan::ColumnBasis &time,
                                            float u, int width) {
  const int x = hs::clamp(static_cast<int>(u * (width / 2)), 0, width / 2 - 1);
  const math::Vector tangent = math::cross(time.azimuth[x], time.up[x]);
  const float fraction =
      Spacetime::surface_time(x, width) +
      math::fast_atan2(math::dot(p, tangent), math::dot(p, time.azimuth[x])) /
          math::PI_F;
  return math::wrap_t(fraction);
}

/**
 * @brief A vertex's column in the frame's moving column basis, packed.
 * @param p Vertex in endpoint camera space.
 * @param time Frame's column basis.
 * @param u Fraction selecting the local column basis.
 * @param width Display width.
 * @return capture_fraction() in 65536ths of an arm, wrapping at 1.
 */
HS_HOT_INLINE inline uint16_t vertex_column(const math::Vector &p,
                                            const Scan::ColumnBasis &time,
                                            float u, int width) {
  return static_cast<uint16_t>(static_cast<uint32_t>(
      capture_fraction(p, time, u, width) * 65536.0f + 0.5f));
}

#if HS_ENABLE_TEST_HOOKS
/** @brief Host work counts for moving mesh capture. */
struct MeshCaptureCounters {
  size_t deformations = 0;   ///< Deformation scans that classified faces.
  size_t vertices = 0;       ///< Moving vertices admitted by the segment clip.
  size_t evaluations = 0;    ///< Vertex callback evaluations.
  size_t captures = 0;       ///< Vertex capture solves.
  size_t shared_faces = 0;   ///< Moving faces scanned once in the shared pass.
  size_t seam_faces = 0;     ///< Seam-face copies in the arms' shared passes.
  size_t cross_faces = 0;    ///< Faces drawn in copies of their own per arm.
  size_t polar_faces = 0;    ///< Faces reaching a pole.
  size_t swept_faces = 0;    ///< Polar faces swept column pair by column pair.
  size_t sweep_columns = 0;  ///< Column pairs those sweeps scanned.
  size_t sweep_rows = 0;     ///< Rows those column pairs scanned, per arm.
  size_t sweep_tests = 0;    ///< Column pairs the sweeps tested.
  size_t sweep_pixels = 0;   ///< Probes of those rows within a pixel of a face.
  size_t sweep_probes = 0;   ///< Probes those sweeps made.
  size_t face_copies = 0;    ///< Faces scanned alone.
  size_t split_vertices = 0; ///< Vertices a face copy draws twice.
  size_t scans = 0;          ///< Mesh scan invocations.
  size_t face_setups = 0;    ///< Faces those scans constructed.
};
/** @brief Instrumentation for the current mesh measurement. */
inline MeshCaptureCounters mesh_capture_counters;
#endif

/**
 * @brief Evaluates one moving vertex at the time its column sees it.
 * @tparam VertexFn Callable (float u) returning the vertex at that fraction.
 * @param start The vertex's column at the frame's start, as vertex_column()
 *        gives it at u = 0.
 * @param end Its column at the frame's end, as vertex_column() gives it at
 *        u = 1.
 * @param centre Arm fraction the capture's branch lies within half a frame of.
 * @param width Display width.
 * @param vertex Closed-form vertex path, continuing the frame's motion past
 *        either endpoint.
 * @return Captured vertex.
 * @details The column a vertex stands in runs linearly between its columns
 * at the frame's two ends; the capture time is that track's arm pass
 * (Spacetime::crossings()) on the branch whose column time lies within half
 * a frame of @p centre, and the path is evaluated once, there. A vertex the
 * neighbouring arm's columns see is captured at that column's time, past 1
 * or before 0.
 */
template <typename VertexFn>
HS_HOT_INLINE math::Vector capture_vertex(uint16_t start, uint16_t end,
                                          float centre, int width,
                                          VertexFn &&vertex) {
#if HS_ENABLE_TEST_HOOKS
  ++mesh_capture_counters.captures;
  ++mesh_capture_counters.evaluations;
#endif
  constexpr float SCALE = 1.0f / 65536.0f;
  const float first = unwrap_capture(start * SCALE, centre);
  const float travel = unwrap_capture(end * SCALE, first) - first;
  // The track over the frames on either side, so a pass past either end of
  // this one is found too; its first-arm branch is the column's own time.
  const float half = width / 2.0f;
  float u = first;
  Spacetime::crossings(width, (first - travel) * half,
                       (first + 2.0f * travel) * half, -1.0f, 2.0f,
                       Spacetime::FrameEnds::NEITHER,
                       [&](const Spacetime::Crossing &pass) HS_HOT_INLINE {
                         if (pass.arm == 0)
                           u = pass.time;
                       });
  return vertex(u);
}

/** @brief Contiguous arm-time interval of a face. */
struct FaceCaptureInterval {
  float lo; ///< First vertex fraction; negative for polar faces.
  /// Last unwrapped vertex fraction, or an upper bound on a polar face's span.
  float hi;
};

/**
 * @brief Bounds a face after unwrapping azimuths about its first vertex.
 * @param angles Signed vertex azimuths in half turns, one per face vertex.
 * @param count Number of face vertices.
 * @param reach Its fastest vertex's reach (vertex_reach()); 0 when none
 *        moves.
 * @param nearest Its nearest vertex's distance from the polar axis.
 * @param sweep The column basis's sweep over the frame.
 * @return Contiguous interval, or a negative lower bound for a face that
 *         winds about a pole, spans a half turn or reaches the pole: its
 *         fastest vertex's reach, with the basis's sweep, reaches half
 *         a right angle about the polar axis from one of its vertices.
 *         For polar faces the upper value bounds the azimuth span in half turns.
 */
HS_HOT_FLASH_MEMBER FaceCaptureInterval face_capture_interval(
    const float *angles, size_t count, float reach, float nearest, float sweep);

namespace detail {

/** @brief Borrowed face-capture callbacks with flash-resident dispatch. */
struct SweptMeshCallbacks {
  void *vertex_context; ///< Borrowed vertex evaluator.
  void *scan_context;   ///< Borrowed mesh scanner.
  math::Vector (*evaluate)(void *, size_t, float); ///< Vertex dispatch.
  void (*draw)(void *);                            ///< Scan dispatch.

  /**
   * @brief Borrows the callbacks for one deformation scan.
   * @tparam VertexFn Vertex path evaluator.
   * @tparam ScanFn Mesh scanner.
   * @param vertex Vertex path evaluator.
   * @param scan Mesh scanner.
   */
  template <typename VertexFn, typename ScanFn>
  SweptMeshCallbacks(VertexFn &vertex, ScanFn &scan)
      : vertex_context(const_cast<void *>(static_cast<const void *>(&vertex))),
        scan_context(const_cast<void *>(static_cast<const void *>(&scan))),
        evaluate(evaluate_vertex<VertexFn>), draw(scan_faces<ScanFn>) {}

  /**
   * @brief Evaluates a borrowed vertex path.
   * @tparam VertexFn Vertex path evaluator.
   * @param context Borrowed evaluator.
   * @param i Vertex index.
   * @param u Capture fraction.
   * @return Evaluated vertex.
   */
  template <typename VertexFn>
  HS_HOT_FLASH_MEMBER static math::Vector evaluate_vertex(void *context,
                                                          size_t i, float u) {
    return (*static_cast<VertexFn *>(context))(i, u);
  }

  /**
   * @brief Invokes a borrowed mesh scanner.
   * @tparam ScanFn Mesh scanner.
   * @param context Borrowed scanner.
   */
  template <typename ScanFn>
  HS_HOT_FLASH_MEMBER static void scan_faces(void *context) {
    (*static_cast<ScanFn *>(context))();
  }
};

/** @brief One deformation scan's frame state, shared by its face passes. */
struct SweptMesh {
  Canvas &canvas;          ///< Canvas carrying the scan's column-time state.
  MeshState &mesh;         ///< Vertex pool, overwritten by the captures.
  const uint16_t *columns; ///< Per-face seeing columns.
  const uint8_t *moving;   ///< vertex_flags() per vertex; 0 held.
  /// Moving vertices' columns at the frame's start (vertex_column()).
  const uint16_t *start;
  uint8_t *selected; ///< The scan's face selection table.
  Arena &scratch;    ///< Arena for the per-vertex and per-face tables.
  float pad; ///< Arc past its vertices that a moving face's scan can reach.
  ClipPieces clip;              ///< Render clip's column runs.
  SweptMeshCallbacks callbacks; ///< Vertex evaluator and mesh scanner.
};

/**
 * @brief Scans a moving mesh with each face's vertices at their column times.
 * @param deformation The scan's frame state, with the vertex pool at the
 *        current endpoint.
 * @details deform_mesh() describes the passes.
 */
HS_HOT_FLASH_MEMBER void scan_swept_faces(const SweptMesh &deformation);

} // namespace detail

/**
 * @brief Scans a moving mesh with every vertex at the time its column is lit.
 * @tparam W Display width.
 * @tparam H Display height.
 * @tparam VertexFn Callable (size_t i, float u) evaluating an endpoint-space vertex.
 * @tparam ScanFn Callable () scanning the mesh into the canvas's open window.
 * @param canvas Canvas with column-time state when deformation is active.
 * @param mesh Endpoint mesh, overwritten by the captured vertices.
 * @param active Whether the mesh deforms during this frame.
 * @param columns Per-face seeing columns; HELD_FACE_COLUMN marks held faces.
 * @param moving vertex_flags() for vertices whose paths are captured; 0 for
 *        those that hold.
 * @param start Each moving vertex's column at the frame's start, as
 *        vertex_column() gives it at u = 0 from the vertex evaluated there.
 * @param scratch Arena for face selection and vertex placement tables.
 * @param reach Bound on vertex motion in radians.
 * @param vertex Closed-form vertex path at a frame fraction, continuing the
 *        frame's motion past either endpoint.
 * @param scan Scans the mesh with its current vertices and selected faces.
 * @details Each moving vertex is captured at the sweep time of the column
 * that sees it (capture_vertex()), on its own arm or across the arm boundary
 * nearer it, once each, and shared by the faces around it. Held faces and
 * the moving faces clear of both arm boundaries draw in one scan. A face across an arm
 * boundary draws once per arm, clipped to that arm's columns, each edge's
 * vertex times continuing across the boundary the edge crosses, so frame
 * N's copy before a boundary and frame N+1's copy after it are one polygon
 * there. A face every column pair's meridian can reach (one around a pole,
 * spanning half a turn, or reaching the pole by its vertices' motion) draws
 * in one scan per half of the arms, each of its vertices at its time nearest
 * the half's middle; a vertex whose motion can carry it half a right angle
 * about the polar axis holds at its endpoint. A polar face whose single
 * placement would show a visible step draws one column pair at a time
 * instead, each pair's face standing as at the pair's sweep time
 * (Mesh::draw_face_sweep()).
 */
template <int W, int H, typename VertexFn, typename ScanFn>
HS_HOT_FLASH_MEMBER void
deform_mesh(Canvas &canvas, MeshState &mesh, bool active,
            const uint16_t *columns, const uint8_t *moving,
            const uint16_t *start, Arena &scratch, float reach,
            VertexFn &&vertex, ScanFn &&scan) {
  const detail::SweptMeshCallbacks callbacks(vertex, scan);
  if (!active) {
    callbacks.draw(callbacks.scan_context);
    return;
  }
  HS_CHECK(canvas.column_basis(),
           "mesh deformation needs the frame's column-time state installed");
  HS_CHECK(canvas.column_basis()->azimuth,
           "mesh deformation needs column azimuth state");
  ScratchScope guard(scratch);
  uint8_t *selected = scratch.allocate_n<uint8_t>(mesh.num_faces());
  Scan::ColumnBasis time = *canvas.column_basis();
  check_mesh_pass(time.mesh, mesh.num_faces(), mesh.vertices.size());
  time.mesh.selection_face_count = mesh.num_faces();
  time.mesh.selection_vertex_count = mesh.vertices.size();
  time.mesh.selected_faces = selected;
  const float sweep = time.sweep_arc;
  const auto pads = FacePadTables::prepare<W, H>(scratch, &time);
  time.mesh.face_pads = pads.full;
  time.mesh.compact_face_pads = pads.compact;
  time.mesh.wide_face_pads = pads.wide;
  ColumnBasisScope scope(canvas, time);
  detail::scan_swept_faces({canvas, mesh, columns, moving, start, selected,
                            scratch, scan_row_pad<W>(sweep) + reach,
                            clip_pieces<W>(canvas.render_x_clip()), callbacks});
}

} // namespace Scan
