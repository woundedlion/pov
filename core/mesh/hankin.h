/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/mesh/mesh.h"

/**
 * @file hankin.h
 * @brief Hankin (polygons-in-contact) star-and-rosette pattern generation on a
 *        spherical mesh.
 *
 * Each edge of the input mesh contributes a midpoint, and from every corner two
 * "contact" rays leave the edge midpoints at a fixed contact angle; where
 * neighbouring rays meet they form the star points whose position varies with
 * the angle. compile_hankin() bakes the angle-independent topology once;
 * update_hankin() then recomputes only the angle-dependent vertices.
 */

/**
 * @brief One angle-dependent (dynamic) vertex's input topology, consumed by
 * MeshOps::solve_hankin_vertices to position a star point.
 */
struct HankinInstruction {
  uint16_t
      v_corner; /**< Index into corner_src (read via corner()) of the corner vertex. */
  uint16_t
      v_prev; /**< Index into corner_src (read via corner()) of the previous corner. */
  uint16_t
      v_next; /**< Index into corner_src (read via corner()) of the next corner. */
  uint16_t
      idx_m1; /**< Index into static_vertices of the first edge midpoint. */
  uint16_t
      idx_m2; /**< Index into static_vertices of the second edge midpoint. */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 *
 * The output mesh concatenates static vertices then dynamic vertices, so an
 * index < static_offset refers to static_vertices and an index >=
 * static_offset refers to star point (index - static_offset).
 */
struct CompiledHankin {
  ArenaVector<math::Vector>
      base_vertices; /**< Owned corner vertices (copy mode). */
  ArenaVector<math::Vector>
      static_vertices; /**< Edge midpoints; angle-independent. */
  ArenaVector<HankinInstruction>
      dynamic_instructions;         /**< One per dynamic vertex. */
  ArenaVector<uint8_t> face_counts; /**< Vertex count of each output face. */
  ArenaVector<uint16_t> faces;      /**< Flat vertex indices for all faces. */
  int static_offset =
      0; /**< static_vertices.size(); base for dynamic indices in faces. */
  /** Resolved corner source: base_vertices in copy mode, the borrowed input
   * vertices in borrow mode. HankinInstruction indices read through corner(). */
  ArenaSpan<math::Vector> corner_src;
  /** MeshOps::connectivity_key of the faces this pattern emits; 0 until
   * compiled. */
  uint32_t topology_key = 0;

  /**
   * @brief Returns the corner vertex a HankinInstruction index refers to.
   * @param i Corner index into `corner_src`.
   * @return The corner vertex.
   */
  const math::Vector &corner(size_t i) const { return corner_src[i]; }

  /**
   * @brief Empties all owned vectors and drops the corner source.
   * @details The storage is reclaimed by the arena, not here.
   */
  void clear() {
    base_vertices.clear();
    static_vertices.clear();
    dynamic_instructions.clear();
    face_counts.clear();
    faces.clear();
    static_offset = 0;
    topology_key = 0;
    corner_src = ArenaSpan<math::Vector>();
  }

  /**
   * @brief Deep-copies all owned data into a target arena.
   * @param src Source instance to copy from; must own its corners, since a
   * borrow-mode source keeps no corner storage the clone could copy.
   * @param dst Destination instance whose vectors are bound and filled.
   * @param arena Arena backing the destination's freshly bound vectors.
   * @details Required by Cloneable; each vector is rebound from @p arena and
   * bulk-copied. Traps if src aliases dst.
   */
  HS_FLASH_MEMBER static void clone(const CompiledHankin &src,
                                    CompiledHankin &dst, Arena &arena) {
    HS_CHECK(&src != &dst, "CompiledHankin::clone src must not alias dst");
    HS_CHECK(src.dynamic_instructions.size() == 0 ||
                 src.corner_src.data() == src.base_vertices.data(),
             "CompiledHankin::clone needs an owned-corner source "
             "(compile_hankin borrow mode has no corners to copy)");
    dst = CompiledHankin();
    MeshOps::copy_vector(dst.base_vertices, src.base_vertices.data(),
                         src.base_vertices.size(), arena);
    MeshOps::copy_vector(dst.static_vertices, src.static_vertices.data(),
                         src.static_vertices.size(), arena);
    MeshOps::copy_vector(dst.dynamic_instructions,
                         src.dynamic_instructions.data(),
                         src.dynamic_instructions.size(), arena);
    MeshOps::copy_vector(dst.face_counts, src.face_counts.data(),
                         src.face_counts.size(), arena);
    MeshOps::copy_vector(dst.faces, src.faces.data(), src.faces.size(), arena);
    dst.static_offset = src.static_offset;
    dst.topology_key = src.topology_key;
    dst.corner_src = ArenaSpan<math::Vector>(dst.base_vertices);
  }
};

namespace MeshOps {

HS_O3_BEGIN

/**
 * @brief How a CompiledHankin sources its corner vertices.
 */
enum class BaseVertices : uint8_t {
  OWNED_COPY, /**< Copy the mesh vertices into the target arena. */
  BORROW,     /**< Alias the source mesh's vertices. */
};

/** @brief Emits star faces and their dynamic vertex instructions.
 * @param he_mesh Closed source topology.
 * @param compiled Receives faces and instructions.
 * @param he_to_midpoint_idx Midpoint index per half edge.
 * @param he_to_dynamic_idx Receives dynamic indices. */
HS_COLD static void emit_hankin_stars(const HalfEdgeMesh &he_mesh,
                                      CompiledHankin &compiled,
                                      const uint16_t *he_to_midpoint_idx,
                                      uint16_t *he_to_dynamic_idx) {
  for (size_t i = 0; i < he_mesh.faces.size(); ++i) {
    const HEFace &face = he_mesh.faces[i];
    int count = 0;
    face_loop(he_mesh, face.half_edge,
              static_cast<int>(he_mesh.half_edges.size()),
              [&](uint16_t he_idx) {
                count += 2;
                const HalfEdge &curr_he = he_mesh.half_edges[he_idx];
                uint16_t prev_idx = curr_he.prev;

                int idx_m1 = he_to_midpoint_idx[prev_idx];
                int idx_m2 = he_to_midpoint_idx[he_idx];

                const HalfEdge &prev_he = he_mesh.half_edges[prev_idx];

                uint16_t i_corner = prev_he.vertex;
                uint16_t i_prev = he_mesh.half_edges[prev_he.prev].vertex;
                uint16_t i_next = curr_he.vertex;

                compiled.dynamic_instructions.push_back(
                    {i_corner, i_prev, i_next, narrow_index(idx_m1),
                     narrow_index(idx_m2)});

                uint16_t dyn_idx =
                    narrow_index(compiled.dynamic_instructions.size() - 1);
                he_to_dynamic_idx[he_idx] = dyn_idx;

                compiled.faces.push_back(narrow_index(idx_m1));
                compiled.faces.push_back(
                    narrow_index(compiled.static_offset + dyn_idx));
              });

    compiled.face_counts.push_back(narrow_face_count(count));
  }
}

/** @brief Emits a rosette around each source vertex.
 * @param he_mesh Closed source topology.
 * @param compiled Receives rosette faces.
 * @param temp_arena Orbit scratch storage.
 * @param v Source vertex count.
 * @param i_value Source half-edge count.
 * @param he_to_midpoint_idx Midpoint index per half edge.
 * @param he_to_dynamic_idx Dynamic index per half edge. */
HS_COLD static void emit_hankin_rosettes(const HalfEdgeMesh &he_mesh,
                                         CompiledHankin &compiled,
                                         Arena &temp_arena, size_t v,
                                         size_t i_value,
                                         const uint16_t *he_to_midpoint_idx,
                                         const uint16_t *he_to_dynamic_idx) {
  bool *visited_verts = temp_arena.allocate_n<bool>(v);
  std::fill_n(visited_verts, v, false);

  // Per-orbit scratch: two indices (midpoint, dynamic vertex) per orbit step.
  uint16_t *face_indices = temp_arena.allocate_n<uint16_t>(2 * i_value);

  for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
    uint16_t he_start_idx = static_cast<uint16_t>(i);
    const HalfEdge &he_start = he_mesh.half_edges[he_start_idx];
    // Rosette origin: he_start's tail vertex.
    uint16_t origin_idx = he_mesh.half_edges[he_start.prev].vertex;
    if (visited_verts[origin_idx])
      continue;
    visited_verts[origin_idx] = true;

    int count = 0;

    vertex_orbit<OrbitDir::PAIR_NEXT>(
        he_mesh, he_start_idx, [&](uint16_t curr_idx) {
          HS_CHECK(count + 2 <= (int)(2 * i_value),
                   "Hankin rosette winding overflow");
          face_indices[count++] = he_to_midpoint_idx[curr_idx];
          const uint16_t next_edge_idx =
              he_mesh.half_edges[he_mesh.half_edges[curr_idx].pair].next;
          face_indices[count++] = narrow_index(
              compiled.static_offset + he_to_dynamic_idx[next_edge_idx]);
        });

    // count = 2 * vertex degree; degree 2 (quad rosette) is legal.
    HS_CHECK(count >= 4, "Hankin rosette winding has degree < 2");
    compiled.face_counts.push_back(narrow_face_count(count));
    for (int k = count - 1; k >= 0; --k) {
      compiled.faces.push_back(face_indices[k]);
    }
  }
}

/**
 * @brief Bakes the angle-independent Hankin topology for a mesh.
 * @param mesh Input closed-manifold mesh to derive the pattern from.
 * @param compiled Output topology, allocated from @p target_arena.
 * @param target_arena Arena that backs the persistent compiled vectors.
 * @param temp_arena Arena holding the transient half-edge structures.
 * @param base How corner reads source vertices. BORROW requires @p mesh to
 * remain alive whenever the compiled topology is read.
 * @details Builds a half-edge mesh, emits one shared midpoint per edge into
 * static_vertices, reserves one dynamic (star-point) slot per half-edge, and
 * records the star and rosette faces.
 */
HS_COLD static void
compile_hankin(const PolyMesh &mesh, CompiledHankin &compiled,
               Arena &target_arena, Arena &temp_arena,
               BaseVertices base = BaseVertices::OWNED_COPY) {
  size_t v = mesh.vertices.size();
  size_t f = mesh.get_face_counts_size();
  size_t i_value = mesh.get_faces_size();

  // Unbind, not clear(): a reused CompiledHankin's bindings may name blocks its
  // arena has already reclaimed.
  compiled = CompiledHankin();

  if (base == BaseVertices::BORROW) {
    compiled.corner_src = ArenaSpan<math::Vector>(mesh.vertices);
  } else {
    compiled.base_vertices.bind(target_arena, v);
    for (size_t i = 0; i < v; ++i) {
      compiled.base_vertices.push_back(mesh.vertices[i]);
    }
    compiled.corner_src = ArenaSpan<math::Vector>(compiled.base_vertices);
  }
  // Static pool is i_value/2 midpoints; largest emitted index is (i_value/2)+(i_value-1). Guard
  // adds 1 to both sides to dodge the unsigned underflow of i_value-1 at i_value == 0.
  HS_CHECK((i_value / 2) + i_value <= MeshLimits::MAX_VERTICES,
           "Hankin output vertex count exceeds int16_t index range "
           "(oversized source mesh?)");
  compiled.static_vertices.bind(target_arena, i_value / 2);
  compiled.dynamic_instructions.bind(target_arena, i_value);
  compiled.face_counts.bind(target_arena, f + v);
  compiled.faces.bind(target_arena, 4 * i_value);

  {
    ScratchScope temp_arena_guard(temp_arena);

    HalfEdgeMesh he_mesh(temp_arena, mesh);

    require_closed_manifold(he_mesh, temp_arena, "compile_hankin");

    uint16_t *he_to_midpoint_idx = temp_arena.allocate_n<uint16_t>(i_value);
    std::fill_n(he_to_midpoint_idx, i_value, HE_NONE);
    uint16_t *he_to_dynamic_idx = temp_arena.allocate_n<uint16_t>(i_value);
    std::fill_n(he_to_dynamic_idx, i_value, HE_NONE);

    emit_edge_midpoints(he_mesh, mesh, compiled.static_vertices,
                        he_to_midpoint_idx, [](const HalfEdge &) {});

    compiled.static_offset = static_cast<int>(compiled.static_vertices.size());

    emit_hankin_stars(he_mesh, compiled, he_to_midpoint_idx, he_to_dynamic_idx);

    emit_hankin_rosettes(he_mesh, compiled, temp_arena, v, i_value,
                         he_to_midpoint_idx, he_to_dynamic_idx);
  }

  compiled.topology_key = MeshOps::connectivity_key(
      compiled.face_counts.data(), compiled.face_counts.size(),
      compiled.faces.data(), compiled.faces.size());
}

/** Squared ratio endpoints of the far-intersection blend. The ratio ramp is
 * scaled by the parallel gate, so the edge-midpoint fallback fully replaces
 * the intersection from the far ratio onward only on near-parallel planes. */
inline constexpr float STAR_FAR_BLEND_START_RATIO_SQ = 2.25f;
/// Squared ratio at the far end of the blend.
inline constexpr float STAR_FAR_RATIO_SQ = 4.0f;
/** Plane-cross squared floor below which fallback is always mixed in. */
inline constexpr float HANKIN_PARALLEL_REGULARIZATION_SQ = 3.0e-4f;
/** Plane-cross squared threshold below which far points are anchored. */
inline constexpr float HANKIN_CONDITIONED_CLEAR_SQ = 2.0e-2f;
/** Raw chord ratio squared at the start of the conditioning ramp. */
inline constexpr float HANKIN_CONDITIONED_NEAR_RATIO_SQ = 1.44f;
/** Raw chord ratio squared at full conditioning. */
inline constexpr float HANKIN_CONDITIONED_FAR_RATIO_SQ = 9.0f;
/** plane_cross_sq (= |cross(n_hankin1, n_hankin2)|^2) window gating the
 * far-star fallback: fully on below LO (near-parallel planes), off above HI. */
inline constexpr float HANKIN_PARALLEL_GATE_LO_SQ = 0.05f;
/// Gate window upper end; the fallback is off above it.
inline constexpr float HANKIN_PARALLEL_GATE_HI_SQ = 0.30f;

/** Branch used to position one dynamic Hankin vertex. */
enum class HankinBranch : uint8_t { INTERSECT, COLLAPSED, FALLBACK, BLENDED };

/** Position and conditioning diagnostics for one dynamic Hankin vertex. */
struct HankinSolve {
  math::Vector pos;    ///< Solved vertex position.
  HankinBranch branch; ///< Conditioning path used for the position.
  /** Squared star-to-corner chord divided by the largest midpoint-to-corner
   * squared chord; zero on degenerate and collapsed branches. */
  float far_ratio = 0;
};

/** @brief Solves one dynamic vertex using the prepared contact half-angle.
 * @param compiled Static vertices and corner lookup data.
 * @param instr Midpoint and corner indices for this vertex.
 * @param is_flat Collapse the vertex to its corner direction.
 * @param cos_ha Cosine of half the contact angle.
 * @param sin_ha Sine of half the contact angle.
 * @return Position and conditioning diagnostics. */
HS_HOT_INLINE inline HankinSolve
solve_hankin_vertex(const CompiledHankin &compiled,
                    const HankinInstruction &instr, bool is_flat, float cos_ha,
                    float sin_ha) {
  const auto blend_above = [](float value, float start, float end) {
    const float t = hs::clamp((value - start) / (end - start), 0.0f, 1.0f);
    return math::quintic_kernel(t);
  };

  math::Vector p_corner = compiled.corner(instr.v_corner);

  if (is_flat) {
    return {math::normalized_or(p_corner, p_corner), HankinBranch::COLLAPSED};
  }

  math::Vector m1 = compiled.static_vertices[instr.idx_m1];
  math::Vector m2 = compiled.static_vertices[instr.idx_m2];
  math::Vector p_prev = compiled.corner(instr.v_prev);
  math::Vector p_next = compiled.corner(instr.v_next);

  math::Vector cross1 = math::cross(p_prev, p_corner);
  math::Vector cross2 = math::cross(p_corner, p_next);

  if (math::dot(cross1, cross1) < math::EPS_CROSS_SQ ||
      math::dot(cross2, cross2) < math::EPS_CROSS_SQ) {
    return {math::normalized_or(p_corner, p_corner), HankinBranch::COLLAPSED};
  }

  math::Vector n_edge1 = cross1.normalized();
  // Opposite-signed contact angles (+ha about m1, -ha about m2) tilt both
  // Hankin planes toward the shared corner; m1/m2 are unit, so q1/q2 are too.
  math::Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
  math::Vector n_hankin1 = math::rotate(n_edge1, q1);

  math::Vector n_edge2 = cross2.normalized();
  math::Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
  math::Vector n_hankin2 = math::rotate(n_edge2, q2);

  math::Vector intersect = math::cross(n_hankin1, n_hankin2);
  const float plane_cross_sq = math::dot(intersect, intersect);
  math::Vector cn = math::normalized_or(p_corner, p_corner);
  math::Vector fallback = math::normalized_or(m1 + m2, cn);
  if (math::dot(fallback, p_corner) < 0.0f)
    fallback = -fallback;
  const math::Vector oriented_intersect =
      intersect * math::dot(intersect, p_corner);
  const math::Vector raw_star =
      math::normalized_or(oriented_intersect, fallback);
  const float local_sq =
      fmaxf(math::distance_squared(m1, cn), math::distance_squared(m2, cn));
  if (!(local_sq > math::EPS_LEN_SQ)) {
    return {fallback, HankinBranch::FALLBACK};
  }
  const float raw_ratio_sq = math::distance_squared(raw_star, cn) / local_sq;
  const float conditioned =
      blend_above(raw_ratio_sq, HANKIN_CONDITIONED_NEAR_RATIO_SQ,
                  HANKIN_CONDITIONED_FAR_RATIO_SQ);
  const float anchor =
      fmaxf(0.0f, HANKIN_PARALLEL_REGULARIZATION_SQ - plane_cross_sq) +
      conditioned * fmaxf(0.0f, HANKIN_CONDITIONED_CLEAR_SQ - plane_cross_sq);
  intersect =
      math::normalized_or(oriented_intersect + fallback * anchor, fallback);
  // Near-parallel contact planes fling the intersection geodesically far
  // from the corner, yielding a sliver face that renders as a long line.
  const float ratio_sq = math::distance_squared(intersect, cn) / local_sq;
  // parallel_gate: 0 above HI (planes well-separated), 1 below LO
  // (near-parallel).
  const float parallel_gate =
      blend_above(-plane_cross_sq, -HANKIN_PARALLEL_GATE_HI_SQ,
                  -HANKIN_PARALLEL_GATE_LO_SQ);
  const float fallback_blend =
      blend_above(ratio_sq, STAR_FAR_BLEND_START_RATIO_SQ, STAR_FAR_RATIO_SQ) *
      parallel_gate;
  if (fallback_blend > 0.0f) {
    return {fallback_blend >= 1.0f
                ? fallback
                : math::normalized_or(intersect * (1.0f - fallback_blend) +
                                          fallback * fallback_blend,
                                      fallback),
            fallback_blend >= 1.0f ? HankinBranch::FALLBACK
                                   : HankinBranch::BLENDED,
            ratio_sq};
  }

  return {intersect, HankinBranch::INTERSECT, ratio_sq};
}

/**
 * @brief Positions static and angle-dependent Hankin vertices.
 * @param compiled Baked topology with corner data.
 * @param vertices Output positions, allocated from target_arena.
 * @param target_arena Arena backing the vertex vector.
 * @param angle Finite contact angle in radians.
 */
HS_COLD static void solve_hankin_vertices(const CompiledHankin &compiled,
                                          ArenaVector<math::Vector> &vertices,
                                          Arena &target_arena, float angle) {
  bool is_flat = std::abs(angle) < math::TOLERANCE;

  vertices.bind(target_arena, compiled.static_vertices.size() +
                                  compiled.dynamic_instructions.size());
  vertices.append_bulk(compiled.static_vertices.data(),
                       compiled.static_vertices.size());

  float cos_ha = cosf(angle * 0.5f);
  float sin_ha = sinf(angle * 0.5f);
  HS_CHECK(compiled.corner_src.data() != nullptr,
           "update_hankin needs a compiled topology (corner_src is null after "
           "CompiledHankin::clear)");
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &instr = compiled.dynamic_instructions[i];
    vertices.push_back(
        solve_hankin_vertex(compiled, instr, is_flat, cos_ha, sin_ha).pos);
  }
}

/**
 * @brief Positions the angle-dependent vertices and writes the final mesh.
 * @tparam MeshT MeshLike output mesh type; may optionally provide
 *   face_offsets.
 * @param compiled Baked angle-independent topology.
 * @param out_mesh Output mesh, allocated from @p target_arena. Its topology
 *   array is retained, not rebuilt; a mesh reused for a different pattern must
 *   be cleared first or this traps.
 * @param target_arena Arena backing @p out_mesh's vertex and face vectors.
 * @param angle Contact angle in radians; domain [0, pi/2]. At ~0 the star points
 *   collapse onto their corners (flat tiling); larger angles form sharper star
 *   points. The construction is pi-periodic. Angles in (pi/2, pi) modulo pi
 *   tilt the contact planes away from the corner and place star points outside
 *   their face.
 * @details A corner whose contact planes are near-parallel has no nearby
 *   intersection; its star point falls back to the edge-midpoint mean instead
 *   of being flung across the sphere (see STAR_FAR_RATIO_SQ).
 */
template <MeshLike MeshT>
HS_FLASH_MEMBER inline void update_hankin(const CompiledHankin &compiled,
                                          MeshT &out_mesh, Arena &target_arena,
                                          float angle) {

  // Read the topology size before set_owned() drops a borrowed view.
  size_t prior_topology_size = out_mesh.topology.size();
  if constexpr (requires { out_mesh.get_topology_size(); }) {
    prior_topology_size = out_mesh.get_topology_size();
  }

  // topology remains unbound; discard any borrowed classification.
  if constexpr (requires { out_mesh.set_owned(); }) {
    out_mesh.set_owned();
  }

  HS_CHECK(prior_topology_size == 0 ||
               out_mesh.topology_key == compiled.topology_key,
           "update_hankin: reused out_mesh carries a topology from a different "
           "compiled pattern (clear it first)");

  // A dropped borrowed classification takes its key with it.
  if (out_mesh.topology.size() == 0)
    out_mesh.topology_key = 0;

  HS_CHECK(std::isfinite(angle), "update_hankin: contact angle must be finite");

  solve_hankin_vertices(compiled, out_mesh.vertices, target_arena, angle);

  out_mesh.face_counts.bind(target_arena, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.bind(target_arena, compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      MeshOps::push_face_offset(out_mesh.face_offsets, current_offset,
                                compiled.face_counts[i]);
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.bind(target_arena, compiled.faces.size());
  out_mesh.faces.append_bulk(compiled.faces.data(), compiled.faces.size());
}

/**
 * @brief One-shot Hankin pattern: compile then update, returning the new mesh.
 * @param mesh Input closed-manifold mesh to derive the pattern from.
 * @param target Arena that backs the returned mesh's persistent data. It also
 *   serves as compile_hankin's working arena, so it must hold max(compile
 *   scratch, output).
 * @param temp Arena for the transient compiled topology, discarded on return.
 * @param angle Contact angle in radians; domain [0, pi/2] (see update_hankin).
 * @return The generated Hankin PolyMesh, allocated from @p target.
 * @details Convenience wrapper for callers that do not vary the angle. To
 * animate the angle, keep a CompiledHankin and call update_hankin repeatedly.
 */
HS_COLD static PolyMesh hankin(const PolyMesh &mesh, Arena &target, Arena &temp,
                               float angle) {
  HS_CHECK(&target != &temp, "hankin: target and temp must differ");
  PolyMesh out;

  {
    ScratchScope temp_guard(temp);
    CompiledHankin compiled;
    // The throwaway CompiledHankin lives in `temp`; `target` is
    // compile_hankin's working arena, then update_hankin builds `out` into it.
    compile_hankin(mesh, compiled, temp, target, BaseVertices::BORROW);
    update_hankin(compiled, out, target, angle);
  }

  return out;
}

HS_O3_END

} // namespace MeshOps
