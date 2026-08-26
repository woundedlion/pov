/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "mesh/mesh.h"

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
 * update_hankin to position a star point.
 */
struct HankinInstruction {
  uint16_t v_corner; /**< Index into base_vertices of the corner vertex. */
  uint16_t v_prev;   /**< Index into base_vertices of the previous corner. */
  uint16_t v_next;   /**< Index into base_vertices of the next corner. */
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
  ArenaVector<Vector> base_vertices; /**< Owned corner vertices (copy mode). */
  ArenaVector<Vector>
      static_vertices; /**< Edge midpoints; angle-independent. */
  ArenaVector<HankinInstruction>
      dynamic_instructions;         /**< One per dynamic vertex. */
  ArenaVector<uint8_t> face_counts; /**< Vertex count of each output face. */
  ArenaVector<uint16_t> faces;      /**< Flat vertex indices for all faces. */
  int static_offset =
      0; /**< static_vertices.size(); base for dynamic indices in faces. */
  /** Resolved corner source: base_vertices in copy mode, the borrowed input
   * vertices in borrow mode. HankinInstruction indices read through corner().
   * A span, not a raw pointer: in debug builds it carries the source's arena
   * stamps, so a corner read after the source arena was reset or rewound
   * faults instead of returning reissued bytes. */
  ArenaSpan<Vector> corner_src;
  /** MeshOps::connectivity_key of the faces this pattern emits; 0 until
   * compiled. Two seeds can share every census figure (a cube and an
   * octahedron both compile to 14 faces over 96 indices), so this is what
   * identifies the pattern a mesh's retained classification belongs to. */
  uint32_t topology_key = 0;

  /** Returns the corner vertex a HankinInstruction index refers to. */
  const Vector &corner(size_t i) const { return corner_src[i]; }

  /**
   * @brief Empties all owned vectors and drops the corner source.
   * @details ArenaVector::clear() zeroes the element count and keeps the arena
   * binding; the storage is reclaimed by the arena, not here.
   */
  void clear() {
    base_vertices.clear();
    static_vertices.clear();
    dynamic_instructions.clear();
    face_counts.clear();
    faces.clear();
    static_offset = 0;
    topology_key = 0;
    corner_src = ArenaSpan<Vector>();
  }

  /**
   * @brief Deep-copies all owned data into a target arena.
   * @param src Source instance to copy from; must own its corners, since a
   * borrow-mode source keeps no corner storage the clone could copy.
   * @param dst Destination instance whose vectors are bound and filled.
   * @param arena Arena backing the destination's freshly bound vectors.
   * @details Required by Cloneable; each vector is rebound from @p arena and
   * bulk-copied (all element types are trivially copyable). Traps if src
   * aliases dst: copy_vector rebinds dst in place, then memcpy's the block onto
   * itself.
   */
  HS_COLD_MEMBER static void clone(const CompiledHankin &src,
                                   CompiledHankin &dst, Arena &arena) {
    HS_CHECK(&src != &dst, "CompiledHankin::clone src must not alias dst");
    HS_CHECK(src.dynamic_instructions.size() == 0 ||
                 src.corner_src.data() == src.base_vertices.data(),
             "CompiledHankin::clone needs an owned-corner source "
             "(compile_hankin borrow mode has no corners to copy)");
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
    dst.corner_src = ArenaSpan<Vector>(dst.base_vertices);
  }
};

namespace MeshOps {

HS_O3_BEGIN

/**
 * @brief Bakes the angle-independent Hankin topology for a mesh.
 * @param mesh Input closed-manifold mesh to derive the pattern from.
 * @param compiled Output topology, allocated from @p target_arena.
 * @param target_arena Arena that backs the persistent compiled vectors.
 * @param temp_arena Arena holding the transient half-edge structures.
 * @param borrow_base_vertices When true, corner reads alias @p mesh's vertices
 * instead of copying them; the caller must keep @p mesh alive for every
 * update_hankin/hankin_at call that reads the compiled topology. Defaults to
 * false (owned copy), which every persisted or re-updated CompiledHankin needs.
 * @details Builds a half-edge mesh, emits one shared midpoint per edge into
 * static_vertices, reserves one dynamic (star-point) slot per half-edge, and
 * records the star and rosette faces. The closed-manifold precondition is
 * enforced up front, so every prev and pair read below is valid and never
 * HE_NONE.
 */
HS_COLD static void compile_hankin(const PolyMesh &mesh,
                                   CompiledHankin &compiled,
                                   Arena &target_arena, Arena &temp_arena,
                                   bool borrow_base_vertices = false) {
  // Topology via accessors (borrowed-mode safe); vertices is always owned.
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();

  // Unbind, not clear(): a reused CompiledHankin's bindings may name blocks its
  // arena has already reclaimed.
  compiled = CompiledHankin();

  if (borrow_base_vertices) {
    compiled.corner_src = ArenaSpan<Vector>(mesh.vertices);
  } else {
    compiled.base_vertices.bind(target_arena, V);
    for (size_t i = 0; i < V; ++i) {
      compiled.base_vertices.push_back(mesh.vertices[i]);
    }
    compiled.corner_src = ArenaSpan<Vector>(compiled.base_vertices);
  }
  // Static pool is I/2 midpoints; largest emitted index is (I/2)+(I-1). Guard
  // adds 1 to both sides to dodge the unsigned underflow of I-1 at I == 0.
  HS_CHECK((I / 2) + I <= static_cast<size_t>(INT16_MAX) + 1,
           "Hankin output vertex count exceeds int16_t index range "
           "(oversized source mesh?)");
  compiled.static_vertices.bind(target_arena, I / 2);
  compiled.dynamic_instructions.bind(target_arena, I);
  compiled.face_counts.bind(target_arena, F + V);
  compiled.faces.bind(target_arena, 4 * I);

  {
    ScratchScope temp_arena_guard(temp_arena);

    HalfEdgeMesh he_mesh(temp_arena, mesh);

    require_closed_manifold(he_mesh, temp_arena, "compile_hankin");

    uint16_t *he_to_midpoint_idx = temp_arena.allocate_n<uint16_t>(I);
    std::fill_n(he_to_midpoint_idx, I, HE_NONE);
    uint16_t *he_to_dynamic_idx = temp_arena.allocate_n<uint16_t>(I);
    std::fill_n(he_to_dynamic_idx, I, HE_NONE);

    // One shared midpoint per undirected edge, cached under both halves, so
    // every entry of he_to_midpoint_idx is populated from here on.
    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      const uint16_t he_idx = static_cast<uint16_t>(i);
      if (he_to_midpoint_idx[he_idx] != HE_NONE)
        continue;
      compiled.static_vertices.push_back(edge_midpoint(he_mesh, mesh, he_idx));
      const uint16_t idx = narrow_index(compiled.static_vertices.size() - 1);
      he_to_midpoint_idx[he_idx] = idx;
      he_to_midpoint_idx[he_mesh.half_edges[he_idx].pair] = idx;
    }

    compiled.static_offset = static_cast<int>(compiled.static_vertices.size());

    // Star faces
    for (size_t i = 0; i < he_mesh.faces.size(); ++i) {
      HEFace &face = he_mesh.faces[i];
      uint16_t he_idx = face.half_edge;
      HS_CHECK(he_idx != HE_NONE, "compile_hankin: empty face");
      uint16_t start_he = he_idx;
      int count = 0;
      int walked = 0;

      do {
        HS_CHECK(walked < (int)he_mesh.half_edges.size(),
                 "hankin star-face walk exceeded half-edge count");
        ++walked;
        count += 2;
        HalfEdge &curr_he = he_mesh.half_edges[he_idx];
        uint16_t prev_idx = curr_he.prev;

        int idx_m1 = he_to_midpoint_idx[prev_idx];
        int idx_m2 = he_to_midpoint_idx[he_idx];

        HalfEdge &prev_he = he_mesh.half_edges[prev_idx];

        uint16_t i_corner = prev_he.vertex;
        uint16_t i_prev = he_mesh.half_edges[prev_he.prev].vertex;
        uint16_t i_next = curr_he.vertex;

        compiled.dynamic_instructions.push_back({i_corner, i_prev, i_next,
                                                 narrow_index(idx_m1),
                                                 narrow_index(idx_m2)});

        // The instruction pushed above owns this star point, so its index is
        // the last instruction slot.
        uint16_t dyn_idx =
            narrow_index(compiled.dynamic_instructions.size() - 1);
        he_to_dynamic_idx[he_idx] = dyn_idx;

        compiled.faces.push_back(narrow_index(idx_m1));
        compiled.faces.push_back(
            narrow_index(compiled.static_offset + dyn_idx));

        he_idx = curr_he.next;
      } while (he_idx != HE_NONE && he_idx != start_he);

      compiled.face_counts.push_back(narrow_face_count(count));
    }

    // Rosette faces
    bool *visited_verts = temp_arena.allocate_n<bool>(V);
    std::fill_n(visited_verts, V, false);

    // Per-orbit scratch buffer. Each orbit step appends two indices (a
    // midpoint and a dynamic vertex), so the absolute upper bound on entries
    // is twice the total half-edge count.
    uint16_t *face_indices = temp_arena.allocate_n<uint16_t>(2 * I);

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
            HS_CHECK(count + 2 <= (int)(2 * I),
                     "Hankin rosette winding overflow");
            face_indices[count++] = he_to_midpoint_idx[curr_idx];
            const uint16_t next_edge_idx =
                he_mesh.half_edges[he_mesh.half_edges[curr_idx].pair].next;
            face_indices[count++] = narrow_index(
                compiled.static_offset + he_to_dynamic_idx[next_edge_idx]);
          });

      // count = 2 * vertex degree. Degree-2 is legal (hankin-of-hankin walks
      // its own degree-2 midpoints -> quad rosette); only degree < 2
      // degenerates.
      HS_CHECK(count >= 4, "Hankin rosette winding has degree < 2");
      compiled.face_counts.push_back(narrow_face_count(count));
      for (int k = count - 1; k >= 0; --k) {
        compiled.faces.push_back(face_indices[k]);
      }
    }
  }

  compiled.topology_key = MeshOps::connectivity_key(
      compiled.face_counts.data(), compiled.face_counts.size(),
      compiled.faces.data(), compiled.faces.size());
}

/** Squared endpoints of the far-intersection blend: the edge-midpoint fallback
 * ramps in at 2.25 and fully replaces the intersection at 4.0. Measured healthy
 * registry intersections peak at 2.16. */
inline constexpr float STAR_FAR_BLEND_START_RATIO_SQ = 2.25f;
inline constexpr float STAR_FAR_RATIO_SQ = 4.0f;
inline constexpr float HANKIN_PARALLEL_REGULARIZATION_SQ = 3.0e-4f;
inline constexpr float HANKIN_CONDITIONED_CLEAR_SQ = 2.0e-2f;
inline constexpr float HANKIN_CONDITIONED_NEAR_RATIO_SQ = 1.44f;
inline constexpr float HANKIN_CONDITIONED_FAR_RATIO_SQ = 9.0f;
/** plane_cross_sq (= |cross(n_hankin1, n_hankin2)|^2) window gating the
 * far-star fallback. A far chord ratio alone is not instability: a large face
 * (hexagon, octagon) legitimately pushes its star point out to raw_ratio_sq ~3
 * while the two contact planes stay well-separated (plane_cross_sq ~0.8).
 * Genuine near-parallel yeets — where the intersection direction is
 * noise-dominated — sit at plane_cross_sq below 0.01; healthy far points stay
 * above 0.5. plane_cross_sq varies smoothly with the sweep angle, so gating on
 * it keeps the transition continuous. */
inline constexpr float HANKIN_PARALLEL_GATE_LO_SQ = 0.05f;
inline constexpr float HANKIN_PARALLEL_GATE_HI_SQ = 0.30f;

/**
 * @brief Positions the angle-dependent vertices and writes the final mesh.
 * @tparam MeshT Output mesh type; may optionally provide face_offsets.
 * @param compiled Baked angle-independent topology.
 * @param out_mesh Output mesh, allocated from @p target_arena. Its topology
 *   array is retained, not rebuilt, so one classification serves every angle
 *   re-solve of the same compiled pattern; a mesh reused for a DIFFERENT
 *   pattern must be cleared first, which the topology-key check below
 *   enforces.
 * @param target_arena Arena backing @p out_mesh's vertex and face vectors.
 * @param angle Contact angle in radians; domain [0, pi/2]. At ~0 the star points
 *   collapse onto their corners (flat tiling); larger angles push the rays out so
 *   the rays from adjacent edges intersect to form sharper star points. Outside
 *   the domain the construction aliases onto an in-domain pattern rather than
 *   failing — it is 2*pi-periodic in the angle and mirrors under negation — so a
 *   caller at an untrusted boundary must range-check it.
 * @details A corner whose contact planes are near-parallel has no nearby
 *   intersection; its star point falls back to the edge-midpoint mean instead
 *   of being flung across the sphere (see STAR_FAR_RATIO_SQ).
 */
template <typename MeshT>
HS_COLD_MEMBER inline void update_hankin(const CompiledHankin &compiled,
                                         MeshT &out_mesh, Arena &target_arena,
                                         float angle) {

  // A borrowed MeshState carries its topology in a view that set_owned() drops,
  // so the reuse check below reads the size through the mode-aware accessor
  // first.
  size_t prior_topology_size = out_mesh.topology.size();
  if constexpr (requires { out_mesh.get_topology_size(); }) {
    prior_topology_size = out_mesh.get_topology_size();
  }

  // Drop any borrowed-mode views a reused MeshState may still carry, so the
  // owned vertex and face arrays bound below are not shadowed by a stale view.
  if constexpr (requires { out_mesh.set_owned(); }) {
    out_mesh.set_owned();
  }

  HS_CHECK(prior_topology_size == 0 ||
               out_mesh.topology_key == compiled.topology_key,
           "update_hankin: reused out_mesh carries a topology from a different "
           "compiled pattern (clear it first)");

  // A borrowed classification is gone with the view dropped above; the key
  // describing it goes with it, so the output leaves here unclassified rather
  // than keyed to a topology it no longer carries.
  if (out_mesh.topology.size() == 0)
    out_mesh.topology_key = 0;

  HS_CHECK(std::isfinite(angle), "update_hankin: contact angle must be finite");

  bool is_flat = std::abs(angle) < math::TOLERANCE;

  // Star points are computed straight into the output: nothing reads them back
  // across calls, so the compiled topology holds no vertex scratch of its own.
  out_mesh.vertices.bind(target_arena,
                         compiled.static_vertices.size() +
                             compiled.dynamic_instructions.size());
  out_mesh.vertices.append_bulk(compiled.static_vertices.data(),
                                compiled.static_vertices.size());

  float cos_ha = cosf(angle * 0.5f);
  float sin_ha = sinf(angle * 0.5f);
  const auto blend_above = [](float value, float start, float end) {
    const float t =
        std::max(0.0f, std::min(1.0f, (value - start) / (end - start)));
    return quintic_kernel(t);
  };
  HS_CHECK(compiled.corner_src.data() != nullptr,
           "update_hankin needs a compiled topology (corner_src is null after "
           "CompiledHankin::clear)");
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &instr = compiled.dynamic_instructions[i];
    Vector p_corner = compiled.corner(instr.v_corner);

    if (is_flat) {
      out_mesh.vertices.push_back(normalized_or(p_corner, p_corner));
      continue;
    }

    Vector m1 = compiled.static_vertices[instr.idx_m1];
    Vector m2 = compiled.static_vertices[instr.idx_m2];
    Vector p_prev = compiled.corner(instr.v_prev);
    Vector p_next = compiled.corner(instr.v_next);

    Vector cross1 = cross(p_prev, p_corner);
    Vector cross2 = cross(p_corner, p_next);

    if (dot(cross1, cross1) < math::EPS_CROSS_SQ ||
        dot(cross2, cross2) < math::EPS_CROSS_SQ) {
      out_mesh.vertices.push_back(normalized_or(p_corner, p_corner));
      continue; // zero length edge
    }

    Vector n_edge1 = cross1.normalized();
    // Sign convention: the two edge normals are rotated by OPPOSITE-signed half
    // contact angles (+ha about m1, -ha about m2) so both Hankin planes tilt
    // toward the shared corner; the dot(intersect, p_corner)<0 flip below then
    // selects the corner-side hemisphere of their intersection.
    // m1/m2 are unit (midpoints normalized at compile time), so (cos_ha,
    // sin_ha*axis) is already a unit quaternion as rotate() requires.
    Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
    Vector n_hankin1 = rotate(n_edge1, q1);

    Vector n_edge2 = cross2.normalized();
    // cos(-x) = cos(x), sin(-x) = -sin(x); m2 unit per the precondition above.
    Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
    Vector n_hankin2 = rotate(n_edge2, q2);

    Vector intersect = cross(n_hankin1, n_hankin2);
    const float plane_cross_sq = dot(intersect, intersect);
    Vector cn = normalized_or(p_corner, p_corner);
    Vector fallback = normalized_or(m1 + m2, cn);
    if (dot(fallback, p_corner) < 0.0f)
      fallback = -fallback;
    const Vector oriented_intersect = intersect * dot(intersect, p_corner);
    const Vector raw_star = normalized_or(oriented_intersect, fallback);
    const float local_sq =
        std::max(distance_squared(m1, cn), distance_squared(m2, cn));
    if (!(local_sq > math::EPS_LEN_SQ)) {
      out_mesh.vertices.push_back(fallback);
      continue;
    }
    const float raw_ratio_sq = distance_squared(raw_star, cn) / local_sq;
    const float conditioned =
        blend_above(raw_ratio_sq, HANKIN_CONDITIONED_NEAR_RATIO_SQ,
                    HANKIN_CONDITIONED_FAR_RATIO_SQ);
    const float anchor =
        std::max(0.0f, HANKIN_PARALLEL_REGULARIZATION_SQ - plane_cross_sq) +
        conditioned *
            std::max(0.0f, HANKIN_CONDITIONED_CLEAR_SQ - plane_cross_sq);
    intersect = normalized_or(oriented_intersect + fallback * anchor, fallback);
    // Near-parallel contact planes fling the intersection geodesically far
    // from the corner, yielding a sliver face that renders as a long line.
    const float ratio_sq = distance_squared(intersect, cn) / local_sq;
    // Gate the fallback on parallelism: blend_above(-plane_cross_sq, ...) rises
    // as plane_cross_sq falls — 0 above the HI threshold (planes well-separated,
    // keep the true intersection), 1 below LO (near-parallel, take the
    // fallback). A far chord ratio then only triggers the fallback in the
    // genuinely unstable regime.
    const float parallel_gate =
        blend_above(-plane_cross_sq, -HANKIN_PARALLEL_GATE_HI_SQ,
                    -HANKIN_PARALLEL_GATE_LO_SQ);
    const float fallback_blend =
        blend_above(ratio_sq, STAR_FAR_BLEND_START_RATIO_SQ,
                    STAR_FAR_RATIO_SQ) *
        parallel_gate;
    if (fallback_blend > 0.0f) {
      out_mesh.vertices.push_back(
          fallback_blend >= 1.0f
              ? fallback
              : normalized_or(intersect * (1.0f - fallback_blend) +
                                  fallback * fallback_blend,
                              fallback));
      continue;
    }

    out_mesh.vertices.push_back(intersect);
  }

  out_mesh.face_counts.bind(target_arena, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.bind(target_arena, compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      HS_CHECK(current_offset + compiled.face_counts[i] <= UINT16_MAX,
               "mesh face_offsets exceeds 16-bit index range");
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.bind(target_arena, compiled.faces.size());
  out_mesh.faces.append_bulk(compiled.faces.data(), compiled.faces.size());
}

/**
 * @brief One-shot Hankin pattern: compile then update, returning the new mesh.
 * @param mesh Input closed-manifold mesh to derive the pattern from.
 * @param target Arena that backs the returned mesh's persistent data. In this
 *   one-shot path it also serves as compile_hankin's working arena (see the
 *   reversed polarity below), so it must hold max(compile scratch, output) —
 *   sizing it for the output mesh alone under-provisions.
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
    // Arena polarity is reversed from the streaming path: the throwaway
    // CompiledHankin is allocated from `temp` while `target` serves as
    // compile_hankin's working arena, then update_hankin builds `out` into it.
    compile_hankin(mesh, compiled, temp, target,
                   /*borrow_base_vertices=*/true);
    update_hankin(compiled, out, target, angle);
  }

  return out;
}

HS_O3_END

} // namespace MeshOps
