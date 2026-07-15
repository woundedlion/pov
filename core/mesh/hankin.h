/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh/mesh.h"

/*
 * Hankin (polygons-in-contact) star-and-rosette pattern generation on a
 * spherical mesh. Each edge of the input mesh contributes a midpoint, and from
 * every corner two "contact" rays leave the edge midpoints at a fixed contact
 * angle; where neighbouring rays meet they form the star points whose position
 * varies with the angle. compile_hankin() bakes the angle-independent topology
 * once; update_hankin() then recomputes only the angle-dependent vertices.
 */

/**
 * @brief One angle-dependent (dynamic) vertex's input topology, consumed by
 * update_hankin to position a star point.
 */
struct HankinInstruction {
  uint16_t v_corner; /**< Index into base_vertices of the corner vertex. */
  uint16_t v_prev;   /**< Index into base_vertices of the previous corner. */
  uint16_t v_next;   /**< Index into base_vertices of the next corner. */
  uint16_t idx_m1;   /**< Index into static_vertices of the first edge midpoint. */
  uint16_t idx_m2;   /**< Index into static_vertices of the second edge midpoint. */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 *
 * The output mesh concatenates static vertices then dynamic vertices, so an
 * index < static_offset refers to static_vertices and an index >=
 * static_offset refers to dynamic_vertices[index - static_offset].
 */
struct CompiledHankin {
  ArenaVector<Vector> base_vertices;  /**< Input mesh corner vertices. */
  ArenaVector<Vector> static_vertices;  /**< Edge midpoints; angle-independent. */
  ArenaVector<Vector> dynamic_vertices; /**< Star points; recomputed per angle. */
  ArenaVector<HankinInstruction> dynamic_instructions; /**< One per dynamic vertex. */
  ArenaVector<uint8_t> face_counts; /**< Vertex count of each output face. */
  ArenaVector<uint16_t> faces;      /**< Flat vertex indices for all faces. */
  int static_offset = 0; /**< static_vertices.size(); base for dynamic indices in faces. */

  /**
   * @brief Resets all owned vectors to empty, releasing their contents.
   */
  void clear() {
    base_vertices.clear();
    static_vertices.clear();
    dynamic_vertices.clear();
    dynamic_instructions.clear();
    face_counts.clear();
    faces.clear();
  }

  /**
   * @brief Deep-copies all owned data into a target arena.
   * @param src Source instance to copy from.
   * @param dst Destination instance whose vectors are bound and filled.
   * @param arena Arena backing the destination's freshly bound vectors.
   * @details Required by Cloneable; each vector is rebound from @p arena and
   * bulk-copied (all element types are trivially copyable).
   */
  HS_COLD_MEMBER static void clone(const CompiledHankin &src, CompiledHankin &dst,
                    Arena &arena) {
    copy_vector(dst.base_vertices, src.base_vertices.data(),
                src.base_vertices.size(), arena);
    copy_vector(dst.static_vertices, src.static_vertices.data(),
                src.static_vertices.size(), arena);
    copy_vector(dst.dynamic_vertices, src.dynamic_vertices.data(),
                src.dynamic_vertices.size(), arena);
    copy_vector(dst.dynamic_instructions, src.dynamic_instructions.data(),
                src.dynamic_instructions.size(), arena);
    copy_vector(dst.face_counts, src.face_counts.data(),
                src.face_counts.size(), arena);
    copy_vector(dst.faces, src.faces.data(), src.faces.size(), arena);
    dst.static_offset = src.static_offset;
  }
};

namespace MeshOps {

/**
 * @brief Bakes the angle-independent Hankin topology for a mesh.
 * @param mesh Input closed-manifold mesh to derive the pattern from.
 * @param compiled Output topology, allocated from @p target_arena.
 * @param target_arena Arena that backs the persistent compiled vectors.
 * @param temp_arena Arena holding the transient half-edge structures.
 * @details Builds a half-edge mesh, emits one shared midpoint per edge into
 * static_vertices, reserves one dynamic (star-point) slot per half-edge, and
 * records the star and rosette faces.
 */
HS_COLD static void compile_hankin(const PolyMesh &mesh, CompiledHankin &compiled,
                                    Arena &target_arena, Arena &temp_arena) {
  // Topology via accessors (borrowed-mode safe); vertices is always owned.
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();

  compiled.base_vertices.bind(target_arena, V);
  for (size_t i = 0; i < V; ++i) {
    compiled.base_vertices.push_back(mesh.vertices[i]);
  }
  // Static pool is I/2 midpoints; largest emitted index is (I/2)+(I-1). Guard
  // adds 1 to both sides to dodge the unsigned underflow of I-1 at I == 0.
  HS_CHECK((I / 2) + I <= static_cast<size_t>(INT16_MAX) + 1,
           "Hankin output vertex count exceeds int16_t index range "
           "(MAX_INDICES raised too high?)");
  compiled.static_vertices.bind(target_arena, I / 2);
  compiled.dynamic_vertices.bind(target_arena, I);
  compiled.dynamic_instructions.bind(target_arena, I);
  compiled.face_counts.bind(target_arena, F + V);
  compiled.faces.bind(target_arena, 4 * I);

  {
    ScratchScope temp_arena_guard(temp_arena);

    HalfEdgeMesh he_mesh(temp_arena, mesh);

    require_closed_manifold(he_mesh, "compile_hankin");

    uint16_t *he_to_midpoint_idx = temp_arena.allocate_n<uint16_t>(I);
    std::fill_n(he_to_midpoint_idx, I, HE_NONE);
    uint16_t *he_to_dynamic_idx = temp_arena.allocate_n<uint16_t>(I);
    std::fill_n(he_to_dynamic_idx, I, HE_NONE);

    // Return the shared midpoint index for a half-edge, lazily creating and
    // caching it (under both the edge and its pair) on first encounter.
    auto get_midpoint_idx = [&](uint16_t he_idx) {
      if (he_to_midpoint_idx[he_idx] != HE_NONE)
        return he_to_midpoint_idx[he_idx];
      HalfEdge &he = he_mesh.half_edges[he_idx];
      if (he_to_midpoint_idx[he.pair] != HE_NONE)
        return he_to_midpoint_idx[he.pair];

      // Closed manifold: prev (the interior face loop) is always valid, so the
      // edge endpoints are simply prev->vertex and he->vertex.
      Vector p_a = mesh.vertices[he_mesh.half_edges[he.prev].vertex];
      Vector p_b = mesh.vertices[he.vertex];
      Vector mid = (p_a + p_b) * 0.5f;
      mid = normalized_or(mid, p_a);

      compiled.static_vertices.push_back(mid);
      uint16_t idx = narrow_index(compiled.static_vertices.size() - 1);
      he_to_midpoint_idx[he_idx] = idx;
      he_to_midpoint_idx[he.pair] = idx;
      return idx;
    };

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      get_midpoint_idx(static_cast<uint16_t>(i));
    }

    compiled.static_offset = static_cast<int>(compiled.static_vertices.size());

    // Star faces
    for (size_t i = 0; i < he_mesh.faces.size(); ++i) {
      HEFace &face = he_mesh.faces[i];
      uint16_t he_idx = face.half_edge;
      uint16_t start_he = he_idx;
      int count = 0;
      int walked = 0;

      if (he_idx == HE_NONE)
        continue;

      do {
        HS_CHECK(walked++ < (int)he_mesh.half_edges.size(),
                 "hankin star-face walk exceeded half-edge count");
        count += 2;
        HalfEdge &curr_he = he_mesh.half_edges[he_idx];
        uint16_t prev_idx = curr_he.prev;

        int idx_m1 = get_midpoint_idx(prev_idx);
        int idx_m2 = get_midpoint_idx(he_idx);

        HalfEdge &prev_he = he_mesh.half_edges[prev_idx];

        // Closed manifold: prev_he->prev is always valid (interior face loop).
        uint16_t i_corner = prev_he.vertex;
        uint16_t i_prev = he_mesh.half_edges[prev_he.prev].vertex;
        uint16_t i_next = curr_he.vertex;

        compiled.dynamic_instructions.push_back({i_corner, i_prev, i_next,
                                                narrow_index(idx_m1),
                                                narrow_index(idx_m2)});

        // size() BEFORE emplace_back is the index the new vertex will occupy.
        uint16_t dyn_idx = narrow_index(compiled.dynamic_vertices.size());
        he_to_dynamic_idx[he_idx] = dyn_idx;
        compiled.dynamic_vertices.emplace_back();

        compiled.faces.push_back(narrow_index(idx_m1));
        compiled.faces.push_back(narrow_index(compiled.static_offset + dyn_idx));

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
    int16_t *face_indices = temp_arena.allocate_n<int16_t>(2 * I);

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      uint16_t he_start_idx = static_cast<uint16_t>(i);
      const HalfEdge &he_start = he_mesh.half_edges[he_start_idx];
      // Closed manifold: the rosette origin is he_start's tail vertex, read via
      // its always-valid prev (interior face loop).
      uint16_t origin_idx = he_mesh.half_edges[he_start.prev].vertex;
      if (visited_verts[origin_idx])
        continue;
      visited_verts[origin_idx] = true;

      uint16_t curr_idx = he_start_idx;
      uint16_t start_orbit = curr_idx;
      int count = 0;

      do {
        const HalfEdge &curr_he = he_mesh.half_edges[curr_idx];
        HS_CHECK(count < (int)(2 * I), "Hankin rosette winding overflow");
        face_indices[count++] = he_to_midpoint_idx[curr_idx];
        // Closed manifold: pair is always valid, so the orbit closes back on
        // start_orbit and never hits HE_NONE.
        uint16_t next_edge_idx = he_mesh.half_edges[curr_he.pair].next;
        HS_CHECK(count < (int)(2 * I), "Hankin rosette winding overflow");
        face_indices[count++] = narrow_index(
            compiled.static_offset + he_to_dynamic_idx[next_edge_idx]);
        curr_idx = next_edge_idx;
      } while (curr_idx != start_orbit);

      // count = 2 * vertex degree. Degree-2 is legal (hankin-of-hankin walks
      // its own degree-2 midpoints -> quad rosette); only degree < 2 degenerates.
      HS_CHECK(count >= 4, "Hankin rosette winding has degree < 2");
      compiled.face_counts.push_back(narrow_face_count(count));
      for (int k = count - 1; k >= 0; --k) {
        compiled.faces.push_back(face_indices[k]);
      }
    }
  }
}

/**
 * @brief Positions the angle-dependent vertices and writes the final mesh.
 * @tparam MeshT Output mesh type; may optionally provide face_offsets.
 * @param compiled Baked topology whose dynamic_vertices are recomputed in place.
 * @param out_mesh Output mesh, allocated from @p target_arena.
 * @param target_arena Arena backing @p out_mesh's vertex and face vectors.
 * @param angle Contact angle in radians. At ~0 the star points collapse onto
 *   their corners (flat tiling); larger angles push the rays out so the rays
 *   from adjacent edges intersect to form sharper star points.
 * @details A corner whose contact planes are near-parallel has no nearby
 *   intersection; its star point falls back to the edge-midpoint mean instead
 *   of being flung across the sphere (see STAR_FAR_RATIO_SQ).
 */
/** Squared ceiling on chord(star point, corner) / chord(corner, midpoint);
 * measured healthy maximum is ~5.5, degenerate intersections start at ~433. */
inline constexpr float STAR_FAR_RATIO_SQ = 36.0f;

template <typename MeshT>
HS_COLD_MEMBER inline void update_hankin(CompiledHankin &compiled, MeshT &out_mesh,
                          Arena &target_arena, float angle) {

  // Drop any borrowed-mode views a reused MeshState may still carry, so the
  // owned topology bound below is unambiguous on reuse.
  if constexpr (requires { out_mesh.set_owned(); }) {
    out_mesh.set_owned();
  }

  bool is_flat = std::abs(angle) < math::TOLERANCE;

  float cos_ha = cosf(angle * 0.5f);
  float sin_ha = sinf(angle * 0.5f);

  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &instr = compiled.dynamic_instructions[i];
    Vector p_corner = compiled.base_vertices[instr.v_corner];

    if (is_flat) {
      compiled.dynamic_vertices[i] = normalized_or(p_corner, p_corner);
      continue;
    }

    Vector m1 = compiled.static_vertices[instr.idx_m1];
    Vector m2 = compiled.static_vertices[instr.idx_m2];
    Vector p_prev = compiled.base_vertices[instr.v_prev];
    Vector p_next = compiled.base_vertices[instr.v_next];

    Vector cross1 = cross(p_prev, p_corner);
    Vector cross2 = cross(p_corner, p_next);

    if (dot(cross1, cross1) < math::EPS_CROSS_SQ ||
        dot(cross2, cross2) < math::EPS_CROSS_SQ) {
      compiled.dynamic_vertices[i] = normalized_or(p_corner, p_corner);
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
    Vector cn = normalized_or(p_corner, p_corner);
    bool degenerate = dot(intersect, intersect) < math::EPS_LEN_SQ;
    if (!degenerate) {
      if (dot(intersect, p_corner) < 0)
        intersect = -intersect;
      intersect = intersect.normalized();
      // Near-parallel contact planes fling the intersection geodesically far
      // from the corner, yielding a sliver face that renders as a long line.
      float local_sq = std::max(distance_squared(m1, cn),
                                distance_squared(m2, cn));
      degenerate =
          distance_squared(intersect, cn) > STAR_FAR_RATIO_SQ * local_sq;
    }
    if (degenerate) {
      Vector fallback = normalized_or(m1 + m2, cn);
      if (dot(fallback, p_corner) < 0)
        fallback = -fallback;
      compiled.dynamic_vertices[i] = fallback;
      continue;
    }

    compiled.dynamic_vertices[i] = intersect;
  }

  out_mesh.vertices.bind(target_arena,
                               compiled.static_vertices.size() +
                                   compiled.dynamic_vertices.size());
  out_mesh.vertices.append_bulk(compiled.static_vertices.data(),
                                compiled.static_vertices.size());
  out_mesh.vertices.append_bulk(compiled.dynamic_vertices.data(),
                                compiled.dynamic_vertices.size());

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
 * @param angle Contact angle in radians (see update_hankin).
 * @return The generated Hankin PolyMesh, allocated from @p target.
 * @details Convenience wrapper for callers that do not vary the angle. To
 * animate the angle, keep a CompiledHankin and call update_hankin repeatedly.
 */
HS_COLD static PolyMesh hankin(const PolyMesh &mesh, Arena &target, Arena &temp,
                                float angle) {
  PolyMesh out;

  {
    ScratchScope temp_guard(temp);
    CompiledHankin compiled;
    // Arena polarity is reversed from the streaming path: the throwaway
    // CompiledHankin is allocated from `temp` while `target` serves as
    // compile_hankin's working arena, then update_hankin builds `out` into it.
    compile_hankin(mesh, compiled, temp, target);
    update_hankin(compiled, out, target, angle);
  }

  return out;
}

} // namespace MeshOps
