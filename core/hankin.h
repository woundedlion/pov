/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh.h"

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

  void clear() {
    base_vertices.clear();
    static_vertices.clear();
    dynamic_vertices.clear();
    dynamic_instructions.clear();
    face_counts.clear();
    faces.clear();
  }

  /// Deep-copy all owned data into a target arena. Required by Cloneable.
  static void clone(const CompiledHankin &src, CompiledHankin &dst,
                    Arena &arena) {
    auto push = [&arena](const auto &s_vec, auto &d_vec) {
      d_vec.bind(arena, s_vec.size());
      for (size_t i = 0; i < s_vec.size(); ++i)
        d_vec.push_back(s_vec[i]);
    };
    push(src.base_vertices, dst.base_vertices);
    push(src.static_vertices, dst.static_vertices);
    push(src.dynamic_vertices, dst.dynamic_vertices);
    push(src.dynamic_instructions, dst.dynamic_instructions);
    push(src.face_counts, dst.face_counts);
    push(src.faces, dst.faces);
    dst.static_offset = src.static_offset;
  }
};

namespace MeshOps {

/**
 * @brief Bakes the angle-independent Hankin topology for a mesh.
 *
 * Builds a half-edge mesh, emits one shared midpoint per edge into
 * static_vertices, reserves one dynamic (star-point) slot per half-edge, and
 * records the star and rosette faces. @p compiled is allocated from
 * @p target_arena; @p temp_arena holds the transient half-edge structures.
 */
template <typename MeshT>
FLASHMEM static void compile_hankin(const MeshT &mesh, CompiledHankin &compiled,
                                    Arena &target_arena, Arena &temp_arena) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  compiled.base_vertices.bind(target_arena, V);
  for (size_t i = 0; i < V; ++i) {
    compiled.base_vertices.push_back(mesh.vertices[i]);
  }
  // compile_hankin requires a closed manifold (every half-edge paired): each
  // undirected edge is shared by 2 half-edges that emit one shared midpoint, so
  // the static pool is exactly I/2 vertices and I is even. The entry check below
  // enforces this; a boundary mesh would otherwise overrun the pool. Combined
  // with the I dynamic vertices, the index (static_offset + dyn_idx, emitted
  // below) must fit the int16_t topology range. narrow_index() would trap an
  // overflow mid-build, but the ceiling is knowable here — fail fast at the
  // binding site if MAX_INDICES is ever raised past what the index width allows.
  HS_CHECK((I / 2) + I <= static_cast<size_t>(INT16_MAX),
           "Hankin output vertex count exceeds int16_t index range "
           "(MAX_INDICES raised too high?)");
  compiled.static_vertices.bind(target_arena, I / 2);
  compiled.dynamic_vertices.bind(target_arena, I);
  compiled.dynamic_instructions.bind(target_arena, I);
  compiled.face_counts.bind(target_arena, F + V);
  compiled.faces.bind(target_arena, 4 * I);

  {
    ScratchScope _(temp_arena);

    HalfEdgeMesh he_mesh(temp_arena, mesh);

    // Closed-manifold invariant: every half-edge must be paired. The pool
    // sizing (I/2 midpoints) and the orbit walks below all assume each edge has
    // its twin; a boundary half-edge would overrun static_vertices and break the
    // rosette traversal. Hankin only ever runs on closed convex solids, so an
    // unpaired edge is a build-invariant violation — trap with a clear message.
    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      HS_CHECK(he_mesh.half_edges[i].pair != HE_NONE,
               "compile_hankin requires a closed manifold (unpaired half-edge)");
    }

    uint16_t *he_to_midpoint_idx = static_cast<uint16_t *>(
        temp_arena.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_midpoint_idx, I, HE_NONE);
    uint16_t *he_to_dynamic_idx = static_cast<uint16_t *>(
        temp_arena.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_dynamic_idx, I, HE_NONE);

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
      mid.normalize();

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

      if (he_idx == HE_NONE)
        continue;

      do {
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
    bool *visited_verts = static_cast<bool *>(
        temp_arena.allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visited_verts, V, false);

    // Per-orbit scratch buffer. Each orbit step appends two indices (a
    // midpoint and a dynamic vertex), so the absolute upper bound on entries
    // is twice the total half-edge count.
    int16_t *face_indices = static_cast<int16_t *>(
        temp_arena.allocate(2 * I * sizeof(int16_t), alignof(int16_t)));

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
        HS_CHECK(count < (int)(2 * I));
        face_indices[count++] = he_to_midpoint_idx[curr_idx];
        // Closed manifold: pair is always valid, so the orbit closes back on
        // start_orbit and never hits HE_NONE.
        uint16_t next_edge_idx = he_mesh.half_edges[curr_he.pair].next;
        HS_CHECK(count < (int)(2 * I));
        face_indices[count++] =
            compiled.static_offset + he_to_dynamic_idx[next_edge_idx];
        curr_idx = next_edge_idx;
      } while (curr_idx != start_orbit);

      if (count > 2) {
        compiled.face_counts.push_back(narrow_face_count(count));
        for (int k = count - 1; k >= 0; --k) {
          compiled.faces.push_back(face_indices[k]);
        }
      }
    }
  }
}

/**
 * @brief Positions the angle-dependent vertices and writes the final mesh.
 *
 * @param angle Contact angle in radians. At ~0 the star points collapse onto
 *   their corners (flat tiling); larger angles push the rays out so the rays
 *   from adjacent edges intersect to form sharper star points. @p out_mesh is
 *   allocated from @p target_arena.
 */
template <typename MeshT>
inline void update_hankin(CompiledHankin &compiled, MeshT &out_mesh,
                          Arena &target_arena, float angle) {

  bool is_flat = std::abs(angle) < math::TOLERANCE;

  // Precompute half-angle trig: one cosf + sinf instead of 2N
  float cos_ha = cosf(angle * 0.5f);
  float sin_ha = sinf(angle * 0.5f);

  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &instr = compiled.dynamic_instructions[i];
    Vector p_corner = compiled.base_vertices[instr.v_corner];

    if (is_flat) {
      compiled.dynamic_vertices[i] = p_corner.normalized();
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
      compiled.dynamic_vertices[i] = p_corner.normalized();
      continue; // zero length edge
    }

    Vector n_edge1 = cross1.normalized();
    // Inline make_rotation using precomputed half-angle trig
    Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
    Vector n_hankin1 = rotate(n_edge1, q1);

    Vector n_edge2 = cross2.normalized();
    // cos(-x) = cos(x), sin(-x) = -sin(x)
    Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
    Vector n_hankin2 = rotate(n_edge2, q2);

    Vector intersect = cross(n_hankin1, n_hankin2);
    float len_sq = dot(intersect, intersect);
    if (len_sq < math::EPS_LEN_SQ)
      intersect = (m1 + m2).normalized();
    if (dot(intersect, p_corner) < 0)
      intersect = -intersect;

    compiled.dynamic_vertices[i] = intersect.normalized();
  }

  out_mesh.vertices.bind(target_arena,
                               compiled.static_vertices.size() +
                                   compiled.dynamic_vertices.size());
  for (size_t i = 0; i < compiled.static_vertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.static_vertices[i]);
  for (size_t i = 0; i < compiled.dynamic_vertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.dynamic_vertices[i]);

  out_mesh.face_counts.bind(target_arena, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.bind(target_arena, compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      // face_offsets is uint16_t; current_offset (cumulative index count) wraps
      // silently past 65535 — trap, mirroring MeshOps::compile.
      HS_CHECK(current_offset <= UINT16_MAX,
               "mesh face_offsets exceeds 16-bit index range");
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.bind(target_arena, compiled.faces.size());
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    out_mesh.faces.push_back(compiled.faces[i]);
}

/**
 * @brief One-shot Hankin pattern: compile then update, returning the new mesh.
 *
 * Convenience wrapper for callers that do not vary the angle. The compiled
 * topology is built in @p temp and discarded; only @p out is allocated from
 * @p target. To animate the angle, keep a CompiledHankin and call
 * update_hankin repeatedly instead.
 */
FLASHMEM static PolyMesh hankin(const PolyMesh &mesh, Arena &target, Arena &temp,
                                float angle) {
  PolyMesh out;

  {
    ScratchScope _(temp);
    CompiledHankin compiled;
    compile_hankin(mesh, compiled, temp, target);
    update_hankin(compiled, out, target, angle);
  }

  return out;
}

} // namespace MeshOps
