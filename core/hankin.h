/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh.h"

/**
 * @brief Structure returned by compile_hankin.
 */
struct HankinInstruction {
  uint16_t v_corner; /**< Index to base_vertices for corner vertex. */
  uint16_t v_prev;   /**< Index to base_vertices for previous vertex. */
  uint16_t v_next;   /**< Index to base_vertices for next vertex. */
  uint16_t idx_m1;   /**< Index of first midpoint (static vertex). */
  uint16_t idx_m2;   /**< Index of second midpoint (static vertex). */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 */
struct CompiledHankin {
  ArenaVector<Vector> base_vertices;
  ArenaVector<Vector> static_vertices;
  ArenaVector<Vector> dynamic_vertices;
  ArenaVector<HankinInstruction> dynamic_instructions;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<uint16_t> faces;
  int static_offset = 0;

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

/**
 * @brief Operations on meshes (Hankin pattern operators).
 */
namespace MeshOps {

/**
 * @brief Compiles the topology for a Hankin pattern.
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
  compiled.static_vertices.bind(target_arena, (I / 2) + 1);
  compiled.dynamic_vertices.bind(target_arena, I);
  compiled.dynamic_instructions.bind(target_arena, I);
  compiled.face_counts.bind(target_arena, F + V);
  compiled.faces.bind(target_arena, 4 * I);

  {
    ScratchScope _(temp_arena);

    HalfEdgeMesh he_mesh(temp_arena, mesh);
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
      if (he.pair != HE_NONE && he_to_midpoint_idx[he.pair] != HE_NONE)
        return he_to_midpoint_idx[he.pair];

      // A closed-manifold half-edge always has a valid prev (interior face
      // loop). If both prev and pair are HE_NONE the edge is degenerate /
      // non-manifold and the he.pair fallback below would read
      // half_edges[HE_NONE] out of bounds — trap on the bad input instead.
      HS_CHECK(he.prev != HE_NONE || he.pair != HE_NONE);
      Vector p_a = he.prev != HE_NONE
                      ? mesh.vertices[he_mesh.half_edges[he.prev].vertex]
                      : mesh.vertices[he_mesh.half_edges[he.pair].vertex];
      Vector p_b = mesh.vertices[he.vertex];
      Vector mid = (p_a + p_b) * 0.5f;
      mid.normalize();

      compiled.static_vertices.push_back(mid);
      uint16_t idx = static_cast<uint16_t>(compiled.static_vertices.size() - 1);
      he_to_midpoint_idx[he_idx] = idx;
      if (he.pair != HE_NONE)
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

        uint16_t i_corner = prev_he.vertex;
        uint16_t i_prev = prev_he.prev != HE_NONE
                             ? he_mesh.half_edges[prev_he.prev].vertex
                             : he_mesh.half_edges[prev_he.pair].vertex;
        uint16_t i_next = curr_he.vertex;

        compiled.dynamic_instructions.push_back({i_corner, i_prev, i_next,
                                                static_cast<uint16_t>(idx_m1),
                                                static_cast<uint16_t>(idx_m2)});

        int16_t dyn_idx = static_cast<int16_t>(compiled.dynamic_vertices.size());
        he_to_dynamic_idx[he_idx] = dyn_idx;
        compiled.dynamic_vertices.emplace_back();

        compiled.faces.push_back(idx_m1);
        compiled.faces.push_back(compiled.static_offset + dyn_idx);

        he_idx = curr_he.next;
      } while (he_idx != HE_NONE && he_idx != start_he);

      compiled.face_counts.push_back(count);
    }

    // Rosette faces
    bool *visited_verts = static_cast<bool *>(
        temp_arena.allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visited_verts, V, false);

    // Per-orbit scratch buffer sized to the absolute upper bound on valence
    // (= total half-edges). Replaces a fixed `face_indices[100]` stack array
    // that silently truncated high-valence orbits.
    int16_t *face_indices = static_cast<int16_t *>(
        temp_arena.allocate(2 * I * sizeof(int16_t), alignof(int16_t)));

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      uint16_t he_start_idx = static_cast<uint16_t>(i);
      const HalfEdge &he_start = he_mesh.half_edges[he_start_idx];
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
        uint16_t next_edge_idx = curr_he.pair != HE_NONE
                                   ? he_mesh.half_edges[curr_he.pair].next
                                   : HE_NONE;
        if (next_edge_idx == HE_NONE)
          break;
        HS_CHECK(count < (int)(2 * I));
        face_indices[count++] =
            compiled.static_offset + he_to_dynamic_idx[next_edge_idx];
        curr_idx = next_edge_idx;
      } while (curr_idx != HE_NONE && curr_idx != start_orbit);

      if (count > 2) {
        compiled.face_counts.push_back(count);
        for (int k = count - 1; k >= 0; --k) {
          compiled.faces.push_back(face_indices[k]);
        }
      }
    }
  }
}

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
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.bind(target_arena, compiled.faces.size());
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    out_mesh.faces.push_back(compiled.faces[i]);
}

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
