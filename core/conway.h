/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh.h"

/**
 * @brief Operations on meshes (Conway operators and utilities).
 */
namespace MeshOps {

// `narrow_index` (the trapping size_t -> uint16_t topology-index cast) lives in
// mesh.h so the Conway and Hankin operators share a single guarded narrowing.

/**
 * @brief Compute the centroid of a face by walking its half-edge loop.
 */
template <typename MeshT>
inline Vector face_centroid(const HalfEdgeMesh &he_mesh,
                            const MeshT &mesh, size_t face_index,
                            int &out_count) {
  const HEFace &face = he_mesh.faces[face_index];
  Vector c(0, 0, 0);
  out_count = 0;
  uint16_t he_idx = face.half_edge;
  uint16_t start = he_idx;
  if (he_idx != HE_NONE) {
    do {
      c = c + mesh.vertices[he_mesh.half_edges[he_idx].vertex];
      out_count++;
      he_idx = he_mesh.half_edges[he_idx].next;
    } while (he_idx != HE_NONE && he_idx != start);
  }
  if (out_count > 0)
    c = c / static_cast<float>(out_count);
  return c;
}

/**
 * @brief Newell's method face normal — robust to non-planar faces and
 * collinear vertex triplets. Walks the half-edge loop, summing the Newell
 * contribution of each edge (consecutive vertex pair). Returns the
 * unnormalized normal; caller normalizes if needed.
 */
template <typename MeshT>
inline Vector face_normal(const HalfEdgeMesh &he_mesh, const MeshT &mesh,
                          size_t face_index) {
  const HEFace &face = he_mesh.faces[face_index];
  Vector n(0, 0, 0);
  uint16_t he_idx = face.half_edge;
  if (he_idx == HE_NONE) return n;
  uint16_t start = he_idx;
  do {
    const HalfEdge &he = he_mesh.half_edges[he_idx];
    const Vector &curr = mesh.vertices[he.vertex];
    const Vector &next = mesh.vertices[he_mesh.half_edges[he.next].vertex];
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
    he_idx = he.next;
  } while (he_idx != start);
  return n;
}

/**
 * @brief Walk all half-edges orbiting a vertex, calling visitor(curr_idx)
 *        for each. Returns the number of visited half-edges.
 *
 * OrbitMode selects the traversal direction:
 *   'P' = prev->pair (dual, ambo, truncate)
 *   'N' = pair->next  (expand, snub)
 */
template <char OrbitMode, typename VisitorFn>
inline int vertex_orbit(const HalfEdgeMesh &he_mesh, uint16_t start_idx,
                        VisitorFn &&visitor) {
  uint16_t curr_idx = start_idx;
  int count = 0;
  // Hard upper bound, independent of asserts (survives NDEBUG): a vertex orbit
  // visits each incident half-edge at most once, so it can never legitimately
  // touch more half-edges than exist. Exceeding that means the half-edge graph
  // is non-manifold/corrupt and the walk would otherwise spin forever — trap so
  // it's caught on the bench instead of hanging the device.
  const int max_orbit = static_cast<int>(he_mesh.half_edges.size());
  do {
    HS_CHECK(count < max_orbit);
    const HalfEdge &curr_he = he_mesh.half_edges[curr_idx];
    if (curr_he.face == HE_NONE)
      break;

    visitor(curr_idx);
    count++;

    if constexpr (OrbitMode == 'P') {
      // prev->pair orbit
      if (he_mesh.half_edges[curr_he.prev].pair == HE_NONE)
        break;
      curr_idx = he_mesh.half_edges[curr_he.prev].pair;
    } else {
      // pair->next orbit
      if (curr_he.pair == HE_NONE)
        break;
      curr_idx = he_mesh.half_edges[curr_he.pair].next;
    }
  } while (curr_idx != HE_NONE && curr_idx != start_idx);
  return count;
}

/**
 * @brief Emit one output face per source-vertex orbit (the shared dual/ambo/
 *        expand/snub scaffold).
 *
 * Walks every half-edge once, and for each not-yet-visited origin vertex builds
 * its orbit via vertex_orbit<DIR>, mapping each visited half-edge to an output
 * vertex index through `value_of(idx)`, then emits the collected indices as one
 * face (skipping degenerate <3-gons). `reverse` flips the winding for the
 * pair->next ('N') orbits whose natural order is opposite the desired front.
 *
 * @param visited_verts caller-allocated bool[V] scratch (filled here).
 * @param orbit_buf      caller-allocated uint16_t[I] scratch.
 */
template <char DIR, typename ValueFn>
inline void emit_vertex_orbit_faces(const HalfEdgeMesh &he_mesh,
                                    PolyMesh &out_mesh, bool *visited_verts,
                                    uint16_t *orbit_buf, size_t V, size_t I,
                                    bool reverse, ValueFn &&value_of) {
  std::fill_n(visited_verts, V, false);
  for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
    uint16_t he_start_idx = static_cast<uint16_t>(i);
    const HalfEdge &he_start = he_mesh.half_edges[he_start_idx];

    uint16_t origin_idx = he_mesh.half_edges[he_start.prev].vertex;
    if (visited_verts[origin_idx])
      continue;
    visited_verts[origin_idx] = true;

    int orbit_count = 0;
    vertex_orbit<DIR>(he_mesh, he_start_idx, [&](uint16_t idx) {
      HS_CHECK(orbit_count < (int)I);
      orbit_buf[orbit_count++] = static_cast<uint16_t>(value_of(idx));
    });

    if (orbit_count >= 3) {
      out_mesh.face_counts.push_back(orbit_count);
      if (reverse) {
        for (int k = orbit_count - 1; k >= 0; --k)
          out_mesh.faces.push_back(orbit_buf[k]);
      } else {
        for (int k = 0; k < orbit_count; ++k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }
  }
}

/**
 * @brief Copies a mesh's topology (views) and vertices into a target mesh.
 * Base case of the variadic overload below: no transformers, so vertices are
 * copied verbatim.
 * @param local_state The source mesh state (centered around origin).
 * @param world_state The destination mesh state to populate.
 * @param arena The memory arena to allocate vertices.
 */
template <typename MeshT>
inline void transform(const MeshT &local_state, MeshT &world_state,
                      Arena& arena) {
  world_state.face_counts_view = ArenaSpan(local_state.face_counts);
  world_state.faces_view = ArenaSpan(local_state.faces);

  if constexpr (requires { world_state.face_offsets_view; }) {
    world_state.face_offsets_view = ArenaSpan(local_state.face_offsets);
  }

  world_state.vertices.bind(arena, local_state.vertices.size());

  for (size_t i = 0; i < local_state.vertices.size(); ++i) {
    world_state.vertices.push_back(local_state.vertices[i]);
  }
}

template <typename MeshT, typename T1, typename... Transformers>
inline void transform(const MeshT &mesh, MeshT &transformed, Arena& arena,
                      const T1 &first_transformer,
                      const Transformers &...transformers) {
  transformed.face_counts_view = ArenaSpan(mesh.face_counts);
  transformed.faces_view = ArenaSpan(mesh.faces);

  if constexpr (requires { transformed.face_offsets_view; }) {
    transformed.face_offsets_view = ArenaSpan(mesh.face_offsets);
  }

  transformed.vertices.bind(arena, mesh.vertices.size());

  for (size_t i = 0; i < mesh.vertices.size(); ++i) {
    Vector v = mesh.vertices[i];

    // Apply first transformer
    v = first_transformer(v);

    // Unroll remaining transformers at compile time using a fold expression
    (..., (v = transformers(v)));

    transformed.vertices.push_back(v);
  }
}

// ---------------------------------------------------------------------------
// Conway operators
//
// All operators take a const MeshT& input (PolyMesh or MeshState — both
// expose the same vertices/face_counts/faces shape) and return a fresh
// PolyMesh in `target`. `temp` is used for scratch (HalfEdgeMesh build,
// per-orbit index buffers); both arenas are checkpointed via ScratchScope.
//
// Per-vertex orbit construction uses an arena scratch buffer sized to the
// maximum possible valence (= total half-edges), so high-valence vertices
// never overflow.
// ---------------------------------------------------------------------------

/**
 * @brief Normalizes all vertices in the mesh to the unit sphere.
 */
template <typename MeshT> static void normalize(MeshT &mesh) {
  for (auto &v : mesh.vertices) {
    v.normalize();
  }
}

/**
 * @brief Computes the dual of a mesh.
 */
template <typename MeshT>
FLASHMEM static PolyMesh dual(const MeshT &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();

  out_mesh.vertices.bind(target, F);
  out_mesh.face_counts.bind(target, V);
  out_mesh.faces.bind(target, I);

  {
    ScratchScope _temp(temp);
    ScratchScope _target(target);

    HalfEdgeMesh he_mesh(temp, mesh);

    for (size_t i = 0; i < he_mesh.faces.size(); ++i) {
      int count;
      Vector c = face_centroid(he_mesh, mesh, i, count);
      out_mesh.vertices.push_back(c.normalized());
    }

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    // Per-orbit scratch buffer sized to the absolute upper bound on valence.
    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    // One face per source vertex: its orbit's incident faces become the dual
    // face's vertices.
    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) { return he_mesh.half_edges[idx].face; });
  }
  return out_mesh;
}

/**
 * @brief Kis operator: Raises a pyramid on each face.
 */
template <typename MeshT>
FLASHMEM static PolyMesh kis(const MeshT &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  const auto *face_counts = mesh.get_face_counts_data();
  const auto *faces = mesh.get_faces_data();

  out_mesh.vertices.bind(target, V + F);
  out_mesh.face_counts.bind(target, I);
  out_mesh.faces.bind(target, 3 * I);

  {
    ScratchScope _(temp);

    for (size_t i = 0; i < V; ++i)
      out_mesh.vertices.push_back(mesh.vertices[i]);

    size_t offset = 0;
    for (size_t fi = 0; fi < F; ++fi) {
      int count = face_counts[fi];
      Vector centroid(0, 0, 0);
      for (int i = 0; i < count; ++i) {
        centroid = centroid + mesh.vertices[faces[offset + i]];
      }
      if (count > 0)
        centroid = centroid / static_cast<float>(count);

      out_mesh.vertices.push_back(centroid);
      int center_idx = narrow_index(out_mesh.vertices.size() - 1);

      for (int i = 0; i < count; ++i) {
        uint16_t vi = faces[offset + i];
        uint16_t vj = faces[offset + (i + 1) % count];
        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(vi);
        out_mesh.faces.push_back(vj);
        out_mesh.faces.push_back(static_cast<uint16_t>(center_idx));
      }
      offset += count;
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Ambo operator: Truncates vertices to edge midpoints.
 */
template <typename MeshT>
FLASHMEM static PolyMesh ambo(const MeshT &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, E);
  out_mesh.face_counts.bind(target, F + V);
  out_mesh.faces.bind(target, 2 * I);

  {
    ScratchScope _temp(temp);
    ScratchScope _target(target);

    HalfEdgeMesh he_mesh(temp, mesh);

    uint16_t *edge_to_vert =
        static_cast<uint16_t *>(target.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(edge_to_vert, I, HE_NONE);

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    // 3. Populate Vertices
    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      if (edge_to_vert[i] == HE_NONE) {
        const HalfEdge &he = he_mesh.half_edges[i];
        uint16_t v1 = he.vertex;
        uint16_t v2 = he_mesh.half_edges[he.prev].vertex;

        Vector mid = (mesh.vertices[v1] + mesh.vertices[v2]) * 0.5f;
        out_mesh.vertices.push_back(mid);
        uint16_t new_idx = narrow_index(out_mesh.vertices.size() - 1);

        edge_to_vert[i] = new_idx;
        if (he.pair != HE_NONE)
          edge_to_vert[he.pair] = new_idx;
      }
    }

    // 4. Reconstruct Original Faces (Shrunk)
    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      const HEFace &face = he_mesh.faces[fi];
      uint16_t he_idx = face.half_edge;
      int count = 0;
      if (he_idx != HE_NONE) {
        uint16_t start = he_idx;
        do {
          count++;
          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);

        if (count >= 3) {
          out_mesh.face_counts.push_back(count);
          he_idx = start;
          do {
            out_mesh.faces.push_back(edge_to_vert[he_idx]);
            he_idx = he_mesh.half_edges[he_idx].next;
          } while (he_idx != HE_NONE && he_idx != start);
        }
      }
    }

    // 5. Build Vertex Orbits (New Faces)
    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) { return edge_to_vert[idx]; });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Truncate operator: Cuts corners off the polyhedron.
 * @param t Truncation depth, the fraction along each edge at which the two cut
 *   points sit, in [0..1]. Each edge `(k1,k2)` yields `k1+(k2-k1)*t` and
 *   `k2+(k1-k2)*t`. For `t<0.5` the cut points stay on their own half; at
 *   exactly `0.5` both reach the midpoint and this short-circuits to `ambo`;
 *   for `t>0.5` the two points cross past each other, producing intentional
 *   self-intersecting cut faces (used by the `*_truncate50d_*` solids). `t`
 *   outside `[0..1]` would place a cut point beyond the edge endpoints, so it
 *   traps per the fail-fast doctrine.
 */
template <typename MeshT>
FLASHMEM static PolyMesh truncate(const MeshT &mesh, Arena &target, Arena &temp,
                                  float t = 0.25f) {
  HS_CHECK(t >= 0.0f && t <= 1.0f);
  if (std::abs(t - 0.5f) < math::TOLERANCE) {
    return ambo(mesh, target, temp);
  }

  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, 2 * E);
  out_mesh.face_counts.bind(target, F + V);
  out_mesh.faces.bind(target, 3 * I);

  {
    ScratchScope _temp(temp);
    ScratchScope _target(target);

    HalfEdgeMesh he_mesh(temp, mesh);

    std::pair<int16_t, int16_t> *edge_to_vert =
        static_cast<std::pair<int16_t, int16_t> *>(
            target.allocate(
                I * sizeof(std::pair<int16_t, int16_t>),
                alignof(std::pair<int16_t, int16_t>)));
    std::fill_n(edge_to_vert, I, std::make_pair<int16_t, int16_t>(-1, -1));

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      if (edge_to_vert[i].first == -1) {
        const HalfEdge &he = he_mesh.half_edges[i];
        uint16_t vi = he_mesh.half_edges[he.prev].vertex;
        uint16_t vj = he.vertex;
        uint16_t k1 = std::min(vi, vj);
        uint16_t k2 = std::max(vi, vj);

        Vector new_u =
            mesh.vertices[k1] + (mesh.vertices[k2] - mesh.vertices[k1]) * t;
        Vector new_v =
            mesh.vertices[k2] + (mesh.vertices[k1] - mesh.vertices[k2]) * t;

        out_mesh.vertices.push_back(new_u);
        uint16_t idx_u = narrow_index(out_mesh.vertices.size() - 1);

        out_mesh.vertices.push_back(new_v);
        uint16_t idx_v = narrow_index(out_mesh.vertices.size() - 1);

        edge_to_vert[i] = {idx_u, idx_v};
        if (he.pair != HE_NONE)
          edge_to_vert[he.pair] = {idx_u, idx_v};
      }
    }

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      const HEFace &face = he_mesh.faces[fi];
      uint16_t he_idx = face.half_edge;
      int count = 0;
      if (he_idx != HE_NONE) {
        uint16_t start = he_idx;
        do {
          count++;
          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);

        if (count >= 3) {
          out_mesh.face_counts.push_back(count * 2);
          he_idx = start;
          do {
            const HalfEdge &he = he_mesh.half_edges[he_idx];
            uint16_t vi = he_mesh.half_edges[he.prev].vertex;
            uint16_t vj = he.vertex;
            uint16_t k1 = std::min(vi, vj);
            std::pair<int16_t, int16_t> new_verts = edge_to_vert[he_idx];

            if (vi == k1) {
              out_mesh.faces.push_back(new_verts.first);
              out_mesh.faces.push_back(new_verts.second);
            } else {
              out_mesh.faces.push_back(new_verts.second);
              out_mesh.faces.push_back(new_verts.first);
            }

            he_idx = he_mesh.half_edges[he_idx].next;
          } while (he_idx != HE_NONE && he_idx != start);
        }
      }
    }

    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) {
          const HalfEdge &curr_he = he_mesh.half_edges[idx];
          uint16_t vi = he_mesh.half_edges[curr_he.prev].vertex;
          uint16_t vj = curr_he.vertex;
          uint16_t k1 = std::min(vi, vj);
          std::pair<int16_t, int16_t> new_verts = edge_to_vert[idx];
          return (vi == k1) ? new_verts.first : new_verts.second;
        });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Expand operator: Separates faces (e = aa).
 * @param t Expansion factor. Default 2-sqrt(2) ~= 0.5857.
 */
template <typename MeshT>
FLASHMEM static PolyMesh expand(const MeshT &mesh, Arena &target, Arena &temp,
                                float t = 2.0f - sqrtf(2.0f)) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, I);
  out_mesh.face_counts.bind(target, F + V + E);
  out_mesh.faces.bind(target, 4 * I);

  {
    ScratchScope _temp(temp);
    ScratchScope _target(target);

    HalfEdgeMesh he_mesh(temp, mesh);
    int16_t *he_to_vert_idx = static_cast<int16_t *>(
        target.allocate(I * sizeof(int16_t), alignof(int16_t)));
    std::fill_n(he_to_vert_idx, I, -1);

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    bool *visited_edges = static_cast<bool *>(
        target.allocate(I * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      out_mesh.face_counts.push_back(count);
      he_idx = start;
      if (he_idx != HE_NONE) {
        do {
          Vector v = mesh.vertices[he_mesh.half_edges[he_idx].vertex];
          Vector new_v = v + (centroid - v) * t;
          out_mesh.vertices.push_back(new_v);
          int idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_vert_idx[he_idx] = idx;

          out_mesh.faces.push_back(idx);
          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    emit_vertex_orbit_faces<'N'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/true,
        [&](uint16_t idx) {
          return he_to_vert_idx[he_mesh.half_edges[idx].prev];
        });

    std::fill_n(visited_edges, I, false);

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      uint16_t he_idx = static_cast<uint16_t>(i);
      const HalfEdge &he = he_mesh.half_edges[he_idx];

      if (visited_edges[he_idx])
        continue;

      visited_edges[he_idx] = true;
      if (he.pair != HE_NONE)
        visited_edges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(4);
        out_mesh.faces.push_back(he_to_vert_idx[he.prev]);
        out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
        out_mesh.faces.push_back(he_to_vert_idx[he_mesh.half_edges[he.pair].prev]);
        out_mesh.faces.push_back(he_to_vert_idx[he_idx]);
      }
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Bitruncate operator: Truncate the rectified mesh.
 */
template <typename MeshT>
FLASHMEM static PolyMesh bitruncate(const MeshT &mesh, Arena &target,
                                    Arena &temp, float t = 1.0f / 3.0f) {
  return truncate(ambo(mesh, target, temp), temp, target, t);
}

/**
 * @brief Chamfer operator: Replaces edges with hexagonal faces.
 * @param t Thickness factor for the new hexagons [0..1].
 */
template <typename MeshT>
FLASHMEM static PolyMesh chamfer(const MeshT &mesh, Arena &target, Arena &temp,
                                 float t = 0.5f) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, V + 2 * E);
  out_mesh.face_counts.bind(target, F + E);
  out_mesh.faces.bind(target, I + 6 * E);

  {
    ScratchScope _(temp);
    HalfEdgeMesh he_mesh(temp, mesh);

    uint16_t *he_to_new_v =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_new_v, I, HE_NONE);

    // 1. Copy original vertices
    for (size_t i = 0; i < V; ++i) {
      out_mesh.vertices.push_back(mesh.vertices[i]);
    }

    // 2. Generate new vertices and shrunk faces
    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      out_mesh.face_counts.push_back(count);
      he_idx = start;

      if (he_idx != HE_NONE) {
        do {
          uint16_t vi =
              he_mesh.half_edges[he_mesh.half_edges[he_idx].prev].vertex;
          Vector v = mesh.vertices[vi];
          Vector new_v = v + (centroid - v) * t;

          out_mesh.vertices.push_back(new_v);
          uint16_t idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_new_v[he_idx] = idx;

          out_mesh.faces.push_back(idx);

          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    // 3. Generate Hexagon faces for edges
    bool *visited_edges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visited_edges, I, false);

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      uint16_t he_idx = static_cast<uint16_t>(i);
      const HalfEdge &he = he_mesh.half_edges[he_idx];

      if (visited_edges[he_idx])
        continue;

      visited_edges[he_idx] = true;
      if (he.pair != HE_NONE)
        visited_edges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(6);

        uint16_t A = he_mesh.half_edges[he.prev].vertex;
        uint16_t B = he.vertex;

        out_mesh.faces.push_back(A);
        out_mesh.faces.push_back(he_to_new_v[he_mesh.half_edges[he.pair].next]);
        out_mesh.faces.push_back(he_to_new_v[he.pair]);
        out_mesh.faces.push_back(B);
        out_mesh.faces.push_back(he_to_new_v[he.next]);
        out_mesh.faces.push_back(he_to_new_v[he_idx]);
      }
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Edge-length relaxation by spring forces on the unit sphere.
 */
template <typename MeshT>
FLASHMEM static PolyMesh relax(const MeshT &mesh, Arena &target, Arena &temp,
                               int iterations = 8) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  const auto *face_counts = mesh.get_face_counts_data();
  const auto *faces = mesh.get_faces_data();

  out_mesh.vertices.bind(target, V);
  out_mesh.face_counts.bind(target, F);
  out_mesh.faces.bind(target, I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);
  for (size_t i = 0; i < F; ++i)
    out_mesh.face_counts.push_back(face_counts[i]);
  for (size_t i = 0; i < I; ++i)
    out_mesh.faces.push_back(faces[i]);

  {
    ScratchScope _(temp);

    ArenaVector<Vector> movements;
    movements.bind(temp, V);
    for (size_t i = 0; i < V; ++i)
      movements.push_back(Vector(0, 0, 0));

    HalfEdgeMesh he_mesh(temp, out_mesh);

    for (int iter = 0; iter < iterations; ++iter) {
      float total_len = 0;
      int edge_count = 0;

      for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
        const HalfEdge &he = he_mesh.half_edges[i];
        int u = he_mesh.half_edges[he.prev].vertex;
        int v = he.vertex;
        if (u < v) {
          total_len +=
              distance_between(out_mesh.vertices[u], out_mesh.vertices[v]);
          edge_count++;
        }
      }

      if (edge_count == 0)
        break;
      float target_len = static_cast<float>(total_len / edge_count);

      for (size_t i = 0; i < V; ++i) {
        Vector force(0, 0, 0);

        HEVertex &hev = he_mesh.vertices[i];
        uint16_t he_idx = hev.half_edge;
        if (he_idx != HE_NONE) {
          uint16_t start = he_idx;
          // A manifold vertex orbit visits at most every half-edge once.
          // Exceeding that means the twin graph is corrupt and this hand-rolled
          // walk would spin forever — trap so it's caught on the bench instead
          // of hanging the device (mirrors vertex_orbit's always-on guard).
          int orbit_count = 0;
          const int max_orbit = static_cast<int>(he_mesh.half_edges.size());
          do {
            HS_CHECK(orbit_count++ < max_orbit);
            const HalfEdge &curr_he = he_mesh.half_edges[he_idx];
            int ni = he_mesh.half_edges[curr_he.prev].vertex;
            Vector vec = out_mesh.vertices[ni] - out_mesh.vertices[i];
            // Compare SQUARED length against the squared tolerance (the
            // codebase idiom), so near-zero edges are skipped before 1/dist
            // can explode into a huge force spike. sqrt only runs for
            // non-degenerate edges.
            float len_sq = dot(vec, vec);
            if (len_sq > math::EPS_LEN_SQ) {
              float dist = sqrtf(len_sq);
              float diff = dist - target_len;
              force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
            }

            if (he_mesh.half_edges[curr_he.next].pair == HE_NONE)
              break;
            he_idx = he_mesh.half_edges[curr_he.next].pair;
          } while (he_idx != HE_NONE && he_idx != start);
        }
        movements[i] = force;
      }

      float max_move_sq = 0.0f;
      for (size_t i = 0; i < V; ++i) {
        max_move_sq = std::max(max_move_sq, dot(movements[i], movements[i]));
        out_mesh.vertices[i] =
            (out_mesh.vertices[i] + movements[i]).normalized();
      }

      // Converged: the largest per-vertex move this pass is below ~1e-4 rad
      // (sub-pixel on the unit sphere), so the remaining iterations would be
      // visual no-ops — stop early. The guard only fires once the spring system
      // has settled; a still-moving mesh runs the full iteration budget.
      constexpr float RELAX_CONVERGE_EPS_SQ = 1e-8f;
      if (max_move_sq < RELAX_CONVERGE_EPS_SQ)
        break;
    }
  }
  return out_mesh;
}

/**
 * @brief Snub operator: Creates a chiral semi-regular polyhedron.
 * Uses Newell's method for face normals — robust to non-planar faces
 * on the unit sphere and to collinear vertex triplets.
 */
template <typename MeshT>
FLASHMEM static PolyMesh snub(const MeshT &mesh, Arena &target, Arena &temp,
                              float t = 0.5f, float twist = 0.0f) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, I);
  out_mesh.face_counts.bind(target, F + V + 2 * E);
  out_mesh.faces.bind(target, 5 * I);

  {
    ScratchScope _(temp);

    HalfEdgeMesh he_mesh(temp, mesh);
    uint16_t *he_to_vert_idx =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_vert_idx, I, HE_NONE);

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        temp.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      // Newell's method face normal — robust for sphere-projected faces.
      Vector normal_raw = face_normal(he_mesh, mesh, fi);
      Vector normal(0, 0, 0);
      if (dot(normal_raw, normal_raw) > math::EPS_NORMAL_SQ) {
        normal = normal_raw.normalized();
      } else if (dot(centroid, centroid) > math::EPS_LEN_SQ) {
        normal = centroid.normalized();
      }

      out_mesh.face_counts.push_back(count);
      he_idx = start;
      if (he_idx != HE_NONE) {
        do {
          Vector v = mesh.vertices[he_mesh.half_edges[he_idx].vertex];
          Vector new_v = v + (centroid - v) * t;

          if (twist != 0.0f) {
            Vector local = new_v - centroid;
            Quaternion q = make_rotation(normal, twist);
            new_v = centroid + rotate(local, q);
          }

          out_mesh.vertices.push_back(new_v);
          int idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_vert_idx[he_idx] = idx;

          out_mesh.faces.push_back(idx);

          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    bool *visited_verts = static_cast<bool *>(
        temp.allocate(V * sizeof(bool), alignof(bool)));
    emit_vertex_orbit_faces<'N'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/true,
        [&](uint16_t idx) {
          return he_to_vert_idx[he_mesh.half_edges[idx].prev];
        });

    bool *visited_edges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visited_edges, I, false);

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      uint16_t he_idx = static_cast<uint16_t>(i);
      const HalfEdge &he = he_mesh.half_edges[he_idx];

      if (visited_edges[he_idx])
        continue;

      visited_edges[he_idx] = true;
      if (he.pair != HE_NONE)
        visited_edges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(he_to_vert_idx[he.prev]);
        out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
        out_mesh.faces.push_back(he_to_vert_idx[he_idx]);

        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
        out_mesh.faces.push_back(he_to_vert_idx[he_mesh.half_edges[he.pair].prev]);
        out_mesh.faces.push_back(he_to_vert_idx[he_idx]);
      }
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

template <typename MeshT>
FLASHMEM static PolyMesh gyro(const MeshT &mesh, Arena &target, Arena &temp) {
  return dual(snub(mesh, target, temp), temp, target);
}

// ---------------------------------------------------------------------------
// Compositional operators (Hart's notation)
//
// These are defined as compositions of primitive operators. Equivalences are
// from Hart's reference implementation. Memory-wise, each composition uses
// the standard ping-pong of (target, temp) arenas.
//   meta   m = kj = kis of ambo (j = a)
//   needle n = kd = kis of dual
//   zip    z = dk = dual of kis (truncated dual)
//   bevel  b = ta = truncate of ambo (rectify-then-truncate)
//
// Standalone:
//   propeller — Hart: each face surrounded by quadrilateral propeller blades
//               and a smaller central face. Implemented standalone.
// ---------------------------------------------------------------------------

/// meta = kj. Hart's `m` operator. (j = a; kj = kis-of-ambo.)
template <typename MeshT>
FLASHMEM static PolyMesh meta(const MeshT &mesh, Arena &target, Arena &temp) {
  return kis(ambo(mesh, target, temp), temp, target);
}

/// needle = kd. Kis of the dual.
template <typename MeshT>
FLASHMEM static PolyMesh needle(const MeshT &mesh, Arena &target, Arena &temp) {
  return kis(dual(mesh, target, temp), temp, target);
}

/// zip = dk. Dual of kis (truncated dual).
template <typename MeshT>
FLASHMEM static PolyMesh zip(const MeshT &mesh, Arena &target, Arena &temp) {
  return dual(kis(mesh, target, temp), temp, target);
}

/// bevel = ta. Truncate of ambo.
template <typename MeshT>
FLASHMEM static PolyMesh bevel(const MeshT &mesh, Arena &target, Arena &temp,
                               float t = 0.25f) {
  return truncate(ambo(mesh, target, temp), temp, target, t);
}

// Propeller (Hart's `p`) and whirl/loft are deliberately not implemented
// here — they are chiral standalone operators whose blade geometry/winding
// is subtle on the sphere, and a half-baked version is worse than none.
// Add them under a separate change with a dedicated test for blade winding.

/**
 * @brief (Re)builds the mesh KDTree into the supplied arena.
 *
 * The tree's nodes live in `arena`. Every caller passes a scratch/temp arena
 * that is reset between frames, so the tree CANNOT be cached across calls: a
 * persistent `cache_valid` flag would leave `mesh.kd_tree.nodes` dangling into
 * reset/overwritten scratch on the next call and return silent garbage
 * nearest-neighbours in release (the stale-binding guard that would catch this
 * is debug-only). We therefore rebuild every call, mirroring
 * closest_point_on_mesh_graph(), which already rebuilds the HalfEdgeMesh — the
 * dominant cost — each call, so the KDTree rebuild adds little.
 */
FLASHMEM static void compute_kdtree(const PolyMesh &mesh, Arena &arena) {
  mesh.kd_tree =
      KDTree(arena, std::span<const Vector>(mesh.vertices.data(),
                                            mesh.vertices.size()));
}

inline Vector closest_point_on_mesh_graph(const Vector &p, const PolyMesh &mesh,
                                          Arena &temp_arena) {
  if (mesh.vertices.empty())
    return Vector(0, 1, 0);

  compute_kdtree(mesh, temp_arena);

  auto nearest_nodes = mesh.kd_tree.nearest(p, 1);
  if (nearest_nodes.size() == 0)
    return mesh.vertices[0];

  const auto &node = nearest_nodes[0];
  int closest_vertex_index = (int)node.original_index;
  Vector closest_vertex_pos = node.point;

  Vector best_point = closest_vertex_pos;
  float max_dot = dot(p, best_point);

  HalfEdgeMesh he_mesh(temp_arena, mesh);
  if (closest_vertex_index < 0 ||
      closest_vertex_index >= (int)he_mesh.vertices.size())
    return best_point;

  HEVertex &hev = he_mesh.vertices[closest_vertex_index];
  uint16_t he_idx = hev.half_edge;
  if (he_idx == HE_NONE)
    return best_point;

  Vector A = closest_vertex_pos;
  uint16_t start = he_idx;

  do {
    const HalfEdge &curr_he = he_mesh.half_edges[he_idx];
    int neighbor_idx = he_mesh.half_edges[curr_he.prev].vertex;
    Vector B = mesh.vertices[neighbor_idx];

    Vector N = cross(A, B);
    float len_sq = dot(N, N);
    if (len_sq >= math::EPS_LEN_SQ) {
      N = N * (1.0f / sqrtf(len_sq));
      float p_dot_n = dot(p, N);
      Vector proj = p + N * (-p_dot_n);
      proj.normalize();

      Vector cross_ac = cross(A, proj);
      Vector cross_cb = cross(proj, B);

      if (dot(cross_ac, N) > 0 && dot(cross_cb, N) > 0) {
        float d = dot(p, proj);
        if (d > max_dot) {
          max_dot = d;
          best_point = proj;
        }
      }
    }

    if (curr_he.pair == HE_NONE)
      break;
    he_idx = he_mesh.half_edges[curr_he.pair].next;
  } while (he_idx != HE_NONE && he_idx != start);

  return best_point;
}

} // namespace MeshOps
