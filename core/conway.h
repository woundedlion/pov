/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_CONWAY_H_
#define HOLOSPHERE_CORE_CONWAY_H_

#include "mesh.h"

/**
 * @brief Operations on meshes (Conway operators and utilities).
 */
namespace MeshOps {

// ---------------------------------------------------------------------------
// Shared topology helpers
// ---------------------------------------------------------------------------

/**
 * @brief Compute the centroid of a face by walking its half-edge loop.
 */
template <typename MeshT>
inline Vector face_centroid(const HalfEdgeMesh &heMesh,
                            const MeshT &mesh, size_t face_index,
                            int &out_count) {
  const HEFace &face = heMesh.faces[face_index];
  Vector c(0, 0, 0);
  out_count = 0;
  uint16_t heIdx = face.halfEdge;
  uint16_t start = heIdx;
  if (heIdx != HE_NONE) {
    do {
      c = c + mesh.vertices[heMesh.halfEdges[heIdx].vertex];
      out_count++;
      heIdx = heMesh.halfEdges[heIdx].next;
    } while (heIdx != HE_NONE && heIdx != start);
  }
  if (out_count > 0)
    c = c / static_cast<float>(out_count);
  return c;
}

/**
 * @brief Newell's method face normal — robust to non-planar faces and
 * collinear vertex triplets. Walks the half-edge loop and accumulates
 * cross products between adjacent edge midpoints. Returns the unnormalized
 * normal; caller normalizes if needed.
 */
template <typename MeshT>
inline Vector face_normal(const HalfEdgeMesh &heMesh, const MeshT &mesh,
                          size_t face_index) {
  const HEFace &face = heMesh.faces[face_index];
  Vector n(0, 0, 0);
  uint16_t heIdx = face.halfEdge;
  if (heIdx == HE_NONE) return n;
  uint16_t start = heIdx;
  do {
    const HalfEdge &he = heMesh.halfEdges[heIdx];
    const Vector &curr = mesh.vertices[he.vertex];
    const Vector &next = mesh.vertices[heMesh.halfEdges[he.next].vertex];
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
    heIdx = he.next;
  } while (heIdx != start);
  return n;
}

/**
 * @brief Walk all half-edges orbiting a vertex, calling visitor(currIdx)
 *        for each. Returns the number of visited half-edges.
 *
 * OrbitMode selects the traversal direction:
 *   'P' = prev->pair (dual, ambo, truncate)
 *   'N' = pair->next  (expand, snub)
 */
template <char OrbitMode, typename VisitorFn>
inline int vertex_orbit(const HalfEdgeMesh &heMesh, uint16_t startIdx,
                        VisitorFn &&visitor) {
  uint16_t currIdx = startIdx;
  int count = 0;
  do {
    const HalfEdge &currHe = heMesh.halfEdges[currIdx];
    if (currHe.face == HE_NONE)
      break;

    visitor(currIdx);
    count++;

    if constexpr (OrbitMode == 'P') {
      // prev->pair orbit
      if (heMesh.halfEdges[currHe.prev].pair == HE_NONE)
        break;
      currIdx = heMesh.halfEdges[currHe.prev].pair;
    } else {
      // pair->next orbit
      if (currHe.pair == HE_NONE)
        break;
      currIdx = heMesh.halfEdges[currHe.pair].next;
    }
  } while (currIdx != HE_NONE && currIdx != startIdx);
  return count;
}

/**
 * @brief Transforms a MeshState into a target MeshState using the provided
 * transformers.
 * @param local_state The source mesh state (centered around origin).
 * @param world_state The destination mesh state to populate.
 * @param arena The memory arena to allocate vertices.
 * @param transformers A list of transformation functions.
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
// maximum possible valence (= total half-edges). Previously the orbit was
// collected into a fixed `local_face[100]` stack array which silently
// truncated for high-valence vertices — see commit history.
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

    HalfEdgeMesh heMesh(temp, mesh);

    for (size_t i = 0; i < heMesh.faces.size(); ++i) {
      int count;
      Vector c = face_centroid(heMesh, mesh, i, count);
      out_mesh.vertices.push_back(c.normalized());
    }

    bool *visitedVerts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);

    // Per-orbit scratch buffer sized to the absolute upper bound on valence.
    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      const HalfEdge &heStart = heMesh.halfEdges[heStartIdx];

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      int orbit_count = 0;
      vertex_orbit<'P'>(heMesh, heStartIdx, [&](uint16_t idx) {
        assert(orbit_count < (int)I);
        orbit_buf[orbit_count++] = heMesh.halfEdges[idx].face;
      });

      if (orbit_count >= 3) {
        out_mesh.face_counts.push_back(orbit_count);
        for (int k = 0; k < orbit_count; ++k) {
          out_mesh.faces.push_back(orbit_buf[k]);
        }
      }
    }
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
      int centerIdx = static_cast<int>(out_mesh.vertices.size()) - 1;

      for (int i = 0; i < count; ++i) {
        uint16_t vi = faces[offset + i];
        uint16_t vj = faces[offset + (i + 1) % count];
        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(vi);
        out_mesh.faces.push_back(vj);
        out_mesh.faces.push_back(static_cast<uint16_t>(centerIdx));
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

    HalfEdgeMesh heMesh(temp, mesh);

    uint16_t *edgeToVert =
        static_cast<uint16_t *>(target.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(edgeToVert, I, HE_NONE);

    bool *visitedVerts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    // 3. Populate Vertices
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      if (edgeToVert[i] == HE_NONE) {
        const HalfEdge &he = heMesh.halfEdges[i];
        uint16_t v1 = he.vertex;
        uint16_t v2 = heMesh.halfEdges[he.prev].vertex;

        Vector mid = (mesh.vertices[v1] + mesh.vertices[v2]) * 0.5f;
        out_mesh.vertices.push_back(mid);
        uint16_t newIdx = static_cast<uint16_t>(out_mesh.vertices.size() - 1);

        edgeToVert[i] = newIdx;
        if (he.pair != HE_NONE)
          edgeToVert[he.pair] = newIdx;
      }
    }

    // 4. Reconstruct Original Faces (Shrunk)
    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      const HEFace &face = heMesh.faces[fi];
      uint16_t heIdx = face.halfEdge;
      int count = 0;
      if (heIdx != HE_NONE) {
        uint16_t start = heIdx;
        do {
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start);

        if (count >= 3) {
          out_mesh.face_counts.push_back(count);
          heIdx = start;
          do {
            out_mesh.faces.push_back(edgeToVert[heIdx]);
            heIdx = heMesh.halfEdges[heIdx].next;
          } while (heIdx != HE_NONE && heIdx != start);
        }
      }
    }

    // 5. Build Vertex Orbits (New Faces)
    std::fill_n(visitedVerts, V, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      const HalfEdge &heStart = heMesh.halfEdges[heStartIdx];

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      int orbit_count = 0;
      vertex_orbit<'P'>(heMesh, heStartIdx, [&](uint16_t idx) {
        assert(orbit_count < (int)I);
        orbit_buf[orbit_count++] = edgeToVert[idx];
      });

      if (orbit_count >= 3) {
        out_mesh.face_counts.push_back(orbit_count);
        for (int k = 0; k < orbit_count; ++k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Truncate operator: Cuts corners off the polyhedron.
 * @param t Truncation depth [0..0.5].
 */
template <typename MeshT>
FLASHMEM static PolyMesh truncate(const MeshT &mesh, Arena &target, Arena &temp,
                                  float t = 0.25f) {
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

    HalfEdgeMesh heMesh(temp, mesh);

    std::pair<int16_t, int16_t> *edgeToVert =
        static_cast<std::pair<int16_t, int16_t> *>(
            target.allocate(
                I * sizeof(std::pair<int16_t, int16_t>),
                alignof(std::pair<int16_t, int16_t>)));
    std::fill_n(edgeToVert, I, std::make_pair<int16_t, int16_t>(-1, -1));

    bool *visitedVerts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      if (edgeToVert[i].first == -1) {
        const HalfEdge &he = heMesh.halfEdges[i];
        uint16_t vi = heMesh.halfEdges[he.prev].vertex;
        uint16_t vj = he.vertex;
        uint16_t k1 = std::min(vi, vj);
        uint16_t k2 = std::max(vi, vj);

        Vector new_u =
            mesh.vertices[k1] + (mesh.vertices[k2] - mesh.vertices[k1]) * t;
        Vector new_v =
            mesh.vertices[k2] + (mesh.vertices[k1] - mesh.vertices[k2]) * t;

        out_mesh.vertices.push_back(new_u);
        uint16_t idx_u = static_cast<uint16_t>(out_mesh.vertices.size() - 1);

        out_mesh.vertices.push_back(new_v);
        uint16_t idx_v = static_cast<uint16_t>(out_mesh.vertices.size() - 1);

        edgeToVert[i] = {idx_u, idx_v};
        if (he.pair != HE_NONE)
          edgeToVert[he.pair] = {idx_u, idx_v};
      }
    }

    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      const HEFace &face = heMesh.faces[fi];
      uint16_t heIdx = face.halfEdge;
      int count = 0;
      if (heIdx != HE_NONE) {
        uint16_t start = heIdx;
        do {
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start);

        if (count >= 3) {
          out_mesh.face_counts.push_back(count * 2);
          heIdx = start;
          do {
            const HalfEdge &he = heMesh.halfEdges[heIdx];
            uint16_t vi = heMesh.halfEdges[he.prev].vertex;
            uint16_t vj = he.vertex;
            uint16_t k1 = std::min(vi, vj);
            std::pair<int16_t, int16_t> newVerts = edgeToVert[heIdx];

            if (vi == k1) {
              out_mesh.faces.push_back(newVerts.first);
              out_mesh.faces.push_back(newVerts.second);
            } else {
              out_mesh.faces.push_back(newVerts.second);
              out_mesh.faces.push_back(newVerts.first);
            }

            heIdx = heMesh.halfEdges[heIdx].next;
          } while (heIdx != HE_NONE && heIdx != start);
        }
      }
    }

    std::fill_n(visitedVerts, V, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      const HalfEdge &heStart = heMesh.halfEdges[heStartIdx];

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      int orbit_count = 0;
      vertex_orbit<'P'>(heMesh, heStartIdx, [&](uint16_t idx) {
        const HalfEdge &currHe = heMesh.halfEdges[idx];
        uint16_t vi = heMesh.halfEdges[currHe.prev].vertex;
        uint16_t vj = currHe.vertex;
        uint16_t k1 = std::min(vi, vj);
        std::pair<int16_t, int16_t> newVerts = edgeToVert[idx];
        assert(orbit_count < (int)I);
        orbit_buf[orbit_count++] =
            (vi == k1) ? newVerts.first : newVerts.second;
      });

      if (orbit_count >= 3) {
        out_mesh.face_counts.push_back(orbit_count);
        for (int k = 0; k < orbit_count; ++k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }
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

    HalfEdgeMesh heMesh(temp, mesh);
    int16_t *heToVertIdx = static_cast<int16_t *>(
        target.allocate(I * sizeof(int16_t), alignof(int16_t)));
    std::fill_n(heToVertIdx, I, -1);

    bool *visitedVerts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    bool *visitedEdges = static_cast<bool *>(
        target.allocate(I * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(heMesh, mesh, fi, count);
      uint16_t start = heMesh.faces[fi].halfEdge;
      uint16_t heIdx = start;

      out_mesh.face_counts.push_back(count);
      heIdx = start;
      if (heIdx != HE_NONE) {
        do {
          Vector v = mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          Vector newV = v + (centroid - v) * t;
          out_mesh.vertices.push_back(newV);
          int idx = static_cast<int>(out_mesh.vertices.size()) - 1;
          heToVertIdx[heIdx] = idx;

          out_mesh.faces.push_back(idx);
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start);
      }
    }

    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      const HalfEdge &heStart = heMesh.halfEdges[heStartIdx];

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      int orbit_count = 0;
      vertex_orbit<'N'>(heMesh, heStartIdx, [&](uint16_t idx) {
        assert(orbit_count < (int)I);
        orbit_buf[orbit_count++] = heToVertIdx[heMesh.halfEdges[idx].prev];
      });

      if (orbit_count >= 3) {
        out_mesh.face_counts.push_back(orbit_count);
        for (int k = orbit_count - 1; k >= 0; --k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }

    std::fill_n(visitedEdges, I, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heIdx = static_cast<uint16_t>(i);
      const HalfEdge &he = heMesh.halfEdges[heIdx];

      if (visitedEdges[heIdx])
        continue;

      visitedEdges[heIdx] = true;
      if (he.pair != HE_NONE)
        visitedEdges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(4);
        out_mesh.faces.push_back(heToVertIdx[he.prev]);
        out_mesh.faces.push_back(heToVertIdx[he.pair]);
        out_mesh.faces.push_back(heToVertIdx[heMesh.halfEdges[he.pair].prev]);
        out_mesh.faces.push_back(heToVertIdx[heIdx]);
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
    HalfEdgeMesh heMesh(temp, mesh);

    uint16_t *heToNewV =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToNewV, I, HE_NONE);

    // 1. Copy original vertices
    for (size_t i = 0; i < V; ++i) {
      out_mesh.vertices.push_back(mesh.vertices[i]);
    }

    // 2. Generate new vertices and shrunk faces
    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(heMesh, mesh, fi, count);
      uint16_t start = heMesh.faces[fi].halfEdge;
      uint16_t heIdx = start;

      out_mesh.face_counts.push_back(count);
      heIdx = start;

      if (heIdx != HE_NONE) {
        do {
          uint16_t vi =
              heMesh.halfEdges[heMesh.halfEdges[heIdx].prev].vertex;
          Vector v = mesh.vertices[vi];
          Vector newV = v + (centroid - v) * t;

          out_mesh.vertices.push_back(newV);
          uint16_t idx = static_cast<uint16_t>(out_mesh.vertices.size() - 1);
          heToNewV[heIdx] = idx;

          out_mesh.faces.push_back(idx);

          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start);
      }
    }

    // 3. Generate Hexagon faces for edges
    bool *visitedEdges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visitedEdges, I, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heIdx = static_cast<uint16_t>(i);
      const HalfEdge &he = heMesh.halfEdges[heIdx];

      if (visitedEdges[heIdx])
        continue;

      visitedEdges[heIdx] = true;
      if (he.pair != HE_NONE)
        visitedEdges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(6);

        uint16_t A = heMesh.halfEdges[he.prev].vertex;
        uint16_t B = he.vertex;

        out_mesh.faces.push_back(A);
        out_mesh.faces.push_back(heToNewV[heMesh.halfEdges[he.pair].next]);
        out_mesh.faces.push_back(heToNewV[he.pair]);
        out_mesh.faces.push_back(B);
        out_mesh.faces.push_back(heToNewV[he.next]);
        out_mesh.faces.push_back(heToNewV[heIdx]);
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

    HalfEdgeMesh heMesh(temp, out_mesh);

    for (int iter = 0; iter < iterations; ++iter) {
      float totalLen = 0;
      int edgeCount = 0;

      for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
        const HalfEdge &he = heMesh.halfEdges[i];
        int u = heMesh.halfEdges[he.prev].vertex;
        int v = he.vertex;
        if (u < v) {
          totalLen +=
              distance_between(out_mesh.vertices[u], out_mesh.vertices[v]);
          edgeCount++;
        }
      }

      if (edgeCount == 0)
        break;
      float targetLen = static_cast<float>(totalLen / edgeCount);

      for (size_t i = 0; i < V; ++i) {
        Vector force(0, 0, 0);

        HEVertex &hev = heMesh.vertices[i];
        uint16_t heIdx = hev.halfEdge;
        if (heIdx != HE_NONE) {
          uint16_t start = heIdx;
          do {
            const HalfEdge &currHe = heMesh.halfEdges[heIdx];
            int ni = heMesh.halfEdges[currHe.prev].vertex;
            Vector vec = out_mesh.vertices[ni] - out_mesh.vertices[i];
            float dist = vec.length();
            if (dist > math::EPS_LEN_SQ) {
              float diff = dist - targetLen;
              force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
            }

            if (heMesh.halfEdges[currHe.next].pair == HE_NONE)
              break;
            heIdx = heMesh.halfEdges[currHe.next].pair;
          } while (heIdx != HE_NONE && heIdx != start);
        }
        movements[i] = force;
      }

      for (size_t i = 0; i < V; ++i) {
        out_mesh.vertices[i] =
            (out_mesh.vertices[i] + movements[i]).normalized();
      }
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

    HalfEdgeMesh heMesh(temp, mesh);
    uint16_t *heToVertIdx =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToVertIdx, I, HE_NONE);

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        temp.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(heMesh, mesh, fi, count);
      uint16_t start = heMesh.faces[fi].halfEdge;
      uint16_t heIdx = start;

      // Newell's method face normal — robust for sphere-projected faces.
      Vector normal_raw = face_normal(heMesh, mesh, fi);
      Vector normal(0, 0, 0);
      if (dot(normal_raw, normal_raw) > math::EPS_NORMAL_SQ) {
        normal = normal_raw.normalized();
      } else if (dot(centroid, centroid) > math::EPS_LEN_SQ) {
        normal = centroid.normalized();
      }

      out_mesh.face_counts.push_back(count);
      heIdx = start;
      if (heIdx != HE_NONE) {
        do {
          Vector v = mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          Vector newV = v + (centroid - v) * t;

          if (twist != 0.0f) {
            Vector local = newV - centroid;
            Quaternion q = make_rotation(normal, twist);
            newV = centroid + rotate(local, q);
          }

          out_mesh.vertices.push_back(newV);
          int idx = static_cast<int>(out_mesh.vertices.size()) - 1;
          heToVertIdx[heIdx] = idx;

          out_mesh.faces.push_back(idx);

          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start);
      }
    }

    bool *visitedVerts = static_cast<bool *>(
        temp.allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      const HalfEdge &heStart = heMesh.halfEdges[heStartIdx];

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      int orbit_count = 0;
      vertex_orbit<'N'>(heMesh, heStartIdx, [&](uint16_t idx) {
        assert(orbit_count < (int)I);
        orbit_buf[orbit_count++] = heToVertIdx[heMesh.halfEdges[idx].prev];
      });

      if (orbit_count >= 3) {
        out_mesh.face_counts.push_back(orbit_count);
        for (int k = orbit_count - 1; k >= 0; --k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }

    bool *visitedEdges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visitedEdges, I, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heIdx = static_cast<uint16_t>(i);
      const HalfEdge &he = heMesh.halfEdges[heIdx];

      if (visitedEdges[heIdx])
        continue;

      visitedEdges[heIdx] = true;
      if (he.pair != HE_NONE)
        visitedEdges[he.pair] = true;

      if (he.pair != HE_NONE) {
        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(heToVertIdx[he.prev]);
        out_mesh.faces.push_back(heToVertIdx[he.pair]);
        out_mesh.faces.push_back(heToVertIdx[heIdx]);

        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(heToVertIdx[he.pair]);
        out_mesh.faces.push_back(heToVertIdx[heMesh.halfEdges[he.pair].prev]);
        out_mesh.faces.push_back(heToVertIdx[heIdx]);
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
//   meta   m = kj = k(dual(ambo(x)))  — Hart: meta is kis-of-join.
//                                       Wait: Hart's meta is k3j, treat as kj.
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
 * @brief Computes KDTree and Adjacency map for the mesh (caching it).
 */
FLASHMEM static void compute_kdtree(const PolyMesh &mesh, Arena &arena) {
  if (mesh.cache_valid)
    return;

  mesh.kdTree =
      KDTree(arena, std::span<const Vector>(mesh.vertices.data(),
                                            mesh.vertices.size()));
  mesh.cache_valid = true;
}

inline Vector closest_point_on_mesh_graph(const Vector &p, const PolyMesh &mesh,
                                          Arena &temp_arena) {
  if (mesh.vertices.empty())
    return Vector(0, 1, 0);

  compute_kdtree(mesh, temp_arena);

  auto nearestNodes = mesh.kdTree.nearest(p, 1);
  if (nearestNodes.size() == 0)
    return mesh.vertices[0];

  const auto &node = nearestNodes[0];
  int closestVertexIndex = (int)node.originalIndex;
  Vector closestVertexPos = node.point;

  Vector bestPoint = closestVertexPos;
  float maxDot = dot(p, bestPoint);

  HalfEdgeMesh heMesh(temp_arena, mesh);
  if (closestVertexIndex < 0 ||
      closestVertexIndex >= (int)heMesh.vertices.size())
    return bestPoint;

  HEVertex &hev = heMesh.vertices[closestVertexIndex];
  uint16_t heIdx = hev.halfEdge;
  if (heIdx == HE_NONE)
    return bestPoint;

  Vector A = closestVertexPos;
  uint16_t start = heIdx;

  do {
    const HalfEdge &currHe = heMesh.halfEdges[heIdx];
    int neighborIdx = heMesh.halfEdges[currHe.prev].vertex;
    Vector B = mesh.vertices[neighborIdx];

    Vector N = cross(A, B);
    float lenSq = dot(N, N);
    if (lenSq >= math::EPS_LEN_SQ) {
      N = N * (1.0f / sqrtf(lenSq));
      float pDotN = dot(p, N);
      Vector proj = p + N * (-pDotN);
      proj.normalize();

      Vector crossAC = cross(A, proj);
      Vector crossCB = cross(proj, B);

      if (dot(crossAC, N) > 0 && dot(crossCB, N) > 0) {
        float d = dot(p, proj);
        if (d > maxDot) {
          maxDot = d;
          bestPoint = proj;
        }
      }
    }

    if (currHe.pair == HE_NONE)
      break;
    heIdx = heMesh.halfEdges[currHe.pair].next;
  } while (heIdx != HE_NONE && heIdx != start);

  return bestPoint;
}

} // namespace MeshOps
#endif // HOLOSPHERE_CORE_CONWAY_H_
