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
                      ScratchFront arena) {
  world_state.face_counts_view = ArenaSpan(local_state.face_counts);
  world_state.faces_view = ArenaSpan(local_state.faces);

  if constexpr (requires { world_state.face_offsets_view; }) {
    world_state.face_offsets_view = ArenaSpan(local_state.face_offsets);
  }

  world_state.vertices.initialize(arena, local_state.vertices.size());

  for (size_t i = 0; i < local_state.vertices.size(); ++i) {
    world_state.vertices.push_back(local_state.vertices[i]);
  }
}

template <typename MeshT, typename T1, typename... Transformers>
inline void transform(const MeshT &mesh, MeshT &transformed, ScratchFront arena,
                      const T1 &first_transformer,
                      const Transformers &...transformers) {
  transformed.face_counts_view = ArenaSpan(mesh.face_counts);
  transformed.faces_view = ArenaSpan(mesh.faces);

  if constexpr (requires { transformed.face_offsets_view; }) {
    transformed.face_offsets_view = ArenaSpan(mesh.face_offsets);
  }

  transformed.vertices.initialize(arena, mesh.vertices.size());

  for (size_t i = 0; i < mesh.vertices.size(); ++i) {
    Vector v = mesh.vertices[i];

    // Apply first transformer
    v = first_transformer(v);

    // Unroll remaining transformers at compile time using a fold expression
    (..., (v = transformers(v)));

    transformed.vertices.push_back(v);
  }
}

/**
 * @brief Computes the dual of a mesh.
 */
FLASHMEM inline PolyMesh dual(const PolyMesh &mesh, MemoryCtx &ctx) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(ctx.get_scratch_front(), F);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), V);
  out_mesh.faces.initialize(ctx.get_scratch_front(), I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    HalfEdgeMesh heMesh(ctx.get_scratch_back(), mesh);

    for (size_t i = 0; i < heMesh.faces.size(); ++i) {
      HEFace &face = heMesh.faces[i];
      Vector c(0, 0, 0);
      int count = 0;
      uint16_t heIdx = face.halfEdge;
      uint16_t start = heIdx;
      if (heIdx != HE_NONE) {
        do {
          c = c + mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start && count < 100);
      }

      if (count > 0)
        c = c / static_cast<float>(count);
      out_mesh.vertices.push_back(c.normalize());
    }

    bool *visitedVerts = static_cast<bool *>(
        ctx.get_scratch_back().allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int new_face_count = 0;
      uint16_t local_face[100];
      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (currHe.face == HE_NONE)
          break;
        if (new_face_count < 100)
          local_face[new_face_count++] =
              currHe.face; // currHe.face matches face index

        if (currHe.prev == HE_NONE ||
            heMesh.halfEdges[currHe.prev].pair == HE_NONE)
          break;
        currIdx = heMesh.halfEdges[currHe.prev].pair;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (new_face_count >= 3) {
        out_mesh.face_counts.push_back(new_face_count);
        for (int k = 0; k < new_face_count; ++k) {
          out_mesh.faces.push_back(local_face[k]);
        }
      }
    }
  }
  return out_mesh;
}

// --- CONWAY OPERATORS ---

/**
 * @brief Normalizes all vertices in the mesh to the unit sphere.
 */
template <typename MeshT> static void normalize(MeshT &mesh) {
  for (auto &v : mesh.vertices) {
    v = v.normalize();
  }
}

/**
 * @brief Kis operator: Raises a pyramid on each face.
 */
FLASHMEM inline PolyMesh kis(const PolyMesh &mesh, MemoryCtx &ctx) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(ctx.get_scratch_front(), V + F);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), I);
  out_mesh.faces.initialize(ctx.get_scratch_front(), 3 * I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    for (size_t i = 0; i < V; ++i)
      out_mesh.vertices.push_back(mesh.vertices[i]);

    size_t offset = 0;
    for (size_t fi = 0; fi < F; ++fi) {
      int count = mesh.face_counts[fi];
      Vector centroid(0, 0, 0);
      for (int i = 0; i < count; ++i) {
        centroid = centroid + mesh.vertices[mesh.faces[offset + i]];
      }
      if (count > 0)
        centroid = centroid / static_cast<float>(count);

      out_mesh.vertices.push_back(centroid);
      int centerIdx = static_cast<int>(out_mesh.vertices.size()) - 1;

      for (int i = 0; i < count; ++i) {
        uint16_t vi = mesh.faces[offset + i];
        uint16_t vj = mesh.faces[offset + (i + 1) % count];
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
FLASHMEM inline PolyMesh ambo(const PolyMesh &mesh, MemoryCtx &ctx) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(ctx.get_scratch_front(), E);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), F + V);
  out_mesh.faces.initialize(ctx.get_scratch_front(), 2 * I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    // 1. Build HE Mesh FIRST
    HalfEdgeMesh heMesh(ctx.get_scratch_back(), mesh);

    uint16_t *edgeToVert =
        static_cast<uint16_t *>(ctx.get_scratch_back().allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(edgeToVert, I, HE_NONE);

    // 3. Populate Vertices
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      if (edgeToVert[i] == HE_NONE) {
        HalfEdge &he = heMesh.halfEdges[i];
        if (he.prev == HE_NONE)
          continue;
        uint16_t v1 = he.vertex;
        uint16_t v2 = heMesh.halfEdges[he.prev].vertex;

        Vector mid = (mesh.vertices[v1] + mesh.vertices[v2]) * 0.5f;
        out_mesh.vertices.push_back(mid);
        uint16_t newIdx = static_cast<int16_t>(out_mesh.vertices.size() - 1);

        edgeToVert[i] = newIdx;
        if (he.pair != HE_NONE)
          edgeToVert[he.pair] = newIdx;
      }
    }

    // 4. Reconstruct Original Faces (Shrunk)
    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      HEFace &face = heMesh.faces[fi];
      uint16_t heIdx = face.halfEdge;
      int count = 0;
      if (heIdx != HE_NONE) {
        uint16_t start = heIdx;
        do {
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start && count < 100);

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
    bool *visitedVerts = static_cast<bool *>(
        ctx.get_scratch_back().allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int count = 0;
      uint16_t local_face[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (currHe.face == HE_NONE)
          break;

        if (count < 100)
          local_face[count++] = edgeToVert[currIdx];

        if (currHe.prev == HE_NONE ||
            heMesh.halfEdges[currHe.prev].pair == HE_NONE)
          break;
        currIdx = heMesh.halfEdges[currHe.prev].pair;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (count >= 3) {
        out_mesh.face_counts.push_back(count);
        for (int k = 0; k < count; ++k)
          out_mesh.faces.push_back(local_face[k]);
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
FLASHMEM inline PolyMesh truncate(const PolyMesh &mesh, MemoryCtx &ctx,
                                  float t = 0.25f) {
  if (std::abs(t - 0.5f) < 1e-4f) {
    return ambo(mesh, ctx);
  }

  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(ctx.get_scratch_front(), 2 * E);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), F + V);
  out_mesh.faces.initialize(ctx.get_scratch_front(), 3 * I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    HalfEdgeMesh heMesh(ctx.get_scratch_back(), mesh);

    std::pair<int16_t, int16_t> *edgeToVert =
        static_cast<std::pair<int16_t, int16_t> *>(
            ctx.get_scratch_back().allocate(
                I * sizeof(std::pair<int16_t, int16_t>),
                alignof(std::pair<int16_t, int16_t>)));
    std::fill_n(edgeToVert, I, std::make_pair<int16_t, int16_t>(-1, -1));

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      if (edgeToVert[i].first == -1) {
        HalfEdge &he = heMesh.halfEdges[i];
        if (he.prev == HE_NONE)
          continue;
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
      HEFace &face = heMesh.faces[fi];
      uint16_t heIdx = face.halfEdge;
      int count = 0;
      if (heIdx != HE_NONE) {
        uint16_t start = heIdx;
        do {
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start && count < 100);

        if (count >= 3) {
          out_mesh.face_counts.push_back(count * 2);
          heIdx = start;
          do {
            HalfEdge &he = heMesh.halfEdges[heIdx];
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

    bool *visitedVerts = static_cast<bool *>(
        ctx.get_scratch_back().allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int count = 0;
      uint16_t local_face[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (currHe.face == HE_NONE)
          break;

        uint16_t vi = heMesh.halfEdges[currHe.prev].vertex;
        uint16_t vj = currHe.vertex;
        uint16_t k1 = std::min(vi, vj);
        std::pair<int16_t, int16_t> newVerts = edgeToVert[currIdx];

        if (count < 100)
          local_face[count++] = (vi == k1) ? newVerts.first : newVerts.second;

        if (currHe.prev == HE_NONE ||
            heMesh.halfEdges[currHe.prev].pair == HE_NONE)
          break;
        currIdx = heMesh.halfEdges[currHe.prev].pair;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (count >= 3) {
        out_mesh.face_counts.push_back(count);
        for (int k = 0; k < count; ++k)
          out_mesh.faces.push_back(local_face[k]);
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
FLASHMEM inline PolyMesh expand(const PolyMesh &mesh, MemoryCtx &ctx,
                                float t = 2.0f - sqrt(2.0f)) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(ctx.get_scratch_front(), I);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), F + V + E);
  out_mesh.faces.initialize(ctx.get_scratch_front(), 4 * I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    HalfEdgeMesh heMesh(ctx.get_scratch_back(), mesh);
    int16_t *heToVertIdx = static_cast<int16_t *>(
        ctx.get_scratch_back().allocate(I * sizeof(int16_t), alignof(int16_t)));
    std::fill_n(heToVertIdx, I, -1);

    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      HEFace &face = heMesh.faces[fi];
      Vector centroid(0, 0, 0);
      int count = 0;
      uint16_t heIdx = face.halfEdge;
      uint16_t start = heIdx;
      if (heIdx != HE_NONE) {
        do {
          centroid = centroid + mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start && count < 100);
      }

      if (count > 0)
        centroid = centroid / static_cast<float>(count);

      out_mesh.face_counts.push_back(count);
      heIdx = start;
      int safety = 0;
      if (heIdx != HE_NONE) {
        do {
          Vector v = mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          Vector newV = v + (centroid - v) * t;
          out_mesh.vertices.push_back(newV);
          int idx = static_cast<int>(out_mesh.vertices.size()) - 1;
          heToVertIdx[heIdx] = idx;

          out_mesh.faces.push_back(idx);
          heIdx = heMesh.halfEdges[heIdx].next;
          safety++;
        } while (heIdx != HE_NONE && heIdx != start && safety < 100);
      }
    }

    bool *visitedVerts = static_cast<bool *>(
        ctx.get_scratch_back().allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int count = 0;
      uint16_t local_face[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (currHe.face == HE_NONE)
          break;
        if (count < 100)
          local_face[count++] = heToVertIdx[currHe.prev];

        if (currHe.pair == HE_NONE)
          break;
        currIdx = heMesh.halfEdges[currHe.pair].next;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (count >= 3) {
        out_mesh.face_counts.push_back(count);
        for (int k = count - 1; k >= 0; --k)
          out_mesh.faces.push_back(local_face[k]);
      }
    }

    bool *visitedEdges = static_cast<bool *>(
        ctx.get_scratch_back().allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visitedEdges, I, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heIdx = static_cast<uint16_t>(i);
      HalfEdge &he = heMesh.halfEdges[heIdx];

      if (he.prev == HE_NONE || visitedEdges[heIdx])
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
FLASHMEM inline PolyMesh bitruncate(const PolyMesh &mesh, MemoryCtx &ctx,
                                    float t = 1.0f / 3.0f) {
  return truncate(ambo(mesh, ctx), ctx, t);
}

FLASHMEM inline PolyMesh canonicalize(const PolyMesh &mesh, MemoryCtx &ctx,
                                      int iterations = 8) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(ctx.get_scratch_front(), V);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), F);
  out_mesh.faces.initialize(ctx.get_scratch_front(), I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);
  for (size_t i = 0; i < F; ++i)
    out_mesh.face_counts.push_back(mesh.face_counts[i]);
  for (size_t i = 0; i < I; ++i)
    out_mesh.faces.push_back(mesh.faces[i]);

  {
    ScopedScratch _(ctx.get_scratch_back());

    ArenaVector<Vector> movements;
    movements.initialize(ctx.get_scratch_back(), V);
    for (size_t i = 0; i < V; ++i)
      movements.push_back(Vector(0, 0, 0));

    HalfEdgeMesh heMesh(ctx.get_scratch_back(),
                        out_mesh); // Use out_mesh as base since we mutate it

    for (int iter = 0; iter < iterations; ++iter) {
      double totalLen = 0;
      int edgeCount = 0;

      for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
        HalfEdge &he = heMesh.halfEdges[i];
        if (he.prev == HE_NONE)
          continue;
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
          int safety = 0;
          do {
            HalfEdge &currHe = heMesh.halfEdges[heIdx];
            if (currHe.prev == HE_NONE)
              break;
            int ni = heMesh.halfEdges[currHe.prev].vertex;
            Vector vec = out_mesh.vertices[ni] - out_mesh.vertices[i];
            float dist = vec.length();
            if (dist > 1e-6f) {
              float diff = dist - targetLen;
              force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
            }

            if (currHe.next == HE_NONE ||
                heMesh.halfEdges[currHe.next].pair == HE_NONE)
              break;
            heIdx = heMesh.halfEdges[currHe.next].pair;
            safety++;
          } while (heIdx != HE_NONE && heIdx != start && safety < 100);
        }
        movements[i] = force;
      }

      for (size_t i = 0; i < V; ++i) {
        out_mesh.vertices[i] =
            (out_mesh.vertices[i] + movements[i]).normalize();
      }
    }
  }
  return out_mesh;
}

/**
 * @brief Snub operator: Creates a chiral semi-regular polyhedron.
 * Updated with twist support.
 */
FLASHMEM inline PolyMesh snub(const PolyMesh &mesh, MemoryCtx &ctx,
                              float t = 0.5f, float twist = 0.0f) {
  ctx.swap_scratch();
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(ctx.get_scratch_front(), I);
  out_mesh.face_counts.initialize(ctx.get_scratch_front(), F + V + 2 * E);
  out_mesh.faces.initialize(ctx.get_scratch_front(), 5 * I);

  {
    ScopedScratch _(ctx.get_scratch_back());

    HalfEdgeMesh heMesh(ctx.get_scratch_back(), mesh);
    uint16_t *heToVertIdx =
        static_cast<uint16_t *>(ctx.get_scratch_back().allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToVertIdx, I, HE_NONE);

    for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
      HEFace &face = heMesh.faces[fi];
      Vector centroid(0, 0, 0);
      int count = 0;
      uint16_t heIdx = face.halfEdge;
      uint16_t start = heIdx;
      if (heIdx != HE_NONE) {
        do {
          centroid = centroid + mesh.vertices[heMesh.halfEdges[heIdx].vertex];
          count++;
          heIdx = heMesh.halfEdges[heIdx].next;
        } while (heIdx != HE_NONE && heIdx != start && count < 100);
      }

      if (count > 0)
        centroid = centroid / static_cast<float>(count);

      Vector normal(0, 0, 0);
      if (count >= 3 && start != HE_NONE) {
        uint16_t nextIdx = heMesh.halfEdges[start].next;
        uint16_t nextNextIdx = heMesh.halfEdges[nextIdx].next;
        Vector ab = mesh.vertices[heMesh.halfEdges[nextIdx].vertex] -
                    mesh.vertices[heMesh.halfEdges[start].vertex];
        Vector ac = mesh.vertices[heMesh.halfEdges[nextNextIdx].vertex] -
                    mesh.vertices[heMesh.halfEdges[start].vertex];
        normal = cross(ab, ac).normalize();
      }
      if (dot(centroid, centroid) > 1e-6f && dot(normal, normal) < 1e-9f) {
        normal = centroid.normalize();
      }

      out_mesh.face_counts.push_back(count);
      heIdx = start;
      int safety = 0;
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
          safety++;
        } while (heIdx != HE_NONE && heIdx != start && safety < 100);
      }
    }

    bool *visitedVerts = static_cast<bool *>(
        ctx.get_scratch_back().allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;

      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int count = 0;
      uint16_t local_face[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (currHe.face == HE_NONE)
          break;
        if (count < 100)
          local_face[count++] = heToVertIdx[currHe.prev];

        if (currHe.pair == HE_NONE)
          break;
        currIdx = heMesh.halfEdges[currHe.pair].next;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (count >= 3) {
        out_mesh.face_counts.push_back(count);
        for (int k = count - 1; k >= 0; --k)
          out_mesh.faces.push_back(local_face[k]);
      }
    }

    bool *visitedEdges = static_cast<bool *>(
        ctx.get_scratch_back().allocate(I * sizeof(bool), alignof(bool)));
    std::fill_n(visitedEdges, I, false);

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heIdx = static_cast<uint16_t>(i);
      HalfEdge &he = heMesh.halfEdges[heIdx];

      if (he.prev == HE_NONE || visitedEdges[heIdx])
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

FLASHMEM inline PolyMesh gyro(const PolyMesh &mesh, MemoryCtx &ctx) {
  return dual(snub(mesh, ctx), ctx);
}

/**
 * @brief Computes KDTree and Adjacency map for the mesh (caching it).
 */
FLASHMEM inline void compute_kdtree(const PolyMesh &mesh, Arena &arena) {
  if (mesh.cache_valid)
    return;

  PolyMesh &m = const_cast<PolyMesh &>(mesh);
  m.kdTree =
      KDTree(arena, std::span<Vector>(m.vertices.data(), m.vertices.size()));
  m.cache_valid = true;
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
  int safety = 0;

  do {
    HalfEdge &currHe = heMesh.halfEdges[heIdx];
    if (currHe.prev == HE_NONE)
      break;
    int neighborIdx = heMesh.halfEdges[currHe.prev].vertex;
    Vector B = mesh.vertices[neighborIdx];

    Vector N = cross(A, B);
    float lenSq = dot(N, N);
    if (lenSq >= 1e-6f) {
      N = N * (1.0f / sqrt(lenSq));
      float pDotN = dot(p, N);
      Vector proj = p + N * (-pDotN);
      proj = proj.normalize();

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
    safety++;
  } while (heIdx != HE_NONE && heIdx != start && safety < 100);

  return bestPoint;
}

} // namespace MeshOps
