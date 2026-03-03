/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "concepts.h"
#include "3dmath.h"
#include "spatial.h"
#include "memory.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <span>
#include <cmath>
#include "concepts.h"

// Forward declarations
struct HEVertex;
struct HEFace;
struct HalfEdge;

/**
 * @brief Forward declaration of MeshState for HalfEdgeMesh.
 */
struct MeshState;

/**
 * @brief A simple dynamic mesh structure compatible with MeshOps templates.
 */
struct PolyMesh {
  ArenaVector<Vector> vertices;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<uint16_t> faces;

  // Cache
  mutable bool cache_valid = false;
  mutable KDTree kdTree;

  PolyMesh() = default;

  void initialize(Arena &arena, size_t num_verts, size_t num_faces,
                  size_t num_indices) {
    vertices.initialize(arena, num_verts);
    face_counts.initialize(arena, num_faces);
    faces.initialize(arena, num_indices);
    cache_valid = false;
  }

  inline void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    cache_valid = false;
  }

  void clear_cache() const {
    cache_valid = false;
    kdTree.clear();
  }
};

constexpr uint16_t HE_NONE = 0xFFFF;

struct HalfEdge {
  uint16_t vertex = HE_NONE; /**< Vertex at the end of this half-edge. */
  uint16_t face = HE_NONE;   /**< Face this half-edge belongs to. */
  uint16_t next = HE_NONE;   /**< Next half-edge in the face loop. */
  uint16_t prev = HE_NONE;   /**< Previous half-edge in the face loop. */
  uint16_t pair = HE_NONE;   /**< Opposite half-edge. */
};

struct HEVertex {
  uint16_t halfEdge =
      HE_NONE; /**< One of the half-edges pointing to this vertex. */
};

struct HEFace {
  uint16_t halfEdge =
      HE_NONE; /**< One of the half-edges bordering this face. */
};

class HalfEdgeMesh {
public:
  ArenaVector<HEVertex> vertices;
  ArenaVector<HEFace> faces;
  ArenaVector<HalfEdge> halfEdges;

  explicit HalfEdgeMesh(Arena &arena, const PolyMesh &mesh) {
    build_from_flat(arena, mesh.vertices, mesh.face_counts, mesh.faces);
  }

  explicit HalfEdgeMesh(Arena &arena, const MeshState &mesh) {
    build_from_flat(arena, mesh.vertices, mesh.face_counts, mesh.faces);
  }

private:
  template <typename Verts, typename Counts, typename Faces>
  void build_from_flat(Arena &arena, const Verts &verts, const Counts &counts,
                       const Faces &faces_arr) {
    size_t num_verts = verts.size();
    size_t num_faces = counts.size();
    size_t total_indices = faces_arr.size();

    vertices.initialize(arena, num_verts);
    for (size_t i = 0; i < num_verts; ++i) {
      vertices.push_back({HE_NONE});
    }

    faces.initialize(arena, num_faces);
    halfEdges.initialize(arena, total_indices);

    struct EdgeRecord {
      uint16_t min_v;
      uint16_t max_v;
      uint16_t he;
    };

    size_t face_offset = 0;
    size_t he_idx = 0;

    {
      ArenaMarker _(arena);
      EdgeRecord *records = static_cast<EdgeRecord *>(arena.allocate(
          total_indices * sizeof(EdgeRecord), alignof(EdgeRecord)));

      for (size_t fi = 0; fi < num_faces; ++fi) {
        int count = counts[fi];

        faces.emplace_back();
        uint16_t currentFaceIdx = static_cast<uint16_t>(faces.size() - 1);
        size_t faceStartHeIdx = he_idx;

        for (int i = 0; i < count; ++i) {
          halfEdges.emplace_back();
        }

        for (int i = 0; i < count; ++i) {
          uint16_t u = faces_arr[face_offset + i];
          uint16_t v = faces_arr[face_offset + (i + 1) % count];

          uint16_t heIdx = static_cast<uint16_t>(faceStartHeIdx + i);
          HalfEdge &he = halfEdges[heIdx];
          he.vertex = v;
          he.face = currentFaceIdx;
          he.next = static_cast<uint16_t>(faceStartHeIdx + (i + 1) % count);
          he.prev =
              static_cast<uint16_t>(faceStartHeIdx + (i - 1 + count) % count);

          vertices[v].halfEdge = heIdx;

          records[he_idx].min_v = std::min(u, v);
          records[he_idx].max_v = std::max(u, v);
          records[he_idx].he = heIdx;

          he_idx++;
        }
        faces[currentFaceIdx].halfEdge = static_cast<uint16_t>(faceStartHeIdx);
        face_offset += count;
      }

      std::sort(records, records + total_indices,
                [](const EdgeRecord &a, const EdgeRecord &b) {
                  if (a.min_v != b.min_v)
                    return a.min_v < b.min_v;
                  return a.max_v < b.max_v;
                });

      for (size_t i = 0; i < total_indices;) {
        if (i + 1 < total_indices && records[i].min_v == records[i + 1].min_v &&
            records[i].max_v == records[i + 1].max_v) {
          halfEdges[records[i].he].pair = records[i + 1].he;
          halfEdges[records[i + 1].he].pair = records[i].he;
          i += 2;
        } else {
          i += 1;
        }
      }
    }
  }
};

/**
 * @brief Structure returned by compile_hankin.
 */
struct HankinInstruction {
  uint16_t vCorner; /**< Index to baseVertices for corner vertex. */
  uint16_t vPrev;   /**< Index to baseVertices for previous vertex. */
  uint16_t vNext;   /**< Index to baseVertices for next vertex. */
  uint16_t idxM1;   /**< Index of first midpoint (static vertex). */
  uint16_t idxM2;   /**< Index of second midpoint (static vertex). */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 */
struct CompiledHankin {
  ArenaVector<Vector> baseVertices;
  ArenaVector<Vector> staticVertices;
  ArenaVector<Vector> dynamicVertices;
  ArenaVector<HankinInstruction> dynamicInstructions;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<uint16_t> faces;
  int staticOffset;
};

/**
 * @brief Operations on meshes (Dual, Hankin, etc.).
 */
namespace MeshOps {

/**
 * @brief Compiles a PolyMesh into a static MeshState.
 * Removes degenerate faces (faces with < 3 vertices) during
 * the process.
 * @param src The source PolyMesh.
 * @return The compiled MeshState.
 */
inline void compile(const PolyMesh &src, MeshState &dst, Arena &geom_arena) {
  dst.clear();

  size_t valid_faces = 0;
  size_t valid_indices = 0;

  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    if (src.face_counts[i] >= 3) {
      valid_faces++;
      valid_indices += src.face_counts[i];
    }
  }

  dst.vertices.initialize(geom_arena, src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    dst.vertices.push_back(src.vertices[i]);
  }

  dst.face_counts.initialize(geom_arena, valid_faces);
  dst.faces.initialize(geom_arena, valid_indices);
  dst.face_offsets.initialize(geom_arena, valid_faces);

  size_t offset = 0;
  int current_offset = 0;
  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    int count = src.face_counts[i];
    if (count >= 3) {
      dst.face_counts.push_back(static_cast<uint8_t>(count));
      dst.face_offsets.push_back(static_cast<uint16_t>(current_offset));
      for (int k = 0; k < count; ++k) {
        dst.faces.push_back(src.faces[offset + k]);
      }
      current_offset += count;
    }
    offset += count;
  }
}

/**
 * @brief Performs a strict deep copy of a MeshState into a target arena.
 * Safe for memory compaction and bouncing between arenas.
 */
template <typename MeshT>
inline void clone(const MeshT &src, MeshT &dst, Arena &arena) {
  dst.vertices.initialize(arena, src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    dst.vertices.push_back(src.vertices[i]);
  }

  dst.face_counts.initialize(arena, src.face_counts.size());
  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    dst.face_counts.push_back(src.face_counts[i]);
  }

  dst.faces.initialize(arena, src.faces.size());
  for (size_t i = 0; i < src.faces.size(); ++i) {
    dst.faces.push_back(src.faces[i]);
  }

  if constexpr (requires { dst.face_offsets; }) {
    dst.face_offsets.initialize(arena, src.face_offsets.size());
    for (size_t i = 0; i < src.face_offsets.size(); ++i) {
      dst.face_offsets.push_back(src.face_offsets[i]);
    }
  }
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
                      Arena &arena) {
  world_state.face_counts.shallow_copy(local_state.face_counts);
  world_state.faces.shallow_copy(local_state.faces);

  if constexpr (requires { world_state.face_offsets; }) {
    world_state.face_offsets.shallow_copy(local_state.face_offsets);
  }

  world_state.vertices.initialize(arena, local_state.vertices.size());

  for (size_t i = 0; i < local_state.vertices.size(); ++i) {
    world_state.vertices.push_back(local_state.vertices[i]);
  }
}

template <typename MeshT, typename T1, typename... Transformers>
inline void transform(const MeshT &mesh, MeshT &transformed, Arena &arena,
                      const T1 &first_transformer,
                      const Transformers &...transformers) {
  transformed.face_counts.shallow_copy(mesh.face_counts);
  transformed.faces.shallow_copy(mesh.faces);

  if constexpr (requires { transformed.face_offsets; }) {
    transformed.face_offsets.shallow_copy(mesh.face_offsets);
  }

  transformed.vertices.initialize(arena, mesh.vertices.size());

  if constexpr (sizeof...(Transformers) > 0) {
    // Convert variadic args to TransformFn to avoid template bloat on transform
    // implementations
    std::array<TransformFn, sizeof...(Transformers)> tf_array = {
        transformers...};

    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
      Vector v = mesh.vertices[i];
      for (const auto &tf : tf_array) {
        if (tf) {
          v = tf(v);
        }
      }
      transformed.vertices.push_back(v);
    }
  } else {
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
      transformed.vertices.push_back(mesh.vertices[i]);
    }
  }
}

/**
 * @brief Computes the dual of a mesh.
 */
inline PolyMesh dual(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(*(ctx.target), F);
  out_mesh.face_counts.initialize(*(ctx.target), V);
  out_mesh.faces.initialize(*(ctx.target), I);

  {
    ArenaMarker _(*(ctx.source));

    HalfEdgeMesh heMesh(*(ctx.source), mesh);

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
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Compiles the topology for a Hankin pattern.
 */
inline void compile_hankin(const PolyMesh &mesh, CompiledHankin &compiled,
                           ScratchContext &ctx) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  compiled.baseVertices.initialize(*(ctx.source), V);
  for (size_t i = 0; i < V; ++i) {
    compiled.baseVertices.push_back(mesh.vertices[i]);
  }
  compiled.staticVertices.initialize(*(ctx.source), I);
  compiled.dynamicVertices.initialize(*(ctx.source), I);
  compiled.dynamicInstructions.initialize(*(ctx.source), I);
  compiled.face_counts.initialize(*(ctx.source), F + V);
  compiled.faces.initialize(*(ctx.source), 4 * I);

  {
    ArenaMarker _(*(ctx.target));

    HalfEdgeMesh heMesh(*(ctx.source), mesh);
    uint16_t *heToMidpointIdx = static_cast<uint16_t *>(
        ctx.source->allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToMidpointIdx, I, HE_NONE);
    uint16_t *heToDynamicIdx = static_cast<uint16_t *>(
        ctx.source->allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToDynamicIdx, I, HE_NONE);

    auto getMidpointIdx = [&](uint16_t heIdx) {
      if (heToMidpointIdx[heIdx] != HE_NONE)
        return heToMidpointIdx[heIdx];
      HalfEdge &he = heMesh.halfEdges[heIdx];
      if (he.pair != HE_NONE && heToMidpointIdx[he.pair] != HE_NONE)
        return heToMidpointIdx[he.pair];

      Vector pA = he.prev != HE_NONE
                      ? mesh.vertices[heMesh.halfEdges[he.prev].vertex]
                      : mesh.vertices[heMesh.halfEdges[he.pair].vertex];
      Vector pB = mesh.vertices[he.vertex];
      Vector mid = (pA + pB) * 0.5f;
      mid = mid.normalize();

      compiled.staticVertices.push_back(mid);
      uint16_t idx = static_cast<uint16_t>(compiled.staticVertices.size() - 1);
      heToMidpointIdx[heIdx] = idx;
      if (he.pair != HE_NONE)
        heToMidpointIdx[he.pair] = idx;
      return idx;
    };

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      getMidpointIdx(static_cast<uint16_t>(i));
    }

    compiled.staticOffset = static_cast<int>(compiled.staticVertices.size());

    // Star faces
    for (size_t i = 0; i < heMesh.faces.size(); ++i) {
      HEFace &face = heMesh.faces[i];
      uint16_t heIdx = face.halfEdge;
      uint16_t startHe = heIdx;
      int count = 0;

      if (heIdx == HE_NONE)
        continue;

      do {
        count += 2;
        HalfEdge &currHe = heMesh.halfEdges[heIdx];
        uint16_t prevIdx = currHe.prev;

        int idxM1 = getMidpointIdx(prevIdx);
        int idxM2 = getMidpointIdx(heIdx);

        HalfEdge &prevHe = heMesh.halfEdges[prevIdx];

        uint16_t iCorner = prevHe.vertex;
        uint16_t iPrev = prevHe.prev != HE_NONE
                             ? heMesh.halfEdges[prevHe.prev].vertex
                             : heMesh.halfEdges[prevHe.pair].vertex;
        uint16_t iNext = currHe.vertex;

        compiled.dynamicInstructions.push_back({iCorner, iPrev, iNext,
                                                static_cast<uint16_t>(idxM1),
                                                static_cast<uint16_t>(idxM2)});

        int16_t dynIdx = static_cast<int16_t>(compiled.dynamicVertices.size());
        heToDynamicIdx[heIdx] = dynIdx;
        compiled.dynamicVertices.emplace_back();

        compiled.faces.push_back(idxM1);
        compiled.faces.push_back(compiled.staticOffset + dynIdx);

        heIdx = currHe.next;
      } while (heIdx != HE_NONE && heIdx != startHe);

      compiled.face_counts.push_back(count);
    }

    // Rosette faces
    bool *visitedVerts = static_cast<bool *>(
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
      int16_t face_indices[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (count < 100)
          face_indices[count++] = heToMidpointIdx[currIdx];
        uint16_t nextEdgeIdx = currHe.pair != HE_NONE
                                   ? heMesh.halfEdges[currHe.pair].next
                                   : HE_NONE;
        if (nextEdgeIdx == HE_NONE)
          break;
        if (count < 100)
          face_indices[count++] =
              compiled.staticOffset + heToDynamicIdx[nextEdgeIdx];
        currIdx = nextEdgeIdx;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

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

  bool is_flat = std::abs(angle) < 1e-4f;

  for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
    const auto &instr = compiled.dynamicInstructions[i];
    Vector pCorner = compiled.baseVertices[instr.vCorner];

    if (is_flat) {
      compiled.dynamicVertices[i] = pCorner.normalize();
      continue;
    }

    Vector m1 = compiled.staticVertices[instr.idxM1];
    Vector m2 = compiled.staticVertices[instr.idxM2];
    Vector pPrev = compiled.baseVertices[instr.vPrev];
    Vector pNext = compiled.baseVertices[instr.vNext];

    Vector cross1 = cross(pPrev, pCorner);
    Vector cross2 = cross(pCorner, pNext);

    if (dot(cross1, cross1) < 1e-8f || dot(cross2, cross2) < 1e-8f) {
      compiled.dynamicVertices[i] = pCorner.normalize();
      continue; // zero length edge
    }

    Vector nEdge1 = cross1.normalize();
    Quaternion q1 = make_rotation(m1, angle);
    Vector nHankin1 = rotate(nEdge1, q1);

    Vector nEdge2 = cross2.normalize();
    Quaternion q2 = make_rotation(m2, -angle);
    Vector nHankin2 = rotate(nEdge2, q2);

    Vector intersect = cross(nHankin1, nHankin2);
    float lenSq = dot(intersect, intersect);
    if (lenSq < 1e-6f)
      intersect = (m1 + m2).normalize();
    if (dot(intersect, pCorner) < 0)
      intersect = -intersect;

    compiled.dynamicVertices[i] = intersect.normalize();
  }

  out_mesh.vertices.initialize(target_arena,
                               compiled.staticVertices.size() +
                                   compiled.dynamicVertices.size());
  for (size_t i = 0; i < compiled.staticVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.staticVertices[i]);
  for (size_t i = 0; i < compiled.dynamicVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.dynamicVertices[i]);

  out_mesh.face_counts.initialize(target_arena, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.initialize(target_arena, compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.initialize(target_arena, compiled.faces.size());
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    out_mesh.faces.push_back(compiled.faces[i]);
}

inline PolyMesh hankin(const PolyMesh &mesh, ScratchContext &ctx, float angle) {
  PolyMesh out;
  CompiledHankin compiled;

  compile_hankin(mesh, compiled, ctx);
  update_hankin(compiled, out, *(ctx.target), angle);

  ctx.swap_and_clear();
  return out;
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
inline PolyMesh kis(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(*(ctx.target), V + F);
  out_mesh.face_counts.initialize(*(ctx.target), I);
  out_mesh.faces.initialize(*(ctx.target), 3 * I);

  {
    ArenaMarker _(*(ctx.source));

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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Ambo operator: Truncates vertices to edge midpoints.
 */
/**
 * @brief Ambo operator: Truncates vertices to edge midpoints.
 */
inline PolyMesh ambo(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(*(ctx.target), E);
  out_mesh.face_counts.initialize(*(ctx.target), F + V);
  out_mesh.faces.initialize(*(ctx.target), 2 * I);

  {
    ArenaMarker _(*(ctx.source));

    // 1. Build HE Mesh FIRST
    HalfEdgeMesh heMesh(*(ctx.source), mesh);

    uint16_t *edgeToVert = static_cast<uint16_t *>(
        ctx.source->allocate(I * sizeof(uint16_t), alignof(uint16_t)));
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
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Truncate operator: Cuts corners off the polyhedron.
 * @param t Truncation depth [0..0.5].
 */
inline PolyMesh truncate(const PolyMesh &mesh, ScratchContext &ctx,
                         float t = 0.25f) {
  if (std::abs(t - 0.5f) < 1e-4f) {
    return ambo(mesh, ctx);
  }

  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(*(ctx.target), 2 * E);
  out_mesh.face_counts.initialize(*(ctx.target), F + V);
  out_mesh.faces.initialize(*(ctx.target), 3 * I);

  {
    ArenaMarker _(*(ctx.source));

    HalfEdgeMesh heMesh(*(ctx.source), mesh);

    std::pair<int16_t, int16_t> *edgeToVert =
        static_cast<std::pair<int16_t, int16_t> *>(
            ctx.source->allocate(I * sizeof(std::pair<int16_t, int16_t>),
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
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Expand operator: Separates faces (e = aa).
 * @param t Expansion factor. Default 2-sqrt(2) ~= 0.5857.
 */
inline PolyMesh expand(const PolyMesh &mesh, ScratchContext &ctx,
                       float t = 2.0f - sqrt(2.0f)) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(*(ctx.target), I);
  out_mesh.face_counts.initialize(*(ctx.target), F + V + E);
  out_mesh.faces.initialize(*(ctx.target), 4 * I);

  {
    ArenaMarker _(*(ctx.source));

    HalfEdgeMesh heMesh(*(ctx.source), mesh);
    int16_t *heToVertIdx = static_cast<int16_t *>(
        ctx.source->allocate(I * sizeof(int16_t), alignof(int16_t)));
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
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
        ctx.source->allocate(I * sizeof(bool), alignof(bool)));
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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Bitruncate operator: Truncate the rectified mesh.
 */
inline PolyMesh bitruncate(const PolyMesh &mesh, ScratchContext &ctx,
                           float t = 1.0f / 3.0f) {
  return truncate(ambo(mesh, ctx), ctx, t);
}

inline PolyMesh canonicalize(const PolyMesh &mesh, ScratchContext &ctx,
                             int iterations = 8) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(*(ctx.target), V);
  out_mesh.face_counts.initialize(*(ctx.target), F);
  out_mesh.faces.initialize(*(ctx.target), I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);
  for (size_t i = 0; i < F; ++i)
    out_mesh.face_counts.push_back(mesh.face_counts[i]);
  for (size_t i = 0; i < I; ++i)
    out_mesh.faces.push_back(mesh.faces[i]);

  {
    ArenaMarker _(*(ctx.source));

    ArenaVector<Vector> movements;
    movements.initialize(*(ctx.source), V);
    for (size_t i = 0; i < V; ++i)
      movements.push_back(Vector(0, 0, 0));

    HalfEdgeMesh heMesh(*(ctx.source),
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
  ctx.swap_and_clear();
  return out_mesh;
}

/**
 * @brief Snub operator: Creates a chiral semi-regular polyhedron.
 * Updated with twist support.
 */
inline PolyMesh snub(const PolyMesh &mesh, ScratchContext &ctx, float t = 0.5f,
                     float twist = 0.0f) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(*(ctx.target), I);
  out_mesh.face_counts.initialize(*(ctx.target), F + V + 2 * E);
  out_mesh.faces.initialize(*(ctx.target), 5 * I);

  {
    ArenaMarker _(*(ctx.source));

    HalfEdgeMesh heMesh(*(ctx.source), mesh);
    uint16_t *heToVertIdx = static_cast<uint16_t *>(
        ctx.source->allocate(I * sizeof(uint16_t), alignof(uint16_t)));
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
        ctx.source->allocate(V * sizeof(bool), alignof(bool)));
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
        ctx.source->allocate(I * sizeof(bool), alignof(bool)));
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
  ctx.swap_and_clear();
  return out_mesh;
}

inline PolyMesh gyro(const PolyMesh &mesh, ScratchContext &ctx) {
  return dual(snub(mesh, ctx), ctx);
}

/**
 * @brief Computes KDTree and Adjacency map for the mesh (caching it).
 */
inline void compute_kdtree(const PolyMesh &mesh, Arena &arena) {
  if (mesh.cache_valid)
    return;

  PolyMesh &m = const_cast<PolyMesh &>(mesh);
  m.kdTree =
      KDTree(arena, std::span<Vector>(m.vertices.data(), m.vertices.size()));
  m.cache_valid = true;
}

inline Vector closest_point_on_mesh_graph(const Vector &p, const PolyMesh &mesh,
                                          Arena &scratch_arena_a) {
  if (mesh.vertices.empty())
    return Vector(0, 1, 0);

  compute_kdtree(mesh, scratch_arena_a);

  auto nearestNodes = mesh.kdTree.nearest(p, 1);
  if (nearestNodes.size() == 0)
    return mesh.vertices[0];

  const auto &node = nearestNodes[0];
  int closestVertexIndex = (int)node.originalIndex;
  Vector closestVertexPos = node.point;

  Vector bestPoint = closestVertexPos;
  float maxDot = dot(p, bestPoint);

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
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

/**
 * @brief Helper to finish hash.
 */
static uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

// A simple hash combine block
static inline void hash_combine(uint32_t &seed, uint32_t v) {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief Colors faces based on their vertex count and neighbor topology.
 */
template <typename MeshT>
static std::vector<int> classify_faces_by_topology(const MeshT &mesh,
                                                   Arena &scratch_arena_a) {
  ArenaMarker _(scratch_arena_a);

  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  ArenaVector<uint32_t> faceHashes;
  faceHashes.initialize(scratch_arena_a, F);

  ArenaVector<uint32_t> finalHashes;
  finalHashes.initialize(scratch_arena_a, F);

  size_t offset = 0;
  for (size_t i = 0; i < F; ++i) {
    int count = mesh.face_counts[i];
    Vector verts[100];
    for (int k = 0; k < count; ++k) {
      if (k < 100)
        verts[k] = mesh.vertices[mesh.faces[offset + k]];
    }

    int vertexCount = std::min(count, 100);
    int angles[100];
    if (vertexCount >= 3) {
      for (int k = 0; k < vertexCount; ++k) {
        const Vector &prev = verts[(k - 1 + vertexCount) % vertexCount];
        const Vector &curr = verts[k];
        const Vector &next = verts[(k + 1) % vertexCount];
        Vector v1 = (prev - curr).normalize();
        Vector v2 = (next - curr).normalize();
        float ang = angle_between(v1, v2);
        angles[k] = (int)std::round(ang * 180.0f / PI_F);
      }
      std::sort(angles, angles + vertexCount);
    }

    uint32_t h = 0x12345678;
    hash_combine(h, static_cast<uint32_t>(vertexCount));
    for (int k = 0; k < vertexCount; ++k) {
      hash_combine(h, static_cast<uint32_t>(angles[k]));
    }
    h = fmix32(h);
    faceHashes.push_back(h);
    finalHashes.push_back(h);
    offset += count;
  }

  {
    ArenaMarker temp_topo(scratch_arena_a);

    uint16_t *heToFace = static_cast<uint16_t *>(
        scratch_arena_a.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *pairArray = static_cast<uint16_t *>(
        scratch_arena_a.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(pairArray, I, HE_NONE);

    {
      ArenaMarker temp_records(scratch_arena_a);
      struct EdgeRecord {
        uint16_t min_v, max_v, he;
      };
      EdgeRecord *records = static_cast<EdgeRecord *>(scratch_arena_a.allocate(
          I * sizeof(EdgeRecord), alignof(EdgeRecord)));

      size_t he_idx = 0;
      size_t face_offset = 0;
      for (size_t fi = 0; fi < F; ++fi) {
        int count = mesh.face_counts[fi];
        for (int k = 0; k < count; ++k) {
          uint16_t u = mesh.faces[face_offset + k];
          uint16_t v = mesh.faces[face_offset + (k + 1) % count];
          records[he_idx].min_v = std::min(u, v);
          records[he_idx].max_v = std::max(u, v);
          records[he_idx].he = static_cast<uint16_t>(he_idx);
          heToFace[he_idx] = static_cast<uint16_t>(fi);
          he_idx++;
        }
        face_offset += count;
      }

      std::sort(records, records + I,
                [](const EdgeRecord &a, const EdgeRecord &b) {
                  if (a.min_v != b.min_v)
                    return a.min_v < b.min_v;
                  return a.max_v < b.max_v;
                });

      for (size_t i = 0; i < I;) {
        if (i + 1 < I && records[i].min_v == records[i + 1].min_v &&
            records[i].max_v == records[i + 1].max_v) {
          pairArray[records[i].he] = records[i + 1].he;
          pairArray[records[i + 1].he] = records[i].he;
          i += 2;
        } else {
          i += 1;
        }
      }
    }

    offset = 0;
    for (size_t fi = 0; fi < F; ++fi) {
      int count = mesh.face_counts[fi];
      uint32_t neighborAcc = 0;
      for (int k = 0; k < count; ++k) {
        uint16_t pIdx = pairArray[offset + k];
        if (pIdx != HE_NONE) {
          uint32_t neigh_h = faceHashes[heToFace[pIdx]];
          hash_combine(neigh_h, 0);
          neighborAcc += fmix32(neigh_h);
        }
      }
      uint32_t final_h = faceHashes[fi];
      hash_combine(final_h, neighborAcc);
      finalHashes[fi] = fmix32(final_h);
      offset += count;
    }
  }

  struct HashNode {
    uint32_t hash;
    int original_face;
  };
  HashNode *nodes = static_cast<HashNode *>(
      scratch_arena_a.allocate(F * sizeof(HashNode), alignof(HashNode)));

  for (size_t fi = 0; fi < F; ++fi) {
    nodes[fi] = {finalHashes[fi], static_cast<int>(fi)};
  }

  std::sort(nodes, nodes + F, [](const HashNode &a, const HashNode &b) {
    return a.hash < b.hash;
  });

  std::vector<int> faceColorIndices(F, 0);
  if (F > 0) {
    int currentID = 0;
    faceColorIndices[nodes[0].original_face] = 0;
    for (size_t i = 1; i < F; ++i) {
      if (nodes[i].hash != nodes[i - 1].hash) {
        currentID++;
      }
      faceColorIndices[nodes[i].original_face] = currentID;
    }
  }

  return faceColorIndices;
}

} // namespace MeshOps
