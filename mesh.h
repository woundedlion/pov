/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "concepts.h"
#include "3dmath.h"
#include "spatial.h"
#include "memory.h"

#include <map>
#include <set>
#include <algorithm>
#include <span>
#include <cmath>

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
  ArenaVector<int> topology;

  // Cache
  mutable bool cache_valid = false;
  mutable KDTree kdTree;

  PolyMesh() = default;

  void initialize(Arena &arena, size_t num_verts, size_t num_faces,
                  size_t num_indices) {
    vertices.bind(arena, num_verts);
    face_counts.bind(arena, num_faces);
    faces.bind(arena, num_indices);
    cache_valid = false;
  }

  inline void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    topology.clear();
    cache_valid = false;
  }

  void clear_cache() const {
    cache_valid = false;
    kdTree.clear();
  }

  // Unified accessors (PolyMesh always owns, so these just forward)
  const uint8_t *get_face_counts_data() const { return face_counts.data(); }
  size_t get_face_counts_size() const { return face_counts.size(); }

  const uint16_t *get_faces_data() const { return faces.data(); }
  size_t get_faces_size() const { return faces.size(); }
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

    vertices.bind(arena, num_verts);
    for (size_t i = 0; i < num_verts; ++i) {
      vertices.push_back({HE_NONE});
    }

    faces.bind(arena, num_faces);
    halfEdges.bind(arena, total_indices);

    struct EdgeRecord {
      uint16_t min_v;
      uint16_t max_v;
      uint16_t he;
    };

    size_t face_offset = 0;
    size_t he_idx = 0;

    {
      ScratchScope _(arena);
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

namespace MeshOps {

/**
 * @brief Compiles a PolyMesh into a static MeshState.
 * Removes degenerate faces (faces with < 3 vertices) during
 * the process.
 */
FLASHMEM inline void compile(const PolyMesh &src, MeshState &dst,
                             Arena &geom_arena) {
  dst.clear();

  size_t valid_faces = 0;
  size_t valid_indices = 0;

  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    if (src.face_counts[i] >= 3) {
      valid_faces++;
      valid_indices += src.face_counts[i];
    }
  }

  dst.vertices.bind(geom_arena, src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    dst.vertices.push_back(src.vertices[i]);
  }

  dst.face_counts.bind(geom_arena, valid_faces);
  dst.faces.bind(geom_arena, valid_indices);
  dst.face_offsets.bind(geom_arena, valid_faces);

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
  dst.vertices.bind(arena, src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    dst.vertices.push_back(src.vertices[i]);
  }

  size_t fc_size = src.get_face_counts_size();
  const uint8_t *fc_data = src.get_face_counts_data();
  dst.face_counts.bind(arena, fc_size);
  for (size_t i = 0; i < fc_size; ++i) {
    dst.face_counts.push_back(fc_data[i]);
  }

  size_t f_size = src.get_faces_size();
  const uint16_t *f_data = src.get_faces_data();
  dst.faces.bind(arena, f_size);
  for (size_t i = 0; i < f_size; ++i) {
    dst.faces.push_back(f_data[i]);
  }

  if constexpr (requires { dst.face_offsets; }) {
    size_t fo_size = src.get_face_offsets_size();
    const uint16_t *fo_data = src.get_face_offsets_data();
    dst.face_offsets.bind(arena, fo_size);
    for (size_t i = 0; i < fo_size; ++i) {
      dst.face_offsets.push_back(fo_data[i]);
    }
  }

  if constexpr (requires { dst.topology; }) {
    dst.topology.bind(arena, src.topology.size());
    for (size_t i = 0; i < src.topology.size(); ++i) {
      dst.topology.push_back(src.topology[i]);
    }
  }
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
FLASHMEM __attribute__((noinline)) static void
classify_faces_by_topology(MeshT &mesh, Arena &scratch_a, Arena &scratch_b,
                           Arena &persistent) {
  ScratchScope _(scratch_a);

  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  ArenaVector<uint32_t> faceHashes;
  faceHashes.bind(scratch_a, F);

  ArenaVector<uint32_t> finalHashes;
  finalHashes.bind(scratch_a, F);

  // Find max face vertex count for scratch allocation
  int maxCount = 0;
  for (size_t i = 0; i < F; ++i) {
    int c = mesh.face_counts[i];
    if (c > maxCount) maxCount = c;
  }

  ArenaVector<Vector> verts;
  verts.bind(scratch_a, maxCount);
  ArenaVector<int> angles;
  angles.bind(scratch_a, maxCount);

  size_t offset = 0;
  for (size_t i = 0; i < F; ++i) {
    int count = mesh.face_counts[i];

    verts.clear();
    for (int k = 0; k < count; ++k) {
      verts.push_back(mesh.vertices[mesh.faces[offset + k]]);
    }

    angles.clear();
    if (count >= 3) {
      for (int k = 0; k < count; ++k) {
        const Vector &prev = verts[(k - 1 + count) % count];
        const Vector &curr = verts[k];
        const Vector &next = verts[(k + 1) % count];
        Vector v1 = (prev - curr).normalize();
        Vector v2 = (next - curr).normalize();
        float ang = angle_between(v1, v2);
        angles.push_back((int)std::round(ang * 180.0f / PI_F));
      }
      std::sort(angles.data(), angles.data() + count);
    }

    uint32_t h = 0x12345678;
    hash_combine(h, static_cast<uint32_t>(count));
    for (int k = 0; k < count; ++k) {
      hash_combine(h, static_cast<uint32_t>(angles[k]));
    }
    h = fmix32(h);
    faceHashes.push_back(h);
    finalHashes.push_back(h);
    offset += count;
  }

  {
    ScratchScope temp_topo(scratch_a);

    uint16_t *heToFace =
        static_cast<uint16_t *>(scratch_a.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *pairArray =
        static_cast<uint16_t *>(scratch_a.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(pairArray, I, HE_NONE);

    {
      ScratchScope temp_records(scratch_b);
      struct EdgeRecord {
        uint16_t min_v, max_v, he;
      };
      EdgeRecord *records =
          static_cast<EdgeRecord *>(scratch_b.allocate(
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
      scratch_a.allocate(F * sizeof(HashNode), alignof(HashNode)));

  for (size_t fi = 0; fi < F; ++fi) {
    nodes[fi] = {finalHashes[fi], static_cast<int>(fi)};
  }

  std::sort(nodes, nodes + F, [](const HashNode &a, const HashNode &b) {
    return a.hash < b.hash;
  });

  if (mesh.topology.capacity() < F) {
    mesh.topology.bind(persistent, F);
  }
  mesh.topology.clear();
  for (size_t i = 0; i < F; ++i) {
    mesh.topology.push_back(0);
  }

  if (F > 0) {
    int currentID = 0;
    mesh.topology[nodes[0].original_face] = 0;
    for (size_t i = 1; i < F; ++i) {
      if (nodes[i].hash != nodes[i - 1].hash) {
        currentID++;
      }
      mesh.topology[nodes[i].original_face] = currentID;
    }
  }
}

} // namespace MeshOps
