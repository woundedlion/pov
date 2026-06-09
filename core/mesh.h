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
#include <type_traits>

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

  PolyMesh() = default;

  inline void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    topology.clear();
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
  uint16_t half_edge =
      HE_NONE; /**< One of the half-edges pointing to this vertex. */
};

struct HEFace {
  uint16_t half_edge =
      HE_NONE; /**< One of the half-edges bordering this face. */
};

/**
 * @brief Record used to pair opposite half-edges by their undirected
 * (min_v, max_v) vertex key.
 */
struct HalfEdgePairRecord {
  uint16_t min_v, max_v, he;
};

/**
 * @brief Sorts half-edge records by (min_v, max_v) and calls set_pair(heA, heB)
 * once per matched opposite-edge pair. Templated + inline, so it lowers to the
 * same code as a hand-inlined sort/scan — no call overhead on the (cold)
 * mesh-build path; this only removes the copy-paste between build_from_flat and
 * classify_faces_by_topology.
 */
template <typename SetPairFn>
inline void pair_half_edges(HalfEdgePairRecord *records, size_t n,
                            SetPairFn set_pair) {
  auto edge_less = [](const HalfEdgePairRecord &a, const HalfEdgePairRecord &b) {
    if (a.min_v != b.min_v) return a.min_v < b.min_v;
    return a.max_v < b.max_v;
  };
  auto edge_greater = [&](const HalfEdgePairRecord &a,
                          const HalfEdgePairRecord &b) {
    return edge_less(b, a);
  };
  std::make_heap(records, records + n, edge_greater);
  std::sort_heap(records, records + n, edge_greater);

  for (size_t i = 0; i < n;) {
    size_t j = i + 1;
    while (j < n && records[j].min_v == records[i].min_v &&
           records[j].max_v == records[i].max_v)
      ++j;
    // A 2-manifold mesh shares each undirected edge between exactly 1 (boundary)
    // or 2 (interior) half-edges. A run of >2 is a non-manifold edge — a build
    // invariant violation — so fail fast rather than pair two arbitrarily and
    // silently drop the rest.
    HS_CHECK(j - i <= 2, "non-manifold edge: >2 half-edges share an edge");
    if (j - i == 2)
      set_pair(records[i].he, records[i + 1].he);
    i = j;
  }
}

class HalfEdgeMesh {
public:
  ArenaVector<HEVertex> vertices;
  ArenaVector<HEFace> faces;
  ArenaVector<HalfEdge> half_edges;

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

    // All vertex/face/half-edge indices are stored as uint16_t below; a mesh
    // past that range would silently truncate. Trap instead (fail-fast).
    HS_CHECK(num_verts <= UINT16_MAX && total_indices <= UINT16_MAX,
             "half-edge mesh exceeds 16-bit index range");

    vertices.bind(arena, num_verts);
    for (size_t i = 0; i < num_verts; ++i) {
      vertices.push_back({HE_NONE});
    }

    faces.bind(arena, num_faces);
    half_edges.bind(arena, total_indices);

    size_t face_offset = 0;
    size_t he_idx = 0;

    {
      ScratchScope _(arena);
      HalfEdgePairRecord *records =
          static_cast<HalfEdgePairRecord *>(arena.allocate(
              total_indices * sizeof(HalfEdgePairRecord),
              alignof(HalfEdgePairRecord)));

      for (size_t fi = 0; fi < num_faces; ++fi) {
        int count = counts[fi];

        faces.emplace_back();
        uint16_t current_face_idx = static_cast<uint16_t>(faces.size() - 1);
        size_t face_start_he_idx = he_idx;

        for (int i = 0; i < count; ++i) {
          half_edges.emplace_back();
        }

        for (int i = 0; i < count; ++i) {
          uint16_t u = faces_arr[face_offset + i];
          uint16_t v = faces_arr[face_offset + (i + 1) % count];

          uint16_t he_index = static_cast<uint16_t>(face_start_he_idx + i);
          HalfEdge &he = half_edges[he_index];
          he.vertex = v;
          he.face = current_face_idx;
          he.next = static_cast<uint16_t>(face_start_he_idx + (i + 1) % count);
          he.prev =
              static_cast<uint16_t>(face_start_he_idx + (i - 1 + count) % count);

          vertices[v].half_edge = he_index;

          records[he_idx].min_v = std::min(u, v);
          records[he_idx].max_v = std::max(u, v);
          records[he_idx].he = he_index;

          he_idx++;
        }
        faces[current_face_idx].half_edge = static_cast<uint16_t>(face_start_he_idx);
        face_offset += count;
      }

      pair_half_edges(records, total_indices, [&](uint16_t a, uint16_t b) {
        half_edges[a].pair = b;
        half_edges[b].pair = a;
      });
    }
  }
};

namespace MeshOps {

// ---------------------------------------------------------------------------
// Shared topology helpers
// ---------------------------------------------------------------------------

/**
 * @brief Narrow a freshly-appended container index to the mesh topology index
 * type, trapping on a budget bump that pushes it past the representable range.
 *
 * Vertex/face indices are stored as uint16_t (faces, with HE_NONE=0xFFFF as the
 * sentinel) and, in some operators' scratch, as int16_t (-1 sentinel). The input
 * range is checked at mesh build, but each operator's OUTPUT index is captured
 * via `static_cast<...>(vertices.size() - 1)`, which silently wraps if a future
 * MAX_VERTS bump makes the count exceed the index type. Routing those casts
 * through here converts that silent geometry corruption into a bench-time trap.
 * Bound is INT16_MAX — the narrowest type these indices land in (matches the
 * `MAX_VERTS <= INT16_MAX` static_assert in solids.h). Cold path: once per
 * emitted vertex during a rebuild, never per pixel.
 */
inline uint16_t narrow_index(size_t i) {
  HS_CHECK(i <= static_cast<size_t>(INT16_MAX),
           "mesh index exceeds int16_t topology range (MAX_VERTS bumped?)");
  return static_cast<uint16_t>(i);
}

/**
 * @brief Trapping int -> uint8_t cast for a face's side count.
 *
 * The Conway/Hankin operators accumulate a face's vertex count in an `int`, then
 * push it into the `uint8_t face_counts` array. A face wider than 255 sides —
 * reachable as a high-valence vertex orbit (`orbit_count`), a doubled edge count
 * (`count * 2` in snub), or a long Hankin rosette walk (bounded only by 2*I) —
 * would silently wrap, corrupting the face topology. Mirror of `narrow_index`:
 * route those casts through here so an over-wide face traps at the bench instead
 * of emitting garbage geometry. Cold path: once per emitted face during a
 * rebuild, never per pixel.
 */
inline uint8_t narrow_face_count(int count) {
  HS_CHECK(count >= 0 && count <= UINT8_MAX,
           "mesh face side count exceeds uint8_t range");
  return static_cast<uint8_t>(count);
}

/**
 * @brief Compiles a PolyMesh into a static MeshState.
 * Removes degenerate faces (faces with < 3 vertices) during
 * the process.
 */
FLASHMEM static void compile(const PolyMesh &src, MeshState &dst,
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
 * @brief Performs a strict deep copy of a mesh into a target arena.
 * Safe for memory compaction and bouncing between arenas.
 *
 * For MeshState this delegates to MeshState::clone (the Cloneable interface
 * used by Persist) so the two implementations can't drift; the generic body
 * below handles PolyMesh (which has no face_offsets).
 */
template <typename MeshT>
inline void clone(const MeshT &src, MeshT &dst, Arena &arena) {
  if constexpr (std::is_same_v<MeshT, MeshState>) {
    MeshState::clone(src, dst, arena);
    return;
  }
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

// Mixes a value into a running hash seed.
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

  ArenaVector<uint32_t> face_hashes;
  face_hashes.bind(scratch_a, F);

  ArenaVector<uint32_t> final_hashes;
  final_hashes.bind(scratch_a, F);

  // Find max face vertex count for scratch allocation
  int max_count = 0;
  for (size_t i = 0; i < F; ++i) {
    int c = mesh.face_counts[i];
    if (c > max_count) max_count = c;
  }

  ArenaVector<Vector> verts;
  verts.bind(scratch_a, max_count);
  ArenaVector<int> angles;
  angles.bind(scratch_a, max_count);

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
        Vector v1 = (prev - curr).normalized();
        Vector v2 = (next - curr).normalized();
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
    face_hashes.push_back(h);
    final_hashes.push_back(h);
    offset += count;
  }

  {
    ScratchScope temp_topo(scratch_a);

    uint16_t *he_to_face =
        static_cast<uint16_t *>(scratch_a.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    uint16_t *pair_array =
        static_cast<uint16_t *>(scratch_a.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(pair_array, I, HE_NONE);

    {
      ScratchScope temp_records(scratch_b);
      HalfEdgePairRecord *records =
          static_cast<HalfEdgePairRecord *>(scratch_b.allocate(
              I * sizeof(HalfEdgePairRecord), alignof(HalfEdgePairRecord)));

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
          he_to_face[he_idx] = static_cast<uint16_t>(fi);
          he_idx++;
        }
        face_offset += count;
      }

      pair_half_edges(records, I, [&](uint16_t a, uint16_t b) {
        pair_array[a] = b;
        pair_array[b] = a;
      });
    }

    offset = 0;
    for (size_t fi = 0; fi < F; ++fi) {
      int count = mesh.face_counts[fi];
      uint32_t neighbor_acc = 0;
      for (int k = 0; k < count; ++k) {
        uint16_t p_idx = pair_array[offset + k];
        if (p_idx != HE_NONE) {
          uint32_t neigh_h = face_hashes[he_to_face[p_idx]];
          hash_combine(neigh_h, 0);
          neighbor_acc += fmix32(neigh_h);
        }
      }
      uint32_t final_h = face_hashes[fi];
      hash_combine(final_h, neighbor_acc);
      final_hashes[fi] = fmix32(final_h);
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
    nodes[fi] = {final_hashes[fi], static_cast<int>(fi)};
  }

  auto hash_greater = [](const HashNode &a, const HashNode &b) {
    return a.hash > b.hash;
  };
  std::make_heap(nodes, nodes + F, hash_greater);
  std::sort_heap(nodes, nodes + F, hash_greater);

  if (mesh.topology.capacity() < F) {
    mesh.topology.bind(persistent, F);
  }
  mesh.topology.clear();
  for (size_t i = 0; i < F; ++i) {
    mesh.topology.push_back(0);
  }

  if (F > 0) {
    int current_id = 0;
    mesh.topology[nodes[0].original_face] = 0;
    for (size_t i = 1; i < F; ++i) {
      if (nodes[i].hash != nodes[i - 1].hash) {
        current_id++;
      }
      mesh.topology[nodes[i].original_face] = current_id;
    }
  }
}

} // namespace MeshOps
