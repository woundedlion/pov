/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "concepts.h"
#include "3dmath.h"
#include "spatial.h"
#include "memory.h"

#include <algorithm>
#include <cmath>
#include <type_traits>

// Forward declarations
struct HEVertex;
struct HEFace;
struct HalfEdge;

struct MeshState;

/**
 * @brief A simple dynamic mesh structure compatible with MeshOps templates.
 */
struct PolyMesh {
  ArenaVector<Vector> vertices;     /**< Vertex positions. */
  ArenaVector<uint8_t> face_counts; /**< Number of sides for each face. */
  ArenaVector<uint16_t> faces;      /**< Flat per-face vertex index list. */
  ArenaVector<int> topology;        /**< Per-face topology class id. */

  /**
   * @brief Constructs an empty mesh with no allocated storage.
   */
  PolyMesh() = default;

  /**
   * @brief Resets all arrays to empty without releasing arena storage.
   */
  inline void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    topology.clear();
  }

  /**
   * @brief Returns a pointer to the face side-count array.
   * @return Pointer to the contiguous uint8_t per-face side counts.
   * @details Unified accessor; PolyMesh always owns its data, so this forwards.
   */
  const uint8_t *get_face_counts_data() const { return face_counts.data(); }

  /**
   * @brief Returns the number of faces.
   * @return Count of entries in the face side-count array.
   */
  size_t get_face_counts_size() const { return face_counts.size(); }

  /**
   * @brief Returns a pointer to the flat face index array.
   * @return Pointer to the contiguous uint16_t per-face vertex indices.
   * @details Unified accessor; PolyMesh always owns its data, so this forwards.
   */
  const uint16_t *get_faces_data() const { return faces.data(); }

  /**
   * @brief Returns the total number of face vertex indices.
   * @return Count of entries in the flat face index array.
   */
  size_t get_faces_size() const { return faces.size(); }
};

constexpr uint16_t HE_NONE = 0xFFFF; /**< Null index sentinel for half-edge connectivity. */

/**
 * @brief One directed half of an undirected edge in the half-edge mesh.
 */
struct HalfEdge {
  uint16_t vertex = HE_NONE; /**< Vertex at the end of this half-edge. */
  uint16_t face = HE_NONE;   /**< Face this half-edge belongs to. */
  uint16_t next = HE_NONE;   /**< Next half-edge in the face loop. */
  uint16_t prev = HE_NONE;   /**< Previous half-edge in the face loop. */
  uint16_t pair = HE_NONE;   /**< Opposite half-edge. */
};

/**
 * @brief Mesh vertex; entry point into the half-edge ring around it.
 */
struct HEVertex {
  uint16_t half_edge =
      HE_NONE; /**< One of the half-edges pointing to this vertex. */
};

/**
 * @brief Mesh face; entry point into its bordering half-edge loop.
 */
struct HEFace {
  uint16_t half_edge =
      HE_NONE; /**< One of the half-edges bordering this face. */
};

/**
 * @brief Record used to pair opposite half-edges by their undirected
 * (min_v, max_v) vertex key.
 */
struct HalfEdgePairRecord {
  uint16_t min_v, max_v, he; /**< Lower vertex index, upper vertex index, and the half-edge index. */
};

/**
 * @brief Fills one half-edge record with the canonical undirected edge key.
 * @param rec Record to populate.
 * @param u One edge endpoint vertex index.
 * @param v The other edge endpoint vertex index.
 * @param he The half-edge index this record represents.
 * @details Centralizes the (min_v, max_v) key construction so the two record-
 * building loops (HalfEdgeMesh::build_from_flat and classify_faces_by_topology)
 * cannot drift if the key convention ever changes. Inline — no call overhead.
 */
inline void fill_edge_record(HalfEdgePairRecord &rec, uint16_t u, uint16_t v,
                             uint16_t he) {
  rec.min_v = std::min(u, v);
  rec.max_v = std::max(u, v);
  rec.he = he;
}

/**
 * @brief Sorts half-edge records by (min_v, max_v) and calls set_pair(heA, heB)
 * once per matched opposite-edge pair.
 * @tparam SetPairFn Callable invoked as set_pair(heA, heB) to link a matched pair.
 * @param records Array of half-edge pair records; sorted in place by edge key.
 * @param n Number of records in the array.
 * @param set_pair Callback linking the two half-edges of each interior edge.
 * @details Templated + inline, so it lowers to the same code as a hand-inlined
 * sort/scan — no call overhead on the cold mesh-build path; this only removes
 * the copy-paste between build_from_flat and classify_faces_by_topology. Traps
 * on a non-manifold edge (>2 half-edges sharing one undirected edge).
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

/**
 * @brief Non-owning (pointer + size) view exposing just the size()/operator[]
 * surface build_from_flat needs.
 * @tparam T Element type of the viewed array.
 * @details Lets the MeshState ctor feed topology through the unified
 * get_*_data()/get_*_size() accessors: a borrowed-mode MeshState serves its
 * topology through the view spans with the owned face_counts/faces left empty,
 * so reading the owned vectors directly would silently build an empty half-edge
 * mesh. POD by value — inlines away, no per-element cost.
 */
template <typename T> struct FlatView {
  const T *ptr; /**< Pointer to the first viewed element. */
  size_t n;     /**< Number of viewed elements. */

  /**
   * @brief Returns the number of viewed elements.
   * @return Element count of the view.
   */
  size_t size() const { return n; }

  /**
   * @brief Returns the element at the given index.
   * @param i Index in [0, size()).
   * @return Const reference to the viewed element.
   */
  const T &operator[](size_t i) const { return ptr[i]; }
};

/**
 * @brief Half-edge connectivity built from a flat face list, giving O(1)
 * adjacency (next/prev/pair) for traversal-based mesh operators.
 */
class HalfEdgeMesh {
public:
  ArenaVector<HEVertex> vertices;   /**< Per-vertex half-edge ring entry points. */
  ArenaVector<HEFace> faces;        /**< Per-face half-edge loop entry points. */
  ArenaVector<HalfEdge> half_edges; /**< Directed half-edges with next/prev/pair links. */

  /**
   * @brief Builds the half-edge connectivity from a PolyMesh.
   * @param arena Arena supplying storage for the half-edge arrays.
   * @param mesh Source mesh whose owned vertices/face_counts/faces are read.
   */
  explicit HalfEdgeMesh(Arena &arena, const PolyMesh &mesh) {
    build_from_flat(arena, mesh.vertices, mesh.face_counts, mesh.faces);
  }

  /**
   * @brief Builds the half-edge connectivity from a MeshState.
   * @param arena Arena supplying storage for the half-edge arrays.
   * @param mesh Source mesh; vertices are owned but topology may be borrowed.
   * @details MeshState's vertices are always owned, but its topology may be
   * borrowed (view spans); route it through the unified accessors so both
   * owned and borrowed modes work.
   */
  explicit HalfEdgeMesh(Arena &arena, const MeshState &mesh) {
    build_from_flat(
        arena, mesh.vertices,
        FlatView<uint8_t>{mesh.get_face_counts_data(),
                          mesh.get_face_counts_size()},
        FlatView<uint16_t>{mesh.get_faces_data(), mesh.get_faces_size()});
  }

private:
  /**
   * @brief Builds the half-edge structure from a flat (vertices, per-face side
   * counts, face index list) representation.
   * @tparam Verts Container of vertex positions (size()/operator[]).
   * @tparam Counts Container of per-face side counts (size()/operator[]).
   * @tparam Faces Container of the flat face vertex indices (size()/operator[]).
   * @param arena Arena supplying storage for the half-edge arrays and scratch.
   * @param verts Vertex positions; only its size() is used here.
   * @param counts Per-face side counts.
   * @param faces_arr Flat per-face vertex index list.
   * @details Emits one half-edge per face index, links each face's next/prev
   * loop, then pairs opposite half-edges by undirected vertex key. Templated
   * over the container types so it accepts both owned ArenaVectors and borrowed
   * FlatViews. Traps on a 16-bit index overflow or a zero-side face.
   */
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

        // A zero-count face emits no half-edges yet still gets its half_edge
        // set below, mis-linking it to the next face's first half-edge (or
        // one-past-the-end for a trailing zero-count face) — corrupt topology.
        // Trap it (fail-fast, matching the 16-bit index trap above). Note we
        // trap only count==0, not all <3 faces: compile() deliberately accepts
        // 1-/2-side degenerate faces and strips them downstream (they still
        // emit linked half-edges here), so a >=3 invariant would be wrong.
        HS_CHECK(count > 0, "half-edge mesh face has zero sides");

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

          fill_edge_record(records[he_idx], u, v, he_index);

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
 * @param i Container index to narrow.
 * @return The index as a uint16_t.
 * @details Vertex/face indices are stored as uint16_t (faces, with
 * HE_NONE=0xFFFF as the sentinel) and, in some operators' scratch, as int16_t
 * (-1 sentinel). The input range is checked at mesh build, but each operator's
 * OUTPUT index is captured via static_cast<...>(vertices.size() - 1), which
 * silently wraps if a future MAX_VERTS bump makes the count exceed the index
 * type. Routing those casts through here converts that silent geometry
 * corruption into a bench-time trap. Bound is INT16_MAX — the narrowest type
 * these indices land in (matches the MAX_VERTS <= INT16_MAX static_assert in
 * solids.h). Cold path: once per emitted vertex during a rebuild, never per
 * pixel.
 */
inline uint16_t narrow_index(size_t i) {
  HS_CHECK(i <= static_cast<size_t>(INT16_MAX),
           "mesh index exceeds int16_t topology range (MAX_VERTS bumped?)");
  return static_cast<uint16_t>(i);
}

/**
 * @brief Trapping int -> uint8_t cast for a face's side count.
 * @param count Face side count to narrow.
 * @return The side count as a uint8_t.
 * @details The Conway/Hankin operators accumulate a face's vertex count in an
 * int, then push it into the uint8_t face_counts array. A face wider than 255
 * sides — reachable as a high-valence vertex orbit (orbit_count), a doubled
 * edge count (count * 2 in snub), or a long Hankin rosette walk (bounded only
 * by 2*I) — would silently wrap, corrupting the face topology. Parallel to
 * narrow_index in intent — route the cast through here so an over-wide value
 * traps at the bench instead of emitting garbage geometry — but a distinct
 * guard: this takes an int and bounds it to [0, UINT8_MAX] (side count), where
 * narrow_index takes a size_t and bounds it to INT16_MAX (vertex index). Cold
 * path: once per emitted face during a rebuild, never per pixel.
 */
inline uint8_t narrow_face_count(int count) {
  HS_CHECK(count >= 0 && count <= UINT8_MAX,
           "mesh face side count exceeds uint8_t range");
  return static_cast<uint8_t>(count);
}

/**
 * @brief Compiles a PolyMesh into a static MeshState.
 * @param src Source dynamic mesh to compile.
 * @param dst Destination MeshState, cleared and populated in place.
 * @param geom_arena Arena supplying storage for the destination arrays.
 * @details Removes degenerate faces (faces with < 3 vertices) during the
 * process. Traps if the cumulative face offset exceeds the 16-bit range.
 */
FLASHMEM static inline void compile(const PolyMesh &src, MeshState &dst,
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

  copy_vector(dst.vertices, src.vertices.data(), src.vertices.size(),
              geom_arena);

  dst.face_counts.bind(geom_arena, valid_faces);
  dst.faces.bind(geom_arena, valid_indices);
  dst.face_offsets.bind(geom_arena, valid_faces);

  size_t offset = 0;
  int current_offset = 0;
  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    int count = src.face_counts[i];
    if (count >= 3) {
      dst.face_counts.push_back(narrow_face_count(count));
      // face_offsets is uint16_t; current_offset is the cumulative index count
      // and wraps silently past 65535. Trap instead, matching build_from_flat's
      // total_indices <= UINT16_MAX gate (not narrow_index's INT16_MAX bound:
      // offsets count half-edges and legitimately exceed the vertex range).
      HS_CHECK(current_offset <= UINT16_MAX,
               "mesh face_offsets exceeds 16-bit index range");
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
 * @tparam MeshT Mesh type (PolyMesh or MeshState).
 * @param src Source mesh to copy from.
 * @param dst Destination mesh, populated in place from the given arena.
 * @param arena Arena supplying storage for the destination arrays.
 * @details Safe for memory compaction and bouncing between arenas. For
 * MeshState this delegates to MeshState::clone (the Cloneable interface used by
 * Persist) so the two implementations can't drift; the generic body below
 * handles PolyMesh (which has no face_offsets).
 */
template <typename MeshT>
inline void clone(const MeshT &src, MeshT &dst, Arena &arena) {
  if constexpr (std::is_same_v<MeshT, MeshState>) {
    MeshState::clone(src, dst, arena);
    return;
  }
  copy_vector(dst.vertices, src.vertices.data(), src.vertices.size(), arena);
  copy_vector(dst.face_counts, src.get_face_counts_data(),
              src.get_face_counts_size(), arena);
  copy_vector(dst.faces, src.get_faces_data(), src.get_faces_size(), arena);

  if constexpr (requires { dst.face_offsets; }) {
    copy_vector(dst.face_offsets, src.get_face_offsets_data(),
                src.get_face_offsets_size(), arena);
  }

  if constexpr (requires { dst.topology; }) {
    copy_vector(dst.topology, src.topology.data(), src.topology.size(), arena);
  }
}

/**
 * @brief MurmurHash3 32-bit finalizer: avalanches the bits of an accumulated
 * hash so small input differences spread across the whole word.
 * @param h Accumulated 32-bit hash value to finalize.
 * @return The avalanched 32-bit hash.
 */
static inline uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

/**
 * @brief Mixes a value into a running hash seed.
 * @param seed Running hash seed, updated in place.
 * @param v Value to fold into the seed.
 */
static inline void hash_combine(uint32_t &seed, uint32_t v) {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief Colors faces based on their vertex count and neighbor topology.
 * @tparam MeshT Mesh type exposing the unified accessors and a topology array.
 * @param mesh Mesh to classify; its topology array is filled with class ids.
 * @param scratch_a Scratch arena for hashes, per-face data, and the node sort.
 * @param scratch_b Scratch arena for the half-edge pairing records.
 * @param persistent Arena backing the mesh topology array when it must grow.
 * @details Hashes each face by side count and interior angles, folds in
 * neighbor hashes via the half-edge pairing, then sorts and assigns a dense
 * topology id per distinct hash. Traps if the half-edge or face count exceeds
 * the 16-bit index range.
 */
template <typename MeshT>
FLASHMEM __attribute__((noinline)) static void
classify_faces_by_topology(MeshT &mesh, Arena &scratch_a, Arena &scratch_b,
                           Arena &persistent) {
  ScratchScope _(scratch_a);

  // Read topology through the unified accessors: a borrowed-mode MeshState
  // serves face_counts/faces via its view spans with the owned vectors empty,
  // so indexing the owned vectors directly would silently classify nothing.
  // Vertices and topology are always owned, so they stay direct.
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  const uint8_t *face_counts = mesh.get_face_counts_data();
  const uint16_t *faces = mesh.get_faces_data();

  ArenaVector<uint32_t> face_hashes;
  face_hashes.bind(scratch_a, F);

  ArenaVector<uint32_t> final_hashes;
  final_hashes.bind(scratch_a, F);

  // Find max face vertex count for scratch allocation
  int max_count = 0;
  for (size_t i = 0; i < F; ++i) {
    int c = face_counts[i];
    if (c > max_count) max_count = c;
  }

  ArenaVector<Vector> verts;
  verts.bind(scratch_a, max_count);
  ArenaVector<int> angles;
  angles.bind(scratch_a, max_count);

  size_t offset = 0;
  for (size_t i = 0; i < F; ++i) {
    int count = face_counts[i];

    verts.clear();
    for (int k = 0; k < count; ++k) {
      verts.push_back(mesh.vertices[faces[offset + k]]);
    }

    // The 32-bit hash (count + sorted rounded interior angles, mixed by fmix32)
    // is used directly as the topology id, so two genuinely different face
    // topologies that collide to the same hash would be merged into one class.
    // No secondary tiebreak is kept: the inputs are the small fixed polyhedron
    // roster (a handful of distinct face shapes per solid), not adversarial, so
    // a 32-bit collision is not reachable in practice. Note the angles below are
    // rounded to whole degrees, which widens the collision target for near-but-
    // distinct angle sets — so the safety margin rests on the fixed roster, not
    // on hash strength.
    //
    // The one untrusted caller is the WASM/JS mesh editor (wasm.cpp's
    // Mesh::classifyFaces over an editor-built mesh). A collision there is still
    // not a memory-safety or invariant breach — it only mis-merges two face
    // *palette* classes (a cosmetic grouping in the tool), and the birthday
    // bound needs ~77k DISTINCT face topologies before a collision is even
    // likely, which no displayable mesh approaches. If this is ever fed meshes
    // dense enough in distinct face classes to matter, add a vertex-key tiebreak
    // here and hash the angles at finer than 1-degree resolution, gated to the
    // untrusted path so the firmware roster keeps its current cheap hashing.
    uint32_t h = 0x12345678;
    hash_combine(h, static_cast<uint32_t>(count));

    // Interior angles only exist for a real polygon. A degenerate face (1 or 2
    // vertices) leaves `angles` empty, so the angle hash-combine must stay
    // inside this guard. Folding angles[k] in unconditionally would read past
    // the empty scratch vector and yield a non-deterministic topology id.
    // Degenerate faces hash on vertex count alone.
    if (count >= 3) {
      angles.clear();
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
      for (int k = 0; k < count; ++k) {
        hash_combine(h, static_cast<uint32_t>(angles[k]));
      }
    }
    h = fmix32(h);
    face_hashes.push_back(h);
    final_hashes.push_back(h);
    offset += count;
  }

  {
    ScratchScope temp_topo(scratch_a);

    // Half-edge ids (he_idx, < I) and face ids (fi, < F) are stored as uint16_t
    // below; an over-range mesh would silently truncate. Trap up front, the same
    // discipline as HalfEdgeMesh::build_from_flat.
    HS_CHECK(I <= UINT16_MAX && F <= UINT16_MAX,
             "classify_faces_by_topology exceeds 16-bit index range");

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
        int count = face_counts[fi];
        for (int k = 0; k < count; ++k) {
          uint16_t u = faces[face_offset + k];
          uint16_t v = faces[face_offset + (k + 1) % count];
          fill_edge_record(records[he_idx], u, v, static_cast<uint16_t>(he_idx));
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
      int count = face_counts[fi];
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

  /**
   * @brief Sortable pairing of a face's final topology hash with its index.
   */
  struct HashNode {
    uint32_t hash;    /**< Final topology hash for the face. */
    int original_face; /**< Index of the face in the original mesh order. */
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
