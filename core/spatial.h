/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include "3dmath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <span>

#include <cfloat>
#include "static_circular_buffer.h"
#include "memory.h"

/**
 * @brief Axis-Aligned Bounding Box.
 * @details Speculative: the implementation (union_with, ray-slab intersect) is
 * sound and exercised by test_spatial.h, but it has no production consumer yet —
 * nothing in the engine constructs or queries an AABB. Wire it to a broad-phase
 * / culling consumer before relying on it as core spatial machinery.
 */
struct AABB {
  Vector min_val; /**< Minimum corner (per-axis lower bound). */
  Vector max_val; /**< Maximum corner (per-axis upper bound). */

  /**
   * @brief Constructs an inverted/empty box (min > max) so the first expand() seeds it.
   */
  AABB()
      : min_val(FLT_MAX, FLT_MAX, FLT_MAX),
        max_val(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

  /**
   * @brief Grows this box to also enclose point p.
   * @param p Point to include, in world units.
   */
  void expand(const Vector &p) {
    if (p.x < min_val.x)
      min_val.x = p.x;
    if (p.y < min_val.y)
      min_val.y = p.y;
    if (p.z < min_val.z)
      min_val.z = p.z;

    if (p.x > max_val.x)
      max_val.x = p.x;
    if (p.y > max_val.y)
      max_val.y = p.y;
    if (p.z > max_val.z)
      max_val.z = p.z;
  }

  /**
   * @brief Grows this box to also enclose another box.
   * @param box Box to merge into this one.
   */
  void union_with(const AABB &box) {
    if (box.min_val.x < min_val.x)
      min_val.x = box.min_val.x;
    if (box.min_val.y < min_val.y)
      min_val.y = box.min_val.y;
    if (box.min_val.z < min_val.z)
      min_val.z = box.min_val.z;

    if (box.max_val.x > max_val.x)
      max_val.x = box.max_val.x;
    if (box.max_val.y > max_val.y)
      max_val.y = box.max_val.y;
    if (box.max_val.z > max_val.z)
      max_val.z = box.max_val.z;
  }

  /**
   * @brief Tests whether the ray from origin along direction hits this box.
   * @param origin Ray start point, in world units.
   * @param direction Ray direction, in world units; need not be normalized.
   * @return True if the ray intersects the box; false for an inverted/empty box
   * or a box entirely behind the origin.
   */
  bool intersect_ray(const Vector &origin, const Vector &direction) const {
    float tmin = -FLT_MAX;
    float tmax = FLT_MAX;

    // Slab test, one axis at a time. A zero direction component means the ray
    // runs parallel to that pair of slab planes: it can only intersect if the
    // origin already lies within the slab, otherwise it misses. Guarding this
    // avoids (plane - origin) / 0, which yields 0/0 = NaN when the origin sits
    // exactly on a slab plane and would silently break every subsequent
    // comparison.
    auto slab = [&](float mn, float mx, float o, float d) -> bool {
      if (d != 0.0f) {
        float inv = 1.0f / d;
        float t1 = (mn - o) * inv;
        float t2 = (mx - o) * inv;
        if (t1 > t2)
          std::swap(t1, t2);
        if (t1 > tmin)
          tmin = t1;
        if (t2 < tmax)
          tmax = t2;
        return tmin <= tmax;
      }
      return o >= mn && o <= mx; // parallel ray: inside the slab?
    };

    if (!slab(min_val.x, max_val.x, origin.x, direction.x))
      return false;
    if (!slab(min_val.y, max_val.y, origin.y, direction.y))
      return false;
    if (!slab(min_val.z, max_val.z, origin.z, direction.z))
      return false;

    // Reject boxes entirely behind the ray origin.
    return tmax >= 0.0f;
  }
};

/**
 * @brief A single node of the KDTree.
 */
struct KDNode {
  /**
   * @brief Copy of the split point (not an index into the source array) to
   * avoid lifetime dependence on that array.
   */
  Vector point;
  uint16_t original_index = 0; /**< Index of this point in the source array. */
  int16_t axis = 0;            /**< Splitting axis: 0=x, 1=y, 2=z. */
  int16_t left = -1;           /**< Left child node index, or -1 if none. */
  int16_t right = -1;          /**< Right child node index, or -1 if none. */
};

/**
 * @brief k-d Tree over 3D points using arena storage; supports k-nearest
 * neighbor search.
 */
class KDTree {
public:
  static constexpr int MAX_K = 5; /**< Maximum number of neighbors a query may request. */

  ArenaVector<KDNode> nodes; /**< Arena-backed node pool, one entry per point. */
  int root_index = -1;       /**< Index of the root node, or -1 if empty. */

  /**
   * @brief Constructs an empty tree with no nodes.
   */
  KDTree() = default;

  /**
   * @brief Builds the tree from a span of points using arena storage.
   * @param arena Arena used for node storage and temporary index sorting.
   * @param points Source points; std::span<const Vector> so callers need not
   * const_cast read-only vertex arrays.
   * @details Allocates one node per point in the arena. The scratch index array
   * is scoped to a ScratchScope so its arena offset rewinds once build() returns.
   */
  KDTree(Arena &arena, std::span<const Vector> points) {
    clear();
    if (points.empty())
      return;

    size_t count = points.size();
    nodes.bind(arena, count);

    // The per-build index array is scratch: build() partitions it in place and
    // never reads it after returning, and `nodes` is already fully reserved
    // above. Scope it to a ScratchScope so the arena offset rewinds past it once
    // build() returns — otherwise ~count*4 bytes (≈34 KB at MAX_VERTS) leak into
    // the arena for the rest of its lifetime.
    ScratchScope scratch(arena);
    // count is narrowed to int here (indices is int*, build() takes int count,
    // and indices[i]=i folds the source index into int). The per-node build()
    // guards bound the *stored* index range (uint16_t original_index), which
    // requires every source index 0..count-1 to fit uint16_t; trap that bound
    // at the seam so an over-large point set fails loud here instead of deep in
    // the recursion. Cold path (once per build).
    HS_CHECK(count <= static_cast<size_t>(UINT16_MAX) + 1,
             "KDTree source point count exceeds uint16_t index range");
    int *indices = (int *)arena.allocate(count * sizeof(int), alignof(int));
    for (size_t i = 0; i < count; ++i)
      indices[i] = (int)i;

    root_index = build(points, indices, count, 0);
  }

  /**
   * @brief Resets the tree to empty by dropping the root reference.
   */
  void clear() { root_index = -1; }

  /**
   * @brief Finds the k nearest neighbors of target, sorted closest-first.
   * @param target Query point, in world units.
   * @param k Number of neighbors to return; capped at MAX_K and at the point count.
   * @return Buffer of whole nodes (index + point), closest first.
   */
  StaticCircularBuffer<KDNode, MAX_K> nearest(const Vector &target,
                                              size_t k = 1) const {
    StaticCircularBuffer<KDNode, MAX_K> result;
    if (root_index == -1 || k == 0) // k is size_t; only k == 0 is the empty case
      return result;

    // The result/candidate buffers are sized MAX_K, so more than MAX_K neighbors
    // cannot be honored. Silently returning fewer than requested would mask the
    // caller's mistake; trap so a k > MAX_K bug surfaces instead. (Requesting
    // more neighbors than points exist is fine — that legitimately caps below.)
    HS_CHECK(k <= static_cast<size_t>(MAX_K));

    // Bounded scratch buffer of the k best candidates so far, kept unsorted.
    // Each entry pairs a squared distance with its node index.
    struct Candidate {
      float d_sq;
      int idx;
    };
    // Named `best`, not `heap`: this is a deliberately unsorted linear array of
    // the k best candidates, not a heap-ordered structure. With small K (1-5) a
    // linear scan for the max is cheaper than maintaining heap order, so there
    // is no O(log k) push / O(1) peek-max invariant to rely on.
    StaticCircularBuffer<Candidate, MAX_K> best;

    // Offer a candidate to the k-best set: append while under k, otherwise
    // displace the current worst entry if this one is closer.
    auto offer_candidate = [&](float d_sq, int idx) {
      if (best.size() < static_cast<size_t>(k)) {
        best.push_back({d_sq, idx});
      } else {
        // Full: replace the current worst candidate only if d_sq beats it.
        float max_d = -1.0f;
        int max_i = -1;
        for (size_t i = 0; i < best.size(); ++i) {
          if (best[i].d_sq > max_d) {
            max_d = best[i].d_sq;
            max_i = i;
          }
        }
        if (d_sq < max_d) {
          best[max_i] = {d_sq, idx};
        }
      }
    };

    // Pruning bound: largest squared distance currently held. Returns FLT_MAX
    // until the set holds k entries so nothing is pruned before it fills.
    auto get_worst_dist = [&]() -> float {
      if (best.size() < k)
        return FLT_MAX;
      float max_d = -1.0f;
      for (size_t i = 0; i < best.size(); ++i) {
        if (best[i].d_sq > max_d)
          max_d = best[i].d_sq;
      }
      return max_d;
    };

    search_k(root_index, target, offer_candidate, get_worst_dist);

    std::sort(
        best.begin(), best.end(),
        [](const Candidate &a, const Candidate &b) { return a.d_sq < b.d_sq; });

    for (const auto &c : best) {
      result.push_back(nodes[c.idx]);
    }
    return result;
  }

private:
  /**
   * @brief Recursively builds a balanced subtree over indices[0..count).
   * @param points Source points indexed by entries of `indices`.
   * @param indices Scratch index array, partitioned in place by this call.
   * @param count Number of indices in this subtree.
   * @param depth Recursion depth; selects the split axis as depth % 3.
   * @return Root node index of the built subtree, or -1 if empty.
   * @details Cycles the split axis by depth%3, partitioning around the median
   * along that axis and reordering `indices` in place.
   */
  int build(std::span<const Vector> points, int *indices, int count, int depth) {
    if (count <= 0)
      return -1; // legitimate empty-subtree sentinel (leaf recursion base case)
    // The node pool is bound to exactly points.size() and build() emits exactly
    // one node per point, so this can never fire. Trap rather than return -1:
    // a -1 here is indistinguishable from the empty-subtree sentinel above, so
    // pool exhaustion would silently drop an entire subtree instead of failing
    // loud. Cold path (once per node during build).
    HS_CHECK(nodes.size() < nodes.capacity(),
             "KDTree node pool exhausted during build");

    int axis = depth % 3;
    int mid = count / 2;

    auto *start = indices;
    auto *end = indices + count;

    std::nth_element(start, start + mid, end, [&](int a, int b) {
      float va = (axis == 0)   ? points[a].x
                 : (axis == 1) ? points[a].y
                               : points[a].z;
      float vb = (axis == 0)   ? points[b].x
                 : (axis == 1) ? points[b].y
                               : points[b].z;
      return va < vb;
    });

    int median_idx = indices[mid];

    int new_node_idx = static_cast<int>(nodes.size());
    // Child links (left/right) are stored as int16_t with -1 as the sentinel,
    // so a node pool bumped past INT16_MAX would silently truncate them into
    // garbage links. The pool-capacity guard above caps nodes.size(), but trap
    // the index range explicitly so a future capacity bump fails loud instead of
    // corrupting the tree. Cold path (once per node during build).
    HS_CHECK(new_node_idx <= INT16_MAX,
             "KDTree node index exceeds int16_t child-link range");
    // original_index is uint16_t, so a source set larger than 65535 points
    // would silently fold a high index back into the wrong vertex. The node
    // guard above bounds nodes.size(), not the source index, so trap the vertex
    // range too. Cold path (once per node during build).
    HS_CHECK(median_idx <= UINT16_MAX,
             "KDTree vertex index exceeds uint16_t original_index range");
    nodes.emplace_back();
    nodes[new_node_idx].point = points[median_idx];
    nodes[new_node_idx].original_index = (uint16_t)median_idx;
    nodes[new_node_idx].axis = (int16_t)axis;

    nodes[new_node_idx].left = (int16_t)build(points, start, mid, depth + 1);
    nodes[new_node_idx].right =
        (int16_t)build(points, start + mid + 1, count - mid - 1, depth + 1);

    return new_node_idx;
  }

  /**
   * @brief Recursive k-NN traversal of the subtree rooted at node_idx.
   * @tparam PushFn Callable offering a (d_sq, idx) candidate to the k-best set.
   * @tparam MaxDistFn Callable returning the current worst squared distance bound.
   * @param node_idx Subtree root node index, or -1 to terminate.
   * @param target Query point, in world units.
   * @param offer_candidate Callback that records a candidate in the k-best set.
   * @param get_worst_dist Callback returning the current pruning bound (squared distance).
   * @details Descends the near child first, then prunes the far child when the
   * splitting plane is farther than the current worst hit. The k-best set and k
   * itself are reached only through the callbacks, which capture them by
   * reference in nearest().
   */
  template <typename PushFn, typename MaxDistFn>
  void search_k(int node_idx, const Vector &target, PushFn &&offer_candidate,
                MaxDistFn &&get_worst_dist) const {
    if (node_idx == -1)
      return;

    const KDNode &node = nodes[node_idx];
    float d_sq = distance_squared(node.point, target);

    if (d_sq < get_worst_dist())
      offer_candidate(d_sq, node_idx);

    float axis_dist = (node.axis == 0)   ? (target.x - node.point.x)
                     : (node.axis == 1) ? (target.y - node.point.y)
                                        : (target.z - node.point.z);

    int near = axis_dist < 0 ? node.left : node.right;
    int far = axis_dist < 0 ? node.right : node.left;

    search_k(near, target, offer_candidate, get_worst_dist);

    // Re-query the worst bound: the near subtree may have filled the k-best set
    // and tightened it, so the far side can be pruned.
    if ((axis_dist * axis_dist) < get_worst_dist()) {
      search_k(far, target, offer_candidate, get_worst_dist);
    }
  }
};

/**
 * @brief Represents the state of a mesh using arena storage to avoid heap
 * allocations.
 */
struct MeshState {
  ArenaVector<Vector> vertices;        /**< Owned vertex positions, in world units. */
  ArenaVector<uint8_t> face_counts;    /**< Owned per-face vertex counts. */
  ArenaVector<uint16_t> faces;         /**< Owned flattened face vertex indices. */
  ArenaVector<uint16_t> face_offsets;  /**< Owned start offset of each face into faces. */
  ArenaVector<int> topology;           /**< Owned adjacency/topology data. */

  ArenaSpan<uint8_t> face_counts_view;   /**< Borrowed face-counts view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t> faces_view;        /**< Borrowed faces view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t> face_offsets_view; /**< Borrowed face-offsets view, populated by MeshOps::transform. */

  /**
   * @brief Constructs an empty, unbound mesh.
   */
  MeshState() = default;

  /**
   * @brief Move-constructs, transferring owned buffers and views.
   * @param other Source mesh; its borrowed views are cleared so it holds no
   * dangling borrows.
   */
  MeshState(MeshState &&other) noexcept
      : vertices(std::move(other.vertices)),
        face_counts(std::move(other.face_counts)),
        faces(std::move(other.faces)),
        face_offsets(std::move(other.face_offsets)),
        topology(std::move(other.topology)),
        face_counts_view(other.face_counts_view),
        faces_view(other.faces_view),
        face_offsets_view(other.face_offsets_view) {
    other.face_counts_view = {};
    other.faces_view = {};
    other.face_offsets_view = {};
  }

  /**
   * @brief Move-assigns, transferring owned buffers and views.
   * @param other Source mesh; its borrowed views are cleared so the moved-from
   * mesh holds no dangling borrows.
   * @return Reference to this mesh.
   */
  MeshState &operator=(MeshState &&other) noexcept {
    if (this != &other) {
      vertices = std::move(other.vertices);
      face_counts = std::move(other.face_counts);
      faces = std::move(other.faces);
      face_offsets = std::move(other.face_offsets);
      topology = std::move(other.topology);
      face_counts_view = other.face_counts_view;
      faces_view = other.faces_view;
      face_offsets_view = other.face_offsets_view;
      other.face_counts_view = {};
      other.faces_view = {};
      other.face_offsets_view = {};
    }
    return *this;
  }

  /**
   * @brief Resets to empty: clears owned buffers and drops the borrowed views.
   */
  void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    face_offsets.clear();
    topology.clear();
    face_counts_view = {};
    faces_view = {};
    face_offsets_view = {};
  }

  /**
   * @brief Checks whether any member vector is bound (has been allocated).
   * @return True if the mesh owns allocated storage.
   */
  bool is_bound() const { return vertices.is_bound(); }

  /**
   * @brief Returns the face-counts pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   * @details Discriminates on is_bound() (which mode this is), NOT on empty():
   * an owned-but-legitimately-empty mesh is bound with size 0, and gating on
   * empty() would wrongly fall through to a stale/unset borrowed view. The same
   * rationale applies to the sibling accessors below.
   */
  const uint8_t *get_face_counts_data() const {
    return face_counts.is_bound() ? face_counts.data() : face_counts_view.data();
  }
  /**
   * @brief Returns the number of face counts for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_face_counts_size() const {
    return face_counts.is_bound() ? face_counts.size() : face_counts_view.size();
  }

  /**
   * @brief Returns the faces pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   */
  const uint16_t *get_faces_data() const {
    return faces.is_bound() ? faces.data() : faces_view.data();
  }
  /**
   * @brief Returns the number of face indices for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_faces_size() const {
    return faces.is_bound() ? faces.size() : faces_view.size();
  }

  /**
   * @brief Returns the face-offsets pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   */
  const uint16_t *get_face_offsets_data() const {
    return face_offsets.is_bound() ? face_offsets.data()
                                   : face_offsets_view.data();
  }
  /**
   * @brief Returns the number of face offsets for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_face_offsets_size() const {
    return face_offsets.is_bound() ? face_offsets.size()
                                   : face_offsets_view.size();
  }

  /**
   * @brief Returns the number of vertices in the mesh.
   * @return Vertex count.
   */
  size_t num_vertices() const { return vertices.size(); }
  /**
   * @brief Returns the number of faces in the mesh.
   * @return Face count (equal to the active face-counts size).
   */
  size_t num_faces() const { return get_face_counts_size(); }

  /**
   * @brief Deep-copies all owned data from src into dst using a target arena.
   * @param src Source mesh to copy from.
   * @param dst Destination mesh to populate.
   * @param arena Arena providing storage for the destination buffers.
   * @details Required by Cloneable.
   */
  static void clone(const MeshState &src, MeshState &dst, Arena &arena) {
    dst.vertices.bind(arena, src.vertices.size());
    for (size_t i = 0; i < src.vertices.size(); ++i)
      dst.vertices.push_back(src.vertices[i]);

    size_t fc_size = src.get_face_counts_size();
    const uint8_t *fc_data = src.get_face_counts_data();
    dst.face_counts.bind(arena, fc_size);
    for (size_t i = 0; i < fc_size; ++i)
      dst.face_counts.push_back(fc_data[i]);

    size_t f_size = src.get_faces_size();
    const uint16_t *f_data = src.get_faces_data();
    dst.faces.bind(arena, f_size);
    for (size_t i = 0; i < f_size; ++i)
      dst.faces.push_back(f_data[i]);

    size_t fo_size = src.get_face_offsets_size();
    if (fo_size > 0) {
      const uint16_t *fo_data = src.get_face_offsets_data();
      dst.face_offsets.bind(arena, fo_size);
      for (size_t i = 0; i < fo_size; ++i)
        dst.face_offsets.push_back(fo_data[i]);
    }

    if (!src.topology.is_empty()) {
      dst.topology.bind(arena, src.topology.size());
      for (size_t i = 0; i < src.topology.size(); ++i)
        dst.topology.push_back(src.topology[i]);
    }
  }
};
