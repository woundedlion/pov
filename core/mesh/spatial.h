/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file spatial.h
 * @brief KDTree and the nearest-point query structures the mesh and effect
 *        code build on.
 */

#include <cstdint>
#include "math/3dmath.h"
#include <algorithm>
#include <span>

#include <cfloat>
#include "engine/static_circular_buffer.h"
#include "engine/memory.h"

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
 * @brief A single neighbor returned by KDTree::nearest().
 * @details Carries the squared distance the search already computed; excludes
 * KDNode's internal tree links, which are meaningless outside the tree.
 */
struct Neighbor {
  Vector point;                /**< Copy of the neighbor's position. */
  uint16_t original_index = 0; /**< Index of this point in the source array. */
  float d_sq = 0.0f;           /**< Squared distance from the query point. */
};

/**
 * @brief k-d Tree over 3D points using arena storage; supports k-nearest
 * neighbor search.
 */
class KDTree {
public:
  static constexpr int MAX_K =
      5; /**< Maximum number of neighbors a query may request. */

  /**
   * @brief Upper bound on source points, set by the int16_t node-link range.
   * @details One node per point, and node indices flow through int16_t
   * left/right links, so the point count must fit in [0, MAX_POINTS]. This also
   * bounds original_index (uint16_t) since indices stay below the count.
   */
  static constexpr size_t MAX_POINTS = static_cast<size_t>(INT16_MAX) + 1;

  ArenaVector<KDNode>
      nodes;           /**< Arena-backed node pool, one entry per point. */
  int root_index = -1; /**< Index of the root node, or -1 if empty. */

  /**
   * @brief Constructs an empty tree with no nodes.
   */
  KDTree() = default;

  /**
   * @brief Builds the tree from a span of points using arena storage.
   * @param arena Arena used for node storage and temporary index sorting.
   * @param points Source points; std::span<const Vector> so callers need not
   * const_cast read-only vertex arrays. Every coordinate must be finite (traps
   * via HS_CHECK otherwise).
   * @details Allocates one node per point in the arena. The scratch index array
   * is scoped to a ScratchScope so its arena offset rewinds once build() returns.
   * Retains N * sizeof(KDNode) bytes and peaks a further N * sizeof(int) over
   * that (N points): 8000 B retained over a 1600 B transient at Voronoi's
   * MAX_SITES = 400, and at most 128 KB transient at MAX_POINTS.
   */
  HS_COLD_MEMBER KDTree(Arena &arena, std::span<const Vector> points) {
    if (points.empty())
      return;

    size_t count = points.size();
    HS_CHECK(
        count <= MAX_POINTS,
        "KDTree source point count exceeds int16_t child-link index range");
    nodes.bind(arena, count);

    // Scope the scratch index array so the arena offset rewinds once build()
    // returns.
    ScratchScope scratch(arena);
    int *indices = arena.allocate_n<int>(count);
    for (size_t i = 0; i < count; ++i) {
      // A NaN coordinate compares equivalent to every value, which breaks the
      // strict weak ordering build()'s comparator owes std::nth_element.
      HS_CHECK(std::isfinite(points[i].x) && std::isfinite(points[i].y) &&
                   std::isfinite(points[i].z),
               "KDTree source point is not finite");
      indices[i] = (int)i;
    }

    root_index = build(points, indices, count, 0);
  }

  /**
   * @brief Finds the k nearest neighbors of target, sorted closest-first.
   * @param target Query point, in world units.
   * @param k Number of neighbors to return; MUST be <= MAX_K (traps via HS_CHECK
   *        otherwise). Soft-capped only at the point count.
   * @return Buffer of neighbors (point + source index + squared distance), closest first.
   */
  StaticCircularBuffer<Neighbor, MAX_K> nearest(const Vector &target,
                                                size_t k = 1) const {
    StaticCircularBuffer<Neighbor, MAX_K> result;
    HS_CHECK(k <= static_cast<size_t>(MAX_K),
             "KDTree::nearest k exceeds MAX_K");
    if (root_index == -1 ||
        k == 0) // k is size_t; only k == 0 is the empty case
      return result;

    // Cached pruning bound: the largest squared distance in `result` and its
    // slot, FLT_MAX until the set fills to k so nothing prunes early.
    // Recomputed only when `result` changes, keeping the per-node prune test O(1).
    float worst_d_sq = FLT_MAX;
    size_t worst_i = 0;
    auto recompute_worst = [&]() {
      worst_d_sq = -1.0f;
      worst_i = 0;
      for (size_t i = 0; i < result.size(); ++i) {
        if (result[i].d_sq > worst_d_sq) {
          worst_d_sq = result[i].d_sq;
          worst_i = i;
        }
      }
    };

    // Offer a node to the k-best set: append while under k, otherwise displace
    // the cached worst entry if this one is closer.
    auto offer_candidate = [&](float d_sq, int idx) {
      if (result.size() < static_cast<size_t>(k)) {
        result.push_back({nodes[idx].point, nodes[idx].original_index, d_sq});
        if (result.size() == static_cast<size_t>(k))
          recompute_worst(); // set just filled: cache its worst for pruning
      } else if (d_sq < worst_d_sq) {
        result[worst_i] = {nodes[idx].point, nodes[idx].original_index, d_sq};
        recompute_worst(); // worst displaced: refresh the cache
      }
    };

    auto get_worst_dist = [&]() -> float { return worst_d_sq; };

    search_k(root_index, target, offer_candidate, get_worst_dist);

    // Total order (distance, then source index): std::sort is unstable, so
    // equidistant neighbors would otherwise come back in an unspecified order.
    std::sort(result.begin(), result.end(),
              [](const Neighbor &a, const Neighbor &b) {
                return a.d_sq != b.d_sq ? a.d_sq < b.d_sq
                                        : a.original_index < b.original_index;
              });
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
  HS_COLD_MEMBER int build(std::span<const Vector> points, int *indices,
                           int count, int depth) {
    if (count <= 0)
      return -1; // legitimate empty-subtree sentinel (leaf recursion base case)
    // Trap rather than return -1: a -1 here is indistinguishable from the
    // empty-subtree sentinel above, so exhaustion would silently drop a subtree.
    HS_CHECK(nodes.size() < nodes.capacity(),
             "KDTree node pool exhausted during build");

    int axis = depth % 3;
    int mid = count / 2;

    auto *start = indices;
    auto *end = indices + count;

    // Total order (axis value, then source index): axis ties are the norm for
    // polyhedron vertices, and an axis-only comparator leaves the tree shape up
    // to the standard library's partition order.
    std::nth_element(start, start + mid, end, [&](int a, int b) {
      float va = (axis == 0)   ? points[a].x
                 : (axis == 1) ? points[a].y
                               : points[a].z;
      float vb = (axis == 0)   ? points[b].x
                 : (axis == 1) ? points[b].y
                               : points[b].z;
      return va != vb ? va < vb : a < b;
    });

    int median_idx = indices[mid];

    int new_node_idx = static_cast<int>(nodes.size());
    // left/right are int16_t (-1 sentinel); original_index is uint16_t. Both
    // ranges are covered by MAX_POINTS, which the source count is capped to.
    HS_CHECK(static_cast<size_t>(new_node_idx) < MAX_POINTS,
             "KDTree node index exceeds int16_t child-link range");
    HS_CHECK(static_cast<size_t>(median_idx) < MAX_POINTS,
             "KDTree vertex index exceeds original_index range");
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
 * @brief Deep-copies n contiguous elements into a freshly-bound ArenaVector.
 * @tparam T Element type.
 * @param dst Destination vector, bound to exactly n elements from arena.
 * @param src Pointer to the first source element (ignored when n == 0).
 * @param n Number of elements to copy.
 * @param arena Arena supplying storage for dst.
 * @details Shared by MeshState::clone, MeshOps::clone and
 * CompiledHankin::clone.
 */
template <typename T>
inline void copy_vector(ArenaVector<T> &dst, const T *src, size_t n,
                        Arena &arena) {
  dst.bind(arena, n);
  dst.append_bulk(src, n);
}

/**
 * @brief Represents the state of a mesh using arena storage to avoid heap
 * allocations.
 */
struct MeshState {
  ArenaVector<Vector> vertices; /**< Owned vertex positions, in world units. */
  ArenaVector<uint8_t> face_counts; /**< Owned per-face vertex counts. */
  ArenaVector<uint16_t> faces;      /**< Owned flattened face vertex indices. */
  ArenaVector<uint16_t>
      face_offsets; /**< Owned start offset of each face into faces. */
  ArenaVector<uint16_t> topology; /**< Optional owned per-face topology class
                   ids from classify_faces_by_topology; either empty or one
                   dense 16-bit-bounded id per face. */

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
        face_counts_view(other.face_counts_view), faces_view(other.faces_view),
        face_offsets_view(other.face_offsets_view),
        topology_view(other.topology_view) {
    other.set_owned();
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
      topology_view = other.topology_view;
      other.set_owned();
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
    set_owned();
  }

  /**
   * @brief Checks whether the vertex buffer is bound (has been allocated).
   * @return True if the mesh owns vertex storage; says nothing about the
   *   topology arrays, which are unbound in borrowed mode.
   */
  bool is_bound() const { return vertices.is_bound(); }

  /**
   * @brief Returns the face-counts pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   * @details Discriminates on is_bound() (which mode this is), NOT on empty():
   * an owned-but-legitimately-empty mesh is bound with size 0, and gating on
   * empty() would wrongly fall through to a stale/unset borrowed view (same for
   * the sibling accessors below).
   */
  const uint8_t *get_face_counts_data() const {
    return face_counts.is_bound() ? face_counts.data()
                                  : face_counts_view.data();
  }
  /**
   * @brief Returns the number of face counts for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_face_counts_size() const {
    return face_counts.is_bound() ? face_counts.size()
                                  : face_counts_view.size();
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
   * @brief Returns the per-face topology-class pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   *   Null when the mesh is unclassified.
   */
  const uint16_t *get_topology_data() const {
    return topology.is_bound() ? topology.data() : topology_view.data();
  }
  /**
   * @brief Returns the number of per-face topology classes for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size; 0 when
   *   the mesh is unclassified.
   */
  size_t get_topology_size() const {
    return topology.is_bound() ? topology.size() : topology_view.size();
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
   * @details Required by Cloneable. Traps if src aliases dst: the set_owned()
   * below would drop a borrowed src's views before they are read, yielding an
   * empty mesh.
   */
  HS_COLD_MEMBER static void clone(const MeshState &src, MeshState &dst,
                                   Arena &arena) {
    HS_CHECK(&src != &dst, "MeshState::clone src must not alias dst");
    // Reused dst may carry stale views; clone produces an owned-mode mesh.
    dst.set_owned();
    copy_vector(dst.vertices, src.vertices.data(), src.vertices.size(), arena);
    copy_vector(dst.face_counts, src.get_face_counts_data(),
                src.get_face_counts_size(), arena);
    copy_vector(dst.faces, src.get_faces_data(), src.get_faces_size(), arena);
    copy_vector(dst.face_offsets, src.get_face_offsets_data(),
                src.get_face_offsets_size(), arena);
    copy_vector(dst.topology, src.get_topology_data(), src.get_topology_size(),
                arena);
  }

  /**
   * @brief Switches to owned mode by dropping the borrowed topology views so the
   *   is_bound() accessors read the owned buffers.
   * @details Call on a possibly-borrowed MeshState after (re)binding its owned
   *   topology, so a leftover view can't shadow the owned data.
   */
  void set_owned() {
    face_counts_view = {};
    faces_view = {};
    face_offsets_view = {};
    topology_view = {};
  }

  /**
   * @brief Reports whether every face offset is the running sum of the face
   *   counts before it.
   * @param face_counts_span Per-face vertex counts.
   * @param face_offsets_span Per-face start offsets into the flat faces list;
   *   must be one entry per face count.
   * @return True when offset[i] equals counts[0] + ... + counts[i-1] for every
   *   face.
   * @details Exact equality at each index, which also establishes that the
   *   offsets are non-decreasing and that no two face spans overlap. Endpoint
   *   agreement alone does not: interior offsets can be scrambled while the
   *   first and last still line up.
   */
  static bool offsets_are_prefix_sum(ArenaSpan<uint8_t> face_counts_span,
                                     ArenaSpan<uint16_t> face_offsets_span) {
    size_t running = 0;
    for (size_t i = 0; i < face_counts_span.size(); ++i) {
      if (static_cast<size_t>(face_offsets_span[i]) != running)
        return false;
      running += face_counts_span[i];
    }
    return true;
  }

  /**
   * @brief Switches to borrowed mode: drops the owned face_counts/faces/
   *   face_offsets/topology arrays so they cannot shadow the views, then points
   *   the four views at the given spans.
   * @param face_counts_span Borrowed per-face vertex counts.
   * @param faces_span Borrowed flattened face vertex indices.
   * @param face_offsets_span Borrowed per-face start offsets into faces. Empty
   *   when the source mesh carries no offsets (only the solid scan path needs
   *   them); otherwise one entry per face.
   * @param topology_span Borrowed per-face topology class ids. Empty when the
   *   source mesh is unclassified; otherwise one entry per face.
   * @details Traps on inconsistent spans: a present offsets array must be one
   *   entry per face, and its last offset plus that face's count must cover the
   *   whole flat faces list. With no offsets the counts must sum to the flat
   *   faces length. A present topology array must be one entry per face. The
   *   interior offsets are audited against the prefix sum under HS_AUDIT_CHECK
   *   only: this runs per frame from MeshOps::transform, so the device pays
   *   nothing for the O(F) walk.
   */
  void set_view(ArenaSpan<uint8_t> face_counts_span,
                ArenaSpan<uint16_t> faces_span,
                ArenaSpan<uint16_t> face_offsets_span,
                ArenaSpan<uint16_t> topology_span = {}) {
    if (!face_offsets_span.is_empty()) {
      HS_CHECK(face_offsets_span.size() == face_counts_span.size(),
               "MeshState::set_view: one face offset per face count required");
      const size_t last = face_counts_span.size() - 1;
      HS_CHECK(static_cast<size_t>(face_offsets_span[last]) +
                       face_counts_span[last] ==
                   faces_span.size(),
               "MeshState::set_view: face offsets do not span faces");
      HS_AUDIT_CHECK(
          offsets_are_prefix_sum(face_counts_span, face_offsets_span),
          "MeshState::set_view: face offsets are not the prefix sum of the "
          "face counts");
    } else {
      size_t counted = 0;
      for (size_t i = 0; i < face_counts_span.size(); ++i)
        counted += face_counts_span[i];
      HS_CHECK(counted == faces_span.size(),
               "MeshState::set_view: face counts do not span faces");
    }
    HS_CHECK(topology_span.is_empty() ||
                 topology_span.size() == face_counts_span.size(),
             "MeshState::set_view: one topology class per face required");
    face_counts = {};
    faces = {};
    face_offsets = {};
    topology = {};
    face_counts_view = face_counts_span;
    faces_view = faces_span;
    face_offsets_view = face_offsets_span;
    topology_view = topology_span;
  }

private:
  // Private so set_view() stays the only way to enter borrowed mode: its
  // consistency traps are what downstream fi + fo[f] indexing rests on.
  ArenaSpan<uint8_t>
      face_counts_view; /**< Borrowed face-counts view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      faces_view; /**< Borrowed faces view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      face_offsets_view; /**< Borrowed face-offsets view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      topology_view; /**< Borrowed per-face topology view, populated by MeshOps::transform. */
};
