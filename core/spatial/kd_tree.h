/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file kd_tree.h
 * @brief KDTree and its nearest-point query structures.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include "math/3dmath.h"
#include <algorithm>
#include <span>

#include <cfloat>
#include "containers/static_circular_buffer.h"
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
   * @brief Whether the tree holds no points.
   * @return True for a default-constructed tree or one built from an empty span.
   */
  bool empty() const { return root_index == -1; }

  /**
   * @brief Finds the k nearest neighbors of target, sorted closest-first.
   * @param target Query point, in world units.
   * @param k Number of neighbors to return; MUST be <= MAX_K (traps via HS_CHECK
   *        otherwise). Soft-capped at the point count.
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
    // A k above the point count can never fill the set, so worst_d_sq would
    // stay FLT_MAX and neither the candidate test nor the subtree prune would
    // ever reject: the query degrades to a full traversal.
    k = std::min(k, nodes.size());

    // Cached pruning bound: the largest squared distance in `result` and its
    // slot, FLT_MAX until the set fills to k so nothing prunes early.
    // Recomputed only when `result` changes, keeping the per-node prune test O(1).
    float worst_d_sq = FLT_MAX;
    size_t worst_i = 0;
    auto recompute_worst = [&]() {
      worst_d_sq = -1.0f;
      worst_i = 0;
      for (size_t i = 0; i < result.size(); ++i) {
        if (result[i].d_sq > worst_d_sq ||
            (result[i].d_sq == worst_d_sq &&
             result[i].original_index > result[worst_i].original_index)) {
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
      } else if (d_sq < worst_d_sq ||
                 (d_sq == worst_d_sq &&
                  nodes[idx].original_index < result[worst_i].original_index)) {
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
  ArenaVector<KDNode>
      nodes;           /**< Arena-backed node pool, one entry per point. */
  int root_index = -1; /**< Index of the root node, or -1 if empty. */

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

    if (d_sq <= get_worst_dist())
      offer_candidate(d_sq, node_idx);

    float axis_dist = (node.axis == 0)   ? (target.x - node.point.x)
                      : (node.axis == 1) ? (target.y - node.point.y)
                                         : (target.z - node.point.z);

    int near = axis_dist < 0 ? node.left : node.right;
    int far = axis_dist < 0 ? node.right : node.left;

    search_k(near, target, offer_candidate, get_worst_dist);

    // Re-query the worst bound: the near subtree may have filled the k-best set
    // and tightened it, so the far side can be pruned.
    if ((axis_dist * axis_dist) <= get_worst_dist()) {
      search_k(far, target, offer_candidate, get_worst_dist);
    }
  }
};
