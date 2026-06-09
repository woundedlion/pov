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
 */
struct AABB {
  Vector min_val;
  Vector max_val;

  AABB()
      : min_val(FLT_MAX, FLT_MAX, FLT_MAX),
        max_val(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

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

  // Grow this box to also enclose box.
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

  bool intersect_ray(const Vector &origin, const Vector &direction) const {
    float tmin = -FLT_MAX;
    float tmax = FLT_MAX;

    // Slab test, one axis at a time. A zero direction component means the ray
    // runs parallel to that pair of slab planes: it can only intersect if the
    // origin already lies within the slab, otherwise it misses. Guarding this
    // avoids the naive (plane - origin) / 0 — which is ±inf for a glancing ray
    // (the old comparisons happened to tolerate that) but 0/0 = NaN when the
    // origin sits exactly on a slab plane, silently breaking every subsequent
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
  // Stores a copy of the point (not an index into the source array) to avoid
  // lifetime dependence on that array.
  Vector point;
  uint16_t original_index = 0; // Index in the source array
  int16_t axis = 0;           // 0=x, 1=y, 2=z
  int16_t left = -1;
  int16_t right = -1;
};

/**
 * @brief k-d Tree over 3D points using arena storage; supports k-nearest
 * neighbor search.
 */
class KDTree {
public:
  static constexpr int MAX_K = 5;

  ArenaVector<KDNode> nodes;
  int root_index = -1;

  KDTree() = default;

  // Build from a Span of vectors, using arena for node storage and temporary
  // index sorting. Accepts std::span<const Vector> so callers don't have to
  // const_cast read-only vertex arrays.
  KDTree(Arena &arena, std::span<const Vector> points) {
    clear();
    if (points.empty())
      return;

    size_t count = points.size();
    nodes.bind(arena, count);

    // Temporary index allocation for building the tree
    int *indices = (int *)arena.allocate(count * sizeof(int), alignof(int));
    for (size_t i = 0; i < count; ++i)
      indices[i] = i;

    root_index = build(points, indices, count, 0);
  }

  void clear() { root_index = -1; }

  // Nearest neighbor search
  // Returns up to K results (Nodes, so we get index + point)
  StaticCircularBuffer<KDNode, MAX_K> nearest(const Vector &target,
                                              size_t k = 1) const {
    StaticCircularBuffer<KDNode, MAX_K> result;
    if (root_index == -1 || k <= 0)
      return result;

    // The result/heap buffers are sized MAX_K, so more than MAX_K neighbors
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
    StaticCircularBuffer<Candidate, MAX_K> heap;

    auto push_heap = [&](float d_sq, int idx) {
      if (heap.size() < static_cast<size_t>(k)) {
        // Kept unsorted: with small K (1-5) a linear scan for the max is
        // cheaper than maintaining real heap order.
        heap.push_back({d_sq, idx});
      } else {
        // Replace worst if better
        float max_d = -1.0f;
        int max_i = -1;
        for (size_t i = 0; i < heap.size(); ++i) {
          if (heap[i].d_sq > max_d) {
            max_d = heap[i].d_sq;
            max_i = i;
          }
        }
        if (d_sq < max_d) {
          heap[max_i] = {d_sq, idx};
        }
      }
    };

    auto get_worst_dist = [&]() -> float {
      if (heap.size() < k)
        return FLT_MAX;
      float max_d = -1.0f;
      for (size_t i = 0; i < heap.size(); ++i) {
        if (heap[i].d_sq > max_d)
          max_d = heap[i].d_sq;
      }
      return max_d;
    };

    search_k(root_index, target, k, heap, push_heap, get_worst_dist);

    // Sort result by distance (closest first)
    std::sort(
        heap.begin(), heap.end(),
        [](const Candidate &a, const Candidate &b) { return a.d_sq < b.d_sq; });

    for (const auto &c : heap) {
      result.push_back(nodes[c.idx]);
    }
    return result;
  }

private:
  int build(std::span<const Vector> points, int *indices, int count, int depth) {
    if (count <= 0)
      return -1;
    if (nodes.size() >= nodes.capacity())
      return -1;

    int axis = depth % 3;
    int mid = count / 2;

    // nth_element / partition
    // We need to partition the indices based on the points they point to
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
    nodes[new_node_idx].original_index = (uint16_t)median_idx; // Store index
    nodes[new_node_idx].axis = (int16_t)axis;

    nodes[new_node_idx].left = (int16_t)build(points, start, mid, depth + 1);
    nodes[new_node_idx].right =
        (int16_t)build(points, start + mid + 1, count - mid - 1, depth + 1);

    return new_node_idx;
  }

  // K-Nearest Search
  template <typename PushFn, typename MaxDistFn>
  void search_k(int node_idx, const Vector &target, int k, const auto &heap,
                PushFn &&push_heap, MaxDistFn &&get_worst_dist) const {
    if (node_idx == -1)
      return;

    const KDNode &node = nodes[node_idx];
    float d_sq = distance_squared(node.point, target);
    float worst_sq = get_worst_dist();

    if (d_sq < worst_sq) {
      push_heap(d_sq, node_idx);
      worst_sq = get_worst_dist(); // Update after push
    }

    float axis_dist = (node.axis == 0)   ? (target.x - node.point.x)
                     : (node.axis == 1) ? (target.y - node.point.y)
                                        : (target.z - node.point.z);

    int near = axis_dist < 0 ? node.left : node.right;
    int far = axis_dist < 0 ? node.right : node.left;

    search_k(near, target, k, heap, push_heap, get_worst_dist);

    // Query worst again
    worst_sq = get_worst_dist();
    if ((axis_dist * axis_dist) < worst_sq) {
      search_k(far, target, k, heap, push_heap, get_worst_dist);
    }
  }
};

/**
 * @brief Represents the state of a mesh using arena storage to avoid heap
 * allocations.
 */
struct MeshState {
  ArenaVector<Vector> vertices;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<uint16_t> faces;
  ArenaVector<uint16_t> face_offsets;
  ArenaVector<int> topology;

  // Borrowed (non-owning) views — populated by MeshOps::transform
  ArenaSpan<uint8_t> face_counts_view;
  ArenaSpan<uint16_t> faces_view;
  ArenaSpan<uint16_t> face_offsets_view;

  MeshState() = default;

  // Explicit move semantics to ensure source is invalidated
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

  /// Check if any member vector is bound (has been allocated).
  bool is_bound() const { return vertices.is_bound(); }

  // Unified accessors: return owned data in owned mode, the borrowed view in
  // borrowed mode. Discriminate on is_bound() (which mode this is), NOT on
  // empty(): an owned-but-legitimately-empty mesh is bound with size 0, and
  // gating on empty() would wrongly fall through to a stale/unset borrowed view.
  const uint8_t *get_face_counts_data() const {
    return face_counts.is_bound() ? face_counts.data() : face_counts_view.data();
  }
  size_t get_face_counts_size() const {
    return face_counts.is_bound() ? face_counts.size() : face_counts_view.size();
  }

  const uint16_t *get_faces_data() const {
    return faces.is_bound() ? faces.data() : faces_view.data();
  }
  size_t get_faces_size() const {
    return faces.is_bound() ? faces.size() : faces_view.size();
  }

  const uint16_t *get_face_offsets_data() const {
    return face_offsets.is_bound() ? face_offsets.data()
                                   : face_offsets_view.data();
  }
  size_t get_face_offsets_size() const {
    return face_offsets.is_bound() ? face_offsets.size()
                                   : face_offsets_view.size();
  }

  // Helper for size accessors if needed, but direct vector access is preferred.
  size_t num_vertices() const { return vertices.size(); }
  size_t num_faces() const { return get_face_counts_size(); }

  /// Deep-copy all owned data into a target arena. Required by Cloneable.
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

    if (!src.topology.empty()) {
      dst.topology.bind(arena, src.topology.size());
      for (size_t i = 0; i < src.topology.size(); ++i)
        dst.topology.push_back(src.topology[i]);
    }
  }
};
