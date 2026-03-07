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

#include <unordered_map>
#include <memory>
#include <cfloat>
#include "static_circular_buffer.h"
#include "memory.h"

/**
 * @brief Axis-Aligned Bounding Box.
 */
struct AABB {
  Vector minVal;
  Vector maxVal;

  AABB()
      : minVal(FLT_MAX, FLT_MAX, FLT_MAX),
        maxVal(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

  void expand(const Vector &p) {
    if (p.x < minVal.x)
      minVal.x = p.x;
    if (p.y < minVal.y)
      minVal.y = p.y;
    if (p.z < minVal.z)
      minVal.z = p.z;

    if (p.x > maxVal.x)
      maxVal.x = p.x;
    if (p.y > maxVal.y)
      maxVal.y = p.y;
    if (p.z > maxVal.z)
      maxVal.z = p.z;
  }

  // Optimize union
  void unionWith(const AABB &box) {
    if (box.minVal.x < minVal.x)
      minVal.x = box.minVal.x;
    if (box.minVal.y < minVal.y)
      minVal.y = box.minVal.y;
    if (box.minVal.z < minVal.z)
      minVal.z = box.minVal.z;

    if (box.maxVal.x > maxVal.x)
      maxVal.x = box.maxVal.x;
    if (box.maxVal.y > maxVal.y)
      maxVal.y = box.maxVal.y;
    if (box.maxVal.z > maxVal.z)
      maxVal.z = box.maxVal.z;
  }

  bool intersectRay(const Vector &origin, const Vector &direction) const {
    float tmin = (minVal.x - origin.x) / direction.x;
    float tmax = (maxVal.x - origin.x) / direction.x;
    if (tmin > tmax)
      std::swap(tmin, tmax);

    float tymin = (minVal.y - origin.y) / direction.y;
    float tymax = (maxVal.y - origin.y) / direction.y;
    if (tymin > tymax)
      std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
      return false;

    if (tymin > tmin)
      tmin = tymin;
    if (tymax < tmax)
      tmax = tymax;

    float tzmin = (minVal.z - origin.z) / direction.z;
    float tzmax = (maxVal.z - origin.z) / direction.z;
    if (tzmin > tzmax)
      std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
      return false;

    return true;
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

  void invalidate() {
    vertices.invalidate();
    face_counts.invalidate();
    faces.invalidate();
    face_offsets.invalidate();
    topology.invalidate();
    face_counts_view = {};
    faces_view = {};
    face_offsets_view = {};
  }

  // Unified accessors: return whichever is populated (owned or borrowed)
  const uint8_t *get_face_counts_data() const {
    return face_counts.empty() ? face_counts_view.data() : face_counts.data();
  }
  size_t get_face_counts_size() const {
    return face_counts.empty() ? face_counts_view.size() : face_counts.size();
  }

  const uint16_t *get_faces_data() const {
    return faces.empty() ? faces_view.data() : faces.data();
  }
  size_t get_faces_size() const {
    return faces.empty() ? faces_view.size() : faces.size();
  }

  const uint16_t *get_face_offsets_data() const {
    return face_offsets.empty() ? face_offsets_view.data()
                                : face_offsets.data();
  }
  size_t get_face_offsets_size() const {
    return face_offsets.empty() ? face_offsets_view.size()
                                : face_offsets.size();
  }

  // Helper for size accessors if needed, but direct vector access is preferred.
  size_t num_vertices() const { return vertices.size(); }
  size_t num_faces() const { return get_face_counts_size(); }
};

/**
 * @brief k-d Tree implementation for 3D points using static memory.
 * @details Stores points and allows nearest neighbor search.
 */
struct KDNode {
  // We store a COPY of the point to avoid lifetime issues,
  // or we could store an index if we had a reference to the point array.
  // Storing copy is safer and simpler for now.
  Vector point;
  int16_t originalIndex = -1; // Index in the source array
  int16_t axis = 0;           // 0=x, 1=y, 2=z
  int16_t left = -1;
  int16_t right = -1;
};

class KDTree {
public:
  static constexpr int MAX_K = 5;

  ArenaVector<KDNode> nodes;
  size_t nodeCount = 0;
  int rootIndex = -1;

  KDTree() = default;

  // Build from a Span of vectors, using arena for node storage and temporary
  // index sorting
  KDTree(Arena &arena, std::span<Vector> points) {
    clear();
    if (points.empty())
      return;

    size_t count = points.size();
    nodes.initialize(arena, count);

    // Temporary index allocation for building the tree
    int *indices = (int *)arena.allocate(count * sizeof(int), alignof(int));
    for (size_t i = 0; i < count; ++i)
      indices[i] = i;

    rootIndex = build(points, indices, count, 0);
  }

  void clear() {
    nodeCount = 0;
    rootIndex = -1;
  }

  // Nearest neighbor search
  // Returns up to K results (Nodes, so we get index + point)
  StaticCircularBuffer<KDNode, MAX_K> nearest(const Vector &target,
                                              size_t k = 1) const {
    StaticCircularBuffer<KDNode, MAX_K> result;
    if (rootIndex == -1 || k <= 0)
      return result;

    // Max-Heap behavior scratch buffer
    // We store pairs of <distanceSq, nodeIdx>
    struct Candidate {
      float dSq;
      int idx;
    };
    StaticCircularBuffer<Candidate, MAX_K> heap;

    auto push_heap = [&](float dSq, int idx) {
      if (heap.size() < static_cast<size_t>(k)) {
        heap.push_back({dSq, idx});
        // Maintain max-heap property: bubble up?
        // With small K (1-5), linear insert/sort is faster than heap complexity
        // Just unsorted push, then find max when needed.
      } else {
        // Replace worst if better
        float maxD = -1.0f;
        int maxI = -1;
        for (size_t i = 0; i < heap.size(); ++i) {
          if (heap[i].dSq > maxD) {
            maxD = heap[i].dSq;
            maxI = i;
          }
        }
        if (dSq < maxD) {
          heap[maxI] = {dSq, idx};
        }
      }
    };

    auto get_worst_dist = [&]() -> float {
      if (heap.size() < k)
        return FLT_MAX;
      float maxD = -1.0f;
      for (size_t i = 0; i < heap.size(); ++i) {
        if (heap[i].dSq > maxD)
          maxD = heap[i].dSq;
      }
      return maxD;
    };

    search_k(rootIndex, target, k, heap, push_heap, get_worst_dist);

    // Sort result by distance (closest first)
    std::sort(
        heap.begin(), heap.end(),
        [](const Candidate &a, const Candidate &b) { return a.dSq < b.dSq; });

    for (const auto &c : heap) {
      result.push_back(nodes[c.idx]);
    }
    return result;
  }

private:
  int build(std::span<Vector> points, int *indices, int count, int depth) {
    if (count <= 0)
      return -1;
    if (nodeCount >= nodes.capacity())
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

    int medianIdx = indices[mid];

    int newNodeIdx = nodeCount++;
    nodes[newNodeIdx].point = points[medianIdx];
    nodes[newNodeIdx].originalIndex = (int16_t)medianIdx; // Store index
    nodes[newNodeIdx].axis = (int16_t)axis;

    nodes[newNodeIdx].left = (int16_t)build(points, start, mid, depth + 1);
    nodes[newNodeIdx].right =
        (int16_t)build(points, start + mid + 1, count - mid - 1, depth + 1);

    return newNodeIdx;
  }

  // K-Nearest Search
  template <typename PushFn, typename MaxDistFn>
  void search_k(int nodeIdx, const Vector &target, int k, const auto &heap,
                PushFn &&push_heap, MaxDistFn &&get_worst_dist) const {
    if (nodeIdx == -1)
      return;

    const KDNode &node = nodes[nodeIdx];
    float dSq = distance_squared(node.point, target);
    float worstSq = get_worst_dist();

    if (dSq < worstSq) {
      push_heap(dSq, nodeIdx);
      worstSq = get_worst_dist(); // Update after push
    }

    float axisDist = (node.axis == 0)   ? (target.x - node.point.x)
                     : (node.axis == 1) ? (target.y - node.point.y)
                                        : (target.z - node.point.z);

    int near = axisDist < 0 ? node.left : node.right;
    int far = axisDist < 0 ? node.right : node.left;

    search_k(near, target, k, heap, push_heap, get_worst_dist);

    // Query worst again
    worstSq = get_worst_dist();
    if ((axisDist * axisDist) < worstSq) {
      search_k(far, target, k, heap, push_heap, get_worst_dist);
    }
  }

  void search(int nodeIdx, const Vector &target, float &bestDistSq,
              int &bestNodeIdx) const {
    if (nodeIdx == -1)
      return;

    const KDNode &node = nodes[nodeIdx];
    float dSq = distance_squared(node.point, target);

    if (dSq < bestDistSq) {
      bestDistSq = dSq;
      bestNodeIdx = nodeIdx;
    }

    float axisDist = (node.axis == 0)   ? (target.x - node.point.x)
                     : (node.axis == 1) ? (target.y - node.point.y)
                                        : (target.z - node.point.z);

    int near = axisDist < 0 ? node.left : node.right;
    int far = axisDist < 0 ? node.right : node.left;

    search(near, target, bestDistSq, bestNodeIdx);

    if (axisDist * axisDist < bestDistSq) {
      search(far, target, bestDistSq, bestNodeIdx);
    }
  }
};

/**
 * @brief Spatial Hashing for fast neighbor queries.
 */
class SpatialHash {
public:
  static constexpr size_t MAX_ENTRIES = 1024;
  static constexpr size_t TABLE_SIZE = 256;

  struct Entry {
    int id;
    int16_t next;
  };

  SpatialHash(float cellSize) : cellSize(cellSize) { clear(); }

  void clear() {
    std::fill(buckets.begin(), buckets.end(), -1);
    poolCount = 0;
  }

  void insert(const Vector &p, int id) {
    if (poolCount >= MAX_ENTRIES)
      return;

    long long h = hash(p);
    int16_t bucketIdx =
        (h % TABLE_SIZE + TABLE_SIZE) % TABLE_SIZE; // Ensure positive

    int idx = poolCount++;
    pool[idx].id = id;
    pool[idx].next = buckets[bucketIdx];
    buckets[bucketIdx] = idx;
  }

  // Returns fixed size buffer to avoid allocations
  StaticCircularBuffer<int, 64> query(const Vector &p) const {
    StaticCircularBuffer<int, 64> result;
    long long h = hash(p);
    int16_t bucketIdx = (h % TABLE_SIZE + TABLE_SIZE) % TABLE_SIZE;

    int16_t curr = buckets[bucketIdx];
    while (curr != -1) {
      result.push_back(pool[curr].id);
      curr = pool[curr].next;
      if (result.is_full())
        break;
    }
    return result;
  }

private:
  float cellSize;
  std::array<int16_t, TABLE_SIZE> buckets;
  std::array<Entry, MAX_ENTRIES> pool;
  size_t poolCount = 0;

  long long hash(const Vector &p) const {
    int x = static_cast<int>(floorf(p.x / cellSize));
    int y = static_cast<int>(floorf(p.y / cellSize));
    int z = static_cast<int>(floorf(p.z / cellSize));
    return ((long long)x * 73856093) ^ ((long long)y * 19349663) ^
           ((long long)z * 83492791);
  }
};

/**
 * @brief Projects a point onto the nearest surface of a mesh.
 * @param p The point to project.
 * @param mesh The mesh configuration.
 * @return The closest point on the mesh surface.
 */
inline Vector project_to_mesh(const Vector &p, const MeshState &mesh) {
  Vector origin(0, 0, 0);

  // Fallback: Closest Vertex
  if (mesh.vertices.empty())
    return p;
  Vector best = mesh.vertices[0];
  float minSq = distance_squared(p, best);
  for (size_t i = 1; i < mesh.vertices.size(); ++i) {
    float d = distance_squared(p, mesh.vertices[i]);
    if (d < minSq) {
      minSq = d;
      best = mesh.vertices[i];
    }
  }
  return best;
}
