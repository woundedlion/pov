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
#include <vector>
#include <unordered_map>
#include <memory>
#include <string>
#include <cfloat>
#include "static_circular_buffer.h"

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
    if (p.i < minVal.i)
      minVal.i = p.i;
    if (p.j < minVal.j)
      minVal.j = p.j;
    if (p.k < minVal.k)
      minVal.k = p.k;

    if (p.i > maxVal.i)
      maxVal.i = p.i;
    if (p.j > maxVal.j)
      maxVal.j = p.j;
    if (p.k > maxVal.k)
      maxVal.k = p.k;
  }

  // Optimize union
  void unionWith(const AABB &box) {
    if (box.minVal.i < minVal.i)
      minVal.i = box.minVal.i;
    if (box.minVal.j < minVal.j)
      minVal.j = box.minVal.j;
    if (box.minVal.k < minVal.k)
      minVal.k = box.minVal.k;

    if (box.maxVal.i > maxVal.i)
      maxVal.i = box.maxVal.i;
    if (box.maxVal.j > maxVal.j)
      maxVal.j = box.maxVal.j;
    if (box.maxVal.k > maxVal.k)
      maxVal.k = box.maxVal.k;
  }

  bool intersectRay(const Vector &origin, const Vector &direction) const {
    float tmin = (minVal.i - origin.i) / direction.i;
    float tmax = (maxVal.i - origin.i) / direction.i;
    if (tmin > tmax)
      std::swap(tmin, tmax);

    float tymin = (minVal.j - origin.j) / direction.j;
    float tymax = (maxVal.j - origin.j) / direction.j;
    if (tymin > tymax)
      std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
      return false;

    if (tymin > tmin)
      tmin = tymin;
    if (tymax < tmax)
      tmax = tymax;

    float tzmin = (minVal.k - origin.k) / direction.k;
    float tzmax = (maxVal.k - origin.k) / direction.k;
    if (tzmin > tzmax)
      std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
      return false;

    return true;
  }
};

/**
 * @brief Represents the state of a mesh using static storage to avoid heap
 * allocations.
 */
struct MeshState {
  // Config constants (kept for compatibility)
  static constexpr size_t MAX_VERTS = 4096;
  static constexpr size_t MAX_FACES = 2048;

  std::vector<Vector> vertices;
  std::vector<uint8_t> face_counts;
  std::vector<int> faces;
  std::vector<uint16_t> face_offsets;

  void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    face_offsets.clear();
  }

  // Helper for size accessors if needed, but direct vector access is preferred.
  size_t num_vertices() const { return vertices.size(); }
  size_t num_faces() const { return face_counts.size(); }
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
  static constexpr int MAX_NODES = 2048;
  static constexpr int MAX_K = 5;

  std::array<KDNode, MAX_NODES> nodes;
  int nodeCount = 0;
  int rootIndex = -1;

  KDTree() = default;

  // Build from a Span of vectors
  KDTree(std::span<Vector> points) {
    clear();
    if (points.empty())
      return;

    // We need a mutable index array to partition
    std::array<int, MAX_NODES> indices;
    size_t count = std::min((size_t)MAX_NODES, points.size());
    for (size_t i = 0; i < count; ++i)
      indices[i] = i;

    rootIndex = build(points, indices.data(), count, 0);
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
    if (nodeCount >= MAX_NODES)
      return -1;

    int axis = depth % 3;
    int mid = count / 2;

    // nth_element / partition
    // We need to partition the indices based on the points they point to
    auto *start = indices;
    auto *end = indices + count;

    std::nth_element(start, start + mid, end, [&](int a, int b) {
      float va = (axis == 0)   ? points[a].i
                 : (axis == 1) ? points[a].j
                               : points[a].k;
      float vb = (axis == 0)   ? points[b].i
                 : (axis == 1) ? points[b].j
                               : points[b].k;
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

    float axisDist = (node.axis == 0)   ? (target.i - node.point.i)
                     : (node.axis == 1) ? (target.j - node.point.j)
                                        : (target.k - node.point.k);

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

    float axisDist = (node.axis == 0)   ? (target.i - node.point.i)
                     : (node.axis == 1) ? (target.j - node.point.j)
                                        : (target.k - node.point.k);

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
    int x = static_cast<int>(floorf(p.i / cellSize));
    int y = static_cast<int>(floorf(p.j / cellSize));
    int z = static_cast<int>(floorf(p.k / cellSize));
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
  Vector dir = Vector(p).normalize();

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
