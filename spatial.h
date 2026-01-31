/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <map>
#include <unordered_map>
#include <string>

/**
 * @brief k-d Tree implementation for 3D points.
 * @details Stores points and allows nearest neighbor search.
 */
struct KDNode {
  Vector point;
  int axis; // 0=x, 1=y, 2=z
  std::unique_ptr<KDNode> left;
  std::unique_ptr<KDNode> right;
  
  KDNode(const Vector& p, int ax) : point(p), axis(ax) {}
};

class KDTree {
public:
  KDTree() = default;
  
  // Build from a list of vectors
  KDTree(const std::vector<Vector>& points) {
      std::vector<Vector> pts = points; // Copy to sort
      root = build(pts, 0);
  }

  // Nearest neighbor search
  std::vector<Vector> nearest(const Vector& target, int k = 1) const {
      std::vector<NodeDist> best;
      search(root.get(), target, k, best);
      
      std::vector<Vector> result;
      for(const auto& b : best) {
          result.push_back(b.node->point);
      }
      return result;
  }

private:
  struct NodeDist {
      float distSq;
      KDNode* node;
      
      bool operator<(const NodeDist& other) const {
          return distSq < other.distSq;
      }
  };

  std::unique_ptr<KDNode> root;

  std::unique_ptr<KDNode> build(std::vector<Vector>& points, int depth) {
      if (points.empty()) return nullptr;
      
      int axis = depth % 3;
      std::sort(points.begin(), points.end(), [axis](const Vector& a, const Vector& b) {
          if (axis == 0) return a.i < b.i;
          if (axis == 1) return a.j < b.j;
          return a.k < b.k;
      });
      
      size_t mid = points.size() / 2;
      auto node = std::make_unique<KDNode>(points[mid], axis);
      
      std::vector<Vector> leftPts(points.begin(), points.begin() + mid);
      std::vector<Vector> rightPts(points.begin() + mid + 1, points.end());
      
      node->left = build(leftPts, depth + 1);
      node->right = build(rightPts, depth + 1);
      
      return node;
  }
  
  void search(KDNode* node, const Vector& target, int k, std::vector<NodeDist>& best) const {
      if (!node) return;
      
      float dSq = distance_squared(node->point, target);
      float axisDist = 0;
      if (node->axis == 0) axisDist = target.i - node->point.i;
      else if (node->axis == 1) axisDist = target.j - node->point.j;
      else axisDist = target.k - node->point.k;
      
      // Update best
      if (best.size() < k || dSq < best.back().distSq) {
          best.push_back({ dSq, node });
          std::sort(best.begin(), best.end());
          if (best.size() > k) best.pop_back();
      }
      
      KDNode* near = axisDist < 0 ? node->left.get() : node->right.get();
      KDNode* far = axisDist < 0 ? node->right.get() : node->left.get();
      
      search(near, target, k, best);
      
      if (best.size() < k || (axisDist * axisDist) < best.back().distSq) {
          search(far, target, k, best);
      }
  }
};

/**
 * @brief Axis-Aligned Bounding Box.
 */
struct AABB {
  Vector minVal;
  Vector maxVal;
  
  AABB() : minVal(INFINITY, INFINITY, INFINITY), maxVal(-INFINITY, -INFINITY, -INFINITY) {}
  
  void expand(const Vector& p) {
      minVal.i = std::min(minVal.i, p.i);
      minVal.j = std::min(minVal.j, p.j);
      minVal.k = std::min(minVal.k, p.k);
      
      maxVal.i = std::max(maxVal.i, p.i);
      maxVal.j = std::max(maxVal.j, p.j);
      maxVal.k = std::max(maxVal.k, p.k);
  }
  
  void unionWith(const AABB& box) {
      minVal.i = std::min(minVal.i, box.minVal.i);
      minVal.j = std::min(minVal.j, box.minVal.j);
      minVal.k = std::min(minVal.k, box.minVal.k);
      
      maxVal.i = std::max(maxVal.i, box.maxVal.i);
      maxVal.j = std::max(maxVal.j, box.maxVal.j);
      maxVal.k = std::max(maxVal.k, box.maxVal.k);
  }
  
  // Ray intersection (Slab method)
  bool intersectRay(const Vector& origin, const Vector& direction) const {
      float tmin = (minVal.i - origin.i) / direction.i;
      float tmax = (maxVal.i - origin.i) / direction.i;
      
      if (tmin > tmax) std::swap(tmin, tmax);
      
      float tymin = (minVal.j - origin.j) / direction.j;
      float tymax = (maxVal.j - origin.j) / direction.j;
      
      if (tymin > tymax) std::swap(tymin, tymax);
      
      if ((tmin > tymax) || (tymin > tmax)) return false;
      
      if (tymin > tmin) tmin = tymin;
      if (tymax < tmax) tmax = tymax;
      
      float tzmin = (minVal.k - origin.k) / direction.k;
      float tzmax = (maxVal.k - origin.k) / direction.k;
      
      if (tzmin > tzmax) std::swap(tzmin, tzmax);
      
      if ((tmin > tzmax) || (tzmin > tmax)) return false;
      
      return true;
  }
};

/**
 * @brief Node in the BVH.
 */
struct BVHNode {
  AABB aabb;
  std::unique_ptr<BVHNode> left;
  std::unique_ptr<BVHNode> right;
  std::vector<int> indices; // Face indices, leaf only
};

struct HitResult {
    float dist = FLT_MAX;
    Vector point;
    int faceIndex = -1;
    bool hit = false;
};

/**
 * @brief Bounding Volume Hierarchy for optimized ray casting.
 */
class BVH {
public:
  BVH(const MeshState& mesh) : mesh(mesh) {}
  
  void build() {
      std::vector<int> indices;
      std::vector<Vector> centroids;
      
      const int* face_ptr = mesh.faces;
      for(size_t i=0; i<mesh.num_faces; ++i) {
          indices.push_back((int)i);
          
          size_t count = mesh.face_counts[i];
          Vector c(0,0,0);
          for(size_t k=0; k<count; ++k) {
              c = c + mesh.vertices[face_ptr[k]];
          }
          c = c / static_cast<float>(count);
          centroids.push_back(c);
          
          face_ptr += count;
      }
      
      root = buildRecursive(indices, centroids);
  }
  
  HitResult intersectRay(const Vector& origin, const Vector& direction) const {
      HitResult bestHit;
      if (!root) return bestHit;
      
      std::vector<BVHNode*> stack;
      stack.push_back(root.get());
      
      while(!stack.empty()) {
          BVHNode* node = stack.back();
          stack.pop_back();
          
          if (!node->aabb.intersectRay(origin, direction)) continue;
          
          if (!node->indices.empty()) {
              // Leaf: Check faces
              for(int idx : node->indices) {
                  // Reconstruct face access from flat arrays
                  // Optimization: Could cache face offsets in build?
                  // For now, scan to find face (O(N) unfortunately in bare MeshState if redundant, likely fast enough for leaf)
                  // Wait, MeshState is array based. We need the specific face pointer.
                  // Since 'idx' is the K-th face, we access via face_counts scanning?
                  // MeshState doesn't store offsets.
                  // BETTER: BVH build should store offsets or direct pointers?
                  // Pointers to int* in MeshState faces array? 
                  // But MeshState faces might move if copied? No, they point to static data usually.
                  // Let's rely on finding it. Or store precomputed offsets in BVH?
                  
                  // For correctness and speed, let's scan.
                  // Since this scanning is O(N) per leaf hit, it defeats BVH purpose if N is large.
                  // MeshState 'faces' pointer typically points to static data in solids.h.
                  // We can pre-calculate offsets in BVH.
                  
                  int offset = faceOffsets[idx];
                  int count = mesh.face_counts[idx];
                  const int* face = &mesh.faces[offset];
                  
                  // Triangulate fan logic
                  Vector v0 = mesh.vertices[face[0]];
                  for(int i=0; i<count-2; ++i) {
                      Vector v1 = mesh.vertices[face[i+1]];
                      Vector v2 = mesh.vertices[face[i+2]];
                      
                      Vector edge1 = v1 - v0;
                      Vector edge2 = v2 - v0;
                      
                      Vector h = cross(direction, edge2);
                      float a = dot(edge1, h);
                      
                      if (a > -1e-6f && a < 1e-6f) continue;
                      
                      float f = 1.0f / a;
                      Vector s = origin - v0;
                      float u = f * dot(s, h);
                      
                      if (u < 0.0f || u > 1.0f) continue;
                      
                      Vector q = cross(s, edge1);
                      float v = f * dot(direction, q);
                      
                      if (v < 0.0f || u + v > 1.0f) continue;
                      
                      float t = f * dot(edge2, q);
                      
                      if (t > 1e-6f) {
                          if (!bestHit.hit || t < bestHit.dist) {
                              bestHit.hit = true;
                              bestHit.dist = t;
                              bestHit.faceIndex = idx;
                              bestHit.point = origin + direction * t;
                          }
                      }
                  }
              }
          } else {
              if (node->left) stack.push_back(node->left.get());
              if (node->right) stack.push_back(node->right.get());
          }
      }
      return bestHit;
  }

private:
  const MeshState& mesh;
  std::unique_ptr<BVHNode> root;
  std::vector<int> faceOffsets; // Cache offsets for faces
  
  std::unique_ptr<BVHNode> buildRecursive(const std::vector<int>& indices, const std::vector<Vector>& centroids) {
      auto node = std::make_unique<BVHNode>();
      
      // Compute cache if root
      if (faceOffsets.empty()) {
          faceOffsets.resize(mesh.num_faces);
          int off = 0;
          for(size_t i=0; i<mesh.num_faces; ++i) {
              faceOffsets[i] = off;
              off += mesh.face_counts[i];
          }
      }
      
      // 1. Compute AABB
      for(int idx : indices) {
          int offset = faceOffsets[idx];
          int count = mesh.face_counts[idx];
          const int* face = &mesh.faces[offset];
          for(int k=0; k<count; ++k) {
              node->aabb.expand(mesh.vertices[face[k]]);
          }
      }
      
      // 2. Leaf
      if (indices.size() <= 4) {
          node->indices = indices;
          return node;
      }
      
      // 3. Split
      Vector size = node->aabb.maxVal - node->aabb.minVal;
      int axis = 0; // x
      if (size.j > size.i && size.j > size.k) axis = 1; // y
      if (size.k > size.i && size.k > size.j) axis = 2; // z
      
      float mid = 0;
      if (axis == 0) mid = (node->aabb.minVal.i + node->aabb.maxVal.i) * 0.5f;
      else if (axis == 1) mid = (node->aabb.minVal.j + node->aabb.maxVal.j) * 0.5f;
      else mid = (node->aabb.minVal.k + node->aabb.maxVal.k) * 0.5f;
      
      std::vector<int> left, right;
      for(int idx : indices) {
          float val = 0;
          if (axis == 0) val = centroids[idx].i;
          else if (axis == 1) val = centroids[idx].j;
          else val = centroids[idx].k;
          
          if (val < mid) left.push_back(idx);
          else right.push_back(idx);
      }
      
      if (left.empty() || right.empty()) {
          node->indices = indices;
          return node;
      }
      
      node->left = buildRecursive(left, centroids);
      node->right = buildRecursive(right, centroids);
      
      return node;
  }
};

/**
 * @brief Spatial Hashing for particles.
 */
class SpatialHash {
public:
  SpatialHash(float cellSize) : cellSize(cellSize) {}

  void clear() {
    grid.clear();
  }

  void insert(const Vector& p, int id) {
    long long key = hash(p);
    grid[key].push_back(id);
  }

  std::vector<int> query(const Vector& p, float radius) {
    std::vector<int> result;
    // Basic implementation: check neighbor cells
    int range = static_cast<int>(ceil(radius / cellSize));
    int cx = static_cast<int>(floor(p.i / cellSize));
    int cy = static_cast<int>(floor(p.j / cellSize));
    int cz = static_cast<int>(floor(p.k / cellSize));
    
    for (int x = cx - range; x <= cx + range; ++x) {
        for (int y = cy - range; y <= cy + range; ++y) {
            for (int z = cz - range; z <= cz + range; ++z) {
                long long k = encode(x,y,z);
                if (grid.count(k)) {
                    const auto& bucket = grid.at(k);
                    result.insert(result.end(), bucket.begin(), bucket.end());
                }
            }
        }
    }
    return result;
  }
  
  // Single bucket query
  std::vector<int> query(const Vector& p) {
    long long key = hash(p);
    if (grid.count(key)) return grid.at(key);
    return {};
  }

private:
  float cellSize;
  std::unordered_map<long long, std::vector<int>> grid;

  long long encode(int x, int y, int z) const {
       // Simple packing for hash
      return ((long long)x * 73856093) ^ ((long long)y * 19349663) ^ ((long long)z * 83492791);
  }
  
  long long hash(const Vector& p) const {
    int x = static_cast<int>(floorf(p.i / cellSize));
    int y = static_cast<int>(floorf(p.j / cellSize));
    int z = static_cast<int>(floorf(p.k / cellSize));
    return encode(x,y,z);
  }
};

/**
 * @brief Raycasts a point onto the mesh surface using BVH if available.
 * @param p The point to project.
 * @param mesh The target MeshState.
 * @return The projected point.
 */
inline Vector project_to_mesh(const Vector& p, MeshState& mesh) {
    if (!mesh.bvh) {
        // Construct BVH if missing (and logic allows modification of MeshState)
        // Since MeshState struct is light and we have shared_ptr, this is safe.
        // We need to modify the mesh to store the BVH cache.
        // Cast away constness or require non-const ref? Reference is non-const.
        mesh.bvh = std::make_shared<BVH>(mesh);
        mesh.bvh->build();
    }
    
    Vector dir = Vector(p).normalize();
    Vector origin(0,0,0);
    
    HitResult hit = mesh.bvh->intersectRay(origin, dir);
    
    if (hit.hit) {
        return hit.point;
    }
    
    // Fallback logic
     if (mesh.num_vertices == 0) return p;
     
     Vector best = mesh.vertices[0];
     float minSq = distance_squared(p, best);
     for(size_t i=1; i<mesh.num_vertices; ++i) {
         float sq = distance_squared(p, mesh.vertices[i]);
         if (sq < minSq) {
             minSq = sq;
             best = mesh.vertices[i];
         }
     }
     return best;
}

// Overload for const MeshState checks BVH existence but cannot build it
inline Vector project_to_mesh(const Vector& p, const MeshState& mesh) {
    if (mesh.bvh) {
        Vector dir = Vector(p).normalize();
        Vector origin(0,0,0);
        HitResult hit = mesh.bvh->intersectRay(origin, dir);
        if (hit.hit) return hit.point;
    }
    // Fallback manual scan
     // ... (Duplicate fallback logic or cast constness if we really want to cache?)
     // For now, simple fallback.
     if (mesh.num_vertices == 0) return p;
     Vector best = mesh.vertices[0];
     float minSq = distance_squared(p, best);
     for(size_t i=1; i<mesh.num_vertices; ++i) {
         float sq = distance_squared(p, mesh.vertices[i]);
         if (sq < minSq) {
             minSq = sq;
             best = mesh.vertices[i];
         }
     }
     return best;
}
