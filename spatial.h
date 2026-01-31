/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include "geometry.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>
#include <unordered_map>



// KDTree commented out to avoid dynamic allocations (std::vector/unique_ptr)
/*
struct KDNode {
  // ...
};
class KDTree {
  // ...
};
*/


// --- SpatialHash (Moved from animation.h) ---
class SpatialHash {
public:
  SpatialHash(float cellSize) : cellSize(cellSize) {}
  void clear() { grid.clear(); }
  void insert(const Vector& p, int id) { grid[hash(p)].push_back(id); }
  std::vector<int> query(const Vector& p) {
    long long key = hash(p);
    return grid.count(key) ? grid.at(key) : std::vector<int>{};
  }
private:
  float cellSize;
  std::unordered_map<long long, std::vector<int>> grid;
  long long hash(const Vector& p) const {
    int x = static_cast<int>(floorf(p.i / cellSize));
    int y = static_cast<int>(floorf(p.j / cellSize));
    int z = static_cast<int>(floorf(p.k / cellSize));
    return ((long long)x * 73856093) ^ ((long long)y * 19349663) ^ ((long long)z * 83492791);
  }
};

// --- BVH Static Implementation ---

// --- BVH Static Implementation ---

namespace BVHImpl {
    // Recursive builder
    // range [start, end) in bvh.faceIndices
    inline int build_recursive(BVH& bvh, const MeshState& mesh, int start, int end, 
                        const int* faceOffsets, 
                        const Vector* centroids) {
        
        if (bvh.nodeCount >= BVH::MAX_NODES) return -1;
        int nodeIdx = bvh.nodeCount++;
        BVHNode& node = bvh.nodes[nodeIdx];
        
        // Compute AABB
        for(int i=start; i<end; ++i) {
            int faceIdx = bvh.faceIndices[i];
            int offset = faceOffsets[faceIdx];
            int count = mesh.face_counts[faceIdx];
            for(int k=0; k<count; ++k) {
                node.aabb.expand(mesh.vertices[mesh.faces[offset+k]]);
            }
        }
        
        int count = end - start;
        if (count <= 4) {
            node.firstFaceIndex = start;
            node.faceCount = count;
            return nodeIdx;
        }
        
        // Split
        Vector size = node.aabb.maxVal - node.aabb.minVal;
        int axis = 0;
        if (size.j > size.i && size.j > size.k) axis = 1;
        if (size.k > size.i && size.k > size.j) axis = 2;
        
        float mid = (axis == 0) ? (node.aabb.minVal.i + node.aabb.maxVal.i) :
                    (axis == 1) ? (node.aabb.minVal.j + node.aabb.maxVal.j) :
                                  (node.aabb.minVal.k + node.aabb.maxVal.k);
        mid *= 0.5f;
        
        // Partition
        // std::partition reorders elements in range
        // We use pointers to the fixed array in BVH
        int16_t* begin = &bvh.faceIndices[start];
        int16_t* endPtr = &bvh.faceIndices[end];
        
        auto it = std::partition(begin, endPtr,
            [&](int16_t faceIdx) {
                float val = (axis == 0) ? centroids[faceIdx].i :
                            (axis == 1) ? centroids[faceIdx].j : centroids[faceIdx].k;
                return val < mid;
            });
        
        int split = start + (int)std::distance(begin, it);
        
        if (split == start || split == end) {
             node.firstFaceIndex = start;
             node.faceCount = count;
             return nodeIdx;
        }
        
        node.left = build_recursive(bvh, mesh, start, split, faceOffsets, centroids);
        node.right = build_recursive(bvh, mesh, split, end, faceOffsets, centroids);
        
        return nodeIdx;
    }
}

inline void build_bvh(MeshState& mesh) {
    mesh.bvh.clear();
    if (mesh.num_faces == 0) return;
    
    // Stack buffers (ensure these fit in stack, 512 * 16 bytes ~ 8KB is fine)
    std::array<int, BVH::MAX_INDICES> faceOffsets;
    std::array<Vector, BVH::MAX_INDICES> centroids;
    
    int off = 0;
    for(size_t i=0; i<mesh.num_faces; ++i) {
        faceOffsets[i] = off;
        int count = mesh.face_counts[i];
        
        Vector c(0,0,0);
        for(int k=0; k<count; ++k) {
            c = c + mesh.vertices[mesh.faces[off+k]];
        }
        if (count > 0) c = c / (float)count;
        centroids[i] = c;
        
        // Init indices
        if (i < BVH::MAX_INDICES) mesh.bvh.faceIndices[i] = (int16_t)i;
        
        off += count;
    }
    
    int numFuncs = std::min((int)mesh.num_faces, BVH::MAX_INDICES);
    mesh.bvh.rootIndex = BVHImpl::build_recursive(mesh.bvh, mesh, 0, numFuncs, faceOffsets.data(), centroids.data());
}

struct HitResult {
    float dist = FLT_MAX;
    Vector point;
    bool hit = false;
};

inline HitResult bvh_intersect(const MeshState& mesh, const Vector& origin, const Vector& direction) {
    HitResult bestHit;
    if (mesh.bvh.rootIndex == -1 || mesh.bvh.nodeCount == 0) return bestHit;

    // Stack
    // Static stack? Depth likely small.
    std::array<int16_t, 64> stack;
    int stackPtr = 0;
    stack[stackPtr++] = mesh.bvh.rootIndex;

    // Precompute offsets if we can? No, we have to scan if we don't store them.
    // Optimization: We could store 'offset' in Cached BVH/Face list.
    // But MeshState doesn't have it.
    // Let's just scan linearly for offsets? Expensive.
    // Better: Build an offset map inside BVH?
    // BVH has indices array.
    // We can't easily change BVH structure now without breaking header, but we can do a quick scan table on stack if needed?
    // Actually, scanning faces array is O(F). 
    // MeshState::faces is a single flat array.
    // We MUST know the offset to access face K.
    // Solution: BVH stores 'firstIndex' into the face list? No, it stores 'faceIndex'.
    // CRITICAL: We need O(1) access to faces.
    // solids.h data doesn't provide offset table.
    // We can compute an offset table into `static` storage?
    // No, MeshState can point to different solids.
    // `build_bvh` computed offsets. We should ideally store them in BVH. 
    // BUT BVH struct is fixed.
    // Wait, `faceIndices` in BVH is int16_t.
    // Maybe we just store the *offset* into the faces array in `faceIndices` instead of the face number?
    // Then we need to know the count. `face_counts` is indexed by face number.
    // So we need both.
    // Alternatively, `BVH` struct in `geometry.h` has spare space?
    // `BVH` has `faceIndices`.
    // Let's assume we scan. It's slow but safe for now.
    // OR: Update geometry.h to include `std::array<int16_t, MAX_INDICES> faceOffsets` in BVH?
    // That requires another edit.
    // Let's stick to scanning for now, assuming N is small (SnubDodeca = 92 faces).
    // Scanning 92 integers is fast.
    
    while(stackPtr > 0) {
        int nodeIdx = stack[--stackPtr];
        const BVHNode& node = mesh.bvh.nodes[nodeIdx];
        
        if (!node.aabb.intersectRay(origin, direction)) continue;
        
        if (node.faceCount > 0) {
            // Leaf
            for(int i=0; i<node.faceCount; ++i) {
                int idxMap = node.firstFaceIndex + i; // Index into faceIndices
                int faceIdx = mesh.bvh.faceIndices[idxMap];
                
                // Find offset
                int offset = 0;
                for(int f=0; f<faceIdx; ++f) offset += mesh.face_counts[f];
                
                int count = mesh.face_counts[faceIdx];
                
                // Intersection
                Vector v0 = mesh.vertices[mesh.faces[offset]];
                for(int k=0; k<count-2; ++k) {
                    Vector v1 = mesh.vertices[mesh.faces[offset+k+1]];
                    Vector v2 = mesh.vertices[mesh.faces[offset+k+2]];
                    
                    // Moller-Trumbore
                    Vector h = cross(direction, v2 - v0); // edge2
                    float a = dot(v1 - v0, h); // edge1
                    
                    if(a > -1e-6f && a < 1e-6f) continue;
                    float f = 1.0f/a;
                    Vector s = origin - v0;
                    float u = f * dot(s, h);
                    if(u<0 || u>1) continue;
                    Vector q = cross(s, v1 - v0);
                    float v = f * dot(direction, q);
                    if(v<0 || u+v>1) continue;
                    float t = f * dot(v2 - v0, q);
                    
                    if(t > 1e-6f) {
                        if(!bestHit.hit || t < bestHit.dist) {
                            bestHit.hit = true;
                            bestHit.dist = t;
                            bestHit.point = origin + direction * t;
                        }
                    }
                }
            }
        } else {
            if (node.left != -1) stack[stackPtr++] = node.left;
            if (node.right != -1) stack[stackPtr++] = node.right;
        }
    }
    return bestHit;
}

inline Vector project_to_mesh(const Vector& p, const MeshState& mesh) {
    Vector origin(0,0,0);
    Vector dir = Vector(p).normalize();
    
    // Attempt BVH
    if (mesh.bvh.nodeCount > 0) {
        HitResult hit = bvh_intersect(mesh, origin, dir);
        if (hit.hit) return hit.point;
    }
    
    // Fallback: Closest Vertex
    if (mesh.num_vertices == 0) return p;
    Vector best = mesh.vertices[0];
    float minSq = distance_squared(p, best);
    for(size_t i=1; i<mesh.num_vertices; ++i) {
        float d = distance_squared(p, mesh.vertices[i]);
        if (d < minSq) {
            minSq = d;
            best = mesh.vertices[i];
        }
    }
    return best;
}

// Non-const overload (in case we want to support building?)
// We implemented build_bvh separately.
inline Vector project_to_mesh(const Vector& p, MeshState& mesh) {
    if (mesh.bvh.nodeCount == 0) {
        build_bvh(mesh); // Auto-build!
    }
    return project_to_mesh(p, static_cast<const MeshState&>(mesh));
}
