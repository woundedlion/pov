/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"
#include "spatial.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <span>
#include <cmath>
#include <iostream>

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
  std::vector<Vector> vertices;
  std::vector<std::vector<int>> faces;

  // Cache
  mutable bool cache_valid = false;
  mutable KDTree kdTree;
  mutable std::vector<std::vector<int>> adjacency;

  void clear_cache() const {
    cache_valid = false;
    kdTree.clear();
    adjacency.clear();
  }
};

struct HalfEdge {
  HEVertex *vertex = nullptr; /**< Vertex at the end of this half-edge. */
  HEFace *face = nullptr;     /**< Face this half-edge belongs to. */
  HalfEdge *next = nullptr;   /**< Next half-edge in the face loop. */
  HalfEdge *prev = nullptr;   /**< Previous half-edge in the face loop. */
  HalfEdge *pair = nullptr;   /**< Opposite half-edge. */
};

struct HEVertex {
  Vector position;
  HalfEdge *halfEdge =
      nullptr; /**< One of the half-edges pointing to this vertex. */
};

struct HEFace {
  HalfEdge *halfEdge =
      nullptr; /**< One of the half-edges bordering this face. */
  int vertexCount = 0;
  uint32_t intrinsicHash = 0;

  void compute_properties() {
    if (!halfEdge)
      return;

    // 1. Collect vertices
    std::vector<Vector> verts;
    HalfEdge *he = halfEdge;
    HalfEdge *start = he;
    int safety = 0;
    do {
      verts.push_back(he->vertex->position);
      he = he->next;
      safety++;
    } while (he && he != start && safety < 100);

    vertexCount = (int)verts.size();

    // 2. Angles
    std::vector<int> angles;
    if (vertexCount >= 3) {
      for (int i = 0; i < vertexCount; ++i) {
        const Vector &prev = verts[(i - 1 + vertexCount) % vertexCount];
        const Vector &curr = verts[i];
        const Vector &next = verts[(i + 1) % vertexCount];

        Vector v1 = (prev - curr).normalize();
        Vector v2 = (next - curr).normalize();

        // angle_between returns radians.
        float ang = angle_between(v1, v2);
        angles.push_back((int)std::round(ang * 180.0f / PI_F));
      }
      std::sort(angles.begin(), angles.end());
    }

    // 3. Hash (JS style hash32)
    auto hash32 = [](uint32_t n, uint32_t seed) -> uint32_t {
      n = (n ^ seed) * 0x5bd1e995;
      n ^= n >> 15;
      return n * 0x97455bcd;
    };

    uint32_t h = hash32(static_cast<uint32_t>(vertexCount), 0x12345678);
    for (int a : angles) {
      h = hash32(static_cast<uint32_t>(a), h);
    }
    intrinsicHash = h;
  }
};

class HalfEdgeMesh {
public:
  std::vector<HEVertex> vertices;
  std::vector<HEFace> faces;
  std::vector<HalfEdge> halfEdges;

  explicit HalfEdgeMesh(const PolyMesh &mesh) {
    vertices.reserve(mesh.vertices.size());
    for (const auto &v : mesh.vertices) {
      vertices.push_back({v, nullptr});
    }

    size_t totalHE = 0;
    for (const auto &f : mesh.faces)
      totalHE += f.size();
    halfEdges.reserve(totalHE);
    faces.reserve(mesh.faces.size());

    std::map<std::pair<int, int>, HalfEdge *> edgeMap;

    for (const auto &f : mesh.faces) {
      faces.emplace_back();
      HEFace *currentFace = &faces.back();
      size_t count = f.size();
      size_t faceStartHeIdx = halfEdges.size();

      // Allocate edges for this face
      for (size_t i = 0; i < count; ++i) {
        halfEdges.emplace_back();
      }

      for (size_t i = 0; i < count; ++i) {
        int u = f[i];
        int v = f[(i + 1) % count];

        HalfEdge *he = &halfEdges[faceStartHeIdx + i];

        // Link basic geometry
        he->vertex = &vertices[v]; // Points TO v
        he->face = currentFace;

        // Circular links
        he->next = &halfEdges[faceStartHeIdx + (i + 1) % count];
        he->prev = &halfEdges[faceStartHeIdx + (i - 1 + count) % count];

        // Vertex ref (just needs one incoming edge)
        vertices[v].halfEdge = he;

        // Pair lookup
        if (edgeMap.count({v, u})) {
          HalfEdge *neighbor = edgeMap[{v, u}];
          he->pair = neighbor;
          neighbor->pair = he;
          edgeMap.erase({v, u});
        } else {
          edgeMap[{u, v}] = he;
        }
      }
      currentFace->halfEdge = &halfEdges[faceStartHeIdx];
    }

    // Compute Properties
    for (auto &f : faces)
      f.compute_properties();
  }

  HalfEdgeMesh(const MeshState &mesh) {
    // Reconstitute from flat arrays if needed
    PolyMesh p;
    p.vertices = mesh.vertices;
    int off = 0;
    for (uint8_t c : mesh.face_counts) {
      std::vector<int> f;
      for (int k = 0; k < c; ++k)
        f.push_back(mesh.faces[off + k]);
      p.faces.push_back(f);
      off += c;
    }
    *this = HalfEdgeMesh(p);
  }
};

/**
 * @brief Structure returned by compile_hankin.
 */
struct HankinInstruction {
  Vector pCorner; /**< Corner vertex position. */
  Vector pPrev;   /**< Previous vertex position. */
  Vector pNext;   /**< Next vertex position. */
  int idxM1;      /**< Index of first midpoint (static vertex). */
  int idxM2;      /**< Index of second midpoint (static vertex). */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 */
struct CompiledHankin {
  std::vector<Vector> staticVertices;  /**< Midpoints that don't move. */
  std::vector<Vector> dynamicVertices; /**< Intersection points that move. */
  std::vector<HankinInstruction>
      dynamicInstructions; /**< Instructions to update dynamic vertices. */
  std::vector<std::vector<int>> faces; /**< Resulting face topology. */
  int staticOffset; /**< Offset where dynamic vertices start. */
};

/**
 * @brief Operations on meshes (Dual, Hankin, etc.).
 */
namespace MeshOps {

/**
 * @brief Computes the dual of a mesh.
 */
template <typename MeshT> static MeshT dual(const MeshT &mesh) {
  HalfEdgeMesh heMesh(mesh);
  MeshT dualMesh;
  std::unordered_map<HEFace *, int> faceToVertIdx;

  // New Vertices (Centroids of original faces)
  for (size_t i = 0; i < heMesh.faces.size(); ++i) {
    HEFace *face = &heMesh.faces[i];
    Vector c(0, 0, 0);
    int count = 0;
    HalfEdge *he = face->halfEdge;
    HalfEdge *start = he;
    do {
      c = c + he->vertex->position;
      count++;
      he = he->next;
    } while (he != start);

    c = c / static_cast<float>(count);
    dualMesh.vertices.push_back(c.normalize());
    faceToVertIdx[face] = static_cast<int>(i);
  }

  // New Faces (Cycles around original vertices)
  std::vector<HEVertex *> visitedVerts; // Naive set replacement

  // Helper to check if visited
  auto isVisited = [&](HEVertex *v) {
    return std::find(visitedVerts.begin(), visitedVerts.end(), v) !=
           visitedVerts.end();
  };

  for (auto &heStart : heMesh.halfEdges) {
    if (!heStart.prev)
      continue; // Safety

    HEVertex *origin = heStart.prev->vertex; // The vertex 'start' comes FROM
    if (isVisited(origin))
      continue;
    visitedVerts.push_back(origin);

    std::vector<int> faceIndices;
    HalfEdge *curr = &heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;

    do {
      if (!curr->face)
        break;
      faceIndices.push_back(faceToVertIdx[curr->face]);

      if (!curr->pair)
        break;
      curr = curr->pair->next;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (faceIndices.size() > 2) {
      std::reverse(faceIndices.begin(), faceIndices.end()); // Maintain CCW
      dualMesh.faces.push_back(faceIndices);
    }
  }

  return dualMesh;
}

/**
 * @brief Compiles the topology for a Hankin pattern.
 */
template <typename MeshT>
static CompiledHankin compile_hankin(const MeshT &mesh) {
  HalfEdgeMesh heMesh(mesh);
  CompiledHankin compiled;

  std::map<HalfEdge *, int> heToMidpointIdx;
  std::map<HalfEdge *, int> heToDynamicIdx;

  // Helper to get/create midpoint index
  auto getMidpointIdx = [&](HalfEdge *he) {
    if (heToMidpointIdx.count(he))
      return heToMidpointIdx[he];
    if (he->pair && heToMidpointIdx.count(he->pair))
      return heToMidpointIdx[he->pair];

    Vector pA =
        he->prev ? he->prev->vertex->position : he->pair->vertex->position;
    Vector pB = he->vertex->position;
    Vector mid = (pA + pB) * 0.5f;
    mid.normalize();

    compiled.staticVertices.push_back(mid);
    int idx = static_cast<int>(compiled.staticVertices.size()) - 1;
    heToMidpointIdx[he] = idx;
    if (he->pair)
      heToMidpointIdx[he->pair] = idx;
    return idx;
  };

  // Ensure all midpoints
  for (auto &he : heMesh.halfEdges) {
    getMidpointIdx(&he);
  }

  compiled.staticOffset = static_cast<int>(compiled.staticVertices.size());

  // Star faces
  for (auto &face : heMesh.faces) {
    std::vector<int> starFaceIndices;
    HalfEdge *he = face.halfEdge;
    HalfEdge *startHe = he;

    do {
      HalfEdge *prev = he->prev;
      HalfEdge *curr = he;

      int idxM1 = getMidpointIdx(prev);
      int idxM2 = getMidpointIdx(curr);

      Vector pCorner = prev->vertex->position;
      Vector pPrev = (prev->prev ? prev->prev->vertex->position
                                 : prev->pair->vertex->position);
      Vector pNext = curr->vertex->position;

      compiled.dynamicInstructions.push_back(
          {pCorner, pPrev, pNext, idxM1, idxM2});

      int dynIdx = static_cast<int>(compiled.dynamicVertices.size());
      heToDynamicIdx[curr] = dynIdx;
      compiled.dynamicVertices.emplace_back(); // Placeholder

      starFaceIndices.push_back(idxM1);
      starFaceIndices.push_back(compiled.staticOffset + dynIdx);

      he = he->next;
    } while (he != startHe);

    compiled.faces.push_back(starFaceIndices);
  }

  // Rosette faces
  std::vector<HEVertex *> visitedVerts;
  auto isVisited = [&](HEVertex *v) {
    return std::find(visitedVerts.begin(), visitedVerts.end(), v) !=
           visitedVerts.end();
  };

  for (auto &heStart : heMesh.halfEdges) {
    if (!heStart.prev)
      continue;
    HEVertex *origin = heStart.prev->vertex;
    if (isVisited(origin))
      continue;
    visitedVerts.push_back(origin);

    std::vector<int> rosetteIndices;
    HalfEdge *curr = &heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;

    do {
      rosetteIndices.push_back(heToMidpointIdx[curr]); // Static
      HalfEdge *nextEdge = curr->pair ? curr->pair->next : nullptr;
      if (!nextEdge)
        break;
      rosetteIndices.push_back(compiled.staticOffset +
                               heToDynamicIdx[nextEdge]); // Dynamic
      curr = nextEdge;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (rosetteIndices.size() > 2) {
      std::reverse(rosetteIndices.begin(),
                   rosetteIndices.end()); // Fix inward normals
      compiled.faces.push_back(rosetteIndices);
    }
  }

  return compiled;
}

/**
 * @brief Updates a compiled Hankin pattern.
 */
template <typename MeshT> // Templated just to match style, though MeshT output
                          // is expected
static MeshT update_hankin(CompiledHankin &compiled, float angle) {
  for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
    const auto &instr = compiled.dynamicInstructions[i];
    Vector m1 = compiled.staticVertices[instr.idxM1];
    Vector m2 = compiled.staticVertices[instr.idxM2];

    Vector nEdge1 = cross(instr.pPrev, instr.pCorner).normalize();
    Quaternion q1 = make_rotation(m1, angle);
    Vector nHankin1 = rotate(nEdge1, q1);

    Vector nEdge2 = cross(instr.pCorner, instr.pNext).normalize();
    Quaternion q2 = make_rotation(m2, -angle);
    Vector nHankin2 = rotate(nEdge2, q2);

    Vector intersect = cross(nHankin1, nHankin2);

    float lenSq = dot(intersect, intersect);
    if (lenSq < 1e-6f) { // Increased epsilon for float precision
      // Degenerate/Parallel.
      // When lines are parallel (e.g. 45 deg on tetrahedron), they often form a
      // straight line passing through the face center. The correct "vertex" for
      // the pattern is the midpoint of the edge connectors m1 and m2.
      intersect = (m1 + m2).normalize();
    }

    // Chirality
    if (dot(intersect, instr.pCorner) < 0)
      intersect = -intersect;

    compiled.dynamicVertices[i] = intersect.normalize();
  }

  MeshT result;
  result.vertices = compiled.staticVertices;
  result.vertices.insert(result.vertices.end(),
                         compiled.dynamicVertices.begin(),
                         compiled.dynamicVertices.end());
  result.faces = compiled.faces;
  return result;
}

/**
 * @brief Helper to do full hankin generation.
 */
template <typename MeshT> static MeshT hankin(const MeshT &mesh, float angle) {
  auto compiled = compile_hankin(mesh);
  return update_hankin<MeshT>(compiled, angle);
}

// --- CONWAY OPERATORS ---

/**
 * @brief Normalizes all vertices in the mesh to the unit sphere.
 */
template <typename MeshT> static void normalize(MeshT &mesh) {
  for (auto &v : mesh.vertices) {
    v = v.normalize();
  }
}

/**
 * @brief Kis operator: Raises a pyramid on each face.
 */
template <typename MeshT> static MeshT kis(const MeshT &mesh) {
  MeshT result;
  result.vertices = mesh.vertices; // Copy existing

  for (const auto &f : mesh.faces) {
    // Add centroid
    Vector centroid(0, 0, 0);
    for (int vi : f) {
      centroid = centroid + mesh.vertices[vi];
    }
    if (!f.empty())
      centroid = centroid / static_cast<float>(f.size());

    result.vertices.push_back(centroid);
    int centerIdx = static_cast<int>(result.vertices.size()) - 1;

    // Create triangles
    for (size_t i = 0; i < f.size(); ++i) {
      int vi = f[i];
      int vj = f[(i + 1) % f.size()];
      result.faces.push_back({vi, vj, centerIdx});
    }
  }

  normalize(result);
  return result;
}

/**
 * @brief Ambo operator: Truncates vertices to edge midpoints.
 */
template <typename MeshT> static MeshT ambo(const MeshT &mesh) {
  MeshT result;
  std::map<std::pair<int, int>, int> edgeMap;

  // 1. Create vertices at edge midpoints
  for (const auto &f : mesh.faces) {
    for (size_t i = 0; i < f.size(); ++i) {
      int vi = f[i];
      int vj = f[(i + 1) % f.size()];
      int u = std::min(vi, vj);
      int v = std::max(vi, vj);

      if (edgeMap.find({u, v}) == edgeMap.end()) {
        Vector mid = (mesh.vertices[vi] + mesh.vertices[vj]) * 0.5f;
        result.vertices.push_back(mid);
        edgeMap[{u, v}] = static_cast<int>(result.vertices.size()) - 1;
      }
    }
  }

  // 2. Create faces
  // A. Shrink old faces
  for (const auto &f : mesh.faces) {
    std::vector<int> faceVerts;
    for (size_t i = 0; i < f.size(); ++i) {
      int vi = f[i];
      int vj = f[(i + 1) % f.size()];
      int u = std::min(vi, vj);
      int v = std::max(vi, vj);
      faceVerts.push_back(edgeMap[{u, v}]);
    }
    result.faces.push_back(faceVerts);
  }

  // B. Create new faces at old vertices
  std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int vi = f[i];
      int vj = f[(i + 1) % f.size()];
      int u = std::min(vi, vj);
      int v = std::max(vi, vj);
      edgeToFaces[{u, v}].push_back(static_cast<int>(fi));
    }
  }

  for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
    std::vector<int> neighborMids;

    // Find a start face
    int startFaceIdx = -1;
    for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
      const auto &f = mesh.faces[fi];
      if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
        startFaceIdx = static_cast<int>(fi);
        break;
      }
    }
    if (startFaceIdx == -1)
      continue;

    int currFaceIdx = startFaceIdx;
    int safety = 0;
    do {
      const auto &face = mesh.faces[currFaceIdx];
      auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
      size_t idxInFace = std::distance(face.begin(), it);
      int nextVi = face[(idxInFace + 1) % face.size()];

      int u = std::min((int)vi, nextVi);
      int v = std::max((int)vi, nextVi);
      neighborMids.push_back(edgeMap[{u, v}]);

      // Find adjacent face
      const auto &adjFaces = edgeToFaces[{u, v}];
      auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(),
                                     [&](int id) { return id != currFaceIdx; });
      if (nextFaceIt == adjFaces.end())
        break;

      currFaceIdx = *nextFaceIt;
      safety++;
    } while (currFaceIdx != startFaceIdx && safety < 20);

    if (neighborMids.size() >= 3) {
      std::reverse(neighborMids.begin(), neighborMids.end());
      result.faces.push_back(neighborMids);
    }
  }

  normalize(result);
  return result;
}

/**
 * @brief Truncate operator: Cuts corners off the polyhedron.
 * @param t Truncation depth [0..0.5].
 */
template <typename MeshT>
static MeshT truncate(const MeshT &mesh, float t = 0.25f) {
  // Singularity check: if t is 0.5, this is geometrically equivalent to ambo
  // (rectification). The standard truncate logic produces degenerate edges at
  // t=0.5, so we redirect to ambo.
  if (std::abs(t - 0.5f) < 1e-4f) {
    return ambo(mesh);
  }

  MeshT result;
  // Map edge (u,v) -> pair of new vertex indices {near_u, near_v}
  // Stored as key {min(u,v), max(u,v)} -> value {idx_near_key_first,
  // idx_near_key_second}
  std::map<std::pair<int, int>, std::pair<int, int>> edgeMap;

  // 1. Create new vertices along edges
  for (const auto &f : mesh.faces) {
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);

      if (edgeMap.find({k1, k2}) == edgeMap.end()) {
        Vector vU = mesh.vertices[u];
        Vector vV = mesh.vertices[v];

        // Near U
        // lerp implementation in 3dmath/geometry? Using manual: a + (b-a)*t
        Vector p1 = vU + (vV - vU) * t;
        result.vertices.push_back(p1);
        int idx1 = static_cast<int>(result.vertices.size()) - 1;

        // Near V
        Vector p2 = vU + (vV - vU) * (1.0f - t);
        result.vertices.push_back(p2);
        int idx2 = static_cast<int>(result.vertices.size()) - 1;

        if (u < v) {
          edgeMap[{k1, k2}] = {idx1, idx2};
        } else {
          edgeMap[{k1, k2}] = {idx2, idx1};
        }
        if (u < v)
          edgeMap[{u, v}] = {idx1, idx2};
        else
          edgeMap[{v, u}] = {idx2, idx1};
      }
    }
  }

  // 2. Modified Faces (internal polygons)
  for (const auto &f : mesh.faces) {
    std::vector<int> faceVerts;
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);

      std::pair<int, int> indices = edgeMap[{k1, k2}];

      if (u < v) {
        faceVerts.push_back(indices.first);
        faceVerts.push_back(indices.second);
      } else {
        faceVerts.push_back(indices.second);
        faceVerts.push_back(indices.first);
      }
    }
    result.faces.push_back(faceVerts);
  }

  // 3. Corner Faces
  std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);
      edgeToFaces[{k1, k2}].push_back(static_cast<int>(fi));
    }
  }

  for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
    int startFaceIdx = -1;
    for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
      const auto &f = mesh.faces[fi];
      if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
        startFaceIdx = static_cast<int>(fi);
        break;
      }
    }
    if (startFaceIdx == -1)
      continue;

    std::vector<int> polyVerts;
    int currFaceIdx = startFaceIdx;
    int safety = 0;
    do {
      const auto &face = mesh.faces[currFaceIdx];
      auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
      size_t idxInFace = std::distance(face.begin(), it);
      int nextVi = face[(idxInFace + 1) % face.size()];

      int k1 = std::min((int)vi, nextVi);
      int k2 = std::max((int)vi, nextVi);
      std::pair<int, int> indices = edgeMap[{k1, k2}];

      int idxNearVi =
          (static_cast<int>(vi) == k1) ? indices.first : indices.second;
      polyVerts.push_back(idxNearVi);

      // Find adjacent face
      const auto &adjFaces = edgeToFaces[{k1, k2}];
      auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(),
                                     [&](int id) { return id != currFaceIdx; });
      if (nextFaceIt == adjFaces.end())
        break;

      currFaceIdx = *nextFaceIt;
      safety++;
    } while (currFaceIdx != startFaceIdx && safety < 20);

    if (polyVerts.size() > 2) {
      std::reverse(polyVerts.begin(), polyVerts.end());
      result.faces.push_back(polyVerts);
    }
  }

  normalize(result);
  return result;
}

/**
 * @brief Expand operator: Separates faces (e = aa).
 * @param t Expansion factor. Default 2-sqrt(2) ~= 0.5857.
 */
template <typename MeshT>
static MeshT expand(const MeshT &mesh, float t = 2.0f - sqrt(2.0f)) {
  MeshT result;
  // Map (faceIdx) -> list of new vertex indices
  std::vector<std::vector<int>> faceVertsMap(mesh.faces.size());

  // 1. Inset Faces
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    Vector centroid(0, 0, 0);
    for (int vi : f)
      centroid = centroid + mesh.vertices[vi];
    if (!f.empty())
      centroid = centroid / static_cast<float>(f.size());

    for (int vi : f) {
      Vector v = mesh.vertices[vi];
      Vector newV = v + (centroid - v) * t;
      result.vertices.push_back(newV);
      faceVertsMap[fi].push_back(static_cast<int>(result.vertices.size()) - 1);
    }
    result.faces.push_back(faceVertsMap[fi]);
  }

  // 2. Vertex Faces
  std::vector<std::vector<size_t>> vertToFaces(mesh.vertices.size());
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    for (int vi : mesh.faces[fi]) {
      vertToFaces[vi].push_back(fi);
    }
  }

  std::map<std::pair<int, int>, std::vector<size_t>> edgeToFacesLookup;
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);
      edgeToFacesLookup[{k1, k2}].push_back(fi);
    }
  }

  for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
    const auto &adjacentFaces = vertToFaces[vi];
    if (adjacentFaces.size() < 3)
      continue;

    std::vector<int> orderedIndices;
    size_t startFace = adjacentFaces[0];
    size_t currFace = startFace;
    int safety = 0;

    do {
      const auto &face = mesh.faces[currFace];
      auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
      int idxInFace = static_cast<int>(std::distance(face.begin(), it));

      orderedIndices.push_back(faceVertsMap[currFace][idxInFace]);

      // Prev edge: prev -> vi
      int prevVi = face[(idxInFace - 1 + face.size()) % face.size()];

      int k1 = std::min((int)vi, prevVi);
      int k2 = std::max((int)vi, prevVi);

      const auto &neighbors = edgeToFacesLookup[{k1, k2}];
      auto nextIt = std::find_if(neighbors.begin(), neighbors.end(),
                                 [&](size_t fid) { return fid != currFace; });

      if (nextIt == neighbors.end())
        break;
      currFace = *nextIt;
      safety++;
    } while (currFace != startFace && safety < 20);

    result.faces.push_back(orderedIndices);
  }

  // 3. Edge Quad Faces
  // Map key -> list of {faceIdx, indexInFace, u, v}
  struct EdgeEntry {
    size_t fi;
    int i;
    int u;
    int v;
  };
  std::map<std::pair<int, int>, std::vector<EdgeEntry>> edgeMap;
  std::vector<std::pair<int, int>>
      edgeOrder; // To preserve insertion order (JS Parity)

  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);

      if (edgeMap.find({k1, k2}) == edgeMap.end()) {
        edgeOrder.push_back({k1, k2});
      }
      edgeMap[{k1, k2}].push_back({fi, (int)i, u, v});
    }
  }

  for (const auto &key : edgeOrder) {
    const auto &entries = edgeMap[key];
    if (entries.size() != 2)
      continue;

    const auto &e1 = entries[0];
    const auto &e2 = entries[1];

    // Vertices in Face A (e1) corresponding to u, v
    int count1 = static_cast<int>(mesh.faces[e1.fi].size());
    int A_u_idx = faceVertsMap[e1.fi][e1.i];
    int A_v_idx = faceVertsMap[e1.fi][(e1.i + 1) % count1];

    // Vertices in Face B (e2) corresponding to u, v
    // Need to find u and v in Face B
    const auto &f2 = mesh.faces[e2.fi];
    auto it_v = std::find(f2.begin(), f2.end(), e1.v);
    auto it_u = std::find(f2.begin(), f2.end(), e1.u);

    int idx_v_in_B = static_cast<int>(std::distance(f2.begin(), it_v));
    int idx_u_in_B = static_cast<int>(std::distance(f2.begin(), it_u));

    int B_v_idx = faceVertsMap[e2.fi][idx_v_in_B];
    int B_u_idx = faceVertsMap[e2.fi][idx_u_in_B];

    // Quad: A_v -> A_u -> B_u -> B_v
    result.faces.push_back({A_v_idx, A_u_idx, B_u_idx, B_v_idx});
  }

  normalize(result);
  return result;
}

/**
 * @brief Bitruncate operator: Truncate the rectified mesh.
 */
template <typename MeshT>
static MeshT bitruncate(const MeshT &mesh, float t = 1.0f / 3.0f) {
  return truncate(ambo(mesh), t);
}

/**
 * @brief Canonicalize operator: Iteratively relaxes the mesh to equalize edge
 * lengths.
 */
template <typename MeshT>
static MeshT canonicalize(const MeshT &mesh, int iterations = 100) {
  MeshT result = mesh; // Copy
  std::vector<Vector> &positions = result.vertices;

  // Build adjacency
  std::vector<std::vector<int>> neighbors(positions.size());
  for (const auto &f : result.faces) {
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];

      bool found = false;
      for (int existing : neighbors[u])
        if (existing == v)
          found = true;
      if (!found)
        neighbors[u].push_back(v);

      found = false;
      for (int existing : neighbors[v])
        if (existing == u)
          found = true;
      if (!found)
        neighbors[v].push_back(u);
    }
  }

  for (int iter = 0; iter < iterations; ++iter) {
    double totalLen = 0;
    int edgeCount = 0;
    for (size_t i = 0; i < positions.size(); ++i) {
      for (int ni : neighbors[i]) {
        if ((int)i < ni) {
          totalLen += distance_between(positions[i], positions[ni]);
          edgeCount++;
        }
      }
    }
    if (edgeCount == 0)
      break;
    float targetLen = static_cast<float>(totalLen / edgeCount);

    std::vector<Vector> movements(positions.size(), Vector(0, 0, 0));

    for (size_t i = 0; i < positions.size(); ++i) {
      Vector force(0, 0, 0);
      for (int ni : neighbors[i]) {
        Vector vec = positions[ni] - positions[i];
        float dist = vec.length();
        float diff = dist - targetLen;
        // Hooke's Law
        if (dist > 1e-6f) {
          force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
        }
      }
      movements[i] = movements[i] + force;
    }

    for (size_t i = 0; i < positions.size(); ++i) {
      positions[i] = positions[i] + movements[i];
      positions[i] = positions[i].normalize();
    }
  }
  return result;
}

/**
 * @brief Snub operator: Creates a chiral semi-regular polyhedron.
 * Updated with twist support.
 */
template <typename MeshT>
static MeshT snub(const MeshT &mesh, float t = 0.5f, float twist = 0.0f) {
  MeshT result;
  std::vector<std::vector<int>> newVertsMap(mesh.faces.size());

  // 1. Create new vertices
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    Vector centroid(0, 0, 0);
    for (int vi : f)
      centroid = centroid + mesh.vertices[vi];
    if (!f.empty())
      centroid = centroid / static_cast<float>(f.size());

    Vector normal(0, 0, 0);
    // Try face normal first
    if (f.size() >= 3) {
      Vector ab = mesh.vertices[f[1]] - mesh.vertices[f[0]];
      Vector ac = mesh.vertices[f[2]] - mesh.vertices[f[0]];
      normal = cross(ab, ac).normalize();
    }
    // Fallback to centroid logic if face is degenerate or we want robust
    // spherical normal
    if (dot(centroid, centroid) > 1e-6f) {
      if (dot(normal, normal) < 1e-9f)
        normal = centroid.normalize();
    }

    newVertsMap[fi].resize(f.size());
    for (size_t i = 0; i < f.size(); ++i) {
      Vector v = mesh.vertices[f[i]];
      // lerp(v, centroid, t) -> v + (centroid - v)*t
      Vector newV = v + (centroid - v) * t;

      if (twist != 0.0f) {
        // Rotate around normal at centroid
        Vector local = newV - centroid;
        Quaternion q = make_rotation(normal, twist);
        newV = centroid + rotate(local, q);
      }

      result.vertices.push_back(newV);
      newVertsMap[fi][i] = static_cast<int>(result.vertices.size()) - 1;
    }
  }

  // 2. Face Faces
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    result.faces.push_back(newVertsMap[fi]);
  }

  // Helper: Edge Map
  std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);
      edgeToFaces[{k1, k2}].push_back(static_cast<int>(fi));
    }
  }

  // 3. Vertex Faces
  for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
    int startFaceIdx = -1;
    for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
      const auto &f = mesh.faces[fi];
      if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
        startFaceIdx = static_cast<int>(fi);
        break;
      }
    }
    if (startFaceIdx == -1)
      continue;

    std::vector<int> orderedFaces;
    int currFaceIdx = startFaceIdx;
    int safety = 0;
    do {
      orderedFaces.push_back(currFaceIdx);
      const auto &face = mesh.faces[currFaceIdx];
      auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
      size_t idxInFace = std::distance(face.begin(), it);
      int prevVi = face[(idxInFace - 1 + face.size()) % face.size()];

      int k1 = std::min((int)vi, prevVi);
      int k2 = std::max((int)vi, prevVi);
      const auto &adjFaces = edgeToFaces[{k1, k2}];
      auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(),
                                     [&](int id) { return id != currFaceIdx; });
      if (nextFaceIt == adjFaces.end())
        break;

      currFaceIdx = *nextFaceIt;
      safety++;
    } while (currFaceIdx != startFaceIdx && safety < 20);

    std::vector<int> faceVerts;
    for (int fi : orderedFaces) {
      const auto &face = mesh.faces[fi];
      auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
      size_t idx = std::distance(face.begin(), it);
      faceVerts.push_back(newVertsMap[fi][idx]);
    }
    result.faces.push_back(faceVerts);
  }

  // 4. Edge Triangles
  std::set<std::pair<int, int>> processedEdges;
  for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
    const auto &f = mesh.faces[fi];
    for (size_t i = 0; i < f.size(); ++i) {
      int vi = f[i];
      int vj = f[(i + 1) % f.size()];
      int k1 = std::min(vi, vj);
      int k2 = std::max(vi, vj);

      if (processedEdges.count({k1, k2}))
        continue;
      processedEdges.insert({k1, k2});

      const auto &adj = edgeToFaces[{k1, k2}];
      if (adj.size() < 2)
        continue;

      int faceA = fi;
      int faceB = (adj[0] == static_cast<int>(fi)) ? adj[1] : adj[0];

      auto itA_u =
          std::find(mesh.faces[faceA].begin(), mesh.faces[faceA].end(), vi);
      int idxA_u =
          static_cast<int>(std::distance(mesh.faces[faceA].begin(), itA_u));
      auto itA_v =
          std::find(mesh.faces[faceA].begin(), mesh.faces[faceA].end(), vj);
      int idxA_v =
          static_cast<int>(std::distance(mesh.faces[faceA].begin(), itA_v));

      auto itB_u =
          std::find(mesh.faces[faceB].begin(), mesh.faces[faceB].end(), vi);
      int idxB_u =
          static_cast<int>(std::distance(mesh.faces[faceB].begin(), itB_u));
      auto itB_v =
          std::find(mesh.faces[faceB].begin(), mesh.faces[faceB].end(), vj);
      int idxB_v =
          static_cast<int>(std::distance(mesh.faces[faceB].begin(), itB_v));

      int A_u = newVertsMap[faceA][idxA_u];
      int A_v = newVertsMap[faceA][idxA_v];
      int B_u = newVertsMap[faceB][idxB_u];
      int B_v = newVertsMap[faceB][idxB_v];

      // Tri 1
      result.faces.push_back({A_v, A_u, B_v});
      // Tri 2
      result.faces.push_back({B_u, B_v, A_u});
    }
  }

  normalize(result);
  return result;
}

/**
 * @brief Gyro operator: dual(snub(mesh)).
 */
template <typename MeshT> static MeshT gyro(const MeshT &mesh) {
  return dual(snub(mesh));
}

/**
 * @brief Computes KDTree and Adjacency map for the mesh (caching it).
 */
template <typename MeshT> static void compute_kdtree(const MeshT &mesh) {
  if (mesh.cache_valid)
    return;

  MeshT &m = const_cast<MeshT &>(mesh);

  // 1. Build Adjacency
  m.adjacency.assign(m.vertices.size(), {});
  for (const auto &f : m.faces) {
    for (size_t i = 0; i < f.size(); ++i) {
      int u = f[i];
      int v = f[(i + 1) % f.size()];

      bool hasV = false;
      for (int x : m.adjacency[u])
        if (x == v)
          hasV = true;
      if (!hasV)
        m.adjacency[u].push_back(v);

      bool hasU = false;
      for (int x : m.adjacency[v])
        if (x == u)
          hasU = true;
      if (!hasU)
        m.adjacency[v].push_back(u);
    }
  }

  // 2. Build KDTree (requires constructing new one)
  // We pass the span of vertices directly.
  // Ensure KDTree constructor calls build() on them.
  // KDTree stores copies of points, so safe even if vector reallocates later
  // (though it shouldn't for static mesh) Note: KDTree constructor is
  // KDTree(std::span<Vector>).
  m.kdTree = KDTree(std::span<Vector>(m.vertices));
  m.cache_valid = true;
}

/**
 * @brief Finds the closest point on the mesh graph (BFS).
 * @param p The query point.
 * @param mesh The mesh to search.
 * @return The position of the closest vertex.
 */
template <typename MeshT>
static Vector closest_point_on_mesh_graph(const Vector &p, const MeshT &mesh) {
  if (mesh.vertices.empty())
    return Vector(0, 1, 0);

  compute_kdtree(mesh);

  // 1. Closest Vertex
  auto nearestNodes = mesh.kdTree.nearest(p, 1);
  if (nearestNodes.size() == 0)
    return mesh.vertices[0];

  const auto &node = nearestNodes[0];
  int closestVertexIndex = (int)node.originalIndex;
  Vector closestVertexPos = node.point;

  Vector bestPoint = closestVertexPos;
  float maxDot = dot(p, bestPoint);

  // 2. Check connected edges
  if (closestVertexIndex < 0 ||
      closestVertexIndex >= (int)mesh.adjacency.size())
    return bestPoint;
  const auto &neighbors = mesh.adjacency[closestVertexIndex];
  if (neighbors.empty())
    return bestPoint;

  Vector A = closestVertexPos;

  for (int neighborIdx : neighbors) {
    if (neighborIdx < 0 || neighborIdx >= (int)mesh.vertices.size())
      continue;
    Vector B = mesh.vertices[neighborIdx];

    // Great circle normal
    Vector N = cross(A, B);
    float lenSq = dot(N, N);
    if (lenSq < 1e-6f)
      continue;
    N = N * (1.0f / sqrt(lenSq));

    // Project P
    float pDotN = dot(p, N);
    Vector proj = p + N * (-pDotN); // P_proj
    proj = proj.normalize();

    // Arc Check using Cross Products
    Vector crossAC = cross(A, proj);
    Vector crossCB = cross(proj, B);

    // Check if winding matches A->B (N)
    if (dot(crossAC, N) > 0 && dot(crossCB, N) > 0) {
      float d = dot(p, proj);
      if (d > maxDot) {
        maxDot = d;
        bestPoint = proj;
      }
    }
  }

  return bestPoint;
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

// JS hash32 for consistency
static uint32_t js_hash32(uint32_t n, uint32_t seed) {
  n = (n ^ seed) * 0x5bd1e995;
  n ^= n >> 15;
  return n * 0x97455bcd;
}

/**
 * @brief Colors faces based on their vertex count and neighbor topology.
 */
template <typename MeshT>
static std::vector<int> classify_faces_by_topology(const MeshT &mesh) {
  HalfEdgeMesh heMesh(mesh);
  std::map<uint32_t, int> signatureToID;

  std::vector<int> faceColorIndices;
  faceColorIndices.resize(heMesh.faces.size());
  int nextID = 0;

  // Compute Contextual Hashes (acc neighbors)
  for (size_t i = 0; i < heMesh.faces.size(); ++i) {
    HEFace *face = &heMesh.faces[i];

    // Start hash with self intrinsic (order independent seed)
    uint32_t neighborAcc = 0;

    HalfEdge *he = face->halfEdge;
    if (he) {
      HalfEdge *start = he;
      int safety = 0;
      do {
        if (he->pair && he->pair->face) {
          uint32_t h = js_hash32(he->pair->face->intrinsicHash, 0);
          neighborAcc += h;
        }
        he = he->next;
        safety++;
      } while (he && he != start && safety < 100);
    }

    // Final Hash: Combine Neighbor Accumulator with Self Intrinsic
    uint32_t finalHash = js_hash32(neighborAcc, face->intrinsicHash);

    if (signatureToID.find(finalHash) == signatureToID.end()) {
      signatureToID[finalHash] = nextID++;
    }
    faceColorIndices[i] = signatureToID[finalHash];
  }

  return faceColorIndices;
}

} // namespace MeshOps
