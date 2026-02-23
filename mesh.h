/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"
#include "spatial.h"
#include "memory.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <span>
#include <cmath>

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
  ArenaVector<Vector> vertices;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<int> faces;

  // Cache
  mutable bool cache_valid = false;
  mutable KDTree kdTree;

  PolyMesh() = default;

  void initialize(Arena &arena, size_t num_verts, size_t num_faces,
                  size_t num_indices) {
    vertices.initialize(arena, num_verts);
    face_counts.initialize(arena, num_faces);
    faces.initialize(arena, num_indices);
    cache_valid = false;
  }

  inline void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    cache_valid = false;
  }

  void clear_cache() const {
    cache_valid = false;
    kdTree.clear();
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

    // 1. Collect vertices on the stack (No malloc!)
    Vector verts[100];
    HalfEdge *he = halfEdge;
    HalfEdge *start = he;
    int safety = 0;
    do {
      if (safety < 100)
        verts[safety] = he->vertex->position;
      he = he->next;
      safety++;
    } while (he && he != start && safety < 100);

    vertexCount = safety;

    // 2. Angles
    int angles[100];
    if (vertexCount >= 3) {
      for (int i = 0; i < vertexCount; ++i) {
        const Vector &prev = verts[(i - 1 + vertexCount) % vertexCount];
        const Vector &curr = verts[i];
        const Vector &next = verts[(i + 1) % vertexCount];

        Vector v1 = (prev - curr).normalize();
        Vector v2 = (next - curr).normalize();

        float ang = angle_between(v1, v2);
        angles[i] = (int)std::round(ang * 180.0f / PI_F);
      }
      // std::sort works perfectly on raw pointers
      std::sort(angles, angles + vertexCount);
    }

    // 3. Hash (JS style hash32)
    auto hash32 = [](uint32_t n, uint32_t seed) -> uint32_t {
      n = (n ^ seed) * 0x5bd1e995;
      n ^= n >> 15;
      return n * 0x97455bcd;
    };

    uint32_t h = hash32(static_cast<uint32_t>(vertexCount), 0x12345678);
    for (int i = 0; i < vertexCount; ++i) {
      h = hash32(static_cast<uint32_t>(angles[i]), h);
    }
    intrinsicHash = h;
  }
};

class HalfEdgeMesh {
public:
  ArenaVector<HEVertex> vertices;
  ArenaVector<HEFace> faces;
  ArenaVector<HalfEdge> halfEdges;

  explicit HalfEdgeMesh(Arena &arena, const PolyMesh &mesh) {
    build_from_flat(arena, mesh.vertices, mesh.face_counts, mesh.faces);
  }

  explicit HalfEdgeMesh(Arena &arena, const MeshState &mesh) {
    build_from_flat(arena, mesh.vertices, mesh.face_counts, mesh.faces);
  }

private:
  template <typename Verts, typename Counts, typename Faces>
  void build_from_flat(Arena &arena, const Verts &verts, const Counts &counts,
                       const Faces &faces_arr) {
    size_t num_verts = verts.size();
    size_t num_faces = counts.size();
    size_t total_indices = faces_arr.size();

    vertices.initialize(arena, num_verts);
    for (size_t i = 0; i < num_verts; ++i) {
      vertices.push_back({verts[i], nullptr});
    }

    faces.initialize(arena, num_faces);
    halfEdges.initialize(arena, total_indices);

    ArenaMap<std::pair<int, int>, HalfEdge *> edgeMap(arena, total_indices);

    size_t face_offset = 0;
    size_t he_idx = 0;

    for (size_t fi = 0; fi < num_faces; ++fi) {
      int count = counts[fi];

      faces.emplace_back();
      HEFace *currentFace = &faces[faces.size() - 1];
      size_t faceStartHeIdx = he_idx;

      for (int i = 0; i < count; ++i) {
        halfEdges.emplace_back();
      }

      for (int i = 0; i < count; ++i) {
        int u = faces_arr[face_offset + i];
        int v = faces_arr[face_offset + (i + 1) % count];

        HalfEdge *he = &halfEdges[faceStartHeIdx + i];
        he->vertex = &vertices[v];
        he->face = currentFace;
        he->next = &halfEdges[faceStartHeIdx + (i + 1) % count];
        he->prev = &halfEdges[faceStartHeIdx + (i - 1 + count) % count];

        vertices[v].halfEdge = he;

        if (edgeMap.contains({v, u})) {
          HalfEdge *neighbor = edgeMap[{v, u}];
          he->pair = neighbor;
          neighbor->pair = he;
        } else {
          edgeMap[{u, v}] = he;
        }
        he_idx++;
      }
      currentFace->halfEdge = &halfEdges[faceStartHeIdx];
      face_offset += count;
    }

    for (size_t fi = 0; fi < num_faces; ++fi)
      faces[fi].compute_properties();
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
  ArenaVector<Vector> staticVertices;
  ArenaVector<Vector> dynamicVertices;
  ArenaVector<HankinInstruction> dynamicInstructions;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<int> faces;
  int staticOffset;
};

/**
 * @brief Operations on meshes (Dual, Hankin, etc.).
 */
namespace MeshOps {

/**
 * @brief Compiles a PolyMesh into a static MeshState.
 * Removes degenerate faces (faces with < 3 vertices) during
 * the process.
 * @param src The source PolyMesh.
 * @return The compiled MeshState.
 */
inline void compile(const PolyMesh &src, MeshState &dst, Arena &geom_arena) {
  dst.clear();

  size_t valid_faces = 0;
  size_t valid_indices = 0;

  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    if (src.face_counts[i] >= 3) {
      valid_faces++;
      valid_indices += src.face_counts[i];
    }
  }

  dst.vertices.initialize(geom_arena, src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    dst.vertices.push_back(src.vertices[i]);
  }

  dst.face_counts.initialize(geom_arena, valid_faces);
  dst.faces.initialize(geom_arena, valid_indices);
  dst.face_offsets.initialize(geom_arena, valid_faces);

  size_t offset = 0;
  int current_offset = 0;
  for (size_t i = 0; i < src.face_counts.size(); ++i) {
    int count = src.face_counts[i];
    if (count >= 3) {
      dst.face_counts.push_back(static_cast<uint8_t>(count));
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
 * @brief Transforms a MeshState into a target MeshState using the provided
 * transformer.
 * @tparam TransformerType The transformer (e.g., RippleTransformer,
 * OrientTransformer).
 * @param local_state The source mesh state (centered around origin).
 * @param world_state The destination mesh state to populate.
 * @param transformer The transformer providing a transform(Vector) method.
 */
template <typename MeshT, typename... TransformerTypes>
inline void transform(const MeshT &local_state, MeshT &world_state,
                      Arena &arena, const TransformerTypes &...transformers) {
  world_state.vertices.initialize(arena, local_state.vertices.size());
  world_state.face_counts.initialize(arena, local_state.face_counts.size());
  world_state.faces.initialize(arena, local_state.faces.size());

  if constexpr (requires { world_state.face_offsets; }) {
    world_state.face_offsets.initialize(arena, local_state.face_counts.size());
  }

  for (size_t i = 0; i < local_state.vertices.size(); ++i) {
    Vector v = local_state.vertices[i];
    if constexpr (sizeof...(transformers) > 0) {
      ((v = transformers.transform(v)), ...);
    }
    world_state.vertices.push_back(v);
  }
  for (size_t i = 0; i < local_state.face_counts.size(); ++i) {
    world_state.face_counts.push_back(local_state.face_counts[i]);
    if constexpr (requires { world_state.face_offsets; }) {
      world_state.face_offsets.push_back(local_state.face_offsets[i]);
    }
  }
  for (size_t i = 0; i < local_state.faces.size(); ++i) {
    world_state.faces.push_back(local_state.faces[i]);
  }
}

/**
 * @brief Computes the dual of a mesh.
 */
inline void dual(const PolyMesh &mesh, PolyMesh &out_mesh,
                 Arena &scratch_arena_a) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(scratch_arena_a, F);
  out_mesh.face_counts.initialize(scratch_arena_a, V);
  out_mesh.faces.initialize(scratch_arena_a, I);

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  ArenaMap<HEFace *, int> faceToVertIdx(scratch_arena_a, F);

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
    } while (he != start && count < 100);

    c = c / static_cast<float>(count);
    out_mesh.vertices.push_back(c.normalize());
    faceToVertIdx[face] = static_cast<int>(i);
  }

  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;

    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    int face_idx = out_mesh.faces.size();
    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int new_face_count = 0;
    int local_face[100];
    do {
      if (!curr->face)
        break;
      if (new_face_count < 100)
        local_face[new_face_count++] = faceToVertIdx[curr->face];

      if (!curr->prev || !curr->prev->pair)
        break;
      curr = curr->prev->pair;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (new_face_count >= 3) {
      out_mesh.face_counts.push_back(new_face_count);
      for (int k = 0; k < new_face_count; ++k) {
        out_mesh.faces.push_back(local_face[k]);
      }
    }
  }
}

inline PolyMesh dual(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out;
  dual(mesh, out, *(ctx.target));
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Compiles the topology for a Hankin pattern.
 */
inline void compile_hankin(const PolyMesh &mesh, CompiledHankin &compiled,
                           Arena &geom_arena, Arena &scratch_arena_a) {
  HalfEdgeMesh heMesh(scratch_arena_a, mesh);

  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  compiled.staticVertices.initialize(geom_arena, I);
  compiled.dynamicVertices.initialize(geom_arena, I);
  compiled.dynamicInstructions.initialize(geom_arena, I);

  compiled.face_counts.initialize(geom_arena, F + V);
  compiled.faces.initialize(geom_arena, 4 * I);

  ArenaMap<HalfEdge *, int> heToMidpointIdx(scratch_arena_a, I);
  ArenaMap<HalfEdge *, int> heToDynamicIdx(scratch_arena_a, I);

  auto getMidpointIdx = [&](HalfEdge *he) {
    if (heToMidpointIdx.contains(he))
      return heToMidpointIdx[he];
    if (he->pair && heToMidpointIdx.contains(he->pair))
      return heToMidpointIdx[he->pair];

    Vector pA =
        he->prev ? he->prev->vertex->position : he->pair->vertex->position;
    Vector pB = he->vertex->position;
    Vector mid = (pA + pB) * 0.5f;
    mid = mid.normalize();

    compiled.staticVertices.push_back(mid);
    int idx = static_cast<int>(compiled.staticVertices.size()) - 1;
    heToMidpointIdx[he] = idx;
    if (he->pair)
      heToMidpointIdx[he->pair] = idx;
    return idx;
  };

  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    getMidpointIdx(&heMesh.halfEdges[i]);
  }

  compiled.staticOffset = static_cast<int>(compiled.staticVertices.size());

  // Star faces
  for (size_t i = 0; i < heMesh.faces.size(); ++i) {
    HEFace &face = heMesh.faces[i];
    HalfEdge *he = face.halfEdge;
    HalfEdge *startHe = he;
    int count = 0;

    do {
      count += 2;
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
      compiled.dynamicVertices.emplace_back();

      compiled.faces.push_back(idxM1);
      compiled.faces.push_back(compiled.staticOffset + dynIdx);

      he = he->next;
    } while (he != startHe);

    compiled.face_counts.push_back(count);
  }

  // Rosette faces
  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;
    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int count = 0;
    int face_indices[100];

    do {
      if (count < 100)
        face_indices[count++] = heToMidpointIdx[curr];
      HalfEdge *nextEdge = curr->pair ? curr->pair->next : nullptr;
      if (!nextEdge)
        break;
      if (count < 100)
        face_indices[count++] =
            compiled.staticOffset + heToDynamicIdx[nextEdge];
      curr = nextEdge;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (count > 2) {
      compiled.face_counts.push_back(count);
      for (int k = count - 1; k >= 0; --k) {
        compiled.faces.push_back(face_indices[k]);
      }
    }
  }
}

template <typename MeshT>
inline void update_hankin(CompiledHankin &compiled, MeshT &out_mesh,
                          Arena &scratch_arena_a, float angle) {
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
    if (lenSq < 1e-6f)
      intersect = (m1 + m2).normalize();
    if (dot(intersect, instr.pCorner) < 0)
      intersect = -intersect;

    compiled.dynamicVertices[i] = intersect.normalize();
  }

  out_mesh.vertices.initialize(scratch_arena_a,
                               compiled.staticVertices.size() +
                                   compiled.dynamicVertices.size());
  for (size_t i = 0; i < compiled.staticVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.staticVertices[i]);
  for (size_t i = 0; i < compiled.dynamicVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.dynamicVertices[i]);

  out_mesh.face_counts.initialize(scratch_arena_a, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.initialize(scratch_arena_a,
                                     compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.initialize(scratch_arena_a, compiled.faces.size());
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    out_mesh.faces.push_back(compiled.faces[i]);
}

inline void hankin(const PolyMesh &mesh, PolyMesh &out_mesh,
                   Arena &scratch_arena_a, float angle) {
  CompiledHankin compiled;
  compile_hankin(mesh, compiled, scratch_arena_a, scratch_arena_a);
  update_hankin(compiled, out_mesh, scratch_arena_a, angle);
}

inline PolyMesh hankin(const PolyMesh &mesh, ScratchContext &ctx, float angle) {
  PolyMesh out;
  hankin(mesh, out, *(ctx.target), angle);
  ctx.swap_and_clear();
  return out;
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
inline void kis(const PolyMesh &mesh, PolyMesh &out_mesh,
                Arena &scratch_arena_a) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(scratch_arena_a, V + F);
  out_mesh.face_counts.initialize(scratch_arena_a, I);
  out_mesh.faces.initialize(scratch_arena_a, 3 * I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);

  size_t offset = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    int count = mesh.face_counts[fi];
    Vector centroid(0, 0, 0);
    for (int i = 0; i < count; ++i) {
      centroid = centroid + mesh.vertices[mesh.faces[offset + i]];
    }
    if (count > 0)
      centroid = centroid / static_cast<float>(count);

    out_mesh.vertices.push_back(centroid);
    int centerIdx = static_cast<int>(out_mesh.vertices.size()) - 1;

    for (int i = 0; i < count; ++i) {
      int vi = mesh.faces[offset + i];
      int vj = mesh.faces[offset + (i + 1) % count];
      out_mesh.face_counts.push_back(3);
      out_mesh.faces.push_back(vi);
      out_mesh.faces.push_back(vj);
      out_mesh.faces.push_back(centerIdx);
    }
    offset += count;
  }
}

inline PolyMesh kis(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out;
  kis(mesh, out, *(ctx.target));
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Ambo operator: Truncates vertices to edge midpoints.
 */
inline void ambo(const PolyMesh &mesh, PolyMesh &out_mesh,
                 Arena &scratch_arena_a) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(scratch_arena_a, E);
  out_mesh.face_counts.initialize(scratch_arena_a, F + V);
  out_mesh.faces.initialize(scratch_arena_a, 2 * I);

  ArenaMap<std::pair<int, int>, int> edgeMap(scratch_arena_a, E);

  size_t offset = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    int count = mesh.face_counts[fi];
    for (int i = 0; i < count; ++i) {
      int vi = mesh.faces[offset + i];
      int vj = mesh.faces[offset + (i + 1) % count];
      int u = std::min(vi, vj);
      int v = std::max(vi, vj);

      if (!edgeMap.contains({u, v})) {
        Vector mid = (mesh.vertices[vi] + mesh.vertices[vj]) * 0.5f;
        out_mesh.vertices.push_back(mid);
        edgeMap[{u, v}] = static_cast<int>(out_mesh.vertices.size()) - 1;
      }
    }
    offset += count;
  }

  offset = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    int count = mesh.face_counts[fi];
    out_mesh.face_counts.push_back(count);
    for (int i = 0; i < count; ++i) {
      int vi = mesh.faces[offset + i];
      int vj = mesh.faces[offset + (i + 1) % count];
      int u = std::min(vi, vj);
      int v = std::max(vi, vj);
      out_mesh.faces.push_back(edgeMap[{u, v}]);
    }
    offset += count;
  }

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);

  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;

    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    int face_idx = out_mesh.faces.size();
    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int count = 0;
    int local_face[100];

    do {
      if (!curr->face)
        break;
      int vi = curr->prev->vertex - heMesh.vertices.data();
      int vj = curr->vertex - heMesh.vertices.data();
      if (count < 100)
        local_face[count++] = edgeMap[{std::min(vi, vj), std::max(vi, vj)}];

      if (!curr->prev || !curr->prev->pair)
        break;
      curr = curr->prev->pair;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (count >= 3) {
      out_mesh.face_counts.push_back(count);
      for (int k = 0; k < count; ++k)
        out_mesh.faces.push_back(local_face[k]);
    }
  }
}

inline PolyMesh ambo(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out;
  ambo(mesh, out, *(ctx.target));
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Truncate operator: Cuts corners off the polyhedron.
 * @param t Truncation depth [0..0.5].
 */
inline void truncate(const PolyMesh &mesh, PolyMesh &out_mesh,
                     Arena &scratch_arena_a, float t = 0.25f) {
  if (std::abs(t - 0.5f) < 1e-4f) {
    ambo(mesh, out_mesh, scratch_arena_a);
    return;
  }

  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(scratch_arena_a, 2 * E);
  out_mesh.face_counts.initialize(scratch_arena_a, F + V);
  out_mesh.faces.initialize(scratch_arena_a, 3 * I);

  ArenaMap<std::pair<int, int>, std::pair<int, int>> edgeMap(scratch_arena_a, E);

  size_t offset = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    int count = mesh.face_counts[fi];
    for (int i = 0; i < count; ++i) {
      int u = mesh.faces[offset + i];
      int v = mesh.faces[offset + (i + 1) % count];
      int k1 = std::min(u, v);
      int k2 = std::max(u, v);

      if (!edgeMap.contains({k1, k2})) {
        Vector new_u =
            mesh.vertices[k1] + (mesh.vertices[k2] - mesh.vertices[k1]) * t;
        Vector new_v =
            mesh.vertices[k2] + (mesh.vertices[k1] - mesh.vertices[k2]) * t;

        out_mesh.vertices.push_back(new_u);
        int idx_u = static_cast<int>(out_mesh.vertices.size()) - 1;

        out_mesh.vertices.push_back(new_v);
        int idx_v = static_cast<int>(out_mesh.vertices.size()) - 1;

        edgeMap[{k1, k2}] = {idx_u, idx_v};
      }
    }
    offset += count;
  }

  offset = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    int count = mesh.face_counts[fi];
    out_mesh.face_counts.push_back(count * 2);

    for (int i = 0; i < count; ++i) {
      int vi = mesh.faces[offset + i];
      int vj = mesh.faces[offset + (i + 1) % count];

      int k1 = std::min(vi, vj);
      int k2 = std::max(vi, vj);

      std::pair<int, int> newVerts = edgeMap[{k1, k2}];

      if (vi == k1) {
        out_mesh.faces.push_back(newVerts.first);
        out_mesh.faces.push_back(newVerts.second);
      } else {
        out_mesh.faces.push_back(newVerts.second);
        out_mesh.faces.push_back(newVerts.first);
      }
    }
    offset += count;
  }

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);

  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;

    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int count = 0;
    int local_face[100];

    do {
      if (!curr->face)
        break;
      int vi = curr->prev->vertex - heMesh.vertices.data();
      int vj = curr->vertex - heMesh.vertices.data();

      int k1 = std::min(vi, vj);
      int k2 = std::max(vi, vj);
      std::pair<int, int> newVerts = edgeMap[{k1, k2}];

      if (count < 100) {
        local_face[count++] = (vi == k1) ? newVerts.first : newVerts.second;
      }

      if (!curr->prev || !curr->prev->pair)
        break;
      curr = curr->prev->pair;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (count >= 3) {
      out_mesh.face_counts.push_back(count);
      for (int k = 0; k < count; ++k)
        out_mesh.faces.push_back(local_face[k]);
    }
  }
}

inline PolyMesh truncate(const PolyMesh &mesh, ScratchContext &ctx,
                         float t = 0.25f) {
  PolyMesh out;
  truncate(mesh, out, *(ctx.target), t);
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Expand operator: Separates faces (e = aa).
 * @param t Expansion factor. Default 2-sqrt(2) ~= 0.5857.
 */
inline void expand(const PolyMesh &mesh, PolyMesh &out_mesh,
                   Arena &scratch_arena_a, float t = 2.0f - sqrt(2.0f)) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(scratch_arena_a, I);
  out_mesh.face_counts.initialize(scratch_arena_a, F + V + E);
  out_mesh.faces.initialize(scratch_arena_a, 4 * I);

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  ArenaMap<HalfEdge *, int> heToVertIdx(scratch_arena_a, I);

  for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
    HEFace *face = &heMesh.faces[fi];
    Vector centroid(0, 0, 0);
    int count = 0;
    HalfEdge *he = face->halfEdge;
    HalfEdge *start = he;
    do {
      centroid = centroid + he->vertex->position;
      count++;
      he = he->next;
    } while (he != start && count < 100);

    centroid = centroid / static_cast<float>(count);

    out_mesh.face_counts.push_back(count);
    he = start;
    int safety = 0;
    do {
      Vector v = he->vertex->position;
      Vector newV = v + (centroid - v) * t;
      out_mesh.vertices.push_back(newV);
      int idx = static_cast<int>(out_mesh.vertices.size()) - 1;
      heToVertIdx[he] = idx;

      out_mesh.faces.push_back(idx);
      he = he->next;
      safety++;
    } while (he != start && safety < 100);
  }

  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;

    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int count = 0;
    int local_face[100];

    do {
      if (!curr->face)
        break;
      if (count < 100)
        local_face[count++] = heToVertIdx[curr->prev];

      if (!curr->pair)
        break;
      curr = curr->pair->next;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (count >= 3) {
      out_mesh.face_counts.push_back(count);
      for (int k = count - 1; k >= 0; --k)
        out_mesh.faces.push_back(local_face[k]);
    }
  }

  ArenaMap<std::pair<int, int>, bool> visitedEdges(scratch_arena_a, E);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *he = &heMesh.halfEdges[i];
    int u = he->prev->vertex - heMesh.vertices.data();
    int v = he->vertex - heMesh.vertices.data();
    int k1 = std::min(u, v);
    int k2 = std::max(u, v);

    if (visitedEdges.contains({k1, k2}))
      continue;
    visitedEdges[{k1, k2}] = true;

    if (he->pair) {
      out_mesh.face_counts.push_back(4);
      out_mesh.faces.push_back(heToVertIdx[he->prev]);
      out_mesh.faces.push_back(heToVertIdx[he->pair]);
      out_mesh.faces.push_back(heToVertIdx[he->pair->prev]);
      out_mesh.faces.push_back(heToVertIdx[he]);
    }
  }
}

inline PolyMesh expand(const PolyMesh &mesh, ScratchContext &ctx,
                       float t = 2.0f - sqrt(2.0f)) {
  PolyMesh out;
  expand(mesh, out, *(ctx.target), t);
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Bitruncate operator: Truncate the rectified mesh.
 */
inline void bitruncate(const PolyMesh &mesh, PolyMesh &out_mesh,
                       Arena &scratch_arena_a, float t = 1.0f / 3.0f) {
  PolyMesh temp;
  ambo(mesh, temp, scratch_arena_a);
  truncate(temp, out_mesh, scratch_arena_a, t);
}

inline PolyMesh bitruncate(const PolyMesh &mesh, ScratchContext &ctx,
                           float t = 1.0f / 3.0f) {
  PolyMesh out;
  bitruncate(mesh, out, *(ctx.target), t);
  ctx.swap_and_clear();
  return out;
}

inline void canonicalize(const PolyMesh &mesh, PolyMesh &out_mesh,
                         Arena &scratch_arena_a, int iterations = 100) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  out_mesh.vertices.initialize(scratch_arena_a, V);
  out_mesh.face_counts.initialize(scratch_arena_a, F);
  out_mesh.faces.initialize(scratch_arena_a, I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);
  for (size_t i = 0; i < F; ++i)
    out_mesh.face_counts.push_back(mesh.face_counts[i]);
  for (size_t i = 0; i < I; ++i)
    out_mesh.faces.push_back(mesh.faces[i]);

  ArenaVector<Vector> movements;
  movements.initialize(scratch_arena_a, V);
  for (size_t i = 0; i < V; ++i)
    movements.push_back(Vector(0, 0, 0));

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);

  for (int iter = 0; iter < iterations; ++iter) {
    double totalLen = 0;
    int edgeCount = 0;

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      HalfEdge *he = &heMesh.halfEdges[i];
      int u = he->prev->vertex - heMesh.vertices.data();
      int v = he->vertex - heMesh.vertices.data();
      if (u < v) {
        totalLen +=
            distance_between(out_mesh.vertices[u], out_mesh.vertices[v]);
        edgeCount++;
      }
    }

    if (edgeCount == 0)
      break;
    float targetLen = static_cast<float>(totalLen / edgeCount);

    for (size_t i = 0; i < V; ++i) {
      Vector force(0, 0, 0);

      HEVertex *hev = &heMesh.vertices[i];
      HalfEdge *he = hev->halfEdge;
      if (he) {
        HalfEdge *start = he;
        int safety = 0;
        do {
          int ni = he->prev->vertex - heMesh.vertices.data();
          Vector vec = out_mesh.vertices[ni] - out_mesh.vertices[i];
          float dist = vec.length();
          if (dist > 1e-6f) {
            float diff = dist - targetLen;
            force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
          }

          if (!he->next || !he->next->pair)
            break;
          he = he->next->pair;
          safety++;
        } while (he != start && safety < 100);
      }
      movements[i] = force;
    }

    for (size_t i = 0; i < V; ++i) {
      out_mesh.vertices[i] = (out_mesh.vertices[i] + movements[i]).normalize();
    }
  }
}

inline PolyMesh canonicalize(const PolyMesh &mesh, ScratchContext &ctx,
                             int iterations = 100) {
  PolyMesh out;
  canonicalize(mesh, out, *(ctx.target), iterations);
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Snub operator: Creates a chiral semi-regular polyhedron.
 * Updated with twist support.
 */
inline void snub(const PolyMesh &mesh, PolyMesh &out_mesh, Arena &scratch_arena_a,
                 float t = 0.5f, float twist = 0.0f) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();
  size_t E = I / 2;

  out_mesh.vertices.initialize(scratch_arena_a, I);
  out_mesh.face_counts.initialize(scratch_arena_a, F + V + 2 * E);
  out_mesh.faces.initialize(scratch_arena_a, 5 * I);

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  ArenaMap<HalfEdge *, int> heToVertIdx(scratch_arena_a, I);

  for (size_t fi = 0; fi < heMesh.faces.size(); ++fi) {
    HEFace *face = &heMesh.faces[fi];
    Vector centroid(0, 0, 0);
    int count = 0;
    HalfEdge *he = face->halfEdge;
    HalfEdge *start = he;
    do {
      centroid = centroid + he->vertex->position;
      count++;
      he = he->next;
    } while (he != start && count < 100);

    centroid = centroid / static_cast<float>(count);

    Vector normal(0, 0, 0);
    if (count >= 3) {
      Vector ab = start->next->vertex->position - start->vertex->position;
      Vector ac = start->next->next->vertex->position - start->vertex->position;
      normal = cross(ab, ac).normalize();
    }
    if (dot(centroid, centroid) > 1e-6f && dot(normal, normal) < 1e-9f) {
      normal = centroid.normalize();
    }

    out_mesh.face_counts.push_back(count);
    he = start;
    int safety = 0;
    do {
      Vector v = he->vertex->position;
      Vector newV = v + (centroid - v) * t;

      if (twist != 0.0f) {
        Vector local = newV - centroid;
        Quaternion q = make_rotation(normal, twist);
        newV = centroid + rotate(local, q);
      }

      out_mesh.vertices.push_back(newV);
      int idx = static_cast<int>(out_mesh.vertices.size()) - 1;
      heToVertIdx[he] = idx;

      out_mesh.faces.push_back(idx);

      he = he->next;
      safety++;
    } while (he != start && safety < 100);
  }

  ArenaMap<HEVertex *, bool> visitedVerts(scratch_arena_a, V);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *heStart = &heMesh.halfEdges[i];
    if (!heStart->prev)
      continue;

    HEVertex *origin = heStart->prev->vertex;
    if (visitedVerts.contains(origin))
      continue;
    visitedVerts[origin] = true;

    HalfEdge *curr = heStart;
    HalfEdge *startOrbit = curr;
    int safety = 0;
    int count = 0;
    int local_face[100];

    do {
      if (!curr->face)
        break;
      if (count < 100)
        local_face[count++] = heToVertIdx[curr->prev];

      if (!curr->pair)
        break;
      curr = curr->pair->next;
      safety++;
    } while (curr != startOrbit && curr && safety < 100);

    if (count >= 3) {
      out_mesh.face_counts.push_back(count);
      for (int k = count - 1; k >= 0; --k)
        out_mesh.faces.push_back(local_face[k]);
    }
  }

  ArenaMap<std::pair<int, int>, bool> visitedEdges(scratch_arena_a, E);
  for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
    HalfEdge *he = &heMesh.halfEdges[i];
    int u = he->prev->vertex - heMesh.vertices.data();
    int v = he->vertex - heMesh.vertices.data();
    int k1 = std::min(u, v);
    int k2 = std::max(u, v);

    if (visitedEdges.contains({k1, k2}))
      continue;
    visitedEdges[{k1, k2}] = true;

    if (he->pair) {
      out_mesh.face_counts.push_back(3);
      out_mesh.faces.push_back(heToVertIdx[he->prev]);
      out_mesh.faces.push_back(heToVertIdx[he->pair]);
      out_mesh.faces.push_back(heToVertIdx[he]);

      out_mesh.face_counts.push_back(3);
      out_mesh.faces.push_back(heToVertIdx[he->pair]);
      out_mesh.faces.push_back(heToVertIdx[he->pair->prev]);
      out_mesh.faces.push_back(heToVertIdx[he]);
    }
  }
}

inline PolyMesh snub(const PolyMesh &mesh, ScratchContext &ctx, float t = 0.5f,
                     float twist = 0.0f) {
  PolyMesh out;
  snub(mesh, out, *(ctx.target), t, twist);
  ctx.swap_and_clear();
  return out;
}

inline void gyro(const PolyMesh &mesh, PolyMesh &out_mesh,
                 Arena &scratch_arena_a) {
  PolyMesh temp;
  snub(mesh, temp, scratch_arena_a);
  dual(temp, out_mesh, scratch_arena_a);
}

inline PolyMesh gyro(const PolyMesh &mesh, ScratchContext &ctx) {
  PolyMesh out;
  gyro(mesh, out, *(ctx.target));
  ctx.swap_and_clear();
  return out;
}

/**
 * @brief Computes KDTree and Adjacency map for the mesh (caching it).
 */
inline void compute_kdtree(const PolyMesh &mesh, Arena &arena) {
  if (mesh.cache_valid)
    return;

  PolyMesh &m = const_cast<PolyMesh &>(mesh);
  m.kdTree =
      KDTree(arena, std::span<Vector>(m.vertices.data(), m.vertices.size()));
  m.cache_valid = true;
}

inline Vector closest_point_on_mesh_graph(const Vector &p, const PolyMesh &mesh,
                                          Arena &scratch_arena_a) {
  if (mesh.vertices.empty())
    return Vector(0, 1, 0);

  compute_kdtree(mesh, scratch_arena_a);

  auto nearestNodes = mesh.kdTree.nearest(p, 1);
  if (nearestNodes.size() == 0)
    return mesh.vertices[0];

  const auto &node = nearestNodes[0];
  int closestVertexIndex = (int)node.originalIndex;
  Vector closestVertexPos = node.point;

  Vector bestPoint = closestVertexPos;
  float maxDot = dot(p, bestPoint);

  HalfEdgeMesh heMesh(scratch_arena_a, mesh);
  if (closestVertexIndex < 0 ||
      closestVertexIndex >= (int)heMesh.vertices.size())
    return bestPoint;

  HEVertex *hev = &heMesh.vertices[closestVertexIndex];
  HalfEdge *he = hev->halfEdge;
  if (!he)
    return bestPoint;

  Vector A = closestVertexPos;
  HalfEdge *start = he;
  int safety = 0;

  do {
    int neighborIdx = he->prev->vertex - heMesh.vertices.data();
    Vector B = mesh.vertices[neighborIdx];

    Vector N = cross(A, B);
    float lenSq = dot(N, N);
    if (lenSq >= 1e-6f) {
      N = N * (1.0f / sqrt(lenSq));
      float pDotN = dot(p, N);
      Vector proj = p + N * (-pDotN);
      proj = proj.normalize();

      Vector crossAC = cross(A, proj);
      Vector crossCB = cross(proj, B);

      if (dot(crossAC, N) > 0 && dot(crossCB, N) > 0) {
        float d = dot(p, proj);
        if (d > maxDot) {
          maxDot = d;
          bestPoint = proj;
        }
      }
    }

    if (!he->pair)
      break;
    he = he->pair->next;
    safety++;
  } while (he != start && safety < 100);

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
static std::vector<int> classify_faces_by_topology(const MeshT &mesh,
                                                   Arena &scratch_arena_a) {
  HalfEdgeMesh heMesh(scratch_arena_a, mesh);

  // Replace std::map with ArenaMap. Max possible signatures = faces.size()
  ArenaMap<uint32_t, int> signatureToID(scratch_arena_a, heMesh.faces.size());

  // Use std::vector so the persistent caller object preserves values out of
  // scope
  std::vector<int> faceColorIndices;
  faceColorIndices.reserve(heMesh.faces.size());
  int nextID = 0;

  for (size_t i = 0; i < heMesh.faces.size(); ++i) {
    HEFace *face = &heMesh.faces[i];
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

    uint32_t finalHash = js_hash32(neighborAcc, face->intrinsicHash);

    if (!signatureToID.contains(finalHash)) {
      signatureToID[finalHash] = nextID++;
    }
    faceColorIndices.push_back(signatureToID[finalHash]);
  }

  return faceColorIndices;
}

} // namespace MeshOps
