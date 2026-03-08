/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh.h"

/**
 * @brief Structure returned by compile_hankin.
 */
struct HankinInstruction {
  uint16_t vCorner; /**< Index to baseVertices for corner vertex. */
  uint16_t vPrev;   /**< Index to baseVertices for previous vertex. */
  uint16_t vNext;   /**< Index to baseVertices for next vertex. */
  uint16_t idxM1;   /**< Index of first midpoint (static vertex). */
  uint16_t idxM2;   /**< Index of second midpoint (static vertex). */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 */
struct CompiledHankin {
  ArenaVector<Vector> baseVertices;
  ArenaVector<Vector> staticVertices;
  ArenaVector<Vector> dynamicVertices;
  ArenaVector<HankinInstruction> dynamicInstructions;
  ArenaVector<uint8_t> face_counts;
  ArenaVector<uint16_t> faces;
  int staticOffset;
};

/**
 * @brief Operations on meshes (Hankin pattern operators).
 */
namespace MeshOps {

/**
 * @brief Compiles the topology for a Hankin pattern.
 */
template <typename MeshT>
FLASHMEM inline void compile_hankin(const MeshT &mesh, CompiledHankin &compiled,
                                    Arena &target_arena, Arena &temp_arena) {
  size_t V = mesh.vertices.size();
  size_t F = mesh.face_counts.size();
  size_t I = mesh.faces.size();

  compiled.baseVertices.initialize(target_arena, V);
  for (size_t i = 0; i < V; ++i) {
    compiled.baseVertices.push_back(mesh.vertices[i]);
  }
  compiled.staticVertices.initialize(target_arena, (I / 2) + 1);
  compiled.dynamicVertices.initialize(target_arena, I);
  compiled.dynamicInstructions.initialize(target_arena, I);
  compiled.face_counts.initialize(target_arena, F + V);
  compiled.faces.initialize(target_arena, 3 * I);

  {
    ScopedScratch _(temp_arena);

    HalfEdgeMesh heMesh(temp_arena, mesh);
    uint16_t *heToMidpointIdx = static_cast<uint16_t *>(
        temp_arena.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToMidpointIdx, I, HE_NONE);
    uint16_t *heToDynamicIdx = static_cast<uint16_t *>(
        temp_arena.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(heToDynamicIdx, I, HE_NONE);

    auto getMidpointIdx = [&](uint16_t heIdx) {
      if (heToMidpointIdx[heIdx] != HE_NONE)
        return heToMidpointIdx[heIdx];
      HalfEdge &he = heMesh.halfEdges[heIdx];
      if (he.pair != HE_NONE && heToMidpointIdx[he.pair] != HE_NONE)
        return heToMidpointIdx[he.pair];

      Vector pA = he.prev != HE_NONE
                      ? mesh.vertices[heMesh.halfEdges[he.prev].vertex]
                      : mesh.vertices[heMesh.halfEdges[he.pair].vertex];
      Vector pB = mesh.vertices[he.vertex];
      Vector mid = (pA + pB) * 0.5f;
      mid = mid.normalize();

      compiled.staticVertices.push_back(mid);
      uint16_t idx = static_cast<uint16_t>(compiled.staticVertices.size() - 1);
      heToMidpointIdx[heIdx] = idx;
      if (he.pair != HE_NONE)
        heToMidpointIdx[he.pair] = idx;
      return idx;
    };

    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      getMidpointIdx(static_cast<uint16_t>(i));
    }

    compiled.staticOffset = static_cast<int>(compiled.staticVertices.size());

    // Star faces
    for (size_t i = 0; i < heMesh.faces.size(); ++i) {
      HEFace &face = heMesh.faces[i];
      uint16_t heIdx = face.halfEdge;
      uint16_t startHe = heIdx;
      int count = 0;

      if (heIdx == HE_NONE)
        continue;

      do {
        count += 2;
        HalfEdge &currHe = heMesh.halfEdges[heIdx];
        uint16_t prevIdx = currHe.prev;

        int idxM1 = getMidpointIdx(prevIdx);
        int idxM2 = getMidpointIdx(heIdx);

        HalfEdge &prevHe = heMesh.halfEdges[prevIdx];

        uint16_t iCorner = prevHe.vertex;
        uint16_t iPrev = prevHe.prev != HE_NONE
                             ? heMesh.halfEdges[prevHe.prev].vertex
                             : heMesh.halfEdges[prevHe.pair].vertex;
        uint16_t iNext = currHe.vertex;

        compiled.dynamicInstructions.push_back({iCorner, iPrev, iNext,
                                                static_cast<uint16_t>(idxM1),
                                                static_cast<uint16_t>(idxM2)});

        int16_t dynIdx = static_cast<int16_t>(compiled.dynamicVertices.size());
        heToDynamicIdx[heIdx] = dynIdx;
        compiled.dynamicVertices.emplace_back();

        compiled.faces.push_back(idxM1);
        compiled.faces.push_back(compiled.staticOffset + dynIdx);

        heIdx = currHe.next;
      } while (heIdx != HE_NONE && heIdx != startHe);

      compiled.face_counts.push_back(count);
    }

    // Rosette faces
    bool *visitedVerts = static_cast<bool *>(
        temp_arena.allocate(V * sizeof(bool), alignof(bool)));
    std::fill_n(visitedVerts, V, false);
    for (size_t i = 0; i < heMesh.halfEdges.size(); ++i) {
      uint16_t heStartIdx = static_cast<uint16_t>(i);
      HalfEdge &heStart = heMesh.halfEdges[heStartIdx];
      if (heStart.prev == HE_NONE)
        continue;
      uint16_t originIdx = heMesh.halfEdges[heStart.prev].vertex;
      if (visitedVerts[originIdx])
        continue;
      visitedVerts[originIdx] = true;

      uint16_t currIdx = heStartIdx;
      uint16_t startOrbit = currIdx;
      int safety = 0;
      int count = 0;
      int16_t face_indices[100];

      do {
        HalfEdge &currHe = heMesh.halfEdges[currIdx];
        if (count < 100)
          face_indices[count++] = heToMidpointIdx[currIdx];
        uint16_t nextEdgeIdx = currHe.pair != HE_NONE
                                   ? heMesh.halfEdges[currHe.pair].next
                                   : HE_NONE;
        if (nextEdgeIdx == HE_NONE)
          break;
        if (count < 100)
          face_indices[count++] =
              compiled.staticOffset + heToDynamicIdx[nextEdgeIdx];
        currIdx = nextEdgeIdx;
        safety++;
      } while (currIdx != HE_NONE && currIdx != startOrbit && safety < 100);

      if (count > 2) {
        compiled.face_counts.push_back(count);
        for (int k = count - 1; k >= 0; --k) {
          compiled.faces.push_back(face_indices[k]);
        }
      }
    }
  }
}

template <typename MeshT>
inline void update_hankin(CompiledHankin &compiled, MeshT &out_mesh,
                          Arena &target_arena, float angle) {

  bool is_flat = std::abs(angle) < 1e-4f;

  // Precompute half-angle trig: one cosf + sinf instead of 2N
  float cos_ha = cosf(angle * 0.5f);
  float sin_ha = sinf(angle * 0.5f);

  for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
    const auto &instr = compiled.dynamicInstructions[i];
    Vector pCorner = compiled.baseVertices[instr.vCorner];

    if (is_flat) {
      compiled.dynamicVertices[i] = pCorner.normalize();
      continue;
    }

    Vector m1 = compiled.staticVertices[instr.idxM1];
    Vector m2 = compiled.staticVertices[instr.idxM2];
    Vector pPrev = compiled.baseVertices[instr.vPrev];
    Vector pNext = compiled.baseVertices[instr.vNext];

    Vector cross1 = cross(pPrev, pCorner);
    Vector cross2 = cross(pCorner, pNext);

    if (dot(cross1, cross1) < 1e-8f || dot(cross2, cross2) < 1e-8f) {
      compiled.dynamicVertices[i] = pCorner.normalize();
      continue; // zero length edge
    }

    Vector nEdge1 = cross1.normalize();
    // Inline make_rotation using precomputed half-angle trig
    Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
    Vector nHankin1 = rotate(nEdge1, q1);

    Vector nEdge2 = cross2.normalize();
    // cos(-x) = cos(x), sin(-x) = -sin(x)
    Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
    Vector nHankin2 = rotate(nEdge2, q2);

    Vector intersect = cross(nHankin1, nHankin2);
    float lenSq = dot(intersect, intersect);
    if (lenSq < 1e-6f)
      intersect = (m1 + m2).normalize();
    if (dot(intersect, pCorner) < 0)
      intersect = -intersect;

    compiled.dynamicVertices[i] = intersect.normalize();
  }

  out_mesh.vertices.initialize(target_arena,
                               compiled.staticVertices.size() +
                                   compiled.dynamicVertices.size());
  for (size_t i = 0; i < compiled.staticVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.staticVertices[i]);
  for (size_t i = 0; i < compiled.dynamicVertices.size(); ++i)
    out_mesh.vertices.push_back(compiled.dynamicVertices[i]);

  out_mesh.face_counts.initialize(target_arena, compiled.face_counts.size());

  if constexpr (requires { out_mesh.face_offsets; }) {
    out_mesh.face_offsets.initialize(target_arena, compiled.face_counts.size());
  }

  int current_offset = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i) {
    out_mesh.face_counts.push_back(compiled.face_counts[i]);
    if constexpr (requires { out_mesh.face_offsets; }) {
      out_mesh.face_offsets.push_back(static_cast<uint16_t>(current_offset));
    }
    current_offset += compiled.face_counts[i];
  }

  out_mesh.faces.initialize(target_arena, compiled.faces.size());
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    out_mesh.faces.push_back(compiled.faces[i]);
}

FLASHMEM inline PolyMesh hankin(const PolyMesh &mesh, MemoryCtx &ctx,
                                float angle) {
  ctx.swap_scratch();
  PolyMesh out;

  {
    ScopedScratch _(ctx.get_scratch_back());
    CompiledHankin compiled;
    compile_hankin(mesh, compiled, ctx.get_scratch_back(),
                   ctx.get_scratch_front());
    update_hankin(compiled, out, ctx.get_scratch_front(), angle);
  }

  return out;
}

} // namespace MeshOps
