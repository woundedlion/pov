/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file mesh_marshal.h
 * @brief Pure (no-Emscripten) validate-and-build layer for the WASM mesh-editor
 *        boundary.
 *
 * The JS mesh editor hands the engine a wholesale mesh as two flat arrays — a
 * Float32 vertex stream [x, y, z, ...] and an Int32 face stream with -1 face
 * delimiters. This is the most input-fragile surface in the bridge: every entry
 * is untrusted, and a malformed mesh that slipped through would later trap a
 * fail-fast HS_CHECK and abort the whole WASM module. The contract is therefore
 * reject-don't-trap — each malformed case returns a distinct status and the
 * caller maps non-Ok to a null result (mirroring setClip/fromSolid).
 *
 * The validation + raw-PolyMesh build is factored out of wasm.cpp (which only
 * adds the emscripten::val array conversion on top) so every reject branch is
 * host-unit-testable without an Emscripten toolchain — see
 * tests/test_mesh_marshal.h.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/memory.h" // Arena, ArenaVector
#include "core/mesh.h"   // PolyMesh, Vector
#include "platform.h"    // hs::log

namespace hs_wasm {

/**
 * @brief Outcome of build_mesh_from_flat(). Every non-kOk value names exactly
 *        one reject branch; the caller maps anything but kOk to a null mesh.
 */
enum class MeshBuildStatus {
  kOk,
  kVerticesNotMultipleOf3, ///< vertex stream length is not divisible by 3
  kFaceIndexOutOfRange,    ///< a face entry is < -1 or >= vertex count
  kFaceTooManyVerts,       ///< a single face has more than UINT8_MAX vertices
  kIndexRangeOverflow,     ///< vertex or half-edge count exceeds UINT16_MAX
  kArenaOverflow,          ///< the well-formed mesh exceeds the build arena
};

/**
 * @brief Validate the flat (vertex, face) streams and, if well-formed, build the
 *        raw PolyMesh into @p arena.
 *
 * @param vData  Flat vertex stream [x, y, z, ...]; length must be a multiple of 3.
 * @param fData  Flat face stream of vertex indices in [0, num_verts), with -1
 *               face delimiters. A trailing face needs no final -1; empty faces
 *               (leading/consecutive -1) emit nothing.
 * @param arena  Build arena; reset and bound into ONLY on success, so a rejected
 *               mesh leaves it untouched.
 * @param out    Receives the built vertices/faces/face_counts on kOk; untouched
 *               otherwise.
 *
 * On any non-kOk return the specific reason is logged (hs::log) here, so the
 * caller need only translate the status to its null/empty result.
 */
inline MeshBuildStatus build_mesh_from_flat(const std::vector<float> &vData,
                                            const std::vector<int> &fData,
                                            Arena &arena, PolyMesh &out) {
  // A vertex stream whose length is not a multiple of 3 would over-read
  // vData[i+1]/[i+2] past the end during the build below.
  if (vData.size() % 3 != 0) {
    hs::log("WASM: fromData vertices length %zu not a multiple of 3 — ignored",
            vData.size());
    return MeshBuildStatus::kVerticesNotMultipleOf3;
  }
  const int num_verts = static_cast<int>(vData.size() / 3);

  // Every entry must be the -1 delimiter or a vertex index in [0, num_verts).
  // An out-of-range index would later dereference a nonexistent vertex; reject
  // the whole mesh up front.
  for (int idx : fData) {
    if (idx < -1 || idx >= num_verts) {
      hs::log("WASM: fromData face index %d out of range [0,%d) — ignored", idx,
              num_verts);
      return MeshBuildStatus::kFaceIndexOutOfRange;
    }
  }

  // Size the arenas by simulating the build pass rather than tallying
  // delimiters: an empty face (a leading or consecutive -1) emits no face_count,
  // while a trailing face with no final -1 still counts. A naive "one face per
  // -1" count would over-allocate face_counts on empty faces, and the "size
  // minus delimiters" trick would under-size faces when the stream doesn't end
  // on a -1. Counting exactly what the loop pushes keeps both arenas right-sized.
  int num_faces = 0;      // non-empty faces => face_counts entries
  int num_face_verts = 0; // total vertex indices => faces entries
  {
    int run = 0;
    for (int idx : fData) {
      if (idx == -1) {
        if (run > 0) {
          num_faces++;
          run = 0;
        }
      } else {
        num_face_verts++;
        // A face_count is stored in a uint8_t, so a face with more than
        // UINT8_MAX vertices would wrap on the (uint8_t) casts below (256 -> 0
        // verts, corrupting topology). Reject it up front rather than routing
        // through a narrow whose HS_CHECK would trap and abort the module.
        if (++run > UINT8_MAX) {
          hs::log("WASM: fromData face has > %d vertices — ignored", UINT8_MAX);
          return MeshBuildStatus::kFaceTooManyVerts;
        }
      }
    }
    if (run > 0)
      num_faces++;
  }

  // The half-edge builder stores every vertex index and half-edge index in a
  // uint16_t and traps (build_from_flat's HS_CHECK) past that range. num_verts
  // bounds the index values pushed into out.faces and num_face_verts is the
  // half-edge count; either exceeding UINT16_MAX would silently truncate at the
  // (uint16_t) out.faces store below and then abort the whole module on the next
  // operator. Reject up front like the other malformed-input cases.
  if (num_verts > UINT16_MAX || num_face_verts > UINT16_MAX) {
    hs::log("WASM: fromData mesh exceeds 16-bit index range "
            "(%d verts, %d face indices) — ignored",
            num_verts, num_face_verts);
    return MeshBuildStatus::kIndexRangeOverflow;
  }

  // A well-formed mesh can still be too large for the build arena. Budget the
  // three allocations (plus one alignment gap each) and soft-reject an oversize
  // mesh — the index-range check above only bounds values, not total size, and
  // an over-allocation would otherwise trip the arena trap and abort the module.
  const size_t needed = (size_t)num_verts * sizeof(Vector) +
                        (size_t)num_face_verts * sizeof(uint16_t) +
                        (size_t)num_faces * sizeof(uint8_t) +
                        3 * alignof(std::max_align_t);
  if (needed > arena.get_capacity()) {
    hs::log("WASM: fromData mesh needs %zu B > tooling arena %zu B — ignored",
            needed, arena.get_capacity());
    return MeshBuildStatus::kArenaOverflow;
  }

  // Validation passed. fromData rebuilds a mesh wholesale and reads no prior
  // arena state, so reset first; this reclaims a long editor session's creep and
  // lets the new mesh use the arena's full capacity.
  arena.reset();

  out.vertices.bind(arena, num_verts);
  for (size_t i = 0; i < vData.size(); i += 3) {
    out.vertices.emplace_back(vData[i], vData[i + 1], vData[i + 2]);
  }

  out.faces.bind(arena, num_face_verts);
  out.face_counts.bind(arena, num_faces);

  // current_count is bounded to UINT8_MAX by the counting pass above, so the
  // (uint8_t) casts below cannot wrap.
  int current_count = 0;
  for (int idx : fData) {
    if (idx == -1) {
      if (current_count > 0) {
        out.face_counts.push_back((uint8_t)current_count);
        current_count = 0;
      }
    } else {
      out.faces.push_back(idx);
      current_count++;
    }
  }
  if (current_count > 0) {
    out.face_counts.push_back((uint8_t)current_count);
  }

  return MeshBuildStatus::kOk;
}

} // namespace hs_wasm
