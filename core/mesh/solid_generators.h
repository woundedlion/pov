/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file solid_generators.h
 * @brief Umbrella for solid tables, SolidBuilder and procedural generators,
 *        with shared truncation/snub constants and finalize_solid.
 */

#include <array>
#include "core/mesh/base_mesh.h"
#include "core/math/geometry.h"
#include "core/mesh/mesh.h" // For MeshOps
#include "core/mesh/hankin.h"
#include "core/mesh/conway.h"
#include "core/mesh/relax_bakes_generated.h"
#include "core/mesh/solid_tables.h"
#include <cmath>

namespace Solids {

// --- Constants for Procedural Generation ---
/** Square root of 2. */
inline constexpr float SQRT2 = 1.414213562373095f;
/** Tribonacci constant t, the real root of t^3 - t^2 - t - 1 = 0 (~1.83928676).
 */
inline constexpr float TRIBONACCI_CONST = 1.839286755214161f;
/** Snub-cube inset parameter. */
inline constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
/** Snub-cube twist. */
inline constexpr float SNUB_CUBE_TWIST = 0.28f;
/** Truncated-dodecahedron / truncated-icosidodecahedron (bevel) truncation parameter. */
inline constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + math::PHI);
/** Truncated-cube/cuboctahedron truncation parameter. */
inline constexpr float T_TRUNC_CUBE = 1.0f / (2.0f + SQRT2);
/** Truncated tetra/octa/icosahedron truncation parameter. */
inline constexpr float T_TRUNC_THIRD = 1.0f / 3.0f;

HS_O3_BEGIN

/**
 * @brief Copies a freshly-generated mesh into the long-lived geometry arena.
 * @param temp Mesh built in the scratch arena pair.
 * @param geom Long-lived arena that backs the returned mesh.
 * @return A PolyMesh owning copies of temp's vertex/face data in geom.
 * @details Copies into geom so the caller can rewind the scratch pair without
 * clobbering the result.
 */
FLASHMEM static PolyMesh finalize_solid(const PolyMesh &temp, Arena &geom) {
  PolyMesh final_mesh;
  MeshOps::clone(temp, final_mesh, geom);
  return final_mesh;
}

/**
 * @brief Materializes a compile-time static mesh into a runtime PolyMesh.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces
 * arrays.
 * @param target Arena that backs the returned mesh's storage.
 * @return A PolyMesh holding copies of the static mesh's data in target.
 * @details Face indices narrow to uint16_t, trapping past MeshLimits::MAX_VERTEX_INDEX.
 */
template <typename StaticMeshT> PolyMesh to_polymesh(Arena &target) {
  PolyMesh mesh;
  mesh.vertices.bind(target, StaticMeshT::vertices.size());
  mesh.vertices.append_bulk(StaticMeshT::vertices.data(),
                            StaticMeshT::vertices.size());
  mesh.face_counts.bind(target, StaticMeshT::face_counts.size());
  mesh.face_counts.append_bulk(StaticMeshT::face_counts.data(),
                               StaticMeshT::face_counts.size());
  mesh.faces.bind(target, StaticMeshT::faces.size());
  for (const auto &f : StaticMeshT::faces)
    mesh.faces.push_back(MeshOps::narrow_index(static_cast<size_t>(f)));
  return mesh;
}

#include "core/mesh/solid_builder.h"
#include "core/mesh/procedural_solids.h"
HS_O3_END

} // namespace Solids
