/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file base_mesh.h
 * @brief Base mesh identities, geometry ceilings, and authoring labels.
 */
#include <cstddef>
#include <cstdint>
#include <iterator>

namespace Solids {

/** @brief Platonic, Archimedean, and Catalan base meshes. */
#define HS_BASE_MESH_LIST(X)                                                   \
  X(TETRAHEDRON, "Tetrahedron")                                                \
  X(CUBE, "Cube")                                                              \
  X(OCTAHEDRON, "Octahedron")                                                  \
  X(DODECAHEDRON, "Dodecahedron")                                              \
  X(ICOSAHEDRON, "Icosahedron")                                                \
  X(TRUNCATED_TETRAHEDRON, "Truncated Tetrahedron")                            \
  X(CUBOCTAHEDRON, "Cuboctahedron")                                            \
  X(TRUNCATED_CUBE, "Truncated Cube")                                          \
  X(TRUNCATED_OCTAHEDRON, "Truncated Octahedron")                              \
  X(RHOMBICUBOCTAHEDRON, "Rhombicuboctahedron")                                \
  X(TRUNCATED_CUBOCTAHEDRON, "Truncated Cuboctahedron")                        \
  X(SNUB_CUBE, "Snub Cube")                                                    \
  X(ICOSIDODECAHEDRON, "Icosidodecahedron")                                    \
  X(TRUNCATED_DODECAHEDRON, "Truncated Dodecahedron")                          \
  X(TRUNCATED_ICOSAHEDRON, "Truncated Icosahedron")                            \
  X(RHOMBICOSIDODECAHEDRON, "Rhombicosidodecahedron")                          \
  X(TRUNCATED_ICOSIDODECAHEDRON, "Truncated Icosidodecahedron")                \
  X(SNUB_DODECAHEDRON, "Snub Dodecahedron")                                    \
  X(TRIAKIS_TETRAHEDRON, "Triakis Tetrahedron")                                \
  X(RHOMBIC_DODECAHEDRON, "Rhombic Dodecahedron")                              \
  X(TRIAKIS_OCTAHEDRON, "Triakis Octahedron")                                  \
  X(TETRAKIS_HEXAHEDRON, "Tetrakis Hexahedron")                                \
  X(DELTOIDAL_ICOSITETRAHEDRON, "Deltoidal Icositetrahedron")                  \
  X(DISDYAKIS_DODECAHEDRON, "Disdyakis Dodecahedron")                          \
  X(PENTAGONAL_ICOSITETRAHEDRON, "Pentagonal Icositetrahedron")                \
  X(RHOMBIC_TRIACONTAHEDRON, "Rhombic Triacontahedron")                        \
  X(TRIAKIS_ICOSAHEDRON, "Triakis Icosahedron")                                \
  X(PENTAKIS_DODECAHEDRON, "Pentakis Dodecahedron")                            \
  X(DELTOIDAL_HEXECONTAHEDRON, "Deltoidal Hexecontahedron")                    \
  X(DISDYAKIS_TRIACONTAHEDRON, "Disdyakis Triacontahedron")                    \
  X(PENTAGONAL_HEXECONTAHEDRON, "Pentagonal Hexecontahedron")

enum class BaseMesh : uint8_t {
#define HS_BASE_MESH_ENUM(name, label) name,
  HS_BASE_MESH_LIST(HS_BASE_MESH_ENUM)
#undef HS_BASE_MESH_ENUM
};

inline constexpr size_t BASE_MESH_COUNT =
    static_cast<size_t>(BaseMesh::PENTAGONAL_HEXECONTAHEDRON) + 1;
inline constexpr size_t PLATONIC_BASE_MESH_COUNT = 5;

/**
 * @brief Geometry bounds every BaseMesh satisfies, for sizing arena storage.
 * @details The generators build at runtime, so these are authored ceilings the
 * loaders check each compiled solid against rather than derived counts. Tight:
 * the truncated icosidodecahedron hits 120 vertices and 180 edges, the
 * disdyakis triacontahedron 120 faces. A mesh drawn through Plot::Mesh must
 * also stay within its DEDUP_CAPACITY vertex ceiling.
 */
inline constexpr size_t MAX_SOLID_VERTICES = 120;
/** @brief Flat face-index slots; each undirected edge is walked twice. */
inline constexpr size_t MAX_SOLID_FACE_SLOTS = 360;
/** @brief Face count ceiling. */
inline constexpr size_t MAX_SOLID_FACES = 120;
/** @brief Unique edges implied by MAX_SOLID_FACE_SLOTS. */
inline constexpr size_t MAX_SOLID_EDGES = MAX_SOLID_FACE_SLOTS / 2;

inline constexpr const char *BASE_MESH_OPTIONS[] = {
#define HS_BASE_MESH_LABEL(name, label) label,
    HS_BASE_MESH_LIST(HS_BASE_MESH_LABEL)
#undef HS_BASE_MESH_LABEL
};

inline constexpr const char *BASE_MESH_EXPORT_OPTIONS[] = {
#define HS_BASE_MESH_EXPORT(name, label) "BaseMesh::" #name,
    HS_BASE_MESH_LIST(HS_BASE_MESH_EXPORT)
#undef HS_BASE_MESH_EXPORT
};
#undef HS_BASE_MESH_LIST

static_assert(BASE_MESH_COUNT == std::size(BASE_MESH_OPTIONS));
static_assert(BASE_MESH_COUNT == std::size(BASE_MESH_EXPORT_OPTIONS));
static_assert(PLATONIC_BASE_MESH_COUNT ==
              static_cast<size_t>(BaseMesh::ICOSAHEDRON) + 1);

} // namespace Solids
