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
#include <string_view>

namespace Solids {

/** @brief Platonic, Archimedean, and Catalan base meshes. */
enum class BaseMesh : uint8_t {
  TETRAHEDRON,
  CUBE,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  TRUNCATED_TETRAHEDRON,
  CUBOCTAHEDRON,
  TRUNCATED_CUBE,
  TRUNCATED_OCTAHEDRON,
  RHOMBICUBOCTAHEDRON,
  TRUNCATED_CUBOCTAHEDRON,
  SNUB_CUBE,
  ICOSIDODECAHEDRON,
  TRUNCATED_DODECAHEDRON,
  TRUNCATED_ICOSAHEDRON,
  RHOMBICOSIDODECAHEDRON,
  TRUNCATED_ICOSIDODECAHEDRON,
  SNUB_DODECAHEDRON,
  TRIAKIS_TETRAHEDRON,
  RHOMBIC_DODECAHEDRON,
  TRIAKIS_OCTAHEDRON,
  TETRAKIS_HEXAHEDRON,
  DELTOIDAL_ICOSITETRAHEDRON,
  DISDYAKIS_DODECAHEDRON,
  PENTAGONAL_ICOSITETRAHEDRON,
  RHOMBIC_TRIACONTAHEDRON,
  TRIAKIS_ICOSAHEDRON,
  PENTAKIS_DODECAHEDRON,
  DELTOIDAL_HEXECONTAHEDRON,
  DISDYAKIS_TRIACONTAHEDRON,
  PENTAGONAL_HEXECONTAHEDRON
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
    "Tetrahedron",
    "Cube",
    "Octahedron",
    "Dodecahedron",
    "Icosahedron",
    "Truncated Tetrahedron",
    "Cuboctahedron",
    "Truncated Cube",
    "Truncated Octahedron",
    "Rhombicuboctahedron",
    "Truncated Cuboctahedron",
    "Snub Cube",
    "Icosidodecahedron",
    "Truncated Dodecahedron",
    "Truncated Icosahedron",
    "Rhombicosidodecahedron",
    "Truncated Icosidodecahedron",
    "Snub Dodecahedron",
    "Triakis Tetrahedron",
    "Rhombic Dodecahedron",
    "Triakis Octahedron",
    "Tetrakis Hexahedron",
    "Deltoidal Icositetrahedron",
    "Disdyakis Dodecahedron",
    "Pentagonal Icositetrahedron",
    "Rhombic Triacontahedron",
    "Triakis Icosahedron",
    "Pentakis Dodecahedron",
    "Deltoidal Hexecontahedron",
    "Disdyakis Triacontahedron",
    "Pentagonal Hexecontahedron"};

inline constexpr const char *BASE_MESH_EXPORT_OPTIONS[] = {
    "BaseMesh::TETRAHEDRON",
    "BaseMesh::CUBE",
    "BaseMesh::OCTAHEDRON",
    "BaseMesh::DODECAHEDRON",
    "BaseMesh::ICOSAHEDRON",
    "BaseMesh::TRUNCATED_TETRAHEDRON",
    "BaseMesh::CUBOCTAHEDRON",
    "BaseMesh::TRUNCATED_CUBE",
    "BaseMesh::TRUNCATED_OCTAHEDRON",
    "BaseMesh::RHOMBICUBOCTAHEDRON",
    "BaseMesh::TRUNCATED_CUBOCTAHEDRON",
    "BaseMesh::SNUB_CUBE",
    "BaseMesh::ICOSIDODECAHEDRON",
    "BaseMesh::TRUNCATED_DODECAHEDRON",
    "BaseMesh::TRUNCATED_ICOSAHEDRON",
    "BaseMesh::RHOMBICOSIDODECAHEDRON",
    "BaseMesh::TRUNCATED_ICOSIDODECAHEDRON",
    "BaseMesh::SNUB_DODECAHEDRON",
    "BaseMesh::TRIAKIS_TETRAHEDRON",
    "BaseMesh::RHOMBIC_DODECAHEDRON",
    "BaseMesh::TRIAKIS_OCTAHEDRON",
    "BaseMesh::TETRAKIS_HEXAHEDRON",
    "BaseMesh::DELTOIDAL_ICOSITETRAHEDRON",
    "BaseMesh::DISDYAKIS_DODECAHEDRON",
    "BaseMesh::PENTAGONAL_ICOSITETRAHEDRON",
    "BaseMesh::RHOMBIC_TRIACONTAHEDRON",
    "BaseMesh::TRIAKIS_ICOSAHEDRON",
    "BaseMesh::PENTAKIS_DODECAHEDRON",
    "BaseMesh::DELTOIDAL_HEXECONTAHEDRON",
    "BaseMesh::DISDYAKIS_TRIACONTAHEDRON",
    "BaseMesh::PENTAGONAL_HEXECONTAHEDRON"};

static_assert(BASE_MESH_COUNT == std::size(BASE_MESH_OPTIONS));
static_assert(BASE_MESH_COUNT == std::size(BASE_MESH_EXPORT_OPTIONS));
static_assert(PLATONIC_BASE_MESH_COUNT ==
              static_cast<size_t>(BaseMesh::ICOSAHEDRON) + 1);

/**
 * @brief Checks both label arrays carry the expected strings for one BaseMesh.
 * @param mesh Enum value indexing both arrays.
 * @param label Expected picker label.
 * @param export_label Expected export spelling.
 * @return True when both arrays match at that index.
 */
inline constexpr bool base_mesh_labelled(BaseMesh mesh, std::string_view label,
                                         std::string_view export_label) {
  const size_t i = static_cast<size_t>(mesh);
  return std::string_view(BASE_MESH_OPTIONS[i]) == label &&
         std::string_view(BASE_MESH_EXPORT_OPTIONS[i]) == export_label;
}

// Pin every entry: both arrays are consumed positionally as BaseMesh labels,
// so any reorder or transposition must fail to compile.
static_assert(base_mesh_labelled(BaseMesh::TETRAHEDRON, "Tetrahedron",
                                 "BaseMesh::TETRAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::CUBE, "Cube", "BaseMesh::CUBE"));
static_assert(base_mesh_labelled(BaseMesh::OCTAHEDRON, "Octahedron",
                                 "BaseMesh::OCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::DODECAHEDRON, "Dodecahedron",
                                 "BaseMesh::DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::ICOSAHEDRON, "Icosahedron",
                                 "BaseMesh::ICOSAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_TETRAHEDRON,
                                 "Truncated Tetrahedron",
                                 "BaseMesh::TRUNCATED_TETRAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::CUBOCTAHEDRON, "Cuboctahedron",
                                 "BaseMesh::CUBOCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_CUBE, "Truncated Cube",
                                 "BaseMesh::TRUNCATED_CUBE"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_OCTAHEDRON,
                                 "Truncated Octahedron",
                                 "BaseMesh::TRUNCATED_OCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::RHOMBICUBOCTAHEDRON,
                                 "Rhombicuboctahedron",
                                 "BaseMesh::RHOMBICUBOCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_CUBOCTAHEDRON,
                                 "Truncated Cuboctahedron",
                                 "BaseMesh::TRUNCATED_CUBOCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::SNUB_CUBE, "Snub Cube",
                                 "BaseMesh::SNUB_CUBE"));
static_assert(base_mesh_labelled(BaseMesh::ICOSIDODECAHEDRON,
                                 "Icosidodecahedron",
                                 "BaseMesh::ICOSIDODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_DODECAHEDRON,
                                 "Truncated Dodecahedron",
                                 "BaseMesh::TRUNCATED_DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_ICOSAHEDRON,
                                 "Truncated Icosahedron",
                                 "BaseMesh::TRUNCATED_ICOSAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::RHOMBICOSIDODECAHEDRON,
                                 "Rhombicosidodecahedron",
                                 "BaseMesh::RHOMBICOSIDODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRUNCATED_ICOSIDODECAHEDRON,
                                 "Truncated Icosidodecahedron",
                                 "BaseMesh::TRUNCATED_ICOSIDODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::SNUB_DODECAHEDRON,
                                 "Snub Dodecahedron",
                                 "BaseMesh::SNUB_DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRIAKIS_TETRAHEDRON,
                                 "Triakis Tetrahedron",
                                 "BaseMesh::TRIAKIS_TETRAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::RHOMBIC_DODECAHEDRON,
                                 "Rhombic Dodecahedron",
                                 "BaseMesh::RHOMBIC_DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRIAKIS_OCTAHEDRON,
                                 "Triakis Octahedron",
                                 "BaseMesh::TRIAKIS_OCTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TETRAKIS_HEXAHEDRON,
                                 "Tetrakis Hexahedron",
                                 "BaseMesh::TETRAKIS_HEXAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::DELTOIDAL_ICOSITETRAHEDRON,
                                 "Deltoidal Icositetrahedron",
                                 "BaseMesh::DELTOIDAL_ICOSITETRAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::DISDYAKIS_DODECAHEDRON,
                                 "Disdyakis Dodecahedron",
                                 "BaseMesh::DISDYAKIS_DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::PENTAGONAL_ICOSITETRAHEDRON,
                                 "Pentagonal Icositetrahedron",
                                 "BaseMesh::PENTAGONAL_ICOSITETRAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::RHOMBIC_TRIACONTAHEDRON,
                                 "Rhombic Triacontahedron",
                                 "BaseMesh::RHOMBIC_TRIACONTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::TRIAKIS_ICOSAHEDRON,
                                 "Triakis Icosahedron",
                                 "BaseMesh::TRIAKIS_ICOSAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::PENTAKIS_DODECAHEDRON,
                                 "Pentakis Dodecahedron",
                                 "BaseMesh::PENTAKIS_DODECAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::DELTOIDAL_HEXECONTAHEDRON,
                                 "Deltoidal Hexecontahedron",
                                 "BaseMesh::DELTOIDAL_HEXECONTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::DISDYAKIS_TRIACONTAHEDRON,
                                 "Disdyakis Triacontahedron",
                                 "BaseMesh::DISDYAKIS_TRIACONTAHEDRON"));
static_assert(base_mesh_labelled(BaseMesh::PENTAGONAL_HEXECONTAHEDRON,
                                 "Pentagonal Hexecontahedron",
                                 "BaseMesh::PENTAGONAL_HEXECONTAHEDRON"));

} // namespace Solids
