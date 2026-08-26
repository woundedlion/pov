/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file solid_generators.h
 * @brief The hardcoded Platonic vertex/face tables, the SolidBuilder that
 *        chains Conway and Hankin operators over them, and the ~50 named
 *        generators the solid registry points at.
 */

#include <array>
#include <string_view>
#include "math/geometry.h"
#include "mesh/mesh.h" // For MeshOps
#include "mesh/hankin.h"
#include "mesh/conway.h"
#include "mesh/relax_bakes_generated.h"
#include <cmath>

// --- Constants for Procedural Generation ---
/** Square root of 2. */
inline constexpr float SQRT2 = 1.414213562373095f;
/** Tribonacci constant t, the real root of t^3 - t^2 - t - 1 = 0 (~1.83928676).
 */
inline constexpr float TRIBONACCI_CONST = 1.839286755214161f;
/** Snub-cube truncation parameter. */
inline constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
/** Snub-cube twist. */
inline constexpr float SNUB_CUBE_TWIST = 0.28f;
/** Truncated-dodecahedron/icosahedron truncation parameter. */
inline constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + PHI);
/** Truncated-cube/cuboctahedron truncation parameter. */
inline constexpr float T_TRUNC_CUBE = 1.0f / (2.0f + SQRT2);
/** Truncated tetra/octa/icosahedron truncation parameter. */
inline constexpr float T_TRUNC_THIRD = 1.0f / 3.0f;

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

HS_O3_BEGIN

/**
 * @brief Copies a freshly-generated mesh into the long-lived geometry arena.
 * @param temp Mesh built in the scratch arena pair.
 * @param geom Long-lived arena that backs the returned mesh.
 * @return A PolyMesh owning copies of temp's vertex/face data in geom.
 * @details Frees the scratch arenas for reuse by the next solid without
 * clobbering the result.
 */
FLASHMEM static PolyMesh finalize_solid(const PolyMesh &temp, Arena &geom) {
  PolyMesh final_mesh;
  MeshOps::clone(temp, final_mesh, geom);
  return final_mesh;
}

// ==========================================================================================
// 1. DATA DEFINITIONS (Hardcoded Platonic Solids)
// ==========================================================================================

/**
 * @brief Tetrahedron geometry data.
 */
struct Tetrahedron {
  static constexpr int NUM_VERTS = 4;
  static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f)};
  static constexpr int NUM_FACES = 4;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {3, 3, 3, 3};
  static constexpr std::array<int, 12> faces = {0, 3, 1, 0, 2, 3,
                                                0, 1, 2, 1, 3, 2};
};

/**
 * @brief Cube geometry data.
 */
struct Cube {
  static constexpr int NUM_VERTS = 8;
  static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(-0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f)};
  static constexpr int NUM_FACES = 6;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {4, 4, 4,
                                                                 4, 4, 4};
  static constexpr std::array<int, 24> faces = {
      0, 3, 2, 1, 0, 1, 5, 4, 0, 4, 7, 3, 6, 5, 1, 2, 6, 2, 3, 7, 6, 7, 4, 5};
};

/**
 * @brief Octahedron geometry data.
 */
struct Octahedron {
  static constexpr int NUM_VERTS = 6;
  static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(1.0000000000000000f, 0.0000000000000000f, 0.0000000000000000f),
      Vector(-1.0000000000000000f, 0.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, 1.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, -1.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, 0.0000000000000000f, 1.0000000000000000f),
      Vector(0.0000000000000000f, 0.0000000000000000f, -1.0000000000000000f)};
  static constexpr int NUM_FACES = 8;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {3, 3, 3, 3,
                                                                 3, 3, 3, 3};
  static constexpr std::array<int, 24> faces = {
      4, 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 5, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3};
};

/**
 * @brief Icosahedron geometry data.
 */
struct Icosahedron {
  static constexpr int NUM_VERTS = 12;
  static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(-0.5257311121191336f, 0.0000000000000000f, 0.8506508083520400f),
      Vector(0.5257311121191336f, 0.0000000000000000f, 0.8506508083520400f),
      Vector(-0.5257311121191336f, 0.0000000000000000f, -0.8506508083520400f),
      Vector(0.5257311121191336f, 0.0000000000000000f, -0.8506508083520400f),
      Vector(0.0000000000000000f, 0.8506508083520400f, 0.5257311121191336f),
      Vector(0.0000000000000000f, 0.8506508083520400f, -0.5257311121191336f),
      Vector(0.0000000000000000f, -0.8506508083520400f, 0.5257311121191336f),
      Vector(0.0000000000000000f, -0.8506508083520400f, -0.5257311121191336f),
      Vector(0.8506508083520400f, 0.5257311121191336f, 0.0000000000000000f),
      Vector(-0.8506508083520400f, 0.5257311121191336f, 0.0000000000000000f),
      Vector(0.8506508083520400f, -0.5257311121191336f, 0.0000000000000000f),
      Vector(-0.8506508083520400f, -0.5257311121191336f, 0.0000000000000000f)};
  static constexpr int NUM_FACES = 20;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  static constexpr std::array<int, 60> faces = {
      0, 1, 4, 0, 4, 9, 9,  4, 5, 4,  8, 5, 4,  1,  8,  8, 1, 10, 8,  10,
      3, 5, 8, 3, 5, 3, 2,  2, 3, 7,  7, 3, 10, 7,  10, 6, 7, 6,  11, 11,
      6, 0, 0, 6, 1, 6, 10, 1, 9, 11, 0, 9, 2,  11, 9,  5, 2, 7,  11, 2};
};

/**
 * @brief Dodecahedron geometry data.
 */
struct Dodecahedron {
  static constexpr int NUM_VERTS = 20;
  static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f),
      Vector(-0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(0.3568220897730897f, 0.9341723589627157f, 0.0000000000000000f),
      Vector(-0.3568220897730897f, 0.9341723589627157f, 0.0000000000000000f),
      Vector(0.3568220897730897f, -0.9341723589627157f, 0.0000000000000000f),
      Vector(-0.3568220897730897f, -0.9341723589627157f, 0.0000000000000000f),
      Vector(0.9341723589627157f, 0.0000000000000000f, 0.3568220897730897f),
      Vector(0.9341723589627157f, 0.0000000000000000f, -0.3568220897730897f),
      Vector(-0.9341723589627157f, 0.0000000000000000f, 0.3568220897730897f),
      Vector(-0.9341723589627157f, 0.0000000000000000f, -0.3568220897730897f),
      Vector(0.0000000000000000f, 0.3568220897730897f, 0.9341723589627157f),
      Vector(0.0000000000000000f, -0.3568220897730897f, 0.9341723589627157f),
      Vector(0.0000000000000000f, 0.3568220897730897f, -0.9341723589627157f),
      Vector(0.0000000000000000f, -0.3568220897730897f, -0.9341723589627157f)};
  static constexpr int NUM_FACES = 12;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  static constexpr std::array<int, 60> faces = {
      0, 8,  9,  4,  16, 0,  12, 13, 1,  8,  0,  16, 17, 2,  12,
      8, 1,  18, 5,  9,  12, 2,  10, 3,  13, 16, 4,  14, 6,  17,
      9, 5,  15, 14, 4,  6,  11, 10, 2,  17, 3,  19, 18, 1,  13,
      7, 15, 5,  18, 19, 7,  11, 6,  14, 15, 7,  19, 3,  10, 11};
};

/**
 * @brief Compile-time consistency check for a hardcoded solid's tables.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces arrays.
 * @return True iff the tables describe a closed orientable surface whose
 * vertices all lie on the unit sphere.
 * @details Checks that face_counts spans the flat face list exactly, every face
 * index addresses a listed vertex, no directed edge repeats and every directed
 * edge has its reverse (so each undirected edge joins exactly two faces in
 * opposite orientation, making E = sum/2), Euler's formula holds, and every
 * vertex is unit length. Face planarity, convexity and non-self-intersection
 * are not checked. Compares squared lengths so the whole check is
 * constant-evaluable.
 */
template <typename StaticMeshT> constexpr bool solid_tables_consistent() {
  size_t total_indices = 0;
  for (uint8_t count : StaticMeshT::face_counts)
    total_indices += count;
  if (total_indices != StaticMeshT::faces.size() || total_indices % 2 != 0)
    return false;

  for (int index : StaticMeshT::faces)
    if (index < 0 || static_cast<size_t>(index) >= StaticMeshT::vertices.size())
      return false;

  // Directed edge incidence; the V*V scratch keeps the pass O(indices + V^2).
  constexpr size_t V = StaticMeshT::vertices.size();
  std::array<bool, V * V> used = {};
  size_t base = 0;
  for (uint8_t count : StaticMeshT::face_counts) {
    for (uint8_t i = 0; i < count; ++i) {
      const size_t a = static_cast<size_t>(StaticMeshT::faces[base + i]);
      const size_t b = static_cast<size_t>(
          StaticMeshT::faces[base + (i + 1 == count ? 0 : i + 1)]);
      if (a == b || used[a * V + b])
        return false;
      used[a * V + b] = true;
    }
    base += count;
  }
  for (size_t a = 0; a < V; ++a)
    for (size_t b = 0; b < a; ++b)
      if (used[a * V + b] != used[b * V + a])
        return false;

  // V - E + F == 2, in a form free of unsigned wrap.
  if (StaticMeshT::vertices.size() + StaticMeshT::face_counts.size() !=
      total_indices / 2 + 2)
    return false;

  constexpr float LEN_SQ_TOL = 1e-4f;
  for (const Vector &v : StaticMeshT::vertices) {
    const float len_sq = v.x * v.x + v.y * v.y + v.z * v.z;
    if (len_sq < 1.0f - LEN_SQ_TOL || len_sq > 1.0f + LEN_SQ_TOL)
      return false;
  }
  return true;
}

static_assert(solid_tables_consistent<Tetrahedron>(),
              "Tetrahedron tables are inconsistent");
static_assert(solid_tables_consistent<Cube>(), "Cube tables are inconsistent");
static_assert(solid_tables_consistent<Octahedron>(),
              "Octahedron tables are inconsistent");
static_assert(solid_tables_consistent<Icosahedron>(),
              "Icosahedron tables are inconsistent");
static_assert(solid_tables_consistent<Dodecahedron>(),
              "Dodecahedron tables are inconsistent");

/**
 * @brief Materializes a compile-time static mesh into a runtime PolyMesh.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces
 * arrays.
 * @param target Arena that backs the returned mesh's storage.
 * @return A PolyMesh holding copies of the static mesh's data in target.
 * @details Face indices widen to the uint16_t the topology path expects.
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

#if defined(HS_RELAX_BAKE_VERIFY)
/** @brief Payloads relax_baked() has re-derived and matched this run. */
inline int relax_bakes_verified = 0;
#endif

/**
 * @brief Fluent builder for chaining Conway operators with automatic arena
 * swapping.
 * @details Each method runs `mesh = op(mesh, output_arena, scratch_arena)` then
 * swaps the two arenas; every operator returns its output in the arena passed
 * as `target` (COMPOSITION POLARITY in conway.h), so after each step the mesh
 * sits in `scratch_arena` and the next step writes into the other arena.
 *
 * The seed may sit in either arena: a base solid builds into `a`, while a
 * nested chain leaves its result in whichever arena its last step wrote. So the
 * first operator may read its input from the arena it writes its output into.
 * That costs peak arena, not correctness: every operator binds its output before
 * opening a ScratchScope, and a bump arena never rewinds below a live
 * allocation.
 *
 * Each step then rewinds the arena the NEXT step writes into back to the offset
 * it held when the chain started, reclaiming that step's spent intermediates
 * (including the ones a composed operator leaves behind in `temp`). The seed
 * sits below both marks and stays for the life of the chain.
 */
class SolidBuilder {
  PolyMesh mesh; /**< Mesh being built; updated in place by each operator. */
  Arena *output_arena;  /**< Current output arena (swapped per op). */
  Arena *scratch_arena; /**< Current scratch arena (swapped per op). */
  size_t output_mark;   /**< output_arena's offset when the chain started. */
  size_t scratch_mark;  /**< scratch_arena's offset when the chain started. */

  /**
   * @brief Swaps the arena roles and reclaims the one the next step writes.
   * @details Every operator returns its output in `output_arena`, so once the
   * roles swap the live mesh is in `scratch_arena` and everything the other
   * arena holds above its start mark is a spent intermediate.
   */
  void advance() {
    std::swap(output_arena, scratch_arena);
    std::swap(output_mark, scratch_mark);
    output_arena->set_offset(output_mark);
  }

public:
  /**
   * @brief Constructs a builder seeded with an initial mesh and arena pair.
   * @param seed Starting mesh, moved into the builder.
   * @param a Initial output arena.
   * @param b Initial scratch arena.
   */
  HS_COLD_MEMBER SolidBuilder(PolyMesh seed, Arena &a, Arena &b)
      : mesh(std::move(seed)), output_arena(&a), scratch_arena(&b),
        output_mark(a.get_offset()), scratch_mark(b.get_offset()) {
    // One arena for both roles would put the live mesh in the arena advance()
    // rewinds.
    HS_CHECK(&a != &b, "SolidBuilder: output and scratch must differ");
  }

  /**
   * @brief Applies the dual operator (faces become vertices and vice versa).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &dual() {
    mesh = MeshOps::dual(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the kis operator (raise a pyramid on every face).
   * @return Reference to this builder for chaining.
   */
  HS_COLD_MEMBER SolidBuilder &kis() {
    mesh = MeshOps::kis(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the ambo operator (rectification: new vertex per edge).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &ambo() {
    mesh = MeshOps::ambo(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the truncate operator (cut corners off each vertex).
   * @param t Truncation depth in [0, 1] along each edge (the fraction at which
   *   each cut point sits). t < 0.5 keeps the cuts on their own half; t == 0.5
   *   is rectification (short-circuits to ambo); t > 0.5 crosses the cuts past
   *   each other for intentional self-intersecting faces (the *_truncate50d_*
   *   recipes pass 50 deg ~= 0.873). See MeshOps::truncate.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &truncate(float t = 0.25f) {
    mesh = MeshOps::truncate(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the expand operator (cantellation: push faces outward).
   * @param t Expansion amount; default places square faces at the canonical
   * gap.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &expand(float t = MeshOps::EXPAND_DEFAULT_T) {
    mesh = MeshOps::expand(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the chamfer operator (replace edges with hexagons).
   * @param t Fraction each face corner moves toward the face centroid, in
   *   [0, 1); at 1 a face collapses to its centroid, so it traps. See
   *   MeshOps::chamfer.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &chamfer(float t = 0.5f) {
    mesh = MeshOps::chamfer(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the snub operator (chiral expansion with a twist).
   * @param t Inset factor of each face toward its centroid, in [0, 1); at 1 a
   *   face collapses to its centroid, so it traps.
   * @param twist Per-face rotation about the face normal, in radians; 0
   *   disables the twist pass. See MeshOps::snub.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &snub(float t = 0.5f, float twist = 0.0f) {
    mesh = MeshOps::snub(mesh, *output_arena, *scratch_arena, t, twist);
    advance();
    return *this;
  }
  /**
   * @brief Applies the gyro operator (dual of snub; pentagonal chiral
   * subdivision).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &gyro() {
    mesh = MeshOps::gyro(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Relaxes vertex positions toward a regular configuration.
   * @param iterations Upper bound on the smoothing passes; relax stops early
   *   once the springs converge, so fewer usually run. Must be non-negative;
   *   0 is a normalize-only pass-through. See MeshOps::relax.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &relax(int iterations = 8) {
    mesh = MeshOps::relax(mesh, *output_arena, *scratch_arena, iterations);
    advance();
    return *this;
  }
  /**
   * @brief Applies a host-generated relaxation payload.
   * @details The same baked vertex bits load on host and device, so both
   * render bit-identically. The payload's own `iterations` count defines it;
   * there is no per-call-site count. In the two host tooling modes the payload
   * is instead reproduced live by running exactly `bake.iterations` smoothing
   * steps: EXTRACT dumps the resulting bits and freshly measured guards, so it
   * can author a payload or recover from a topology change (generation), while
   * VERIFY asserts both against the committed payload (the native
   * re-derivation test).
   * @param bake Payload whose guarded topology must match the current mesh.
   * @return Reference to this builder for chaining.
   */
  HS_COLD_MEMBER SolidBuilder &relax_baked(const MeshOps::RelaxBake &bake) {
#if defined(HS_RELAX_BAKE_EXTRACT) || defined(HS_RELAX_BAKE_VERIFY)
#if defined(HS_RELAX_BAKE_VERIFY)
    HS_CHECK(MeshOps::relax_topology_hash(mesh) == bake.topology_hash,
             "relax bake verify: source topology differs");
    MeshOps::check_relax_bake_source(mesh, bake);
#else
    const uint32_t topology_hash = MeshOps::relax_topology_hash(mesh);
    const uint32_t source_hash = MeshOps::relax_source_hash(mesh);
    const float source_margin = MeshOps::relax_source_quantization_margin(mesh);
#endif
    mesh = MeshOps::relax(mesh, *output_arena, *scratch_arena, bake.iterations);
    uint32_t output_hash = MeshOps::FNV1A_BASIS;
    for (const Vector &v : mesh.vertices) {
      output_hash =
          MeshOps::fnv1a_step(output_hash, std::bit_cast<uint32_t>(v.x));
      output_hash =
          MeshOps::fnv1a_step(output_hash, std::bit_cast<uint32_t>(v.y));
      output_hash =
          MeshOps::fnv1a_step(output_hash, std::bit_cast<uint32_t>(v.z));
    }
#if defined(HS_RELAX_BAKE_VERIFY)
    HS_CHECK(mesh.vertices.size() == bake.vertex_count &&
                 mesh.face_counts.size() == bake.face_count &&
                 mesh.faces.size() == bake.index_count,
             "relax bake verify: dimensions differ");
    HS_CHECK(output_hash == bake.output_hash,
             "relax bake verify: output differs");
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
      HS_CHECK(std::bit_cast<uint32_t>(mesh.vertices[i].x) ==
                       bake.vertex_bits[3 * i] &&
                   std::bit_cast<uint32_t>(mesh.vertices[i].y) ==
                       bake.vertex_bits[3 * i + 1] &&
                   std::bit_cast<uint32_t>(mesh.vertices[i].z) ==
                       bake.vertex_bits[3 * i + 2],
               "relax bake verify: vertex differs");
    }
    ++relax_bakes_verified;
#else // HS_RELAX_BAKE_EXTRACT: emit the payload for the generated header.
    hs::log("RELAX_BAKE_BEGIN %s %d %lu %lu %lu %08lx %08lx %08lx %d %08lx "
            "%08lx %08lx",
            bake.name, static_cast<int>(bake.iterations),
            static_cast<unsigned long>(mesh.vertices.size()),
            static_cast<unsigned long>(mesh.face_counts.size()),
            static_cast<unsigned long>(mesh.faces.size()),
            static_cast<unsigned long>(topology_hash),
            static_cast<unsigned long>(source_hash),
            static_cast<unsigned long>(output_hash),
            static_cast<int>(MeshOps::RELAX_SOURCE_SCALE),
            static_cast<unsigned long>(
                std::bit_cast<uint32_t>(MeshOps::RELAX_SOURCE_BIAS)),
            static_cast<unsigned long>(
                std::bit_cast<uint32_t>(MeshOps::RELAX_SOURCE_MIN_MARGIN)),
            static_cast<unsigned long>(std::bit_cast<uint32_t>(source_margin)));
    for (const Vector &v : mesh.vertices)
      hs::log("RELAX_BAKE_DATA %08lx %08lx %08lx",
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.x)),
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.y)),
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.z)));
    hs::log("RELAX_BAKE_END");
#endif
    advance();
    return *this;
#else
    mesh = MeshOps::relax_baked(mesh, *output_arena, bake);
    advance();
    return *this;
#endif
  }
  /**
   * @brief Applies the meta operator (kis composed with join).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &meta() {
    mesh = MeshOps::meta(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the needle operator (kis of the dual; n = kd).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &needle() {
    mesh = MeshOps::needle(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the zip operator (dual of kis).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &zip() {
    mesh = MeshOps::zip(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the bevel operator (truncate composed with ambo).
   * @param t Truncation depth forwarded to the truncate step, in [0, 1]. At
   *   exactly 0.5 truncate short-circuits to ambo, so the chain is ambo(ambo)
   *   with that census rather than a bevel. See MeshOps::bevel.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &bevel(float t = 0.25f) {
    mesh = MeshOps::bevel(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the Hankin star-pattern operator to each face.
   * @param angle Contact angle of the star pattern, in radians.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &hankin(float angle) {
    // *scratch_arena may hold mesh itself; hankin's ScratchScope marks above
    // the input, so compiling the topology into it leaves the input intact.
    mesh = MeshOps::hankin(mesh, *output_arena, *scratch_arena, angle);
    advance();
    return *this;
  }

  /**
   * @brief Finalizes the chain and yields the built mesh.
   * @return The accumulated PolyMesh, moved out of the builder.
   */
  HS_COLD_MEMBER PolyMesh build() { return std::move(mesh); }
};

// ==========================================================================================
// 2. PROCEDURAL GENERATORS
// ==========================================================================================

namespace Platonic {
/**
 * @brief Builds a tetrahedron (V=4, F=4, I=12).
 * @param a Arena that backs the returned mesh.
 * @return The tetrahedron mesh.
 */
FLASHMEM static PolyMesh tetrahedron(Arena &a, Arena &) {
  return to_polymesh<Tetrahedron>(a);
}
/**
 * @brief Builds a cube (V=8, F=6, I=24).
 * @param a Arena that backs the returned mesh.
 * @return The cube mesh.
 */
FLASHMEM static PolyMesh cube(Arena &a, Arena &) {
  return to_polymesh<Cube>(a);
}
/**
 * @brief Builds an octahedron (V=6, F=8, I=24).
 * @param a Arena that backs the returned mesh.
 * @return The octahedron mesh.
 */
FLASHMEM static PolyMesh octahedron(Arena &a, Arena &) {
  return to_polymesh<Octahedron>(a);
}
/**
 * @brief Builds a dodecahedron (V=20, F=12, I=60).
 * @param a Arena that backs the returned mesh.
 * @return The dodecahedron mesh.
 */
FLASHMEM static PolyMesh dodecahedron(Arena &a, Arena &) {
  return to_polymesh<Dodecahedron>(a);
}
/**
 * @brief Builds an icosahedron (V=12, F=20, I=60).
 * @param a Arena that backs the returned mesh.
 * @return The icosahedron mesh.
 */
FLASHMEM static PolyMesh icosahedron(Arena &a, Arena &) {
  return to_polymesh<Icosahedron>(a);
}
} // namespace Platonic

namespace Archimedean {
using namespace Platonic;

/**
 * @brief Builds a truncated tetrahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated tetrahedron mesh.
 */
FLASHMEM static PolyMesh truncatedTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(tetrahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a cuboctahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The cuboctahedron mesh.
 */
FLASHMEM static PolyMesh cuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).ambo().build();
}
/**
 * @brief Builds a truncated cube.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated cube mesh.
 */
FLASHMEM static PolyMesh truncatedCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).truncate(T_TRUNC_CUBE).build();
}
/**
 * @brief Builds a truncated octahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated octahedron mesh.
 */
FLASHMEM static PolyMesh truncatedOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a rhombicuboctahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rhombicuboctahedron mesh.
 */
FLASHMEM static PolyMesh rhombicuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).expand().build();
}
/**
 * @brief Builds a truncated cuboctahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated cuboctahedron mesh.
 */
FLASHMEM static PolyMesh truncatedCuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b)
      .bevel(T_TRUNC_CUBE)
      .relax_baked(RelaxBakes::truncated_cuboctahedron_converged)
      .build();
}
/**
 * @brief Builds a snub cube.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The snub cube mesh.
 */
FLASHMEM static PolyMesh snubCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b)
      .snub(T_SNUB_CUBE, SNUB_CUBE_TWIST)
      .relax_baked(RelaxBakes::snub_cube_converged)
      .build();
}
/**
 * @brief Builds an icosidodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The icosidodecahedron mesh.
 */
FLASHMEM static PolyMesh icosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).ambo().build();
}
/**
 * @brief Builds a truncated dodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated dodecahedron mesh.
 */
FLASHMEM static PolyMesh truncatedDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).truncate(T_TRUNC_ICOS).build();
}
/**
 * @brief Builds a truncated icosahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated icosahedron mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a rhombicosidodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rhombicosidodecahedron mesh.
 */
FLASHMEM static PolyMesh rhombicosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .expand()
      .relax_baked(RelaxBakes::rhombicosidodecahedron_converged)
      .build();
}
/**
 * @brief Builds a truncated icosidodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated icosidodecahedron mesh.
 */
FLASHMEM static PolyMesh truncatedIcosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .bevel(T_TRUNC_ICOS)
      .relax_baked(RelaxBakes::truncated_icosidodecahedron_converged)
      .build();
}
/**
 * @brief Builds a snub dodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The snub dodecahedron mesh.
 */
FLASHMEM static PolyMesh snubDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .snub(0.5f)
      .relax_baked(RelaxBakes::snub_dodecahedron_converged)
      .build();
}
} // namespace Archimedean

namespace Catalan {
using namespace Archimedean;

/**
 * @brief Builds a triakis tetrahedron (dual of the truncated tetrahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The triakis tetrahedron mesh.
 */
FLASHMEM static PolyMesh triakisTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedTetrahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a rhombic dodecahedron (dual of the cuboctahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rhombic dodecahedron mesh.
 */
FLASHMEM static PolyMesh rhombicDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(cuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a triakis octahedron (dual of the truncated cube).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The triakis octahedron mesh.
 */
FLASHMEM static PolyMesh triakisOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCube(a, b), a, b).dual().build();
}
/**
 * @brief Builds a tetrakis hexahedron (dual of the truncated octahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The tetrakis hexahedron mesh.
 */
FLASHMEM static PolyMesh tetrakisHexahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedOctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a deltoidal icositetrahedron (dual of the rhombicuboctahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The deltoidal icositetrahedron mesh.
 */
FLASHMEM static PolyMesh deltoidalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a disdyakis dodecahedron (dual of the truncated cuboctahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The disdyakis dodecahedron mesh.
 */
FLASHMEM static PolyMesh disdyakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentagonal icositetrahedron (dual of the snub cube).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The pentagonal icositetrahedron mesh.
 */
FLASHMEM static PolyMesh pentagonalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubCube(a, b), a, b).dual().build();
}
/**
 * @brief Builds a rhombic triacontahedron (dual of the icosidodecahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rhombic triacontahedron mesh.
 */
FLASHMEM static PolyMesh rhombicTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a triakis icosahedron (dual of the truncated dodecahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The triakis icosahedron mesh.
 */
FLASHMEM static PolyMesh triakisIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedDodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentakis dodecahedron (dual of the truncated icosahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The pentakis dodecahedron mesh.
 */
FLASHMEM static PolyMesh pentakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a deltoidal hexecontahedron (dual of the
 * rhombicosidodecahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The deltoidal hexecontahedron mesh.
 */
FLASHMEM static PolyMesh deltoidalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a disdyakis triacontahedron (dual of the truncated
 * icosidodecahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The disdyakis triacontahedron mesh.
 */
FLASHMEM static PolyMesh disdyakisTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentagonal hexecontahedron (dual of the snub dodecahedron).
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The pentagonal hexecontahedron mesh.
 */
FLASHMEM static PolyMesh pentagonalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubDodecahedron(a, b), a, b).dual().build();
}
} // namespace Catalan

namespace IslamicStarPatterns {

/** Degrees-to-radians conversion factor. */
inline constexpr float D2R = PI_F / 180.0f;

/** Truncation depth of the `*_truncate5d_*` recipes, bit-exactly 5.0f * D2R
 * and named for it, consumed by truncate as a dimensionless edge fraction
 * short of the ambo pinch at t = 0.5. */
inline constexpr float TRUNCATE_T_NEAR = 0.0872664601f;
/** Truncation depth of the `*_truncate50d_*` recipes, bit-exactly 50.0f * D2R
 * and named for it, consumed by truncate as a dimensionless edge fraction past
 * the ambo pinch, where the cut faces self-intersect by design. */
inline constexpr float TRUNCATE_T_FAR = 0.87266463f;

/**
 * @brief Builds the truncatedIcosahedron_hk58_chamfer63 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_hk58_chamfer63(Arena &a,
                                                             Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .hankin(58.0f * D2R)
      .chamfer(0.63f)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk62_ambo_hk62 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk62_ambo_hk62(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(62.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .build();
}
/**
 * @brief Builds the octahedron_hk17_ambo_hk73 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh octahedron_hk17_ambo_hk73(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::octahedron(a, b), a, b)
      .hankin(17.0f * D2R)
      .ambo()
      .hankin(73.0f * D2R)
      .build();
}
/**
 * @brief Builds the icosahedron_kis_gyro star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_kis_gyro(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b).kis().gyro().build();
}
/**
 * @brief Builds the truncatedIcosidodecahedron_truncate50d_ambo_dual star
 * pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosidodecahedron_truncate50d_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosidodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_FAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the icosidodecahedron_truncate5d_ambo_dual star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosidodecahedron_truncate5d_ambo_dual(Arena &a,
                                                                Arena &b) {
  return SolidBuilder(Archimedean::icosidodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the snubDodecahedron_truncate5d_ambo_dual star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh snubDodecahedron_truncate5d_ambo_dual(Arena &a,
                                                               Arena &b) {
  return SolidBuilder(Archimedean::snubDodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the octahedron_hk34_ambo_hk72 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh octahedron_hk34_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::octahedron(a, b), a, b)
      .hankin(34.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the rhombicuboctahedron_hk63_ambo_hk63 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh rhombicuboctahedron_hk63_ambo_hk63(Arena &a,
                                                            Arena &b) {
  return SolidBuilder(Archimedean::rhombicuboctahedron(a, b), a, b)
      .hankin(63.0f * D2R)
      .ambo()
      .hankin(63.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_hk54_ambo_hk72 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_hk54_ambo_hk72(Arena &a,
                                                             Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk54_ambo_hk72 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk54_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk72_ambo_dual_hk20 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk72_ambo_dual_hk20(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(72.0f * D2R)
      .ambo()
      .dual()
      .hankin(20.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_truncate50d_ambo_dual star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_truncate50d_ambo_dual(Arena &a,
                                                                    Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .truncate(TRUNCATE_T_FAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the icosahedron_snub_relax_truncate033_hankin62 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_snub_relax_truncate033_hankin62(Arena &a,
                                                                     Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .snub()
      .relax_baked(RelaxBakes::icosahedron_snub_converged)
      .truncate(0.33f)
      .hankin(62.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 * @details The final contact angle sits clear of the ~43-degree resonance
 * where one corner class's contact planes go near-parallel.
 */
FLASHMEM static PolyMesh dodecahedron_hk35_ambo_hk62_ambo_relax_hk42(Arena &a,
                                                                     Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .ambo()
      .relax_baked(RelaxBakes::dodecahedron_hankin_ambo_hankin_ambo_converged)
      .hankin(42.0f * D2R)
      .build();
}
/**
 * @brief Builds the icosahedron_ambo_truncate033_hankin59 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_ambo_truncate033_hankin59(Arena &a,
                                                               Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .ambo()
      .truncate(0.33f)
      .hankin(59.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin59 star
 * pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin59(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.01f)
      .hankin(59.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin73 star
 * pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin73(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.01f)
      .hankin(73.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedOctahedron_gyro_kis_hk17 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedOctahedron_gyro_kis_hk17(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedOctahedron(a, b), a, b)
      .gyro()
      .kis()
      .hankin(17.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosidodecahedron_bevel5_relax_hk77 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosidodecahedron_bevel5_relax_hk77(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosidodecahedron(a, b), a, b)
      .bevel(0.5f)
      .relax_baked(RelaxBakes::truncated_icosidodecahedron_bevel50_relax100)
      .hankin(77.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_bevel2_relax_gyro star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_bevel2_relax_gyro(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .bevel(0.2f)
      .relax_baked(RelaxBakes::dodecahedron_bevel20_converged)
      .gyro()
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate33_hk64 star
 * pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate33_hk64(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.33f)
      .hankin(64.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_ambo_bevel33_relax_hk66 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_ambo_bevel33_relax_hk66(Arena &a,
                                                              Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .ambo()
      .bevel(0.33f)
      .relax_baked(RelaxBakes::dodecahedron_ambo_bevel33_converged)
      .hankin(66.0f * D2R)
      .build();
}
} // namespace IslamicStarPatterns

HS_O3_END

} // namespace Solids
