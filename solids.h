/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "geometry.h"
#include "mesh.h" // For MeshOps
#include <cmath>
#include <string>
#include <vector>
#include <map>

// --- Constants for Procedural Generation ---
static constexpr float SQRT2 = 1.414213562373095f;
static constexpr float TRIBONACCI =
    (1.0f + 1.839286755214161f + 1.839286755214161f) / 3.0f; // Approx
// Real tribonacci constant t: t^3 - t^2 - t - 1 = 0. ~1.839286755214161
static constexpr float TRIBONACCI_CONST = 1.839286755214161f;
static constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
static constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + PHI);

namespace Solids {

inline PolyMesh finalize_solid(const PolyMesh &temp, Arena &geom) {
  PolyMesh final_mesh;
  final_mesh.vertices.initialize(geom, temp.vertices.size());
  for (size_t i = 0; i < temp.vertices.size(); ++i)
    final_mesh.vertices.push_back(temp.vertices[i]);
  final_mesh.face_counts.initialize(geom, temp.face_counts.size());
  for (size_t i = 0; i < temp.face_counts.size(); ++i)
    final_mesh.face_counts.push_back(temp.face_counts[i]);
  final_mesh.faces.initialize(geom, temp.faces.size());
  for (size_t i = 0; i < temp.faces.size(); ++i)
    final_mesh.faces.push_back(temp.faces[i]);
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

// Helper to convert Static Mesh Data to Dynamic PolyMesh
template <typename StaticMeshT> PolyMesh to_polymesh(Arena &scratch) {
  PolyMesh mesh;
  mesh.vertices.initialize(scratch, StaticMeshT::vertices.size());
  for (const auto &v : StaticMeshT::vertices)
    mesh.vertices.push_back(v);
  mesh.face_counts.initialize(scratch, StaticMeshT::face_counts.size());
  for (const auto &c : StaticMeshT::face_counts)
    mesh.face_counts.push_back(c);
  mesh.faces.initialize(scratch, StaticMeshT::faces.size());
  for (const auto &f : StaticMeshT::faces)
    mesh.faces.push_back(f);
  return mesh;
}

// ==========================================================================================
// 2. PROCEDURAL GENERATORS
// ==========================================================================================

namespace Platonic {
inline PolyMesh tetrahedron(Arena &scratch) {
  return to_polymesh<Tetrahedron>(scratch);
}
inline PolyMesh cube(Arena &scratch) { return to_polymesh<Cube>(scratch); }
inline PolyMesh octahedron(Arena &scratch) {
  return to_polymesh<Octahedron>(scratch);
}
inline PolyMesh dodecahedron(Arena &scratch) {
  return to_polymesh<Dodecahedron>(scratch);
}
inline PolyMesh icosahedron(Arena &scratch) {
  return to_polymesh<Icosahedron>(scratch);
}
} // namespace Platonic

namespace Archimedean {
using namespace Platonic;
using namespace MeshOps;

inline PolyMesh truncatedTetrahedron(Arena &scratch) {
  return truncate(tetrahedron(scratch), scratch, 1.0f / 3.0f);
}
inline PolyMesh cuboctahedron(Arena &scratch) {
  return ambo(cube(scratch), scratch);
}
inline PolyMesh truncatedCube(Arena &scratch) {
  return truncate(cube(scratch), scratch, 1.0f / (2.0f + SQRT2));
}
inline PolyMesh truncatedOctahedron(Arena &scratch) {
  return truncate(octahedron(scratch), scratch, 1.0f / 3.0f);
}
inline PolyMesh rhombicuboctahedron(Arena &scratch) {
  return expand(cube(scratch), scratch);
}
inline PolyMesh truncatedCuboctahedron(Arena &scratch) {
  return canonicalize(bitruncate(cube(scratch), scratch, 1.0f / (2.0f + SQRT2)),
                      scratch, 50);
}
inline PolyMesh snubCube(Arena &scratch) {
  return canonicalize(snub(cube(scratch), scratch, T_SNUB_CUBE, 0.28f), scratch,
                      50);
}
inline PolyMesh icosidodecahedron(Arena &scratch) {
  return ambo(dodecahedron(scratch), scratch);
}
inline PolyMesh truncatedDodecahedron(Arena &scratch) {
  return truncate(dodecahedron(scratch), scratch, T_TRUNC_ICOS);
}
inline PolyMesh truncatedIcosahedron(Arena &scratch) {
  return truncate(icosahedron(scratch), scratch, 1.0f / 3.0f);
}
inline PolyMesh rhombicosidodecahedron(Arena &scratch) {
  return canonicalize(expand(dodecahedron(scratch), scratch), scratch, 50);
}
inline PolyMesh truncatedIcosidodecahedron(Arena &scratch) {
  return canonicalize(
      bitruncate(dodecahedron(scratch), scratch, 1.0f / (2.0f + PHI)), scratch,
      50);
}
inline PolyMesh snubDodecahedron(Arena &scratch) {
  return canonicalize(snub(dodecahedron(scratch), scratch, 0.5f), scratch, 50);
}
} // namespace Archimedean

namespace IslamicStarPatterns {
using namespace Platonic;
using namespace Archimedean;
using namespace MeshOps;

static constexpr float D2R = PI_F / 180.0f;

inline PolyMesh bitruncate(const PolyMesh &m, Arena &scratch, float t) {
  return truncate(ambo(m, scratch), scratch, t);
}

inline PolyMesh icosahedron_hk59_bitruncate033(Arena &scratch) {
  return hankin(bitruncate(icosahedron(scratch), scratch, 0.33f), scratch,
                59.0f * D2R);
}

inline PolyMesh octahedron_hk17_ambo_hk72(Arena &scratch) {
  return hankin(
      ambo(hankin(octahedron(scratch), scratch, 17.0f * D2R), scratch), scratch,
      73.0f * D2R);
}
inline PolyMesh icosahedron_kis_gyro(Arena &scratch) {
  return gyro(kis(icosahedron(scratch), scratch), scratch);
}

inline PolyMesh
truncatedIcosidodecahedron_truncate05_ambo_dual(Arena &scratch) {
  return dual(
      ambo(truncate(truncatedIcosidodecahedron(scratch), scratch, 50.0f * D2R),
           scratch),
      scratch);
}

inline PolyMesh icosidodecahedron_truncate05_ambo_dual(Arena &scratch) {
  return dual(
      ambo(truncate(icosidodecahedron(scratch), scratch, 5.0f * D2R), scratch),
      scratch);
}
inline PolyMesh snubDodecahedron_truncate05_ambo_dual(Arena &scratch) {
  return dual(
      ambo(truncate(snubDodecahedron(scratch), scratch, 5.0f * D2R), scratch),
      scratch);
}

inline PolyMesh octahedron_hk34_ambo_hk72(Arena &scratch) {
  return hankin(
      ambo(hankin(octahedron(scratch), scratch, 34.0f * D2R), scratch), scratch,
      72.0f * D2R);
}

inline PolyMesh rhombicuboctahedron_hk63_ambo_hk63(Arena &scratch) {
  return hankin(
      ambo(hankin(rhombicuboctahedron(scratch), scratch, 63.0f * D2R), scratch),
      scratch, 63.0f * D2R);
}

inline PolyMesh truncatedIcosahedron_hk54_ambo_hk72(Arena &scratch) {
  return hankin(
      ambo(hankin(truncatedIcosahedron(scratch), scratch, 54.0f * D2R),
           scratch),
      scratch, 72.0f * D2R);
}

inline PolyMesh dodecahedron_hk54_ambo_hk72(Arena &scratch) {
  return hankin(
      ambo(hankin(dodecahedron(scratch), scratch, 54.0f * D2R), scratch),
      scratch, 72.0f * D2R);
}

inline PolyMesh dodecahedron_hk72_ambo_dual_hk20(Arena &scratch) {
  return hankin(
      dual(ambo(hankin(dodecahedron(scratch), scratch, 72.0f * D2R), scratch),
           scratch),
      scratch, 20.0f * D2R);
}

inline PolyMesh truncatedIcosahedron_truncate05_ambo_dual(Arena &scratch) {
  return dual(
      ambo(truncate(truncatedIcosahedron(scratch), scratch, 50.0f * D2R),
           scratch),
      scratch);
}

inline PolyMesh
icosahedron_snub_canonicalize_truncate033_hankin62(Arena &scratch) {
  return hankin(
      truncate(canonicalize(snub(icosahedron(scratch), scratch), scratch),
               scratch, 0.33f),
      scratch, 62.0f * D2R);
}

inline PolyMesh dodecahedron_hk35_ambo_100_hk62_ambo_100_hk43(Arena &scratch) {
  return hankin(
      canonicalize(ambo(hankin(canonicalize(ambo(hankin(dodecahedron(scratch),
                                                        scratch, 3.0f * D2R),
                                                 scratch),
                                            scratch),
                               scratch, 62.0f * D2R),
                        scratch),
                   scratch),
      scratch, 43.0f * D2R);
}

inline PolyMesh icosahedron_ambo_truncate033_hankin59(Arena &scratch) {
  return hankin(truncate(ambo(icosahedron(scratch), scratch), scratch, 0.33f),
                scratch, 59.0f * D2R);
}

inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate001_hankin59(Arena &scratch) {
  return hankin(
      truncate(
          canonicalize(ambo(truncatedIcosahedron(scratch), scratch), scratch),
          scratch, 0.01f),
      scratch, 59.0f * D2R);
}

inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate001_hankin73(Arena &scratch) {
  return hankin(
      truncate(
          canonicalize(ambo(truncatedIcosahedron(scratch), scratch), scratch),
          scratch, 0.01f),
      scratch, 73.0f * D2R);
}

inline PolyMesh truncatedOctahedron_gyro_kis_hk17(Arena &scratch) {
  return hankin(kis(gyro(truncatedOctahedron(scratch), scratch), scratch),
                scratch, 17.0f * D2R);
}

inline PolyMesh
truncatedIcosidodecahedron_bitruncate5_canonicalize_hk77(Arena &scratch) {
  return hankin(canonicalize(bitruncate(truncatedIcosidodecahedron(scratch),
                                        scratch, 0.5f),
                             scratch, 100),
                scratch, 77.0f * D2R);
}

inline PolyMesh dodecahedron_bitruncate2_canonicalize_gyro(Arena &scratch) {
  return gyro(canonicalize(bitruncate(dodecahedron(scratch), scratch, 0.2f),
                           scratch, 100),
              scratch);
}

inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate33_hk64(Arena &scratch) {
  return hankin(
      truncate(canonicalize(ambo(truncatedIcosahedron(scratch), scratch),
                            scratch, 217),
               scratch, 0.33f),
      scratch, 64.0f * D2R);
}

} // namespace IslamicStarPatterns

enum class Category { Simple, Complex };

struct Entry {
  const char *name;
  PolyMesh (*generate)(Arena &scratch);
  Category category;
};

// Centralized Registry
// Order: Platonic (0-4), Archimedean (5-17)
static constexpr Entry simple_registry[] = {

    // Platonic
    {"tetrahedron", Platonic::tetrahedron, Category::Simple},
    {"cube", Platonic::cube, Category::Simple},
    {"octahedron", Platonic::octahedron, Category::Simple},
    {"dodecahedron", Platonic::dodecahedron, Category::Simple},
    {"icosahedron", Platonic::icosahedron, Category::Simple},

    // Archimedean
    {"truncatedTetrahedron", Archimedean::truncatedTetrahedron,
     Category::Simple},
    {"cuboctahedron", Archimedean::cuboctahedron, Category::Simple},
    {"truncatedCube", Archimedean::truncatedCube, Category::Simple},
    {"truncatedOctahedron", Archimedean::truncatedOctahedron, Category::Simple},
    {"rhombicuboctahedron", Archimedean::rhombicuboctahedron, Category::Simple},
    {"truncatedCuboctahedron", Archimedean::truncatedCuboctahedron,
     Category::Simple},
    {"snubCube", Archimedean::snubCube, Category::Simple},
    {"icosidodecahedron", Archimedean::icosidodecahedron, Category::Simple},
    {"truncatedDodecahedron", Archimedean::truncatedDodecahedron,
     Category::Simple},
    {"truncatedIcosahedron", Archimedean::truncatedIcosahedron,
     Category::Simple},
    {"rhombicosidodecahedron", Archimedean::rhombicosidodecahedron,
     Category::Simple},
    {"truncatedIcosidodecahedron", Archimedean::truncatedIcosidodecahedron,
     Category::Simple},
    {"snubDodecahedron", Archimedean::snubDodecahedron, Category::Simple}};

// Islamic Star Patterns
static constexpr Entry islamic_registry[] = {
    {"truncatedIcosahedron_ambo_canonicalize_truncate33_hk64",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_canonicalize_truncate33_hk64,
     Category::Complex},
    {"dodecahedron_bitruncate2_canonicalize_gyro",
     IslamicStarPatterns::dodecahedron_bitruncate2_canonicalize_gyro,
     Category::Complex},
    {"truncatedIcosidodecahedron_bitruncate5_canonicalize_hk77",
     IslamicStarPatterns::
         truncatedIcosidodecahedron_bitruncate5_canonicalize_hk77,
     Category::Complex},
    {"truncatedOctahedron_gyro_kis_hk17",
     IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17, Category::Complex},
    {"truncatedIcosahedron_ambo_canonicalize_truncate001_hankin59",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_canonicalize_truncate001_hankin59,
     Category::Complex},
    {"truncatedIcosahedron_ambo_canonicalize_truncate001_hankin73",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_canonicalize_truncate001_hankin73,
     Category::Complex},
    {"icosahedron_ambo_truncate033_hankin59",
     IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59,
     Category::Complex},
    {"dodecahedron_hk35_ambo_100_hk62_ambo_100_hk43",
     IslamicStarPatterns::dodecahedron_hk35_ambo_100_hk62_ambo_100_hk43,
     Category::Complex},
    {"icosahedron_hk59_bitruncate033",
     IslamicStarPatterns::icosahedron_hk59_bitruncate033, Category::Complex},
    {"octahedron_hk17_ambo_hk72",
     IslamicStarPatterns::octahedron_hk17_ambo_hk72, Category::Complex},
    {"icosahedron_kis_gyro", IslamicStarPatterns::icosahedron_kis_gyro,
     Category::Complex},
    {"truncatedIcosidodecahedron_truncate05_ambo_dual",
     IslamicStarPatterns::truncatedIcosidodecahedron_truncate05_ambo_dual,
     Category::Complex},
    {"icosidodecahedron_truncate05_ambo_dual",
     IslamicStarPatterns::icosidodecahedron_truncate05_ambo_dual,
     Category::Complex},
    {"snubDodecahedron_truncate05_ambo_dual",
     IslamicStarPatterns::snubDodecahedron_truncate05_ambo_dual,
     Category::Complex},
    {"octahedron_hk34_ambo_hk72",
     IslamicStarPatterns::octahedron_hk34_ambo_hk72, Category::Complex},
    {"rhombicuboctahedron_hk63_ambo_hk63",
     IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63,
     Category::Complex},
    {"truncatedIcosahedron_hk54_ambo_hk72",
     IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72,
     Category::Complex},
    {"dodecahedron_hk54_ambo_hk72",
     IslamicStarPatterns::dodecahedron_hk54_ambo_hk72, Category::Complex},
    {"dodecahedron_hk72_ambo_dual_hk20",
     IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20, Category::Complex},
    {"truncatedIcosahedron_truncate05_ambo_dual",
     IslamicStarPatterns::truncatedIcosahedron_truncate05_ambo_dual,
     Category::Complex},
    {"icosahedron_snub_canonicalize_truncate033_hankin62",
     IslamicStarPatterns::icosahedron_snub_canonicalize_truncate033_hankin62,
     Category::Complex}};

static constexpr int NUM_ENTRIES =
    sizeof(simple_registry) / sizeof(simple_registry[0]) +
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);

namespace Collections {
static constexpr const Entry *simple_solids = simple_registry;
static constexpr int num_simple_solids =
    sizeof(simple_registry) / sizeof(simple_registry[0]);
static constexpr const Entry *islamic_solids = islamic_registry;
static constexpr int num_islamic_solids =
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);
} // namespace Collections

inline const Entry &get_entry(int index) {
  if (index < 0 || index >= NUM_ENTRIES) {
    // Fallback if out of bounds (should not happen with correct logic)
    return simple_registry[3]; // Dodecahedron
  }

  if (index < Collections::num_simple_solids) {
    return simple_registry[index];
  } else {
    return islamic_registry[index - Collections::num_simple_solids];
  }
}

inline PolyMesh get(Arena &geom, Arena &scratch, int index) {
  return finalize_solid(get_entry(index).generate(scratch), geom);
}

inline PolyMesh get_by_name(Arena &geom, Arena &scratch,
                            const std::string &name) {
  for (const auto &entry : simple_registry) {
    if (name == entry.name) {
      return finalize_solid(entry.generate(scratch), geom);
    }
  }
  for (const auto &entry : islamic_registry) {
    if (name == entry.name) {
      return finalize_solid(entry.generate(scratch), geom);
    }
  }
  return finalize_solid(Platonic::cube(scratch), geom); // Fallback
}

} // namespace Solids
