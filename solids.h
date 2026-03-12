/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "geometry.h"
#include "mesh.h" // For MeshOps
#include "hankin.h"
#include "conway.h"
#include <cmath>
#include <string_view>
#include <span>
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

static constexpr int MAX_VERTS = 8700;
static constexpr int MAX_FACES = 1000;
static constexpr int MAX_INDICES = 20000;

FLASHMEM inline PolyMesh finalize_solid(const PolyMesh &temp, Arena &geom) {
  PolyMesh final_mesh;
  final_mesh.vertices.bind(geom, temp.vertices.size());
  for (size_t i = 0; i < temp.vertices.size(); ++i)
    final_mesh.vertices.push_back(temp.vertices[i]);
  final_mesh.face_counts.bind(geom, temp.face_counts.size());
  for (size_t i = 0; i < temp.face_counts.size(); ++i)
    final_mesh.face_counts.push_back(temp.face_counts[i]);
  final_mesh.faces.bind(geom, temp.faces.size());
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
template <typename StaticMeshT> PolyMesh to_polymesh(Arena &target) {
  PolyMesh mesh;
  mesh.vertices.bind(target, StaticMeshT::vertices.size());
  for (const auto &v : StaticMeshT::vertices)
    mesh.vertices.push_back(v);
  mesh.face_counts.bind(target, StaticMeshT::face_counts.size());
  for (const auto &c : StaticMeshT::face_counts)
    mesh.face_counts.push_back(c);
  mesh.faces.bind(target, StaticMeshT::faces.size());
  for (const auto &f : StaticMeshT::faces)
    mesh.faces.push_back(f);
  return mesh;
}

/// Fluent builder for chaining Conway operators with automatic arena swapping.
class SolidBuilder {
  PolyMesh mesh_;
  Arena* a_;
  Arena* b_;
public:
  SolidBuilder(PolyMesh seed, Arena& a, Arena& b)
    : mesh_(std::move(seed)), a_(&a), b_(&b) {}

  SolidBuilder& dual() { mesh_ = MeshOps::dual(mesh_, *a_, *b_); std::swap(a_, b_); return *this; }
  SolidBuilder& kis() { mesh_ = MeshOps::kis(mesh_, *a_, *b_); std::swap(a_, b_); return *this; }
  SolidBuilder& ambo() { mesh_ = MeshOps::ambo(mesh_, *a_, *b_); std::swap(a_, b_); return *this; }
  SolidBuilder& truncate(float t = 0.25f) { mesh_ = MeshOps::truncate(mesh_, *a_, *b_, t); std::swap(a_, b_); return *this; }
  SolidBuilder& expand(float t = 2.0f - 1.414213562373095f) { mesh_ = MeshOps::expand(mesh_, *a_, *b_, t); std::swap(a_, b_); return *this; }
  SolidBuilder& chamfer(float t = 0.5f) { mesh_ = MeshOps::chamfer(mesh_, *a_, *b_, t); std::swap(a_, b_); return *this; }
  SolidBuilder& snub(float t = 0.5f, float twist = 0.0f) { mesh_ = MeshOps::snub(mesh_, *a_, *b_, t, twist); std::swap(a_, b_); return *this; }
  SolidBuilder& gyro() { mesh_ = MeshOps::gyro(mesh_, *a_, *b_); std::swap(a_, b_); return *this; }
  SolidBuilder& bitruncate(float t = 1.0f/3.0f) { mesh_ = MeshOps::bitruncate(mesh_, *a_, *b_, t); std::swap(a_, b_); return *this; }
  SolidBuilder& canonicalize(int iterations = 8) { mesh_ = MeshOps::canonicalize(mesh_, *a_, *b_, iterations); std::swap(a_, b_); return *this; }
  SolidBuilder& hankin(float angle) { mesh_ = MeshOps::hankin(mesh_, *a_, *b_, angle); std::swap(a_, b_); return *this; }

  PolyMesh build() { return std::move(mesh_); }
};

// ==========================================================================================
// 2. PROCEDURAL GENERATORS
// ==========================================================================================

namespace Platonic {
// V=4, F=4, I=12
FLASHMEM inline PolyMesh tetrahedron(Arena &a, Arena &b) {
  return to_polymesh<Tetrahedron>(a);
}
// V=8, F=6, I=24
FLASHMEM inline PolyMesh cube(Arena &a, Arena &b) { return to_polymesh<Cube>(a); }
// V=6, F=8, I=24
FLASHMEM inline PolyMesh octahedron(Arena &a, Arena &b) {
  return to_polymesh<Octahedron>(a);
}
// V=20, F=12, I=60
FLASHMEM inline PolyMesh dodecahedron(Arena &a, Arena &b) {
  return to_polymesh<Dodecahedron>(a);
}
// V=12, F=20, I=60
FLASHMEM inline PolyMesh icosahedron(Arena &a, Arena &b) {
  return to_polymesh<Icosahedron>(a);
}
} // namespace Platonic

namespace Archimedean {
using namespace Platonic;

FLASHMEM inline PolyMesh truncatedTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(tetrahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
}
FLASHMEM inline PolyMesh cuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).ambo().build();
}
FLASHMEM inline PolyMesh truncatedCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).truncate(1.0f / (2.0f + SQRT2)).build();
}
FLASHMEM inline PolyMesh truncatedOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
}
FLASHMEM inline PolyMesh rhombicuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).expand().build();
}
FLASHMEM inline PolyMesh truncatedCuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).bitruncate(1.0f / (2.0f + SQRT2)).canonicalize(50).build();
}
FLASHMEM inline PolyMesh snubCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).snub(T_SNUB_CUBE, 0.28f).canonicalize(50).build();
}
FLASHMEM inline PolyMesh icosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).ambo().build();
}
FLASHMEM inline PolyMesh truncatedDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).truncate(T_TRUNC_ICOS).build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
}
FLASHMEM inline PolyMesh rhombicosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).expand().canonicalize(50).build();
}
FLASHMEM inline PolyMesh truncatedIcosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).bitruncate(1.0f / (2.0f + PHI)).canonicalize(50).build();
}
FLASHMEM inline PolyMesh snubDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).snub(0.5f).canonicalize(50).build();
}
} // namespace Archimedean

namespace Catalan {
using namespace Archimedean;

FLASHMEM inline PolyMesh triakisTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedTetrahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh rhombicDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(cuboctahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh triakisOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCube(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh tetrakisHexahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedOctahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh deltoidalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicuboctahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh disdyakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCuboctahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh pentagonalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubCube(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh rhombicTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosidodecahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh triakisIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedDodecahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh pentakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh deltoidalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicosidodecahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh disdyakisTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosidodecahedron(a, b), a, b).dual().build();
}
FLASHMEM inline PolyMesh pentagonalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubDodecahedron(a, b), a, b).dual().build();
}
} // namespace Catalan

namespace IslamicStarPatterns {
using namespace Platonic;
using namespace Archimedean;

static constexpr float D2R = PI_F / 180.0f;

FLASHMEM inline PolyMesh cube_canonicalize_bitruncate33_canonicalize_hk68_expand5(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).canonicalize(100).bitruncate(0.33f).canonicalize(100).hankin(67.5f * D2R).expand(0.5f).build();
}
FLASHMEM inline PolyMesh icosahedron_hk59_bitruncate033(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).bitruncate(0.33f).hankin(59.0f * D2R).build();
}
FLASHMEM inline PolyMesh octahedron_hk17_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).hankin(17.0f * D2R).ambo().hankin(73.0f * D2R).build();
}
FLASHMEM inline PolyMesh icosahedron_kis_gyro(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).kis().gyro().build();
}
FLASHMEM inline PolyMesh truncatedIcosidodecahedron_truncate05_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosidodecahedron(a, b), a, b).truncate(50.0f * D2R).ambo().dual().build();
}
FLASHMEM inline PolyMesh icosidodecahedron_truncate05_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(icosidodecahedron(a, b), a, b).truncate(5.0f * D2R).ambo().dual().build();
}
FLASHMEM inline PolyMesh snubDodecahedron_truncate05_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(snubDodecahedron(a, b), a, b).truncate(5.0f * D2R).ambo().dual().build();
}
FLASHMEM inline PolyMesh octahedron_hk34_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).hankin(34.0f * D2R).ambo().hankin(72.0f * D2R).build();
}
FLASHMEM inline PolyMesh rhombicuboctahedron_hk63_ambo_hk63(Arena &a, Arena &b) {
  return SolidBuilder(rhombicuboctahedron(a, b), a, b).hankin(63.0f * D2R).ambo().hankin(63.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron_hk54_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).hankin(54.0f * D2R).ambo().hankin(72.0f * D2R).build();
}
FLASHMEM inline PolyMesh dodecahedron_hk54_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).hankin(54.0f * D2R).ambo().hankin(72.0f * D2R).build();
}
FLASHMEM inline PolyMesh dodecahedron_hk72_ambo_dual_hk20(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).hankin(72.0f * D2R).ambo().dual().hankin(20.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron_truncate05_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).truncate(50.0f * D2R).ambo().dual().build();
}
FLASHMEM inline PolyMesh icosahedron_snub_canonicalize_truncate033_hankin62(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).snub().canonicalize().truncate(0.33f).hankin(62.0f * D2R).build();
}
FLASHMEM inline PolyMesh dodecahedron_hk35_ambo_hk62_ambo_canonicalize_hk43(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).hankin(35.0f * D2R).ambo().hankin(62.0f * D2R).ambo().canonicalize(100).hankin(43.0f * D2R).build();
}
FLASHMEM inline PolyMesh icosahedron_ambo_truncate033_hankin59(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).ambo().truncate(0.33f).hankin(59.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron_ambo_canonicalize_truncate001_hankin59(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).ambo().canonicalize().truncate(0.01f).hankin(59.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron_ambo_canonicalize_truncate001_hankin73(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).ambo().canonicalize().truncate(0.01f).hankin(73.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedOctahedron_gyro_kis_hk17(Arena &a, Arena &b) {
  return SolidBuilder(truncatedOctahedron(a, b), a, b).gyro().kis().hankin(17.0f * D2R).build();
}
FLASHMEM inline PolyMesh truncatedIcosidodecahedron_bitruncate5_canonicalize_hk77(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosidodecahedron(a, b), a, b).bitruncate(0.5f).canonicalize(100).hankin(77.0f * D2R).build();
}
FLASHMEM inline PolyMesh dodecahedron_bitruncate2_canonicalize_gyro(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).bitruncate(0.2f).canonicalize(100).gyro().build();
}
FLASHMEM inline PolyMesh truncatedIcosahedron_ambo_canonicalize_truncate33_hk64(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).ambo().canonicalize(217).truncate(0.33f).hankin(64.0f * D2R).build();
}
FLASHMEM inline PolyMesh dodecahedron_ambo_bitruncate33_canonicalize_hk66(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).ambo().bitruncate(0.33f).canonicalize(100).hankin(66.0f * D2R).build();
}
} // namespace IslamicStarPatterns

enum class Category { Simple, Complex };

struct Entry {
  const char *name;
  PolyMesh (*generate)(Arena &a, Arena &b);
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
     Category::Simple}};

static constexpr Entry catalan_registry[] = {
    // Catalan
    {"triakisTetrahedron", Catalan::triakisTetrahedron, Category::Simple},
    {"rhombicDodecahedron", Catalan::rhombicDodecahedron, Category::Simple},
    {"triakisOctahedron", Catalan::triakisOctahedron, Category::Simple},
    {"tetrakisHexahedron", Catalan::tetrakisHexahedron, Category::Simple},
    {"deltoidalIcositetrahedron", Catalan::deltoidalIcositetrahedron,
     Category::Simple},
    {"disdyakisDodecahedron", Catalan::disdyakisDodecahedron, Category::Simple},
    {"pentagonalIcositetrahedron", Catalan::pentagonalIcositetrahedron,
     Category::Simple},
    {"rhombicTriacontahedron", Catalan::rhombicTriacontahedron,
     Category::Simple},
    {"triakisIcosahedron", Catalan::triakisIcosahedron, Category::Simple},
    {"pentakisDodecahedron", Catalan::pentakisDodecahedron, Category::Simple},
    {"deltoidalHexecontahedron", Catalan::deltoidalHexecontahedron,
     Category::Simple},
    {"disdyakisTriacontahedron", Catalan::disdyakisTriacontahedron,
     Category::Simple},
    {"pentagonalHexecontahedron", Catalan::pentagonalHexecontahedron,
     Category::Simple}};

// Islamic Star Patterns
static constexpr Entry islamic_registry[] = {
    {"cube_canonicalize_bitruncate33_canonicalize_hk68_expand5",
     IslamicStarPatterns::
         cube_canonicalize_bitruncate33_canonicalize_hk68_expand5,
     Category::Complex},
    {"dodecahedron_ambo_bitruncate33_canonicalize_hk66",
     IslamicStarPatterns::dodecahedron_ambo_bitruncate33_canonicalize_hk66,
     Category::Complex},
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
    {"dodecahedron_hk35_ambo_hk62_ambo_canonicalize_hk43",
     IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_canonicalize_hk43,
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
    sizeof(catalan_registry) / sizeof(catalan_registry[0]) +
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);

namespace Collections {
inline std::span<const Entry> get_platonic_solids() {
  return std::span<const Entry>(simple_registry, 5);
}
inline std::span<const Entry> get_archimedean_solids() {
  return std::span<const Entry>(simple_registry + 5, 11);
}
inline std::span<const Entry> get_simple_solids() {
  return std::span<const Entry>(simple_registry);
}
inline std::span<const Entry> get_catalan_solids() {
  return std::span<const Entry>(catalan_registry);
}
inline std::span<const Entry> get_islamic_solids() {
  return std::span<const Entry>(islamic_registry);
}
} // namespace Collections

inline const Entry &get_entry(size_t index) {
  if (index < 0 || index >= NUM_ENTRIES)
    return simple_registry[3];

  if (index < std::size(simple_registry))
    return simple_registry[index];

  if (index < std::size(simple_registry) + std::size(catalan_registry))
    return catalan_registry[index - std::size(simple_registry)];

  return islamic_registry[index - (std::size(simple_registry) +
                                   std::size(catalan_registry))];
}

FLASHMEM inline PolyMesh get(Arena &geom, Arena &a, Arena &b, int index) {
  return finalize_solid(get_entry(index).generate(a, b), geom);
}

FLASHMEM inline PolyMesh get_by_name(Arena &geom, Arena &a, Arena &b,
                                     std::string_view name) {
  for (const auto &entry : simple_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(a, b), geom);
  }
  for (const auto &entry : catalan_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(a, b), geom);
  }
  for (const auto &entry : islamic_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(a, b), geom);
  }
  return finalize_solid(Platonic::cube(a, b), geom);
}

} // namespace Solids
