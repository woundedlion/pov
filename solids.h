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
template <typename StaticMeshT> PolyMesh to_polymesh(MemoryCtx &ctx) {
  ctx.swap_scratch();
  PolyMesh mesh;
  mesh.vertices.initialize(ctx.get_scratch_front(),
                           StaticMeshT::vertices.size());
  for (const auto &v : StaticMeshT::vertices)
    mesh.vertices.push_back(v);
  mesh.face_counts.initialize(ctx.get_scratch_front(),
                              StaticMeshT::face_counts.size());
  for (const auto &c : StaticMeshT::face_counts)
    mesh.face_counts.push_back(c);
  mesh.faces.initialize(ctx.get_scratch_front(), StaticMeshT::faces.size());
  for (const auto &f : StaticMeshT::faces)
    mesh.faces.push_back(f);
  return mesh;
}

// ==========================================================================================
// 2. PROCEDURAL GENERATORS
// ==========================================================================================

namespace Platonic {
// V=4, F=4, I=12
FLASHMEM inline PolyMesh tetrahedron(MemoryCtx &ctx) {
  return to_polymesh<Tetrahedron>(ctx);
}
// V=8, F=6, I=24
FLASHMEM inline PolyMesh cube(MemoryCtx &ctx) { return to_polymesh<Cube>(ctx); }
// V=6, F=8, I=24
FLASHMEM inline PolyMesh octahedron(MemoryCtx &ctx) {
  return to_polymesh<Octahedron>(ctx);
}
// V=20, F=12, I=60
FLASHMEM inline PolyMesh dodecahedron(MemoryCtx &ctx) {
  return to_polymesh<Dodecahedron>(ctx);
}
// V=12, F=20, I=60
FLASHMEM inline PolyMesh icosahedron(MemoryCtx &ctx) {
  return to_polymesh<Icosahedron>(ctx);
}
} // namespace Platonic

namespace Archimedean {
using namespace Platonic;
using namespace MeshOps;

// V=12, F=8, I=36
FLASHMEM inline PolyMesh truncatedTetrahedron(MemoryCtx &ctx) {
  return truncate(tetrahedron(ctx), ctx, 1.0f / 3.0f);
}
// V=12, F=14, I=48
FLASHMEM inline PolyMesh cuboctahedron(MemoryCtx &ctx) {
  return ambo(cube(ctx), ctx);
}
// V=24, F=14, I=72
FLASHMEM inline PolyMesh truncatedCube(MemoryCtx &ctx) {
  return truncate(cube(ctx), ctx, 1.0f / (2.0f + SQRT2));
}
// V=24, F=14, I=72
FLASHMEM inline PolyMesh truncatedOctahedron(MemoryCtx &ctx) {
  return truncate(octahedron(ctx), ctx, 1.0f / 3.0f);
}
// V=24, F=26, I=96
FLASHMEM inline PolyMesh rhombicuboctahedron(MemoryCtx &ctx) {
  return expand(cube(ctx), ctx);
}
// V=24, F=14, I=72
FLASHMEM inline PolyMesh truncatedCuboctahedron(MemoryCtx &ctx) {
  return canonicalize(bitruncate(cube(ctx), ctx, 1.0f / (2.0f + SQRT2)), ctx,
                      50);
}
// V=24, F=38, I=120
FLASHMEM inline PolyMesh snubCube(MemoryCtx &ctx) {
  return canonicalize(snub(cube(ctx), ctx, T_SNUB_CUBE, 0.28f), ctx, 50);
}
// V=30, F=32, I=120
FLASHMEM inline PolyMesh icosidodecahedron(MemoryCtx &ctx) {
  return ambo(dodecahedron(ctx), ctx);
}
// V=60, F=32, I=180
FLASHMEM inline PolyMesh truncatedDodecahedron(MemoryCtx &ctx) {
  return truncate(dodecahedron(ctx), ctx, T_TRUNC_ICOS);
}
// V=60, F=32, I=180
FLASHMEM inline PolyMesh truncatedIcosahedron(MemoryCtx &ctx) {
  return truncate(icosahedron(ctx), ctx, 1.0f / 3.0f);
}
// V=60, F=62, I=240
FLASHMEM inline PolyMesh rhombicosidodecahedron(MemoryCtx &ctx) {
  return canonicalize(expand(dodecahedron(ctx), ctx), ctx, 50);
}
// V=60, F=32, I=180
FLASHMEM inline PolyMesh truncatedIcosidodecahedron(MemoryCtx &ctx) {
  return canonicalize(bitruncate(dodecahedron(ctx), ctx, 1.0f / (2.0f + PHI)),
                      ctx, 50);
}
// V=60, F=92, I=300
FLASHMEM inline PolyMesh snubDodecahedron(MemoryCtx &ctx) {
  return canonicalize(snub(dodecahedron(ctx), ctx, 0.5f), ctx, 50);
}
} // namespace Archimedean

namespace Catalan {
using namespace Archimedean;
using namespace MeshOps;

// Dual of Truncated Tetrahedron
FLASHMEM inline PolyMesh triakisTetrahedron(MemoryCtx &ctx) {
  return dual(truncatedTetrahedron(ctx), ctx);
}
// Dual of Cuboctahedron
FLASHMEM inline PolyMesh rhombicDodecahedron(MemoryCtx &ctx) {
  return dual(cuboctahedron(ctx), ctx);
}
// Dual of Truncated Cube
FLASHMEM inline PolyMesh triakisOctahedron(MemoryCtx &ctx) {
  return dual(truncatedCube(ctx), ctx);
}
// Dual of Truncated Octahedron
FLASHMEM inline PolyMesh tetrakisHexahedron(MemoryCtx &ctx) {
  return dual(truncatedOctahedron(ctx), ctx);
}
// Dual of Rhombicuboctahedron
FLASHMEM inline PolyMesh deltoidalIcositetrahedron(MemoryCtx &ctx) {
  return dual(rhombicuboctahedron(ctx), ctx);
}
// Dual of Truncated Cuboctahedron
FLASHMEM inline PolyMesh disdyakisDodecahedron(MemoryCtx &ctx) {
  return dual(truncatedCuboctahedron(ctx), ctx);
}
// Dual of Snub Cube
FLASHMEM inline PolyMesh pentagonalIcositetrahedron(MemoryCtx &ctx) {
  return dual(snubCube(ctx), ctx);
}
// Dual of Icosidodecahedron
FLASHMEM inline PolyMesh rhombicTriacontahedron(MemoryCtx &ctx) {
  return dual(icosidodecahedron(ctx), ctx);
}
// Dual of Truncated Dodecahedron
FLASHMEM inline PolyMesh triakisIcosahedron(MemoryCtx &ctx) {
  return dual(truncatedDodecahedron(ctx), ctx);
}
// Dual of Truncated Icosahedron
FLASHMEM inline PolyMesh pentakisDodecahedron(MemoryCtx &ctx) {
  return dual(truncatedIcosahedron(ctx), ctx);
}
// Dual of Rhombicosidodecahedron
FLASHMEM inline PolyMesh deltoidalHexecontahedron(MemoryCtx &ctx) {
  return dual(rhombicosidodecahedron(ctx), ctx);
}
// Dual of Truncated Icosidodecahedron
FLASHMEM inline PolyMesh disdyakisTriacontahedron(MemoryCtx &ctx) {
  return dual(truncatedIcosidodecahedron(ctx), ctx);
}
// Dual of Snub Dodecahedron
FLASHMEM inline PolyMesh pentagonalHexecontahedron(MemoryCtx &ctx) {
  return dual(snubDodecahedron(ctx), ctx);
}
} // namespace Catalan

namespace IslamicStarPatterns {
using namespace Platonic;
using namespace Archimedean;
using namespace MeshOps;

static constexpr float D2R = PI_F / 180.0f;

// V=576, F=434, I=2016
FLASHMEM inline PolyMesh
cube_canonicalize_bitruncate33_canonicalize_hk68_expand5(MemoryCtx &ctx) {
  return expand(
      hankin(canonicalize(
                 bitruncate(canonicalize(cube(ctx), ctx, 100), ctx, 0.33f), ctx,
                 100),
             ctx, 67.5f * D2R),
      ctx, 0.5f);
}

// V=240, F=146, I=768
FLASHMEM inline PolyMesh icosahedron_hk59_bitruncate033(MemoryCtx &ctx) {
  return hankin(bitruncate(icosahedron(ctx), ctx, 0.33f), ctx, 59.0f * D2R);
}

// V=432, F=50, I=960
FLASHMEM inline PolyMesh octahedron_hk17_ambo_hk72(MemoryCtx &ctx) {
  return hankin(ambo(hankin(octahedron(ctx), ctx, 17.0f * D2R), ctx), ctx,
                73.0f * D2R);
}
// V=272, F=452, I=1444
FLASHMEM inline PolyMesh icosahedron_kis_gyro(MemoryCtx &ctx) {
  return gyro(kis(icosahedron(ctx), ctx), ctx);
}

// V=542, F=540, I=2160
FLASHMEM inline PolyMesh
truncatedIcosidodecahedron_truncate05_ambo_dual(MemoryCtx &ctx) {
  return dual(
      ambo(truncate(truncatedIcosidodecahedron(ctx), ctx, 50.0f * D2R), ctx),
      ctx);
}

// V=182, F=180, I=720
FLASHMEM inline PolyMesh
icosidodecahedron_truncate05_ambo_dual(MemoryCtx &ctx) {
  return dual(ambo(truncate(icosidodecahedron(ctx), ctx, 5.0f * D2R), ctx),
              ctx);
}
// V=452, F=450, I=1800
FLASHMEM inline PolyMesh snubDodecahedron_truncate05_ambo_dual(MemoryCtx &ctx) {
  return dual(ambo(truncate(snubDodecahedron(ctx), ctx, 5.0f * D2R), ctx), ctx);
}

// V=432, F=50, I=960
FLASHMEM inline PolyMesh octahedron_hk34_ambo_hk72(MemoryCtx &ctx) {
  return hankin(ambo(hankin(octahedron(ctx), ctx, 34.0f * D2R), ctx), ctx,
                72.0f * D2R);
}

// V=1728, F=194, I=3840
FLASHMEM inline PolyMesh rhombicuboctahedron_hk63_ambo_hk63(MemoryCtx &ctx) {
  return hankin(ambo(hankin(rhombicuboctahedron(ctx), ctx, 63.0f * D2R), ctx),
                ctx, 63.0f * D2R);
}

// V=3240, F=362, I=7200
FLASHMEM inline PolyMesh truncatedIcosahedron_hk54_ambo_hk72(MemoryCtx &ctx) {
  return hankin(ambo(hankin(truncatedIcosahedron(ctx), ctx, 54.0f * D2R), ctx),
                ctx, 72.0f * D2R);
}

// V=1080, F=122, I=2400
FLASHMEM inline PolyMesh dodecahedron_hk54_ambo_hk72(MemoryCtx &ctx) {
  return hankin(ambo(hankin(dodecahedron(ctx), ctx, 54.0f * D2R), ctx), ctx,
                72.0f * D2R);
}

// V=1122, F=164, I=2568
FLASHMEM inline PolyMesh dodecahedron_hk72_ambo_dual_hk20(MemoryCtx &ctx) {
  return hankin(
      dual(ambo(hankin(dodecahedron(ctx), ctx, 72.0f * D2R), ctx), ctx), ctx,
      20.0f * D2R);
}

// V=272, F=270, I=1080
FLASHMEM inline PolyMesh
truncatedIcosahedron_truncate05_ambo_dual(MemoryCtx &ctx) {
  return dual(ambo(truncate(truncatedIcosahedron(ctx), ctx, 50.0f * D2R), ctx),
              ctx);
}

// V=2100, F=302, I=4800
FLASHMEM inline PolyMesh
icosahedron_snub_canonicalize_truncate033_hankin62(MemoryCtx &ctx) {
  return hankin(
      truncate(canonicalize(snub(icosahedron(ctx), ctx), ctx), ctx, 0.33f), ctx,
      62.0f * D2R);
}

// V=3240, F=1082, I=8640
FLASHMEM inline PolyMesh
dodecahedron_hk35_ambo_hk62_ambo_canonicalize_hk43(MemoryCtx &ctx) {
  return hankin(
      canonicalize(
          ambo(hankin(ambo(hankin(dodecahedron(ctx), ctx, 35.0f * D2R), ctx),
                      ctx, 62.0f * D2R),
               ctx),
          ctx, 100),
      ctx, 43.0f * D2R);
}

// V=840, F=122, I=1920
FLASHMEM inline PolyMesh icosahedron_ambo_truncate033_hankin59(MemoryCtx &ctx) {
  return hankin(truncate(ambo(icosahedron(ctx), ctx), ctx, 0.33f), ctx,
                59.0f * D2R);
}

// V=2520, F=362, I=5760
FLASHMEM inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate001_hankin59(MemoryCtx &ctx) {
  return hankin(
      truncate(canonicalize(ambo(truncatedIcosahedron(ctx), ctx), ctx), ctx,
               0.01f),
      ctx, 59.0f * D2R);
}

// V=2520, F=362, I=5760
FLASHMEM inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate001_hankin73(MemoryCtx &ctx) {
  return hankin(
      truncate(canonicalize(ambo(truncatedIcosahedron(ctx), ctx), ctx), ctx,
               0.01f),
      ctx, 73.0f * D2R);
}

// V=2452, F=294, I=5488
FLASHMEM inline PolyMesh truncatedOctahedron_gyro_kis_hk17(MemoryCtx &ctx) {
  return hankin(kis(gyro(truncatedOctahedron(ctx), ctx), ctx), ctx,
                17.0f * D2R);
}

// V=2520, F=362, I=5760
FLASHMEM inline PolyMesh
truncatedIcosidodecahedron_bitruncate5_canonicalize_hk77(MemoryCtx &ctx) {
  return hankin(
      canonicalize(bitruncate(truncatedIcosidodecahedron(ctx), ctx, 0.5f), ctx,
                   100),
      ctx, 77.0f * D2R);
}

// V=272, F=452, I=1444
FLASHMEM inline PolyMesh
dodecahedron_bitruncate2_canonicalize_gyro(MemoryCtx &ctx) {
  return gyro(canonicalize(bitruncate(dodecahedron(ctx), ctx, 0.2f), ctx, 100),
              ctx);
}

// V=2520, F=362, I=5760
FLASHMEM inline PolyMesh
truncatedIcosahedron_ambo_canonicalize_truncate33_hk64(MemoryCtx &ctx) {
  return hankin(
      truncate(canonicalize(ambo(truncatedIcosahedron(ctx), ctx), ctx, 217),
               ctx, 0.33f),
      ctx, 64.0f * D2R);
}

// V=1080, F=362, I=1440
FLASHMEM inline PolyMesh
dodecahedron_ambo_bitruncate33_canonicalize_hk66(MemoryCtx &ctx) {
  return hankin(
      canonicalize(bitruncate(ambo(dodecahedron(ctx), ctx), ctx, 0.33f), ctx,
                   100),
      ctx, 66.0f * D2R);
}

} // namespace IslamicStarPatterns

enum class Category { Simple, Complex };

struct Entry {
  const char *name;
  PolyMesh (*generate)(MemoryCtx &ctx);
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

FLASHMEM inline PolyMesh get(Arena &geom, MemoryCtx &ctx, int index) {
  return finalize_solid(get_entry(index).generate(ctx), geom);
}

FLASHMEM inline PolyMesh get_by_name(Arena &geom, MemoryCtx &ctx,
                                     std::string_view name) {
  for (const auto &entry : simple_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(ctx), geom);
  }
  for (const auto &entry : catalan_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(ctx), geom);
  }
  for (const auto &entry : islamic_registry) {
    if (name == entry.name)
      return finalize_solid(entry.generate(ctx), geom);
  }
  return finalize_solid(Platonic::cube(ctx), geom);
}

} // namespace Solids
