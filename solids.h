/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include <cmath>
#include <string>
#include <vector>
#include <map>

// --- Constants for Procedural Generation ---
static constexpr float SQRT2 = 1.414213562373095f;
static constexpr float TRIBONACCI = (1.0f + 1.839286755214161f + 1.839286755214161f) / 3.0f; // Approx
// Real tribonacci constant t: t^3 - t^2 - t - 1 = 0. ~1.839286755214161
static constexpr float TRIBONACCI_CONST = 1.839286755214161f; 
static constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
static constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + PHI);

namespace Solids {

  // ==========================================================================================
  // 1. DATA DEFINITIONS (Hardcoded Platonic Solids)
  // ==========================================================================================

  struct Tetrahedron {
    static constexpr int NUM_VERTS = 4;
    static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f),
      Vector(0.5773502691896258f, -0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, 0.5773502691896258f, -0.5773502691896258f),
      Vector(-0.5773502691896258f, -0.5773502691896258f, 0.5773502691896258f)
    };
    static constexpr int NUM_FACES = 4;
    static constexpr std::array<uint8_t, NUM_FACES> face_counts = { 3, 3, 3, 3 };
    static constexpr std::array<int, 12> faces = { 0, 3, 1, 0, 2, 3, 0, 1, 2, 1, 3, 2 };
  };

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
      Vector(-0.5773502691896258f, 0.5773502691896258f, 0.5773502691896258f)
    };
    static constexpr int NUM_FACES = 6;
    static constexpr std::array<uint8_t, NUM_FACES> face_counts = { 4, 4, 4, 4, 4, 4 };
    static constexpr std::array<int, 24> faces = {
      0, 3, 2, 1, 0, 1, 5, 4, 0, 4, 7, 3, 6, 5, 1, 2, 6, 2, 3, 7, 6, 7, 4, 5
    };
  };

  struct Octahedron {
    static constexpr int NUM_VERTS = 6;
    static constexpr std::array<Vector, NUM_VERTS> vertices = {
      Vector(1.0000000000000000f, 0.0000000000000000f, 0.0000000000000000f),
      Vector(-1.0000000000000000f, 0.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, 1.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, -1.0000000000000000f, 0.0000000000000000f),
      Vector(0.0000000000000000f, 0.0000000000000000f, 1.0000000000000000f),
      Vector(0.0000000000000000f, 0.0000000000000000f, -1.0000000000000000f)
    };
    static constexpr int NUM_FACES = 8;
    static constexpr std::array<uint8_t, NUM_FACES> face_counts = { 3, 3, 3, 3, 3, 3, 3, 3 };
    static constexpr std::array<int, 24> faces = {
      4, 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 5, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3
    };
  };

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
      Vector(-0.8506508083520400f, -0.5257311121191336f, 0.0000000000000000f)
    };
    static constexpr int NUM_FACES = 20;
    static constexpr std::array<uint8_t, NUM_FACES> face_counts = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
    static constexpr std::array<int, 60> faces = {
      0, 1, 4, 0, 4, 9, 9, 4, 5, 4, 8, 5, 4, 1, 8, 8, 1, 10, 8, 10, 3, 5, 8, 3, 5, 3, 2, 2, 3, 7, 7, 3, 10, 7, 10, 6, 7, 6, 11, 11, 6, 0, 0, 6, 1, 6, 10, 1, 9, 11, 0, 9, 2, 11, 9, 5, 2, 7, 11, 2
    };
  };

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
      Vector(0.0000000000000000f, -0.3568220897730897f, -0.9341723589627157f)
    };
    static constexpr int NUM_FACES = 12;
    static constexpr std::array<uint8_t, NUM_FACES> face_counts = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };
    static constexpr std::array<int, 60> faces = {
      0, 8, 9, 4, 16, 0, 12, 13, 1, 8, 0, 16, 17, 2, 12, 8, 1, 18, 5, 9,
      12, 2, 10, 3, 13, 16, 4, 14, 6, 17, 9, 5, 15, 14, 4, 6, 11, 10, 2, 17,
      3, 19, 18, 1, 13, 7, 15, 5, 18, 19, 7, 11, 6, 14, 15, 7, 19, 3, 10, 11
    };
  };

  // Helper to convert Static Mesh Data to Dynamic PolyMesh
  template<typename StaticMeshT>
  PolyMesh to_polymesh() {
    PolyMesh mesh;
    mesh.vertices.assign(StaticMeshT::vertices.begin(), StaticMeshT::vertices.end());
    
    int offset = 0;
    for (size_t i = 0; i < StaticMeshT::face_counts.size(); ++i) {
      int count = StaticMeshT::face_counts[i];
      std::vector<int> face;
      for (int k = 0; k < count; ++k) {
        face.push_back(StaticMeshT::faces[offset + k]);
      }
      mesh.faces.push_back(face);
      offset += count;
    }
    return mesh;
  }

  // ==========================================================================================
  // 2. PROCEDURAL GENERATORS
  // ==========================================================================================

  namespace Platonic {
    inline PolyMesh tetrahedron() { return to_polymesh<Tetrahedron>(); }
    inline PolyMesh cube() { return to_polymesh<Cube>(); }
    inline PolyMesh octahedron() { return to_polymesh<Octahedron>(); }
    inline PolyMesh dodecahedron() { return to_polymesh<Dodecahedron>(); }
    inline PolyMesh icosahedron() { return to_polymesh<Icosahedron>(); }
  }

  namespace Archimedean {
    using namespace MeshOps;
    using namespace Platonic;

    inline PolyMesh truncatedTetrahedron() { return truncate(tetrahedron(), 1.0f/3.0f); }
    inline PolyMesh cuboctahedron() { return ambo(cube()); }
    inline PolyMesh truncatedCube() { return truncate(cube(), 1.0f/(2.0f + SQRT2)); }
    inline PolyMesh truncatedOctahedron() { return truncate(octahedron(), 1.0f/3.0f); }
    inline PolyMesh rhombicuboctahedron() { return expand(cube()); }
    inline PolyMesh truncatedCuboctahedron() { return canonicalize(bitruncate(cube(), 1.0f/(2.0f + SQRT2)), 50); }
    inline PolyMesh snubCube() { return canonicalize(snub(cube(), T_SNUB_CUBE, 0.28f), 50); }
    inline PolyMesh icosidodecahedron() { return ambo(dodecahedron()); }
    inline PolyMesh truncatedDodecahedron() { return truncate(dodecahedron(), T_TRUNC_ICOS); }
    inline PolyMesh truncatedIcosahedron() { return truncate(icosahedron(), 1.0f/3.0f); }
    inline PolyMesh rhombicosidodecahedron() { return canonicalize(expand(dodecahedron()), 50); }
    inline PolyMesh truncatedIcosidodecahedron() { return canonicalize(bitruncate(dodecahedron(), 1.0f/(2.0f + PHI)), 50); }
    inline PolyMesh snubDodecahedron() { return canonicalize(snub(dodecahedron(), 0.5f), 50); }
  }

  namespace IslamicStarPatterns {
    using namespace MeshOps;
    using namespace Platonic;
    using namespace Archimedean;

    static constexpr float D2R = PI_F / 180.0f;

    // Helper for bitruncate (truncate(ambo(m), t))
    inline PolyMesh bitruncate(const PolyMesh& m, float t) {
      return truncate(ambo(m), t);
    }

    // Base: icosahedron, Ops: Hk(59), Bitruncate(0.33)
    // JS: icosahedron_hk59_bitruncate033: () => MeshOps.hankin(MeshOps.bitruncate(Solids.PlatonicSolids.icosahedron(), 0.33), 59 * (Math.PI / 180))
    // Note: JS order is hankin(bitruncate(...), 59).
    inline PolyMesh icosahedron_hk59_bitruncate033() {
        return hankin(bitruncate(icosahedron(), 0.33f), 59.0f * D2R);
    }

    // Base: octahedron, Ops: Hk(17), Ambo, Hk(73) (Note: 73 deg in JS)
    // JS: octahedron_hk17_ambo_hk72: ... hankin(..., 73 * Math.PI / 180)
    // Name says 72, code says 73. We match code.
    inline PolyMesh octahedron_hk17_ambo_hk72() {
        return hankin(ambo(hankin(octahedron(), 17.0f * D2R)), 73.0f * D2R);
    }

    // Base: icosahedron, Ops: Kis, Gyro
    inline PolyMesh icosahedron_kis_gyro() {
        return gyro(kis(icosahedron()));
    }

    // Base: truncatedIcosidodecahedron, Ops: Tr(50), Ambo, Dual
    inline PolyMesh truncatedIcosidodecahedron_truncate05_ambo_dual() {
        return dual(ambo(truncate(truncatedIcosidodecahedron(), 50.0f * D2R)));
    }

    // Base: icosidodecahedron, Ops: Tr(5), Ambo, Dual
    inline PolyMesh icosidodecahedron_truncate05_ambo_dual() {
        return dual(ambo(truncate(icosidodecahedron(), 5.0f * D2R)));
    }

    // Base: snubDodecahedron, Ops: Tr(5), Ambo, Dual
    inline PolyMesh snubDodecahedron_truncate05_ambo_dual() {
        return dual(ambo(truncate(snubDodecahedron(), 5.0f * D2R)));
    }

    // Base: octahedron, Ops: Hk(34), Ambo, Hk(72)
    inline PolyMesh octahedron_hk34_ambo_hk72() {
        return hankin(ambo(hankin(octahedron(), 34.0f * D2R)), 72.0f * D2R);
    }

    // Base: rhombicuboctahedron, Ops: Hk(63), Ambo, Hk(63)
    // JS: rhombicuboctahedron_hk63_ambo_hk63
    inline PolyMesh rhombicuboctahedron_hk63_ambo_hk63() {
        return hankin(ambo(hankin(rhombicuboctahedron(), 63.0f * D2R)), 63.0f * D2R);
    }

    // Base: truncatedIcosahedron, Ops: Hk(54), Ambo, Hk(72)
    inline PolyMesh truncatedIcosahedron_hk54_ambo_hk72() {
        return hankin(ambo(hankin(truncatedIcosahedron(), 54.0f * D2R)), 72.0f * D2R);
    }

    // Base: dodecahedron, Ops: Hk(54), Ambo, Hk(72)
    inline PolyMesh dodecahedron_hk54_ambo_hk72() {
        return hankin(ambo(hankin(dodecahedron(), 54.0f * D2R)), 72.0f * D2R);
    }

    // Base: dodecahedron, Ops: Hk(72), Ambo, Dual, Hk(20)
    inline PolyMesh dodecahedron_hk72_ambo_dual_hk20() {
        return hankin(dual(ambo(hankin(dodecahedron(), 72.0f * D2R))), 20.0f * D2R);
    }

    // Base: truncatedIcosahedron, Ops: Tr(50), Ambo, Dual
    inline PolyMesh truncatedIcosahedron_truncate05_ambo_dual() {
        return dual(ambo(truncate(truncatedIcosahedron(), 50.0f * D2R)));
    }
  }
}
