/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "math/geometry.h"
#include "mesh/mesh.h" // For MeshOps
#include "mesh/hankin.h"
#include "mesh/conway.h"
#include <cmath>
#include <string_view>
#include <span>

// --- Constants for Procedural Generation ---
/** Square root of 2. */
static constexpr float SQRT2 = 1.414213562373095f;
/** Tribonacci constant t, the real root of t^3 - t^2 - t - 1 = 0 (~1.83928676). */
static constexpr float TRIBONACCI_CONST = 1.839286755214161f;
/** Snub-cube truncation parameter. */
static constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
/** Truncated-dodecahedron/icosahedron truncation parameter. */
static constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + PHI);

namespace Solids {

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
  final_mesh.vertices.bind(geom, temp.vertices.size());
  final_mesh.vertices.append_bulk(temp.vertices.data(), temp.vertices.size());
  final_mesh.face_counts.bind(geom, temp.face_counts.size());
  final_mesh.face_counts.append_bulk(temp.face_counts.data(),
                                     temp.face_counts.size());
  final_mesh.faces.bind(geom, temp.faces.size());
  final_mesh.faces.append_bulk(temp.faces.data(), temp.faces.size());
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
 * @brief Materializes a compile-time static mesh into a runtime PolyMesh.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces arrays.
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

/**
 * @brief Fluent builder for chaining Conway operators with automatic arena
 * swapping.
 * @details Each method runs `mesh_ = op(mesh_, a_, b_)` then swaps a_/b_ (one
 * swap per call, PRIMITIVE polarity). Even-length composed ops
 * (gyro/needle/zip/bevel) return in `temp` (see COMPOSITION POLARITY in
 * conway.h), shifting the swap parity for one step; odd-length ones (meta)
 * behave like a primitive.
 */
class SolidBuilder {
  PolyMesh mesh_;     /**< Mesh being built; updated in place by each operator. */
  Arena *a_;          /**< Current output arena (swapped with b_ per op). */
  Arena *b_;          /**< Current scratch arena (swapped with a_ per op). */

public:
  /**
   * @brief Constructs a builder seeded with an initial mesh and arena pair.
   * @param seed Starting mesh, moved into the builder.
   * @param a Initial output arena.
   * @param b Initial scratch arena.
   */
  SolidBuilder(PolyMesh seed, Arena &a, Arena &b)
      : mesh_(std::move(seed)), a_(&a), b_(&b) {}

  /**
   * @brief Applies the dual operator (faces become vertices and vice versa).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &dual() {
    mesh_ = MeshOps::dual(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the kis operator (raise a pyramid on every face).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &kis() {
    mesh_ = MeshOps::kis(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the ambo operator (rectification: new vertex per edge).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &ambo() {
    mesh_ = MeshOps::ambo(mesh_, *a_, *b_);
    std::swap(a_, b_);
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
    mesh_ = MeshOps::truncate(mesh_, *a_, *b_, t);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the expand operator (cantellation: push faces outward).
   * @param t Expansion amount; default places square faces at the canonical gap.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &expand(float t = MeshOps::EXPAND_DEFAULT_T) {
    mesh_ = MeshOps::expand(mesh_, *a_, *b_, t);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the chamfer operator (replace edges with hexagons).
   * @param t Chamfer width as a fraction of the edge.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &chamfer(float t = 0.5f) {
    mesh_ = MeshOps::chamfer(mesh_, *a_, *b_, t);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the snub operator (chiral expansion with a twist).
   * @param t Expansion amount.
   * @param twist Rotation applied to each face, in radians.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &snub(float t = 0.5f, float twist = 0.0f) {
    mesh_ = MeshOps::snub(mesh_, *a_, *b_, t, twist);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the gyro operator (dual of snub; pentagonal chiral subdivision).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &gyro() {
    mesh_ = MeshOps::gyro(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Relaxes vertex positions toward a regular configuration.
   * @param iterations Number of smoothing iterations to run.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &relax(int iterations = 8) {
    mesh_ = MeshOps::relax(mesh_, *a_, *b_, iterations);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the meta operator (kis composed with join).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &meta() {
    mesh_ = MeshOps::meta(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the needle operator (kis of the dual; n = kd).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &needle() {
    mesh_ = MeshOps::needle(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the zip operator (dual of kis).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &zip() {
    mesh_ = MeshOps::zip(mesh_, *a_, *b_);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the bevel operator (truncate composed with ambo).
   * @param t Bevel depth along each edge.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &bevel(float t = 0.25f) {
    mesh_ = MeshOps::bevel(mesh_, *a_, *b_, t);
    std::swap(a_, b_);
    return *this;
  }
  /**
   * @brief Applies the Hankin star-pattern operator to each face.
   * @param angle Contact angle of the star pattern, in radians.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &hankin(float angle) {
    mesh_ = MeshOps::hankin(mesh_, *a_, *b_, angle);
    std::swap(a_, b_);
    return *this;
  }

  /**
   * @brief Finalizes the chain and yields the built mesh.
   * @return The accumulated PolyMesh, moved out of the builder.
   */
  PolyMesh build() { return std::move(mesh_); }
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
  return SolidBuilder(tetrahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
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
  return SolidBuilder(cube(a, b), a, b).truncate(1.0f / (2.0f + SQRT2)).build();
}
/**
 * @brief Builds a truncated octahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The truncated octahedron mesh.
 */
FLASHMEM static PolyMesh truncatedOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
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
      .bevel(1.0f / (2.0f + SQRT2))
      .relax(50)
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
      .snub(T_SNUB_CUBE, 0.28f)
      .relax(50)
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
  return SolidBuilder(icosahedron(a, b), a, b).truncate(1.0f / 3.0f).build();
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
      .relax(50)
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
      .bevel(1.0f / (2.0f + PHI))
      .relax(50)
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
      .relax(50)
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
 * @brief Builds a deltoidal hexecontahedron (dual of the rhombicosidodecahedron).
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
static constexpr float D2R = PI_F / 180.0f;

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
 * @brief Builds the cube_relax_bevel33_relax_hk675_expand5 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
cube_relax_bevel33_relax_hk675_expand5(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::cube(a, b), a, b)
      .relax(100)
      .bevel(0.33f)
      .relax(100)
      .hankin(67.5f * D2R)
      .expand(0.5f)
      .build();
}
/**
 * @brief Builds the icosahedron_bevel033_hk59 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_bevel033_hk59(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .bevel(0.33f)
      .hankin(59.0f * D2R)
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
 * @brief Builds the truncatedIcosidodecahedron_truncate50d_ambo_dual star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosidodecahedron_truncate50d_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosidodecahedron(a, b), a, b)
      .truncate(50.0f * D2R)
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
      .truncate(5.0f * D2R)
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
      .truncate(5.0f * D2R)
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
      .truncate(50.0f * D2R)
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
FLASHMEM static PolyMesh
icosahedron_snub_relax_truncate033_hankin62(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .snub()
      .relax()
      .truncate(0.33f)
      .hankin(62.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk35_ambo_hk62_ambo_relax_hk43 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
dodecahedron_hk35_ambo_hk62_ambo_relax_hk43(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .ambo()
      .relax(100)
      .hankin(43.0f * D2R)
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
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin59 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin59(Arena &a,
                                                            Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax()
      .truncate(0.01f)
      .hankin(59.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin73 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin73(Arena &a,
                                                            Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax()
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
      .relax(100)
      .hankin(77.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_bevel2_relax_gyro star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_bevel2_relax_gyro(Arena &a,
                                                                    Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .bevel(0.2f)
      .relax(100)
      .gyro()
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate33_hk64 star pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate33_hk64(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax(217)
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
FLASHMEM static PolyMesh
dodecahedron_ambo_bevel33_relax_hk66(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .ambo()
      .bevel(0.33f)
      .relax(100)
      .hankin(66.0f * D2R)
      .build();
}
} // namespace IslamicStarPatterns

/**
 * @brief Coarse generator-cost hint surfaced to the picker UI.
 * @details Display label only (Complex = Islamic star-pattern registry, Simple =
 * everything else); no runtime path gates on it.
 */
enum class Category { Simple, Complex };

/**
 * @brief One named solid in a registry.
 */
struct Entry {
  const char *name;                        /**< Registry key / display name. */
  PolyMesh (*generate)(Arena &a, Arena &b); /**< Generator building into arena pair (a, b). */
  Category category;                       /**< Simple or Complex cost class. */
};

/**
 * @brief Registry of the Platonic and the 13 Archimedean solids.
 * @details Order is load-bearing:
 * Collections::get_platonic/archimedean_solids() slice this array by fixed
 * offsets (Platonic 0-4, Archimedean 5-17).
 */
inline constexpr Entry simple_registry[] = {

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

/**
 * @brief Registry of Catalan solids (duals of the Archimedean solids).
 */
inline constexpr Entry catalan_registry[] = {
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

/**
 * @brief Registry of Islamic star-pattern solids.
 */
inline constexpr Entry islamic_registry[] = {
    {"cube_relax_bevel33_relax_hk675_expand5",
     IslamicStarPatterns::
         cube_relax_bevel33_relax_hk675_expand5,
     Category::Complex},
    {"truncatedIcosahedron_hk58_chamfer63",
     IslamicStarPatterns::truncatedIcosahedron_hk58_chamfer63,
     Category::Complex},
    {"dodecahedron_ambo_bevel33_relax_hk66",
     IslamicStarPatterns::dodecahedron_ambo_bevel33_relax_hk66,
     Category::Complex},
    {"truncatedIcosahedron_ambo_relax_truncate33_hk64",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_relax_truncate33_hk64,
     Category::Complex},
    {"dodecahedron_bevel2_relax_gyro",
     IslamicStarPatterns::dodecahedron_bevel2_relax_gyro,
     Category::Complex},
    {"truncatedIcosidodecahedron_bevel5_relax_hk77",
     IslamicStarPatterns::
         truncatedIcosidodecahedron_bevel5_relax_hk77,
     Category::Complex},
    {"truncatedOctahedron_gyro_kis_hk17",
     IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17, Category::Complex},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin59",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_relax_truncate001_hankin59,
     Category::Complex},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin73",
     IslamicStarPatterns::
         truncatedIcosahedron_ambo_relax_truncate001_hankin73,
     Category::Complex},
    {"icosahedron_ambo_truncate033_hankin59",
     IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59,
     Category::Complex},
    {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk43",
     IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_relax_hk43,
     Category::Complex},
    {"icosahedron_bevel033_hk59",
     IslamicStarPatterns::icosahedron_bevel033_hk59, Category::Complex},
    {"octahedron_hk17_ambo_hk73",
     IslamicStarPatterns::octahedron_hk17_ambo_hk73, Category::Complex},
    {"icosahedron_kis_gyro", IslamicStarPatterns::icosahedron_kis_gyro,
     Category::Complex},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosidodecahedron_truncate50d_ambo_dual,
     Category::Complex},
    {"icosidodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::icosidodecahedron_truncate5d_ambo_dual,
     Category::Complex},
    {"snubDodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::snubDodecahedron_truncate5d_ambo_dual,
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
    {"truncatedIcosahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosahedron_truncate50d_ambo_dual,
     Category::Complex},
    {"icosahedron_snub_relax_truncate033_hankin62",
     IslamicStarPatterns::icosahedron_snub_relax_truncate033_hankin62,
     Category::Complex}};

/** Total number of solids across all three registries. */
inline constexpr int NUM_ENTRIES =
    sizeof(simple_registry) / sizeof(simple_registry[0]) +
    sizeof(catalan_registry) / sizeof(catalan_registry[0]) +
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);

// simple_registry is [Platonic | Archimedean]; the static_assert checks the two
// counts exactly tile it so a boundary move can't silently mis-slice.
inline constexpr size_t PLATONIC_COUNT = 5;
inline constexpr size_t ARCHIMEDEAN_COUNT = 13;
static_assert(PLATONIC_COUNT + ARCHIMEDEAN_COUNT == std::size(simple_registry),
              "PLATONIC_COUNT + ARCHIMEDEAN_COUNT must equal simple_registry "
              "size; update the counts if the registry layout changes");

namespace Collections {
/**
 * @brief Returns the five Platonic solids.
 * @return Span over the Platonic entries (offset 0, count 5) of simple_registry.
 */
inline std::span<const Entry> get_platonic_solids() {
  return std::span<const Entry>(simple_registry, PLATONIC_COUNT);
}
/**
 * @brief Returns the 13 Archimedean solids.
 * @return Span over the Archimedean entries (offset 5, count 13) of
 * simple_registry.
 */
inline std::span<const Entry> get_archimedean_solids() {
  return std::span<const Entry>(simple_registry + PLATONIC_COUNT,
                                ARCHIMEDEAN_COUNT);
}
/**
 * @brief Returns all simple (Platonic and Archimedean) solids.
 * @return Span over the entire simple_registry.
 */
inline std::span<const Entry> get_simple_solids() {
  return std::span<const Entry>(simple_registry);
}
/**
 * @brief Returns all Catalan solids.
 * @return Span over the entire catalan_registry.
 */
inline std::span<const Entry> get_catalan_solids() {
  return std::span<const Entry>(catalan_registry);
}
/**
 * @brief Returns all Islamic star-pattern solids.
 * @return Span over the entire islamic_registry.
 */
inline std::span<const Entry> get_islamic_solids() {
  return std::span<const Entry>(islamic_registry);
}
} // namespace Collections

/**
 * @brief The three solid registries in flat global-index order.
 * @return Spans over the simple, then Catalan, then Islamic registries.
 * @details Single source of truth for the registry enumeration order; the
 * index- and name-based lookups below all derive from it.
 */
inline std::array<std::span<const Entry>, 3> all_registries() {
  return {std::span<const Entry>(simple_registry),
          std::span<const Entry>(catalan_registry),
          std::span<const Entry>(islamic_registry)};
}

/**
 * @brief Finds a registry entry by name across all registries.
 * @param name Candidate solid name.
 * @return Pointer to the matching entry, or nullptr when no name matches.
 */
inline const Entry *find_entry(std::string_view name) {
  for (std::span<const Entry> reg : all_registries())
    for (const Entry &entry : reg)
      if (name == entry.name)
        return &entry;
  return nullptr;
}

/**
 * @brief Looks up a registry entry by global index across all three registries.
 * @param index Zero-based index in [0, NUM_ENTRIES); traps if out of range.
 * @return Reference to the entry at that index.
 */
inline const Entry &get_entry(size_t index) {
  HS_CHECK(index < static_cast<size_t>(NUM_ENTRIES),
           "Solids::get_entry: index out of range");

  for (std::span<const Entry> reg : all_registries()) {
    if (index < reg.size())
      return reg[index];
    index -= reg.size();
  }
  __builtin_unreachable();
}

/**
 * @brief Builds the solid with the given name into the geometry arena.
 * @param geom Long-lived arena that backs the returned mesh.
 * @param a Output arena for even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @param name Registry name of the solid to build; traps if unknown.
 * @return The finalized solid mesh owned by geom.
 * @details For trusted (firmware) callers; an unknown name fails fast.
 */
FLASHMEM static PolyMesh get_by_name(Arena &geom, Arena &a, Arena &b,
                                     std::string_view name) {
  const Entry *entry = find_entry(name);
  HS_CHECK(entry, "Solids::get_by_name: unknown solid name");
  return finalize_solid(entry->generate(a, b), geom);
}

/**
 * @brief Tests whether a solid name exists in any registry, without trapping.
 * @param name Candidate registry name to look up.
 * @return True if the name matches a registered solid, false otherwise.
 * @details Trusted (firmware) callers use get_by_name() with its fail-fast
 * contract; untrusted boundaries (e.g. the WASM/JS mesh editor) must validate
 * with this first and reject unknown names rather than abort.
 */
inline bool has_name(std::string_view name) {
  return find_entry(name) != nullptr;
}

/**
 * @brief Builds a registry solid's unit vertex directions plus per-vertex
 *        orientation quaternions and nearest-neighbour gaps.
 * @param geom Arena backing the intermediate mesh (typically a caller scratch
 *        scope; nothing is retained after return).
 * @param scratch Scratch arena for odd pipeline stages.
 * @param name Registry name of the solid; traps if unknown.
 * @param max_points Capacity of the output arrays; traps if exceeded.
 * @param points Out: vertex directions projected onto the unit sphere (Catalan
 *        vertices sit at multiple radii).
 * @param quats Out: per-vertex Y-axis-to-direction rotations.
 * @param nn_angle Out: per-vertex nearest-neighbour angle (radians), for
 *        sizing per-vertex geometry to its local gap.
 * @return The vertex count written.
 * @details Shared by the volume-scatter effects (Raymarch). HS_COLD:
 * setup-only, keeps the build loops out of ITCM.
 */
HS_COLD static int build_vertex_directions(Arena &geom, Arena &scratch,
                                           std::string_view name,
                                           int max_points, Vector *points,
                                           Quaternion *quats,
                                           float *nn_angle) {
  PolyMesh mesh = get_by_name(geom, geom, scratch, name);
  int count = static_cast<int>(mesh.vertices.size());
  HS_CHECK(count <= max_points,
           "build_vertex_directions: vertex count exceeds capacity");
  for (int i = 0; i < count; ++i) {
    points[i] = mesh.vertices[i].normalized();
    quats[i] = make_rotation(Y_AXIS, points[i]);
  }
  // O(n^2) nearest-neighbor scan is intentional: cold setup path, vertex counts
  // are small, so the KD-tree's build overhead is not worth it here.
  for (int i = 0; i < count; ++i) {
    float max_dot = -1.0f;
    for (int j = 0; j < count; ++j)
      if (j != i)
        max_dot = std::max(max_dot, dot(points[i], points[j]));
    nn_angle[i] = acosf(hs::clamp(max_dot, -1.0f, 1.0f));
  }
  return count;
}

} // namespace Solids
