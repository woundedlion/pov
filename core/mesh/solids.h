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
#include "mesh/relax_bakes_generated.h"
#include <cmath>
#include <string_view>
#include <span>

// --- Constants for Procedural Generation ---
/** Square root of 2. */
static constexpr float SQRT2 = 1.414213562373095f;
/** Tribonacci constant t, the real root of t^3 - t^2 - t - 1 = 0 (~1.83928676).
 */
static constexpr float TRIBONACCI_CONST = 1.839286755214161f;
/** Snub-cube truncation parameter. */
static constexpr float T_SNUB_CUBE = 1.0f / (1.0f + TRIBONACCI_CONST);
/** Snub-cube twist. */
static constexpr float SNUB_CUBE_TWIST = 0.28f;
/** Truncated-dodecahedron/icosahedron truncation parameter. */
static constexpr float T_TRUNC_ICOS = 1.0f / (2.0f + PHI);
/** Truncated-cube/cuboctahedron truncation parameter. */
static constexpr float T_TRUNC_CUBE = 1.0f / (2.0f + SQRT2);
/** Truncated tetra/octa/icosahedron truncation parameter. */
static constexpr float T_TRUNC_THIRD = 1.0f / 3.0f;

namespace Solids {

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
 * @return True iff the tables form a closed unit-sphere polyhedron.
 * @details Checks that face_counts spans the flat face list exactly, every face
 * index addresses a listed vertex, Euler's formula holds for E = sum/2, and
 * every vertex is unit length. Compares squared lengths so the whole check is
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

/**
 * @brief Fluent builder for chaining Conway operators with automatic arena
 * swapping.
 * @details Each method runs `mesh = op(mesh, output_arena, scratch_arena)` then
 * swaps the two arenas; every operator returns its output in the arena passed
 * as `target` (COMPOSITION POLARITY in conway.h), so after each step the mesh
 * sits in `scratch_arena` and the next step writes into the other arena.
 *
 * Callers build the seed into `a`, so the first operator reads its input from
 * the arena it writes its output into. That costs peak arena, not correctness:
 * every operator binds its output before opening a ScratchScope, and a bump
 * arena never rewinds below a live allocation.
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
   * steps: EXTRACT dumps the resulting bits (generation), VERIFY asserts they
   * match the committed payload (the native re-derivation test).
   * @param bake Payload whose guarded topology must match the current mesh.
   * @return Reference to this builder for chaining.
   */
  HS_COLD_MEMBER SolidBuilder &relax_baked(const MeshOps::RelaxBake &bake) {
#if defined(HS_RELAX_BAKE_EXTRACT) || defined(HS_RELAX_BAKE_VERIFY)
    HS_CHECK(MeshOps::relax_topology_hash(mesh) == bake.topology_hash,
             "relax bake: source topology differs");
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
#else // HS_RELAX_BAKE_EXTRACT: emit the payload for the generated header.
    hs::log("RELAX_BAKE_BEGIN %s %d %lu %lu %lu %08lx %08lx", bake.name,
            static_cast<int>(bake.iterations),
            static_cast<unsigned long>(mesh.vertices.size()),
            static_cast<unsigned long>(mesh.face_counts.size()),
            static_cast<unsigned long>(mesh.faces.size()),
            static_cast<unsigned long>(bake.topology_hash),
            static_cast<unsigned long>(output_hash));
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
  return SolidBuilder(cube(a, b), a, b).bevel(T_TRUNC_CUBE).relax(50).build();
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
  return SolidBuilder(icosahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a rhombicosidodecahedron.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rhombicosidodecahedron mesh.
 */
FLASHMEM static PolyMesh rhombicosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).expand().relax(50).build();
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
static constexpr float D2R = PI_F / 180.0f;

/** Truncation depth of the `*_truncate5d_*` recipes: a dimensionless edge
 * fraction short of the ambo pinch at t = 0.5. The `5d` in those names is the
 * degree slot of the parameter sweep that found the value, not a unit of it. */
static constexpr float TRUNCATE_T_NEAR = 0.0872664601f;
/** Truncation depth of the `*_truncate50d_*` recipes: a dimensionless edge
 * fraction past the ambo pinch, where the cut faces self-intersect by design.
 * The `50d` in those names is a sweep slot, not a unit of the value. */
static constexpr float TRUNCATE_T_FAR = 0.87266463f;

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
 * @brief Builds the truncatedIcosahedron_ambo_relax100_hk54_needle star
 * pattern.
 * @param a Output arena for the result and even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax100_hk54_needle(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .hankin(54.0f * D2R)
      .needle()
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

/**
 * @brief Coarse generator-cost hint surfaced to the picker UI.
 * @details Display label only (Complex = Islamic star-pattern registry, Simple
 * = everything else); no runtime path gates on it.
 */
enum class Category { Simple, Complex };

/**
 * @brief Authored Conway operator, including the composite ops.
 * @details expand_to_primitives() (mesh/recipe.h) lowers composites to the
 * primitives plus AMBO; authored recipes keep the composite.
 */
enum class Op : uint8_t {
  TRUNCATE,
  EXPAND,
  SNUB,
  CHAMFER,
  HANKIN,
  RELAX,
  KIS,
  DUAL,
  AMBO,
  BEVEL,
  GYRO,
  META,
  NEEDLE,
  ZIP
};

/**
 * @brief One authored step in a recipe's op chain.
 */
struct OpStep {
  Op op; /**< Operator applied at this step. */
  /**
   * @brief t / contact angle (radians) / RELAX iterations.
   * @details Unread on a RELAX step carrying a `bake`; such steps leave it at
   * zero rather than naming a count the replay never runs.
   */
  float param = 0.0f;
  float twist = 0.0f; /**< SNUB face rotation, radians. */
  /**
   * @brief RELAX bake this step lands on, mirroring the generator's
   * relax_baked() call; null replays `param` live iterations.
   * @details When set, every replay of the step (bitwise gate, on-screen build
   * leg, and eager clean endpoint) resolves through the baked converged mesh
   * instead of `param` smoothing steps, so the recipe lands on the shipped
   * geometry rather than a mid-convergence freeze. Left null for a mid-chain
   * relax on a mesh with no bake, which keeps iterating.
   */
  const MeshOps::RelaxBake *bake = nullptr;
};

/**
 * @brief Declarative op chain mirroring a registry generator.
 * @details Parallel declaration of a generator's chain, proven bitwise-equal
 * to it in tests/test_solids.h. The generators stay the source of truth for
 * shipping geometry.
 */
struct Recipe {
  uint8_t seed;        /**< simple_registry index of the base solid. */
  const OpStep *steps; /**< Authored op chain, applied left to right. */
  uint8_t count;       /**< Number of steps. */
};

/**
 * @brief One named solid in a registry.
 */
struct Entry {
  const char *name; /**< Registry key / display name. */
  PolyMesh (*generate)(
      Arena &a, Arena &b); /**< Generator building into arena pair (a, b). */
  Category category;       /**< Simple or Complex cost class. */
  const Recipe *recipe = nullptr; /**< Declarative chain mirror; null = none. */
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

// Recipe seed indices into simple_registry, pinned to the registry order so a
// reorder fails to compile.
inline constexpr uint8_t SEED_OCTAHEDRON = 2;
inline constexpr uint8_t SEED_DODECAHEDRON = 3;
inline constexpr uint8_t SEED_ICOSAHEDRON = 4;
inline constexpr uint8_t SEED_TRUNCATED_OCTAHEDRON = 8;
inline constexpr uint8_t SEED_RHOMBICUBOCTAHEDRON = 9;
inline constexpr uint8_t SEED_ICOSIDODECAHEDRON = 12;
inline constexpr uint8_t SEED_TRUNCATED_ICOSAHEDRON = 14;
inline constexpr uint8_t SEED_TRUNCATED_ICOSIDODECAHEDRON = 16;
inline constexpr uint8_t SEED_SNUB_DODECAHEDRON = 17;

static_assert(std::string_view(simple_registry[SEED_OCTAHEDRON].name) ==
              "octahedron");
static_assert(std::string_view(simple_registry[SEED_DODECAHEDRON].name) ==
              "dodecahedron");
static_assert(
    std::string_view(simple_registry[SEED_RHOMBICUBOCTAHEDRON].name) ==
    "rhombicuboctahedron");
static_assert(std::string_view(simple_registry[SEED_ICOSIDODECAHEDRON].name) ==
              "icosidodecahedron");
static_assert(std::string_view(simple_registry[SEED_SNUB_DODECAHEDRON].name) ==
              "snubDodecahedron");
static_assert(std::string_view(simple_registry[SEED_ICOSAHEDRON].name) ==
              "icosahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_OCTAHEDRON].name) ==
    "truncatedOctahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_ICOSAHEDRON].name) ==
    "truncatedIcosahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_ICOSIDODECAHEDRON].name) ==
    "truncatedIcosidodecahedron");

/** Step table for dodecahedron_hk62_ambo_hk62. */
inline constexpr OpStep DODECAHEDRON_HK62_AMBO_HK62_STEPS[] = {
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk62_ambo_hk62. */
inline constexpr Recipe DODECAHEDRON_HK62_AMBO_HK62_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK62_AMBO_HK62_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK62_AMBO_HK62_STEPS))};

/** Step table for octahedron_hk17_ambo_hk73. */
inline constexpr OpStep OCTAHEDRON_HK17_AMBO_HK73_STEPS[] = {
    {Op::HANKIN, 17.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 73.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::octahedron_hk17_ambo_hk73. */
inline constexpr Recipe OCTAHEDRON_HK17_AMBO_HK73_RECIPE = {
    SEED_OCTAHEDRON, OCTAHEDRON_HK17_AMBO_HK73_STEPS,
    static_cast<uint8_t>(std::size(OCTAHEDRON_HK17_AMBO_HK73_STEPS))};

/** Step table for octahedron_hk34_ambo_hk72. */
inline constexpr OpStep OCTAHEDRON_HK34_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 34.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::octahedron_hk34_ambo_hk72. */
inline constexpr Recipe OCTAHEDRON_HK34_AMBO_HK72_RECIPE = {
    SEED_OCTAHEDRON, OCTAHEDRON_HK34_AMBO_HK72_STEPS,
    static_cast<uint8_t>(std::size(OCTAHEDRON_HK34_AMBO_HK72_STEPS))};

/** Step table for dodecahedron_hk54_ambo_hk72. */
inline constexpr OpStep DODECAHEDRON_HK54_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 54.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk54_ambo_hk72. */
inline constexpr Recipe DODECAHEDRON_HK54_AMBO_HK72_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK54_AMBO_HK72_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK54_AMBO_HK72_STEPS))};

/** Step table for icosahedron_ambo_truncate033_hankin59. */
inline constexpr OpStep ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS[] = {
    {Op::AMBO},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 59.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59.
 */
inline constexpr Recipe ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS))};

/** Step table for rhombicuboctahedron_hk63_ambo_hk63. */
inline constexpr OpStep RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS[] = {
    {Op::HANKIN, 63.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 63.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63. */
inline constexpr Recipe RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE = {
    SEED_RHOMBICUBOCTAHEDRON, RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS,
    static_cast<uint8_t>(std::size(RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS))};

/** Step table for icosahedron_kis_gyro. */
inline constexpr OpStep ICOSAHEDRON_KIS_GYRO_STEPS[] = {{Op::KIS}, {Op::GYRO}};
/** Recipe mirror of IslamicStarPatterns::icosahedron_kis_gyro. */
inline constexpr Recipe ICOSAHEDRON_KIS_GYRO_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_KIS_GYRO_STEPS,
    static_cast<uint8_t>(std::size(ICOSAHEDRON_KIS_GYRO_STEPS))};

/** Step table for dodecahedron_hk72_ambo_dual_hk20. */
inline constexpr OpStep DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS[] = {
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::DUAL},
    {Op::HANKIN, 20.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20. */
inline constexpr Recipe DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS))};

/** Step table for icosidodecahedron_truncate5d_ambo_dual. */
inline constexpr OpStep ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_NEAR},
    {Op::AMBO},
    {Op::DUAL}};
/** Recipe mirror of IslamicStarPatterns::icosidodecahedron_truncate5d_ambo_dual. */
inline constexpr Recipe ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE = {
    SEED_ICOSIDODECAHEDRON, ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS))};

/** Step table for snubDodecahedron_truncate5d_ambo_dual. */
inline constexpr OpStep SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_NEAR},
    {Op::AMBO},
    {Op::DUAL}};
/** Recipe mirror of IslamicStarPatterns::snubDodecahedron_truncate5d_ambo_dual. */
inline constexpr Recipe SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE = {
    SEED_SNUB_DODECAHEDRON, SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS))};

/** Step table for dodecahedron_ambo_bevel33_relax_hk66. */
inline constexpr OpStep DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS[] = {
    {Op::AMBO},
    {Op::BEVEL, 0.33f},
    {.op = Op::RELAX, .bake = &RelaxBakes::dodecahedron_ambo_bevel33_converged},
    {Op::HANKIN, 66.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_ambo_bevel33_relax_hk66. */
inline constexpr Recipe DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS,
    static_cast<uint8_t>(
        std::size(DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax100_hk54_needle. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::HANKIN, 54.0f * IslamicStarPatterns::D2R},
        {Op::NEEDLE}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax100_hk54_needle.
 */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_RECIPE =
    {SEED_TRUNCATED_ICOSAHEDRON,
     TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_STEPS,
     static_cast<uint8_t>(
         std::size(TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_STEPS))};

/** Step table for truncatedIcosahedron_hk58_chamfer63. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS[] = {
    {Op::HANKIN, 58.0f * IslamicStarPatterns::D2R}, {Op::CHAMFER, 0.63f}};
/** Recipe mirror of IslamicStarPatterns::truncatedIcosahedron_hk58_chamfer63. */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON, TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate33_hk64. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.33f},
        {Op::HANKIN, 64.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate33_hk64.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS,
        static_cast<uint8_t>(
            std::size(TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS))};

/** Step table for dodecahedron_bevel2_relax_gyro. */
inline constexpr OpStep DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS[] = {
    {Op::BEVEL, 0.2f},
    {.op = Op::RELAX, .bake = &RelaxBakes::dodecahedron_bevel20_converged},
    {Op::GYRO}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_bevel2_relax_gyro. */
inline constexpr Recipe DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS))};

/** Step table for truncatedIcosidodecahedron_bevel5_relax_hk77. */
inline constexpr OpStep TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS[] =
    {{Op::BEVEL, 0.5f},
     {.op = Op::RELAX,
      .bake = &RelaxBakes::truncated_icosidodecahedron_bevel50_relax100},
     {Op::HANKIN, 77.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosidodecahedron_bevel5_relax_hk77.
 */
inline constexpr Recipe TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_RECIPE = {
    SEED_TRUNCATED_ICOSIDODECAHEDRON,
    TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS))};

/** Step table for truncatedOctahedron_gyro_kis_hk17. */
inline constexpr OpStep TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS[] = {
    {Op::GYRO}, {Op::KIS}, {Op::HANKIN, 17.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17. */
inline constexpr Recipe TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE = {
    SEED_TRUNCATED_OCTAHEDRON, TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS,
    static_cast<uint8_t>(std::size(TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate001_hankin59. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.01f},
        {Op::HANKIN, 59.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin59.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate001_hankin73. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.01f},
        {Op::HANKIN, 73.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin73.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS))};

/** Step table for dodecahedron_hk35_ambo_hk62_ambo_relax_hk42. */
inline constexpr OpStep DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS[] = {
    {Op::HANKIN, 35.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {.op = Op::RELAX,
     .bake = &RelaxBakes::dodecahedron_hankin_ambo_hankin_ambo_converged},
    {Op::HANKIN, 42.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_relax_hk42.
 */
inline constexpr Recipe DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS,
    static_cast<uint8_t>(
        std::size(DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS))};

/** Step table for truncatedIcosidodecahedron_truncate50d_ambo_dual. */
inline constexpr OpStep
    TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
        {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_FAR},
        {Op::AMBO},
        {Op::DUAL}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosidodecahedron_truncate50d_ambo_dual.
 */
inline constexpr Recipe
    TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE = {
        SEED_TRUNCATED_ICOSIDODECAHEDRON,
        TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS))};

/** Step table for truncatedIcosahedron_hk54_ambo_hk72. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 54.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72. */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON, TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS))};

/** Step table for truncatedIcosahedron_truncate50d_ambo_dual. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_FAR},
    {Op::AMBO},
    {Op::DUAL}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_truncate50d_ambo_dual.
 */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON,
    TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS))};

/** Step table for icosahedron_snub_relax_truncate033_hankin62. */
inline constexpr OpStep ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS[] = {
    {Op::SNUB, 0.5f, 0.0f},
    {.op = Op::RELAX, .bake = &RelaxBakes::icosahedron_snub_converged},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::icosahedron_snub_relax_truncate033_hankin62.
 */
inline constexpr Recipe ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS))};

/**
 * @brief Registry of Islamic star-pattern solids.
 */
inline constexpr Entry islamic_registry[] = {
    {"truncatedIcosahedron_ambo_relax100_hk54_needle",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax100_hk54_needle,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_RECIPE},
    {"dodecahedron_hk62_ambo_hk62",
     IslamicStarPatterns::dodecahedron_hk62_ambo_hk62, Category::Complex,
     &DODECAHEDRON_HK62_AMBO_HK62_RECIPE},
    {"truncatedIcosahedron_hk58_chamfer63",
     IslamicStarPatterns::truncatedIcosahedron_hk58_chamfer63,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_RECIPE},
    {"dodecahedron_ambo_bevel33_relax_hk66",
     IslamicStarPatterns::dodecahedron_ambo_bevel33_relax_hk66,
     Category::Complex, &DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate33_hk64",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate33_hk64,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE},
    {"dodecahedron_bevel2_relax_gyro",
     IslamicStarPatterns::dodecahedron_bevel2_relax_gyro, Category::Complex,
     &DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE},
    {"truncatedIcosidodecahedron_bevel5_relax_hk77",
     IslamicStarPatterns::truncatedIcosidodecahedron_bevel5_relax_hk77,
     Category::Complex, &TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_RECIPE},
    {"truncatedOctahedron_gyro_kis_hk17",
     IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17, Category::Complex,
     &TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin59",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin59,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin73",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin73,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_RECIPE},
    {"icosahedron_ambo_truncate033_hankin59",
     IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59,
     Category::Complex, &ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_RECIPE},
    {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42",
     IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_relax_hk42,
     Category::Complex, &DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE},
    {"octahedron_hk17_ambo_hk73",
     IslamicStarPatterns::octahedron_hk17_ambo_hk73, Category::Complex,
     &OCTAHEDRON_HK17_AMBO_HK73_RECIPE},
    {"icosahedron_kis_gyro", IslamicStarPatterns::icosahedron_kis_gyro,
     Category::Complex, &ICOSAHEDRON_KIS_GYRO_RECIPE},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosidodecahedron_truncate50d_ambo_dual,
     Category::Complex,
     &TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE},
    {"icosidodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::icosidodecahedron_truncate5d_ambo_dual,
     Category::Complex, &ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE},
    {"snubDodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::snubDodecahedron_truncate5d_ambo_dual,
     Category::Complex, &SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE},
    {"octahedron_hk34_ambo_hk72",
     IslamicStarPatterns::octahedron_hk34_ambo_hk72, Category::Complex,
     &OCTAHEDRON_HK34_AMBO_HK72_RECIPE},
    {"rhombicuboctahedron_hk63_ambo_hk63",
     IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63, Category::Complex,
     &RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE},
    {"truncatedIcosahedron_hk54_ambo_hk72",
     IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE},
    {"dodecahedron_hk54_ambo_hk72",
     IslamicStarPatterns::dodecahedron_hk54_ambo_hk72, Category::Complex,
     &DODECAHEDRON_HK54_AMBO_HK72_RECIPE},
    {"dodecahedron_hk72_ambo_dual_hk20",
     IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20, Category::Complex,
     &DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE},
    {"truncatedIcosahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosahedron_truncate50d_ambo_dual,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE},
    {"icosahedron_snub_relax_truncate033_hankin62",
     IslamicStarPatterns::icosahedron_snub_relax_truncate033_hankin62,
     Category::Complex, &ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE}};

/** Total number of solids across all three registries. */
inline constexpr int NUM_ENTRIES =
    sizeof(simple_registry) / sizeof(simple_registry[0]) +
    sizeof(catalan_registry) / sizeof(catalan_registry[0]) +
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);

// simple_registry is [Platonic | Archimedean]; the static_asserts check the two
// counts exactly tile it and name the entries either side of the slice
// boundary, so a boundary move can't silently mis-slice.
inline constexpr size_t PLATONIC_COUNT = 5;
inline constexpr size_t ARCHIMEDEAN_COUNT = 13;
static_assert(PLATONIC_COUNT + ARCHIMEDEAN_COUNT == std::size(simple_registry),
              "PLATONIC_COUNT + ARCHIMEDEAN_COUNT must equal simple_registry "
              "size; update the counts if the registry layout changes");
static_assert(std::string_view(simple_registry[PLATONIC_COUNT - 1].name) ==
                  "icosahedron",
              "PLATONIC_COUNT must end the Platonic run");
static_assert(std::string_view(simple_registry[PLATONIC_COUNT].name) ==
                  "truncatedTetrahedron",
              "PLATONIC_COUNT must start the Archimedean run");

inline constexpr size_t CATALAN_COUNT = 13;
inline constexpr size_t ISLAMIC_COUNT = 24;
static_assert(CATALAN_COUNT == std::size(catalan_registry),
              "catalan_registry size changed; update CATALAN_COUNT and the "
              "README registry table");
static_assert(ISLAMIC_COUNT == std::size(islamic_registry),
              "islamic_registry size changed; update ISLAMIC_COUNT and the "
              "README registry table");
static_assert(PLATONIC_COUNT + ARCHIMEDEAN_COUNT + CATALAN_COUNT +
                      ISLAMIC_COUNT ==
                  static_cast<size_t>(NUM_ENTRIES),
              "NUM_ENTRIES must equal the sum of the per-registry counts");

namespace Collections {
/**
 * @brief Returns the five Platonic solids.
 * @return Span over the Platonic entries (offset 0, count 5) of
 * simple_registry.
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
inline constexpr std::array<std::span<const Entry>, 3> all_registries() {
  return {std::span<const Entry>(simple_registry),
          std::span<const Entry>(catalan_registry),
          std::span<const Entry>(islamic_registry)};
}

/**
 * @brief Tests whether every registered recipe's seed indexes simple_registry.
 * @return True when no entry carries an out-of-range Recipe::seed.
 */
inline constexpr bool all_recipe_seeds_in_range() {
  for (std::span<const Entry> reg : all_registries())
    for (const Entry &entry : reg)
      if (entry.recipe && entry.recipe->seed >= std::size(simple_registry))
        return false;
  return true;
}
static_assert(all_recipe_seeds_in_range(),
              "Recipe::seed must index simple_registry; the replay sites index "
              "it without a runtime bound");

/**
 * @brief Tests whether every entry's Category matches the registry it sits in.
 * @return True when islamic_registry is uniformly Complex and the other two are
 * uniformly Simple.
 */
inline constexpr bool all_categories_match_registry() {
  for (const Entry &entry : simple_registry)
    if (entry.category != Category::Simple)
      return false;
  for (const Entry &entry : catalan_registry)
    if (entry.category != Category::Simple)
      return false;
  for (const Entry &entry : islamic_registry)
    if (entry.category != Category::Complex)
      return false;
  return true;
}
static_assert(all_categories_match_registry(),
              "Category must be Complex on islamic_registry and Simple "
              "everywhere else");

/**
 * @brief Tests whether every Islamic star-pattern entry carries a Recipe mirror.
 * @return True when no islamic_registry entry has a null recipe.
 */
inline constexpr bool all_islamic_entries_have_recipes() {
  for (const Entry &entry : islamic_registry)
    if (!entry.recipe)
      return false;
  return true;
}
static_assert(all_islamic_entries_have_recipes(),
              "every islamic_registry entry needs its Recipe mirror; the "
              "lowering and morph-feasibility sweeps skip recipe-less entries");

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
[[maybe_unused]] FLASHMEM static PolyMesh
get_by_name(Arena &geom, Arena &a, Arena &b, std::string_view name) {
  const Entry *entry = find_entry(name);
  HS_CHECK(entry, "Solids::get_by_name: unknown solid name");
  return finalize_solid(entry->generate(a, b), geom);
}

/**
 * @brief Builds a registry solid's unit vertex directions plus per-vertex
 *        orientation quaternions and nearest-neighbour gaps.
 * @param scratch Arena backing the intermediate mesh; nothing is retained
 *        after return.
 * @param temp Alternate scratch arena for odd pipeline stages.
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
[[maybe_unused]] HS_COLD static int
build_vertex_directions(Arena &scratch, Arena &temp, std::string_view name,
                        int max_points, Vector *points, Quaternion *quats,
                        float *nn_angle) {
  const Entry *entry = find_entry(name);
  HS_CHECK(entry, "build_vertex_directions: unknown solid name");
  // Read straight out of the generator's arena pair; nothing outlives the call,
  // so finalize_solid's long-lived copy would buy nothing.
  PolyMesh mesh = entry->generate(scratch, temp);
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

HS_O3_END

} // namespace Solids
