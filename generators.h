/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include "mesh.h"
#include "solids.h" // For access to Solids::* functions

/**
 * @brief A generic concept/interface for generating geometry
 * @tparam GeometryType The target buffer type to populate (e.g., PolyMesh,
 * std::array)
 */
template <typename GeometryType> struct IGenerator {
  /**
   * @brief Generates and returns the target geometry buffer in local space.
   * @param geom The arena to allocate persistent geometry data.
   * @param a Scratch arena A.
   * @param b Scratch arena B.
   * @return The generated geometry box.
   */
  virtual GeometryType generate(Arena &geom, Arena &a, Arena &b) const = 0;

  virtual ~IGenerator() = default;
};

/**
 * @brief Specialized interface for PolyMesh generation
 */
struct IMeshGenerator : public IGenerator<PolyMesh> {
  virtual PolyMesh generate(Arena &geom, Arena &a, Arena &b) const = 0;
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by ID.
 * Defaults to avoiding allocations by generating into the provided reference.
 */
struct SolidGenerator : public IMeshGenerator {
  int solid_id;

  SolidGenerator(int id) : solid_id(id) {}

  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::get(geom, a, b, solid_id);
  }
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by name.
 */
struct SolidNameGenerator : public IMeshGenerator {
  std::string_view solid_name;

  SolidNameGenerator(std::string_view name) : solid_name(name) {}

  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::get_by_name(geom, a, b, solid_name);
  }
};

/**
 * @brief Helper generator for specific platonic solids
 */
struct IcosahedronGenerator : public IMeshGenerator {
  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::finalize_solid(Solids::Platonic::icosahedron(a, b), geom);
  }
};

struct DodecahedronGenerator : public IMeshGenerator {
  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), geom);
  }
};

struct CubeGenerator : public IMeshGenerator {
  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::finalize_solid(Solids::Platonic::cube(a, b), geom);
  }
};

struct OctahedronGenerator : public IMeshGenerator {
  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::finalize_solid(Solids::Platonic::octahedron(a, b), geom);
  }
};

struct TetrahedronGenerator : public IMeshGenerator {
  FLASHMEM PolyMesh generate(Arena &geom, Arena &a, Arena &b) const override {
    ScopedScratch _(b);
    return Solids::finalize_solid(Solids::Platonic::tetrahedron(a, b), geom);
  }
};

/**
 * @brief Generates a mesh using scratch arenas.
 * Encapsulates ScopedScratch boilerplate for any generator
 * that implements generate(Arena&, Arena&, Arena&).
 *
 * Usage: mesh = generate_mesh<DodecahedronGenerator>(persistent_arena);
 *        mesh = generate_mesh<SolidGenerator>(persistent_arena, solid_id);
 */
template <typename Gen, typename... Args>
inline auto generate_mesh(Arena &geom, Args &&...args) {
  Gen gen(std::forward<Args>(args)...);
  scratch_arena_a.reset();
  scratch_arena_b.reset();
  ScopedScratch _(scratch_arena_a);
  return gen.generate(geom, scratch_arena_a, scratch_arena_b);
}
