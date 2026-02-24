/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include "mesh.h"
#include "solids.h" // For access to Solids::* functions
#include <vector>

/**
 * @brief A generic concept/interface for generating geometry
 * @tparam GeometryType The target buffer type to populate (e.g., PolyMesh,
 * std::array)
 */
template <typename GeometryType> struct IGenerator {
  /**
   * @brief Generates and returns the target geometry buffer in local space.
   * @param geom The arena to allocate persistent geometry data.
   * @param scratch The arena to allocate temporary structures.
   * @return The generated geometry box.
   */
  virtual GeometryType generate(Arena &geom, Arena &scratch) const = 0;

  virtual ~IGenerator() = default;
};

/**
 * @brief Specialized interface for PolyMesh generation
 */
struct IMeshGenerator : public IGenerator<PolyMesh> {
  virtual PolyMesh generate(Arena &geom, Arena &scratch) const = 0;
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by ID.
 * Defaults to avoiding allocations by generating into the provided reference.
 */
struct SolidGenerator : public IMeshGenerator {
  int solid_id;

  SolidGenerator(int id) : solid_id(id) {}

  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    // Generate fresh mesh from solids registry
    return Solids::get(geom, scratch, scratch_arena_b, solid_id);
  }
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by name.
 */
struct SolidNameGenerator : public IMeshGenerator {
  std::string solid_name;

  SolidNameGenerator(const std::string &name) : solid_name(name) {}

  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    return Solids::get_by_name(geom, scratch, scratch_arena_b, solid_name);
  }
};

/**
 * @brief Helper generator for specific platonic solids
 */
struct IcosahedronGenerator : public IMeshGenerator {
  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    ScratchContext ctx(scratch, scratch_arena_b);
    return Solids::finalize_solid(Solids::Platonic::icosahedron(ctx), geom);
  }
};

struct DodecahedronGenerator : public IMeshGenerator {
  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    ScratchContext ctx(scratch, scratch_arena_b);
    return Solids::finalize_solid(Solids::Platonic::dodecahedron(ctx), geom);
  }
};

struct CubeGenerator : public IMeshGenerator {
  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    ScratchContext ctx(scratch, scratch_arena_b);
    return Solids::finalize_solid(Solids::Platonic::cube(ctx), geom);
  }
};

struct OctahedronGenerator : public IMeshGenerator {
  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    ScratchContext ctx(scratch, scratch_arena_b);
    return Solids::finalize_solid(Solids::Platonic::octahedron(ctx), geom);
  }
};

struct TetrahedronGenerator : public IMeshGenerator {
  PolyMesh generate(Arena &geom, Arena &scratch) const override {
    ScratchContext ctx(scratch, scratch_arena_b);
    return Solids::finalize_solid(Solids::Platonic::tetrahedron(ctx), geom);
  }
};
