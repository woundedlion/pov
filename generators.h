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
   * @return The generated geometry box.
   */
  virtual GeometryType generate() const = 0;

  virtual ~IGenerator() = default;
};

/**
 * @brief Specialized interface for PolyMesh generation
 */
struct IMeshGenerator : public IGenerator<PolyMesh> {
  virtual PolyMesh generate() const = 0;
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by ID.
 * Defaults to avoiding allocations by generating into the provided reference.
 */
struct SolidGenerator : public IMeshGenerator {
  int solid_id;

  SolidGenerator(int id) : solid_id(id) {}

  PolyMesh generate() const override {
    // Generate fresh mesh from solids registry
    return Solids::get(solid_id);
  }
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by name.
 */
struct SolidNameGenerator : public IMeshGenerator {
  std::string solid_name;

  SolidNameGenerator(const std::string &name) : solid_name(name) {}

  PolyMesh generate() const override { return Solids::get_by_name(solid_name); }
};

/**
 * @brief Helper generator for specific platonic solids
 */
struct IcosahedronGenerator : public IMeshGenerator {
  PolyMesh generate() const override { return Solids::Platonic::icosahedron(); }
};

struct DodecahedronGenerator : public IMeshGenerator {
  PolyMesh generate() const override {
    return Solids::Platonic::dodecahedron();
  }
};

struct CubeGenerator : public IMeshGenerator {
  PolyMesh generate() const override { return Solids::Platonic::cube(); }
};

struct OctahedronGenerator : public IMeshGenerator {
  PolyMesh generate() const override { return Solids::Platonic::octahedron(); }
};

struct TetrahedronGenerator : public IMeshGenerator {
  PolyMesh generate() const override { return Solids::Platonic::tetrahedron(); }
};
