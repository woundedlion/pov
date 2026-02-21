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
   * @brief Populates the target geometry buffer in local space.
   * @param out_geometry The buffer to populate.
   */
  virtual void generate(GeometryType &out_geometry) const = 0;

  virtual ~IGenerator() = default;
};

/**
 * @brief Specialized interface for PolyMesh generation
 */
struct IMeshGenerator : public IGenerator<PolyMesh> {
  virtual void generate(PolyMesh &out_mesh) const = 0;
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by ID.
 * Defaults to avoiding allocations by generating into the provided reference.
 */
struct SolidGenerator : public IMeshGenerator {
  int solid_id;

  SolidGenerator(int id) : solid_id(id) {}

  void generate(PolyMesh &out_mesh) const override {
    // Generate fresh mesh from solids registry
    out_mesh = Solids::get(solid_id);
  }
};

/**
 * @brief A generic solid generator wrapping the Solids:: registry by name.
 */
struct SolidNameGenerator : public IMeshGenerator {
  std::string solid_name;

  SolidNameGenerator(const std::string &name) : solid_name(name) {}

  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::get_by_name(solid_name);
  }
};

/**
 * @brief Helper generator for specific platonic solids
 */
struct IcosahedronGenerator : public IMeshGenerator {
  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::Platonic::icosahedron();
  }
};

struct DodecahedronGenerator : public IMeshGenerator {
  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::Platonic::dodecahedron();
  }
};

struct CubeGenerator : public IMeshGenerator {
  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::Platonic::cube();
  }
};

struct OctahedronGenerator : public IMeshGenerator {
  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::Platonic::octahedron();
  }
};

struct TetrahedronGenerator : public IMeshGenerator {
  void generate(PolyMesh &out_mesh) const override {
    out_mesh = Solids::Platonic::tetrahedron();
  }
};
