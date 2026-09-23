/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file relax_bake.h
 * @brief Relax bake payloads and platform-independent mesh identity checks.
 */
#include "mesh/mesh.h"

namespace MeshOps {
HS_O3_BEGIN

/**
 * @brief Host-generated exact vertex payload for one deterministic relax input.
 * @details `iterations` is the relax count this payload was baked at
 * (early-stop on convergence applies, so any count past convergence yields the
 * same converged mesh; a count short of convergence deliberately freezes a
 * pre-converged configuration). The same bits are loaded on host and device, so
 * the two platforms render bit-identically. `source_hash` and `topology_hash`
 * guard that the live source mesh's quantized vertices and connectivity still
 * match what was baked against.
 */
struct RelaxBake {
  const char *name;
  const uint32_t *vertex_bits;
  uint16_t vertex_count;
  uint16_t face_count;
  uint16_t index_count;
  uint16_t iterations;
  uint32_t source_hash;
  uint32_t topology_hash;
  uint32_t output_hash;
};

/**
 * @brief Fixed-point scale used to identify relax source coordinates.
 * @details Together with RELAX_SOURCE_BIAS, the committed bake sources remain
 * at least RELAX_SOURCE_MIN_MARGIN from every quantization boundary while
 * topology-and-dimension peers retain distinct identities.
 */
inline constexpr float RELAX_SOURCE_SCALE = 2013.0f;

/** @brief Offset placing baked sources away from quantization boundaries. */
inline constexpr float RELAX_SOURCE_BIAS = 0.7361977398f;

/** @brief Required distance between a baked source and the nearest boundary. */
inline constexpr float RELAX_SOURCE_MIN_MARGIN = 1.0e-5f;

/** @brief Quantizes one relax source coordinate onto its identity grid. */
inline int32_t relax_source_coordinate(float coordinate) {
  const float magnitude =
      (coordinate < 0.0f ? -coordinate : coordinate) * RELAX_SOURCE_SCALE;
  const int32_t quantized = static_cast<int32_t>(magnitude + RELAX_SOURCE_BIAS);
  return coordinate < 0.0f ? -quantized : quantized;
}

/** @brief FNV-1a 32-bit offset basis, seeding every relax hash. */
inline constexpr uint32_t FNV1A_BASIS = 2166136261u;

/**
 * @brief One FNV-1a 32-bit round.
 * @param hash Accumulator, seeded from FNV1A_BASIS.
 * @param word Word mixed into the accumulator.
 * @return The updated accumulator.
 * @details Shared by the relax topology hash, the baked-payload load check and
 *   the bake extract/verify tooling, which must agree bit-for-bit.
 */
inline uint32_t fnv1a_step(uint32_t hash, uint32_t word) {
  return (hash ^ word) * 16777619u;
}

/** @brief Hashes platform-independent relax topology and dimensions. */
inline uint32_t relax_topology_hash(const PolyMesh &mesh) {
  uint32_t hash = FNV1A_BASIS;
  auto mix = [&](uint32_t word) { hash = fnv1a_step(hash, word); };
  mix(static_cast<uint32_t>(mesh.vertices.size()));
  mix(static_cast<uint32_t>(mesh.get_face_counts_size()));
  mix(static_cast<uint32_t>(mesh.get_faces_size()));
  for (size_t i = 0; i < mesh.get_face_counts_size(); ++i)
    mix(mesh.get_face_counts_data()[i]);
  for (size_t i = 0; i < mesh.get_faces_size(); ++i)
    mix(mesh.get_faces_data()[i]);
  return hash;
}

/**
 * @brief Hashes fixed-point-rounded vertices identifying a relax source mesh.
 * @details A fixed-point grid absorbs endpoint reconstruction roundoff while
 * retaining vertex order and parameterized geometry in the identity. The grid
 * phase keeps generated sources clear of rounding boundaries.
 */
inline uint32_t relax_source_hash(const PolyMesh &mesh) {
  uint32_t hash = FNV1A_BASIS;
  auto mix = [&](float coordinate) {
    hash = fnv1a_step(
        hash, static_cast<uint32_t>(relax_source_coordinate(coordinate)));
  };
  for (const math::Vector &v : mesh.vertices) {
    mix(v.x);
    mix(v.y);
    mix(v.z);
  }
  return hash;
}

/**
 * @brief Finds the nearest source-identity quantization boundary.
 * @return Coordinate distance to the nearest boundary across the mesh.
 */
inline float relax_source_quantization_margin(const PolyMesh &mesh) {
  float minimum = 1.0f;
  auto measure = [&](float coordinate) {
    const float magnitude =
        (coordinate < 0.0f ? -coordinate : coordinate) * RELAX_SOURCE_SCALE;
    const float biased = magnitude + RELAX_SOURCE_BIAS;
    const float fraction = biased - static_cast<int32_t>(biased);
    const float grid_distance =
        (fraction < 1.0f - fraction ? fraction : 1.0f - fraction) /
        RELAX_SOURCE_SCALE;
    if (grid_distance < minimum)
      minimum = grid_distance;
  };
  for (const math::Vector &v : mesh.vertices) {
    measure(v.x);
    measure(v.y);
    measure(v.z);
  }
  return minimum;
}

/**
 * @brief Checks that a relax bake belongs to the source mesh's vertices.
 * @param mesh Source mesh to identify.
 * @param bake Bake carrying the expected source identity.
 */
inline void check_relax_bake_source(const PolyMesh &mesh,
                                    const RelaxBake &bake) {
  HS_CHECK(relax_source_hash(mesh) == bake.source_hash,
           "relax_baked: source vertices differ");
}

HS_O3_END
} // namespace MeshOps
