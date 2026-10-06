/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#pragma once

/**
 * @file reaction_graph.h
 * @brief Fibonacci lattice K-NN graph and nearest-node lookup.
 *
 * node() mirrors the lattice math of scripts/generate_reaction_graph.py, which
 * emits neighbors[].
 */

#include "platform/platform.h"
#include "math/3dmath.h"
#include "memory.h"
#include <cassert>
#include <cmath>

namespace ReactionGraph {

// Changing RD_N requires regenerating neighbors[] and updating D_AVG.
inline constexpr int RD_N = 7680;
inline constexpr int RD_K = 6;

/**
 * @brief Characteristic spacing sqrt(4π / RD_N) for an RD_N-point unit-sphere
 *        lattice, not a measured mean nearest-neighbor distance.
 */
inline constexpr float D_AVG = 0.0404505398f; // sqrt(4π / 7680)

// node()'s `(RD_N - 1)` divisor degenerates (divide-by-zero) when RD_N <= 1.
static_assert(RD_N >= 2, "node() lattice mapping degenerates when RD_N <= 1");

// neighbors[] elements and stored node indices are int16_t.
static_assert(RD_N <= INT16_MAX, "node index must fit int16_t");

// D_AVG = sqrt(4π / RD_N), so D_AVG*D_AVG*RD_N must equal 4π (12.566...). sqrtf
// isn't constexpr here, so check the squared, multiply-only form.
static_assert(D_AVG * D_AVG * RD_N - 12.566370614f < 0.0006f &&
                  12.566370614f - D_AVG * D_AVG * RD_N < 0.0006f,
              "D_AVG out of sync with RD_N (must stay sqrt(4*pi / RD_N))");

/**
 * @brief Computes Fibonacci-lattice node i as a unit vector on the sphere.
 * @param i Node index in [0, RD_N), ordered from north pole (i=0) southward;
 *        out of range traps.
 * @return Unit-length direction for lattice point i.
 * @details Analytic reference for the generated node_positions array.
 */
HS_COLD_MEMBER inline math::Vector node(int i) {
  HS_CHECK(i >= 0 && i < RD_N, "node() index outside the lattice");
  // Folds y, radius, and theta in double to reproduce neighbors[] bit-for-bit:
  // float32 flips near-tie sort order, and theta = golden_angle*i reaches ~18,400
  // rad at i=RD_N-1. Bit-exactness is a pinned-toolchain provenance contract;
  // the runtime tolerates ULP drift.
  constexpr double golden_angle = 2.399963229728653;
  constexpr double two_pi = 6.283185307179586;
  double y = 1.0 - (static_cast<double>(i) / (RD_N - 1)) * 2.0;
  // y=±1 at the poles (i=0, i=RD_N-1) collapses radius to 0, placing the node
  // exactly on the axis regardless of theta.
  double radius = std::sqrt(1.0 - y * y);
  double theta = std::fmod(golden_angle * i, two_pi);
  return math::Vector(static_cast<float>(std::cos(theta) * radius),
                      static_cast<float>(y),
                      static_cast<float>(std::sin(theta) * radius));
}

/** @brief Generated Fibonacci-lattice positions in flash. */
extern HS_PROGMEM_UNIQUE(node_positions) const math::Vector
    node_positions[RD_N];

/**
 * @brief Precomputed K-nearest-neighbor indices for every lattice node.
 * @details neighbors[i][k] is the node index of the k-th nearest neighbor of node
 *          i. Every slot is a node index in [0, RD_N), with no vacant-slot
 *          sentinel; validate_neighbors() checks this.
 */
extern HS_PROGMEM_UNIQUE(neighbors) const int16_t neighbors[RD_N][RD_K];

/** @brief Consecutive nodes sharing an ordered tuple of neighbor offsets. */
struct NeighborRun {
  uint16_t end;        /**< Exclusive node index; first run starts at zero. */
  int16_t delta[RD_K]; /**< neighbors[i][k] - i, in nearest-first order. */
};

/** @brief Generated lossless run encoding for sequential neighbor sweeps. */
extern HS_PROGMEM_UNIQUE(neighbor_runs) const NeighborRun neighbor_runs[];
extern HS_PROGMEM_UNIQUE(neighbor_run_count) const unsigned NEIGHBOR_RUN_COUNT;

/** @brief Node-to-run map for random access to ordered neighbor offsets. */
extern HS_PROGMEM_UNIQUE(neighbor_run_index) const uint8_t
    neighbor_run_index[RD_N];

/**
 * @brief Traps unless every slot of a neighbor table is a lattice node index.
 * @param table Neighbor rows to check, RD_N rows of RD_K indices.
 */
HS_COLD_MEMBER inline void
validate_neighbors(const int16_t (&table)[RD_N][RD_K]) {
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k)
      HS_CHECK(table[i][k] >= 0 && table[i][k] < RD_N,
               "neighbors[] slot is not a lattice node index");
}

// ---------------------------------------------------------------------------
// CubemapLUT: O(1) direction → nearest Fibonacci node lookup (no runtime trig)
// ---------------------------------------------------------------------------

/**
 * @brief 6-face cubemap LUT mapping unit vectors to nearest lattice node
 *        indices, giving O(1) direction → node lookup with no runtime trig.
 * @details Memory: 6 × RES² × 2B = 48 KB at RES=64.
 */
struct CubemapLUT {
  static constexpr int RES = 64;

  struct Projection {
    int face;
    float u;
    float v;
  };

  /**
   * @brief Allocates and populates the LUT from the given arena (48 KB
   *        persistent + ~90 KB transient).
   * @param arena Arena providing backing storage for the 6×RES² table (48 KB,
   *        retained) plus a transient ~90 KB lattice scratch, scoped to build() and
   *        rewound on return. A caller must provision for the peak (~138 KB), not
   *        the 48 KB persistent table alone, or this traps mid-build.
   * @details Texels are filled in lookup()'s (face*RES+y)*RES+x index order.
   */
  HS_COLD_MEMBER void build(Arena &arena) {
    validate_neighbors(neighbors);
    data.bind(arena, 6 * RES * RES);
    // Precompute every lattice point once into scratch for the hill-climb.
    ScratchScope lattice_guard(arena);
    math::Vector *lattice = arena.allocate_n<math::Vector>(RD_N);
    for (int i = 0; i < RD_N; ++i)
      lattice[i] = node(i);
    fill(lattice);
  }

  /** @brief Builds from a resident lattice whose neighbors are validated. */
  HS_COLD_MEMBER void build(Arena &arena, const math::Vector *lattice) {
    data.bind(arena, 6 * RES * RES);
    fill(lattice);
  }

private:
  /** @brief Fills storage allocated by build(). */
  HS_COLD_MEMBER void fill(const math::Vector *lattice) {
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < RES; ++y) {
        int seed = -1;
        for (int x = 0; x < RES; ++x) {
          float u = (x + 0.5f) / RES * 2.0f - 1.0f;
          float v = (y + 0.5f) / RES * 2.0f - 1.0f;
          seed = find_nearest_node(texel_direction(face, u, v), lattice, seed);
          data.push_back(static_cast<uint16_t>(seed));
        }
      }
    }
  }

public:
  /**
   * @brief O(1) cubemap lookup projecting a unit vector to an approximately
   *        nearest lattice node.
   * @param p Query direction; MUST be unit-length. The dominant-axis magnitude
   *        divides with no zero-guard; a unit `p` keeps it >= 1/sqrt(3).
   * @return A seed lattice node index in [0, RD_N) close to p. The table holds
   *         hill-climb local minima at quantized face cells, so the true
   *         nearest node can lie beyond the seed's one-ring.
   */
  int lookup(const math::Vector &p) const { return lookup(project(p)); }

  /** @brief Projects a unit direction into the lattice cubemap face coordinates. */
  static __attribute__((always_inline)) Projection
  project(const math::Vector &p) {
    assert(std::fabs(p.x * p.x + p.y * p.y + p.z * p.z - 1.0f) < 1e-3f);
    float ax = fabsf(p.x), ay = fabsf(p.y), az = fabsf(p.z);
    int face = 0;
    float u = 0, v = 0;

    if (ax >= ay && ax >= az) {
      float inv = 1.0f / ax;
      if (p.x >= 0) {
        face = 0;
        u = -p.z * inv;
        v = p.y * inv;
      } else {
        face = 1;
        u = p.z * inv;
        v = p.y * inv;
      }
    } else if (ay >= ax && ay >= az) {
      float inv = 1.0f / ay;
      if (p.y >= 0) {
        face = 2;
        u = p.x * inv;
        v = -p.z * inv;
      } else {
        face = 3;
        u = p.x * inv;
        v = p.z * inv;
      }
    } else {
      float inv = 1.0f / az;
      if (p.z >= 0) {
        face = 4;
        u = p.x * inv;
        v = p.y * inv;
      } else {
        face = 5;
        u = -p.x * inv;
        v = p.y * inv;
      }
    }

    return {face, u, v};
  }

  /** @brief Looks up a previously projected lattice cubemap coordinate. */
  __attribute__((always_inline)) int
  lookup(const Projection &projection) const {
    const int FACE = projection.face;
    const float U = projection.u, V = projection.v;
    // Clamp nonfinite coordinates before converting them to integer texels.
    int ui = static_cast<int>(
        hs::clamp((U + 1.0f) * 0.5f * RES, 0.0f, static_cast<float>(RES - 1)));
    int vi = static_cast<int>(
        hs::clamp((V + 1.0f) * 0.5f * RES, 0.0f, static_cast<float>(RES - 1)));
    return data[static_cast<size_t>((FACE * RES + vi) * RES + ui)];
  }

private:
  /**
   * @brief Texel-to-node table, indexed (face*RES+y)*RES+x.
   */
  ArenaVector<uint16_t> data;

  /**
   * @brief Converts a cubemap face plus texel coordinates to a unit direction.
   * @param face Cube face index in [0, 6): +X,-X,+Y,-Y,+Z,-Z.
   * @param u Horizontal texel coordinate in [-1, 1].
   * @param v Vertical texel coordinate in [-1, 1].
   * @return Unit-length direction vector for the texel.
   * @details One axis is always ±1, so the length is >= 1.
   */
  static math::Vector texel_direction(int face, float u, float v) {
    math::Vector dir;
    if (face == 0)
      dir = math::Vector(1.0f, v, -u); // +X
    else if (face == 1)
      dir = math::Vector(-1.0f, v, u); // -X
    else if (face == 2)
      dir = math::Vector(u, 1.0f, -v); // +Y
    else if (face == 3)
      dir = math::Vector(u, -1.0f, v); // -Y
    else if (face == 4)
      dir = math::Vector(u, v, 1.0f); // +Z
    else
      dir = math::Vector(-u, v, -1.0f); // -Z
    return dir.normalized();
  }

  /**
   * @brief Finds the near-nearest Fibonacci node to p via greedy K-NN descent
   *        from a supplied node or a latitude seed.
   * @param p Query direction (expected unit-length) on the sphere.
   * @param lattice Precomputed node() positions for all RD_N points, indexed by
   *        node id.
   * @param seed Starting node, or -1 to choose by latitude.
   * @return Lattice node index in [0, RD_N) at a local distance minimum.
   * @details Hill-climbs toward closer neighbors and stops at a local minimum; on
   *          the near-uniform Fibonacci sphere this lands on the true nearest node
   *          in practice but is not guaranteed to (not a global argmin).
   *
   *          The inner loop moves `cur` the instant it sees a closer neighbor, so
   *          one `iter` can chain several hops; the 64-iteration cap depends on
   *          this, since the latitude-only seed can start dozens of hops from an
   *          equatorial query's node.
   */
  HS_COLD_MEMBER static int find_nearest_node(const math::Vector &p,
                                              const math::Vector *lattice,
                                              int seed) {
    int cur = seed >= 0 ? seed
                        : static_cast<int>(
                              hs::clamp((1.0f - p.y) * 0.5f * (RD_N - 1) + 0.5f,
                                        0.0f, static_cast<float>(RD_N - 1)));
    float best_d = math::distance_squared(p, lattice[cur]);
    bool converged = false;
    for (int iter = 0; iter < 64; ++iter) {
      bool improved = false;
      for (int k = 0; k < RD_K; ++k) {
        int ni = neighbors[cur][k];
        float d = math::distance_squared(p, lattice[ni]);
        if (d < best_d) {
          best_d = d;
          cur = ni;
          improved = true;
        }
      }
      if (!improved) {
        converged = true;
        break;
      }
    }
    HS_CHECK(converged,
             "find_nearest_node hit the 64-iter cap without converging "
             "(RD_N grew past the calibrated hop bound)");
    return cur;
  }
};

} // namespace ReactionGraph
