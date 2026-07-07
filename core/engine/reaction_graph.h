// Fibonacci lattice K-NN graph. The neighbors[] table (reaction_graph.cpp) is
// emitted by scripts/generate_reaction_graph.py; node() below MUST stay in
// lockstep with that generator's lattice math (CI: reaction-graph-provenance).
#pragma once
#include "engine/platform.h"
#include "math/3dmath.h"
#include "engine/memory.h"
#include <cassert>
#include <cmath>

namespace ReactionGraph {

// Changing RD_N requires regenerating neighbors[] (generate_reaction_graph.py)
// and re-pasting D_AVG below — both are guarded, neither is auto-derived.
static constexpr int RD_N = 7680;
static constexpr int RD_K = 6;

/**
 * @brief Mean nearest-neighbor spacing on the unit sphere for an RD_N-point
 *        lattice: sqrt(4π / RD_N), in unit-sphere radians.
 * @details Used as the base radius for reaction-diffusion interpolation kernels
 *          (BZ / GS). Precomputed because std::sqrt isn't constexpr here.
 */
static constexpr float D_AVG = 0.0404505398f; // sqrt(4π / 7680)

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
 * @param i Node index in [0, RD_N), ordered from north pole (i=0) southward.
 * @return Unit-length direction for lattice point i.
 * @details Recomputed on demand to avoid storing 90KB in flash. Only ever called
 *          at init, so the double-precision folding below is off the per-frame
 *          render path (zero per-frame cost); see the implementation note for why
 *          the wider math is needed.
 */
inline Vector node(int i) {
  // Must fold y, radius, and theta in double to reproduce neighbors[] bit-for-bit:
  // float32 flips near-tie sort order, and theta = golden_angle*i reaches ~18,400
  // rad at i=RD_N-1 where a float holds only ~1e-3 rad of azimuth.
  //
  // Bit-exact reproduction of the generated neighbors[] is a CI-pinned-toolchain
  // provenance contract, not a portable runtime guarantee: it needs host x86,
  // device ARM, and the numpy table generator to agree to the last ULP on
  // cos/sin/fmod, which is not promised for transcendentals across libms. The
  // runtime tolerates disagreement — find_nearest_node is a nearest-node search,
  // not a table-index equality — so a ULP drift degrades seeding quality at
  // worst, never correctness.
  constexpr double golden_angle = 2.399963229728653;
  constexpr double two_pi = 6.283185307179586;
  double y = 1.0 - (static_cast<double>(i) / (RD_N - 1)) * 2.0;
  // y=±1 at the poles (i=0, i=RD_N-1) collapses radius to 0, placing the node
  // exactly on the axis regardless of theta.
  double radius = std::sqrt(1.0 - y * y);
  double theta = std::fmod(golden_angle * i, two_pi);
  return Vector(static_cast<float>(std::cos(theta) * radius),
                static_cast<float>(y),
                static_cast<float>(std::sin(theta) * radius));
}

/**
 * @brief Precomputed K-nearest-neighbor indices for every lattice node.
 * @details neighbors[i][k] is the node index of the k-th nearest neighbor of
 *          node i. Every entry is a valid index in [0, RD_N): the RD_N=7680
 *          lattice always yields a full RD_K-neighbor ring, so the table is fully
 *          populated and the `ni < 0` guard in find_nearest_node is a defensive
 *          belt-and-suspenders, not a live table state. Total size 92160 bytes
 *          (RD_N × RD_K × 2B). PROGMEM matches the generated definition; it is a
 *          no-op on the supported flat-address targets, where the direct
 *          `neighbors[i][k]` subscripting below is intentional.
 */
extern PROGMEM const int16_t neighbors[RD_N][RD_K];

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

  /**
   * @brief Allocates and populates the LUT from the given arena (48 KB
   *        persistent + ~90 KB transient).
   * @param arena Arena providing backing storage for the 6×RES² table (48 KB,
   *        retained) plus a transient RD_N×sizeof(Vector) (~90 KB) lattice
   *        scratch, scoped to build() and rewound on return. A caller must
   *        provision for the peak (~138 KB), not the 48 KB persistent table
   *        alone, or this traps mid-build.
   * @details The triple loop visits texels in linear (face, y, x) order — exactly
   *          the index (face*RES+y)*RES+x used by lookup() — so sequential
   *          push_back fills the table in place with no random-access writes.
   */
  void build(Arena &arena) {
    data.bind(arena, 6 * RES * RES);
    // node() is double-precision trig; precompute every lattice point once into
    // scratch so the hill-climb reads a table instead of recomputing per hop.
    ScratchScope lattice_guard(arena);
    Vector *lattice = static_cast<Vector *>(
        arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    for (int i = 0; i < RD_N; ++i)
      lattice[i] = node(i);
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < RES; ++y) {
        for (int x = 0; x < RES; ++x) {
          float u = (x + 0.5f) / RES * 2.0f - 1.0f;
          float v = (y + 0.5f) / RES * 2.0f - 1.0f;
          data.push_back(static_cast<uint16_t>(
              find_nearest_node(texel_direction(face, u, v), lattice)));
        }
      }
    }
  }

  /**
   * @brief O(1) cubemap lookup projecting a unit vector to an approximately
   *        nearest lattice node.
   * @param p Query direction; MUST be unit-length. This precondition is
   *        load-bearing for the divide below: the dominant-axis magnitude is
   *        used as the divisor with no zero-guard, and only a unit `p`
   *        guarantees that axis is `>= 1/sqrt(3)` (~0.577), so the reciprocal
   *        never blows up. A zero or denormal `p` would divide by ~0.
   * @return A seed lattice node index in [0, RD_N) close to p — the texel's
   *         precomputed node. Two approximations stack: the table is built from
   *         find_nearest_node (a hill-climb local minimum, not a global argmin)
   *         and the query is quantized to a face cell, so the result can be off
   *         by a neighbor. Callers needing the true nearest node must refine
   *         among the seed and its neighbors (see
   *         ReactionDiffusionBase::refine_nearest_node in effects/).
   */
  int lookup(const Vector &p) const {
    // debug-only: device hot path stays a single load (assert compiles out).
    assert(std::fabs(p.x * p.x + p.y * p.y + p.z * p.z - 1.0f) < 1e-3f);
    float ax = fabsf(p.x), ay = fabsf(p.y), az = fabsf(p.z);
    int face = 0;
    float u = 0, v = 0;

    if (ax >= ay && ax >= az) {
      float inv = 1.0f / ax;
      if (p.x >= 0) { face = 0; u = -p.z * inv; v = p.y * inv; }
      else          { face = 1; u =  p.z * inv; v = p.y * inv; }
    } else if (ay >= ax && ay >= az) {
      float inv = 1.0f / ay;
      if (p.y >= 0) { face = 2; u = p.x * inv; v = -p.z * inv; }
      else          { face = 3; u = p.x * inv; v =  p.z * inv; }
    } else {
      float inv = 1.0f / az;
      if (p.z >= 0) { face = 4; u =  p.x * inv; v = p.y * inv; }
      else          { face = 5; u = -p.x * inv; v = p.y * inv; }
    }

    int ui = hs::clamp(static_cast<int>((u + 1.0f) * 0.5f * RES), 0, RES - 1);
    int vi = hs::clamp(static_cast<int>((v + 1.0f) * 0.5f * RES), 0, RES - 1);
    return data[static_cast<size_t>((face * RES + vi) * RES + ui)];
  }

private:
  /**
   * @brief Texel-to-node table, indexed (face*RES+y)*RES+x.
   * @details Arena-backed rather than a bare uint16_t*: inherits ArenaVector's
   *          debug generation tracking (use-after-free if the arena is reset out
   *          from under the LUT) and bound-check (lookup() before build() traps
   *          via operator[]'s check_bound, instead of a silent null/garbage
   *          read). On device (NDEBUG) the checks compile out, so lookup() is a
   *          single load.
   */
  ArenaVector<uint16_t> data;

  /**
   * @brief Converts a cubemap face plus texel coordinates to a unit direction.
   * @param face Cube face index in [0, 6): +X,-X,+Y,-Y,+Z,-Z.
   * @param u Horizontal texel coordinate in [-1, 1].
   * @param v Vertical texel coordinate in [-1, 1].
   * @return Unit-length direction vector for the texel.
   * @details One axis is always ±1, so the length is >= 1; the strict
   *          normalized() traps if that invariant is ever broken instead of
   *          dividing by zero.
   */
  static Vector texel_direction(int face, float u, float v) {
    Vector dir;
    if (face == 0)      dir = Vector( 1.0f, v, -u);  // +X
    else if (face == 1) dir = Vector(-1.0f, v,  u);  // -X
    else if (face == 2) dir = Vector(u,  1.0f, -v);  // +Y
    else if (face == 3) dir = Vector(u, -1.0f,  v);  // -Y
    else if (face == 4) dir = Vector(u, v,  1.0f);   // +Z
    else                dir = Vector(-u, v, -1.0f);  // -Z
    return dir.normalized();
  }

  /**
   * @brief Finds the near-nearest Fibonacci node to p via greedy K-NN descent
   *        from a latitude seed.
   * @param p Query direction (expected unit-length) on the sphere.
   * @param lattice Precomputed node() positions for all RD_N points, indexed by
   *        node id; built once by the caller to avoid recomputing node()'s trig
   *        per hop.
   * @return Lattice node index in [0, RD_N) at a local distance minimum.
   * @details Hill-climbs toward closer neighbors and stops at a local minimum, so
   *          on the near-uniform Fibonacci sphere this lands on the true nearest
   *          node in practice but is not guaranteed to — it is not a global
   *          argmin.
   *
   *          LOAD-BEARING: the inner loop reassigns `cur` (and `best_d`) the
   *          instant it sees a closer neighbor, so each later `k` reads
   *          `neighbors[cur]` of the just-updated node — one `iter` therefore
   *          chains through several nodes, not one. This is not a textbook
   *          best-of-neighbors step and the 64-iteration cap depends on it: the
   *          seed is latitude-only, so for an equatorial query the true node can
   *          sit dozens of hops away around the longitude circle (the RD_N=7680
   *          lattice has long iso-latitude rings). A refactor that scanned all
   *          RD_K neighbors of a fixed `cur` before moving would advance one hop
   *          per iter, exceed 64 hops near the equator, and silently return a
   *          wrong node — which the round-trip test cannot catch because it seeds
   *          at the answer. Do not convert this to best-of-neighbors without also
   *          raising the cap.
   */
  static int find_nearest_node(const Vector &p, const Vector *lattice) {
    auto dist2 = [](const Vector &a, const Vector &b) {
      float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      return dx * dx + dy * dy + dz * dz;
    };
    int cur = static_cast<int>(hs::clamp(
        (1.0f - p.y) * 0.5f * (RD_N - 1) + 0.5f, 0.0f,
        static_cast<float>(RD_N - 1)));
    float best_d = dist2(p, lattice[cur]);
    bool converged = false;
    for (int iter = 0; iter < 64; ++iter) {
      bool improved = false;
      for (int k = 0; k < RD_K; ++k) {
        int ni = neighbors[cur][k];
        if (ni < 0) continue;
        float d = dist2(p, lattice[ni]);
        if (d < best_d) { best_d = d; cur = ni; improved = true; }
      }
      if (!improved) { converged = true; break; }
    }
    HS_CHECK(converged, "find_nearest_node hit the 64-iter cap without converging "
                        "(RD_N grew past the calibrated hop bound)");
    return cur;
  }
};

} // namespace ReactionGraph
