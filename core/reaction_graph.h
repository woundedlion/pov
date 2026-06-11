// Auto-generated -- do not edit
// Fibonacci lattice K-NN graph: RD_N=7680, RD_K=6
#pragma once
#include "platform.h"
#include "3dmath.h"
#include "memory.h"
#include <cmath>

namespace ReactionGraph {

static constexpr int RD_N = 7680;
static constexpr int RD_K = 6;

/// Mean nearest-neighbor spacing on the unit sphere for an RD_N-point lattice:
/// sqrt(4π / RD_N). Used as the base radius for reaction-diffusion interpolation
/// kernels (BZ / GS). Precomputed because std::sqrt isn't constexpr here.
static constexpr float D_AVG = 0.04045f; // sqrt(4π / 7680)

/** Compute Fibonacci lattice node i (avoids storing 90KB in flash). */
inline Vector node(int i) {
  // theta = golden_angle·i reaches ~18,400 rad at i=RD_N-1, where a float holds
  // only ~1e-3 rad of azimuth — enough to disagree with the offline
  // double-precision generator that built neighbors[]. Accumulate and fold to
  // [0, 2π) in double, then narrow once for the trig. node() runs only at init
  // (CubemapLUT::build, build_nodes), so the wider math costs nothing per frame.
  constexpr double golden_angle = 2.399963229728653;
  constexpr double two_pi = 6.283185307179586;
  float y = 1.0f - (static_cast<float>(i) / (RD_N - 1)) * 2.0f;
  float radius = sqrtf(1.0f - y * y);
  float theta = static_cast<float>(std::fmod(golden_angle * i, two_pi));
  return Vector(cosf(theta) * radius, y, sinf(theta) * radius);
}

/** Precomputed K-NN neighbor indices (92160 bytes). */
extern const int16_t neighbors[RD_N][RD_K];

// ---------------------------------------------------------------------------
// CubemapLUT: O(1) direction → nearest Fibonacci node lookup (no runtime trig)
// ---------------------------------------------------------------------------

/** 6-face cubemap LUT mapping unit vectors to nearest lattice node indices.
 *  Memory: 6 × RES² × 2B = 48 KB at RES=64. */
struct CubemapLUT {
  static constexpr int RES = 64;
  // Arena-backed rather than a bare uint16_t*: inherits ArenaVector's debug
  // generation tracking (use-after-free if the arena is reset out from under the
  // LUT) and bound-check (lookup() before build() traps via operator[]'s
  // check_bound, instead of a silent null/garbage read). On device (NDEBUG) the
  // checks compile out, so lookup() is the same single load as before.
  ArenaVector<uint16_t> data;

  /** Allocate and populate the LUT from the given arena (48 KB). */
  void build(Arena &arena) {
    // The triple loop visits texels in linear (face, y, x) order — exactly the
    // index (face*RES+y)*RES+x used by lookup() — so sequential push_back fills
    // the table in place with no random-access writes.
    data.bind(arena, 6 * RES * RES);
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < RES; ++y) {
        for (int x = 0; x < RES; ++x) {
          float u = (x + 0.5f) / RES * 2.0f - 1.0f;
          float v = (y + 0.5f) / RES * 2.0f - 1.0f;
          data.push_back(
              static_cast<uint16_t>(find_nearest_node(texel_direction(face, u, v))));
        }
      }
    }
  }

  /** O(1) cubemap lookup: project a unit vector to the nearest node index. */
  int lookup(const Vector &p) const {
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
  /** Convert a cubemap face + texel to a unit direction vector. */
  static Vector texel_direction(int face, float u, float v) {
    Vector dir;
    if (face == 0)      dir = Vector( 1.0f, v, -u);  // +X
    else if (face == 1) dir = Vector(-1.0f, v,  u);  // -X
    else if (face == 2) dir = Vector(u,  1.0f, -v);  // +Y
    else if (face == 3) dir = Vector(u, -1.0f,  v);  // -Y
    else if (face == 4) dir = Vector(u, v,  1.0f);   // +Z
    else                dir = Vector(-u, v, -1.0f);   // -Z
    // One axis is always ±1, so the length is >= 1; the strict normalized()
    // traps if that invariant is ever broken instead of dividing by zero.
    return dir.normalized();
  }

  /** Near-nearest Fibonacci node via greedy K-NN descent from a latitude seed.
   *  Hill-climbs toward closer neighbors and stops at a local minimum, so on
   *  the near-uniform Fibonacci sphere this lands on the true nearest node in
   *  practice but is not guaranteed to — it is not a global argmin.
   *
   *  LOAD-BEARING: the inner loop reassigns `cur` (and `best_d`) the instant it
   *  sees a closer neighbor, so each later `k` reads `neighbors[cur]` of the
   *  just-updated node — one `iter` therefore chains through *several* nodes,
   *  not one. This is not a textbook best-of-neighbors step and the 64-iteration
   *  cap depends on it: the seed is latitude-only, so for an equatorial query
   *  the true node can sit dozens of hops away around the longitude circle (the
   *  RD_N=7680 lattice has long iso-latitude rings). A refactor that scanned all
   *  RD_K neighbors of a *fixed* `cur` before moving would advance one hop per
   *  iter, exceed 64 hops near the equator, and silently return a wrong node —
   *  which the round-trip test cannot catch because it seeds at the answer. Do
   *  not convert this to best-of-neighbors without also raising the cap. */
  static int find_nearest_node(const Vector &p) {
    auto dist2 = [](const Vector &a, const Vector &b) {
      float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
      return dx * dx + dy * dy + dz * dz;
    };
    int cur = hs::clamp(
        static_cast<int>((1.0f - p.y) * 0.5f * (RD_N - 1) + 0.5f), 0, RD_N - 1);
    float best_d = dist2(p, node(cur));
    for (int iter = 0; iter < 64; ++iter) {
      bool improved = false;
      for (int k = 0; k < RD_K; ++k) {
        int ni = neighbors[cur][k];
        if (ni < 0) continue;
        float d = dist2(p, node(ni));
        // Mid-scan reassignment is deliberate (see LOAD-BEARING note above):
        // subsequent k read neighbors of `cur`, chaining multiple hops per iter.
        if (d < best_d) { best_d = d; cur = ni; improved = true; }
      }
      if (!improved) break;
    }
    return cur;
  }
};

} // namespace ReactionGraph
