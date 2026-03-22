// Auto-generated -- do not edit
// Fibonacci lattice K-NN graph: RD_N=7680, RD_K=6
#pragma once
#ifndef HOLOSPHERE_CORE_REACTION_GRAPH_H_
#define HOLOSPHERE_CORE_REACTION_GRAPH_H_
#include "platform.h"
#include "3dmath.h"
#include "memory.h"
#include <cmath>

namespace ReactionGraph {

static constexpr int RD_N = 7680;
static constexpr int RD_K = 6;

/** Compute Fibonacci lattice node i (avoids storing 180KB in flash). */
inline Vector node(int i) {
  constexpr float phi = 2.399963229728653f;
  float y = 1.0f - (static_cast<float>(i) / (RD_N - 1)) * 2.0f;
  float radius = sqrtf(1.0f - y * y);
  float theta = phi * i;
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
  uint16_t *data = nullptr;

  /** Allocate and populate the LUT from the given arena (48 KB). */
  void build(Arena &arena) {
    data = static_cast<uint16_t *>(
        arena.allocate(6 * RES * RES * sizeof(uint16_t), alignof(uint16_t)));
    for (int face = 0; face < 6; ++face) {
      for (int y = 0; y < RES; ++y) {
        for (int x = 0; x < RES; ++x) {
          float u = (x + 0.5f) / RES * 2.0f - 1.0f;
          float v = (y + 0.5f) / RES * 2.0f - 1.0f;
          data[(face * RES + y) * RES + x] =
              static_cast<uint16_t>(find_nearest_node(texel_direction(face, u, v)));
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
    return data[(face * RES + vi) * RES + ui];
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
    float len = sqrtf(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    dir.x /= len; dir.y /= len; dir.z /= len;
    return dir;
  }

  /** Greedy walk from a latitude seed to the true nearest Fibonacci node. */
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
        if (d < best_d) { best_d = d; cur = ni; improved = true; }
      }
      if (!improved) break;
    }
    return cur;
  }
};

} // namespace ReactionGraph
#endif // HOLOSPHERE_CORE_REACTION_GRAPH_H_
