// Auto-generated -- do not edit
// Fibonacci lattice K-NN graph: RD_N=3840, RD_K=6
#pragma once
#include "platform.h"
#include "3dmath.h"
#include <cmath>

namespace ReactionGraph {

static constexpr int RD_N = 3840;
static constexpr int RD_K = 6;

/** Compute Fibonacci lattice node i (avoids storing 46KB in flash). */
inline Vector node(int i) {
  constexpr float phi = 2.399963229728653f;
  float y = 1.0f - (static_cast<float>(i) / (RD_N - 1)) * 2.0f;
  float radius = sqrtf(1.0f - y * y);
  float theta = phi * i;
  return Vector(cosf(theta) * radius, y, sinf(theta) * radius);
}

/** Precomputed K-NN neighbor indices (46080 bytes). */
extern const int16_t neighbors[RD_N][RD_K];

} // namespace ReactionGraph
