#pragma once

#include "3dmath.h"

namespace hs {
  float rand_dbl() {
    return static_cast<float>(::random(0, std::numeric_limits<int32_t>::max()))
      / std::numeric_limits<int32_t>::max();
  }

  int rand_int(int min, int max) {
    return ::random(min, max);
  }
}

template <typename T, typename U>
inline T wrap(T x, U m) {
  if (std::abs(x) < TOLERANCE) {
    return 0;
  }

  T r = std::fmod(x, m);
  if (r < 0) {
    r += m;
  }
  return (r >= m) ? 0 : r;
}

float shortest_distance(float a, float b, float m) {
  float d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

constexpr float fwd_distance(float a, float b, float m) {
  auto d = b - a;
  if (d < 0) {
    d += m;
  }
  return d;
}
