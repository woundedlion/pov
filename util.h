#pragma once

#include "3dmath.h"

namespace hs {
  double rand_dbl() {
    return static_cast<double>(::random(0, std::numeric_limits<int32_t>::max()))
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

double shortest_distance(double a, double b, double m) {
  double d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

constexpr double fwd_distance(double a, double b, double m) {
  auto d = b - a;
  if (d < 0) {
    d += m;
  }
  return d;
}
