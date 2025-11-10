#pragma once

namespace hs {
  float_t rand_dbl() {
    return static_cast<float_t>(::random(0, std::numeric_limits<int32_t>::max()))
      / std::numeric_limits<int32_t>::max();
  }

  int rand_int(int min, int max) {
    return ::random(min, max);
  }
}

int wrap(int x, int m) {
  return (x >= 0 ?
    x % m :
    ((x % m) + m) % m);
}

float_t wrap(float_t x, int m) {
  return x >= 0 ?
    fmod(x, m) :
    fmod(fmod(x, m) + m, m);
}

float_t wrap(float_t x, float_t m) {
  return x >= 0 ?
    fmod(x, m) :
    fmod(fmod(x, m) + m, m);
}

float_t shortest_distance(float_t a, float_t b, float_t m) {
  float_t d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

float_t fwd_distance(float_t a, float_t b, float_t m) {
  auto d = b - a;
  if (d < 0) {
    d += m;
  }
  return d;
}