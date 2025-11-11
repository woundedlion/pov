#pragma once

namespace hs {
  double rand_dbl() {
    return static_cast<double>(::random(0, std::numeric_limits<int32_t>::max()))
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

double wrap(double x, int m) {
  return x >= 0 ?
    fmod(x, m) :
    fmod(fmod(x, m) + m, m);
}

double wrap(double x, double m) {
  return x >= 0 ?
    fmod(x, m) :
    fmod(fmod(x, m) + m, m);
}

double shortest_distance(double a, double b, double m) {
  double d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

double fwd_distance(double a, double b, double m) {
  auto d = b - a;
  if (d < 0) {
    d += m;
  }
  return d;
}