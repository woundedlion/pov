#pragma once

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

double distance(double a, double b, double m) {
  return std::min(abs(b - a), m - abs(a - b));
}