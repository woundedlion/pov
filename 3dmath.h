#pragma once

#include <cmath>
#include <cfloat>

static const float pi = 3.1415926535897f;
static const float tau = 2 * pi;

inline float mod_tau(float n) {
  if (n > tau)
    return n - tau;
  else if (n < 0)
    return n + tau;
  return n;
}

inline float modf(float x, int m) {
  return (static_cast<int>(x) % m) + (x - static_cast<int>(x));
}

float degrees_to_radians(uint32_t deg) {
  return (deg % 360) * pi / 180;
}

float radians_to_degrees(float rad) {
  return mod_tau(rad) * 180 / pi;
}

struct Vector {
  Vector(float i, float j, float k) : i(i), j(j), k(k) {}
  Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  Vector& operator=(const Vector& v) { 
    i = v.i; 
    j = v.j, 
    k = v.k;
    return *this;
  }

  bool operator==(const Vector& v) const {
    return fabs(i - v.i) < FLT_EPSILON
      && fabs(j - v.j) < FLT_EPSILON
      && fabs(k - v.k) < FLT_EPSILON;
  }

  Vector operator-() const {
    return Vector(-i, -j, -k);
  }

  float length() const { return sqrt(i * i + j * j + k * k); }

  float i = 0;
  float j = 0;
  float k = 0;
};

Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
}

Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
}

Vector operator*(const Vector& v, float s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

Vector operator*(float s, const Vector& v) {
  return v * s;
}

Vector operator/(const Vector& v, float s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

float dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

struct Quaternion;
Quaternion operator*(const Quaternion& q1, const Quaternion& q2);

struct Quaternion {
  Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  Quaternion(float a, float b, float c, float d) : r(a), v(b, c, d) {}
  Quaternion(float a, const Vector& v) : r(a), v(v) {}
  
  Quaternion& operator=(const Quaternion& q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  bool operator==(const Quaternion& q) {
    return fabs(q.r - r) < FLT_EPSILON
      && q.v == v;
  }

  void operator*=(const Quaternion& q) {
    *this = *this * q;
  }

  Quaternion inverse() const { 
    float n = (r * r) + (v.i * v.i) + (v.j * v.j) + (v.k * v.k);
    return Quaternion(r / n, (-v) / n);
  }

  float r = 1;
  Vector v;
};

Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r + q2.r, q1.v + q2.v);
}

Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r - q2.r, q1.v - q2.v);
}

Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(
    q1.r * q2.r - dot(q1.v, q2.v),
    q1.r * q2.v + q2.r * q1.v + cross(q1.v, q2.v));
}

Quaternion operator*(const Quaternion& q, float s) {
  return Quaternion(q.r * s, q.v * s);
}

Quaternion operator*(float s, const Quaternion& q) {
  return q * s;
}

float dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

float angle_between(const Quaternion& q1, const Quaternion& q2) {
  return acosf(dot(q1, q2));
}

Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t) {
  float theta = angle_between(q1, q2);
  return(((sinf(1 - t) * theta) / sinf(theta)) * q1)
    + (((sinf(t) * theta) / sinf(theta)) * q2);
}

