#pragma once

#include <cmath>
#include <cfloat>
#include <limits>
#include <cassert>

static constexpr float PHI = 1.6180339887498948482045868343656;
static constexpr float G = 1 / PHI;
static constexpr float TOLERANCE = 0.0001f;
static constexpr float PI_F = static_cast<float>(PI);

struct Vector;
struct Spherical {
  Spherical(float theta, float phi) : theta(theta), phi(phi) {}
  Spherical(const Vector& v);

  float theta; // azimuthal angle
  float phi; // polar angle
};

struct Vector {
  constexpr Vector() {}
  constexpr Vector(float i, float j, float k) : i(i), j(j), k(k) {}
  constexpr Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  Vector(const Spherical& s) :
    i(sinf(s.phi)* cosf(s.theta)),
    j(sinf(s.phi)* sinf(s.theta)),
    k(cosf(s.phi))
  {}
  Vector(float theta, float phi) : Vector(Spherical(theta, phi)) {}

  constexpr Vector& operator=(const Vector& v) {
    i = v.i; 
    j = v.j;
    k = v.k;
    return *this;
  }

  bool operator==(const Vector& v) const {
    return std::abs(i - v.i) <= TOLERANCE
      && std::abs(j - v.j) <= TOLERANCE
      && std::abs(k - v.k) <= TOLERANCE;
  }

  bool operator!=(const Vector& v) const {
    return !(*this == v);
  }

  constexpr Vector operator-() const {
    return Vector(-i, -j, -k);
  }

  constexpr float length() const { return sqrt(i * i + j * j + k * k); }
  
  Vector& normalize() {
    float m = length();
    if (m < std::numeric_limits<float>::epsilon()) {
      Serial.println("Can't normalize a zero vector!");
      assert(false);
      i = 1;
      j = k = 0;
    } else {
      i = i / m;
      j = j / m;
      k = k / m;
    }
    return *this;
  }

  float i = 0;
  float j = 0;
  float k = 0;
};

Spherical::Spherical(const Vector& v)
{
  Vector n(v);
  n.normalize();
  theta = atan2f(n.j, n.i);
  phi = acosf(std::clamp(n.k, -1.0f, 1.0f));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
constexpr Quaternion operator/(const Quaternion& q, float s);
constexpr Vector operator/(const Vector& v, float s);

struct Quaternion {
  constexpr Quaternion() {}
  constexpr Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  constexpr  Quaternion(float a, float b, float c, float d) : r(a), v(b, c, d) {}
  constexpr Quaternion(float a, const Vector& v) : r(a), v(v) {}
  
  Quaternion& operator=(const Quaternion& q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  bool operator==(const Quaternion& q) const {
    return std::abs(q.r - r) < TOLERANCE
      && q.v == v;
  }

  void operator*=(const Quaternion& q) {
    *this = *this * q;
  }

  constexpr float squared_magnitude() const {
    return r * r + v.i * v.i + v.j * v.j + v.k * v.k;
  }

  Quaternion inverse() const {
    float sq_mag = squared_magnitude();
    return Quaternion(r, -v) / sq_mag;
  }

  constexpr Quaternion operator-() const {
    return Quaternion(-r, -v);
  }

  constexpr float magnitude() const {
    return sqrt(r * r + v.i * v.i + v.j * v.j + v.k * v.k);
  }

  Quaternion& normalize() {
    auto m = magnitude();
    if (m <= std::numeric_limits<float>::epsilon()) {
      Serial.println("Can't normalize a zaro Quaternion!");
      assert(false);
      r = 1;
      v = Vector(0, 0, 0);
    } else {
      r = r / m;
      v = v / m;
    }
    return *this;
  }

  float r = 1;
  Vector v;
};

Vector rotate(const Vector& v, const Quaternion& q) {
  assert(std::abs(q.magnitude() - 1) < TOLERANCE);
  Quaternion p(0, v);
  auto r = q * p * q.inverse();
  return r.v;
}

constexpr Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
}

constexpr Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
}

constexpr Vector operator*(const Vector& v, float s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

constexpr Vector operator*(float s, const Vector& v) {
  return v * s;
}

constexpr Vector operator/(const Vector& v, float s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

constexpr float dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

constexpr Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

constexpr float distance_between(const Vector& a, const Vector& b) {
  float di = b.i - a.i;
  float dj = b.j - a.j;
  float dk = b.k - a.k;
  return sqrt(di * di + dj * dj + dk * dk);
}

constexpr float angle_between(const Vector& v1, const Vector& v2) {
  float len_product = v1.length() * v2.length();
  if (len_product <= std::numeric_limits<float>::epsilon()) {
    assert(false);
  }
  float d = dot(v1, v2) / len_product;
  return acosf(std::clamp(d, -1.0f, 1.0f));
}

constexpr Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r + q2.r, q1.v + q2.v);
}

constexpr Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r - q2.r, q1.v - q2.v);
}

constexpr Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(
    q1.r * q2.r - dot(q1.v, q2.v),
    q1.r * q2.v + q2.r * q1.v + cross(q1.v, q2.v));
}

constexpr Quaternion operator*(const Quaternion& q, float s) {
  return Quaternion(q.r * s, q.v * s);
}

constexpr Quaternion operator*(float s, const Quaternion& q) {
  return q * s;
}

constexpr Quaternion operator/(const Quaternion& q, float s) {
  assert(std::abs(s) > std::numeric_limits<float>::epsilon());
  return Quaternion(q.r / s, q.v / s);
}

constexpr float dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

constexpr float angle_between(const Quaternion& q1, const Quaternion& q2) {
  assert(std::abs(q1.magnitude() - 1.0f) < TOLERANCE);
  assert(std::abs(q2.magnitude() - 1.0f) < TOLERANCE);
  return acosf(std::clamp(dot(q1, q2), -1.0f, 1.0f));
}

Quaternion make_rotation(const Vector& axis, float theta) {
  assert(std::abs(axis.length() - 1.0f) < TOLERANCE);
  return Quaternion(cosf(theta / 2), sinf(theta / 2) * axis).normalize();
}

Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t, bool long_way = false) {
  assert(std::abs(q1.magnitude() - 1.0f) < TOLERANCE);
  assert(std::abs(q2.magnitude() - 1.0f) < TOLERANCE);
  float d = dot(q1, q2);
  Quaternion p(q1);
  Quaternion q(q2);

  if ((long_way && d > 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d > (1 - TOLERANCE)) {
    Quaternion r = p + t * (q - p);
    return r.normalize();
  }

  float theta = acosf(std::clamp(d, -1.0f, 1.0f));
  float sin_theta = sinf(theta);
  float s1 = sinf((1 - t) * theta) / sin_theta;
  float s2 = sinf(t * theta) / sin_theta;
  return ((s1 * p) + (s2 * q)).normalize();
}

bool is_over(const Vector& v, const Vector& normal) {
  return dot(normal, v) >= 0;
}

bool intersects_plane(const Vector& v1, const Vector& v2, const Vector& normal) {
  return (is_over(v1, normal) && !is_over(v2, normal))
    || (!is_over(v1, normal) && is_over(v2, normal));
}

Vector intersection(const Vector& u, const Vector& v, const Vector& normal) {
  assert(std::abs(u.length() - 1.0) < TOLERANCE);
  assert(std::abs(v.length() - 1.0) < TOLERANCE);
  assert(std::abs(normal.length() - 1.0) < TOLERANCE);
  Vector w = cross(v, u).normalize();
  Vector i1 = cross(w, normal).normalize();
  Vector i2 = cross(normal, w).normalize();

  float a1 = angle_between(u, v);
  float a2 = angle_between(i1, u);
  float a3 = angle_between(i1, v);
  if (std::abs(a2 + a3 - a1) < .0001) {
    return i1;
  }

  a1 = angle_between(u, v);
  a2 = angle_between(i2, u);
  a3 = angle_between(i2, v);
  if (std::abs(a2 + a3 - a1) < TOLERANCE) {
    return i2;
  }

  assert(false);
  return Vector(0, 0, 0);
}

static constexpr float SHIFT_DIV = 96;

std::vector<Vector> split_point(const Vector& v, const Vector& normal) {
  assert(std::abs(v.length() - 1.0f) < TOLERANCE);
  assert(std::abs(normal.length() - 1.0f) < TOLERANCE);
  float shift = sinf(PI_F / SHIFT_DIV);
  return {
    {(v + (normal * shift)).normalize()},
    {(v + (-normal * shift)).normalize()},
  };
}
