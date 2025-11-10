#pragma once

#include <cmath>
#include <cfloat>
#include <limits>
#include <cassert>

static constexpr float PI_F = 3.1415926535f;
static constexpr float_t PHI = 1.6180339887498948482045868343656f;
static constexpr float_t G = 1 / PHI;
struct Vector;
struct Spherical {
  Spherical(float_t theta, float_t phi) : theta(theta), phi(phi) {}
  Spherical(const Vector& v);

  float_t theta; // azimuthal angle
  float_t phi; // polar angle
};

struct Vector {
  Vector() {}
  Vector(float_t i, float_t j, float_t k) : i(i), j(j), k(k) {}
  Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  Vector(const Spherical& s) :
    i(sin(s.phi)* cos(s.theta)),
    j(sin(s.phi)* sin(s.theta)),
    k(cos(s.phi))
  {}
  Vector(float_t theta, float_t phi) : Vector(Spherical(theta, phi)) {}

  Vector& operator=(const Vector& v) { 
    i = v.i; 
    j = v.j, 
    k = v.k;
    return *this;
  }

  bool operator==(const Vector& v) const {
    return fabs(i - v.i) <= 0.0001
      && fabs(j - v.j) <= 0.0001
      && fabs(k - v.k) <= 0.0001;
  }

  bool operator!=(const Vector& v) const {
    return !(*this == v);
  }

  Vector operator-() const {
    return Vector(-i, -j, -k);
  }

  float_t length() const { return sqrt(i * i + j * j + k * k); }
  
  Vector& normalize() {
    float_t m = length();
    if (m == 0) {
      Serial.println("Can't normalize a zero vector!");
      assert(false);
    } else {
      i = i / m;
      j = j / m;
      k = k / m;
    }
    return *this;
  }

  float_t i = 0;
  float_t j = 0;
  float_t k = 0;
};

Spherical::Spherical(const Vector& v) :
  theta(atan2(v.j, v.i)),
  phi(acos(v.k))
{}

struct Quaternion;
Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
Vector operator/(const Vector& v, float_t s);

struct Quaternion {
  Quaternion() {}
  Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  Quaternion(float_t a, float_t b, float_t c, float_t d) : r(a), v(b, c, d) {}
  Quaternion(float_t a, const Vector& v) : r(a), v(v) {}
  
  Quaternion& operator=(const Quaternion& q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  bool operator==(const Quaternion& q) {
    return fabs(q.r - r) < 0.0001
      && q.v == v;
  }

  void operator*=(const Quaternion& q) {
    *this = *this * q;
  }

  Quaternion inverse() const { 
    float_t n = (r * r) + (v.i * v.i) + (v.j * v.j) + (v.k * v.k);
    return Quaternion(r / n, (-v) / n);
  }

  Quaternion operator-() const {
    return Quaternion(-r, -v);
  }

  float_t magnitude() const {
    return sqrt(r * r + v.i * v.i + v.j * v.j + v.k * v.k);
  }

  Quaternion& normalize() {
    auto m = magnitude();
    if (m > 0) {
      r = r / m;
      v = v / m;
    }
    return *this;
  }

  float_t r = 1;
  Vector v;
};

Vector rotate(const Vector& v, const Quaternion& q) {
  assert(q.magnitude() - 1 < 0.00001);
  Quaternion p(0, v);
  auto r = q * p * q.inverse();
  return r.v;
}

Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
}

Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
}

Vector operator*(const Vector& v, float_t s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

Vector operator*(float_t s, const Vector& v) {
  return v * s;
}

Vector operator/(const Vector& v, float_t s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

float_t dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

float_t distance_between(const Vector& a, const Vector& b) {
  return sqrt(pow(b.i - a.i, 2)
    + pow(b.j - a.j, 2)
    + pow(b.k - a.k, 2));
}

float_t angle_between(const Vector& v1, const Vector& v2) {
  return acos(std::clamp(dot(v1, v2), -1.0f, 1.0f));
}

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

Quaternion operator*(const Quaternion& q, float_t s) {
  return Quaternion(q.r * s, q.v * s);
}

Quaternion operator*(float_t s, const Quaternion& q) {
  return q * s;
}

float_t dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

float_t angle_between(const Quaternion& q1, const Quaternion& q2) {
  return acos(dot(q1, q2));
}

Quaternion make_rotation(const Vector& axis, float_t theta) {
  return Quaternion(cos(theta / 2), sin(theta / 2) * axis);
}

Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float_t t, bool long_way = false) {
  float_t d = dot(q1, q2);
  Quaternion p(q1);
  Quaternion q(q2);

  if ((long_way && d > 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d > 0.99995) {
    Quaternion r = p + t * (q - p);
    return r.normalize();
  }

  float_t theta = acos(d);
  float_t sin_theta = sin(theta);
  float_t s1 = sin((1 - t) * theta) / sin_theta;
  float_t s2 = sin(t * theta) / sin_theta;
  return (s1 * p) + (s2 * q);
}

bool is_over(const Vector& v, const Vector& normal) {
  return dot(normal, v) >= 0;
}

bool intersects_plane(const Vector& v1, const Vector& v2, const Vector& normal) {
  return (is_over(v1, normal) && !is_over(v2, normal))
    || (!is_over(v1, normal) && is_over(v2, normal));
}

Vector intersection(const Vector& u, const Vector& v, const Vector& normal) {
  Vector w = cross(v, u).normalize();
  Vector i1 = cross(w, normal).normalize();
  Vector i2 = cross(normal, w).normalize();

  float_t a1 = angle_between(u, v);
  float_t a2 = angle_between(i1, u);
  float_t a3 = angle_between(i1, v);
  if (fabs(a2 + a3 - a1) < .0001) {
    return i1;
  }

  a1 = angle_between(u, v);
  a2 = angle_between(i2, u);
  a3 = angle_between(i2, v);
  if (fabs(a2 + a3 - a1) < 0.0001) {
    return i2;
  }

  assert(false);
  return Vector(0, 0, 0);
}

std::vector<Vector> split_point(const Vector& v, const Vector& normal) {
  float_t shift = sin(PI / MAX_W);
  return {
    {(v + (normal * shift)).normalize()},
    {(v + (-normal * shift)).normalize()},
  };
}

Vector lissajous(float_t m1, float_t m2, float_t a, float_t t) {

  Vector v(
    sin(m2 * t) * cos(m1 * t - a * PI),
    sin(m2 * t) * sin(m1 * t - a * PI),
    cos(m2 * t)
  );
  return v;
}