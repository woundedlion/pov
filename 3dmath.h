#pragma once

#include <cmath>
#include <cfloat>
#include <limits>
#include <cassert>

static constexpr double PHI = 1.6180339887498948482045868343656;
static constexpr double G = 1 / PHI;
static constexpr double TOLERANCE = 0.0001f;

struct Vector;
struct Spherical {
  Spherical(double theta, double phi) : theta(theta), phi(phi) {}
  Spherical(const Vector& v);

  double theta; // azimuthal angle
  double phi; // polar angle
};

struct Vector {
  constexpr Vector() {}
  constexpr Vector(double i, double j, double k) : i(i), j(j), k(k) {}
  constexpr Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  Vector(const Spherical& s) :
    i(sin(s.phi)* cos(s.theta)),
    j(sin(s.phi)* sin(s.theta)),
    k(cos(s.phi))
  {}
  Vector(double theta, double phi) : Vector(Spherical(theta, phi)) {}

  constexpr Vector& operator=(const Vector& v) {
    i = v.i; 
    j = v.j;
    k = v.k;
    return *this;
  }

  bool operator==(const Vector& v) const {
    return fabs(i - v.i) <= TOLERANCE
      && fabs(j - v.j) <= TOLERANCE
      && fabs(k - v.k) <= TOLERANCE;
  }

  bool operator!=(const Vector& v) const {
    return !(*this == v);
  }

  constexpr Vector operator-() const {
    return Vector(-i, -j, -k);
  }

  constexpr double length() const { return sqrt(i * i + j * j + k * k); }
  
  Vector& normalize() {
    double m = length();
    if (m < std::numeric_limits<double>::epsilon()) {
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

  double i = 0;
  double j = 0;
  double k = 0;
};

Spherical::Spherical(const Vector& v)
{
  Vector n(v);
  n.normalize();
  theta = atan2(n.j, n.i);
  phi = acos(std::clamp(n.k, -1.0, 1.0));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
constexpr Quaternion operator/(const Quaternion& q, double s);
constexpr Vector operator/(const Vector& v, double s);

struct Quaternion {
  constexpr Quaternion() {}
  constexpr Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  constexpr  Quaternion(double a, double b, double c, double d) : r(a), v(b, c, d) {}
  constexpr Quaternion(double a, const Vector& v) : r(a), v(v) {}
  
  Quaternion& operator=(const Quaternion& q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  bool operator==(const Quaternion& q) const {
    return fabs(q.r - r) < TOLERANCE
      && q.v == v;
  }

  void operator*=(const Quaternion& q) {
    *this = *this * q;
  }

  constexpr double squared_magnitude() const {
    return r * r + v.i * v.i + v.j * v.j + v.k * v.k;
  }

  Quaternion inverse() const {
    double sq_mag = squared_magnitude();
    return Quaternion(r, -v) / sq_mag;
  }

  constexpr Quaternion operator-() const {
    return Quaternion(-r, -v);
  }

  constexpr double magnitude() const {
    return sqrt(r * r + v.i * v.i + v.j * v.j + v.k * v.k);
  }

  Quaternion& normalize() {
    auto m = magnitude();
    if (m <= std::numeric_limits<double>::epsilon()) {
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

  double r = 1;
  Vector v;
};

Vector rotate(const Vector& v, const Quaternion& q) {
  assert(fabs(q.magnitude() - 1) < TOLERANCE);
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

constexpr Vector operator*(const Vector& v, double s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

constexpr Vector operator*(double s, const Vector& v) {
  return v * s;
}

constexpr Vector operator/(const Vector& v, double s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

constexpr double dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

constexpr Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

constexpr double distance_between(const Vector& a, const Vector& b) {
  double di = b.i - a.i;
  double dj = b.j - a.j;
  double dk = b.k - a.k;
  return sqrt(di * di + dj * dj + dk * dk);
}

constexpr double angle_between(const Vector& v1, const Vector& v2) {
  double len_product = v1.length() * v2.length();
  if (len_product <= std::numeric_limits<double>::epsilon()) {
    assert(false);
  }
  double d = dot(v1, v2) / len_product;
  return acos(std::clamp(d, -1.0, 1.0));
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

constexpr Quaternion operator*(const Quaternion& q, double s) {
  return Quaternion(q.r * s, q.v * s);
}

constexpr Quaternion operator*(double s, const Quaternion& q) {
  return q * s;
}

constexpr Quaternion operator/(const Quaternion& q, double s) {
  assert(fabs(s) > std::numeric_limits<double>::epsilon());
  return Quaternion(q.r / s, q.v / s);
}

constexpr double dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

constexpr double angle_between(const Quaternion& q1, const Quaternion& q2) {
  assert(fabs(q1.magnitude() - 1.0) < TOLERANCE);
  assert(fabs(q2.magnitude() - 1.0) < TOLERANCE);
  return acos(std::clamp(dot(q1, q2), -1.0, 1.0));
}

Quaternion make_rotation(const Vector& axis, double theta) {
  assert(fabs(axis.length() - 1.0) < TOLERANCE);
  return Quaternion(cos(theta / 2), sin(theta / 2) * axis).normalize();
}

Quaternion slerp(const Quaternion& q1, const Quaternion& q2, double t, bool long_way = false) {
  assert(fabs(q1.magnitude() - 1.0) < TOLERANCE);
  assert(fabs(q2.magnitude() - 1.0) < TOLERANCE);
  double d = dot(q1, q2);
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

  double theta = acos(std::clamp(d, -1.0, 1.0));
  double sin_theta = sin(theta);
  double s1 = sin((1 - t) * theta) / sin_theta;
  double s2 = sin(t * theta) / sin_theta;
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
  assert(fabs(u.length() - 1.0) < TOLERANCE);
  assert(fabs(v.length() - 1.0) < TOLERANCE);
  assert(fabs(normal.length() - 1.0) < TOLERANCE);
  Vector w = cross(v, u).normalize();
  Vector i1 = cross(w, normal).normalize();
  Vector i2 = cross(normal, w).normalize();

  double a1 = angle_between(u, v);
  double a2 = angle_between(i1, u);
  double a3 = angle_between(i1, v);
  if (fabs(a2 + a3 - a1) < .0001) {
    return i1;
  }

  a1 = angle_between(u, v);
  a2 = angle_between(i2, u);
  a3 = angle_between(i2, v);
  if (fabs(a2 + a3 - a1) < TOLERANCE) {
    return i2;
  }

  assert(false);
  return Vector(0, 0, 0);
}

std::vector<Vector> split_point(const Vector& v, const Vector& normal) {
  assert(fabs(v.length() - 1.0) < TOLERANCE);
  assert(fabs(normal.length() - 1.0) < TOLERANCE);
  double shift = sin(PI / MAX_W);
  return {
    {(v + (normal * shift)).normalize()},
    {(v + (-normal * shift)).normalize()},
  };
}

Vector lissajous(double m1, double m2, double a, double t) {

  Vector v(
    sin(m2 * t) * cos(m1 * t - a * PI),
    sin(m2 * t) * sin(m1 * t - a * PI),
    cos(m2 * t)
  );
  return v.normalize();
}