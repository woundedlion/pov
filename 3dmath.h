#pragma once

#include <cmath>
#include <cfloat>
#include <limits>
#include <cassert>

static constexpr double PHI = 1.6180339887498948482045868343656;
static constexpr double G = 1 / PHI;

struct Vector;
struct Spherical {
  Spherical(double theta, double phi) : theta(theta), phi(phi) {}
  Spherical(const Vector& v);

  double theta; // azimuthal angle
  double phi; // polar angle
};

struct Vector {
  Vector() {}
  Vector(double i, double j, double k) : i(i), j(j), k(k) {}
  Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  Vector(const Spherical& s) :
    i(sin(s.phi)* cos(s.theta)),
    j(sin(s.phi)* sin(s.theta)),
    k(cos(s.phi))
  {}
  Vector(double theta, double phi) : Vector(Spherical(theta, phi)) {}

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

  double length() const { return sqrt(i * i + j * j + k * k); }
  
  Vector& normalize() {
    double m = length();
    i = i / m;
    j = j / m;
    k = k / m;
    return *this;
  }

  Vector inverse() {
    return Vector(-i, -j, -k);
  }

  double i = 0;
  double j = 0;
  double k = 0;
};

Spherical::Spherical(const Vector& v) :
  theta(atan2(v.j, v.i)),
  phi(acos(v.k))
{}

struct Quaternion;
Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
Vector operator/(const Vector& v, double s);

struct Quaternion {
  Quaternion() {}
  Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  Quaternion(double a, double b, double c, double d) : r(a), v(b, c, d) {}
  Quaternion(double a, const Vector& v) : r(a), v(v) {}
  
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
    double n = (r * r) + (v.i * v.i) + (v.j * v.j) + (v.k * v.k);
    return Quaternion(r / n, (-v) / n);
  }

  double magnitude() const {
    return sqrt(r * r + v.i * v.i + v.j * v.j + v.k * v.k);
  }

  Quaternion& normalize() {
    auto m = magnitude();
    r = r / m;
    v = v / m;
    return *this;
  }

  double r = 1;
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

Vector operator*(const Vector& v, double s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

Vector operator*(double s, const Vector& v) {
  return v * s;
}

Vector operator/(const Vector& v, double s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

double dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

double distance_between(const Vector& a, const Vector& b) {
  return sqrt(pow(b.i - a.i, 2)
    + pow(b.j - a.j, 2)
    + pow(b.k - a.k, 2);
}

double angle_between(const Vector& v1, const Vector& v2) {
  return acos(std::min(1.0, dot(v1, v2)));
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

Quaternion operator*(const Quaternion& q, double s) {
  return Quaternion(q.r * s, q.v * s);
}

Quaternion operator*(double s, const Quaternion& q) {
  return q * s;
}

double dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

double angle_between(const Quaternion& q1, const Quaternion& q2) {
  return acos(dot(q1, q2));
}

Quaternion make_rotation(const Vector& axis, double theta) {
  return Quaternion(cos(theta / 2), sin(theta / 2) * axis);
}

Quaternion slerp(const Quaternion& q1, const Quaternion& q2, double t) {
  double theta = angle_between(q1, q2);
  return((sin((1 - t) * theta) / sin(theta)) * q1)
    + ((sin(t * theta) / sin(theta)) * q2);
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

  double a1 = angle_between(u, v);
  double a2 = angle_between(i1, u);
  double a3 = angle_between(i1, v);
  if (fabs(a2 + a3 - a1) < .0001) {
    return i1;
  }

  a1 = angle_between(u, v);
  a2 = angle_between(i2, u);
  a3 = angle_between(i2, v);
  if (fabs(a2 + a3 - a1) < 0.0001) {
    return i2;
  }

  return Vector(0, 0, 0);
}

std::vector<Vector> split_point(const Vector& v, const Vector& normal) {
  // TODO: optimize shift
  double shift = sin(PI / 96);
  return {
    {(v + (normal * shift)).normalize()},
    {(v + (-normal * shift)).normalize()},
  };
}

///////////////////////////////////////////////////////////////////////////////

double ease_in_out_bicubic(double t) {
  return t < 0.5 ? 4 * pow(t, 3) : 1 - pow(-2 * t + 2, 3) / 2;
}

double ease_in_out_sin(double t) {
  return -(cos(PI * t) - 1) / 2;
}

double ease_in_sin(double t) {
  return 1 - cos((t * PI) / 2);
}

double ease_out_sin(double t) {
  return sin((t * PI) / 2);
}

double ease_in_cubic(double t) {
  return pow(t, 3);
}

double ease_in_circ(double t) {
  return 1 - sqrt(1 - pow(t, 2));
}

double ease_mid(double t) {
  return t;
}

double ease_out_expo(double t) {
  return t == 1 ? 1 : 1 - pow(2, -10 * t);
}

double ease_out_circ(double t) {
  return sqrt(1 - pow(t - 1, 2));
}

Vector lissajous(double m1, double m2, double a, double t) {

  Vector v(
    sin(m2 * t) * cos(m1 * t - a * PI),
    sin(m2 * t) * sin(m1 * t - a * PI),
    cos(m2 * t)
  );
  return v;
}