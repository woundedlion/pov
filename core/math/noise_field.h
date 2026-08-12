/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/math/3dmath.h"
#include "core/vendor/FastNoiseLite.h"

enum class NoiseDomain : uint8_t { SPHERE_3D, PROJECTED_2D };
enum class NoiseBasis : uint8_t { SIMPLEX, FBM3, RIDGED3 };
enum class NoiseChannelLayout : uint8_t {
  SCALAR_V1,
  DIRECT_V1,
  CURL_V1,
  DIRECT_VECTOR_V2
};

struct NoiseFieldSpec {
  NoiseDomain domain = NoiseDomain::SPHERE_3D;
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  int32_t seed = 1337;
  float scale = 1.0f;
  float rate = 0.0f;
  float phase = 0.0f;
  NoiseChannelLayout channel_layout = NoiseChannelLayout::SCALAR_V1;
  uint8_t octave_layout = 1;
  uint8_t loop_layout = 1;
  uint8_t stencil_layout = 1;

  constexpr bool operator==(const NoiseFieldSpec &) const = default;
};

struct NoiseFieldKey {
  NoiseDomain domain = NoiseDomain::SPHERE_3D;
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  int32_t seed = 1337;
  NoiseChannelLayout channel_layout = NoiseChannelLayout::SCALAR_V1;
  uint8_t octave_layout = 1;
  uint8_t loop_layout = 1;
  uint8_t stencil_layout = 1;
  uint8_t noise_type = FastNoiseLite::NoiseType_OpenSimplex2;
  float generator_frequency = 1.0f;

  HS_COLD_MEMBER constexpr bool
  operator==(const NoiseFieldKey &) const = default;
};

constexpr NoiseFieldKey noise_field_key(const NoiseFieldSpec &spec) {
  return {spec.domain,
          spec.basis,
          spec.seed,
          spec.channel_layout,
          spec.octave_layout,
          spec.loop_layout,
          spec.stencil_layout,
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

constexpr std::array<Vector, 3> NOISE_CHANNEL_OFFSETS = {
    Vector(0.0f, 0.0f, 0.0f), Vector(17.0f, -29.0f, 43.0f),
    Vector(-47.0f, 11.0f, -23.0f)};
constexpr std::array<Vector, 3> NOISE_RIDGED_OFFSETS = {
    Vector(-73.271f, 19.119f, 5.0f), Vector(61.731f, 89.417f, -7.0f),
    Vector(13.571f, -59.213f, 97.331f)};
constexpr float NOISE_LOOP_RADIUS = 32.0f;
constexpr float NOISE_STENCIL_RADIUS = 1.0f / 64.0f;

inline Vector noise_sphere_coordinate(const Vector &v, float scale,
                                      float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  return scale * v + Vector(NOISE_LOOP_RADIUS * cosf(angle),
                            NOISE_LOOP_RADIUS * sinf(angle), 0.0f);
}

HS_FLASH_INLINE inline Vector
noise_sphere_coordinate(const Vector &v, float scale,
                        const Vector &loop_offset) {
  return scale * v + loop_offset;
}

HS_FLASH_INLINE inline Vector
noise_projected_coordinate(const Complex &p, float scale, float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  const float diagonal = NOISE_LOOP_RADIUS * sinf(angle) * 0.7071067811865475f;
  return Vector(scale * p.re + diagonal, scale * p.im + diagonal,
                NOISE_LOOP_RADIUS * cosf(angle));
}

HS_FLASH_INLINE inline float sample_noise_octaves(const FastNoiseLite &noise,
                                                  NoiseBasis basis,
                                                  const Vector &q) {
  const float first = noise.GetNoiseSingle(q.x, q.y, q.z);
  if (basis == NoiseBasis::SIMPLEX)
    return first;
  const float second = noise.GetNoiseSingle(2.0f * q.x, 2.0f * q.y, 2.0f * q.z);
  const float third = noise.GetNoiseSingle(4.0f * q.x, 4.0f * q.y, 4.0f * q.z);
  if (basis == NoiseBasis::FBM3)
    return (first + 0.5f * second + 0.25f * third) / 1.75f;
  const float ridge = (1.0f - std::abs(first)) +
                      0.5f * (1.0f - std::abs(second)) +
                      0.25f * (1.0f - std::abs(third));
  return 2.0f * (ridge / 1.75f) - 1.0f;
}

HS_FLASH_INLINE inline float
sample_noise_vector_channel(const FastNoiseLite &noise, NoiseBasis basis,
                            const Vector &q, size_t channel) {
  if (basis != NoiseBasis::RIDGED3)
    return sample_noise_octaves(noise, basis,
                                q + NOISE_CHANNEL_OFFSETS[channel]);
  const float a =
      sample_noise_octaves(noise, basis, q + NOISE_CHANNEL_OFFSETS[channel]);
  const float b =
      sample_noise_octaves(noise, basis, q + NOISE_RIDGED_OFFSETS[channel]);
  return 0.5f * (a - b);
}

__attribute__((always_inline)) inline Vector
sample_simplex_vector(const FastNoiseLite &noise, const Vector &q) {
  float x = q.x;
  float y = q.y;
  float z = q.z;
  noise.GetVectorNoiseSingle(x, y, z);
  return Vector(x - q.x, y - q.y, z - q.z);
}

HS_FLASH_INLINE inline Vector
sample_direct_tangent(const FastNoiseLite &noise, NoiseBasis basis,
                      const Vector &q, const Vector &v, float direction) {
  Vector u;
  if (basis == NoiseBasis::SIMPLEX)
    u = sample_simplex_vector(noise, q);
  else
    u = Vector(sample_noise_vector_channel(noise, basis, q, 0),
               sample_noise_vector_channel(noise, basis, q, 1),
               sample_noise_vector_channel(noise, basis, q, 2));
  u -= dot(u, v) * v;
  const float length_sq = u.x * u.x + u.y * u.y + u.z * u.z;
  if (length_sq > 1.0f)
    u *= fast_rsqrt(length_sq);
  const float angle = TWO_PI_F * direction;
  return cosf(angle) * u + sinf(angle) * cross(v, u);
}

HS_FLASH_INLINE inline Vector
sample_direct_tangent(const FastNoiseLite &noise, NoiseBasis basis,
                      const Vector &q, const Vector &v, float direction_cos,
                      float direction_sin) {
  Vector u;
  if (basis == NoiseBasis::SIMPLEX)
    u = sample_simplex_vector(noise, q);
  else
    u = Vector(sample_noise_vector_channel(noise, basis, q, 0),
               sample_noise_vector_channel(noise, basis, q, 1),
               sample_noise_vector_channel(noise, basis, q, 2));
  u -= dot(u, v) * v;
  const float length_sq = u.x * u.x + u.y * u.y + u.z * u.z;
  if (length_sq > 1.0f)
    u *= fast_rsqrt(length_sq);
  return direction_cos * u + direction_sin * cross(v, u);
}

__attribute__((always_inline)) inline Vector
sample_direct_simplex_tangent(const FastNoiseLite &noise, const Vector &q,
                              const Vector &v) {
  const Vector u0 = sample_simplex_vector(noise, q);
  Vector u = u0 - dot(u0, v) * v;
  const float length_sq = dot(u, u);
  if (length_sq > 1.0f)
    u *= fast_rsqrt(length_sq);
  return u;
}

template <typename Sample>
inline Vector tetrahedral_gradient(const Vector &q, Sample sample) {
  constexpr float INV_SQRT_3 = 0.5773502691896258f;
  constexpr std::array<Vector, 4> DIRECTIONS = {
      Vector(INV_SQRT_3, INV_SQRT_3, INV_SQRT_3),
      Vector(INV_SQRT_3, -INV_SQRT_3, -INV_SQRT_3),
      Vector(-INV_SQRT_3, INV_SQRT_3, -INV_SQRT_3),
      Vector(-INV_SQRT_3, -INV_SQRT_3, INV_SQRT_3)};
  Vector gradient;
  for (const Vector &direction : DIRECTIONS)
    gradient += sample(q + NOISE_STENCIL_RADIUS * direction) * direction;
  return (3.0f / (4.0f * NOISE_STENCIL_RADIUS)) * gradient;
}

HS_FLASH_INLINE inline Vector sample_curl_tangent(const FastNoiseLite &noise,
                                                  NoiseBasis basis,
                                                  const Vector &q,
                                                  const Vector &v) {
  const Vector gradient = tetrahedral_gradient(q, [&](const Vector &point) {
    return sample_noise_octaves(noise, basis, point);
  });
  const Vector tangent_gradient = gradient - dot(gradient, v) * v;
  Vector u = cross(v, tangent_gradient);
  const float length = u.length();
  if (length > 1.0f)
    u /= length;
  return u;
}

inline Vector sphere_exp_map(const Vector &v, const Vector &tangent) {
  const float distance = tangent.length();
  if (distance == 0.0f)
    return v;
  return cosf(distance) * v + (sinf(distance) / distance) * tangent;
}

HS_FLASH_INLINE inline Vector
sphere_exp_map_half_radian(const Vector &v, const Vector &tangent) {
  const float distance_sq = dot(tangent, tangent);
  const float cosine =
      1.0f + distance_sq *
                 (-0.5f + distance_sq *
                              (1.0f / 24.0f + distance_sq * (-1.0f / 720.0f)));
  const float sinc =
      1.0f + distance_sq * (-1.0f / 6.0f +
                            distance_sq * (1.0f / 120.0f +
                                           distance_sq * (-1.0f / 5040.0f)));
  return cosine * v + sinc * tangent;
}

inline Vector transport_tangent(const Vector &from, const Vector &to,
                                const Vector &tangent_at_from) {
  const float denominator = 1.0f + dot(from, to);
  HS_CHECK(denominator > 0.0f, "tangent transport crosses antipodes");
  return tangent_at_from -
         (dot(tangent_at_from, to) / denominator) * (from + to);
}
