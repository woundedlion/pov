/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file noise_field.h
 * @brief Sphere-domain noise sampling: field specs and their cache keys, the
 *        looped sample coordinates, the octave/basis kernels, the scalar,
 *        direct and curl tangent-field kernels, and the sphere exponential
 *        map and parallel transport that advect a point along one.
 */

#include <array>

#include "math/geometry.h"
#include "math/3dmath.h"
#include "vendor/FastNoiseLite.h"

/** @brief Which coordinate the caller hands the lattice. */
enum class NoiseDomain : uint8_t {
  /** A direction on the sphere. */
  SPHERE_3D,
  /** A point in a projection's plane. */
  PROJECTED_2D
};

/** @brief Octave structure of one scalar sample. */
enum class NoiseBasis : uint8_t {
  /** One generator sample. */
  SIMPLEX,
  /** Three octaves at 1x, 2x, 4x, weighted 1, 1/2, 1/4. */
  FBM3,
  /** The same three octaves folded through 1 - |n|. */
  RIDGED3
};

/** @brief How a sample's channels are generated. */
enum class NoiseChannelLayout : uint8_t {
  /** One scalar channel. */
  SCALAR_V1,
  /** Three scalar channels read at NOISE_CHANNEL_OFFSETS. */
  DIRECT_V1,
  /** Gradient by finite difference on the tetrahedral stencil. */
  CURL_V1,
  /** Three channels from the generator's one vector-noise call. */
  DIRECT_VECTOR_V2,
  /** Gradient from the generator's one analytic-gradient call. */
  CURL_ANALYTIC_V2
};

/** @brief Everything one noise field is sampled with. */
struct NoiseFieldSpec {
  /** Coordinate the caller hands the lattice. */
  NoiseDomain domain = NoiseDomain::SPHERE_3D;
  /** Octave structure of a scalar sample. */
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  /** Generator seed. */
  int32_t seed = 1337;
  /** Lattice units per unit of the caller's coordinate. */
  float scale = 1.0f;
  /** Turns of the time loop per unit time. */
  float rate = 0.0f;
  /** Position on the time loop, in turns. */
  float phase = 0.0f;
  /** How a sample's channels are generated. */
  NoiseChannelLayout channel_layout = NoiseChannelLayout::SCALAR_V1;
  /** Which octave stack the consumer walks. */
  uint8_t octave_layout = 1;
  /** Which time-loop construction the consumer applies. */
  uint8_t loop_layout = 1;
  /** Which gradient stencil the consumer walks; 0 when it takes none. */
  uint8_t stencil_layout = 1;

  constexpr bool operator==(const NoiseFieldSpec &) const = default;
};

/**
 * @brief Identity of a shared generator: a spec without its per-frame drive.
 * @details Only seed, noise_type and generator_frequency configure the
 * generator; the rest of the key separates consumers that walk one generator
 * differently, so each keeps its own cached resource.
 */
struct NoiseFieldKey {
  /** Coordinate the caller hands the lattice. */
  NoiseDomain domain = NoiseDomain::SPHERE_3D;
  /** Octave structure of a scalar sample. */
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  /** Generator seed. */
  int32_t seed = 1337;
  /** How a sample's channels are generated. */
  NoiseChannelLayout channel_layout = NoiseChannelLayout::SCALAR_V1;
  /** Which octave stack the consumer walks. */
  uint8_t octave_layout = 1;
  /** Which time-loop construction the consumer applies. */
  uint8_t loop_layout = 1;
  /** Which gradient stencil the consumer walks; 0 when it takes none. */
  uint8_t stencil_layout = 1;
  /** FastNoiseLite::NoiseType the generator is set to. */
  uint8_t noise_type = FastNoiseLite::NoiseType_OpenSimplex2;
  /** Frequency the generator is set to; the spec's scale is applied by the
   *  caller instead, so one generator serves every scale. */
  float generator_frequency = 1.0f;

  HS_COLD_MEMBER constexpr bool
  operator==(const NoiseFieldKey &) const = default;
};

/**
 * @brief Projects a spec onto the generator identity it needs.
 * @param spec Field being sampled.
 * @return The key two specs share exactly when one generator serves both.
 */
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

/** Lattice offsets that decorrelate the three channels of a DIRECT_V1 vector. */
constexpr std::array<Vector, 3> NOISE_CHANNEL_OFFSETS = {
    Vector(0.0f, 0.0f, 0.0f), Vector(17.0f, -29.0f, 43.0f),
    Vector(-47.0f, 11.0f, -23.0f)};
/** Second set of channel offsets, subtracted to recentre a ridged basis. */
constexpr std::array<Vector, 3> NOISE_RIDGED_OFFSETS = {
    Vector(-73.271f, 19.119f, 5.0f), Vector(61.731f, 89.417f, -7.0f),
    Vector(13.571f, -59.213f, 97.331f)};
/** Radius of the lattice circle phase traverses, in lattice units. */
constexpr float NOISE_LOOP_RADIUS = 32.0f;
/** Arm length of the tetrahedral gradient stencil, in lattice units. */
constexpr float NOISE_STENCIL_RADIUS = 1.0f / 64.0f;

/**
 * @brief Lattice coordinate for a sphere direction at a point on the loop.
 * @param v Direction on the sphere.
 * @param scale Lattice units per unit direction.
 * @param phase Position on the time loop, in turns.
 * @return The sampling coordinate.
 * @details Time traverses a NOISE_LOOP_RADIUS circle in the lattice's xy
 * plane, so the field returns to itself after one turn of phase.
 */
inline Vector noise_sphere_coordinate(const Vector &v, float scale,
                                      float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  return scale * v + Vector(NOISE_LOOP_RADIUS * cosf(angle),
                            NOISE_LOOP_RADIUS * sinf(angle), 0.0f);
}

/**
 * @brief noise_sphere_coordinate with the loop point precomputed.
 * @param v Direction on the sphere.
 * @param scale Lattice units per unit direction.
 * @param loop_offset Point on the loop, hoisted out of a per-pixel walk.
 * @return The sampling coordinate.
 */
HS_FLASH_INLINE inline Vector
noise_sphere_coordinate(const Vector &v, float scale,
                        const Vector &loop_offset) {
  return scale * v + loop_offset;
}

/**
 * @brief Lattice coordinate for a plane point at a point on the loop.
 * @param p Point in a projection's plane.
 * @param scale Lattice units per plane unit.
 * @param phase Position on the time loop, in turns.
 * @return The sampling coordinate.
 * @details The loop circle spans z and the plane's x = y diagonal, so neither
 * plane axis carries the whole of it.
 */
HS_FLASH_INLINE inline Vector
noise_projected_coordinate(const Complex &p, float scale, float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  const float diagonal = NOISE_LOOP_RADIUS * sinf(angle) * 0.7071067811865475f;
  return Vector(scale * p.re + diagonal, scale * p.im + diagonal,
                NOISE_LOOP_RADIUS * cosf(angle));
}

/**
 * @brief One scalar sample of the basis's octave stack.
 * @param noise Prepared generator.
 * @param basis Octave structure to apply.
 * @param q Lattice coordinate.
 * @return Noise value in [-1, 1].
 * @details Costs one generator sample for SIMPLEX and three for the others.
 */
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

/**
 * @brief One component of a DIRECT_V1 vector sample.
 * @param noise Prepared generator.
 * @param basis Octave structure to apply.
 * @param q Lattice coordinate.
 * @param channel Component index in [0, 3), selecting the lattice offset.
 * @return Noise value in [-1, 1].
 * @details RIDGED3 is one-sided, so its component is the difference of two
 * ridged stacks read at independent offsets — six generator samples against
 * the other bases' one or three.
 */
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

/**
 * @brief Three correlated noise components from one generator call.
 * @param noise Prepared generator.
 * @param q Lattice coordinate.
 * @return The displacement the generator applies to @p q.
 * @details The DIRECT_VECTOR_V2 path: one simplex vector-noise call in place
 * of the three or six scalar samples DIRECT_V1 costs. Simplex only.
 */
__attribute__((always_inline)) inline Vector
sample_simplex_vector(const FastNoiseLite &noise, const Vector &q) {
  float x = q.x;
  float y = q.y;
  float z = q.z;
  noise.GetVectorNoiseSingle(x, y, z);
  return Vector(x - q.x, y - q.y, z - q.z);
}

/**
 * @brief Tangent field from a 3D noise vector projected onto the sphere.
 * @param noise Prepared generator.
 * @param basis Octave structure to apply.
 * @param q Lattice coordinate.
 * @param v Unit point the tangent is taken at.
 * @param direction_cos Cosine of the in-plane rotation applied to the tangent.
 * @param direction_sin Sine of that rotation.
 * @return A tangent at @p v of length at most 1.
 * @details Takes the DIRECT_VECTOR_V2 path for SIMPLEX and DIRECT_V1
 * otherwise. Only lengths above 1 are rescaled, so the field keeps its own
 * magnitude where the noise is quiet.
 */
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

/**
 * @brief sample_direct_tangent with the rotation given in turns.
 * @param noise Prepared generator.
 * @param basis Octave structure to apply.
 * @param q Lattice coordinate.
 * @param v Unit point the tangent is taken at.
 * @param direction In-plane rotation applied to the tangent, in turns.
 * @return A tangent at @p v of length at most 1.
 */
HS_FLASH_INLINE inline Vector
sample_direct_tangent(const FastNoiseLite &noise, NoiseBasis basis,
                      const Vector &q, const Vector &v, float direction) {
  const float angle = TWO_PI_F * direction;
  return sample_direct_tangent(noise, basis, q, v, cosf(angle), sinf(angle));
}

/**
 * @brief sample_direct_tangent's SIMPLEX path with no in-plane rotation.
 * @param noise Prepared generator.
 * @param q Lattice coordinate.
 * @param v Unit point the tangent is taken at.
 * @return A tangent at @p v of length at most 1.
 */
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

/**
 * @brief Gradient of a scalar field from a four-probe tetrahedral stencil.
 * @tparam Sample Callable taking a Vector and returning the field value.
 * @param q Lattice coordinate the gradient is taken at.
 * @param sample The scalar field.
 * @return The estimated gradient.
 * @details Four probes on a regular tetrahedron of arm NOISE_STENCIL_RADIUS,
 * which is the fewest that spans three dimensions. The arms sum to zero, so the
 * constant term cancels, but they are not antipodal pairs and the curvature
 * term survives — the estimate is first-order in the arm, where a paired
 * central difference would be second-order. Each probe costs whatever
 * @p sample costs, so a three-octave basis pays twelve generator samples.
 */
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

/**
 * @brief Divergence-free tangent field from the generator's analytic gradient.
 * @param noise Prepared generator.
 * @param q Lattice coordinate.
 * @param v Unit point the tangent is taken at.
 * @return A tangent at @p v of length at most 1.
 * @details The CURL_ANALYTIC_V2 path: one generator call, against the twelve
 * samples the stencil costs. Simplex only.
 */
HS_FLASH_INLINE inline Vector
sample_simplex_curl_tangent(const FastNoiseLite &noise, const Vector &q,
                            const Vector &v) {
  Vector gradient;
  noise.GetNoiseGradientSingle(q.x, q.y, q.z, gradient.x, gradient.y,
                               gradient.z);
  const Vector tangent_gradient = gradient - dot(gradient, v) * v;
  Vector u = cross(v, tangent_gradient);
  const float length = u.length();
  if (length > 1.0f)
    u /= length;
  return u;
}

/**
 * @brief Divergence-free tangent field, by whichever gradient the basis has.
 * @param noise Prepared generator.
 * @param basis Octave structure to apply.
 * @param q Lattice coordinate.
 * @param v Unit point the tangent is taken at.
 * @return A tangent at @p v of length at most 1.
 * @details SIMPLEX takes the analytic CURL_ANALYTIC_V2 path at one generator
 * sample; the three-octave bases take the CURL_V1 stencil at twelve.
 */
HS_FLASH_INLINE inline Vector sample_curl_tangent(const FastNoiseLite &noise,
                                                  NoiseBasis basis,
                                                  const Vector &q,
                                                  const Vector &v) {
  if (basis == NoiseBasis::SIMPLEX)
    return sample_simplex_curl_tangent(noise, q, v);
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

/**
 * @brief Advances a unit point along the great circle its tangent selects.
 * @param v Unit point on the sphere.
 * @param tangent Tangent at @p v; its length is the arc travelled in radians.
 * @return The unit point reached after that arc.
 * @details Exact at any arc length; costs a length, a sinf and a cosf.
 */
inline Vector sphere_exp_map(const Vector &v, const Vector &tangent) {
  const float distance = tangent.length();
  if (distance == 0.0f)
    return v;
  return cosf(distance) * v + (sinf(distance) / distance) * tangent;
}

/**
 * @brief sphere_exp_map for short arcs, with the transcendentals unrolled.
 * @param v Unit point on the sphere.
 * @param tangent Tangent at @p v; its length is the arc travelled in radians.
 * @return The unit point reached after that arc.
 * @details cos and sin/x carried to their degree-6 Taylor terms, so the arc
 *   must stay within about half a radian: truncation error is under ~1e-7
 *   there, and grows as the eighth power of the arc beyond it. Takes the
 *   squared length, so it needs no sqrt and no branch at a zero tangent.
 */
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

/** Conditioning bound on |cross(from, to)|^2 for transport_tangent; below it
 *  the great circle is ill-determined. */
inline constexpr float MIN_TRANSPORT_CROSS_SQ = 1e-4f;

/**
 * @brief Parallel-transports a tangent vector along the great circle from
 *   @p from to @p to.
 * @param from Unit source point.
 * @param to Unit destination point.
 * @param tangent_at_from Tangent vector at @p from.
 * @return The transported vector, tangent at @p to.
 * @details Traps once the pair is inside MIN_TRANSPORT_CROSS_SQ of antipodal
 *   rather than only at the exact antipode: any non-tangent component of
 *   @p tangent_at_from is amplified by 2/|cross(from, to)| there, capped at
 *   200x by the bound. Gated on |cross|^2, which stays accurate where 1 + dot
 *   cancels; the dot > 0 half-space short-circuits before the cross product.
 */
inline Vector transport_tangent(const Vector &from, const Vector &to,
                                const Vector &tangent_at_from) {
  const float denominator = 1.0f + dot(from, to);
  HS_CHECK(denominator > 1.0f ||
               dot(cross(from, to), cross(from, to)) > MIN_TRANSPORT_CROSS_SQ,
           "tangent transport is within the antipodal conditioning bound");
  return tangent_at_from -
         (dot(tangent_at_from, to) / denominator) * (from + to);
}
