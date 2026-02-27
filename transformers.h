/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"
#include "concepts.h"
#include <array>
#include "FastNoiseLite.h"
#include "animation.h"

using Animation::NoiseParams;
using Animation::RippleParams;

/**
 * @brief A generic manager for state-based geometry transformations.
 * @tparam W Canvas width (for Timeline).
 * @tparam ParamsT The configuration struct (e.g., RippleParams, MobiusParams).
 * @tparam AnimT The animation class (e.g., Animation::Ripple).
 * @tparam TransformFunc The static function to apply the transformation.
 * @tparam CAPACITY Max number of active transformations.
 */
template <int W, typename ParamsT, typename AnimT, auto TransformFunc,
          int CAPACITY = 32>
  requires TransformerFunction<TransformFunc, ParamsT>
class Transformer {
public:
  struct Entity {
    ParamsT params;
    bool active = false;
  };

  ParamsT params;
  std::array<Entity, CAPACITY> entities;
  Timeline<W> &timeline;

  Transformer(Timeline<W> &tl) : timeline(tl) {}

  /**
   * @brief Spawns a new transformation animation.
   * @param in_frames Delay before starting.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   */
  template <typename... Args> void spawn(int in_frames, Args &&...args) {
    // Linear scan for a free slot
    for (auto &e : entities) {
      if (!e.active) {
        e.params = params;
        e.active = true;
        // Create animation with reference to the stable entity params
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        timeline.add(in_frames, anim.then([&e]() { e.active = false; }));
        return;
      }
    }
    // If we get here, no free slots. Drop the spawn (safe failure).
  }

  /**
   * @brief Applies all active transformations to a vector.
   */
  Vector transform(Vector v) const {
    for (const auto &e : entities) {
      if (e.active) {
        v = TransformFunc(v, e.params);
      }
    }
    return v;
  }
};

/**
 * @brief A transformer adapter for an Orientation object.
 */
template <int W> struct OrientTransformer {
  const Orientation<W> &orientation;

  explicit OrientTransformer(const Orientation<W> &ori) : orientation(ori) {}

  Vector transform(const Vector &v) const { return orientation.orient(v); }
};

/**
 * @brief Applies Mobius Transformation to a vector.
 * Wraps the complex math version from 3dmath.h.
 */
inline Vector mobius_transform(const Vector &v, const MobiusParams &params) {
  return inv_stereo(mobius(stereo(v), params));
}

/**
 * @brief Applies Gnomonic Mobius Transformation to a vector.
 * Projects to gnomonic plane, applies Mobius, projects back.
 */
inline Vector gnomonic_mobius_transform(const Vector &v,
                                        const MobiusParams &params) {
  Complex z = gnomonic(v);
  Complex w = mobius(z, params);
  return inv_gnomonic(w, (v.j >= 0 ? 1.0f : -1.0f));
}

/**
 * @brief Parameters for noise transformation.
 */

/**
 * @brief Applies 3D noise distortion to a vector.
 * @param v The vector to transform.
 * @param params The noise parameters.
 * @return The distorted vector.
 */
inline Vector ripple_transform(const Vector &v, const RippleParams &params) {
  float d = angle_between(v, params.center);
  float dist_from_peak = d - params.phase;

  // Defines the width of the single pulse
  float half_width = params.thickness * 0.5f;
  if (half_width < 0.001f)
    half_width = 0.001f; // Prevent division by zero

  if (std::abs(dist_from_peak) < half_width * 2.0f) {

    // Normalize distance (-2 to 2 range covers the whole wavelet)
    float t = (dist_from_peak / half_width) * 2.0f;

    // Ricker Wavelet: (1 - t^2) * e^(-t^2/2)
    float ricker = (1.0f - t * t) * expf(-0.5f * t * t);

    // Distance attenuation (ripples get smaller as they spread)
    float attenuation = expf(-params.decay * d);

    float theta = params.amplitude * ricker * attenuation;

    // Apply rotation
    Vector axis = cross(params.center, v);
    float lenSq = dot(axis, axis);
    if (lenSq > 1e-6f) {
      axis = axis * (1.0f / sqrtf(lenSq));
      Quaternion q = make_rotation(axis, theta);
      return rotate(v, q);
    }
  }

  return v;
}
inline Vector noise_transform(const Vector &v, const NoiseParams &params) {
  if (params.amplitude <= 0.001f)
    return v;

  float scale = 4.0f;
  float time_val = params.time * params.speed;

  // Use member noise object
  float noiseValX =
      params.noise.GetNoise(v.i * scale, v.j * scale, v.k * scale + time_val);
  float noiseValY = params.noise.GetNoise(
      v.i * scale + 100.0f, v.j * scale + 100.0f, v.k * scale + time_val);

  // Simple perturbation
  float dist = 0.2f * params.amplitude; // Arbitrary scale factor

  Vector res = v;
  res.i += noiseValX * dist;
  res.j += noiseValY * dist;

  return res;
}

/**
 * @brief Generates ripples that warp the sphere
 */
template <int W, int CAPACITY>
using RippleTransformer =
    Transformer<W, RippleParams, Animation::Ripple, ripple_transform, CAPACITY>;

/**
 * @brief Performs Mobius warps that return to the identity.
 */
template <int W, int CAPACITY>
using MobiusWarpTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarp, mobius_transform,
                CAPACITY>;

/**
 * @brief Performs circular Mobius warps stay warped throughout sutiable for
 * repeating animations
 */
template <int W, int CAPACITY>
using MobiusWarpCircularTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarpCircular,
                mobius_transform, CAPACITY>;

/**
 * @brief Performs a changing Mobius warp using gnomonic projection.
 */
template <int W, int CAPACITY>
using MobiusWarpGnomonicTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarpEvolving,
                gnomonic_mobius_transform, CAPACITY>;

/**
 * @brief Applies 3D noise distortion to vectors.
 */
template <int W, int CAPACITY>
using NoiseTransformer =
    Transformer<W, NoiseParams, Animation::Noise, noise_transform, CAPACITY>;
