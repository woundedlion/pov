/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_TRANSFORMERS_H_
#define HOLOSPHERE_CORE_TRANSFORMERS_H_

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
template <int W, typename ParamsT, typename AnimT,
          Vector (*TransformFunc)(const Vector &, const ParamsT &),
          int CAPACITY = 32>
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
   * @return Pointer to the spawned animation, or nullptr if no slots.
   */
  template <typename... Args> AnimT *spawn(int in_frames, Args &&...args) {
    // Linear scan for a free slot
    for (auto &e : entities) {
      if (!e.active) {
        e.params = params;
        e.active = true;
        // Create animation with reference to the stable entity params
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        bool repeats = anim.repeats();
        return timeline.add_get(in_frames,
                                std::move(anim).then([&e, repeats]() {
                                  if (!repeats)
                                    e.active = false;
                                }));
      }
    }
    // If we get here, no free slots. Drop the spawn (safe failure).
    return nullptr;
  }

  /**
   * @brief Prepares per-frame cached state for all active entities.
   * Call once per frame before any transform calls.
   */
  void prepare_frame() {
    for (auto &e : entities) {
      if (e.active) {
        if constexpr (requires { e.params.prepare_thresholds(); }) {
          e.params.prepare_thresholds();
        }
      }
    }
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

  Vector operator()(const Vector &v) const { return transform(v); }
};

/**
 * @brief A transformer adapter for an Orientation object.
 */
template <int W> struct OrientTransformer {
  const Orientation<W> &orientation;

  explicit OrientTransformer(const Orientation<W> &ori) : orientation(ori) {}

  Vector transform(const Vector &v) const { return orientation.orient(v); }

  Vector operator()(const Vector &v) const { return transform(v); }
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
  return inv_gnomonic(w, (v.y >= 0 ? 1.0f : -1.0f));
}

/**
 * @brief Applies 3D noise distortion to a vector.
 * @param v The vector to transform.
 * @param params The noise parameters.
 * @return The distorted vector.
 */
inline Vector ripple_transform(const Vector &v, const RippleParams &params) {
  // fast reject
  float cos_d = dot(v, params.center);
  if (cos_d > params.cos_threshold_min || cos_d < params.cos_threshold_max) {
    return v;
  }

  float d = acosf(hs::clamp(cos_d, -1.0f, 1.0f));
  float dist_from_peak = d - params.phase;

  float half_width = params.thickness * 0.5f;
  if (half_width < 0.001f)
    half_width = 0.001f;

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

  return v;
}
inline Vector noise_transform(const Vector &v, const NoiseParams &params) {
  if (params.amplitude <= 0.001f)
    return v;

  float scale = params.scale;
  float time_val = params.time * params.speed;

  // Sample 3D noise field
  float nx =
      params.noise.GetNoise(v.x * scale, v.y * scale, v.z * scale + time_val);
  float ny = params.noise.GetNoise(v.x * scale + 100.0f, v.y * scale + 100.0f,
                                   v.z * scale + time_val);
  float nz = params.noise.GetNoise(v.x * scale + 200.0f, v.y * scale + 200.0f,
                                   v.z * scale + time_val);

  Vector raw_noise = Vector(nx, ny, nz) * (params.amplitude * 0.05f);

  // Project noise onto the tangent plane
  float inward_pull = dot(raw_noise, v);
  Vector surface_distortion = raw_noise - (v * inward_pull);

  // Soft-cap the slide distance to prevent cross-hemisphere grabs
  constexpr float max_slide = 0.5f;
  float sd_len_sq = dot(surface_distortion, surface_distortion);
  if (sd_len_sq > max_slide * max_slide) {
    surface_distortion = surface_distortion * (max_slide / sqrtf(sd_len_sq));
  }

  return (v + surface_distortion).normalized();
}

/**
 * @brief Result of a stereographic-space noise warp.
 */
struct StereoWarpResult {
  Complex coords;       ///< Displaced coordinate
  float displacement;   ///< Magnitude of the displacement vector
};

/**
 * @brief Applies noise-based displacement in stereographic space.
 * Displacement is attenuated near the projection pole to prevent
 * singularity blowup. Returns both the warped coordinate and the
 * scalar displacement for downstream effects (e.g., hue shifting).
 *
 * @param z         Stereographic coordinate to warp.
 * @param r_sq      Pre-computed |z|² (z.re² + z.im²).
 * @param noise     FastNoiseLite instance for sampling.
 * @param scale     Noise frequency scale.
 * @param strength  Maximum displacement amplitude.
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @param time      Noise time coordinate (for animation).
 */
inline StereoWarpResult stereo_noise_warp(const Complex &z, float r_sq,
                                          const FastNoiseLite &noise,
                                          float scale, float strength,
                                          float pole_fade, float time) {
  float atten = 1.0f / (1.0f + (r_sq / (pole_fade * pole_fade)));
  float s = strength * atten;
  float dx = noise.GetNoise(z.re * scale, z.im * scale, time) * s;
  float dy =
      noise.GetNoise(z.re * scale + 100.0f, z.im * scale + 100.0f, time) * s;
  return {Complex(z.re + dx, z.im + dy), sqrtf(dx * dx + dy * dy)};
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
#endif // HOLOSPHERE_CORE_TRANSFORMERS_H_
