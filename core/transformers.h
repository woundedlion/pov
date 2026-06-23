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
 * @tparam ParamsT The configuration struct (e.g., RippleParams, MobiusParams).
 * @tparam AnimT The animation class (e.g., Animation::Ripple).
 * @tparam TransformFunc The static function to apply the transformation.
 * @tparam CAPACITY Max number of active transformations.
 */
template <typename ParamsT, typename AnimT,
          Vector (*TransformFunc)(const Vector &, const ParamsT &),
          int CAPACITY = 32>
class Transformer {
public:
  /**
   * @brief Per-slot storage for one active transformation.
   */
  struct Entity {
    ParamsT params;         /**< Per-entity configuration the transform reads. */
    bool active = false;    /**< Whether this slot currently holds a live animation. */
  };

  ParamsT template_params;              /**< Template params copied into each new entity on spawn. */
  std::array<Entity, CAPACITY> entities; /**< Fixed-capacity pool of transformation slots. */
  Timeline &timeline;                   /**< Timeline that schedules and steps the spawned animations. */

  /**
   * @brief Constructs a transformer bound to a timeline.
   * @param tl Timeline used to schedule spawned animations; retained by reference.
   */
  Transformer(Timeline &tl) : timeline(tl) {}

  /**
   * @brief Spawns a new transformation animation.
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no slots are free.
   * @details pin=false: the returned pointer is transient (used at the call
   * site, not retained across frames). These animations are often finite and
   * are compacted normally; pinning them would trap on routine completion.
   */
  template <typename... Args> AnimT *spawn(int in_frames, Args &&...args) {
    return spawn_impl(/*pin=*/false, in_frames, std::forward<Args>(args)...);
  }

  /**
   * @brief Like spawn(), but pins the event so the returned pointer may be
   * retained across frames (e.g. registered as a live GUI param).
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no slots are free.
   * @details Only valid when the spawned animation is infinite and added before
   * any finite timeline event, so compaction never shifts it — the standard
   * retained-handle contract (see Timeline::add_get). If that invariant is ever
   * broken, step()'s compaction traps loudly instead of dangling the pointer.
   */
  template <typename... Args>
  AnimT *spawn_pinned(int in_frames, Args &&...args) {
    return spawn_impl(/*pin=*/true, in_frames, std::forward<Args>(args)...);
  }

private:
  /**
   * @brief Compact, ascending list of the slots currently active.
   * @details transform() and prepare_frame() are hot (transform runs per
   * pixel), so they iterate only the active slots — O(active) instead of
   * O(CAPACITY). Kept sorted so the transform composition order matches a
   * low-to-high slot scan; the warps are not all commutative, so the order is
   * load-bearing. Maintenance is O(active) but cold: once per spawn and once
   * per finished animation, never per pixel.
   */
  std::array<int, CAPACITY> active_slots_{};
  int active_count_ = 0; /**< Number of valid entries at the front of active_slots_. */

  /**
   * @brief Inserts a slot index into active_slots_ keeping it ascending.
   * @param idx Slot index to insert; idx is the lowest free slot, so this
   * matches the order a 0..CAPACITY scan would have applied it in.
   */
  void add_active(int idx) {
    int pos = active_count_;
    while (pos > 0 && active_slots_[pos - 1] > idx) {
      active_slots_[pos] = active_slots_[pos - 1];
      --pos;
    }
    active_slots_[pos] = idx;
    ++active_count_;
  }

  /**
   * @brief Removes a slot index from active_slots_, preserving order.
   * @param idx Slot index to drop; a no-op if it is not present.
   */
  void remove_active(int idx) {
    int pos = 0;
    while (pos < active_count_ && active_slots_[pos] != idx)
      ++pos;
    if (pos == active_count_)
      return; // already gone
    for (int k = pos; k + 1 < active_count_; ++k)
      active_slots_[k] = active_slots_[k + 1];
    --active_count_;
  }

  /**
   * @brief Allocates a free slot and schedules its animation on the timeline.
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param pin Whether the timeline event is pinned (retained-handle contract).
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no slots are free.
   */
  template <typename... Args>
  AnimT *spawn_impl(bool pin, int in_frames, Args &&...args) {
    // Linear scan for a free slot (cold path).
    for (int idx = 0; idx < CAPACITY; ++idx) {
      Entity &e = entities[idx];
      if (!e.active) {
        e.params = template_params;
        e.active = true;
        add_active(idx);
        // Create animation with reference to the stable entity params. The
        // completion callback captures `this` + the slot index (both stable:
        // the Transformer is effect-owned and holds a Timeline& so it never
        // moves) to deactivate the slot and drop it from the active list.
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        AnimT *p = timeline.add_get(in_frames, std::move(anim), pin);
        if (p) {
          // A non-pinned spawn returns a transient handle the caller does not
          // retain, so the slot is reclaimed only when the animation reaches
          // done() (the then() callback below fires once and frees it) or a
          // later cancel() routes it through the same path. A perpetual,
          // non-repeating animation (duration<0, repeat=false) reaches neither
          // and has no retained handle to cancel it, so its slot would leak;
          // after CAPACITY such spawns the pool silently returns nullptr. That
          // combination must instead use spawn_pinned and retain the handle.
          // Trap the misuse here rather than leak. (Pinned perpetual handles —
          // the documented retained-handle path — are exempt.)
          if (!pin)
            HS_CHECK(p->is_finite() || p->repeats(),
                     "Transformer::spawn of a perpetual non-repeating "
                     "animation leaks its pool slot; use spawn_pinned");
          // Attach the completion callback to the *stored* animation so it can
          // capture its pointer and re-query repeats() at callback time, not at
          // spawn time. A repeating animation that is later cancel()ed reports
          // repeats()==false and is routed through the removal branch (firing
          // this once); a spawn-time snapshot of repeats==true would skip the
          // free and leak the pool slot forever. The handle p is stable —
          // add_get traps on relocating a handled animation.
          p->then([this, idx, p]() {
            if (!p->repeats()) {
              entities[idx].active = false;
              remove_active(idx);
            }
          });
        }
        return p;
      }
    }
    // If we get here, no free slots. Drop the spawn (safe failure).
    return nullptr;
  }

public:

  /**
   * @brief Prepares per-frame cached state for all active entities.
   * @details ORDERING CONTRACT: call once per frame, before transform(), in any
   * frame where an active entity's params were *live-updated* since they were
   * last prepared (e.g. a GUI slider moved this frame). It refreshes each active
   * entity's derived state (RippleParams::prepare_thresholds /
   * NoiseParams::sync). transform() reads that state but cannot verify it is
   * current: the dependency is value-dependent (did the params change?), not
   * structural, so it is a caller contract, not a compile-time/assert invariant.
   * It is intentionally NOT required when there are no active entities (transform
   * is then the identity) or when an entity's params have not changed since spawn
   * (the spawn-time copy already prepared them) — both are exercised by the unit
   * tests, so a blanket "prepared this frame" assert would be a false positive.
   */
  void prepare_frame() {
    for (int k = 0; k < active_count_; ++k) {
      Entity &e = entities[active_slots_[k]];
      if constexpr (requires { e.params.prepare_thresholds(); }) {
        e.params.prepare_thresholds();
      }
      // Symmetric to prepare_thresholds for RippleParams: NoiseParams::sync()
      // pushes the (possibly live-updated) frequency into the embedded
      // FastNoiseLite. Centralizing it here means a noise effect cannot render a
      // stale frequency by forgetting a per-entity sync loop.
      if constexpr (requires { e.params.sync(); }) {
        e.params.sync();
      }
    }
  }

  /**
   * @brief Applies all active transformations to a vector, in slot order.
   * @param v Vector to transform.
   * @return The vector after every active transform has been composed onto it.
   * @note Reads each active entity's prepared state; see prepare_frame() for the
   * ordering contract (must run first in any frame whose params were live-updated
   * since last prepared). Per-pixel hot path — no guard here by design.
   */
  Vector transform(Vector v) const {
    for (int k = 0; k < active_count_; ++k) {
      v = TransformFunc(v, entities[active_slots_[k]].params);
    }
    return v;
  }

  /**
   * @brief Function-call alias for transform().
   * @param v Vector to transform.
   * @return The transformed vector.
   */
  Vector operator()(const Vector &v) const { return transform(v); }
};

/**
 * @brief A transformer adapter for an Orientation object.
 */
struct OrientTransformer {
  const Orientation<> &orientation; /**< Orientation applied by each transform; retained by reference. */

  /**
   * @brief Constructs an adapter wrapping an orientation.
   * @param ori Orientation to apply; retained by reference.
   */
  explicit OrientTransformer(const Orientation<> &ori) : orientation(ori) {}

  /**
   * @brief Orients a vector through the wrapped orientation.
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  Vector transform(const Vector &v) const { return orientation.orient(v); }

  /**
   * @brief Function-call alias for transform().
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  Vector operator()(const Vector &v) const { return transform(v); }
};

/**
 * @brief Applies a Mobius transformation to a vector.
 * @param v Unit vector to transform.
 * @param params Mobius transformation coefficients.
 * @return The transformed vector.
 * @details Wraps the complex math version from 3dmath.h via stereographic
 * projection: project to the plane, apply mobius, project back.
 */
inline Vector mobius_transform(const Vector &v, const MobiusParams &params) {
  return inv_stereo(mobius(stereo(v), params));
}

/**
 * @brief Applies a gnomonic Mobius transformation to a vector.
 * @param v Unit vector to transform.
 * @param params Mobius transformation coefficients.
 * @return The transformed vector.
 * @details Projects to the gnomonic plane, applies the Mobius map, then
 * projects back to the hemisphere selected by the sign of v.y.
 */
inline Vector gnomonic_mobius_transform(const Vector &v,
                                        const MobiusParams &params) {
  Complex z = gnomonic(v);
  Complex w = mobius(z, params);
  return inv_gnomonic(w, (v.y >= 0 ? 1.0f : -1.0f));
}

/**
 * @brief Rotates a point along a Ricker-wavelet ripple radiating from a center.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @return The displaced vector.
 */
inline Vector ripple_transform(const Vector &v, const RippleParams &params) {
  // Between ripples the envelope drives amplitude to 0; skip the whole
  // per-pixel wavelet (fast_acos + two expf) when there is nothing to displace.
  if (params.amplitude <= 0.001f)
    return v;

  // fast reject
  float cos_d = dot(v, params.center);
  if (cos_d > params.cos_threshold_min || cos_d < params.cos_threshold_max) {
    return v;
  }

  float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
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
/**
 * @brief Slides a point along the sphere surface by a 3D-noise field.
 * @param v The unit vector to transform.
 * @param params Noise field, scale, amplitude, speed and time.
 * @return The displaced unit vector.
 * @details Samples three decorrelated noise channels (each offset by 100/200 on
 * all three axes) to
 * build a displacement, projects it onto the tangent plane at v so the point
 * stays on the sphere, soft-caps the slide to avoid cross-hemisphere jumps,
 * then renormalizes. No-op when amplitude is negligible.
 */
inline Vector noise_transform(const Vector &v, const NoiseParams &params) {
  if (params.amplitude <= 0.001f)
    return v;

  float scale = params.scale;
  float time_val = params.time * params.speed;

  // Sample 3D noise field. Offset all three axes per channel (not just X/Y): a
  // shared Z/time input would leave the three displacement components correlated
  // along Z, so the slide direction is biased rather than isotropic. Distinct
  // per-axis offsets decorrelate every component — each channel reads a region
  // of the field far enough from the others (>> the unit-sphere coordinate span)
  // that their samples are independent.
  constexpr float kChannelYOffset = 100.0f; // channel 2 (ny) field shift
  constexpr float kChannelZOffset = 200.0f; // channel 3 (nz) field shift
  float nx =
      params.noise.GetNoise(v.x * scale, v.y * scale, v.z * scale + time_val);
  float ny = params.noise.GetNoise(v.x * scale + kChannelYOffset,
                                   v.y * scale + kChannelYOffset,
                                   v.z * scale + time_val + kChannelYOffset);
  float nz = params.noise.GetNoise(v.x * scale + kChannelZOffset,
                                   v.y * scale + kChannelZOffset,
                                   v.z * scale + time_val + kChannelZOffset);

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
  Complex coords;     /**< Displaced coordinate. */
  float displacement; /**< Magnitude of the displacement vector. */
};

/**
 * @brief Smooth pole attenuation for stereographic-space effects.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @return Falloff factor 1/(1 + r²/pole_fade²) in (0, 1].
 * @details Stereographic projection sends the far pole to infinity, so |z|²
 * grows without bound near it; this falloff is 1 at the projection origin and
 * decays toward 0 with distance, taming that singularity. Shared by
 * stereo_noise_warp and the stereo pattern effects (Flyby/Liquid2D) so the
 * falloff stays identical across warp and shading.
 */
inline float pole_attenuation(float r_sq, float pole_fade) {
  return 1.0f / (1.0f + (r_sq / (pole_fade * pole_fade)));
}

/**
 * @brief Applies noise-based displacement in stereographic space.
 * @param z Stereographic coordinate to warp.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param noise FastNoiseLite instance for sampling.
 * @param scale Noise frequency scale.
 * @param strength Maximum displacement amplitude.
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @param time Noise time coordinate (for animation).
 * @return Warped coordinate plus the scalar displacement magnitude.
 * @details Displacement is attenuated near the projection pole to prevent
 * singularity blowup. The scalar displacement is returned alongside the warped
 * coordinate for downstream effects (e.g., hue shifting).
 */
inline StereoWarpResult stereo_noise_warp(const Complex &z, float r_sq,
                                          const FastNoiseLite &noise,
                                          float scale, float strength,
                                          float pole_fade, float time) {
  float atten = pole_attenuation(r_sq, pole_fade);
  float s = strength * atten;
  float dx = noise.GetNoise(z.re * scale, z.im * scale, time) * s;
  float dy =
      noise.GetNoise(z.re * scale + 100.0f, z.im * scale + 100.0f, time) * s;
  return {Complex(z.re + dx, z.im + dy), sqrtf(dx * dx + dy * dy)};
}

/**
 * @brief Generates ripples that warp the sphere.
 * @tparam CAPACITY Maximum number of concurrent ripple transformations.
 */
template <int CAPACITY>
using RippleTransformer =
    Transformer<RippleParams, Animation::Ripple, ripple_transform, CAPACITY>;

/**
 * @brief Performs Mobius warps that return to the identity.
 * @tparam CAPACITY Maximum number of concurrent Mobius warp transformations.
 */
template <int CAPACITY>
using MobiusWarpTransformer =
    Transformer<MobiusParams, Animation::MobiusWarp, mobius_transform, CAPACITY>;

/**
 * @brief Performs circular Mobius warps that stay warped throughout, suitable
 * for repeating animations.
 * @tparam CAPACITY Maximum number of concurrent circular Mobius warps.
 * @warning This variant deliberately never returns to identity — unlike
 * MobiusWarpTransformer (which uses the same `mobius_transform` but an animation
 * that eases back to identity), `Animation::MobiusWarpCircular` traces a closed
 * loop in parameter space that holds the warp at full strength. That makes it
 * correct ONLY in a repeating slot, where the loop re-enters seamlessly. In a
 * non-repeating slot it freezes off-identity on the final composed frame, so the
 * teardown is a one-frame discontinuity (the warp snaps away when the slot ends)
 * rather than a smooth relax-to-identity. Use MobiusWarpTransformer for one-shot
 * slots that must land back on the unwarped sphere.
 */
template <int CAPACITY>
using MobiusWarpCircularTransformer =
    Transformer<MobiusParams, Animation::MobiusWarpCircular, mobius_transform,
                CAPACITY>;

/**
 * @brief Performs a changing Mobius warp using gnomonic projection.
 * @tparam CAPACITY Maximum number of concurrent gnomonic Mobius warps.
 */
template <int CAPACITY>
using MobiusWarpGnomonicTransformer =
    Transformer<MobiusParams, Animation::MobiusWarpEvolving,
                gnomonic_mobius_transform, CAPACITY>;

/**
 * @brief Applies 3D noise distortion to vectors.
 * @tparam CAPACITY Maximum number of concurrent noise transformations.
 */
template <int CAPACITY>
using NoiseTransformer =
    Transformer<NoiseParams, Animation::Noise, noise_transform, CAPACITY>;
