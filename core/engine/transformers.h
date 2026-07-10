/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "math/3dmath.h"
#include "engine/concepts.h"
#include <array>
#include "vendor/FastNoiseLite.h"
#include "animation/animation.h"

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
  Timeline &timeline;                   /**< Timeline that schedules and steps the spawned animations. */

  /**
   * @brief Constructs a transformer bound to a timeline.
   * @param tl Timeline used to schedule spawned animations; retained by reference.
   */
  Transformer(Timeline &tl) : timeline(tl) {}

  // spawn_impl's one-shot callbacks capture this+slot index; relocation would
  // dangle them, so the object is fixed in place.
  Transformer(const Transformer &) = delete;
  Transformer(Transformer &&) = delete;

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
   * @details Only valid when the spawned animation never completes on its own —
   * infinite, or repeating (it rewinds rather than reaching done()) — and is
   * added before any finite timeline event, so compaction never shifts it: the
   * standard retained-handle contract (see Timeline::add_get). If that invariant
   * is ever broken, step()'s compaction traps loudly instead of dangling it.
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

  std::array<Entity, CAPACITY> entities; /**< Fixed-capacity pool of transformation slots. */

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
          // A non-pinned spawn keeps no retained handle, so the slot is
          // reclaimed only by the one-shot then() below — which fires only when
          // the animation reaches done() once and is removed. That happens only
          // for a finite, non-repeating animation: an infinite one never reaches
          // done(), and a repeating one rewinds instead of being removed, so
          // either would hold its slot for the effect's life with no handle to
          // cancel it (and after CAPACITY such spawns the pool returns nullptr).
          // Those must use spawn_pinned and the retained-handle contract.
          if (!pin)
            HS_CHECK(p->is_finite() && !p->repeats(),
                     "Transformer::spawn needs a finite, non-repeating "
                     "animation; infinite or repeating spawns leak their pool "
                     "slot — use spawn_pinned");
          // Recycle the pool slot at the animation's final removal. The callback
          // must tell a removal apart from a mid-repeat post fire, which it reads
          // from repeats(); how it reads repeats() differs by pin because the two
          // paths have opposite handle stability:
          //   - Pinned: the event never relocates (add_get traps on moving a
          //     handled event) and may be cancel()ed through the retained
          //     pointer, flipping repeats() to false. Re-query live through p so
          //     a post-cancel teardown routes through the free branch.
          //   - Non-pinned: the event is compacted (relocated) by step() when an
          //     earlier event completes, so p must not be retained. The transient
          //     handle is never cancel()ed, so repeats() is fixed at spawn —
          //     snapshot it rather than dereference a possibly-relocated slot.
          if (pin) {
            p->then([this, idx, p]() {
              if (!p->repeats()) {
                entities[idx].active = false;
                remove_active(idx);
              }
            });
          } else {
            const bool repeats = p->repeats();
            p->then([this, idx, repeats]() {
              if (!repeats) {
                entities[idx].active = false;
                remove_active(idx);
              }
            });
          }
        } else {
          // Timeline pool full: undo the activation so the slot is not leaked
          // (its reclaim callback above never registered).
          e.active = false;
          remove_active(idx);
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
   * last prepared (e.g. a GUI slider moved this frame). It re-reads live config
   * from template_params (NoiseParams::refresh_from, when provided) and refreshes
   * each active entity's derived state (NoiseParams::sync). transform() reads
   * that state but cannot verify it is
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
      // Pull live-tunable config from template_params into the spawned entity:
      // slider edits land on template_params, not the spawn-time copy. Runs
      // before sync() so the refreshed frequency reaches the generator.
      if constexpr (requires { e.params.refresh_from(template_params); }) {
        e.params.refresh_from(template_params);
      }
      // NoiseParams::sync() pushes the (possibly live-updated) frequency into
      // the embedded FastNoiseLite. Centralizing it here means a noise effect
      // cannot render a stale frequency by forgetting a per-entity sync loop.
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
  // per-pixel wavelet (fast_acos + two fast_expf) when there is nothing to
  // displace.
  if (params.amplitude <= 0.001f)
    return v;

  // fast reject
  float cos_d = dot(v, params.center);
  if (cos_d > params.cos_threshold_min || cos_d < params.cos_threshold_max) {
    return v;
  }

  float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  float dist_from_peak = d - params.phase;

  float half_width = params.half_width();

  // Normalize distance; within the accept band |t| reaches 4 (−2..2 is the primary lobe)
  float t = (dist_from_peak / half_width) * 2.0f;

  // Ricker Wavelet: (1 - t^2) * e^(-t^2/2)
  float ricker = (1.0f - t * t) * fast_expf(-0.5f * t * t);

  // Distance attenuation (ripples get smaller as they spread)
  float attenuation = fast_expf(-params.decay * d);

  float theta = params.amplitude * ricker * attenuation;

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

  // ny/nz read the same field as nx under a constant spatial translation
  // (100/200 on every axis). The channels decorrelate because that translation
  // exceeds the noise correlation length, not because the per-axis offsets
  // differ. The same time_val drives the z input of all three, so they animate
  // together in time.
  constexpr float CHANNEL_Y_OFFSET = 100.0f; // channel 2 (ny) field shift
  constexpr float CHANNEL_Z_OFFSET = 200.0f; // channel 3 (nz) field shift
  float nx =
      params.noise.GetNoise(v.x * scale, v.y * scale, v.z * scale + time_val);
  float ny = params.noise.GetNoise(v.x * scale + CHANNEL_Y_OFFSET,
                                   v.y * scale + CHANNEL_Y_OFFSET,
                                   v.z * scale + time_val + CHANNEL_Y_OFFSET);
  float nz = params.noise.GetNoise(v.x * scale + CHANNEL_Z_OFFSET,
                                   v.y * scale + CHANNEL_Z_OFFSET,
                                   v.z * scale + time_val + CHANNEL_Z_OFFSET);

  Vector raw_noise = Vector(nx, ny, nz) * (params.amplitude * 0.05f);

  // Project noise onto the tangent plane at v.
  float inward_pull = dot(raw_noise, v);
  Vector surface_distortion = raw_noise - (v * inward_pull);

  // Soft-cap the slide distance to prevent cross-hemisphere grabs.
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
  // Floor the radius so a 0 pole_fade can't divide by zero and poison the warp.
  const float pf = pole_fade > 1e-3f ? pole_fade : 1e-3f;
  return 1.0f / (1.0f + (r_sq / (pf * pf)));
}

/**
 * @brief Pole-attenuates a stereographic pattern value and maps it to [0, 1].
 * @param pattern Raw pattern value in [-1, 1].
 * @param r_sq Pre-computed |z|² driving the pole fade.
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @return Pole-attenuated value normalized to [0, 1].
 * @details Shared by the stereo pattern effects (Flyby/Liquid2D) so the
 * attenuation and normalization stay identical across them.
 */
inline float pole_normalize_pattern(float pattern, float r_sq, float pole_fade) {
  return (pattern * pole_attenuation(r_sq, pole_fade) + 1.0f) * 0.5f;
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
  constexpr float CHANNEL_OFFSET = 100.0f; // decorrelates dy's field from dx's
  float dx = noise.GetNoise(z.re * scale, z.im * scale, time) * s;
  float dy = noise.GetNoise(z.re * scale + CHANNEL_OFFSET,
                            z.im * scale + CHANNEL_OFFSET, time) *
             s;
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
