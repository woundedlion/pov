/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "math/3dmath.h"
#include "engine/concepts.h"
#include <new>
#include "vendor/FastNoiseLite.h"
#include "animation/animation.h"

using Animation::BumpParams;
using Animation::NoiseParams;
using Animation::NoiseProductParams;
using Animation::RippleParams;

/**
 * @brief Fixed-capacity pool of animation-driven parameter entities.
 * @tparam ParamsT The configuration struct (e.g., RippleParams, MobiusParams).
 * @tparam AnimT The animation class (e.g., Animation::Ripple).
 * @tparam CAPACITY Max number of active entities.
 * @details Owns the entity slots, their timeline lifecycle (spawn, completion
 * reclaim), and the per-frame param refresh. Derived classes add the hot-path
 * composition over the active entities (Transformer composes Vector warps,
 * FieldTransformer sums scalar fields).
 */
template <typename ParamsT, typename AnimT, int CAPACITY = 32>
class TransformerPool {
public:
  /**
   * @brief Per-slot storage for one active entity.
   */
  struct Entity {
    ParamsT params;         /**< Per-entity configuration the composition reads. */
    bool active = false;    /**< Whether this slot currently holds a live animation. */
  };

  ParamsT template_params;              /**< Template params copied into each new entity on spawn. */
  Timeline &timeline;                   /**< Timeline that schedules and steps the spawned animations. */

  /**
   * @brief Constructs a pool bound to a timeline.
   * @param tl Timeline used to schedule spawned animations; retained by reference.
   */
  TransformerPool(Timeline &tl) : timeline(tl) {}

  // spawn_impl's one-shot callbacks capture this+slot index; relocation would
  // dangle them, so the object is fixed in place.
  TransformerPool(const TransformerPool &) = delete;
  TransformerPool(TransformerPool &&) = delete;

  /**
   * @brief Allocates the entity pool from the persistent arena.
   * @param arena Persistent arena supplying CAPACITY entity slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet), after any configure_arenas() and before the first spawn.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    entities = arena.allocate_n<Entity>(CAPACITY);
    for (int i = 0; i < CAPACITY; ++i)
      new (&entities[i]) Entity();
    active_slots_ = arena.allocate_n<int>(CAPACITY);
    active_count_ = 0;
  }

  /**
   * @brief Re-claims the pool's storage after its arena was reset, preserving
   * live entities.
   * @param arena The arena init_storage() allocated from, freshly reset.
   * @details For arenas that are compacted mid-effect (e.g. a mesh carousel's
   * after-reset callback). Spawned animations hold Params references into the
   * slots, so the blocks must re-land at their original addresses — the caller
   * must replay the same allocation order after the reset as after
   * init_storage() (asserted). The bytes are left untouched: an arena reset
   * only rewinds the offset, so live entities carry through.
   */
  HS_COLD_MEMBER void reclaim_storage(Arena &arena) {
    Entity *e = arena.allocate_n<Entity>(CAPACITY);
    int *s = arena.allocate_n<int>(CAPACITY);
    HS_CHECK(e == entities && s == active_slots_,
             "TransformerPool: reclaimed storage moved");
  }

  /**
   * @brief Number of currently active entities.
   * @return Count of live pool slots.
   */
  int active_count() const { return active_count_; }

  /**
   * @brief Params of the k-th active entity, in spawn order.
   * @param k Active index in [0, active_count()).
   * @return The entity's live params.
   */
  const ParamsT &active_params(int k) const {
    return entities[active_slots_[k]].params;
  }

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

protected:
  /**
   * @brief Compact list of the active slots, in spawn order.
   * @details The derived compositions and prepare_frame() are hot (they run per
   * pixel), so they iterate only the active slots — O(active) instead of O(CAPACITY).
   * Held in spawn order (append on activation, order-preserving removal) so the
   * composition order follows spawn order; the warps are not all commutative, so the
   * order is load-bearing and must not depend on which freed slot was recycled.
   */
  int *active_slots_ = nullptr;
  int active_count_ = 0; /**< Number of valid entries at the front of active_slots_. */

  Entity *entities = nullptr; /**< CAPACITY-slot pool, allocated by init_storage(). */

private:
  /**
   * @brief Appends a slot index to active_slots_ in spawn order.
   * @param idx Slot index to insert; appended at the end so composition order
   * follows spawn order regardless of which freed slot was recycled.
   */
  void add_active(int idx) { active_slots_[active_count_++] = idx; }

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
    HS_CHECK(entities, "TransformerPool: call init_storage() before spawn");
    // Linear scan for a free slot (cold path).
    for (int idx = 0; idx < CAPACITY; ++idx) {
      Entity &e = entities[idx];
      if (!e.active) {
        e.params = template_params;
        e.active = true;
        add_active(idx);
        // Completion callback captures `this` + the slot index (both stable) to
        // deactivate the slot and drop it from the active list.
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        AnimT *p = timeline.add_get(in_frames, std::move(anim), pin);
        if (p) {
          // A non-pinned spawn keeps no retained handle, so the slot is reclaimed
          // only by the one-shot then() below, which fires when the animation
          // reaches done() once. A finite, non-repeating animation is required: an
          // infinite one never reaches done() and a repeating one rewinds instead
          // of being removed, so either would hold its slot for the effect's life
          // (nullptr after CAPACITY spawns) — use spawn_pinned.
          if (!pin)
            HS_CHECK(p->is_finite() && !p->repeats(),
                     "Transformer::spawn needs a finite, non-repeating "
                     "animation; infinite or repeating spawns leak their pool "
                     "slot — use spawn_pinned");
          // Recycle the pool slot at final removal. The paths differ by handle
          // stability:
          //   - Pinned: the event never relocates and may be cancel()ed through
          //     the retained pointer (flipping repeats() to false), so re-query
          //     live through p to tell a removal from a mid-repeat post fire.
          //   - Non-pinned: p must not be retained (step() compacts the event),
          //     and the HS_CHECK above pins repeats() false, so fire once at done().
          if (pin) {
            p->then([this, idx, p]() {
              if (!p->repeats()) {
                entities[idx].active = false;
                remove_active(idx);
              }
            });
          } else {
            p->then([this, idx]() {
              entities[idx].active = false;
              remove_active(idx);
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
   * @details ORDERING CONTRACT: call once per frame, before the derived
   * composition (transform()/field()), in any frame where an active entity's
   * params were *live-updated* (e.g. a GUI slider moved). It re-reads live config
   * from template_params (NoiseParams::refresh_from) and refreshes each active
   * entity's derived state (NoiseParams::sync). The composition reads that state
   * but cannot verify it is current: the dependency is
   * value-dependent, so it is a caller contract, not an assert. NOT required when
   * there are no active entities (transform is then the identity) or when params
   * have not changed since spawn.
   */
  void prepare_frame() {
    for (int k = 0; k < active_count_; ++k) {
      Entity &e = entities[active_slots_[k]];
      // Pull live-tunable config from template_params into the spawned entity
      // (slider edits land on template_params, not the spawn-time copy); runs
      // before sync() so the refreshed frequency reaches the generator.
      if constexpr (requires { e.params.refresh_from(template_params); }) {
        e.params.refresh_from(template_params);
      }
      // sync() pushes the (possibly live-updated) frequency into the embedded
      // FastNoiseLite. Centralizing it here means a noise effect cannot render a
      // stale frequency by forgetting a per-entity sync loop.
      if constexpr (requires { e.params.sync(); }) {
        e.params.sync();
      }
    }
  }

};

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
class Transformer : public TransformerPool<ParamsT, AnimT, CAPACITY> {
public:
  using TransformerPool<ParamsT, AnimT, CAPACITY>::TransformerPool;

  /**
   * @brief Applies all active transformations to a vector, in slot order.
   * @param v Vector to transform.
   * @return The vector after every active transform has been composed onto it.
   * @note Reads each active entity's prepared state; see prepare_frame() for the
   * ordering contract (must run first in any frame whose params were live-updated
   * since last prepared). Per-pixel hot path — no guard here by design.
   */
  Vector transform(Vector v) const {
    for (int k = 0; k < this->active_count_; ++k) {
      v = TransformFunc(v, this->entities[this->active_slots_[k]].params);
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
 * @brief A generic manager for animation-driven scalar displacement fields.
 * @tparam ParamsT The configuration struct (e.g., BumpParams).
 * @tparam AnimT The animation class (e.g., Animation::BallDrop).
 * @tparam FieldFunc The static function evaluating one entity's field.
 * @tparam CAPACITY Max number of active fields.
 * @details The scalar counterpart of Transformer: entities superpose by
 * summation instead of composing as warps, so a caller can feed the summed
 * field into a displacement path (e.g. a DistortedRing shift LUT).
 */
template <typename ParamsT, typename AnimT,
          float (*FieldFunc)(const Vector &, const ParamsT &),
          int CAPACITY = 32>
class FieldTransformer : public TransformerPool<ParamsT, AnimT, CAPACITY> {
public:
  using TransformerPool<ParamsT, AnimT, CAPACITY>::TransformerPool;

  /**
   * @brief Sums every active entity's field at a point.
   * @param p Sample point (unit vector).
   * @return The superposed field value; 0 with no active entities.
   * @note Reads each active entity's prepared state; see prepare_frame() for the
   * ordering contract. Per-sample hot path — no guard here by design.
   */
  float field(const Vector &p) const {
    float s = 0.0f;
    for (int k = 0; k < this->active_count_; ++k) {
      s += FieldFunc(p, this->entities[this->active_slots_[k]].params);
    }
    return s;
  }

  /**
   * @brief Function-call alias for field().
   * @param p Sample point (unit vector).
   * @return The superposed field value.
   */
  float operator()(const Vector &p) const { return field(p); }

  /**
   * @brief Magnitude-weighted blend of the active entities' fields: the
   * strongest contribution dominates without stacking.
   * @param p Sample point (unit vector).
   * @return sum(s_i^3) / sum(s_i^2); 0 with no active entities.
   * @details Use instead of field() when overlapping entities must not add
   * (e.g. solid bodies displacing a shared sheet). Unlike a hard max by
   * magnitude — which jumps discontinuously where opposite-signed fields cross
   * in strength — the blend is continuous everywhere: a single entity passes
   * through exactly, equal same-signed overlaps yield the shared value, and
   * opposite-signed overlaps cancel smoothly.
   */
  float field_dominant(const Vector &p) const {
    float num = 0.0f;
    float den = 0.0f;
    for (int k = 0; k < this->active_count_; ++k) {
      float f = FieldFunc(p, this->entities[this->active_slots_[k]].params);
      num += f * f * f;
      den += f * f;
    }
    return den > 1e-9f ? num / den : 0.0f;
  }

  /**
   * @brief field_dominant() restricted to a subset of the active entities.
   * @param p Sample point (unit vector).
   * @param ks Active indices (as accepted by active_params()).
   * @param n Number of indices.
   * @return sum(s_i^3) / sum(s_i^2) over the subset; 0 when empty.
   * @details Exact when every excluded entity's field is 0 at @p p (zero terms
   * do not move the blend), so callers can prefilter with per-entity bounds.
   */
  float field_dominant(const Vector &p, const int *ks, int n) const {
    float num = 0.0f;
    float den = 0.0f;
    for (int j = 0; j < n; ++j) {
      float f = FieldFunc(p, this->active_params(ks[j]));
      num += f * f * f;
      den += f * f;
    }
    return den > 1e-9f ? num / den : 0.0f;
  }

  /**
   * @brief Upper bound on |field()| over the sphere this frame.
   * @return Sum of the active entities' per-entity bounds.
   * @details Requires ParamsT::field_bound() (a true upper bound on
   * |FieldFunc|); callers use it to size conservative culls.
   */
  float field_bound() const {
    float b = 0.0f;
    for (int k = 0; k < this->active_count_; ++k) {
      b += this->entities[this->active_slots_[k]].params.field_bound();
    }
    return b;
  }
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
  // Between ripples the envelope drives amplitude to 0; skip the whole per-pixel
  // wavelet (fast_acos + two fast_expf) when there is nothing to displace.
  if (params.amplitude <= 0.001f)
    return v;

  // Fast reject outside the [d_min, d_max] angular band. cos decreases with
  // angle, so cos_threshold_min holds the LARGER cosine (nearest angle d_min)
  // and cos_threshold_max the smaller (farthest d_max) — the ordering reads
  // inverted but is correct; do not "fix" it.
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
 * @brief Evaluates a spherical-cap drape push: rings bow away from the cap
 * center as if draping over a ball beneath them.
 * @param v Sample point (unit vector).
 * @param params Bump center, stack axis, footprint and lifecycle envelope
 * (the envelope scales the effective footprint, inflating/deflating the cap).
 * @return The signed polar displacement (radians): the depth inside the cap's
 * boundary arc, weighted by a drape factor that is zero for the ring through
 * the center (it rides straight over the top) and zero at the footprint edge
 * (the ball's equator), peaking between. The amplitude gain scales the drape
 * weight, saturating at 1 (full clearance to the boundary arc), so the gain
 * morphs the look from a soft drape toward a solid punch-through. Positive
 * pushes toward larger colatitude about the axis; points outside the cap are
 * untouched.
 */
inline float bump_field(const Vector &v, const BumpParams &params) {
  float r_eff = params.radius * params.envelope;
  if (r_eff <= 1e-3f || params.amplitude <= 0.001f)
    return 0.0f;
  float cos_d = dot(v, params.center);
  if (cos_d <= params.cos_radius)
    return 0.0f;
  float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  if (d >= r_eff)
    return 0.0f;

  // Local cap coords: signed polar offset y from the center (positive toward
  // larger colatitude) and azimuthal offset x, with x^2 = d^2 - y^2 by the
  // small-cap approximation. The boundary arc at this azimuth sits at
  // +-sqrt(r_eff^2 - x_sq), so (arc - |y|) is the polar depth inside the cap —
  // an arc-shaped profile along the ring, which keeps the bulge round.
  float y = fast_acos(hs::clamp(dot(params.axis, v), -1.0f, 1.0f)) -
            fast_acos(hs::clamp(dot(params.axis, params.center), -1.0f, 1.0f));
  float abs_y = std::fabs(y);
  float x_sq = std::max(d * d - y * y, 0.0f);
  float depth = sqrtf(std::max(r_eff * r_eff - x_sq, 0.0f)) - abs_y;
  float drape =
      std::min(params.amplitude * sinf(PI_F * abs_y / r_eff), 1.0f);
  return copysignf(depth * drape, y);
}

/**
 * @brief Evaluates a two-octave product noise field at a point.
 * @param v Sample point (unit vector).
 * @param params Octave scales, amplitude, and field time.
 * @return The field value at @p v.
 * @details Octave 1 envelopes octave 2, so perturbations bunch where the
 * envelope is strong and vanish where it crosses zero.
 */
inline float noise_product_field(const Vector &v,
                                 const NoiseProductParams &params) {
  if (std::fabs(params.amplitude) <= 0.001f)
    return 0.0f;
  float n1 = params.noise.GetNoise(v.x * params.scale1, v.y * params.scale1,
                                   v.z * params.scale1 + params.time);
  float n2 = params.noise.GetNoise(
      v.x * params.scale2 + NoiseProductParams::OCTAVE2_OFFSET,
      v.y * params.scale2, v.z * params.scale2 + params.time);
  return params.amplitude * n1 * n2;
}

/**
 * @brief Generates ripples that warp the sphere.
 * @tparam CAPACITY Maximum number of concurrent ripple transformations.
 */
template <int CAPACITY>
using RippleTransformer =
    Transformer<RippleParams, Animation::Ripple, ripple_transform, CAPACITY>;

/**
 * @brief Bump displacement fields that fall pole-to-pole through a frame.
 * @tparam CAPACITY Maximum number of concurrent falling bumps.
 */
template <int CAPACITY>
using BallDropTransformer =
    FieldTransformer<BumpParams, Animation::BallDrop, bump_field, CAPACITY>;

/**
 * @brief A two-octave product noise displacement field.
 * @tparam CAPACITY Maximum number of concurrent noise fields.
 */
template <int CAPACITY>
using NoiseProductTransformer =
    FieldTransformer<NoiseProductParams, Animation::NoiseProduct,
                     noise_product_field, CAPACITY>;

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
 * @warning This variant never returns to identity — unlike MobiusWarpTransformer
 * (same `mobius_transform` but an animation that eases back to identity),
 * `Animation::MobiusWarpCircular` traces a closed loop that holds the warp at full
 * strength. Correct ONLY in a repeating slot, where the loop re-enters seamlessly;
 * in a non-repeating slot it freezes off-identity on the final composed frame (a
 * one-frame teardown discontinuity). Use MobiusWarpTransformer for one-shot slots
 * that must land back on the unwarped sphere.
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
