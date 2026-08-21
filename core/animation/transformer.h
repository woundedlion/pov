/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file transformer.h
 * @brief TransformerPool and its Transformer, FieldTransformer and
 *        OrientTransformer specializations: animation-driven warp and field
 *        composition.
 * @details Also holds the free warp and field functions those pools compose
 *          over Animation params. The Mobius sphere maps the Mobius pools
 *          compose live in math/stereographic.h alongside their coefficients.
 */

#include "math/3dmath.h"
#include "math/stereographic.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include <new>
#include <type_traits>
#include "vendor/FastNoiseLite.h"
#include "animation/animation.h"

/**
 * @brief A params type's declared refresh_from-hook intent, false when
 *        undeclared.
 * @tparam T Candidate params type.
 */
template <typename T> constexpr bool declared_needs_refresh_from() {
  if constexpr (requires { T::NEEDS_REFRESH_FROM; })
    return T::NEEDS_REFRESH_FROM;
  else
    return false;
}

/**
 * @brief A params type's declared sync-hook intent, false when undeclared.
 * @tparam T Candidate params type.
 */
template <typename T> constexpr bool declared_needs_sync() {
  if constexpr (requires { T::NEEDS_SYNC; })
    return T::NEEDS_SYNC;
  else
    return false;
}

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
    ParamsT params; /**< Per-entity configuration the composition reads. */
    bool active =
        false; /**< Whether this slot currently holds a live animation. */
  };

  static_assert(CAPACITY > 0, "TransformerPool requires CAPACITY >= 1");

  static_assert(std::is_trivially_destructible_v<ParamsT>,
                "TransformerPool placement-news CAPACITY entities into the "
                "arena and never destroys them, so ParamsT must own no state "
                "outside its slot.");

  /** @brief Whether ParamsT exposes prepare_frame()'s live-config hook. */
  static constexpr bool HAS_REFRESH_FROM =
      requires(ParamsT &p, ParamsT &t) { p.refresh_from(t); };
  /** @brief Whether ParamsT exposes prepare_frame()'s derived-state hook. */
  static constexpr bool HAS_SYNC = requires(ParamsT &p) { p.sync(); };

  static_assert(
      requires { ParamsT::NEEDS_REFRESH_FROM; },
      "ParamsT must declare `static constexpr bool NEEDS_REFRESH_FROM`: true "
      "if it carries prepare_frame()'s refresh_from() hook, false if it does "
      "not.");
  static_assert(
      requires { ParamsT::NEEDS_SYNC; },
      "ParamsT must declare `static constexpr bool NEEDS_SYNC`: true if it "
      "carries prepare_frame()'s sync() hook, false if it does not.");
  static_assert(
      declared_needs_refresh_from<ParamsT>() == HAS_REFRESH_FROM,
      "ParamsT::NEEDS_REFRESH_FROM disagrees with whether ParamsT exposes "
      "refresh_from(). prepare_frame() finds each hook by detection "
      "independently, so a renamed hook or a drifted signature would leave the "
      "entity unrefreshed with no other signal.");
  static_assert(
      declared_needs_sync<ParamsT>() == HAS_SYNC,
      "ParamsT::NEEDS_SYNC disagrees with whether ParamsT exposes sync(). "
      "prepare_frame() finds each hook by detection independently, so a "
      "renamed hook or a drifted signature would leave the entity's derived "
      "state stale with no other signal.");

  ParamsT
      template_params; /**< Template params copied into each new entity on spawn. */
  Timeline &
      timeline; /**< Timeline that schedules and steps the spawned animations. */

  /**
   * @brief Constructs a pool bound to a timeline.
   * @param tl Timeline used to schedule spawned animations; retained by reference.
   */
  HS_COLD_MEMBER TransformerPool(Timeline &tl)
      : timeline(tl), pool_serial(next_pool_serial++) {
    next_live = live_head;
    live_head = this;
  }

  /**
   * @brief Drops the pool's clear hook and its liveness record.
   * @details A spawned animation's completion callback outlives the pool
   * whenever the timeline does; dropping the liveness record is what turns
   * those callbacks into no-ops instead of writes through a dead pool.
   */
  HS_COLD_MEMBER ~TransformerPool() {
    unlink_live();
    if (clear_hook_registered) {
      // The timeline reference is used below, so an effect that declares its
      // Timeline after its pools is a use-after-free, not a contract to bend.
      HS_CHECK(global_timeline_live,
               "TransformerPool outlived its Timeline: declare the Timeline "
               "before the pools that schedule on it");
      timeline.remove_clear_hook(this);
    }
  }

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
    HS_CHECK(!entities, "TransformerPool: init_storage() called twice");
    entities = arena.make_n<Entity>(CAPACITY);
    active_slots = arena.allocate_n<int>(CAPACITY);
    active_slot_count = 0;
#ifndef NDEBUG
    stamp.record(arena);
#endif
    if (!clear_hook_registered) {
      timeline.add_clear_hook(this, [](void *self) {
        static_cast<TransformerPool *>(self)->release_all();
      });
      clear_hook_registered = true;
    }
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
    HS_CHECK(entities,
             "TransformerPool: call init_storage() before reclaim_storage");
    Entity *e = arena.allocate_n<Entity>(CAPACITY);
    int *s = arena.allocate_n<int>(CAPACITY);
    HS_CHECK(e == entities && s == active_slots,
             "TransformerPool: reclaimed storage moved");
#ifndef NDEBUG
    assert(stamp.source_arena == &arena &&
           "TransformerPool::reclaim_storage() on a different arena than "
           "init_storage() allocated from");
    stamp.record(arena);
#endif
  }

  /**
   * @brief Number of currently active entities.
   * @return Count of live pool slots.
   */
  int active_count() const {
    check_storage_alive();
    return active_slot_count;
  }

  /**
   * @brief Params of the k-th active entity, in spawn order.
   * @param k Active index in [0, active_count()).
   * @return The entity's live params.
   */
  const ParamsT &active_params(int k) const {
    check_storage_alive();
    HS_CHECK(k >= 0 && k < active_slot_count,
             "TransformerPool: active index out of range");
    return entities[active_slots[k]].params;
  }

  /**
   * @brief Spawns a new transformation animation.
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot or
   * timeline event is available.
   * @details pin=false: the returned pointer is transient (used at the call
   * site, not retained across frames). These animations are often finite and
   * are compacted normally; pinning them would trap on routine completion.
   * The pool claims the animation's single then() slot to recycle the entity,
   * so the caller must not attach one (Animation::then() traps on a second).
   */
  template <typename... Args> AnimT *spawn(int in_frames, Args &&...args) {
    return spawn_impl(Timeline::Pin::UNPINNED, in_frames,
                      std::forward<Args>(args)...);
  }

  /**
   * @brief Like spawn(), but pins the event so the returned pointer may be
   * retained across frames (e.g. registered as a live GUI param).
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot is
   * available. A full timeline traps here instead of returning null, which a
   * retained handle would not check.
   * @details Only valid when the spawned animation never completes on its own —
   * infinite, or repeating (it rewinds rather than reaching done()) — and is
   * added before any finite timeline event, so compaction never shifts it: the
   * standard retained-handle contract (see Timeline::add_get). If that invariant
   * is ever broken, step()'s compaction traps loudly instead of dangling it.
   * The pool claims the animation's single then() slot to recycle the entity,
   * so the retained handle must not attach one (Animation::then() traps on a
   * second).
   */
  template <typename... Args>
  AnimT *spawn_pinned(int in_frames, Args &&...args) {
    return spawn_impl(Timeline::Pin::PINNED, in_frames,
                      std::forward<Args>(args)...);
  }

  /**
   * @brief Prepares per-frame cached state for all active entities.
   * @details ORDERING CONTRACT: call before the derived composition
   * (transform()/field()) whenever active params changed since their previous
   * preparation, whether through animation or live config. It re-reads live
   * config from template_params and refreshes each active entity's derived
   * state. The composition reads that state but cannot verify it is current.
   * NOT required when there are no active entities or when params are unchanged.
   */
  void prepare_frame() {
    check_storage_alive();
    for (int k = 0; k < active_slot_count; ++k) {
      Entity &e = entities[active_slots[k]];
      // Pull live-tunable config from template_params into the spawned entity.
      if constexpr (HAS_REFRESH_FROM) {
        e.params.refresh_from(template_params);
      }
      // Refresh derived state after copying live values.
      if constexpr (HAS_SYNC) {
        e.params.sync();
      }
    }
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
  int *active_slots = nullptr;
  int active_slot_count =
      0; /**< Number of valid entries at the front of active_slots. */

  Entity *entities =
      nullptr; /**< CAPACITY-slot pool, allocated by init_storage(). */

private:
  bool clear_hook_registered =
      false; /**< Whether init_storage() registered the timeline clear hook. */

  static inline TransformerPool *live_head =
      nullptr; /**< Head of the live-pool list. */
  static inline uint32_t next_pool_serial = 0; /**< Source of pool_serial. */
  TransformerPool *next_live = nullptr; /**< Next link in the live list. */
  uint32_t pool_serial; /**< Identity a completion callback tests against. */

  /**
   * @brief Whether @p pool is a live instance still carrying @p serial.
   * @param pool Pool pointer captured by a completion callback; only compared,
   *        never dereferenced unless it is found live.
   * @param serial The pool_serial captured alongside it.
   * @return True iff that exact pool is still constructed.
   * @details Timeline events outlive a pool held in a narrower scope than the
   * timeline, and their completion callbacks reach back into it. The serial
   * separates a destroyed pool from a later one built at the same address.
   */
  HS_COLD_MEMBER static bool is_live(const TransformerPool *pool,
                                     uint32_t serial) {
    for (const TransformerPool *p = live_head; p; p = p->next_live) {
      if (p == pool)
        return p->pool_serial == serial;
    }
    return false;
  }

  /** @brief Unlinks this pool from the live list. */
  HS_COLD_MEMBER void unlink_live() {
    for (TransformerPool **p = &live_head; *p; p = &(*p)->next_live) {
      if (*p == this) {
        *p = next_live;
        return;
      }
    }
    HS_CHECK(false, "TransformerPool: destroyed pool not in the live list");
  }

#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state at the last
                              init_storage()/reclaim_storage(). */

  /**
   * @brief Debug-only use-after-free check on the pool's arena blocks.
   * @details Asserts if the arena was reset or rebound after init_storage() —
   * the slots then alias reissued bytes — or rewound below either block, either
   * while the bytes are still uncovered or after a later allocation reissued
   * them. The raw slot pointers stay non-null through all of it, so nothing else
   * detects it.
   */
  void check_storage_alive() const {
    constexpr size_t ENTITY_BYTES = CAPACITY * sizeof(Entity);
    constexpr size_t SLOT_BYTES = CAPACITY * sizeof(int);
    assert(!stamp.arena_reset() && "TransformerPool use-after-free!");
    assert(!stamp.block_uncovered(entities, ENTITY_BYTES) &&
           !stamp.block_uncovered(active_slots, SLOT_BYTES) &&
           "TransformerPool use-after-free (arena rewound below block)!");
    assert(!stamp.block_reissued(entities, ENTITY_BYTES) &&
           !stamp.block_reissued(active_slots, SLOT_BYTES) &&
           "TransformerPool use-after-free (slots reclaimed by a rewind and "
           "reissued)!");
  }
#else
  /**
   * @brief No-op use-after-free check in release builds.
   */
  void check_storage_alive() const {}
#endif

  /**
   * @brief Frees every slot without touching the timeline.
   * @details Reachable only through the clear hook init_storage() registers:
   * slots are normally reclaimed by each animation's completion callback, and
   * Timeline::clear() destroys events without running those, so spawning past
   * CAPACITY would otherwise return nullptr forever. Freeing slots while their
   * animations are still live would hand a recycled slot's ParamsT to a second
   * entity and let the stale completion callback deactivate it.
   */
  HS_COLD_MEMBER void release_all() {
    HS_CHECK(entities,
             "TransformerPool: call init_storage() before release_all");
    check_storage_alive();
    for (int i = 0; i < CAPACITY; ++i)
      entities[i].active = false;
    active_slot_count = 0;
  }

  /**
   * @brief Appends a slot index to active_slots in spawn order.
   * @param idx Slot index to insert; appended at the end so composition order
   * follows spawn order regardless of which freed slot was recycled.
   */
  void add_active(int idx) { active_slots[active_slot_count++] = idx; }

  /**
   * @brief Removes a slot index from active_slots, preserving order.
   * @param idx Slot index to drop; a no-op if it is not present.
   */
  void remove_active(int idx) {
    int pos = 0;
    while (pos < active_slot_count && active_slots[pos] != idx)
      ++pos;
    if (pos == active_slot_count)
      return; // already gone
    for (int k = pos; k + 1 < active_slot_count; ++k)
      active_slots[k] = active_slots[k + 1];
    --active_slot_count;
  }

  /**
   * @brief Allocates a free slot and schedules its animation on the timeline.
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param pin Whether the timeline event is pinned (retained-handle contract).
   * @param in_frames Delay in frames before the animation starts.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot or
   * timeline event is available.
   */
  template <typename... Args>
  AnimT *spawn_impl(Timeline::Pin pin, int in_frames, Args &&...args) {
    HS_CHECK(entities, "TransformerPool: call init_storage() before spawn");
    check_storage_alive();
    // Linear scan for a free slot (cold path).
    for (int idx = 0; idx < CAPACITY; ++idx) {
      Entity &e = entities[idx];
      if (!e.active) {
        e.params = template_params;
        e.active = true;
        add_active(idx);
        // Completion callback captures `this` + the slot index (both stable) to
        // deactivate the slot and drop it from the active list, plus the serial
        // is_live() needs to reject a pool destroyed before the callback fires.
        const uint32_t serial = pool_serial;
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        // The slot composes from here on, possibly before any prepare_frame(),
        // so derive its cached state from the fields the constructor seeded.
        if constexpr (HAS_SYNC) {
          e.params.sync();
        }
        AnimT *p = timeline.add_get(in_frames, std::move(anim), pin);
        if (p) {
          // A non-pinned spawn keeps no retained handle, so the slot is reclaimed
          // only by the one-shot then() below, which fires when the animation
          // reaches done() once. A finite, non-repeating animation is required: an
          // infinite one never reaches done() and a repeating one rewinds instead
          // of being removed, so either would hold its slot for the effect's life
          // (nullptr after CAPACITY spawns) — use spawn_pinned.
          if (pin == Timeline::Pin::UNPINNED)
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
          if (pin == Timeline::Pin::PINNED) {
            // Capture order sizes the callable: the two pointers lead so the
            // pair of 32-bit fields packs into Fn's inline storage.
            p->then([this, p, idx, serial]() {
              if (!is_live(this, serial))
                return;
              if (!p->repeats()) {
                entities[idx].active = false;
                remove_active(idx);
              }
            });
          } else {
            p->then([this, idx, serial]() {
              if (!is_live(this, serial))
                return;
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
   * ordering contract. Per-pixel hot path — no guard here by design.
   */
  HS_O3_FN Vector transform(Vector v) const {
    for (int k = 0; k < this->active_slot_count; ++k) {
      v = TransformFunc(v, this->entities[this->active_slots[k]].params);
    }
    return v;
  }

  /**
   * @brief Function-call alias for transform().
   * @param v Vector to transform.
   * @return The transformed vector.
   */
  HS_O3_FN Vector operator()(const Vector &v) const { return transform(v); }
};

/** @brief Denominator floor below which DominantFieldAccumulator::value()
 * reports 0 instead of dividing. */
constexpr float FIELD_DOMINANT_DEN_EPS = 1e-9f;

/**
 * @brief Accumulates a magnitude-weighted blend of scalar fields: the strongest
 * contribution dominates without stacking.
 * @details Use instead of summation when overlapping entities must not add
 * (e.g. solid bodies displacing a shared sheet). Unlike a hard max by
 * magnitude — which jumps discontinuously where opposite-signed fields cross in
 * strength — the blend is continuous everywhere: a single field passes through
 * exactly, equal same-signed overlaps yield the shared value, and
 * opposite-signed overlaps cancel smoothly.
 */
struct DominantFieldAccumulator {
  /** @brief Folds one field sample into the blend. */
  void add(float field) {
    numerator += field * field * field;
    denominator += field * field;
  }

  /** @brief The blend so far: sum(s_i^3) / sum(s_i^2); 0 with nothing added. */
  float value() const {
    return denominator > FIELD_DOMINANT_DEN_EPS ? numerator / denominator
                                                : 0.0f;
  }

private:
  float numerator = 0.0f;
  float denominator = 0.0f;
};

/**
 * @brief A generic manager for animation-driven scalar displacement fields.
 * @tparam ParamsT The configuration struct (e.g., BumpParams).
 * @tparam AnimT The animation class (e.g., Animation::BallDrop).
 * @tparam FieldFunc The static function evaluating one entity's field.
 * @tparam CAPACITY Max number of active fields.
 * @details The scalar counterpart of Transformer: entities compose as scalars
 * instead of as warps, so a caller can feed the composed field into a
 * displacement path (e.g. a DistortedRing shift LUT). field() offers
 * superposition by summation; a caller whose entities must not stack composes
 * its own way over active_count()/active_params() (see
 * DominantFieldAccumulator). field_bound() bounds either composition, since the
 * dominant blend never exceeds the largest contribution in magnitude.
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
    for (int k = 0; k < this->active_slot_count; ++k) {
      s += FieldFunc(p, this->entities[this->active_slots[k]].params);
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
   * @brief Upper bound on |field()| over the sphere this frame.
   * @return Sum of the active entities' per-entity bounds.
   * @details Requires ParamsT::field_bound() (a true upper bound on
   * |FieldFunc|); callers use it to size conservative culls.
   */
  float field_bound() const {
    float b = 0.0f;
    for (int k = 0; k < this->active_slot_count; ++k) {
      b += this->entities[this->active_slots[k]].params.field_bound();
    }
    return b;
  }
};

/**
 * @brief A transformer adapter for an Orientation object.
 * @tparam CAPACITY History capacity of the wrapped Orientation.
 */
template <int CAPACITY = 4> struct OrientTransformer {
  const Orientation<CAPACITY> &
      orientation; /**< Orientation applied by each transform; retained by reference. */

  /**
   * @brief Constructs an adapter wrapping an orientation.
   * @param ori Orientation to apply; retained by reference.
   */
  explicit OrientTransformer(const Orientation<CAPACITY> &ori)
      : orientation(ori) {}

  /**
   * @brief Deleted constructor from a temporary Orientation.
   * @details The adapter retains its argument by reference, so binding a
   * temporary would leave every later transform() reading a dead object.
   */
  explicit OrientTransformer(const Orientation<CAPACITY> &&) = delete;

  /**
   * @brief Orients a vector through the wrapped orientation.
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  HS_O3_FN Vector transform(const Vector &v) const {
    return orientation.orient(v);
  }

  /**
   * @brief Function-call alias for transform().
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  HS_O3_FN Vector operator()(const Vector &v) const { return transform(v); }
};

template <int CAPACITY>
OrientTransformer(const Orientation<CAPACITY> &) -> OrientTransformer<CAPACITY>;

/** @brief Largest ripple rotation the series-form quaternion may take; at
 * theta/2 <= 0.075 the truncated sin/cos series err under 3 float ulps. */
constexpr float RIPPLE_SMALL_ANGLE_MAX = 0.15f;

/**
 * @brief Applies one ripple wavelet at a known angular distance and phase.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @param distance Angular distance from the ripple center.
 * @param phase Angular position of the wavelet peak.
 * @return The displaced vector.
 */
HS_O3_FN inline Vector
ripple_transform_at_distance(const Vector &v,
                             const Animation::RippleParams &params,
                             float distance, float phase) {
  float dist_from_peak = distance - phase;
  float half_width = params.half_width();
  float t = (dist_from_peak / half_width) * 2.0f;
  float theta = params.amplitude * (1.0f - t * t) *
                fast_expf(-0.5f * t * t - params.decay * distance);

  Vector axis = cross(params.center, v);
  float len_sq = dot(axis, axis);
  if (len_sq > 1e-6f) {
    axis = axis * (1.0f / sqrtf(len_sq));
    if (fabsf(theta) <= RIPPLE_SMALL_ANGLE_MAX) {
      float h = 0.5f * theta;
      float h2 = h * h;
      float s = h * (1.0f - h2 * (1.0f / 6.0f));
      float c = 1.0f - h2 * (0.5f - h2 * (1.0f / 24.0f));
      return rotate(v, Quaternion(c, s * axis));
    }
    Quaternion q = make_rotation(axis, theta);
    return rotate(v, q);
  }

  return v;
}

/**
 * @brief Rotates a point along a Ricker-wavelet ripple radiating from a center.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @return The displaced vector.
 */
HS_O3_FN inline Vector ripple_transform(const Vector &v,
                                        const Animation::RippleParams &params) {
  // Between ripples the envelope drives amplitude to 0; skip the whole per-pixel
  // wavelet (fast_acos + fast_expf) when there is nothing to displace.
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
  return ripple_transform_at_distance(v, params, d, params.phase);
}

/**
 * @brief Applies a spatially periodic train of spherical ripple wavelets.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @param period Angular spacing between wavelet peaks.
 * @return The displaced vector.
 */
HS_O3_FN inline Vector
periodic_ripple_transform(const Vector &v,
                          const Animation::RippleParams &params, float period) {
  if (params.amplitude <= 0.001f)
    return v;

  float cos_d = dot(v, params.center);
  float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  float winding = floorf((d - params.phase) / period + 0.5f);
  float phase = params.phase + winding * period;
  if (fabsf(d - phase) > 2.0f * params.half_width())
    return v;
  return ripple_transform_at_distance(v, params, d, phase);
}
/**
 * @brief Slides a point along the sphere surface by a 3D-noise field.
 * @param v The unit vector to transform.
 * @param params Noise field, scale, amplitude and time.
 * @return The displaced unit vector.
 * @details Samples three decorrelated noise channels (the second and third
 * field-shifted by 100 and 200 on all three axes) to build a displacement,
 * projects it onto the tangent plane at v so the point stays on the sphere,
 * soft-caps the slide to avoid cross-hemisphere jumps, then renormalizes. No-op
 * when amplitude is negligible.
 */
inline Vector noise_transform(const Vector &v,
                              const Animation::NoiseParams &params) {
  if (params.amplitude <= 0.001f)
    return v;

  float scale = params.scale;
  float time_val = params.time;

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
 * @brief Computes the bump profile from prepared local cap geometry.
 * @param params Bump field geometry and gain.
 * @param r_eff Envelope-scaled footprint radius.
 * @param d Angular distance from the cap center, < @p r_eff.
 * @param y Signed polar offset from the cap center, |y| <= @p d.
 * @return The signed polar displacement (radians).
 * @details Both bounds are required: past |y| = r_eff the drape term's sine
 * turns negative while the depth term does too, so the profile grows with |y|
 * instead of decaying and overruns BumpParams::field_bound().
 */
inline float bump_field_profile(const Animation::BumpParams &params,
                                float r_eff, float d, float y) {
  float abs_y = std::fabs(y);
  float x_sq = std::max(d * d - y * y, 0.0f);
  float depth = sqrtf(std::max(r_eff * r_eff - x_sq, 0.0f)) - abs_y;
  float drape = std::min(params.amplitude * sinf(PI_F * abs_y / r_eff), 1.0f);
  return copysignf(depth * drape, y);
}

/**
 * @brief Tests a sample against the bump's effective cap.
 * @param v Sample point (unit vector).
 * @param params Bump field geometry and gain.
 * @param r_eff Receives the envelope-scaled footprint radius.
 * @param d Receives the angular distance from the cap center; only meaningful
 * on a true return.
 * @return Whether @p v lies inside the cap and the gain is non-negligible.
 */
__attribute__((always_inline)) inline bool
bump_cap_hit(const Vector &v, const Animation::BumpParams &params, float &r_eff,
             float &d) {
  r_eff = params.radius * params.envelope;
  if (r_eff <= 1e-3f || params.amplitude <= 0.001f)
    return false;
  float cos_d = dot(v, params.center);
  if (cos_d <= params.cos_radius)
    return false;
  d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  return d < r_eff;
}

/**
 * @brief Evaluates a bump using a caller-provided signed ring offset.
 * @param v Sample point (unit vector).
 * @param params Bump field geometry and gain.
 * @param y Signed polar offset of @p v from the bump center about
 * params.axis — the same quantity bump_field() derives itself, so it must agree
 * with @p v; |y| never exceeds the angular distance between them.
 * @return The signed polar displacement (radians).
 * @details For callers that already hold the offset (a ring stack sharing the
 * bump axis). A @p y disagreeing with @p v breaks bump_field_profile()'s bound
 * and with it BumpParams::field_bound().
 */
inline float bump_field_with_y(const Vector &v,
                               const Animation::BumpParams &params, float y) {
  float r_eff, d;
  if (!bump_cap_hit(v, params, r_eff, d))
    return 0.0f;

  return bump_field_profile(params, r_eff, d, y);
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
inline float bump_field(const Vector &v, const Animation::BumpParams &params) {
  float r_eff, d;
  if (!bump_cap_hit(v, params, r_eff, d))
    return 0.0f;

  // Local cap coords: signed polar offset y from the center (positive toward
  // larger colatitude) and azimuthal offset x, with x^2 = d^2 - y^2 by the
  // small-cap approximation. The boundary arc at this azimuth sits at
  // +-sqrt(r_eff^2 - x_sq), so (arc - |y|) is the polar depth inside the cap —
  // an arc-shaped profile along the ring, which keeps the bulge round.
  float y = fast_acos(hs::clamp(dot(params.axis, v), -1.0f, 1.0f)) -
            fast_acos(hs::clamp(dot(params.axis, params.center), -1.0f, 1.0f));
  return bump_field_profile(params, r_eff, d, y);
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
                                 const Animation::NoiseProductParams &params) {
  if (std::fabs(params.amplitude) <= 0.001f)
    return 0.0f;
  float n1 = params.noise.GetNoise(v.x * params.scale1, v.y * params.scale1,
                                   v.z * params.scale1 + params.time);
  float n2 = params.noise.GetNoise(
      v.x * params.scale2 + Animation::NoiseProductParams::OCTAVE2_OFFSET,
      v.y * params.scale2, v.z * params.scale2 + params.time);
  return params.amplitude * n1 * n2;
}

/**
 * @brief Generates ripples that warp the sphere.
 * @tparam CAPACITY Maximum number of concurrent ripple transformations.
 */
template <int CAPACITY>
using RippleTransformer =
    Transformer<Animation::RippleParams, Animation::Ripple, ripple_transform,
                CAPACITY>;

/**
 * @brief Bump displacement fields that fall pole-to-pole through a frame.
 * @tparam CAPACITY Maximum number of concurrent falling bumps.
 * @tparam ORIENT_CAP Sub-frame capacity of the orientation the bumps push
 * along.
 */
template <int CAPACITY, int ORIENT_CAP = 4>
using BallDropTransformer =
    FieldTransformer<Animation::BumpParams, Animation::BallDrop<ORIENT_CAP>,
                     bump_field, CAPACITY>;

/**
 * @brief A two-octave product noise displacement field.
 * @tparam CAPACITY Maximum number of concurrent noise fields.
 */
template <int CAPACITY>
using NoiseProductTransformer =
    FieldTransformer<Animation::NoiseProductParams, Animation::NoiseProduct,
                     noise_product_field, CAPACITY>;

/**
 * @brief Performs Mobius warps that return to the identity.
 * @tparam CAPACITY Maximum number of concurrent Mobius warp transformations.
 */
template <int CAPACITY>
using MobiusWarpTransformer = Transformer<MobiusParams, Animation::MobiusWarp,
                                          mobius_transform, CAPACITY>;

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
using NoiseTransformer = Transformer<Animation::NoiseParams, Animation::Noise,
                                     noise_transform, CAPACITY>;
