/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file transformer.h
 * @brief TransformerPool, its Transformer and FieldTransformer derived pools,
 *        and the standalone OrientTransformer adapter.
 * @details Also holds the free warp and field functions those pools compose
 *          over Animation params.
 */

#include "animation/orientation.h"
#include "math/3dmath.h"
#include "math/mobius.h"
#include "engine/concepts.h"
#include "memory.h"
#include <new>
#include <type_traits>
#include "vendor/FastNoiseLite.h"
#include "animation/animation.h"

namespace transformer_detail {

/**
 * @brief Out-of-class declaration of a params type's prepare_frame() hook
 *        intent.
 * @tparam T Params type.
 * @details Params types owned by this layer declare `NEEDS_REFRESH_FROM` and
 *          `NEEDS_SYNC` as members. A params type owned by a layer below
 *          animation specializes this instead. A specialization must declare
 *          both constants.
 */
template <typename T> struct ExternalParamsHooks {};

/** @brief MobiusParams (math/mobius.h) carries neither hook. */
template <> struct ExternalParamsHooks<math::MobiusParams> {
  static constexpr bool NEEDS_REFRESH_FROM = false;
  static constexpr bool NEEDS_SYNC = false;
};

/**
 * @brief Whether a params type declares its hook intent at all.
 * @tparam T Candidate params type.
 */
template <typename T> constexpr bool declares_hooks() {
  constexpr bool REFRESH_DECLARED = requires {
    T::NEEDS_REFRESH_FROM;
  } || requires { ExternalParamsHooks<T>::NEEDS_REFRESH_FROM; };
  constexpr bool SYNC_DECLARED = requires { T::NEEDS_SYNC; } || requires {
    ExternalParamsHooks<T>::NEEDS_SYNC;
  };
  return REFRESH_DECLARED && SYNC_DECLARED;
}

/**
 * @brief A params type's declared refresh_from-hook intent, false when
 *        undeclared.
 * @tparam T Candidate params type.
 */
template <typename T> constexpr bool declared_needs_refresh_from() {
  if constexpr (requires { T::NEEDS_REFRESH_FROM; })
    return T::NEEDS_REFRESH_FROM;
  else if constexpr (requires { ExternalParamsHooks<T>::NEEDS_REFRESH_FROM; })
    return ExternalParamsHooks<T>::NEEDS_REFRESH_FROM;
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
  else if constexpr (requires { ExternalParamsHooks<T>::NEEDS_SYNC; })
    return ExternalParamsHooks<T>::NEEDS_SYNC;
  else
    return false;
}

} // namespace transformer_detail

/**
 * @brief Fixed-capacity pool of animation-driven parameter entities.
 * @tparam ParamsT The configuration struct (e.g., RippleParams, MobiusParams).
 * @tparam AnimT The animation class (e.g., Animation::Ripple).
 * @tparam CAPACITY Max number of active entities.
 * @details Owns the entity slots, their timeline lifecycle, and the per-frame
 * param refresh; derived classes compose over the active entities. Each pool
 * claims one Timeline clear-hook slot at init_storage().
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
      transformer_detail::declares_hooks<ParamsT>(),
      "ParamsT must declare `static constexpr bool NEEDS_REFRESH_FROM` and "
      "`NEEDS_SYNC` as members - true if it carries the matching "
      "prepare_frame() hook, false if it does not - or, when ParamsT belongs to "
      "a layer below animation, specialize "
      "transformer_detail::ExternalParamsHooks for it.");
  static_assert(
      transformer_detail::declared_needs_refresh_from<ParamsT>() ==
          HAS_REFRESH_FROM,
      "ParamsT::NEEDS_REFRESH_FROM disagrees with whether ParamsT exposes "
      "refresh_from(). prepare_frame() finds each hook by detection "
      "independently, so a renamed hook or a drifted signature would leave the "
      "entity unrefreshed with no other signal.");
  static_assert(
      transformer_detail::declared_needs_sync<ParamsT>() == HAS_SYNC,
      "ParamsT::NEEDS_SYNC disagrees with whether ParamsT exposes sync(). "
      "prepare_frame() finds each hook by detection independently, so a "
      "renamed hook or a drifted signature would leave the entity's derived "
      "state stale with no other signal.");

  ParamsT
      template_params; /**< Template params copied into each new entity on spawn. */

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
   * @brief Cancels owned events and drops the pool's liveness record.
   * @pre Later pinned animation owners must retire first.
   * @details Spawned animations' completion callbacks can outlive the pool;
   * dropping the liveness record turns them into no-ops.
   */
  HS_COLD_MEMBER ~TransformerPool() {
    unlink_live();
    if (clear_hook_registered) {
      HS_CHECK(global_timeline_live,
               "TransformerPool outlived its Timeline: declare the Timeline "
               "before the pools that schedule on it");
      timeline.cancel_owner(this);
      timeline.remove_clear_hook(this);
    }
  }

  // Completion callbacks capture this and a slot index; the pool must not move.
  TransformerPool(const TransformerPool &) = delete;
  TransformerPool(TransformerPool &&) = delete;

  /**
   * @brief Allocates the entity pool from the persistent arena.
   * @param arena Persistent arena supplying CAPACITY entity slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet), after any configure_arenas() and before the first spawn.
   * Registers one Timeline clear hook.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    HS_CHECK(!entities, "TransformerPool: init_storage() called twice");
    entities = arena.make_n<Entity>(CAPACITY);
    active_slots = arena.allocate_n<int>(CAPACITY);
    active_slot_count = 0;
    storage_arena = &arena;
    storage_end = arena.get_offset();
#ifndef NDEBUG
    stamp.record(arena);
#endif
    timeline.add_clear_hook(this, [](void *self) {
      static_cast<TransformerPool *>(self)->release_all();
    });
    clear_hook_registered = true;
  }

  /**
   * @brief Re-claims the pool's storage after its arena was reset, preserving
   * live entities.
   * @param arena The arena init_storage() allocated from, freshly reset.
   * @details For arenas that are compacted mid-effect. Spawned animations hold
   * Params references into the slots, so the caller must replay the same
   * allocation order after the reset as after init_storage() (asserted). The
   * bytes are left untouched, so live entities carry through.
   */
  HS_COLD_MEMBER void reclaim_storage(Arena &arena) {
    HS_CHECK(entities,
             "TransformerPool: call init_storage() before reclaim_storage");
    Entity *e = arena.allocate_n<Entity>(CAPACITY);
    int *s = arena.allocate_n<int>(CAPACITY);
    HS_CHECK(e == entities && s == active_slots,
             "TransformerPool: reclaimed storage moved");
    storage_end = arena.get_offset();
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
   * @param in_frames Delay before the first animation step. Constructor-seeded
   * params compose during the delay; delayed effects must seed a neutral state.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot or
   * timeline event is available.
   * @details Pin::UNPINNED: the returned pointer is transient (used at the call
   * site, not retained across frames).
   * The pool claims the animation's single then() slot to recycle the entity,
   * so the caller must not attach one (Animation::then() traps on a second).
   */
  template <typename... Args> AnimT *spawn(int in_frames, Args &&...args) {
    return spawn_impl(Timeline::Pin::UNPINNED, nullptr, in_frames,
                      std::forward<Args>(args)...);
  }

  /**
   * @brief Like spawn(), but gates the timeline event on a pause flag.
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param paused Pause flag that must outlive the event.
   * @param in_frames Delay before the first animation step. Constructor-seeded
   * params compose during the delay; delayed effects must seed a neutral state.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot or
   * timeline event is available.
   * @details Freezes the whole event, delay included, the way
   * Timeline::add_pausable does.
   */
  template <typename... Args>
  AnimT *spawn_pausable(const bool *paused, int in_frames, Args &&...args) {
    HS_CHECK(paused != nullptr, "pausable spawn needs a pause flag");
    return spawn_impl(Timeline::Pin::UNPINNED, paused, in_frames,
                      std::forward<Args>(args)...);
  }

  /**
   * @brief Like spawn(), but pins the event so the returned pointer may be
   * retained across frames (e.g. registered as a live GUI param).
   * @tparam Args Constructor argument types forwarded to the Animation.
   * @param in_frames Delay before the first animation step. Constructor-seeded
   * params compose during the delay; delayed effects must seed a neutral state.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot is
   * available. A full timeline traps.
   * @details Only valid when the spawned animation is infinite or repeating and
   * is added before any finite, non-repeating timeline event (see
   * Timeline::add_get). The pool claims the animation's single then() slot to
   * recycle the entity, so the retained handle must not attach one
   * (Animation::then() traps on a second).
   */
  template <typename... Args>
  AnimT *spawn_pinned(int in_frames, Args &&...args) {
    return spawn_impl(Timeline::Pin::PINNED, nullptr, in_frames,
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
    HS_CHECK(entities,
             "TransformerPool: call init_storage() before prepare_frame");
    check_storage_watermark();
    check_storage_alive();
    for (int k = 0; k < active_slot_count; ++k) {
      Entity &e = entities[active_slots[k]];
      // Pull live-tunable config from template_params into the spawned entity.
      if constexpr (HAS_REFRESH_FROM) {
        e.params.refresh_from(template_params);
      }
      if constexpr (HAS_SYNC) {
        e.params.sync();
      }
    }
  }

protected:
  /**
   * @brief Compact list of the active slots, in spawn order.
   * @details Held in spawn order: the warps do not all commute, so composition
   * order must follow spawn order, not slot recycling.
   */
  int *active_slots = nullptr;
  int active_slot_count =
      0; /**< Number of valid entries at the front of active_slots. */

  Entity *entities =
      nullptr; /**< CAPACITY-slot pool, allocated by init_storage(). */

private:
  Timeline &
      timeline; /**< Timeline that schedules and steps the spawned animations. */
  bool clear_hook_registered =
      false; /**< Whether init_storage() registered the timeline clear hook. */
  Arena *storage_arena = nullptr; /**< Arena init_storage() allocated from. */
  size_t storage_end = 0; /**< Arena offset just past the pool's blocks. */

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
   * @details The serial separates a destroyed pool from a later one built at
   * the same address.
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

  /**
   * @brief Traps if the pool's arena was reclaimed under the live slots.
   * @details Detects reclamation only while the arena offset remains below
   * the storage watermark; intervening allocations can mask it.
   */
  HS_COLD_MEMBER void check_storage_watermark() const {
    HS_CHECK(storage_arena->get_offset() >= storage_end,
             "TransformerPool: arena reclaimed under a live pool; "
             "init_storage() runs after configure_arenas()");
  }

#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state at the last
                              init_storage()/reclaim_storage(). */
#endif

  /**
   * @brief Debug-only use-after-free check on the pool's arena blocks.
   * @details Asserts if the arena was reset, rebound, or rewound below either
   * block since the last init_storage()/reclaim_storage().
   */
  void check_storage_alive() const {
    HS_ASSERT_BLOCK_ALIVE(stamp, entities, CAPACITY * sizeof(Entity),
                          "TransformerPool entities");
    HS_ASSERT_BLOCK_ALIVE(stamp, active_slots, CAPACITY * sizeof(int),
                          "TransformerPool active slots");
  }

  /**
   * @brief Frees every slot without touching the timeline.
   * @details The slots' animations must not step afterwards. Timeline::clear()
   * destroys them without completion callbacks after this hook returns.
   */
  HS_COLD_MEMBER void release_all() {
    HS_CHECK(entities,
             "TransformerPool: call init_storage() before release_all");
    check_storage_watermark();
    check_storage_alive();
    for (int i = 0; i < CAPACITY; ++i)
      entities[i].active = false;
    active_slot_count = 0;
  }

  /**
   * @brief Appends a slot index to active_slots in spawn order.
   * @param idx Slot index to append.
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
   * @param paused Event-level pause gate, or nullptr for none.
   * @param in_frames Delay before the first animation step. Constructor-seeded
   * params compose during the delay; delayed effects must seed a neutral state.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   * @return Pointer to the spawned animation, or nullptr if no pool slot or
   * timeline event is available.
   */
  template <typename... Args>
  AnimT *spawn_impl(Timeline::Pin pin, const bool *paused, int in_frames,
                    Args &&...args) {
    HS_CHECK(entities, "TransformerPool: call init_storage() before spawn");
    check_storage_watermark();
    check_storage_alive();
    // Linear scan for a free slot (cold path).
    for (int idx = 0; idx < CAPACITY; ++idx) {
      Entity &e = entities[idx];
      if (!e.active) {
        e.params = template_params;
        e.active = true;
        add_active(idx);
        // The serial lets the completion callback reject a destroyed pool.
        const uint32_t serial = pool_serial;
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        // The slot composes from here on, possibly before any prepare_frame(),
        // so derive its cached state from the fields the constructor seeded.
        if constexpr (HAS_SYNC) {
          e.params.sync();
        }
        HS_CHECK(pin == Timeline::Pin::PINNED ||
                     (anim.is_finite() && !anim.repeats()),
                 "Transformer::spawn needs a finite, non-repeating animation; "
                 "use spawn_pinned for infinite or repeating animations");
        AnimT *p =
            timeline.add_get(in_frames, std::move(anim), pin, paused, this);
        if (p) {
          // Pinned callbacks re-query live repeats(); non-pinned pointers expire
          // when step() compacts events.
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
          // Timeline full: undo the activation so the slot is not leaked.
          e.active = false;
          remove_active(idx);
        }
        return p;
      }
    }
    // No free slot: drop the spawn.
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
          math::Vector (*TransformFunc)(const math::Vector &, const ParamsT &),
          int CAPACITY = 32>
class Transformer : public TransformerPool<ParamsT, AnimT, CAPACITY> {
public:
  using TransformerPool<ParamsT, AnimT, CAPACITY>::TransformerPool;

  /**
   * @brief Applies all active transformations to a vector, in spawn order.
   * @param v Vector to transform.
   * @return The vector after every active transform has been composed onto it.
   * @note Reads each active entity's prepared state; see prepare_frame() for the
   * ordering contract.
   */
  HS_O3_FN math::Vector transform(math::Vector v) const {
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
  HS_O3_FN math::Vector operator()(const math::Vector &v) const {
    return transform(v);
  }
};

/** @brief Denominator floor below which DominantFieldAccumulator::value()
 * reports 0 instead of dividing. */
constexpr float FIELD_DOMINANT_DEN_EPS = 1e-9f;

/**
 * @brief Accumulates a magnitude-weighted blend of scalar fields: the strongest
 * contribution dominates without stacking.
 * @details Use instead of summation when overlapping entities must not add.
 * At or below FIELD_DOMINANT_DEN_EPS the result is zero, a discontinuity.
 */
struct DominantFieldAccumulator {
  /** @brief Folds one field sample into the blend. */
  void add(float field) { accumulate(numerator, denominator, field); }

  /** @brief Adds one sample to caller-owned blend terms. */
  static void accumulate(float &num, float &den, float field) {
    num += field * field * field;
    den += field * field;
  }

  /** @brief The blend so far: sum(s_i^3) / sum(s_i^2); 0 with nothing added. */
  float value() const { return resolve(numerator, denominator); }

  /** @brief Resolves caller-owned blend terms. */
  static float resolve(float num, float den) {
    return den > FIELD_DOMINANT_DEN_EPS ? num / den : 0.0f;
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
 * @details field() sums the entities; a caller whose entities must not stack
 * composes over active_count()/active_params() (see DominantFieldAccumulator).
 * field_bound() bounds either composition.
 */
template <typename ParamsT, typename AnimT,
          float (*FieldFunc)(const math::Vector &, const ParamsT &),
          int CAPACITY = 32>
class FieldTransformer : public TransformerPool<ParamsT, AnimT, CAPACITY> {
public:
  using TransformerPool<ParamsT, AnimT, CAPACITY>::TransformerPool;

  /**
   * @brief Sums every active entity's field at a point.
   * @param p Sample point (unit vector).
   * @return The superposed field value; 0 with no active entities.
   * @note Reads each active entity's prepared state; see prepare_frame() for the
   * ordering contract.
   */
  float field(const math::Vector &p) const {
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
  float operator()(const math::Vector &p) const { return field(p); }

  /**
   * @brief Upper bound on |field()| over the sphere this frame.
   * @return Sum of the active entities' per-entity bounds.
   * @details Requires ParamsT::field_bound() (a true upper bound on
   * |FieldFunc|).
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
  const math::Orientation<CAPACITY> &
      orientation; /**< Orientation applied by each transform; retained by reference. */

  /**
   * @brief Constructs an adapter wrapping an orientation.
   * @param ori Orientation to apply; retained by reference.
   */
  explicit OrientTransformer(const math::Orientation<CAPACITY> &ori)
      : orientation(ori) {}

  /**
   * @brief Deleted constructor from a temporary Orientation.
   * @details The adapter retains its argument by reference, so binding a
   * temporary would leave every later transform() reading a dead object.
   */
  explicit OrientTransformer(const math::Orientation<CAPACITY> &&) = delete;

  /**
   * @brief Orients a vector through the wrapped orientation.
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  HS_O3_FN math::Vector transform(const math::Vector &v) const {
    return orientation.orient(v);
  }

  /**
   * @brief Function-call alias for transform().
   * @param v Vector to transform.
   * @return The oriented vector.
   */
  HS_O3_FN math::Vector operator()(const math::Vector &v) const {
    return transform(v);
  }
};

template <int CAPACITY>
OrientTransformer(const math::Orientation<CAPACITY> &)
    -> OrientTransformer<CAPACITY>;

/** @brief Largest ripple rotation the series-form quaternion may take, within
 * which the truncated sin/cos series stay at float rounding. */
constexpr float RIPPLE_SMALL_ANGLE_MAX = 0.15f;

/**
 * @brief Applies one ripple wavelet at a known angular distance and phase.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @param distance Angular distance from the ripple center.
 * @param phase Angular position of the wavelet peak.
 * @return The displaced vector.
 */
HS_O3_FN inline math::Vector
ripple_transform_at_distance(const math::Vector &v,
                             const Animation::RippleParams &params,
                             float distance, float phase) {
  float dist_from_peak = distance - phase;
  float half_width = params.half_width();
  float t = (dist_from_peak / half_width) * 2.0f;
  float theta = params.amplitude * (1.0f - t * t) *
                math::fast_expf(-0.5f * t * t - params.decay * distance);

  math::Vector axis = math::cross(params.center, v);
  float len_sq = math::dot(axis, axis);
  if (len_sq > 1e-6f) {
    axis = axis * (1.0f / sqrtf(len_sq));
    if (fabsf(theta) <= RIPPLE_SMALL_ANGLE_MAX) {
      float h = 0.5f * theta;
      float h2 = h * h;
      float s = h * (1.0f - h2 * (1.0f / 6.0f));
      float c = 1.0f - h2 * (0.5f - h2 * (1.0f / 24.0f));
      return math::rotate(v, math::Quaternion(c, s * axis));
    }
    math::Quaternion q = math::make_rotation(axis, theta);
    return math::rotate(v, q);
  }

  return v;
}

/**
 * @brief Rotates a point along a Ricker-wavelet ripple radiating from a center.
 * @param v The vector to transform.
 * @param params The ripple parameters.
 * @return The displaced vector.
 */
HS_O3_FN inline math::Vector
ripple_transform(const math::Vector &v, const Animation::RippleParams &params) {
  // Between ripples the envelope drives amplitude to 0.
  if (params.amplitude <= 0.001f)
    return v;

  // Reject outside the [d_min, d_max] angular band; see RippleParams.
  float cos_d = math::dot(v, params.center);
  if (cos_d > params.cos_threshold_min || cos_d < params.cos_threshold_max) {
    return v;
  }

  float d = math::fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  return ripple_transform_at_distance(v, params, d, params.phase);
}

/**
 * @brief Slides a point along the sphere surface by a 3D-noise field.
 * @param v The unit vector to transform.
 * @param params Noise field, scale, amplitude and time.
 * @return The displaced unit vector.
 * @details Projects a three-channel noise displacement onto the tangent plane
 * at v, soft-caps the slide, and renormalizes.
 */
inline math::Vector noise_transform(const math::Vector &v,
                                    const Animation::NoiseParams &params) {
  if (params.amplitude <= 0.001f)
    return v;
  HS_AUDIT_CHECK(std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z),
                 "noise_transform: non-finite direction");

  float scale = params.scale;
  float time_val = params.time;

  // Constant spatial shifts decorrelate channels; all share the same time input.
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

  math::Vector raw_noise =
      math::Vector(nx, ny, nz) * (params.amplitude * 0.05f);

  float inward_pull = math::dot(raw_noise, v);
  math::Vector surface_distortion = raw_noise - (v * inward_pull);

  // Soft-cap the slide distance to prevent cross-hemisphere grabs.
  constexpr float MAX_SLIDE = 0.5f;
  float sd_len_sq = math::dot(surface_distortion, surface_distortion);
  if (sd_len_sq > MAX_SLIDE * MAX_SLIDE) {
    surface_distortion = surface_distortion * (MAX_SLIDE / sqrtf(sd_len_sq));
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
  float x_sq = fmaxf(d * d - y * y, 0.0f);
  float depth = sqrtf(fmaxf(r_eff * r_eff - x_sq, 0.0f)) - abs_y;
  float drape =
      fminf(params.amplitude * sinf(math::PI_F * abs_y / r_eff), 1.0f);
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
bump_cap_hit(const math::Vector &v, const Animation::BumpParams &params,
             float &r_eff, float &d) {
  r_eff = params.radius * params.envelope;
  if (r_eff <= 1e-3f || params.amplitude <= 0.001f)
    return false;
  float cos_d = math::dot(v, params.center);
  if (cos_d <= params.cos_radius)
    return false;
  d = math::fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
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
inline float bump_field_with_y(const math::Vector &v,
                               const Animation::BumpParams &params, float y) {
  float r_eff, d;
  if (!bump_cap_hit(v, params, r_eff, d))
    return 0.0f;

  // fast_acos and near-unit dot-product rounding perturb angular distances.
  HS_AUDIT_CHECK(std::abs(y) <= d + 1e-3f,
                 "bump offset exceeds the angular distance to its center");

  return bump_field_profile(params, r_eff, d, y);
}

/**
 * @brief Evaluates a spherical-cap drape push: rings bow away from the cap
 * center as if draping over a ball beneath them.
 * @param v Sample point (unit vector).
 * @param params Bump center, stack axis, footprint and lifecycle envelope
 * (the envelope scales the effective footprint, inflating/deflating the cap).
 * @return The signed polar displacement (radians): the depth inside the cap's
 * boundary arc, weighted by a drape factor that is zero at the center ring and
 * the footprint edge; the amplitude gain scales that weight, saturating at 1.
 * Positive pushes toward larger colatitude about the axis; 0 outside the cap.
 */
inline float bump_field(const math::Vector &v,
                        const Animation::BumpParams &params) {
  float r_eff, d;
  if (!bump_cap_hit(v, params, r_eff, d))
    return 0.0f;

  // Signed polar offset from the center, positive toward larger colatitude.
  float y = math::fast_acos(hs::clamp(math::dot(params.axis, v), -1.0f, 1.0f)) -
            math::fast_acos(
                hs::clamp(math::dot(params.axis, params.center), -1.0f, 1.0f));
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
inline float noise_product_field(const math::Vector &v,
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
 * @note Spawn through spawn_pinned(): Animation::NoiseProduct is perpetual,
 * which spawn()/spawn_pausable() reject.
 */
template <int CAPACITY>
using NoiseProductTransformer =
    FieldTransformer<Animation::NoiseProductParams, Animation::NoiseProduct,
                     noise_product_field, CAPACITY>;

/**
 * @brief Performs Mobius warps that return to the identity.
 * @tparam CAPACITY Maximum number of concurrent Mobius warp transformations.
 * @note Repeats by default: use spawn_pinned(), or pass repeat=false to
 * spawn()/spawn_pausable().
 */
template <int CAPACITY>
using MobiusWarpTransformer =
    Transformer<math::MobiusParams, Animation::MobiusWarp,
                math::mobius_transform, CAPACITY>;

/**
 * @brief Performs circular Mobius warps at constant strength, suitable
 * for repeating animations.
 * @tparam CAPACITY Maximum number of concurrent circular Mobius warps.
 * @warning With nonzero scale this never returns to identity:
 * `Animation::MobiusWarpCircular` traces a closed loop at full strength. Use it
 * only in a repeating slot; a non-repeating slot freezes off-identity on its
 * final frame. Use MobiusWarpTransformer for one-shot slots.
 * @note Spawn through spawn_pinned(): spawn()/spawn_pausable() reject a
 * repeating animation. A zero-scale warp remains the identity.
 */
template <int CAPACITY>
using MobiusWarpCircularTransformer =
    Transformer<math::MobiusParams, Animation::MobiusWarpCircular,
                math::mobius_transform, CAPACITY>;

/**
 * @brief Performs a changing Mobius warp using gnomonic projection.
 * @tparam CAPACITY Maximum number of concurrent gnomonic Mobius warps.
 * @note Spawn through spawn_pinned(): Animation::MobiusWarpEvolving is
 * perpetual, which spawn()/spawn_pausable() reject.
 */
template <int CAPACITY>
using MobiusWarpGnomonicTransformer =
    Transformer<math::MobiusParams, Animation::MobiusWarpEvolving,
                math::gnomonic_mobius_transform, CAPACITY>;

/**
 * @brief Applies 3D noise distortion to vectors.
 * @tparam CAPACITY Maximum number of concurrent noise transformations.
 * @note Animation::Noise defaults to an indefinite duration, which
 * spawn()/spawn_pausable() reject; pass a finite duration or spawn_pinned().
 */
template <int CAPACITY>
using NoiseTransformer = Transformer<Animation::NoiseParams, Animation::Noise,
                                     noise_transform, CAPACITY>;
