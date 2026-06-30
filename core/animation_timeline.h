/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "animation_core.h"
#include "animation_orientation.h"
#include "animation_scalars.h"
#include "animation_motion.h"
#include "animation_mesh.h"

/**
 * @brief Structure linking an animation with its starting time.
 * @details Stores the animation inline to avoid arena allocation (survives
 * compaction).
 */
struct TimelineEvent {
  // Inline storage budget for a type-erased animation: 112 B for device/WASM.
  // The native test build's std::function is wider, so HS_TEST_BUILD widens the
  // budget for that build only, leaving the device footprint unchanged.
#ifdef HS_TEST_BUILD
  static constexpr size_t MAX_ANIM_SIZE = 256;
#else
  static constexpr size_t MAX_ANIM_SIZE = 112;
#endif

  uint32_t start = 0; /**< Global frame at which the animation begins stepping. */
  /**
   * @brief Whether this event's inline animation pointer was handed out via
   * Timeline::add_get().
   * @details Compaction must never relocate such an event — doing so dangles the
   * caller's cached pointer. (Occupies existing alignment padding before
   * `storage`, so it costs no extra bytes.)
   */
  bool handled = false;
  alignas(std::max_align_t) uint8_t storage[MAX_ANIM_SIZE]; /**< Inline
                                  type-erased animation storage. */

  /**
   * @brief Type-erased move/destroy operation for the stored animation.
   * @details dst != nullptr → move src into dst and destroy src.
   *          dst == nullptr → just destroy src.
   */
  void (*manager)(TimelineEvent &src, TimelineEvent *dst) = nullptr;

  /**
   * @brief The animation's IAnimation* view, captured by static_cast at
   * construction (and re-captured by the manager on every move).
   * @details reinterpret_cast<IAnimation*> on the raw storage would be formally
   * UB: animation types are non-standard-layout (virtual functions), so the
   * standard does not guarantee the IAnimation base subobject sits at offset 0
   * — only a properly-typed upcast performs the (here-zero, but not
   * portably-zero) base adjustment. Trails `manager` so it lands in the existing
   * tail padding (no extra bytes).
   */
  IAnimation *iface = nullptr;

  /**
   * @brief Gets the live animation interface for this slot.
   * @return The IAnimation* view, or nullptr if the slot is empty.
   */
  IAnimation *animation() { return manager ? iface : nullptr; }

  /**
   * @brief Relocates this event's animation into a destination slot.
   * @param dst The destination slot to move into.
   */
  void move_into(TimelineEvent &dst) {
    // Relocating a handled event would dangle the caller's cached animation
    // pointer (handled animations are never meant to move); trap instead.
    HS_CHECK(!handled);
    dst.start = start;
    dst.handled = handled; // always false past the check; kept for symmetry
    dst.manager = manager;
    if (manager) {
      manager(*this, &dst);
      manager = nullptr;
    }
  }

  /**
   * @brief Destroys the stored animation and empties the slot.
   */
  void destroy() {
    if (manager) {
      manager(*this, nullptr);
      manager = nullptr;
    }
  }
};

// Device inline-storage budget audit. The per-type static_assert in
// Timeline::add only fires for types actually add()ed in this build, so pin the
// largest concrete animation type here so it is checked in every build that
// includes this header. Compared against MAX_ANIM_SIZE (not a literal 112) so it
// is the real 112 B device budget on the 32-bit WASM/device build and a no-op
// headroom check on the wider native host. (Templated animations stay covered at
// their add() sites.)
/** @brief Largest sizeof over a pack of types. */
template <typename... Ts> constexpr size_t largest_sizeof() {
  return std::max({sizeof(Ts)...});
}

// Add every new non-templated Animation type to this pack; the audit folds over
// it, so the static_assert tracks the list automatically.
constexpr size_t kLargestConcreteAnimSize = largest_sizeof<
    Animation::RandomTimer, Animation::PeriodicTimer, Animation::Transition,
    Animation::Mutation, Animation::Driver, Animation::Lerp, Animation::Sprite,
    Animation::ColorWipe, Animation::MobiusFlow, Animation::MobiusWarp,
    Animation::MobiusWarpCircular, Animation::MeshMorph,
    Animation::MobiusWarpEvolving, Animation::Ripple, Animation::Noise>();
static_assert(kLargestConcreteAnimSize <= TimelineEvent::MAX_ANIM_SIZE,
              "A concrete animation type exceeds the TimelineEvent inline-storage "
              "budget (on the 32-bit WASM/device build MAX_ANIM_SIZE is the "
              "112-byte device budget); shrink the type or raise the budget.");

/**
 * @brief Global storage for the timeline to prevent template instantiation
 * bloat.
 */
static constexpr int TIMELINE_MAX_EVENTS = 64;
extern DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
// True while a Timeline instance is alive (guards the single-singleton invariant).
extern bool global_timeline_live;
extern uint32_t global_timeline_t;     // current global frame count
extern int global_timeline_num_events; // current number of active events

/**
 * @brief Manages all active animations and their execution over time.
 *
 * Not templated: the entire instance is backed by one global event array
 * (`global_timeline_events`) sized to hold the largest effect, and the
 * live-guard permits exactly one instance at a time, so the capacity is a
 * single process-wide budget (`MAX_EVENTS`), not a per-instance knob.
 */
class Timeline {
public:
  /**
   * @brief Constructs a Timeline.
   *
   * Traps if one is already alive: all Timelines share global_timeline_events,
   * so two live instances would corrupt each other's events. The real app keeps
   * exactly one (destroys the old effect before building the next); this makes
   * the latent footgun a bench-time crash instead of silent stomping.
   */
  Timeline() {
    HS_CHECK(!global_timeline_live);
    global_timeline_live = true;
    clear();
  }

  /**
   * @brief Cleans up remaining animations, invoked on effect destruction.
   */
  ~Timeline() {
    clear();
    global_timeline_live = false;
  }

  // Singleton over global state — not copyable/movable; call clear() to reset.
  Timeline(const Timeline &) = delete;
  Timeline &operator=(const Timeline &) = delete;
  Timeline(Timeline &&) = delete;
  Timeline &operator=(Timeline &&) = delete;

  /**
   * @brief Destroys all events and resets the global frame counter.
   */
  void clear() {
    for (int i = 0; i < global_timeline_num_events; ++i) {
      global_timeline_events[i].destroy();
    }
    global_timeline_num_events = 0;
    global_timeline_t = 0;
  }

  /**
   * @brief Adds a new animation event to the timeline.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting.
   * @param animation The animation object.
   * @return Reference to the Timeline object.
   */
  template <typename A> Timeline &add(int in_frames, A animation) {
    // pin=false: an add() caller keeps no handle, so the event compacts normally.
    add_get(in_frames, std::move(animation), /*pin=*/false);
    return *this;
  }

  /**
   * @brief Like add(), but returns the typed pointer to the inline-stored
   * animation. Use when you need to hold a reference for later mutation.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting.
   * @param animation The animation object.
   * @param pin If true (default), the caller intends to RETAIN this pointer
   * across frames, so the event is marked handled and step()'s compaction traps
   * (move_into) rather than relocating it out from under the cached pointer.
   * Such a retained handle is only safe when the animation is infinite and added
   * before any finite one (so no earlier event is ever removed to shift it) —
   * the contract the direct callers rely on; the trap enforces it. Pass
   * `pin=false` for a TRANSIENT pointer used only at the call site and not kept
   * across frames (e.g. TransformerPool::spawn, whose finite animations are
   * compacted normally and whose return is typically discarded).
   * @return Typed pointer to the inline-stored animation, or nullptr if full.
   */
  template <typename A> A *add_get(int in_frames, A animation, bool pin = true) {
    static_assert(sizeof(A) <= TimelineEvent::MAX_ANIM_SIZE,
                  "Animation type exceeds TimelineEvent inline storage");
    static_assert(alignof(A) <= alignof(std::max_align_t),
                  "Animation type is over-aligned for TimelineEvent inline "
                  "storage (placement-new would be misaligned)");
    if (global_timeline_num_events >= MAX_EVENTS) {
      hs::log("Timeline full, failed to add animation!");
      return nullptr;
    }
    if (pin) {
      // A finite, non-repeating predecessor is removed on completion and would
      // relocate this pinned event (move_into traps). Repeating or infinite
      // predecessors are never removed, so they are safe to precede a pin.
      for (int i = 0; i < global_timeline_num_events; ++i) {
        IAnimation *prev = global_timeline_events[i].animation();
        HS_CHECK(!prev || !prev->is_finite() || prev->repeats(),
                 "pinned animation added after a finite non-repeating one");
      }
    }
    auto &e = global_timeline_events[global_timeline_num_events++];
    e.start = global_timeline_t + in_frames;
    e.handled = pin;
    auto *ptr = new (e.storage) A(std::move(animation));
    e.iface = static_cast<IAnimation *>(ptr);
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      // std::launder to recover a usable A* from the placement-new'd storage.
      A *obj = std::launder(reinterpret_cast<A *>(src.storage));
      if (dst) {
        dst->iface = static_cast<IAnimation *>(new (dst->storage) A(std::move(*obj)));
      }
      obj->~A();
    };
    return ptr;
  }

  /**
   * @brief Advances the timeline by one frame, stepping all active or starting
   * animations.
   * @param canvas The current canvas buffer.
   */
  void step(Canvas &canvas) {
    ++global_timeline_t;

    int write_idx = 0;
    int active_cnt = global_timeline_num_events; // Snapshot count before callbacks
                                 // potentially add more

    // Collapse each distinct Orientation exactly once, before any animation steps
    // it: collapsing per-animation would discard the sub-frame motion-blur history
    // a sharing animation already built this frame. Track the already-collapsed
    // ids in a fixed-size stack scratch set.
    const void *collapsed_ids[MAX_EVENTS];
    int collapsed_n = 0;
    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];
      if (global_timeline_t < e.start)
        continue;
      IAnimation *anim = e.animation();
      if (!anim)
        continue;
      const void *id = anim->orientation_id();
      if (!id)
        continue;
      bool already_collapsed = false;
      for (int j = 0; j < collapsed_n; ++j) {
        if (collapsed_ids[j] == id) {
          already_collapsed = true;
          break;
        }
      }
      if (!already_collapsed) {
        anim->collapse_orientation();
        collapsed_ids[collapsed_n++] = id;
      }
    }

    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];

      // 1. Check start time
      if (global_timeline_t < e.start) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
        continue;
      }

      // 2. Step (Orientation already collapsed once-per-frame above)
      IAnimation *anim = e.animation();
      HS_CHECK(anim);
      anim->step(canvas);

      // 4. Completion & Cleanup
      bool is_done = anim->done();
      bool keep = true;

      if (is_done) {
        bool does_repeat = anim->repeats();
        if (does_repeat) {
          anim->rewind();
          anim->post_callback();
        } else {
          keep = false;
          anim->post_callback();
        }
      }

      if (keep) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
      } else {
        // A pinned event should never reach natural completion (pinned ⇒
        // infinite); destroying one dangles the caller's pointer. cancel() is the
        // one sanctioned teardown, so it is exempt.
        HS_CHECK(!e.handled || anim->is_canceled());
        e.destroy();
      }
    }

    // 5. Move new events (added during callbacks) to fill the gap left by
    //    completed ones. A pinned event spawned inside a callback would trap in
    //    move_into here; callback-spawners use pin=false, so this is safe.
    int new_vals_count = global_timeline_num_events - active_cnt;
    // Each kept event advances write_idx by one, so it never outruns the events
    // scanned: the source span [active_cnt, ...) and the gap-fill destination
    // span [write_idx, ...) stay disjoint.
    HS_CHECK(write_idx <= active_cnt);
    if (new_vals_count > 0 && write_idx < active_cnt) {
      for (int i = 0; i < new_vals_count; ++i) {
        global_timeline_events[active_cnt + i].move_into(
            global_timeline_events[write_idx + i]);
      }
    }

    global_timeline_num_events = write_idx + new_vals_count;
  }

  /**
   * @brief Current global frame count (number of step() calls since clear()).
   * @return The shared timeline frame counter, advanced once per step().
   * @details Lets an effect derive a phase from the timeline's own clock rather
   *          than a parallel per-effect counter that can silently desync if it
   *          isn't advanced in exact lockstep with step().
   */
  uint32_t frame() const { return global_timeline_t; }

  static constexpr int MAX_EVENTS =
      TIMELINE_MAX_EVENTS; /**< Must match global_timeline_events array size. */
};

/**
 * @brief A double-buffered pair of persistent MeshState slots plus the
 *        arena-compaction primitives effects need to swap between them.
 * @details Holds two MeshState slots in `persistent_arena` and a front/back
 * index. Effects own the transition policy themselves — they generate into a
 * slot, schedule whatever animation they want (an opacity-crossfade `Sprite`, a
 * geometric `MeshMorph`, ...), flip the front index, and reclaim the old slot's
 * space. This class only provides the slot storage, index management, and the
 * two compaction helpers (`compact`, `compact_keep_front`) that those policies
 * share; it intentionally does not encode any single transition shape, because
 * the real clients (IslamicStars, HankinSolids, MeshFeedback) diverge on
 * animation type and on the per-slot side-band state they carry.
 *
 * Usage:
 *   MeshCarousel carousel;  // in effect members
 *
 *   // Build the initial shape directly into the front slot:
 *   carousel.current().clear();
 *   MeshOps::compile(mesh, carousel.current(), persistent_arena);
 *
 *   // To transition: generate into the back slot, schedule an animation,
 *   // then flip and compact (see IslamicStars::spawn_shape for the pattern).
 */
class MeshCarousel {
public:
  /**
   * @brief Constructs an empty carousel with front slot index 0.
   */
  MeshCarousel() {}

  /**
   * @brief Gets the currently visible (front) mesh.
   * @return Const reference to the front MeshState slot.
   */
  const MeshState &current() const { return slots_[front_]; }

  /**
   * @brief Gets the currently visible (front) mesh (mutable).
   * @return Mutable reference to the front MeshState slot.
   */
  MeshState &current() { return slots_[front_]; }

  /**
   * @brief Direct slot access by index (for effects that need both).
   * @param i Slot index (0 or 1).
   * @return Const reference to the requested MeshState slot.
   */
  const MeshState &slot(int i) const { return slots_[i]; }

  /**
   * @brief Direct slot access by index (mutable).
   * @param i Slot index (0 or 1).
   * @return Mutable reference to the requested MeshState slot.
   */
  MeshState &slot(int i) { return slots_[i]; }

  /**
   * @brief Gets the front slot index (for capture in lambdas).
   * @return The index (0 or 1) of the front slot.
   */
  int front_index() const { return front_; }

  /**
   * @brief Manually sets the front index (for effects that manage transitions
   * themselves).
   * @param idx The new front slot index (0 or 1).
   */
  void set_front(int idx) { front_ = idx; }

  /**
   * @brief Compacts the persistent arena, reclaiming fragmented space.
   * @details Evacuates tracked MeshStates and resets the arena. Call before
   * allocating new persistent data.
   */
  void compact() {
    // Both evacuations share scratch_arena_a, which must hold both populated slots.
    Persist<MeshState> p0(slots_[0], scratch_arena_a, persistent_arena);
    Persist<MeshState> p1(slots_[1], scratch_arena_a, persistent_arena);
    persistent_arena.reset();
  }

  /**
   * @brief Frees the back slot and compacts, preserving only the front slot.
   * @tparam AfterReset Callable type invoked as `void(Arena&)`.
   * @param after_reset Callback run immediately after the reset, while the front
   * slot is still evacuated.
   * @details Runs `after_reset(persistent_arena)` immediately after the reset —
   * while the front slot is still evacuated — so the caller can re-bake
   * effect-owned persistent data (e.g. a palette bank) into the fresh arena
   * *before* the front mesh is restored on top of it. Use when only the visible
   * (front) shape must survive a regeneration of the back slot.
   */
  template <typename AfterReset> void compact_keep_front(AfterReset after_reset) {
    int back = 1 - front_;
    slots_[back] = MeshState();
    Persist<MeshState> p(slots_[front_], scratch_arena_b, persistent_arena);
    persistent_arena.reset();
    after_reset(persistent_arena);
  }

private:
  MeshState slots_[2]; /**< Front/back double-buffered mesh slots. */
  int front_ = 0;      /**< Index (0 or 1) of the visible front slot. */
};

/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @tparam CAP Orientation sub-frame capacity.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&,
 * float t)`.
 */
template <int CAP>
void tween(const Orientation<CAP> &o, TweenFn callback) {
  int len = o.length();
  int start = (len > 1) ? 1 : 0;
  for (int i = start; i < len; ++i) {
    // A lone snapshot is the newest sub-position → t = 1 (age-neutral).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(o.get(i), t);
  }
}

/**
 * @brief Helper to iterate over a VectorTrail's historical frames.
 * @tparam CAPACITY Trail capacity.
 * @param trail The trail to iterate.
 * @param callback The function to call for each frame: `void(const Vector&,
 * float t)`.
 */
template <int CAPACITY>
void tween(const Animation::VectorTrail<CAPACITY> &trail,
           VectorTweenFn callback) {
  size_t len = trail.length();
  if (len == 0)
    return;

  for (size_t i = 0; i < len; ++i) {
    // A lone sample is the newest (head) position → t = 1 (age-neutral).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(trail.get(i), t);
  }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or
 * Animation::OrientationTrail).
 * @param trail The Tweenable object to iterate.
 * @param callback The function to call for each step: `void(const T&, float
 * t)`.
 */
void deep_tween(const Tweenable auto &trail, TweenFn callback) {
  size_t trail_len = trail.length();
  if (trail_len == 0)
    return;

  // The age ramp must end at 1.0 on the newest *plotted* orientation. A motionless
  // tail frame collapses to a single sub-frame that emits nothing, so normalizing
  // by trail_len would strand the ramp below 1.0; normalize against the newest
  // frame that actually contributes instead.
  size_t last = trail_len - 1;
  while (last > 0 && trail.get(last).length() <= 1)
    --last;

  // A fully motionless trail collapses to frame 0's lone sub-position, which must
  // read t = 1.0 (age-neutral) or quintic_kernel(0) renders a static head invisible.
  if (last == 0 && trail.get(0).length() == 1) {
    const auto &frame = trail.get(0);
    callback(frame.get(0), 1.0f);
    return;
  }

  float span = static_cast<float>(last + 1);

  for (size_t i = 0; i <= last; ++i) {
    const auto &frame = trail.get(i);
    size_t frame_size = frame.length();
    size_t start_j = (i == 0) ? 0 : 1;

    for (size_t j = start_j; j < frame_size; ++j) {
      const auto &q = frame.get(j);
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      float global_t = (static_cast<float>(i) + sub_t) / span;
      callback(q, global_t);
    }
  }
}
