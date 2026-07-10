/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

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
      // relocate this pinned event, so reject it up front. A repeating/infinite
      // predecessor can still be removed if cancel()ed later; the move_into
      // HS_CHECK is the real backstop for that case.
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

      // Check start time
      if (global_timeline_t < e.start) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
        continue;
      }

      // Step (Orientation already collapsed once-per-frame above)
      IAnimation *anim = e.animation();
      HS_CHECK(anim);
      anim->step(canvas);

      // Completion & Cleanup
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

    // Move new events (added during callbacks) to fill the gap left by
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
