/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file timeline.h
 * @brief Animation fragment: TimelineEvent inline storage and the Timeline
 * scheduler.
 */

/**
 * @brief Structure linking an animation with its starting time.
 * @details Stores the animation inline; no arena allocation.
 */
struct TimelineEvent {
  // Inline storage budget for a type-erased animation: 112 B for device/WASM,
  // 256 B for 64-bit host effect harnesses.
  static constexpr size_t MAX_ANIM_SIZE = HS_TIMELINE_MAX_ANIM_BYTES;

  uint32_t start =
      0; /**< Global frame at which the animation becomes eligible to step. */
  /**
   * @brief Whether the caller retains this event's animation pointer across
   * frames via Timeline::add_get(..., Pin::PINNED).
   * @details Compaction must never relocate such an event — doing so dangles the
   * caller's cached pointer.
   */
  bool pinned = false;
  const bool *paused = nullptr; /**< Optional event-level pause gate. */
  const void *owner = nullptr;  /**< Optional lifetime owner. */
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
   * @details The IAnimation base need not sit at offset 0 of the storage, so a
   * reinterpret_cast of the storage would be UB.
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
    // Relocating a pinned event would dangle the caller's cached pointer.
    HS_CHECK(!pinned, "move_into would dangle a pinned animation's retained "
                      "pointer");
    HS_CHECK(!dst.manager,
             "move_into would leak the destination's live animation");
    dst.start = start;
    // Clears a stale flag: destroy() leaves pinned set on a canceled pinned
    // event, and this slot may be recycling one.
    dst.pinned = false;
    dst.paused = paused;
    dst.owner = owner;
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

/** @brief Event capacity of the process-wide timeline. */
inline constexpr int TIMELINE_MAX_EVENTS = 64;
/** @brief Process-wide timeline storage shared by all template instances. */
extern DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
// True while a Timeline instance is alive.
extern bool global_timeline_live;
extern uint32_t global_timeline_t;       // current global frame count
extern int global_timeline_num_events;   // current number of active events
extern uint32_t global_timeline_dropped; // wrapping full-timeline drop count
// Set once per saturation episode, cleared whenever the event table empties.
extern bool global_timeline_drop_logged;

/**
 * @brief Manages all active animations and their execution over time.
 *
 * Backed by the process-wide `global_timeline_events` array; at most one
 * instance may be live, and `MAX_EVENTS` is a process-wide budget.
 */
class Timeline {
public:
  /**
   * @brief Whether an add_get() caller retains the returned pointer across
   * frames (see add_get()).
   */
  enum class Pin { UNPINNED, PINNED };

  /**
   * @brief Constructs a Timeline.
   *
   * Traps if one is already alive: all Timelines share global_timeline_events.
   */
  Timeline() {
    HS_CHECK(!global_timeline_live,
             "a second live Timeline would stomp the shared global events");
    global_timeline_live = true;
    reset_storage();
  }

  /**
   * @brief Cleans up remaining animations, invoked on effect destruction.
   */
  ~Timeline() {
    reset_storage();
    global_timeline_live = false;
  }

  // Singleton over global state — not copyable/movable; call clear() to reset.
  Timeline(const Timeline &) = delete;
  Timeline &operator=(const Timeline &) = delete;
  Timeline(Timeline &&) = delete;
  Timeline &operator=(Timeline &&) = delete;

  /**
   * @brief Destroys all events, leaving the timeline empty and reusable.
   * @details Traps if called from inside step() (it would free a running
   * completion callback) or while a pinned event (add_get(Pin::PINNED)) is
   * live. Runs the clear hooks first. Does not rewind the global frame cursor;
   * only construction and destruction do.
   */
  void clear() {
    HS_CHECK(!stepping, "clear() from inside step() would destroy the "
                        "animation whose callback is running");
    for (int i = 0; i < global_timeline_num_events; ++i) {
      HS_CHECK(!global_timeline_events[i].pinned,
               "clear() would destroy a pinned animation");
    }
    const int event_count = global_timeline_num_events;
    for (int i = 0; i < clear_hook_count; ++i)
      clear_hooks[i].fn(clear_hooks[i].ctx);
    HS_CHECK(global_timeline_num_events == event_count,
             "clear hook added or removed timeline events");
    destroy_events();
  }

  /**
   * @brief Registers a callback clear() runs before it destroys the events.
   * @param ctx Opaque pointer handed back to @p fn; also the removal key.
   * @param fn Callback; must not add or remove timeline events.
   * @details For state an owner reclaims from a completion callback, which
   * clear() never runs.
   */
  HS_COLD_MEMBER void add_clear_hook(void *ctx, void (*fn)(void *)) {
    HS_CHECK(clear_hook_count < MAX_CLEAR_HOOKS,
             "Timeline clear-hook table full");
    clear_hooks[clear_hook_count++] = ClearHook{ctx, fn};
  }

  /**
   * @brief Unregisters the hook added under @p ctx; a no-op if absent.
   * @param ctx The registration key passed to add_clear_hook().
   * @details Hook order is not preserved.
   */
  HS_COLD_MEMBER void remove_clear_hook(void *ctx) {
    for (int i = 0; i < clear_hook_count; ++i) {
      if (clear_hooks[i].ctx == ctx) {
        clear_hooks[i] = clear_hooks[--clear_hook_count];
        return;
      }
    }
  }

  /**
   * @brief Adds a new animation event to the timeline.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting; 0 and 1 both
   * start on the next step().
   * @param animation The animation object.
   * @return Reference to the Timeline object.
   */
  template <typename A> Timeline &add(int in_frames, A animation) {
    add_get(in_frames, std::move(animation), Pin::UNPINNED);
    return *this;
  }

  /**
   * @brief Adds an event whose delay, animation, and callbacks freeze together.
   * @tparam A The animation type.
   * @param in_frames Active frames to wait before starting.
   * @param animation The animation object.
   * @param paused Pause flag that must outlive the event.
   * @return Reference to the Timeline object.
   * @note The gate a paused effect wants: a pending start delay is preserved,
   * step() advancing e.start in lockstep. An animation's own `paused` pointer
   * instead only early-returns from step(), so its event's delay keeps
   * elapsing while the flag is set.
   */
  template <typename A>
  Timeline &add_pausable(int in_frames, A animation, const bool *paused) {
    HS_CHECK(paused != nullptr, "pausable timeline event needs a pause flag");
    add_get(in_frames, std::move(animation), Pin::UNPINNED, paused);
    return *this;
  }

  /**
   * @brief Like add(), but returns the typed pointer to the inline-stored
   * animation. Use when you need to hold a reference for later mutation.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting; 0 and 1 both
   * start on the next step().
   * @param animation The animation object.
   * @param pin Pin::PINNED: the caller retains the pointer across frames, so
   * the event must never move: the animation must be infinite or repeating and
   * no finite, non-repeating event may precede it (both trap). Pin::UNPINNED:
   * the pointer is valid only at the call site.
   * @param paused Optional event-level pause gate.
   * @param owner Optional lifetime owner used by cancel_owner().
   * @return Typed pointer to the inline-stored animation, or nullptr if full
   * for an UNPINNED add. A PINNED add traps when full.
   */
  template <typename A>
  A *add_get(int in_frames, A animation, Pin pin, const bool *paused = nullptr,
             const void *owner = nullptr) {
    static_assert(sizeof(A) <= TimelineEvent::MAX_ANIM_SIZE,
                  "Animation type exceeds TimelineEvent inline storage");
    static_assert(alignof(A) <= alignof(std::max_align_t),
                  "Animation type is over-aligned for TimelineEvent inline "
                  "storage (placement-new would be misaligned)");
    HS_CHECK(in_frames >= 0, "Timeline delay must be non-negative");
    const uint32_t delay = static_cast<uint32_t>(in_frames);
    HS_CHECK(delay <= UINT32_MAX - global_timeline_t,
             "Timeline start frame overflow");
    if (global_timeline_num_events >= MAX_EVENTS) {
      HS_CHECK(pin == Pin::UNPINNED,
               "Timeline full, dropped a pinned animation");
      // Only the saturation episode's first drop logs.
      ++global_timeline_dropped;
      if (!global_timeline_drop_logged) {
        global_timeline_drop_logged = true;
        hs::log("Timeline full, failed to add animation!");
      }
      return nullptr;
    }
    if (pin == Pin::PINNED) {
      // The pinned animation itself must never complete: it would be destroyed
      // under the caller's retained pointer.
      HS_CHECK(!animation.is_finite() || animation.repeats(),
               "pinned animation must be infinite or repeating");
      // A finite, non-repeating predecessor is removed on completion and would
      // relocate this pinned event, so reject it up front. A repeating/infinite
      // predecessor can still be removed if cancel()ed later; move_into traps
      // then.
      for (int i = 0; i < global_timeline_num_events; ++i) {
        IAnimation *prev = global_timeline_events[i].animation();
        HS_CHECK(!prev || !prev->is_finite() || prev->repeats(),
                 "pinned animation added after a finite non-repeating one");
      }
    }
    auto &e = global_timeline_events[global_timeline_num_events++];
    HS_CHECK(!e.manager, "add_get would overwrite a live animation");
    e.start = global_timeline_t + delay;
    e.pinned = (pin == Pin::PINNED);
    e.paused = paused;
    e.owner = owner;
    auto *ptr = new (e.storage) A(std::move(animation));
    e.iface = static_cast<IAnimation *>(ptr);
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      // std::launder to recover a usable A* from the placement-new'd storage.
      A *obj = std::launder(reinterpret_cast<A *>(src.storage));
      if (dst) {
        dst->iface =
            static_cast<IAnimation *>(new (dst->storage) A(std::move(*obj)));
      }
      obj->~A();
    };
    return ptr;
  }

  /**
   * @brief Cancels every event bound to the retiring lifetime owner.
   * @param owner Lifetime token shared by the events to cancel.
   */
  HS_COLD_MEMBER void cancel_owner(const void *owner) {
    bool retiring_predecessor = false;
    for (int i = 0; i < global_timeline_num_events; ++i) {
      const auto &event = global_timeline_events[i];
      if (event.owner == owner)
        retiring_predecessor = true;
      else
        HS_CHECK(!retiring_predecessor || event.owner == nullptr ||
                     !event.pinned || event.iface->is_canceled(),
                 "retire later pinned owners before their predecessors");
    }
    for (int i = 0; i < global_timeline_num_events; ++i) {
      auto &event = global_timeline_events[i];
      if (event.owner == owner && event.animation()) {
        static_cast<Animation::AnimationCommon *>(event.animation())->cancel();
        event.paused = nullptr;
        event.owner = nullptr;
      }
    }
  }

  /**
   * @brief Animations rejected so far because the timeline was full.
   * @details Process-wide and wraps modulo 2^32; clear() and a new Timeline do
   * not reset it. Only the first drop of each saturation episode logs. A drop
   * permanently ends any chain that re-arms itself from a .then() callback.
   * @return Number of dropped add()/add_get() calls modulo 2^32.
   */
  static uint32_t dropped_events() { return global_timeline_dropped; }

  /**
   * @brief Event slots still free before add()/add_get() starts dropping.
   * @return MAX_EVENTS minus the current event count.
   * @details step() runs a completing event's post_callback() before it destroys
   * that event and recomputes the count, so a .then() re-arm issued from the
   * callback is appended while the completing event still holds its slot. A
   * chain that re-arms itself must budget against this count, not against the
   * post-compaction one.
   */
  static int remaining() { return MAX_EVENTS - global_timeline_num_events; }

  /**
   * @brief Event slots currently held.
   * @return Number of live events, including one whose post_callback() is
   * running (see remaining()).
   */
  static int event_count() { return global_timeline_num_events; }

  /**
   * @brief Advances the timeline by one frame, stepping all active or starting
   * animations.
   * @param canvas The current canvas buffer.
   */
  HS_COLD_MEMBER void step(Canvas &canvas) {
    struct StepScope {
      bool &flag;
      explicit StepScope(bool &f) : flag(f) { flag = true; }
      ~StepScope() { flag = false; }
    } step_scope(stepping);

    ++global_timeline_t;

    int write_idx = 0;
    int active_cnt =
        global_timeline_num_events; // Snapshot count before callbacks
                                    // potentially add more

    // Collapse each shared Orientation once before any animation steps it.
    const void *collapsed_ids[MAX_COLLAPSE_IDS];
    int collapsed_cnt = 0;
    for (int i = 0; i < active_cnt; ++i) {
      const void *id = event_orientation_id(global_timeline_events[i]);
      if (!id)
        continue;
      bool already_collapsed = false;
      for (int j = 0; j < collapsed_cnt; ++j) {
        if (collapsed_ids[j] == id) {
          already_collapsed = true;
          break;
        }
      }
      if (!already_collapsed && collapsed_cnt == MAX_COLLAPSE_IDS) {
        for (int j = 0; j < i; ++j) {
          if (event_orientation_id(global_timeline_events[j]) == id) {
            already_collapsed = true;
            break;
          }
        }
      }
      if (already_collapsed)
        continue;
      if (collapsed_cnt < MAX_COLLAPSE_IDS)
        collapsed_ids[collapsed_cnt++] = id;
      global_timeline_events[i].animation()->collapse_orientation();
    }

    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];

      if (global_timeline_t < e.start && e.animation()->is_canceled()) {
        e.animation()->post_callback();
        e.destroy();
        continue;
      }

      if (event_paused(e)) {
        const bool started = global_timeline_t >= e.start;
        HS_CHECK(e.start < UINT32_MAX, "paused timeline start frame overflow");
        ++e.start;
        if (started) {
          IAnimation *anim = e.animation();
          HS_CHECK(anim, "paused timeline event holds no animation");
          if (!anim->is_canceled())
            anim->step_paused(canvas);
          // step_paused() never advances the animation, so done() here means
          // cancel() (or finish()); complete the event now.
          if (anim->done() && !anim->repeats()) {
            anim->post_callback();
            HS_CHECK(!e.pinned || anim->is_canceled(),
                     "pinned animation completed while paused; only cancel() "
                     "may destroy a pinned event");
            e.destroy();
            continue;
          }
        }
        if (i != write_idx)
          e.move_into(global_timeline_events[write_idx]);
        ++write_idx;
        continue;
      }

      if (global_timeline_t < e.start) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
        continue;
      }

      IAnimation *anim = e.animation();
      HS_CHECK(anim, "timeline event holds no animation");
      if (!anim->is_canceled())
        anim->step(canvas);

      // Completion & Cleanup
      bool is_done = anim->done();
      bool keep = true;

      if (is_done) {
        bool does_repeat = anim->repeats();
        if (does_repeat) {
          anim->rewind();
          anim->post_callback();
          // A callback that cancels (or finishes) its own repeating animation
          // leaves it done and non-repeating; completing it here keeps the
          // removal branch from firing .then() a second time next frame.
          if (anim->done() && !anim->repeats())
            keep = false;
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
        HS_CHECK(!e.pinned || anim->is_canceled(),
                 "pinned animation completed; only cancel() may destroy a "
                 "pinned event");
        e.destroy();
      }
    }

    // Move events added during callbacks into the gap left by completed ones;
    // a pinned one would trap in move_into.
    HS_CHECK(global_timeline_num_events >= active_cnt,
             "callback shrank the timeline mid-step; "
             "new_vals_count would go negative");
    int new_vals_count = global_timeline_num_events - active_cnt;
    HS_CHECK(write_idx <= active_cnt,
             "timeline compaction wrote past the events it scanned");
    if (new_vals_count > 0 && write_idx < active_cnt) {
      // The source span [active_cnt, ...) and the destination span
      // [write_idx, ...) can overlap, but write_idx + i < active_cnt + i for
      // every i, so this forward loop reads each source slot before a write
      // reaches it.
      for (int i = 0; i < new_vals_count; ++i) {
        global_timeline_events[active_cnt + i].move_into(
            global_timeline_events[write_idx + i]);
      }
    }

    global_timeline_num_events = write_idx + new_vals_count;
    if (global_timeline_num_events == 0)
      global_timeline_drop_logged = false;
  }

  /**
   * @brief Current global frame count (number of step() calls since the
   * Timeline was constructed).
   * @return The shared timeline frame counter, advanced once per step().
   */
  static uint32_t frame() { return global_timeline_t; }

  static constexpr int MAX_EVENTS =
      TIMELINE_MAX_EVENTS; /**< Must match global_timeline_events array size. */

  /** @brief clear_hooks capacity. */
  static constexpr int MAX_CLEAR_HOOKS = 4;

  /**
   * @brief Distinct Orientation ids step()'s collapse pass caches per frame.
   * @details Exceeding it costs a rescan, not a wrong collapse.
   */
  static constexpr int MAX_COLLAPSE_IDS = 16;

private:
  static bool event_paused(const TimelineEvent &event) {
    return event.paused && *event.paused;
  }

  /// Orientation id of an event due to step this frame, nullptr otherwise.
  static const void *event_orientation_id(TimelineEvent &event) {
    if (event_paused(event) || global_timeline_t < event.start)
      return nullptr;
    IAnimation *anim = event.animation();
    return anim ? anim->orientation_id() : nullptr;
  }

  /**
   * @brief One registered clear() observer.
   */
  struct ClearHook {
    void *ctx = nullptr;
    void (*fn)(void *) = nullptr;
  };

  ClearHook clear_hooks[MAX_CLEAR_HOOKS] = {};
  int clear_hook_count = 0;
  bool stepping = false; /**< True for the duration of step(); gates clear(). */

  /**
   * @brief Destroys every stored animation and empties the slot table.
   */
  void destroy_events() {
    for (int i = 0; i < global_timeline_num_events; ++i) {
      global_timeline_events[i].destroy();
    }
    global_timeline_num_events = 0;
    global_timeline_drop_logged = false;
  }

  /**
   * @brief Unguarded teardown for construction/destruction: destroys every event
   * and rewinds the global frame cursor.
   * @details Skips clear()'s pin check: no add_get() handle spans an instance
   * boundary.
   */
  void reset_storage() {
    destroy_events();
    for (auto &event : global_timeline_events) {
      event.start = 0;
      event.pinned = false;
      event.paused = nullptr;
      event.owner = nullptr;
      event.manager = nullptr;
      event.iface = nullptr;
    }
    global_timeline_t = 0;
  }
};
