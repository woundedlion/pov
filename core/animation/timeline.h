/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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
 * @details Stores the animation inline to avoid arena allocation (survives
 * compaction).
 */
struct TimelineEvent {
  // Inline storage budget for a type-erased animation: 112 B for device/WASM.
  // The 64-bit host inflates every embedded pointer, pushing pointer-bearing
  // animations past that budget (ColorWipe is 128 B there, 100 B on wasm32), so
  // HS_TEST_BUILD widens it for that build only, leaving the device footprint
  // unchanged.
#ifdef HS_TEST_BUILD
  static constexpr size_t MAX_ANIM_SIZE = 256;
#else
  static constexpr size_t MAX_ANIM_SIZE = 112;
#endif

  uint32_t start =
      0; /**< Global frame at which the animation becomes eligible to step. */
  /**
   * @brief Whether this event's inline animation pointer was handed out via
   * Timeline::add_get().
   * @details Compaction must never relocate such an event — doing so dangles the
   * caller's cached pointer.
   */
  bool handled = false;
  const bool *paused = nullptr; /**< Optional event-level pause gate. */
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
   * — only a properly-typed upcast performs the base adjustment.
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
    HS_CHECK(!dst.manager,
             "move_into would leak the destination's live animation");
    dst.start = start;
    // Clears a stale flag: destroy() leaves handled set on a canceled pinned
    // event, and this slot may be recycling one.
    dst.handled = handled;
    dst.paused = paused;
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
extern uint32_t global_timeline_t;       // current global frame count
extern int global_timeline_num_events;   // current number of active events
extern uint32_t global_timeline_dropped; // monotonic full-timeline drop count

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
   * @brief Whether an add_get() caller retains the returned pointer across
   * frames (see add_get()).
   */
  enum class Pin { UNPINNED, PINNED };

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
   * @details Traps if called from inside step() — a completion callback runs
   * before its event is destroyed, so destroy_events() would free the callable
   * whose frame is still executing. Traps on a pinned event
   * (add_get(Pin::PINNED)) too: destroying one dangles the caller's retained
   * animation pointer, the same contract move_into() and step()'s destroy branch
   * enforce. The global frame cursor is deliberately NOT rewound — it is
   * process-global and effects derive a phase from frame(), so a mid-run rewind
   * desynchronises them. Only construction/destruction, which no retained handle
   * spans, resets the cursor.
   */
  void clear() {
    HS_CHECK(!stepping, "clear() from inside step() would destroy the "
                        "animation whose callback is running");
    for (int i = 0; i < global_timeline_num_events; ++i) {
      HS_CHECK(!global_timeline_events[i].handled,
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
   * clear() never runs (it destroys events outright) — TransformerPool's slot
   * table is the case. Cold path: registration happens at effect init.
   */
  HS_COLD_MEMBER void add_clear_hook(void *ctx, void (*fn)(void *)) {
    HS_CHECK(clear_hook_count < MAX_CLEAR_HOOKS,
             "Timeline clear-hook table full");
    clear_hooks[clear_hook_count++] = ClearHook{ctx, fn};
  }

  /**
   * @brief Unregisters the hook added under @p ctx; a no-op if absent.
   * @param ctx The registration key passed to add_clear_hook().
   * @details Hook order is not preserved (the last entry backfills the hole):
   * the hooks are independent.
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
    // An add() caller keeps no handle, so the event compacts normally.
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
   * start on the next step(), and every existing schedule is tuned to that.
   * @param animation The animation object.
   * @param pin Pin::PINNED (default): the caller intends to RETAIN this pointer
   * across frames, so the event is marked handled and step()'s compaction traps
   * (move_into) rather than relocating it out from under the cached pointer.
   * Such a retained handle is only safe when the animation is infinite and added
   * before any finite one (so no earlier event is ever removed to shift it) —
   * the contract the direct callers rely on; the trap enforces it. Pass
   * Pin::UNPINNED for a TRANSIENT pointer used only at the call site and not
   * kept across frames (e.g. TransformerPool::spawn, whose finite animations are
   * compacted normally and whose return is typically discarded).
   * @param paused Optional event-level pause gate.
   * @return Typed pointer to the inline-stored animation, or nullptr if full.
   */
  template <typename A>
  A *add_get(int in_frames, A animation, Pin pin = Pin::PINNED,
             const bool *paused = nullptr) {
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
      // A pinned caller retains the return value; dropping hands back a nullptr
      // no call site null-checks.
      HS_CHECK(pin == Pin::UNPINNED, "Timeline full, dropped a pinned animation");
      ++global_timeline_dropped;
      hs::log("Timeline full, failed to add animation! (%lu dropped)",
              (unsigned long)global_timeline_dropped);
      return nullptr;
    }
    if (pin == Pin::PINNED) {
      // The pinned animation itself must never complete: it would be destroyed
      // under the caller's retained pointer.
      HS_CHECK(!animation.is_finite() || animation.repeats(),
               "pinned animation must be infinite or repeating");
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
    e.start = global_timeline_t + delay;
    e.handled = (pin == Pin::PINNED);
    e.paused = paused;
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
   * @brief Animations rejected so far because the timeline was full.
   * @details Monotonic and process-wide: never reset, not even by clear() or a
   * new Timeline. A drop is silent apart from its log line and permanently
   * ends any chain that re-arms itself from a .then() callback, so a nonzero
   * count is the only lasting evidence.
   * @return Total number of dropped add()/add_get() calls.
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
  int remaining() const { return MAX_EVENTS - global_timeline_num_events; }

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

    // Collapse each distinct Orientation exactly once, before any animation steps
    // it: collapsing per-animation would discard the sub-frame motion-blur history
    // a sharing animation already built this frame. Track the already-collapsed
    // ids in a fixed-size stack scratch set.
    const void *collapsed_ids[MAX_EVENTS];
    int collapsed_n = 0;
    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];
      if (event_paused(e))
        continue;
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

      if (event_paused(e)) {
        const bool started = global_timeline_t >= e.start;
        HS_CHECK(e.start < UINT32_MAX,
                 "paused timeline start frame overflow");
        ++e.start;
        if (started) {
          IAnimation *anim = e.animation();
          HS_CHECK(anim);
          anim->step_paused(canvas);
        }
        if (i != write_idx)
          e.move_into(global_timeline_events[write_idx]);
        ++write_idx;
        continue;
      }

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
    HS_CHECK(global_timeline_num_events >= active_cnt,
             "callback shrank the timeline mid-step; "
             "new_vals_count would go negative");
    int new_vals_count = global_timeline_num_events - active_cnt;
    // Each kept event advances write_idx by one, so it never outruns the events
    // scanned.
    HS_CHECK(write_idx <= active_cnt);
    if (new_vals_count > 0 && write_idx < active_cnt) {
      // The source span [active_cnt, ...) and the destination span
      // [write_idx, ...) can overlap, but write_idx + i < active_cnt + i for
      // every i, so this forward loop reads each source slot before a write
      // reaches it. Reversing the loop would overwrite unread sources.
      for (int i = 0; i < new_vals_count; ++i) {
        global_timeline_events[active_cnt + i].move_into(
            global_timeline_events[write_idx + i]);
      }
    }

    global_timeline_num_events = write_idx + new_vals_count;
  }

  /**
   * @brief Current global frame count (number of step() calls since the
   * Timeline was constructed).
   * @return The shared timeline frame counter, advanced once per step().
   * @details Lets an effect derive a phase from the timeline's own clock rather
   *          than a parallel per-effect counter that can silently desync if it
   *          isn't advanced in exact lockstep with step().
   */
  static uint32_t frame() { return global_timeline_t; }

  static constexpr int MAX_EVENTS =
      TIMELINE_MAX_EVENTS; /**< Must match global_timeline_events array size. */

  static constexpr int MAX_CLEAR_HOOKS = 4; /**< clear_hooks capacity. */

private:
  static bool event_paused(const TimelineEvent &event) {
    return event.paused && *event.paused;
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
  }

  /**
   * @brief Unguarded teardown for construction/destruction: destroys every event
   * and rewinds the global frame cursor.
   * @details The instance boundary is the one point no add_get() handle can span
   * (the live-guard permits a single instance), so this skips clear()'s pin
   * check and owns the cursor rewind that starts each effect at frame 0.
   */
  void reset_storage() {
    destroy_events();
    global_timeline_t = 0;
  }
};
