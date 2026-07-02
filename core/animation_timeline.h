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
 * @brief Compile-time segue policies for MeshCarousel: the styles by which one
 * mesh hands the sphere to the next.
 * @details A segue owns the scheduling shape of a mesh-to-mesh transition (its
 * schedule() hook, whose return value is the delay until the next transition
 * begins) and the meaning of the phase value the scheduled animation feeds the
 * draw functor: phase ramps 0→1 over the incoming window, holds 1, and falls
 * back to 0 over the outgoing window. The shading hooks translate that phase
 * into pixels:
 *
 *   bool   visible(phase)   — whether drawing at this phase is worthwhile
 *   float  opacity(phase)   — global alpha multiplier
 *   float  fill(&t, phase)  — coverage mask; may remap the edge-distance t
 *   Color4 grade(c, phase)  — color regrade after the palette lookup
 *
 * Optional hooks, detected per policy with `requires` so unused paths compile
 * out of the render loop:
 *
 *   void   retarget(v)               — re-randomize per-transition state
 *   Vector warp(v, phase)            — pre-ripple unit-sphere vertex warp
 *   float  face_offset(center, i)    — per-face sweep ordering in [0, 1]
 *   float  face_phase(phase, offset) — face-local phase from the sweep front
 *
 * Policies are resolved at compile time (no virtuals); Base's identity hooks
 * inline to nothing.
 */
namespace Segue {

/**
 * @brief Schedules one sequential fade-in/fade-out sprite: consecutive sprites
 * never overlap, so a transition renders a single mesh per frame.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to duration/2
 * so the in/out windows never collide.
 * @return duration — the next transition starts as this sprite ends.
 */
inline int schedule_sequential(Timeline &timeline, SpriteFn draw_fn,
                               int duration, int window) {
  int fade = std::min(window, duration / 2);
  timeline.add(0, Animation::Sprite(std::move(draw_fn), duration, fade,
                                    ease_linear, fade, ease_linear));
  return duration;
}

/**
 * @brief Soft sweep front shared by the per-face segues.
 * @param phase Global segue phase in [0, 1].
 * @param offset Face's sweep ordering in [0, 1]; larger extinguishes earlier.
 * @param band Softness of the front, in phase units.
 * @return The face-local phase in [0, 1]: 1 everywhere at phase 1, 0
 * everywhere at phase 0, with faces crossing the front in offset order.
 */
inline float sweep_phase(float phase, float offset, float band) {
  return hs::clamp((phase - offset * (1.0f - band)) / band, 0.0f, 1.0f);
}

/**
 * @brief Identity hooks every segue inherits; a policy shadows only the hooks
 * its transition uses.
 */
struct Base {
  /** @brief Default scheduling: one sequential sprite (see schedule_sequential). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    return schedule_sequential(timeline, std::move(draw_fn), duration, window);
  }
  /** @brief Whether drawing at this phase can produce visible output. */
  bool visible(float phase) const { return phase > 0.005f; }
  /** @brief Global alpha at this phase. */
  float opacity(float) const { return 1.0f; }
  /**
   * @brief Coverage mask over the face interior.
   * @param t Fragment edge distance in [0, 1] (0 at the edge, ~1 at the face
   * center); may be remapped in place for the palette lookup.
   * @return Coverage alpha in [0, 1]; 0 culls the fragment.
   */
  float fill(float &t, float phase) const {
    (void)t;
    (void)phase;
    return 1.0f;
  }
  /** @brief Color regrade applied after the palette lookup. */
  Color4 grade(Color4 c, float) const { return c; }
};

/**
 * @brief Opacity cross-fade between consecutive meshes.
 * @details Phase is opacity. Each transition is one fade-in/fade-out Sprite;
 * the returned delay starts the next transition one fade window before this
 * sprite ends, so the outgoing and incoming sprites overlap and both meshes
 * render during the fade — the cost of this segue is two rasterized meshes
 * per overlap frame. Every other segue is sequential (single mesh per frame).
 */
struct Crossfade : Base {
  /**
   * @brief Schedules the incoming mesh's fading sprite.
   * @param timeline Timeline receiving the sprite.
   * @param draw_fn Draws the incoming mesh at the given opacity.
   * @param duration Total frames the mesh is on screen.
   * @param window Requested fade length in frames; clamped to duration/2 so
   * the fade windows never overlap and sprites cannot pile up beyond two.
   * @return Frames after which the next transition should begin.
   */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    int fade = std::min(window, duration / 2);
    timeline.add(0, Animation::Sprite(std::move(draw_fn), duration, fade,
                                      ease_linear, fade, ease_linear));
    return duration - fade;
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief Faces contract to glowing points at their centers, then the next
 * tessellation blooms back out of the point field.
 * @details Only fragments deeper than the phase-driven inset survive, so the
 * pattern dissolves into a constellation of face-center glints at the swap.
 * The surviving core's edge distance is renormalized so it keeps the full
 * palette gradient as it shrinks.
 */
struct IrisBloom : Base {
  static constexpr float kSoft = 0.08f; /**< Soft rim width, in edge-distance units. */
  float fill(float &t, float phase) const {
    float inset = 1.0f - phase;
    if (t < inset - kSoft)
      return 0.0f;
    float cover = hs::clamp((t - (inset - kSoft)) / kSoft, 0.0f, 1.0f);
    t = hs::clamp((t - inset) / std::max(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief The fill drains until only a glowing band along the edges survives,
 * the meshes swap as lace, then the new fill floods back in.
 * @details The inverse mask of IrisBloom: fragments within the phase-driven
 * band of an edge survive. A line network changing shape reads far less
 * jarring than filled regions changing, which hides the swap.
 */
struct Lace : Base {
  static constexpr float kSoft = 0.08f; /**< Soft band-edge width, in edge-distance units. */
  float fill(float &t, float phase) const {
    if (t > phase + kSoft)
      return 0.0f;
    float cover = hs::clamp((phase + kSoft - t) / kSoft, 0.0f, 1.0f);
    t = hs::clamp(t / std::max(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief A world-fixed day/night line sweeps the sphere, extinguishing the old
 * pattern face by face; the return sweep ignites the new one.
 * @details Face offsets are recomputed from world-space centers each frame, so
 * the terminator stays fixed in the room while the mesh rotates through it.
 */
struct TerminatorSweep : Base {
  static constexpr float kBand = 0.25f; /**< Sweep-front softness, in phase units. */
  Vector axis = Y_AXIS; /**< World-space sweep axis. */
  void retarget(const Vector &v) { axis = v; }
  float face_offset(const Vector &center, int) const {
    return 0.5f * (1.0f + dot(center, axis));
  }
  float face_phase(float phase, float offset) const {
    return sweep_phase(phase, offset, kBand);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief An expanding shockwave erases the pattern outward from a point; its
 * echo redraws the new one.
 * @details Faces nearest the origin extinguish first, so the wave visibly
 * expands. Pairs naturally with the effect's ripple bursts sharing the origin.
 */
struct Shockwave : Base {
  static constexpr float kBand = 0.3f; /**< Wave-front softness, in phase units. */
  Vector origin = Y_AXIS; /**< World-space wave origin. */
  void retarget(const Vector &v) { origin = v; }
  float face_offset(const Vector &center, int) const {
    float angle = fast_acos(hs::clamp(dot(center, origin), -1.0f, 1.0f));
    return 1.0f - angle * (1.0f / PI_F);
  }
  float face_phase(float phase, float offset) const {
    return sweep_phase(phase, offset, kBand);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief Faces wink out individually in a hashed order and the new ones wink
 * back in — a glitter dissolve.
 * @details The order comes from a face-index hash, not geometry, so it is
 * stable across frames (a world-space source would shimmer as the mesh
 * rotates).
 */
struct Sparkle : Base {
  static constexpr float kBand = 0.2f; /**< Per-face fade width, in phase units. */
  float face_offset(const Vector &, int face) const {
    uint32_t h = static_cast<uint32_t>(face) * 2654435761u;
    h ^= h >> 16;
    return static_cast<float>(h & 0xFFFF) * (1.0f / 65535.0f);
  }
  float face_phase(float phase, float offset) const {
    return sweep_phase(phase, offset, kBand);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The polyhedron drains into a vanishing point like water down a sink,
 * then the next one unfurls from it.
 * @details The pull is capped and the mesh fades out well before full
 * collapse: converging faces stack along the view ray exactly like the ripple
 * self-fold, so the deepest pile-up regime is never rendered.
 */
struct Drain : Base {
  static constexpr float kMaxPull = 0.85f; /**< Collapse ceiling; never reaches the focus. */
  Vector focus = Y_AXIS; /**< Vanishing point on the unit sphere. */
  void retarget(const Vector &v) { focus = v; }
  Vector warp(const Vector &v, float phase) const {
    float pull = 1.0f - phase;
    return slerp(v, focus, pull * pull * kMaxPull);
  }
  float opacity(float phase) const {
    return hs::clamp(phase * 3.0f, 0.0f, 1.0f);
  }
  bool visible(float phase) const { return phase > 0.02f; }
};

/**
 * @brief The pattern winds into a spiral around an axis, shears apart, swaps,
 * and unwinds as the new tessellation.
 * @details Twist angle scales with latitude (dot with the axis), so the mesh
 * shears rather than rotating rigidly. Turns are capped and the mesh fades
 * before peak wind to bound the shear-fold overdraw.
 */
struct Vortex : Base {
  static constexpr float kMaxTurns = 1.5f; /**< Full-wind twist, in revolutions. */
  Vector axis = Y_AXIS; /**< Twist axis. */
  void retarget(const Vector &v) { axis = v; }
  Vector warp(const Vector &v, float phase) const {
    float wind = 1.0f - phase;
    float theta = wind * wind * kMaxTurns * 2.0f * PI_F * dot(v, axis);
    return rotate(v, make_rotation(axis, theta));
  }
  float opacity(float phase) const {
    return hs::clamp(phase * 2.0f, 0.0f, 1.0f);
  }
};

/**
 * @brief The whole mesh spins up around an axis until the POV display smears
 * it into bands, swaps at peak speed, and spins back down — a coin flip.
 * @details The warp is rigid, so there is no fold/overdraw hazard and the mesh
 * never fades; the swap hides in the motion blur.
 */
struct SpinFlip : Base {
  static constexpr float kRevs = 3.0f; /**< Extra revolutions at peak spin. */
  Vector axis = Y_AXIS; /**< Spin axis. */
  void retarget(const Vector &v) { axis = v; }
  Vector warp(const Vector &v, float phase) const {
    float wind = 1.0f - phase;
    return rotate(v, make_rotation(axis, wind * wind * kRevs * 2.0f * PI_F));
  }
};

/**
 * @brief Both palettes converge to molten gold around the swap, then the new
 * mesh blooms back into color.
 * @details Purely palette-domain: geometry never moves. A mild opacity dip at
 * the swap softens the topology pop while both meshes are monochrome.
 */
struct GoldConvergence : Base {
  Pixel gold = Color4(uint8_t{255}, uint8_t{196}, uint8_t{64}).color; /**< Linear-space convergence color. */
  Color4 grade(Color4 c, float phase) const {
    return c.lerp(Color4(gold, c.alpha), 1.0f - phase);
  }
  float opacity(float phase) const { return 0.4f + 0.6f * phase; }
};

} // namespace Segue

/**
 * @brief A double-buffered pair of persistent MeshState slots, the
 *        arena-compaction primitives effects need to swap between them, and a
 *        pluggable compile-time segue.
 * @tparam SegueT Segue policy (see namespace Segue) behind schedule_segue().
 * Clients that run their own transition animations (e.g. a `MeshMorph`) keep
 * the default and simply never call it.
 * @details Holds two MeshState slots in `persistent_arena` and a front/back
 * index. Effects still own generation and drawing — they generate into a
 * slot, flip the front index, and reclaim the old slot's space — while the
 * segue owns the transition's animation scheduling.
 *
 * Usage:
 *   MeshCarousel<Segue::Crossfade> carousel;  // in effect members
 *
 *   // Build the initial shape directly into the front slot:
 *   carousel.current().clear();
 *   MeshOps::compile(mesh, carousel.current(), persistent_arena);
 *
 *   // To transition: generate into the back slot, flip, then let the segue
 *   // schedule the animation via schedule_segue (see
 *   // IslamicStars::spawn_shape for the pattern).
 */
template <typename SegueT = Segue::Crossfade> class MeshCarousel {
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
   * @brief Schedules the segue's transition animation for the (already
   * front-flipped) incoming mesh.
   * @param timeline Timeline receiving the segue's animation.
   * @param draw_fn Draws the incoming mesh; the float argument is the segue's
   * phase (opacity for Segue::Crossfade).
   * @param duration Total frames the mesh is on screen.
   * @param window Transition window length in frames, segue-interpreted.
   * @return Frames after which the effect should begin the next transition.
   */
  int schedule_segue(Timeline &timeline, SpriteFn draw_fn, int duration,
                     int window) {
    return segue_.schedule(timeline, std::move(draw_fn), duration, window);
  }

  /**
   * @brief The carousel's segue policy instance (holds per-transition state
   * such as a sweep axis or wave origin).
   */
  SegueT &segue() { return segue_; }
  /** @brief Const view of the segue policy instance. */
  const SegueT &segue() const { return segue_; }

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
  SegueT segue_;       /**< Segue policy instance; per-transition state lives here. */
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

  // An interior length-1 frame holds only the boundary it shares with the prior
  // frame, so it emits nothing; counting it in the span would leave a hole in the
  // age ramp. Normalize against, and index by, the frames that actually contribute.
  size_t num_active = 0;
  for (size_t i = 0; i <= last; ++i)
    if (i == 0 || trail.get(i).length() > 1)
      ++num_active;
  float span = static_cast<float>(num_active);

  size_t active_idx = 0;
  for (size_t i = 0; i <= last; ++i) {
    const auto &frame = trail.get(i);
    size_t frame_size = frame.length();
    if (i != 0 && frame_size <= 1)
      continue;
    size_t start_j = (i == 0) ? 0 : 1;

    for (size_t j = start_j; j < frame_size; ++j) {
      const auto &q = frame.get(j);
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      float global_t = (static_cast<float>(active_idx) + sub_t) / span;
      callback(q, global_t);
    }
    ++active_idx;
  }
}
