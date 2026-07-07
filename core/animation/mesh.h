/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

namespace Animation {

/**
 * @brief Animates a vertex-interpolated crossfade between two meshes.
 * @details Owns its transient state (cloned meshes + SLERP buffers) via a
 * pointer to arena-allocated storage — keeps inline size small for
 * TimelineEvent. The caller provides source/dest MeshState references and an
 * Arena; MeshMorph clones both, builds nearest-vertex correspondence, and
 * interpolates each frame. Destruction runs only member destructors — there is
 * no destructor that reclaims the arena, so the transient bytes persist in the
 * arena until the caller compacts or resets it.
 *
 * Crossfade contract (intentional): only the incoming mesh (mesh_B, carrying
 * dest topology) morphs — each frame its vertices SLERP from their nearest
 * source vertex toward their dest position. The outgoing mesh (mesh_A, a clone
 * of source) holds its geometry and fades out via opacity (op_A = 1 - alpha)
 * while the incoming fades in; the opacities sum to 1 for constant total
 * brightness across the topology swap. The source is cloned (not borrowed) so
 * the animation is self-contained against the caller recycling/mutating source
 * mid-morph, consistent with the borrow contract enforced on the draw callbacks
 * below. The separate draw_outgoing/draw_incoming callbacks exist precisely to
 * shade the two halves independently (see HankinSolids).
 */
class MeshMorph : public AnimationBase<MeshMorph> {
public:
  /**
   * @brief Non-owning per-half draw callback: `void(Canvas&, const MeshState&,
   * float opacity)`. A StoredFunctionRef (not a plain FunctionRef) because it is
   * held as a member and invoked across many frames: the type rejects rvalue
   * temporaries so a dangling inline lambda is a compile error, not a silent
   * use-after-free (see the borrow contract below).
   */
  using MorphDrawFn = StoredFunctionRef<void(Canvas &, const MeshState &, float)>;

  /**
   * @brief Constructs a MeshMorph with separate shading for the two halves.
   * @param source The outgoing mesh (cloned, not borrowed).
   * @param dest The incoming mesh whose topology the morph targets.
   * @param arena Arena providing backing storage for cloned meshes and buffers.
   * @param draw_outgoing Draw callback for the fading-out source clone.
   * @param draw_incoming Draw callback for the fading-in morphing mesh.
   * @param duration The crossfade duration in frames.
   * @param easing_fn The easing function applied to crossfade progress.
   * @note The cloned meshes and position buffers are arena-allocated with no
   *   per-instance reclamation; the caller must compact the arena between
   *   successive morphs or it grows unbounded (see HankinSolids/MeshFeedback).
   */
  MeshMorph(const MeshState &source, const MeshState &dest, Arena &arena,
            MorphDrawFn draw_outgoing, MorphDrawFn draw_incoming, int duration,
            EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(duration, false), easing_fn(easing_fn),
        draw_outgoing(draw_outgoing), draw_incoming(draw_incoming) {
    HS_CHECK(!source.vertices.is_empty());
    HS_CHECK(!dest.vertices.is_empty());
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();

    MeshOps::clone(source, buf_->mesh_A, arena);
    MeshOps::clone(dest, buf_->mesh_B, arena);

    // step() indexes mesh_B.vertices by dest's vertex count, so trap here if
    // clone() ever welds/dedups rather than writing OOB per-frame.
    HS_CHECK(buf_->mesh_B.vertices.size() == dest.vertices.size());

    buf_->start_pos.bind(arena, dest.vertices.size());
    buf_->end_pos.bind(arena, dest.vertices.size());

    // Symmetry-breaking twist to avoid degenerate nearest-vertex mapping
    Vector twist_axis = Vector(0.0f, 0.0f, 1.0f);
    bool has_poles = false;
    for (const auto &v : source.vertices) {
      if (std::abs(v.z) > 0.99f && std::abs(v.x) < 0.01f && std::abs(v.y) < 0.01f)
        has_poles = true;
    }
    if (has_poles) {
      twist_axis = Vector(1.0f, 1.0f, 1.0f).normalized();
    }
    Quaternion twist = make_rotation(twist_axis, 0.05f);

    // Nearest-vertex matching (greatest dot) and per-frame slerp both require
    // unit-length inputs; all mesh sources sit on the unit sphere.
    for (const auto &v : source.vertices)
      HS_CHECK(std::abs(dot(v, v) - 1.0f) < 1e-3f, "MeshMorph source vertex not unit-length");
    for (const auto &v : dest.vertices)
      HS_CHECK(std::abs(dot(v, v) - 1.0f) < 1e-3f, "MeshMorph dest vertex not unit-length");

    // Build nearest-vertex correspondence: an O(V_dest * V_source) brute force,
    // run once at construction. Matched by greatest dot product against the
    // twist-biased dest vertex (the twist breaks ties on symmetric meshes).
    for (size_t i = 0; i < dest.vertices.size(); ++i) {
      Vector v_biased = rotate(dest.vertices[i], twist);
      int best_idx = 0;
      float max_dot = -9999.0f;
      for (size_t j = 0; j < source.vertices.size(); ++j) {
        float d = dot(v_biased, source.vertices[j]);
        if (d > max_dot) {
          max_dot = d;
          best_idx = static_cast<int>(j);
        }
      }
      buf_->start_pos.push_back(source.vertices[best_idx]);
      buf_->end_pos.push_back(dest.vertices[i]);
    }
  }

  // Borrow contract: the draw callbacks are non-owning StoredFunctionRefs read
  // every frame, so they must outlive the timeline; StoredFunctionRef rejects a
  // temporary at the MorphDrawFn parameter, so no `= delete` overload is needed.

  /**
   * @brief Steps the crossfade: interpolates vertices and renders both halves.
   * @param canvas The canvas buffer passed to the draw callbacks.
   */
  void step(Canvas &canvas) override {
    // Increment-first so the final frame lands exactly on the destination mesh
    // (alpha == 1); the skipped progress==0 frame is immaterial.
    AnimationBase::step(canvas);

    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    float alpha = easing_fn(progress);

    for (size_t i = 0; i < buf_->end_pos.size(); ++i) {
      buf_->mesh_B.vertices[i] =
          slerp(buf_->start_pos[i], buf_->end_pos[i], alpha);
    }

    float op_A = 1.0f - alpha;
    if (op_A > 0.01f)
      draw_outgoing(canvas, buf_->mesh_A, op_A);
    if (alpha > 0.01f)
      draw_incoming(canvas, buf_->mesh_B, alpha);
  }

private:
  /**
   * @brief Arena-allocated transient data — keeps MeshMorph inline size small.
   */
  struct Transients {
    MeshState mesh_A;            /**< Outgoing mesh clone. */
    MeshState mesh_B;           /**< Incoming morphing mesh clone. */
    ArenaVector<Vector> start_pos; /**< Per-vertex nearest-source start points. */
    ArenaVector<Vector> end_pos;   /**< Per-vertex dest end points. */
  };

  Transients *buf_;          /**< Pointer to arena-allocated transient state. */
  EasingFn easing_fn;        /**< Easing curve applied to crossfade progress. */
  MorphDrawFn draw_outgoing; /**< Draw callback for the outgoing half. */
  MorphDrawFn draw_incoming; /**< Draw callback for the incoming half. */
};

} // namespace Animation

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
 *   float  face_offset(center, i, cls) — per-face sweep ordering in [0, 1]
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
 * @brief Soft sweep front used by Shockwave.
 * @param phase Global segue phase in [0, 1].
 * @param offset Face's sweep ordering in [0, 1]; larger extinguishes earlier.
 * @param band Softness of the front, in phase units.
 * @return The face-local phase in [0, 1]: 1 everywhere at phase 1, 0
 * everywhere at phase 0, with faces crossing the front in offset order.
 * @details The sqrt ease keeps the hand-off out of black: sweeps schedule
 * sequentially (single mesh per frame), so both meshes sit at low phase
 * around the swap and a linear front would leave the sphere mostly dark for
 * a large slice of each fade window. Accelerating the front through the
 * low-phase end compresses that all-dark valley to a blink, while the
 * endpoints stay exact (phase 1 remains the identity plateau).
 */
inline float sweep_phase(float phase, float offset, float band) {
  float p = std::sqrt(phase);
  return hs::clamp((p - offset * (1.0f - band)) / band, 0.0f, 1.0f);
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
 * @brief A world-fixed day/night line sweeps the sphere at constant speed; the
 * moment it reaches a face, that face fades over kFadeFrames frames. The
 * extinguishing sweep takes the old pattern, the return sweep ignites the new.
 * @details Face offsets are recomputed from world-space centers each frame, so
 * the terminator stays fixed in the room while the mesh rotates through it.
 * The front crosses the sphere over the fade window minus one per-face fade,
 * so the last-touched face still completes: face phases are exactly 1 at
 * phase 1 and 0 at phase 0.
 */
struct TerminatorSweep : Base {
  static constexpr int kFadeFrames = 16; /**< Per-face fade length, in frames, from the front's touch (1 s at 16 fps). */
  Vector axis = Y_AXIS;   /**< World-space sweep axis. */
  float fade_frac = 0.5f; /**< kFadeFrames over the scheduled fade window; set by schedule(). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    int fade = std::min(window, duration / 2);
    fade_frac = std::min(1.0f, static_cast<float>(kFadeFrames) /
                                   static_cast<float>(std::max(fade, 1)));
    return schedule_sequential(timeline, std::move(draw_fn), duration, window);
  }
  void retarget(const Vector &v) { axis = v; }
  float face_offset(const Vector &center, int, int) const {
    return 0.5f * (1.0f + dot(center, axis));
  }
  float face_phase(float phase, float offset) const {
    return hs::clamp((phase - offset * (1.0f - fade_frac)) / fade_frac, 0.0f,
                     1.0f);
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
  float face_offset(const Vector &center, int, int) const {
    float angle = fast_acos(hs::clamp(dot(center, origin), -1.0f, 1.0f));
    return 1.0f - angle * (1.0f / PI_F);
  }
  float face_phase(float phase, float offset) const {
    return sweep_phase(phase, offset, kBand);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The pattern breaks down one topology class at a time: all faces of a
 * class fade together, each class fully gone before the next starts, in a
 * random class order reshuffled per transition; the new tessellation
 * reassembles class by class the same way.
 * @details Faces group by their palette-slot class (the same topology→slot
 * mapping the fragment shader uses), so each color family vanishes as a unit.
 * The class windows are exactly abutting equal slices of the phase range,
 * linear in time — deliberately not sweep_phase's eased front, so every class
 * gets an equal share of the window. The kBlackDwell slice nearest the swap
 * is held fully black: without it the last class's fade runs to the sprite's
 * final frame and the incoming mesh appears one frame later, so the class
 * visibly pops instead of completing. The client supplies the class count
 * and triggers the reshuffle through reorder(), once per transition.
 */
struct Breakdown : Base {
  static constexpr int kMaxClasses = 16; /**< rank[] capacity. */
  static constexpr float kBlackDwell = 0.1f; /**< Phase slice held all-black at the swap end. */
  int num_classes = 1;            /**< Live class count, set by reorder(). */
  uint8_t rank[kMaxClasses] = {}; /**< rank[class]: fade position; 0 vanishes first. */
  /** @brief Re-randomizes the class fade order for the next transition. */
  void reorder(int classes) {
    num_classes = hs::clamp(classes, 1, kMaxClasses);
    for (int i = 0; i < num_classes; ++i)
      rank[i] = static_cast<uint8_t>(i);
    std::shuffle(rank, rank + num_classes, hs::random());
  }
  float face_offset(const Vector &, int, int cls) const {
    if (num_classes <= 1)
      return 0.0f;
    int r = rank[(cls >= 0 && cls < num_classes) ? cls : 0];
    return static_cast<float>(num_classes - 1 - r) /
           static_cast<float>(num_classes - 1);
  }
  float face_phase(float phase, float offset) const {
    // Class windows tile [kBlackDwell, 1]; phase 1 stays the identity plateau.
    float band = (1.0f - kBlackDwell) / static_cast<float>(num_classes);
    return hs::clamp(
        (phase - kBlackDwell - offset * (1.0f - kBlackDwell - band)) / band,
        0.0f, 1.0f);
  }
  float opacity(float phase) const { return phase; }
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
