/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file segue.h
 * @brief Animation fragment: the Segue transition library.
 */

#include "color/color.h"

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
 *   float  face_phase(phase, offset, fade_frac) — face-local phase from the
 *                                                 front
 *   void   reorder(face_classes)     — re-derive the per-transition class
 *                                      ordering (Segue::NeedsClasses)
 *   MaskPair mask_pair(phase, frame) — complementary edge masks the effect
 *                                      hands its two draws (Segue::Masked)
 *
 * A policy defining face_offset must define face_phase; face_fade_frac is
 * Base's (1 = fade over the whole window) unless the policy shadows it.
 *
 * reorder and mask_pair are contracts on the effect, not the draw path: a
 * NeedsClasses policy left un-reordered fades every class as one, and a Masked
 * policy drawn without its masks rasterizes both meshes.
 *
 * The per-face hooks and the fragment hooks are mutually exclusive: a per-face
 * draw path resolves phase and opacity once per face and shades through a
 * palette pointer, so it never calls fill or grade. MeshCarousel
 * static_asserts the combination away rather than dropping the hooks silently.
 *
 * A per-face policy may also declare `static constexpr bool LOCAL_SWEEP =
 * true` to order faces by the untransformed mesh instead of world-space
 * centers: the front then rides the mesh's rotation rather than staying
 * fixed in the room.
 *
 * Policies are resolved at compile time (no virtuals); Base's identity hooks
 * inline to nothing.
 */
namespace Segue {

/**
 * @brief Schedules one sprite with symmetric linear fade-in/fade-out.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to duration/2
 * so the in/out windows never collide.
 * @param paused Optional event-level pause gate.
 * @return The clamped fade length, from which each policy derives its own
 * return offset.
 */
inline int schedule_faded_sprite(Timeline &timeline, SpriteFn draw_fn,
                                 int duration, int window,
                                 const bool *paused = nullptr) {
  int fade = std::min(window, duration / 2);
  Animation::Sprite sprite(std::move(draw_fn), duration, fade, ease_linear,
                           fade, ease_linear);
  if (paused)
    timeline.add_pausable(0, std::move(sprite), paused);
  else
    timeline.add(0, std::move(sprite));
  return fade;
}

/**
 * @brief Schedules one sequential fade-in/fade-out sprite: consecutive sprites
 * never overlap, so a transition renders a single mesh per frame.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to duration/2
 * so the in/out windows never collide.
 * @param paused Optional event-level pause gate.
 * @return duration — the next transition starts as this sprite ends.
 */
inline int schedule_sequential(Timeline &timeline, SpriteFn draw_fn,
                               int duration, int window,
                               const bool *paused = nullptr) {
  schedule_faded_sprite(timeline, std::move(draw_fn), duration, window, paused);
  return duration;
}

/**
 * @brief Schedules one fade-in/fade-out sprite whose successor starts before it
 * ends, so consecutive sprites coexist and both meshes render during the
 * overlap.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to duration/2
 * so the in/out windows never collide.
 * @param overlap Frames consecutive sprites coexist, clamped to the fade
 * window; negative selects the full window. At 0 the schedule is sequential.
 * @param paused Optional event-level pause gate.
 * @return duration minus the effective overlap.
 */
inline int schedule_overlapped(Timeline &timeline, SpriteFn draw_fn,
                               int duration, int window, int overlap,
                               const bool *paused = nullptr) {
  int fade = schedule_faded_sprite(timeline, std::move(draw_fn), duration,
                                   window, paused);
  return duration - (overlap < 0 ? fade : std::min(overlap, fade));
}

/**
 * @brief Soft sweep front used by Shockwave.
 * @param phase Global segue phase in [0, 1].
 * @param offset Face's sweep ordering in [0, 1]; larger extinguishes earlier.
 * @param band Softness of the front, in phase units.
 * @return The face-local phase in [0, 1]: 1 everywhere at phase 1, 0
 * everywhere at phase 0, with faces crossing the front in offset order.
 * @details The sqrt ease keeps the hand-off out of black: both meshes sit at
 * low phase around the swap, so a linear front would leave the sphere mostly
 * dark; accelerating the front through the low-phase end compresses that to a
 * blink. Endpoints stay exact (phase 1 remains the identity plateau).
 */
inline float sweep_phase(float phase, float offset, float band) {
  float p = std::sqrt(phase);
  return hs::clamp((p - offset * (1.0f - band)) / band, 0.0f, 1.0f);
}

/**
 * @brief Cull gate for a policy whose shading vanishes as phase falls to 0.
 * @param phase Global segue phase in [0, 1].
 * @return Whether drawing at this phase can produce visible output.
 */
inline bool fades_to_black(float phase) { return phase > 0.005f; }

/**
 * @brief Identity hooks every segue inherits; a policy shadows only the hooks
 * its transition uses.
 */
struct Base {
  /** @brief Default scheduling: one sequential sprite (see
   * schedule_sequential). @p paused is the optional event-level pause gate
   * every policy's schedule() takes and forwards. */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    return schedule_sequential(timeline, std::move(draw_fn), duration, window,
                               paused);
  }
  /** @brief Whether drawing at this phase can produce visible output. The
   * identity policy never culls; only a policy whose shading actually vanishes
   * with phase shadows this with fades_to_black(). */
  bool visible(float) const { return true; }
  /** @brief Global alpha at this phase. */
  float opacity(float) const { return 1.0f; }
  /**
   * @brief Coverage mask over the face interior.
   * @param t Fragment edge distance in [0, 1] (0 at the edge, ~1 at the face
   * center); may be remapped in place for the palette lookup.
   * @param phase Transition phase; unused by the identity policy.
   * @return Coverage alpha in [0, 1]; 0 culls the fragment.
   */
  float fill(float &t, float phase) const {
    (void)t;
    (void)phase;
    return 1.0f;
  }
  /** @brief Color regrade applied after the palette lookup. */
  Color4 grade(Color4 c, float) const { return c; }
  /**
   * @brief Per-face fade length as a fraction of the transition window.
   * @return 1 — the face fades over the whole window, so face_phase reduces to
   * the global phase for a policy with no per-face fade of its own.
   */
  float face_fade_frac(int) const { return 1.0f; }
};

/** @brief Whether a policy's schedule() hook takes the full argument set,
 * including the trailing pause gate the carousel forwards. A policy shadowing
 * schedule() with a shorter signature hides Base's and is rejected here rather
 * than silently losing the gate. */
template <typename S>
concept Schedulable =
    requires(S &s, Timeline &timeline, SpriteFn draw_fn, const bool *paused) {
      {
        s.schedule(timeline, std::move(draw_fn), 0, 0, paused)
      } -> std::same_as<int>;
    };

/** @brief Whether a policy defines the per-face hook set. */
template <typename S>
concept PerFace = requires(const S &s, const Vector &c) {
  s.face_offset(c, 0, 0);
  s.face_phase(0.5f, 0.5f, 0.1f);
};

/** @brief Whether a policy orders faces by topology class, so the effect must
 * hand it the per-face classes before each transition. */
template <typename S>
concept NeedsClasses = requires(S &s, const ArenaVector<uint16_t> &classes) {
  s.reorder(classes);
};

/** @brief Whether a policy splits one frame's rasterizer work between the two
 * meshes with complementary pixel masks, which the effect passes to the two
 * draws itself. */
template <typename S>
concept Masked = requires(const S &s) { s.mask_pair(0.5f, 0u); };

/**
 * @brief Whether a policy shadows Base's fragment hooks (fill/grade).
 * @details Pointer-to-member type identity: an inherited hook keeps Base as
 * its class type, a shadowing one takes the policy's.
 */
template <typename S>
inline constexpr bool SHADOWS_FRAGMENT_HOOKS =
    !std::is_same_v<decltype(&S::fill), decltype(&Base::fill)> ||
    !std::is_same_v<decltype(&S::grade), decltype(&Base::grade)>;

/**
 * @brief Opacity cross-fade between consecutive meshes.
 * @details Phase is opacity. Each transition is one fade-in/fade-out Sprite;
 * the returned delay starts the next transition `overlap` frames before this
 * sprite ends, so the outgoing and incoming sprites coexist and both meshes
 * render during those frames — the cost of this segue is two rasterized
 * meshes per overlap frame. At overlap 0 the schedule is sequential (a fade
 * through black) and a single mesh renders per frame.
 */
struct Crossfade : Base {
  int overlap = -1; /**< Frames consecutive sprites coexist; clamped to the
                         fade window, negative selects the full window. */
  /**
   * @brief Schedules the incoming mesh's fading sprite.
   * @param timeline Timeline receiving the sprite.
   * @param draw_fn Draws the incoming mesh at the given opacity.
   * @param duration Total frames the mesh is on screen.
   * @param window Requested fade length in frames; clamped to duration/2 so
   * the fade windows never overlap and sprites cannot pile up beyond two.
   * @param paused Optional event-level pause gate.
   * @return Frames after which the next transition should begin: duration
   * minus the effective overlap.
   */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    return schedule_overlapped(timeline, std::move(draw_fn), duration, window,
                               overlap, paused);
  }
  bool visible(float phase) const { return fades_to_black(phase); }
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
  static constexpr float SOFT =
      0.08f; /**< Soft rim width, in edge-distance units. */
  float fill(float &t, float phase) const {
    float inset = 1.0f - phase;
    if (t < inset - SOFT)
      return 0.0f;
    float cover = hs::clamp((t - (inset - SOFT)) / SOFT, 0.0f, 1.0f);
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
  static constexpr float SOFT =
      0.08f; /**< Soft band-edge width, in edge-distance units. */
  float fill(float &t, float phase) const {
    if (t > phase + SOFT)
      return 0.0f;
    float cover = hs::clamp((phase + SOFT - t) / SOFT, 0.0f, 1.0f);
    t = hs::clamp(t / std::max(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief A day/night line pinned to the mesh sweeps across it; when it reaches
 * a face, that face fades over a per-face random length in [fade_frames_min,
 * fade_frames_max] frames.
 * @details LOCAL_SWEEP anchors the line to the untransformed mesh. Each face's
 * fade length is a stable per-transition hash of its index, so the front frays
 * into an irregular edge. The front crosses over the fade window minus one face
 * fade, so face phases are exactly 1 at phase 1 and 0 at phase 0 for every fade
 * length.
 */
struct TerminatorSweep : Base {
  static constexpr bool LOCAL_SWEEP = true; /**< Sweep in mesh-local space. */
  Vector axis = Y_AXIS;                     /**< Mesh-local sweep axis. */
  float fade_frames_min =
      4.0f; /**< Shortest per-face fade length, in frames. */
  float fade_frames_max =
      12.0f; /**< Longest per-face fade length, in frames. */
  uint32_t fade_seed = 0x9e3779b9u; /**< Per-transition seed for the per-face
                                       fade hash; rolled by retarget(). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    int fade = schedule_faded_sprite(timeline, std::move(draw_fn), duration,
                                     window, paused);
    inv_window = 1.0f / static_cast<float>(std::max(fade, 1));
    return duration;
  }
  void retarget(const Vector &v) {
    axis = v.normalized();
    fade_seed = static_cast<uint32_t>(hs::random()());
  }
  float face_offset(const Vector &center, int, int) const {
    return 0.5f * (1.0f + dot(center, axis));
  }
  /** @brief Per-face fade length as a window fraction: a stable hash of the
   * face index into the frame range, divided by the scheduled window. Computed
   * once per face (not per fragment), so it must stay a pure function of the
   * index, the seed and the sliders — the frame bounds are read live so a
   * mid-transition slider move takes effect on the next frame. */
  float face_fade_frac(int i) const {
    float t = hash01(static_cast<uint32_t>(i), fade_seed);
    float lo = std::min(1.0f, fade_frames_min * inv_window);
    float hi = std::min(1.0f, fade_frames_max * inv_window);
    return lo + (hi - lo) * t;
  }
  float face_phase(float phase, float offset, float fade_frac) const {
    float ff = std::max(fade_frac, 1e-4f);
    return hs::clamp((phase - offset * (1.0f - ff)) / ff, 0.0f, 1.0f);
  }
  bool visible(float phase) const { return fades_to_black(phase); }
  /** @brief Squared: alpha scales linear-light color, where a linear ramp
   * reads mostly-bright almost immediately; the square spreads the perceived
   * fade across the face's fade window. */
  float opacity(float phase) const { return phase * phase; }

private:
  /** @brief Reciprocal of the scheduled fade window in frames; set by
   * schedule(). The initializer only covers a query before the first
   * schedule(). */
  float inv_window = 1.0f / 64.0f;
};

/**
 * @brief An expanding shockwave erases the pattern outward from a point; its
 * echo redraws the new one.
 * @details Faces nearest the origin extinguish first, so the wave visibly
 * expands. Pairs naturally with the effect's ripple bursts sharing the origin.
 */
struct Shockwave : Base {
  static constexpr float BAND =
      0.3f;               /**< Wave-front softness, in phase units. */
  Vector origin = Y_AXIS; /**< World-space wave origin. */
  void retarget(const Vector &v) { origin = v; }
  float face_offset(const Vector &center, int, int) const {
    float angle = fast_acos(hs::clamp(dot(center, origin), -1.0f, 1.0f));
    return 1.0f - angle * (1.0f / PI_F);
  }
  float face_phase(float phase, float offset, float = 0.0f) const {
    return sweep_phase(phase, offset, BAND);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The pattern breaks down one topology class at a time: all faces of a
 * class fade together, each class fully gone before the next starts, in a
 * random class order reshuffled per transition; the new tessellation
 * reassembles class by class the same way.
 * @details Faces group by palette-slot class, so each color family vanishes as
 * a unit. Class windows are abutting equal slices of the phase range (linear,
 * not sweep_phase's eased front). The BLACK_DWELL slice nearest the swap is
 * held fully black so the last class completes before the incoming mesh
 * appears, instead of popping. reorder() derives the class count from the
 * per-face classes, so it can never disagree with the mesh.
 */
struct Breakdown : Base {
  static constexpr int MAX_CLASSES = 16; /**< rank[] capacity. */
  static constexpr float BLACK_DWELL =
      0.1f;            /**< Phase slice held all-black at the swap end. */
  int num_classes = 1; /**< Live class count, derived by reorder(). */
  uint8_t rank[MAX_CLASSES] =
      {}; /**< rank[class]: fade position; 0 vanishes first. */
  /**
   * @brief Derives the class count from the per-face classes and re-randomizes
   *        the fade order for the next transition.
   * @param face_classes Per-face palette-slot class ids (dense [0, k), the same
   *        values face_offset receives). num_classes is set to max+1, so a
   *        caller can never mis-declare it. Traps past MAX_CLASSES: a slot id
   *        comes from MeshPaletteBank, whose bank is far smaller, so an
   *        over-range id means raw topology classes were handed over instead.
   */
  template <typename Classes> void reorder(const Classes &face_classes) {
    int detected = 1;
    for (size_t i = 0; i < face_classes.size(); ++i) {
      int c = static_cast<int>(face_classes[i]) + 1;
      if (c > detected)
        detected = c;
    }
    HS_CHECK(detected <= MAX_CLASSES, "Breakdown: class count over rank[]");
    num_classes = detected;
    for (int i = 0; i < num_classes; ++i)
      rank[i] = static_cast<uint8_t>(i);
    hs::shuffle(rank, rank + num_classes);
  }
  float face_offset(const Vector &, int, int cls) const {
    if (num_classes <= 1)
      return 0.0f;
    int r = rank[(cls >= 0 && cls < num_classes) ? cls : 0];
    return static_cast<float>(num_classes - 1 - r) /
           static_cast<float>(num_classes - 1);
  }
  float face_phase(float phase, float offset, float = 0.0f) const {
    // Class windows tile [BLACK_DWELL, 1]; phase 1 stays the identity plateau.
    float band = (1.0f - BLACK_DWELL) / static_cast<float>(num_classes);
    return hs::clamp(
        (phase - BLACK_DWELL - offset * (1.0f - BLACK_DWELL - band)) / band,
        0.0f, 1.0f);
  }
  bool visible(float phase) const { return fades_to_black(phase); }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The whole mesh spins up around an axis until the POV display smears
 * it into bands, swaps at peak speed, and spins back down — a coin flip.
 * @details The warp is rigid, so there is no fold/overdraw hazard and the mesh
 * never fades; the swap hides in the motion blur.
 */
struct SpinFlip : Base {
  static constexpr float REVS = 3.0f; /**< Extra revolutions at peak spin. */
  Vector axis = Y_AXIS;               /**< Spin axis. */
  void retarget(const Vector &v) { axis = v.normalized(); }
  Vector warp(const Vector &v, float phase) const {
    float wind = 1.0f - phase;
    return rotate(v, make_rotation(axis, wind * wind * REVS * 2.0f * PI_F));
  }
};

/**
 * @brief Both palettes converge to molten gold around the swap, then the new
 * mesh blooms back into color.
 * @details Purely palette-domain: geometry never moves. A mild opacity dip at
 * the swap softens the topology pop while both meshes are monochrome.
 */
struct GoldConvergence : Base {
  Pixel gold = Color4(uint8_t{255}, uint8_t{196}, uint8_t{64})
                   .color; /**< Linear-space convergence color. */
  Color4 grade(Color4 c, float phase) const {
    return c.lerp(Color4(gold, c.alpha), 1.0f - phase);
  }
  float opacity(float phase) const { return 0.4f + 0.6f * phase; }
};

/**
 * @brief Stochastic ownership dissolve: each wireframe edge shows exactly one
 * of the two meshes, with the owned fraction tracking the phase.
 * @details The two draws receive complementary DissolveMasks (same threshold
 * and salt, opposite invert), so together they rasterize each edge once — a
 * two-mesh transition costs one wireframe's draw per frame instead of two,
 * which is what keeps heavy-pair crossfades inside one display window. Owned
 * edges draw at full opacity; the dissolve percept is the spatial mix ratio,
 * blurred by POV persistence. The salt folds a frame counter into the
 * per-transition seed so the pattern re-rolls every frame (temporal dither).
 * Unlike the other policies this one partitions rasterizer work (see
 * DissolveMask) rather than fragments in the shader; effects pass the masks to
 * Plot::Mesh::draw's edge-list overload themselves. Only that path takes a
 * mask, so a solid-mesh pair cannot dissolve.
 */
struct Dissolve : Base {
  uint32_t seed =
      0x9e3779b9u; /**< Per-transition seed; rolled by retarget(). */
  void retarget(const Vector &) {
    seed = static_cast<uint32_t>(hs::random()());
  }
  /**
   * @brief The two halves of one frame's ownership split.
   */
  struct MaskPair {
    DissolveMask incoming; /**< Mask the incoming mesh's draw takes. */
    DissolveMask outgoing; /**< Its complement, for the outgoing mesh's draw. */
  };

  /**
   * @brief Builds both halves of one frame's ownership split.
   * @param phase Transition phase in [0, 1]; the incoming mesh owns this
   *        fraction of the edges, the outgoing mesh the complement.
   * @param frame Monotonic frame counter (temporal dither; never wall time).
   * @return The complementary pair.
   * @details The masks partition only when they share a threshold, and
   * schedule() puts the two draws on independent sprites carrying independent
   * phases — so the pair is derived from one phase here and split by the
   * caller, rather than each half deriving its own.
   */
  MaskPair mask_pair(float phase, uint32_t frame) const {
    const uint32_t thr =
        static_cast<uint32_t>(hs::clamp(phase, 0.0f, 1.0f) * 65536.0f);
    const uint32_t salt = frame * 0x9E3779B9u ^ seed;
    return {DissolveMask{thr, salt, false}, DissolveMask{thr, salt, true}};
  }
  /** @brief Overlapping schedule, fixed at the full fade window: the two masks
   * partition the edges only while both meshes are on the timeline, so a
   * shorter overlap would leave the complement unlit (the masks keep the cost
   * at one mesh). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    return schedule_overlapped(timeline, std::move(draw_fn), duration, window,
                               -1, paused);
  }
};

} // namespace Segue
