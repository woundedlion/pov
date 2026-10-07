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

#include "animation/preset_departure.h"

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
 * Declaring an optional hook at the wrong signature is a compile error.
 *
 * reorder and mask_pair are contracts on the effect, not the draw path: a
 * NeedsClasses policy left un-reordered fades every class as one, and a Masked
 * policy drawn without its masks rasterizes both meshes.
 *
 * The per-face hooks and the fragment hooks are mutually exclusive: the
 * per-face draw path never calls fill or grade.
 *
 * A per-face policy may also declare `static constexpr bool LOCAL_SWEEP =
 * true` to order faces by the untransformed mesh instead of world-space
 * centers, so the front rides the mesh's rotation.
 */
namespace Segue {

/**
 * @brief Schedules one sprite with symmetric linear fade-in/fade-out.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to [0,
 * duration/2] so the in/out windows never collide.
 * @param paused Optional event-level pause gate.
 * @return The clamped fade length, from which each policy derives its own
 * return offset.
 * @details Traps when the timeline has no free slot.
 */
inline int schedule_faded_sprite(Timeline &timeline, SpriteFn draw_fn,
                                 int duration, int window,
                                 const bool *paused = nullptr) {
  HS_CHECK(Timeline::remaining() >= 1,
           "segue: the transition sprite needs a free timeline slot");
  int fade = hs::clamp(window, 0, std::max(duration / 2, 0));
  Animation::Sprite sprite(std::move(draw_fn), duration,
                           {.fade_in = {fade}, .fade_out = {fade}});
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
 * @param window Requested transition window in frames; clamped to [0,
 * duration/2] so the in/out windows never collide.
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
 * @param window Requested transition window in frames; clamped to [0,
 * duration/2] so the in/out windows never collide.
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
 * @brief Soft, eased sweep front.
 * @param phase Global segue phase in [0, 1].
 * @param offset Face's sweep ordering in [0, 1]; larger extinguishes earlier.
 * @param band Softness of the front, in phase units.
 * @return The face-local phase in [0, 1]: 1 everywhere at phase 1, 0
 * everywhere at phase 0, with faces crossing the front in offset order.
 * @details The sqrt ease accelerates the front through low phase, shortening
 * the dark interval around the swap. Endpoints stay exact.
 */
inline float sweep_phase(float phase, float offset, float band) {
  float p = std::sqrt(phase);
  return hs::clamp((p - offset * (1.0f - band)) / band, 0.0f, 1.0f);
}

/**
 * @brief Cull gate for a policy whose shading vanishes as phase falls to 0.
 * @param phase The phase the policy's opacity reads: the global segue phase, or
 * — for a per-face policy — the face-local phase of the brightest face.
 * @return Whether drawing at this phase can produce visible output.
 */
inline bool fades_to_black(float phase) { return phase > 0.005f; }

/**
 * @brief Identity hooks every segue inherits; a policy shadows only the hooks
 * its transition uses.
 */
struct Base {
  /**
   * @brief Whether consecutive sprites coexist, so the previous transition is
   * still drawing when the next one is scheduled.
   * @details A policy scheduling through schedule_overlapped() must set this
   * true. One policy instance serves both coexisting sprites: a per-face
   * policy may not overlap, a retarget() under overlap requires
   * RETARGET_SAFE_UNDER_OVERLAP, and an overlapping schedule() must not write
   * state its draw hooks read (unchecked).
   */
  static constexpr bool OVERLAPS = false;
  /**
   * @brief Whether the per-transition state retarget() rolls survives being
   * rewritten while an overlapping predecessor's sprite still reads it.
   * @details Read only for an OVERLAPS policy that declares retarget().
   */
  static constexpr bool RETARGET_SAFE_UNDER_OVERLAP = false;
  /** @brief Default scheduling: one sequential sprite (see
   * schedule_sequential). @p paused is the optional event-level pause gate. */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    return schedule_sequential(timeline, std::move(draw_fn), duration, window,
                               paused);
  }
  /** @brief Whether drawing at this phase can produce visible output. The
   * identity policy never culls. */
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
 * including the trailing pause gate. A policy shadowing schedule() with a
 * shorter signature hides Base's and fails this concept. */
template <typename S>
concept Schedulable =
    requires(S &s, Timeline &timeline, SpriteFn draw_fn, const bool *paused) {
      {
        s.schedule(timeline, std::move(draw_fn), 0, 0, paused)
      } -> std::same_as<int>;
    };

namespace detail {

/** @brief The function type of a const pointer-to-member, dropping the class
 * the hook was found on so an inherited and a shadowing declaration of the same
 * signature compare equal. */
template <typename M> struct HookSignature {};
template <typename C, typename R, typename... A>
struct HookSignature<R (C::*)(A...) const> {
  using type = R(A...);
};

} // namespace detail

/** @brief Whether a policy keeps Base's always-present hook signatures; fails
 * on a hook shadowed at a drifted signature (e.g. a `fill` taking its edge
 * distance by value, which still compiles at the call site). */
template <typename S>
concept HasPhaseHooks = requires(const S &s) {
  { s.visible(0.5f) } -> std::same_as<bool>;
  { s.opacity(0.5f) } -> std::same_as<float>;
  { s.face_fade_frac(0) } -> std::same_as<float>;
  requires std::same_as<
      typename detail::HookSignature<decltype(&S::fill)>::type,
      float(float &, float)>;
  requires std::same_as<
      typename detail::HookSignature<decltype(&S::grade)>::type,
      Color4(Color4, float)>;
};

/** @brief Whether a policy defines the per-face ordering hook, with or without
 * the face_phase that makes the set usable. */
template <typename S>
concept HasFaceOffset =
    requires(const S &s, const math::Vector &c) { s.face_offset(c, 0, 0); };

/** @brief Whether a policy defines the face-local phase hook at the arity the
 * per-face draw path calls. */
template <typename S>
concept HasFacePhase = requires(const S &s) { s.face_phase(0.5f, 0.5f, 0.1f); };

/** @brief Whether a policy defines the per-face hook set. */
template <typename S>
concept PerFace = HasFaceOffset<S> && HasFacePhase<S>;

/** @brief Whether a policy orders faces by topology class, so the effect must
 * hand it the per-face classes before each transition. */
template <typename S>
concept NeedsClasses = requires(S &s, const ArenaVector<uint16_t> &classes) {
  { s.reorder(classes) } -> std::same_as<void>;
};

/** @brief Whether a policy splits one frame's rasterizer work between the two
 * meshes with complementary edge-key masks. The effect calls mask_pair() itself
 * and hands the two halves to Plot::Mesh::draw's edge-list overload. */
template <typename S>
concept Masked = requires(const S &s) {
  typename S::MaskPair;
  { s.mask_pair(0.5f, 0u) } -> std::same_as<typename S::MaskPair>;
};

/** @brief Whether a policy warps vertices before the effect's ripple. */
template <typename S>
concept HasWarp = requires(const S &s, const math::Vector &v) {
  { s.warp(v, 0.5f) } -> std::same_as<math::Vector>;
};

/** @brief Whether a policy re-randomizes its per-transition state on a
 * direction. */
template <typename S>
concept HasRetarget = requires(S &s, const math::Vector &v) {
  { s.retarget(v) } -> std::same_as<void>;
};

/** @brief Whether a policy's LOCAL_SWEEP is a constant-usable bool. */
template <typename S>
concept LocalSweeps = requires {
  requires std::same_as<std::remove_cv_t<decltype(S::LOCAL_SWEEP)>, bool>;
  typename std::bool_constant<S::LOCAL_SWEEP>;
};

/**
 * @brief Signature-agnostic "the policy declares this hook name" probes.
 * @details Each carrier declares one hook or trait name; merged into a policy,
 * the name is ambiguous exactly when the policy declares it too — whatever
 * signature or type it carries, and including a template member the call-shaped
 * concepts cannot see. A final or non-class policy cannot be merged into and
 * reports false.
 */
namespace detail {

struct WarpName {
  void warp();
};
struct RetargetName {
  void retarget();
};
struct ReorderName {
  void reorder();
};
struct MaskPairName {
  void mask_pair();
};
struct FaceOffsetName {
  void face_offset();
};
struct FacePhaseName {
  void face_phase();
};
struct LocalSweepName {
  static constexpr int LOCAL_SWEEP = 0;
};

template <typename S, typename Name> struct Merged : S, Name {};

template <typename S>
concept Mergeable = std::is_class_v<S> && !std::is_final_v<S>;

} // namespace detail

/** @brief Whether a policy declares a `warp` hook of any signature. */
template <typename S>
concept DeclaresWarp = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::WarpName>::warp;
};

/** @brief Whether a policy declares a `retarget` hook of any signature. */
template <typename S>
concept DeclaresRetarget = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::RetargetName>::retarget;
};

/** @brief Whether a policy declares a `reorder` hook of any signature. */
template <typename S>
concept DeclaresReorder = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::ReorderName>::reorder;
};

/** @brief Whether a policy declares a `mask_pair` hook of any signature. */
template <typename S>
concept DeclaresMaskPair = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::MaskPairName>::mask_pair;
};

/** @brief Whether a policy declares a `face_offset` hook of any signature. */
template <typename S>
concept DeclaresFaceOffset = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::FaceOffsetName>::face_offset;
};

/** @brief Whether a policy declares a `face_phase` hook of any signature. */
template <typename S>
concept DeclaresFacePhase = detail::Mergeable<S> && !requires {
  &detail::Merged<S, detail::FacePhaseName>::face_phase;
};

/** @brief Whether a policy declares a `LOCAL_SWEEP` member of any type. */
template <typename S>
concept DeclaresLocalSweep = detail::Mergeable<S> && !requires {
  detail::Merged<S, detail::LocalSweepName>::LOCAL_SWEEP;
};

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
 * @details Phase is opacity. Consecutive sprites coexist for `overlap` frames,
 * rasterizing both meshes; at overlap 0 the schedule is sequential (a fade
 * through black).
 */
struct Crossfade : Base {
  static constexpr bool OVERLAPS = true;
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
  /** @brief Culls once the fade has taken the mesh to black. */
  bool visible(float phase) const { return fades_to_black(phase); }
  /** @brief Global alpha: the fade envelope itself. */
  float opacity(float phase) const { return phase; }
};

/**
 * @brief Faces contract to glowing points at their centers, then the next
 * tessellation blooms back out of the point field.
 */
struct IrisBloom : Base {
  static constexpr float SOFT =
      0.08f; /**< Soft rim width, in edge-distance units. */
  /**
   * @brief Keeps only the core deeper than the phase-driven inset, renormalized
   * to the full palette gradient.
   * @param t Fragment edge distance in [0, 1]; remapped in place across the
   * surviving core.
   * @param phase Transition phase; the core shrinks as it falls.
   * @return Coverage alpha in [0, 1]; 0 culls the fragment.
   */
  float fill(float &t, float phase) const {
    float inset = 1.0f - phase;
    if (t < inset - SOFT)
      return 0.0f;
    float cover = hs::clamp((t - (inset - SOFT)) / SOFT, 0.0f, 1.0f);
    t = hs::clamp((t - inset) / fmaxf(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief The fill drains until only a glowing band along the edges survives,
 * the meshes swap as lace, then the new fill floods back in.
 * @details The inverse mask of IrisBloom: fragments within the phase-driven
 * band of an edge survive.
 */
struct Lace : Base {
  static constexpr float SOFT =
      0.08f; /**< Soft band-edge width, in edge-distance units. */
  /**
   * @brief Keeps only the band within the phase-driven distance of an edge,
   * renormalized to the full palette gradient.
   * @param t Fragment edge distance in [0, 1]; remapped in place across the
   * surviving band.
   * @param phase Transition phase; the band narrows as it falls.
   * @return Coverage alpha in [0, 1]; 0 culls the fragment.
   */
  float fill(float &t, float phase) const {
    if (t > phase + SOFT)
      return 0.0f;
    float cover = hs::clamp((phase + SOFT - t) / SOFT, 0.0f, 1.0f);
    t = hs::clamp(t / fmaxf(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief A day/night line pinned to the mesh sweeps across it; when it reaches
 * a face, that face fades over a per-face random length in [fade_frames_min,
 * fade_frames_max] frames, capped by the scheduled fade window.
 * @details The front crosses over the fade window minus one face fade, so
 * face phases are exactly 1 at phase 1 and 0 at phase 0 for every fade length.
 */
struct TerminatorSweep : Base {
  static constexpr bool LOCAL_SWEEP = true; /**< Sweep in mesh-local space. */
  math::Vector axis = math::Y_AXIS;         /**< Mesh-local sweep axis. */
  float fade_frames_min =
      4.0f; /**< Shortest per-face fade length, in frames. */
  float fade_frames_max =
      12.0f; /**< Longest per-face fade length, in frames. */
  uint32_t fade_seed = 0x9e3779b9u; /**< Per-transition seed for the per-face
                                       fade hash; rolled by retarget(). */
  /**
   * @brief Schedules the sprite and caches the fade window the per-face phases
   * are measured against.
   * @param timeline Timeline receiving the sprite.
   * @param draw_fn Draws the mesh at the envelope phase.
   * @param duration Total frames the mesh is on screen.
   * @param window Requested fade length in frames; clamped to duration/2.
   * @param paused Optional event-level pause gate.
   * @return duration — consecutive sprites do not overlap.
   */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    int fade = schedule_faded_sprite(timeline, std::move(draw_fn), duration,
                                     window, paused);
    inv_window = 1.0f / static_cast<float>(std::max(fade, 1));
    return duration;
  }
  /**
   * @brief Aims the sweep and re-rolls the per-face fade hash for the next
   * transition.
   * @param v Sweep direction; normalized into the mesh-local axis.
   */
  void retarget(const math::Vector &v) {
    axis = v.normalized();
    fade_seed = static_cast<uint32_t>(hs::random()());
  }
  /**
   * @brief Orders faces along the sweep axis.
   * @return Position in [0, 1]: 1 at the axis's positive pole, which the front
   * reaches first.
   */
  float face_offset(const math::Vector &center, int, int) const {
    return 0.5f * (1.0f + math::dot(center, axis));
  }
  /** @brief Per-face fade length as a window fraction: a stable hash of the
   * face index into the frame range; the frame bounds may be unordered. */
  float face_fade_frac(int i) const {
    float t = math::hash01(static_cast<uint32_t>(i), fade_seed);
    float lo = min_fade_frac();
    float hi =
        fminf(1.0f, fmaxf(fade_frames_min, fade_frames_max) * inv_window);
    return lo + (hi - lo) * t;
  }
  /**
   * @brief Face-local phase behind the linear front.
   * @param phase Global segue phase in [0, 1].
   * @param offset The face's position along the sweep axis.
   * @param fade_frac The face's fade length as a window fraction.
   * @return The face's own phase in [0, 1]: 0 before the front arrives, 1 once
   * its fade has completed.
   */
  float face_phase(float phase, float offset, float fade_frac) const {
    float ff = fmaxf(fade_frac, 1e-4f);
    return hs::clamp((phase - offset * (1.0f - ff)) / ff, 0.0f, 1.0f);
  }
  /** @brief Culls once the last face the front reaches, at the shortest
   * fade, has gone black. */
  bool visible(float phase) const {
    return fades_to_black(face_phase(phase, 0.0f, min_fade_frac()));
  }
  /** @brief Squared phase, for a perceptually even fade in linear light. */
  float opacity(float phase) const { return phase * phase; }

private:
  /** @brief Shortest per-face fade length as a window fraction. */
  float min_fade_frac() const {
    return fminf(1.0f, fminf(fade_frames_min, fade_frames_max) * inv_window);
  }
  /** @brief Reciprocal of the scheduled fade window in frames; set by
   * schedule(). The initializer only covers a query before the first
   * schedule(). */
  float inv_window = 1.0f / 64.0f;
};

/**
 * @brief An expanding shockwave erases the pattern outward from a point; its
 * echo redraws the new one.
 */
struct Shockwave : Base {
  static constexpr float BAND =
      0.3f; /**< Wave-front softness, in phase units. */
  math::Vector origin = math::Y_AXIS; /**< Unit-length world-space wave origin;
                               face_offset's acos orders faces only for a unit
                               vector. */
  /**
   * @brief Moves the wave origin for the next transition.
   * @param v World-space direction the wave expands from; normalized into the
   * origin.
   */
  void retarget(const math::Vector &v) { origin = v.normalized(); }
  /**
   * @brief Orders faces by angular distance from the origin.
   * @return Position in [0, 1]: 1 at the origin, which extinguishes first.
   */
  float face_offset(const math::Vector &center, int, int) const {
    float angle =
        math::fast_acos(hs::clamp(math::dot(center, origin), -1.0f, 1.0f));
    return 1.0f - angle * (1.0f / math::PI_F);
  }
  /**
   * @brief Face-local phase behind the eased wave front (see sweep_phase); the
   * fade fraction is unused, every face crossing over the same BAND.
   * @return The face's own phase in [0, 1].
   */
  float face_phase(float phase, float offset, float = 0.0f) const {
    return sweep_phase(phase, offset, BAND);
  }
  /** @brief Culls once the last face the front reaches has gone black. */
  bool visible(float phase) const {
    return fades_to_black(face_phase(phase, 0.0f));
  }
  /** @brief Global alpha, applied to the face-local phase: a linear fade. */
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The pattern breaks down one topology class at a time: all faces of a
 * class fade together, each class fully gone before the next starts, in a
 * random class order reshuffled per transition; the new tessellation
 * reassembles class by class the same way.
 * @details Faces group by palette-slot class. The BLACK_DWELL slice nearest
 * the swap is held fully black.
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
   *        values face_offset receives). num_classes is set to max+1. Traps
   *        past MAX_CLASSES, which signals raw topology classes rather than
   *        palette slots.
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
  /**
   * @brief Orders faces by their class's draw of the shuffled fade order; an
   * out-of-range class shares class 0's rank.
   * @return Position in [0, 1]: 1 for the class that vanishes first.
   */
  float face_offset(const math::Vector &, int, int cls) const {
    if (num_classes <= 1)
      return 0.0f;
    int r = rank[(cls >= 0 && cls < num_classes) ? cls : 0];
    return static_cast<float>(num_classes - 1 - r) /
           static_cast<float>(num_classes - 1);
  }
  /**
   * @brief Face-local phase within the class's own slice of the phase range;
   * the fade fraction is unused, a class fading over its whole slice.
   * @return The face's own phase in [0, 1].
   */
  float face_phase(float phase, float offset, float = 0.0f) const {
    // Class windows tile [BLACK_DWELL, 1]; phase 1 stays the identity plateau.
    float band = (1.0f - BLACK_DWELL) / static_cast<float>(num_classes);
    return hs::clamp(
        (phase - BLACK_DWELL - offset * (1.0f - BLACK_DWELL - band)) / band,
        0.0f, 1.0f);
  }
  /** @brief Culls once the last class has faded, including the BLACK_DWELL
   * slice. */
  bool visible(float phase) const {
    return fades_to_black(face_phase(phase, 0.0f));
  }
  /** @brief Global alpha, applied to the face-local phase: a linear fade. */
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The whole mesh spins up around an axis until the POV display smears
 * it into bands, swaps at peak speed, and spins back down — a coin flip.
 * @details The warp is rigid and the mesh never fades.
 */
struct SpinFlip : Base {
  static constexpr float REVS = 3.0f; /**< Extra revolutions at peak spin. */
  math::Vector axis = math::Y_AXIS;   /**< Spin axis. */
  /**
   * @brief Aims the spin for the next transition.
   * @param v Spin direction; normalized into the axis.
   */
  void retarget(const math::Vector &v) { axis = v.normalized(); }
  /**
   * @brief Winds a vertex around the axis, fastest at phase 0 and stationary
   * at phase 1.
   * @param v Unit-sphere vertex, before the effect's ripple.
   * @param phase Transition phase in [0, 1].
   * @return The rotated vertex.
   */
  math::Vector warp(const math::Vector &v, float phase) const {
    float wind = 1.0f - phase;
    return math::rotate(
        v, math::make_rotation(axis, wind * wind * REVS * 2.0f * math::PI_F));
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
  /**
   * @brief Pulls the palette color toward gold as phase falls.
   * @param c Color from the palette lookup.
   * @param phase Transition phase in [0, 1]; 0 is fully gold.
   * @return The regraded color, alpha untouched.
   */
  Color4 grade(Color4 c, float phase) const {
    return c.lerp(Color4(gold, c.alpha), 1.0f - phase);
  }
  /** @brief Global alpha: dips at the swap, never to black. */
  float opacity(float phase) const { return 0.4f + 0.6f * phase; }
};

/**
 * @brief Stochastic wireframe dissolve with complementary edge-key ownership
 * and an incoming share tracking the phase.
 * @details The two draws receive complementary DissolveMasks (same threshold
 * and salt, opposite invert), which partition shared vertex-index keys. Owned
 * edges draw at full opacity. The salt folds a frame counter into the
 * per-transition seed so the pattern re-rolls every frame. A solid-mesh pair
 * cannot dissolve.
 */
struct Dissolve : Base {
  static constexpr bool OVERLAPS = true;
  /** @brief One mask_pair() call feeds both halves of a frame, so a
   * mid-overlap seed rewrite keeps the split complementary. */
  static constexpr bool RETARGET_SAFE_UNDER_OVERLAP = true;
  uint32_t seed =
      0x9e3779b9u; /**< Per-transition seed; rolled by retarget(). */
  /** @brief Re-rolls the ownership pattern for the next transition; the
   * direction is unused. */
  void retarget(const math::Vector &) {
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
   * @param phase Transition phase in [0, 1], scaled into the incoming mask
   *        hash threshold; the outgoing mask takes the complementary keys.
   * @param frame Monotonic frame counter (temporal dither; never wall time).
   * @return The complementary pair.
   * @details The masks partition only when they share a threshold, so both
   * halves derive from the one phase.
   */
  MaskPair mask_pair(float phase, uint32_t frame) const {
    const uint32_t thr =
        static_cast<uint32_t>(hs::clamp(phase, 0.0f, 1.0f) * 65536.0f);
    const uint32_t salt = frame * 0x9E3779B9u ^ seed;
    return {DissolveMask{thr, salt, false}, DissolveMask{thr, salt, true}};
  }
  /** @brief Overlapping schedule, fixed at the full fade window: the masks
   * partition the edges only while both meshes are on the timeline. */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window,
               const bool *paused = nullptr) {
    return schedule_overlapped(timeline, std::move(draw_fn), duration, window,
                               -1, paused);
  }
};

/**
 * @brief A pack of policies with roster-wide conformance folds.
 * @tparam Ts The policies.
 */
template <typename... Ts> struct PolicyList {
  /** @brief Whether every listed policy keeps the expected hook signatures. */
  static constexpr bool CONFORMING =
      ((Schedulable<Ts> && HasPhaseHooks<Ts>) && ...);

  /** @brief Whether every listed policy that is per-face also schedules
   * sequentially. */
  static constexpr bool SEQUENTIAL_PER_FACE =
      ((!PerFace<Ts> || !Ts::OVERLAPS) && ...);

  /** @brief Whether every listed policy that retargets per-transition state
   * either schedules sequentially or declares the state safe to rewrite
   * mid-overlap. */
  static constexpr bool RETARGET_SURVIVES_OVERLAP =
      ((!DeclaresRetarget<Ts> || !Ts::OVERLAPS ||
        Ts::RETARGET_SAFE_UNDER_OVERLAP) &&
       ...);

  /** @brief Whether every listed policy that declares LOCAL_SWEEP declares it
   * as a constant-usable bool. */
  static constexpr bool LOCAL_SWEEPS_TYPED =
      ((!DeclaresLocalSweep<Ts> || LocalSweeps<Ts>) && ...);

  /** @brief Whether every listed policy is a non-final class; a final one
   * answers false to every Declares* probe and passes them vacuously. */
  static constexpr bool MERGEABLE = (detail::Mergeable<Ts> && ...);

  /**
   * @brief Invokes @p fn once per policy, on a default-constructed instance.
   * @param fn Callable taking any policy.
   */
  template <typename F> static void for_each(F &&fn) { (fn(Ts{}), ...); }
};

/** @brief Every library policy, Base first. */
using AllPolicies =
    PolicyList<Base, Crossfade, IrisBloom, Lace, TerminatorSweep, Shockwave,
               Breakdown, SpinFlip, GoldConvergence, Dissolve>;

static_assert(AllPolicies::CONFORMING,
              "a segue policy shadows schedule(), visible(), opacity() or "
              "face_fade_frac(), fill() or grade() at a drifted signature");

static_assert(AllPolicies::SEQUENTIAL_PER_FACE,
              "a per-face segue must schedule sequentially: schedule() and "
              "retarget() rewrite the single policy instance's per-transition "
              "state, which an overlapping predecessor's sprite is still "
              "reading");

static_assert(AllPolicies::RETARGET_SURVIVES_OVERLAP,
              "an overlapping segue's retarget() rewrites the single policy "
              "instance's per-transition state, which the outgoing sprite is "
              "still reading: schedule sequentially, or set "
              "RETARGET_SAFE_UNDER_OVERLAP once the rewrite is shown harmless");

static_assert(AllPolicies::LOCAL_SWEEPS_TYPED,
              "a segue's LOCAL_SWEEP must be static constexpr bool: the "
              "per-face draw path reads it as a constant expression");

static_assert(AllPolicies::MERGEABLE,
              "a segue must be a non-final class: the Declares* probes merge a "
              "name carrier into it, and a final policy answers false to every "
              "one of them, passing them vacuously");

} // namespace Segue
