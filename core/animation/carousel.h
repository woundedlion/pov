/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file carousel.h
 * @brief Animation fragment: MeshCarousel, the double-buffered mesh slot pair.
 */

#include "color/color.h"

/**
 * @brief A double-buffered pair of persistent MeshState slots, the
 *        arena-compaction primitives effects need to swap between them, and a
 *        pluggable compile-time segue.
 * @tparam SegueT Segue policy (see namespace Segue) behind schedule_segue().
 * Clients that run their own transition animations (e.g. an `OpLeg`) keep
 * the default and simply never call it.
 * @details Holds two MeshState slots in `persistent_arena` and a front/back
 * index. Effects own generation and drawing (generate into a slot, flip the
 * front index, reclaim the old slot); the segue owns transition scheduling.
 *
 * Usage:
 *   MeshCarousel<Segue::Crossfade> carousel;  // in effect members
 *
 *   // Build the initial shape directly into the front slot:
 *   carousel.current().clear();
 *   MeshOps::compile(mesh, carousel.current(), persistent_arena,
 * scratch_arena_a);
 *
 *   // To transition: generate into the back slot, flip, then let the segue
 *   // schedule the animation via schedule_segue (see
 *   // IslamicStars::spawn_shape for the pattern).
 */
template <typename SegueT = Segue::Crossfade> class MeshCarousel {
  static_assert(!Segue::PerFace<SegueT> ||
                    !Segue::SHADOWS_FRAGMENT_HOOKS<SegueT>,
                "a per-face segue's draw path never calls fill/grade");
  static_assert(Segue::Schedulable<SegueT>,
                "a segue's schedule() must take (timeline, draw_fn, duration, "
                "window, paused)");

public:
  /**
   * @brief Constructs an empty carousel with front slot index 0.
   */
  MeshCarousel() {}

  /**
   * @brief Gets the currently visible (front) mesh.
   * @return Const reference to the front MeshState slot.
   */
  const MeshState &current() const { return slots[front]; }

  /**
   * @brief Gets the currently visible (front) mesh (mutable).
   * @return Mutable reference to the front MeshState slot.
   */
  MeshState &current() { return slots[front]; }

  /**
   * @brief Direct slot access by index (for effects that need both).
   * @param i Slot index (0 or 1).
   * @return Const reference to the requested MeshState slot.
   */
  const MeshState &slot(int i) const {
    HS_CHECK(i == 0 || i == 1, "MeshCarousel slot index out of range");
    return slots[i];
  }

  /**
   * @brief Direct slot access by index (mutable).
   * @param i Slot index (0 or 1).
   * @return Mutable reference to the requested MeshState slot.
   */
  MeshState &slot(int i) {
    HS_CHECK(i == 0 || i == 1, "MeshCarousel slot index out of range");
    return slots[i];
  }

  /**
   * @brief Gets the front slot index (for capture in lambdas).
   * @return The index (0 or 1) of the front slot.
   */
  int front_index() const { return front; }

  /**
   * @brief Manually sets the front index (for effects that manage transitions
   * themselves).
   * @param idx The new front slot index (0 or 1).
   */
  void set_front(int idx) {
    HS_CHECK(idx == 0 || idx == 1, "MeshCarousel front index out of range");
    front = idx;
  }

  /**
   * @brief Schedules the segue's transition animation for the (already
   * front-flipped) incoming mesh.
   * @param timeline Timeline receiving the segue's animation.
   * @param slot Incoming mesh slot captured by @p draw_fn.
   * @param draw_fn Draws the incoming mesh; the float argument is the segue's
   * phase (opacity for Segue::Crossfade).
   * @param duration Total frames the mesh is on screen.
   * @param window Transition window length in frames, segue-interpreted.
   * @param paused Optional event-level pause gate handed to the scheduled
   * sprite; null leaves the transition unpausable.
   * @return Frames after which the effect should begin the next transition.
   */
  int schedule_segue(Timeline &timeline, int slot, SpriteFn draw_fn,
                     int duration, int window, const bool *paused = nullptr) {
    HS_CHECK(slot == front,
             "MeshCarousel segue scheduled before incoming slot flip");
    return segue_policy.schedule(timeline, std::move(draw_fn), duration, window,
                                 paused);
  }

  /**
   * @brief The carousel's segue policy instance (holds per-transition state
   * such as a sweep axis or wave origin).
   */
  SegueT &segue() { return segue_policy; }
  /** @brief Const view of the segue policy instance. */
  const SegueT &segue() const { return segue_policy; }

  /**
   * @brief Frees the back slot and compacts, preserving only the front slot.
   * @tparam AfterReset Callable type invoked as `void(Arena&)`.
   * @param after_reset Callback run immediately after the reset, while the
   * front slot is still evacuated.
   * @details Runs `after_reset(persistent_arena)` immediately after the reset —
   * while the front slot is still evacuated — so the caller can re-bake
   * effect-owned persistent data (e.g. a palette bank) into the fresh arena
   * *before* the front mesh is restored on top of it. Use when only the visible
   * (front) shape must survive a regeneration of the back slot. Legal while no
   * sprite still draws the back slot — an overlapping segue's outgoing sprite
   * draws the slot it was spawned into, which set_front() already made the
   * front one, so it is the slot this preserves.
   */
  template <typename AfterReset>
  void compact_keep_front(AfterReset after_reset) {
    int back = 1 - front;
    slots[back] = MeshState();
    Persist<MeshState> p(slots[front], scratch_arena_b, persistent_arena);
    release_gamut_lut();
    persistent_arena.reset();
    after_reset(persistent_arena);
  }

  /**
   * @brief Frees both slots and compacts, evacuating nothing.
   * @tparam AfterReset Callable type invoked as `void(Arena&)`.
   * @param after_reset Callback run immediately after the reset.
   * @details Only legal while no sprite is still drawing either slot — with a
   * sequential segue the outgoing sprite has finished by the time the next
   * transition is scheduled. Callers that regenerate both slots before the next
   * draw reclaim a whole MeshState over compact_keep_front.
   */
  template <typename AfterReset> void compact_drop_all(AfterReset after_reset) {
    slots[0] = MeshState();
    slots[1] = MeshState();
    release_gamut_lut();
    persistent_arena.reset();
    after_reset(persistent_arena);
  }

private:
  MeshState slots[2]; /**< Front/back double-buffered mesh slots. */
  int front = 0;      /**< Index (0 or 1) of the visible front slot. */
  /** Segue policy instance; per-transition state lives here. */
  SegueT segue_policy;
};
