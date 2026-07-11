/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>

#include <new> // std::launder
#include <type_traits>
#include "math/3dmath.h"
#include "engine/platform.h"
#include "vendor/FastNoiseLite.h"
#include "engine/generators.h"
#include "math/geometry.h"
#include "engine/concepts.h" // Canvas, PlotFn/ScalarFn/TimerFn
#include "mesh/mesh.h"     // MeshOps
#include "engine/memory.h"
#include "mesh/spatial.h"
#include "engine/static_circular_buffer.h"
#include "math/rotate.h"
#include "engine/util.h" // wrap_t
#include "math/easing.h"

/**
 * @brief Non-templated interface for all animations.
 * @details Enables virtual dispatch in Timeline without std::visit.
 */
class IAnimation {
public:
  /**
   * @brief Virtual destructor for safe polymorphic deletion.
   */
  virtual ~IAnimation() = default;

  /**
   * @brief Advances the animation by one frame.
   * @param canvas The canvas buffer to draw into.
   */
  virtual void step(Canvas &canvas) = 0;

  /**
   * @brief Reports whether the animation has finished.
   * @return True if the animation is done.
   */
  virtual bool done() const = 0;

  /**
   * @brief Reports whether the animation repeats after finishing.
   * @return True if the animation repeats.
   */
  virtual bool repeats() const = 0;

  /**
   * @brief Reports whether the animation has a finite duration.
   * @return True if it can reach done() on its own; false for an indefinite
   * (perpetual) animation that only ends via cancel(). Defaults to true.
   */
  virtual bool is_finite() const { return true; }

  /**
   * @brief Reports whether the animation was explicitly canceled.
   * @return True if cancel() drove it to done(); false for a natural end.
   * @details Lets Timeline distinguish a deliberate teardown from a natural
   * completion (the pin-completion guard exempts cancellation). Defaults to
   * false for animations with no cancel concept.
   */
  virtual bool is_canceled() const { return false; }

  /**
   * @brief Resets the animation to its start.
   */
  virtual void rewind() = 0;

  /**
   * @brief Fires the registered completion callback.
   */
  virtual void post_callback() = 0;

  /**
   * @brief Collapses the owned Orientation's motion-blur history.
   * @details Override in types that own an Orientation to collapse it.
   */
  virtual void collapse_orientation() {}

  /**
   * @brief Identity of the owned Orientation.
   * @return Pointer identifying the owned Orientation, or nullptr if none.
   * @details Used by Timeline to collapse each distinct Orientation exactly
   * once per frame, so animations sharing one Orientation compose their
   * motion-blur history instead of clobbering it.
   */
  virtual const void *orientation_id() const { return nullptr; }
};
namespace Animation {

/**
 * @brief Defines the coordinate space for rotations.
 */
enum class Space {
  World, /**< Rotations are applied relative to the fixed world axes
            (Pre-multiply). */
  Local  /**< Rotations are applied relative to the object's current axes
            (Post-multiply). */
};

/**
 * @brief Base class for all animations, providing core timing and state
 * management.
 * @tparam Derived The class inheriting from Animation (Curiously Recurring
 * Template Pattern).
 */
template <typename Derived> class AnimationBase : public IAnimation {
public:
  /**
   * @brief Cancels the animation on the next step.
   */
  void cancel() { canceled = true; }
  /**
   * @brief Checks if the animation is finished or canceled.
   * @return True if done.
   */
  bool done() const override {
    return canceled || (duration >= 0 && t >= static_cast<uint32_t>(duration));
  }
  /**
   * @brief Checks if the animation should repeat after finishing.
   * @return True if repeating.
   * @details A canceled animation must not repeat: done() is permanently true
   * once canceled, so without the !canceled guard Timeline::step would rewind it
   * and re-fire its .then() every frame, never removing it — a per-frame zombie.
   */
  bool repeats() const override { return repeat && !canceled; }
  /**
   * @brief Reports whether this animation was canceled (vs. ran to its end).
   * @return True once cancel() has been called.
   */
  bool is_canceled() const override { return canceled; }
  /**
   * @brief Reports whether the animation has a finite duration.
   * @return True if it can reach done() on its own (duration >= 0); false for an
   * indefinite (perpetual) animation that only ends via cancel().
   */
  bool is_finite() const override { return duration >= 0; }
  /**
   * @brief Advances the animation state by one frame.
   * @details The canvas buffer is unused by the base class and passed through
   * to derived classes.
   */
  void step(Canvas &) override { t++; }

  /**
   * @brief Resets the internal timer to zero.
   */
  void rewind() override { t = 0; }

  /**
   * @brief Sets a callback fired at the end of each completion cycle.
   *
   * Fires when the animation reaches done(): once for a one-shot, per cycle for
   * a repeating one, once per frame for a Driver. RandomTimer/PeriodicTimer
   * never reach done() (duration=-1), so they fire it from their own step().
   * Do not attach a one-shot callback to a repeating target.
   *
   * Single post slot: then() traps (HS_CHECK) rather than overwrite an existing
   * callback.
   *
   * @param callback The function to execute at each completion.
   * @return LValue Reference to the derived animation object.
   */
  Derived &then(Fn<void(), 24> callback) & {
    HS_CHECK(!post);
    post = std::move(callback);
    return static_cast<Derived &>(*this);
  }

  /**
   * @brief Sets a per-cycle completion callback (RValue overload).
   *
   * See the lvalue overload for the per-cycle semantics across one-shot,
   * repeating, and Driver targets.
   * @param callback The function to execute at each completion.
   * @return RValue Reference to the derived animation object.
   */
  Derived &&then(Fn<void(), 24> callback) && {
    HS_CHECK(!post);
    post = std::move(callback);
    return static_cast<Derived &&>(*this);
  }

  /**
   * @brief Executes the registered post-completion callback.
   */
  void post_callback() override {
    if (post)
      post();
  }

protected:
  /**
   * @brief Constructor for the base animation class.
   * @param duration Total number of frames the animation should run (-1 for
   * indefinite).
   * @param repeat If true, the animation rewinds and restarts when finished.
   */
  AnimationBase(int duration, bool repeat)
      : duration(duration == 0 ? 1 : duration), repeat(repeat),
        canceled(false) {
    // -1 is the sole perpetual sentinel; any other negative would make done()
    // permanently false so a one-shot never completes or fires .then().
    HS_CHECK(duration >= 0 || duration == -1,
             "AnimationBase duration must be >= 0 or -1 (perpetual)");
  }

  /**
   * @brief Default constructor: an indefinite, non-repeating animation.
   */
  AnimationBase() : duration(-1), repeat(false), canceled(false) {}

  /**
   * @brief Evaluates a wired pause flag.
   * @param flag Optional pointer to an effect's pause bool (may be null).
   * @return True when a non-null flag points to a set bool.
   * @details Shared gate for the pausable animations (Mutation/Driver/Lerp/Sprite)
   * so pause semantics live in one place.
   */
  static bool is_paused(const bool *flag) { return flag && *flag; }

  int duration; /**< Total length of the animation in frames. */
  bool repeat;  /**< Flag indicating if the animation should repeat. */
  uint32_t t = 0; /**< Internal frame counter. Finite animations bound it by
                     `duration`; perpetual ones (duration == -1) increment it
                     forever. Unsigned so `t++` and the `t >= next` trigger
                     comparisons wrap with defined behavior past 2^32 frames
                     (~552 days at 90 fps) rather than incurring signed-overflow
                     UB. Restart before then if uptime can approach that bound. */

private:
  bool canceled;       /**< Flag to signal immediate cancellation. */
  Fn<void(), 24> post; /**< Callback executed when the animation finishes. */
};

} // namespace Animation

// Internal fragments, in dependency order: every animation type derives from
// AnimationBase above; sprites builds on trails, mesh on timeline and sprites.
#define HS_ANIMATION_INTERNAL
#include "animation/timers.h"
#include "animation/params.h"
#include "animation/motion.h"
#include "animation/trails.h"
#include "animation/sprites.h"
#include "animation/timeline.h"
#include "animation/mesh.h"
#undef HS_ANIMATION_INTERNAL

// Device inline-storage budget audit. The per-type static_assert in
// Timeline::add only fires for types actually add()ed in this build, so pin the
// largest concrete animation type here to check it in every build that includes
// this header. Compared against MAX_ANIM_SIZE (the real 112 B device budget on
// the 32-bit WASM/device build, a no-op on the wider native host).
/** @brief Largest sizeof over a pack of types. */
template <typename... Ts> constexpr size_t largest_sizeof() {
  return std::max({sizeof(Ts)...});
}

// Add every new non-templated Animation type to this pack; the audit folds over
// it, so the static_assert tracks the list automatically.
constexpr size_t LARGEST_CONCRETE_ANIM_SIZE = largest_sizeof<
    Animation::RandomTimer, Animation::PeriodicTimer, Animation::Transition,
    Animation::Mutation, Animation::Driver, Animation::Lerp, Animation::Sprite,
    Animation::ColorWipe, Animation::MobiusFlow, Animation::MobiusWarp,
    Animation::MobiusWarpCircular, Animation::MeshMorph,
    Animation::MobiusWarpEvolving, Animation::Ripple, Animation::Noise>();
static_assert(LARGEST_CONCRETE_ANIM_SIZE <= TimelineEvent::MAX_ANIM_SIZE,
              "A concrete animation type exceeds the TimelineEvent inline-storage "
              "budget (on the 32-bit WASM/device build MAX_ANIM_SIZE is the "
              "112-byte device budget); shrink the type or raise the budget.");
