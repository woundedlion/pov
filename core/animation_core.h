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
#include "3dmath.h"
#include "platform.h"
#include "FastNoiseLite.h"
#include "generators.h"
#include "geometry.h"
#include "concepts.h" // Canvas, PlotFn/ScalarFn/TimerFn
#include "mesh.h"     // MeshOps
#include "memory.h"
#include "spatial.h"
#include "static_circular_buffer.h"
#include "rotate.h"
#include "util.h" // wrap_t
namespace Animation {
class Noise;
} // namespace Animation

/**
 * @brief Represents a customizable path.
 * @tparam W Display width (carried for downstream sizing).
 * @tparam RESOLUTION Capacity of the internal point ring buffer.
 * @details Retains an internal buffer for state, but draws to the pipeline.
 */
template <int W, int RESOLUTION = 1024> class Path {
public:
  /**
   * @brief Constructs an empty path.
   */
  Path() {}

  /**
   * @brief Appends a procedurally generated segment to the path.
   * @param plot The function generating points.
   * @param domain The input domain scale.
   * @param samples The number of sample intervals to take (>= 1); the loop
   *        emits samples + 1 points so the endpoint easing(1.0) is always hit.
   * @param easing The easing function for sample distribution.
   * @return Reference to self for chaining.
   */
  Path &append_segment(PlotFn plot, float domain, int samples,
                       ScalarFn easing) {
    // samples >= 1 also keeps the t / samples divide below non-zero.
    HS_CHECK(samples >= 1);
    // Account for the pop_back below: a non-empty path drops its last point
    // before appending samples + 1, so the final size is (size - 1) + samples + 1.
    size_t retained = points.is_empty() ? points.size() : points.size() - 1;
    HS_CHECK(retained + static_cast<size_t>(samples) + 1 <= RESOLUTION);
    if (!points.is_empty())
      points.pop_back();
    for (int t = 0; t <= samples; t++) {
      points.push_back(plot(easing(static_cast<float>(t) / samples) * domain));
    }
    return *this;
  }

  /**
   * @brief Retrieves a point along the path by interpolation.
   * @param t Progress along the path (0.0 to 1.0).
   * @return The interpolated vector.
   */
  Vector get_point(float t) const {
    if (points.is_empty())
      return Vector(0, 0, 0);
    // Clamp: a negative t makes raw_index negative, and casting that to size_t
    // is UB (t > 1 is caught by the i >= size-1 guard below, t < 0 is not).
    t = hs::clamp(t, 0.0f, 1.0f);
    float raw_index = t * (points.size() - 1);
    size_t i = static_cast<size_t>(raw_index);
    float f = raw_index - i;
    if (i >= points.size() - 1)
      return points.back();
    const Vector &p1 = points[i];
    const Vector &p2 = points[i + 1];
    return p1 * (1.0f - f) + p2 * f;
  }

  /**
   * @brief Collapses the path to its newest point only.
   */
  void collapse() {
    // Clear in place rather than assigning a fresh buffer: a full
    // StaticCircularBuffer temporary (~12.3 KB) overflows the WASM 8 KB stack.
    if (points.size() > 1) {
      Vector last = points.back();
      points.clear();
      points.push_back(last);
    }
  }

private:
  StaticCircularBuffer<Vector, RESOLUTION> points;
};

/**
 * @brief Represents a path defined by a single procedural function.
 * @details Matches the interface of Path for use in Motion animations.
 */
struct ProceduralPath {
  PlotFn f; /**< The procedural function mapping t to a point. */

  /**
   * @brief Constructs an empty procedural path.
   */
  ProceduralPath() = default;

  /**
   * @brief Constructs a procedural path from a plotting function.
   * @param path_fn Function mapping a parameter to a point on the path.
   */
  ProceduralPath(PlotFn path_fn) : f(std::move(path_fn)) {}

  /**
   * @brief Evaluates the path at a given parameter.
   * @param t Path parameter.
   * @return The point produced by the procedural function at t.
   */
  Vector get_point(float t) const { return f(t); }
};

#include "easing.h"

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
   * @details Lets Timeline distinguish a deliberate teardown (a held handle
   * calling cancel()) from a natural completion. Only the latter is misuse for a
   * pinned event, so the pin-completion guard exempts cancellation. Defaults to
   * false for animations that have no cancel concept.
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
   * once canceled, so without the !canceled guard Timeline::step would see
   * done()&&repeats(), rewind it, fire its .then() every frame, and never remove
   * it — a per-frame zombie. Dropping repeat here routes cancel through the
   * removal branch (fires post_callback once, then erased), exactly as a
   * canceled one-shot already behaves.
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
   * "Completion" is per-cycle, not strictly once. The callback runs every time
   * the animation reaches done():
   *   - One-shot animation (repeat=false): fires exactly once, then the
   *     animation is removed.
   *   - Repeating animation (repeat=true): fires at the end of every cycle,
   *     just after the timer rewinds. RandomTimer/PeriodicTimer never reach
   *     done() (duration=-1, self-resetting), so they fire this hook directly
   *     from their own step() on each trigger to keep the contract.
   *   - Driver (perpetual, one-frame cycle): fires once per frame, since each
   *     step completes a cycle. See the Driver class doc.
   * This is one polymorphic hook by design — per-cycle is the only framing that
   * generalizes across all three, since a driver's "cycle" is a single
   * increment. Do not attach a one-shot callback to a repeating target.
   *
   * There is a single post slot: then() refuses (HS_CHECK) to overwrite an
   * already-registered callback rather than silently dropping it — e.g. the
   * Transformer attaches a slot-recycling callback to every spawned animation,
   * so a second then() on that handle would leak the slot.
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
                     `duration`; perpetual ones (duration == -1:
                     RandomTimer/PeriodicTimer/Driver/RandomWalk) increment it
                     forever. Unsigned so `t++` and the `t >= next` trigger
                     comparisons wrap with defined behavior past 2^32 frames
                     (~552 days at 90 fps) rather than incurring signed-overflow
                     UB. Restart before then if uptime can approach that bound. */

private:
  bool canceled;       /**< Flag to signal immediate cancellation. */
  Fn<void(), 24> post; /**< Callback executed when the animation finishes. */
};

} // namespace Animation
