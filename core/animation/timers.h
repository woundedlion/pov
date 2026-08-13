/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file timers.h
 * @brief Animation fragment: the RandomTimer and PeriodicTimer callback
 * schedulers.
 */

namespace Animation {

/**
 * @brief Trigger loop shared by the callback schedulers.
 * @tparam Derived Timer supplying reset(), which schedules the next trigger
 * from the current frame counter.
 */
template <typename Derived> class TimerBase : public AnimationBase<Derived> {
public:
  /**
   * @brief Steps the timer, calling the function if the delay has elapsed.
   * @param canvas The canvas buffer (forwarded to the timer callback).
   */
  void step(Canvas &canvas) override {
    AnimationBase<Derived>::step(canvas);
    if (this->t < next)
      return;
    f(canvas);
    if (this->repeat) {
      static_cast<Derived *>(this)->reset();
      // A repeating timer never reaches done(), so fire the per-cycle .then()
      // directly to honor the contract.
      this->post_callback();
    } else {
      // finish(), not cancel(): Timeline's pin guard exempts cancellation, so
      // a canceled one-shot would be destroyed under a retained pointer with no
      // diagnostic.
      this->finish();
    }
  }

  /** @brief A one-shot timer is removed on its single trigger, so it is finite.
   */
  bool is_finite() const override { return !this->repeat; }

protected:
  /**
   * @brief Constructs the perpetual timer body.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  TimerBase(TimerFn f, bool repeat)
      : AnimationBase<Derived>(-1, repeat), f(std::move(f)) {}

  TimerFn f;         /**< The callback function. */
  uint32_t next = 0; /**< The target frame count for the next trigger. */
};

/**
 * @brief An animation that triggers a callback after a random delay.
 */
class RandomTimer : public TimerBase<RandomTimer> {
public:
  /**
   * @brief Constructs a RandomTimer.
   * @param min Minimum delay in frames.
   * @param max Maximum delay in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  RandomTimer(int min, int max, TimerFn f, bool repeat = false)
      : TimerBase(std::move(f), repeat), min(min), max(max) {
    HS_CHECK(min >= 0 && min <= max);
    HS_CHECK(max < std::numeric_limits<int>::max(),
             "RandomTimer max must be < INT_MAX (reset adds 1)");
    reset();
  }

  /**
   * @brief Calculates the next random trigger time.
   */
  HS_COLD_MEMBER void reset() {
    // +1 because hs::rand_int is half-open [min, max); the documented maximum
    // delay is inclusive, and RandomTimer(n, n) must yield exactly n.
    next = t + hs::rand_int(min, max + 1);
  }

private:
  int min; /**< Minimum frame delay. */
  int max; /**< Maximum frame delay. */
};

/**
 * @brief An animation that triggers a callback at regular intervals.
 */
class PeriodicTimer : public TimerBase<PeriodicTimer> {
public:
  /**
   * @brief Constructs a PeriodicTimer.
   * @param period The interval between calls, in frames; clamped to >= 1, so
   * period 0 fires on the next step.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  PeriodicTimer(int period, TimerFn f, bool repeat = false)
      : TimerBase(std::move(f), repeat), period(clamp_period(period)) {
    reset();
  }

  /**
   * @brief Calculates the next periodic trigger time.
   */
  void reset() { next = t + period; }

  /**
   * @brief Live-updates the trigger interval; reschedules the next trigger from
   * now when the clamped interval actually changes.
   * @param new_period New interval in frames; clamped to >= 1.
   * @details Clamps to >= 1: a 0/negative period makes `next = t + period <= t`,
   * which fires the callback every frame (and re-triggers on its own reset).
   * An unchanged period returns without rescheduling: a PeriodicTimer never
   * reaches done(), so a per-frame call that always reset() would defer the
   * callback forever with nothing to record it.
   */
  void set_period(int new_period) {
    int clamped = clamp_period(new_period);
    if (clamped == period)
      return;
    period = clamped;
    reset();
  }

private:
  static int clamp_period(int p) { return p < 1 ? 1 : p; }

  int period; /**< The interval in frames. */
};

} // namespace Animation
