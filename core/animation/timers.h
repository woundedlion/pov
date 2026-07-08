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
 * @brief An animation that triggers a callback after a random delay.
 */
class RandomTimer : public AnimationBase<RandomTimer> {
public:
  /**
   * @brief Constructs a RandomTimer.
   * @param min Minimum delay in frames.
   * @param max Maximum delay in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  RandomTimer(int min, int max, TimerFn f, bool repeat = false)
      : AnimationBase(-1, repeat), min(min), max(max), f(std::move(f)),
        next(0) {
    HS_CHECK(min >= 0 && min <= max);
    HS_CHECK(max < std::numeric_limits<int>::max(),
             "RandomTimer max must be < INT_MAX (reset adds 1)");
    reset();
  }

  /**
   * @brief Calculates the next random trigger time.
   */
  void reset() {
    // +1 because hs::rand_int is half-open [min, max); the documented maximum
    // delay is inclusive, and RandomTimer(n, n) must yield exactly n.
    next = t + hs::rand_int(min, max + 1);
  }

  /**
   * @brief Steps the timer, calling the function if the delay has elapsed.
   * @param canvas The canvas buffer (forwarded to the timer callback).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
        // A repeating timer never reaches done(), so fire the per-cycle .then()
        // directly to honor the contract.
        this->post_callback();
      } else {
        cancel();
      }
    }
  }

  /** @brief A one-shot timer is removed on its single trigger, so it is finite. */
  bool is_finite() const override { return !repeat; }

private:
  int min;   /**< Minimum frame delay. */
  int max;   /**< Maximum frame delay. */
  TimerFn f;      /**< The callback function. */
  uint32_t next;  /**< The target frame count for the next trigger. */
};

/**
 * @brief An animation that triggers a callback at regular intervals.
 */
class PeriodicTimer : public AnimationBase<PeriodicTimer> {
public:
  /**
   * @brief Constructs a PeriodicTimer.
   * @param period The interval between calls, in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  PeriodicTimer(int period, TimerFn f, bool repeat = false)
      : AnimationBase(-1, repeat), period(clamp_period(period)),
        f(std::move(f)) {
    reset();
  }

  /**
   * @brief Calculates the next periodic trigger time.
   */
  void reset() { next = t + period; }

  /**
   * @brief Live-updates the trigger interval; reschedules the next trigger from
   * now.
   * @param new_period New interval in frames; clamped to >= 1.
   * @details Clamps to >= 1: a 0/negative period makes `next = t + period <= t`,
   * which fires the callback every frame (and re-triggers on its own reset).
   */
  void set_period(int new_period) {
    period = clamp_period(new_period);
    reset();
  }

  /**
   * @brief Steps the timer, calling the function if the period has elapsed.
   * @param canvas The canvas buffer (forwarded to the timer callback).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
        this->post_callback(); // see RandomTimer::step
      } else {
        cancel();
      }
    }
  }

  /** @brief A one-shot timer is removed on its single trigger, so it is finite. */
  bool is_finite() const override { return !repeat; }

private:
  static int clamp_period(int p) { return p < 1 ? 1 : p; }

  int period;    /**< The interval in frames. */
  TimerFn f;     /**< The callback function. */
  uint32_t next; /**< The target frame count for the next trigger. */
};

} // namespace Animation
