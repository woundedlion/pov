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
 * as an active-frame countdown.
 */
template <typename Derived> class TimerBase : public AnimationBase<Derived> {
public:
  /**
   * @brief Steps the timer, calling the function if the delay has elapsed.
   * @param canvas The canvas buffer (forwarded to the timer callback).
   */
  void step(Canvas &canvas) override {
    AnimationBase<Derived>::step(canvas);
    if (remaining_delay > 1) {
      --remaining_delay;
      return;
    }
    remaining_delay = 0;
    f(canvas);
    if (this->repeats()) {
      static_cast<Derived *>(this)->reset();
      // A repeating timer never reaches done(), so fire the per-cycle .then()
      // directly to honor the contract.
      this->post_callback();
    } else {
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

  TimerFn f;                    /**< The callback function. */
  uint32_t remaining_delay = 0; /**< Active frames until the next trigger. */
};

/** @brief Delay range and repeat behaviour of a RandomTimer. */
struct RandomTimerOptions {
  /** @brief Minimum sampled delay in frames; zero fires on the next step, like one. */
  int min = 0;
  /** @brief Maximum delay in frames, inclusive. */
  int max = 0;
  /** @brief Resets the timer after each call instead of ending it. */
  bool repeat = false;
};

/**
 * @brief An animation that triggers a callback after a random delay.
 */
class RandomTimer : public TimerBase<RandomTimer> {
public:
  using Options = RandomTimerOptions;

  /**
   * @brief Constructs a RandomTimer.
   * @param options The delay range and repeat behaviour.
   * @param f The function to call when the timer elapses.
   */
  RandomTimer(const Options &options, TimerFn f)
      : TimerBase(std::move(f), options.repeat), min(options.min),
        max(options.max) {
    HS_CHECK(min >= 0 && min <= max, "RandomTimer: invalid frame range");
    HS_CHECK(max < std::numeric_limits<int>::max(),
             "RandomTimer max must be < INT_MAX (reset adds 1)");
    reset();
  }

  /**
   * @brief Schedules the next random active-frame delay.
   */
  HS_COLD_MEMBER void reset() {
    // +1 because hs::rand_int is half-open [min, max); the documented maximum
    // delay is inclusive. A sampled zero fires on the next step, like one.
    remaining_delay = hs::rand_int(min, max + 1);
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
   * @brief Schedules the next periodic active-frame delay.
   */
  void reset() { remaining_delay = period; }

  /**
   * @brief Live-updates the trigger interval; reschedules the next trigger from
   * now when the clamped interval actually changes.
   * @param new_period New interval in frames; clamped to >= 1.
   * @details An unchanged period does not reschedule, so calling this every
   * frame cannot defer the callback forever.
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
