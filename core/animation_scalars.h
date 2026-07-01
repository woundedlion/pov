/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "animation_core.h"

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
      : AnimationBase(-1, repeat), period(period < 1 ? 1 : period),
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
    period = (new_period < 1 ? 1 : new_period);
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

private:
  int period;    /**< The interval in frames. */
  TimerFn f;     /**< The callback function. */
  uint32_t next; /**< The target frame count for the next trigger. */
};

/**
 * @brief An animation that smoothly transitions a float variable over time.
 */
class Transition : public AnimationBase<Transition> {
public:
  /**
   * @brief Constructs a Transition animation.
   * @param mutant The float variable to modify.
   * @param to The target value.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to use.
   * @param quantized If true, the final value is rounded down to an integer.
   * @param repeat If true, the transition repeats indefinitely.
   */
  Transition(float &mutant, float to, int duration, EasingFn easing_fn,
             bool quantized = false, bool repeat = false)
      : AnimationBase(duration, repeat), mutant(mutant), from(mutant), to(to),
        easing_fn(std::move(easing_fn)), quantized(quantized) {
    // Reject the perpetual -1 the base permits: step() would clamp t_norm to 0
    // forever, silently freezing the tween instead of driving it.
    HS_CHECK(duration >= 0, "Transition duration must be >= 0");
  }

  /**
   * @brief Performs one step of the transition.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    // Snapshot the start once, on the first step. A repeating Transition rewinds
    // to t==0 with mutant == to, so re-snapshotting would freeze from = to.
    if (!captured) {
      from = mutant;
      captured = true;
    }
    AnimationBase::step(canvas);
    // Clamp both ends so a negative duration can't feed easing a negative arg.
    auto t_norm = hs::clamp(static_cast<float>(this->t) / duration, 0.0f, 1.0f);
    auto n = easing_fn(t_norm) * (to - from) + from;
    if (quantized) {
      n = std::floor(n);
    }
    mutant.get() = n;
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   * @details Rebinds only the target reference; the `from` snapshot and the
   * `captured` flag are deliberately retained. This is therefore safe for
   * *relocating* the binding to a variable holding the same value (e.g. a moved
   * params struct), but NOT for *retargeting* mid-flight to a variable with a
   * different start value — the in-flight interpolation would still run from the
   * original `from` toward `to`. Retargeting requires a fresh Transition.
   */
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

private:
  std::reference_wrapper<float>
      mutant;         /**< Reference to the float variable being animated. */
  float from;         /**< Starting value, captured on the first step. */
  float to;           /**< Target value. */
  EasingFn easing_fn; /**< Easing curve. */
  bool quantized;     /**< Flag to round result to integer. */
  bool captured = false; /**< True once `from` has been snapshotted. */
};

/**
 * @brief An animation that applies a custom function to a float variable over
 * time.
 */
class Mutation : public AnimationBase<Mutation> {
public:
  /**
   * @brief Constructs a Mutation animation.
   * @param mutant The float variable to modify.
   * @param f The custom mutation function: `float f(float eased_t)`.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to apply to the time factor.
   * @param repeat If true, the mutation repeats indefinitely.
   * @param paused Optional pause gate; null = always runs.
   */
  Mutation(float &mutant, ScalarFn f, int duration, EasingFn easing_fn,
           bool repeat = false, const bool *paused = nullptr)
      : AnimationBase(duration, repeat), mutant(mutant), f(std::move(f)),
        easing_fn(std::move(easing_fn)), paused_(paused) {
    // Reject the perpetual -1 the base permits: step() would clamp t_norm to 0
    // forever, silently freezing the curve instead of driving it.
    HS_CHECK(duration >= 0, "Mutation duration must be >= 0");
  }

  /**
   * @brief Performs one step of the mutation.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details When wired to a pause flag (an effect's `anims_paused_`), the
   * mutation freezes while paused: it neither advances its timer nor writes the
   * mutant, so a GUI slider bound to the same member is the sole writer and the
   * user's value holds. Resuming hands the member back to the curve.
   */
  void step(Canvas &canvas) override {
    if (is_paused(paused_))
      return;
    AnimationBase::step(canvas);
    auto t_norm = hs::clamp(static_cast<float>(this->t) / duration, 0.0f, 1.0f);
    mutant.get() = f(easing_fn(t_norm));
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   */
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

private:
  std::reference_wrapper<float>
      mutant; /**< Reference to the float variable being modified. */
  ScalarFn f; /**< The custom function to apply. */
  EasingFn easing_fn; /**< Easing curve. */
  const bool *paused_; /**< Optional pause gate; freezes the mutation when set
                          and true. Null = always runs. */
};

/**
 * @brief An animation that continuously increments a float variable over time.
 *
 * A Driver is perpetual with a one-frame cycle (duration 1, repeat true): it has
 * no natural completion, so a `.then()` callback attached to it fires once PER
 * FRAME (after each step + rewind). That is intentional — for a per-frame
 * incrementer, a per-frame hook is the only meaningful interpretation of
 * `.then()`. Do not attach a one-shot callback expecting a single fire.
 */
class Driver : public AnimationBase<Driver> {
public:
  /**
   * @brief Constructs a Driver animation for continuous progression.
   * @param mutant The float variable to modify.
   * @param speed The amount to add per frame.
   * @param wrap If true, wraps the value back into the [0,1) range each step.
   * @param paused Optional pause gate; null = always runs.
   */
  Driver(float &mutant, float speed, bool wrap = true,
         const bool *paused = nullptr)
      : AnimationBase(1, true), mutant(mutant), speed(speed), wrap_(wrap),
        paused_(paused) {
    // A non-finite speed permanently poisons `mutant` (wrap_t can't recover NaN).
    HS_CHECK(std::isfinite(speed), "Driver: fixed speed must be finite");
  }

  /**
   * @brief Constructs a Driver whose speed tracks a live param every step.
   * @param mutant The float variable to modify.
   * @param speed_src Pointer to a live value (e.g. a registered slider); the
   *   effective speed is `*speed_src * scale`, re-read each step.
   * @param scale Per-unit multiplier applied to *speed_src.
   * @param wrap If true, wraps the value back into the [0,1) range each step.
   * @param paused Optional pause gate; null = always runs.
   * @details The Driver pulls the slider itself each step, so callers can
   *   `timeline.add(...)` and drop the handle rather than re-syncing the speed
   *   via set_speed() every frame.
   */
  // Borrow contract: speed_src is stored and dereferenced every step, so the
  // referent must outlive the Driver (e.g. an effect's registered param). A
  // pointer parameter already rejects a temporary at compile time, so no deleted
  // overload is needed.
  Driver(float &mutant, const float *speed_src, float scale, bool wrap = true,
         const bool *paused = nullptr)
      : AnimationBase(1, true), mutant(mutant), speed(0.0f), wrap_(wrap),
        paused_(paused), speed_src_(speed_src), scale_(scale) {
    HS_CHECK(speed_src != nullptr, "Driver: live speed_src is null");
    // On a non-finite initial read keep the 0.0f seed (see step()).
    float s = *speed_src * scale;
    if (std::isfinite(s))
      speed = s;
  }

  /**
   * @brief Performs one step by adding the speed to the mutant.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details Freezes while a wired pause flag is set (see Mutation::step).
   */
  void step(Canvas &canvas) override {
    if (is_paused(paused_))
      return;
    AnimationBase::step(canvas);
    // Re-read the live source, keeping the last good speed on a non-finite read:
    // a one-frame NaN/Inf would otherwise poison `mutant` permanently.
    if (speed_src_) {
      float s = *speed_src_ * scale_;
      if (std::isfinite(s))
        speed = s;
    }
    mutant.get() += speed;
    if (wrap_) {
      mutant.get() = wrap_t(mutant.get());
    }
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   */
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

  /**
   * @brief Gets the current speed.
   * @return The amount added to the mutant per frame.
   */
  float get_speed() const { return speed; }

  /**
   * @brief Sets the current speed.
   * @param new_speed The new amount to add to the mutant per frame.
   * @pre No live speed_src is bound; step() would otherwise clobber the manual
   *   speed each frame from the source.
   */
  void set_speed(float new_speed) {
    HS_CHECK(!speed_src_, "Driver::set_speed with a live speed_src bound");
    speed = new_speed;
  }

  /**
   * @brief Gets the current mutant value.
   * @return Const reference to the driven float variable.
   */
  const float &get_mutant() const { return mutant.get(); }

private:
  std::reference_wrapper<float> mutant; /**< Reference to the float variable. */
  float speed;                          /**< Amount added per frame. */
  bool wrap_;                           /**< If true, wraps value to 0-1 range. */
  const bool *paused_; /**< Optional pause gate; freezes the driver when set and
                          true. Null = always runs. */
  const float *speed_src_ =
      nullptr;            /**< Live speed source; null = fixed speed. */
  float scale_ = 1.0f;   /**< Multiplier applied to *speed_src_. */
};

/**
 * @brief An animation that interpolates between states.
 * @details The caller owns the start, subject, and target data. Lerp just holds
 * pointers and a type-erased lerp function. Supports any type T that implements
 * lerp(start, target, t).
 */
class Lerp : public AnimationBase<Lerp> {
public:
  /**
   * @brief Constructs a Lerp animation.
   * @tparam T Interpolated type implementing lerp(start, target, t).
   * @tparam Easing Easing-function type.
   * @param subject The caller-owned object to write each frame.
   * @param start The caller-owned start state.
   * @param target The caller-owned target state.
   * @param duration The duration in frames.
   * @param easing_fn The easing function applied to progress.
   * @param paused Optional pause gate; null = always runs.
   */
  template <typename T, typename Easing>
  Lerp(T &subject, const T &start, const T &target, int duration,
       Easing easing_fn, const bool *paused = nullptr)
      : AnimationBase(duration, false), subject_ptr(&subject),
        start_ptr(&start), target_ptr(&target), easing(easing_fn),
        paused_(paused) {
    do_lerp = [](void *subj, const void *s, const void *tgt, float t) {
      static_cast<T *>(subj)->lerp(*static_cast<const T *>(s),
                                   *static_cast<const T *>(tgt), t);
    };
  }

  // Borrow contract: subject/start/target are stored as pointers and read across
  // many frames, so they must outlive the Timeline. These deleted overloads
  // reject a temporary start/target at compile time; the defaulted `paused`
  // mirrors the real ctor so the 6-arg paused call is rejected too.
  template <typename T, typename Easing>
  Lerp(T &subject, const T &&start, const T &target, int duration,
       Easing easing_fn, const bool *paused = nullptr) = delete;
  template <typename T, typename Easing>
  Lerp(T &subject, const T &start, const T &&target, int duration,
       Easing easing_fn, const bool *paused = nullptr) = delete;
  template <typename T, typename Easing>
  Lerp(T &subject, const T &&start, const T &&target, int duration,
       Easing easing_fn, const bool *paused = nullptr) = delete;

  /**
   * @brief Performs one step of the interpolation.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details Freezes while a wired pause flag is set (see Mutation::step).
   * Preset-cycling effects pass an Effect's `anims_paused_` so a paused lerp
   * neither advances nor writes its subject, and (for `.then`-chained lerps)
   * never completes — so the whole preset chain halts and the slider-bound
   * subject holds.
   */
  void step(Canvas &canvas) override {
    if (is_paused(paused_))
      return;
    AnimationBase::step(canvas);
    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    do_lerp(subject_ptr, start_ptr, target_ptr, easing(progress));
  }

private:
  void *subject_ptr;       /**< Pointer to the caller-owned subject. */
  const void *start_ptr;   /**< Pointer to the caller-owned start state. */
  const void *target_ptr;  /**< Pointer to the caller-owned target state. */
  EasingFn easing;         /**< Easing curve applied to progress. */
  /** @brief Type-erased lerp thunk. */
  void (*do_lerp)(void *, const void *, const void *, float);
  const bool *paused_ = nullptr; /**< Optional pause gate; null = always runs. */
};

/**
 * @brief An animation that draws a sprite while managing its fade-in/out
 * effects.
 * @details Computes opacity inline rather than embedding Transition objects.
 * An indefinite sprite (duration -1) never completes, so a `.then()` callback
 * attached to it never fires.
 */
class Sprite : public AnimationBase<Sprite> {
public:
  /**
   * @brief Constructs a Sprite animation.
   * @param draw_fn The function to call each frame for drawing (`void
   * draw_fn(Canvas&, float opacity)`).
   * @param duration The duration the sprite is fully visible (-1 for
   * indefinite).
   * @param fade_in_duration Frames for fading in.
   * @param fade_in_easing_fn Easing for fade-in.
   * @param fade_out_duration Frames for fading out.
   * @param fade_out_easing_fn Easing for fade-out.
   * @param paused Optional pause gate; null = always runs.
   */
  Sprite(SpriteFn draw_fn, int duration, int fade_in_duration = 0,
         EasingFn fade_in_easing_fn = ease_linear, int fade_out_duration = 0,
         EasingFn fade_out_easing_fn = ease_linear, const bool *paused = nullptr)
      : AnimationBase(duration, false), draw_fn(std::move(draw_fn)),
        fade_in_duration(fade_in_duration),
        fade_out_duration(fade_out_duration),
        fade_in_easing(std::move(fade_in_easing_fn)),
        fade_out_easing(std::move(fade_out_easing_fn)), paused_(paused) {
    HS_CHECK(fade_in_duration >= 0, "Sprite fade-in duration must be >= 0");
    HS_CHECK(fade_out_duration >= 0, "Sprite fade-out duration must be >= 0");
    // Overlapping windows (durations are independent GUI sliders): scale both
    // fades proportionally to fit the visible duration, so the envelope still
    // peaks at full opacity and stays a continuous triangle (definite sprites
    // only; indefinite ones skip fade-out).
    int fade_total = this->fade_in_duration + this->fade_out_duration;
    if (duration >= 0 && fade_total > duration) {
      this->fade_in_duration =
          static_cast<long long>(duration) * this->fade_in_duration / fade_total;
      this->fade_out_duration = duration - this->fade_in_duration;
    }
  }

  /**
   * @brief Updates the drawing function used by the sprite.
   * @param new_draw_fn The new per-frame drawing functor.
   */
  void rebind_draw(SpriteFn new_draw_fn) { draw_fn = std::move(new_draw_fn); }

  /**
   * @brief Steps the animation, computes the current opacity inline, and calls
   * the draw function.
   * @param canvas The canvas buffer passed to the draw function.
   */
  void step(Canvas &canvas) override {
    // Paused: hold the frame (don't advance the timer) but keep drawing at the
    // current opacity.
    if (!is_paused(paused_))
      AnimationBase::step(canvas);

    // Trapezoid envelope as the MIN of an independent fade-in and fade-out ramp.
    // Computing both keeps opacity continuous when the windows overlap (the
    // durations are independent GUI sliders), degrading to a triangle.
    float fade_in = 1.0f;
    if (fade_in_duration > 0 && t < static_cast<uint32_t>(fade_in_duration)) {
      float progress = static_cast<float>(t) / fade_in_duration;
      fade_in = fade_in_easing(hs::clamp(progress, 0.0f, 1.0f));
    }

    // An indefinite sprite (duration < 0) has no end frame to fade toward, so
    // the `duration >= 0` guard skips fade-out for it (load-bearing, not just
    // overflow safety).
    float fade_out = 1.0f;
    if (duration >= 0 && fade_out_duration > 0 &&
        t + static_cast<uint32_t>(fade_out_duration) >=
            static_cast<uint32_t>(duration)) {
      float elapsed = static_cast<float>(t + static_cast<uint32_t>(fade_out_duration) -
                                         static_cast<uint32_t>(duration));
      float progress = elapsed / fade_out_duration;
      fade_out = 1.0f - fade_out_easing(hs::clamp(progress, 0.0f, 1.0f));
    }

    draw_fn(canvas, std::min(fade_in, fade_out));
  }

private:
  SpriteFn draw_fn;         /**< The drawing function functor. */
  int fade_in_duration;     /**< Duration of fade-in phase in frames. */
  int fade_out_duration;    /**< Duration of fade-out phase in frames. */
  EasingFn fade_in_easing;  /**< Easing curve for fade-in. */
  EasingFn fade_out_easing; /**< Easing curve for fade-out. */
  const bool *paused_ = nullptr; /**< Optional pause gate; holds the frame (no
                                    timer advance, still draws) when set. */
};

} // namespace Animation
