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
      : AnimationBase(duration, repeat), mutant(mutant), from(0.0f), to(to),
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
   * @details Retains the `from` snapshot and `captured` flag: safe to relocate
   * the binding to a same-value variable, but not to retarget mid-flight (the
   * tween still runs from the original `from`). Retargeting needs a fresh
   * Transition.
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
        easing_fn(std::move(easing_fn)), paused(paused) {
    // Reject the perpetual -1 the base permits: step() would clamp t_norm to 0
    // forever, silently freezing the curve instead of driving it.
    HS_CHECK(duration >= 0, "Mutation duration must be >= 0");
  }

  /**
   * @brief Performs one step of the mutation.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details While the wired pause flag is set, freezes: neither advances the
   * timer nor writes the mutant, so a GUI slider bound to the same member holds.
   */
  void step(Canvas &canvas) override {
    if (is_paused(paused))
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
  const bool *paused; /**< Optional pause gate; freezes the mutation when set
                          and true. Null = always runs. */
};

/**
 * @brief An animation that continuously increments a float variable over time.
 *
 * Perpetual, one-frame cycle (duration 1, repeat true): a `.then()` callback
 * fires once PER FRAME. Do not attach a one-shot callback expecting a single
 * fire.
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
      : AnimationBase(1, true), mutant(mutant), speed(speed), wrap(wrap),
        paused(paused) {
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
   * @details Re-reads speed_src each step, so callers need not set_speed() per
   *   frame.
   */
  // Borrow contract: speed_src is dereferenced every step, so the referent must
  // outlive the Driver (e.g. an effect's registered param).
  Driver(float &mutant, const float *speed_src, float scale, bool wrap = true,
         const bool *paused = nullptr)
      : AnimationBase(1, true), mutant(mutant), speed(0.0f), wrap(wrap),
        paused(paused), speed_src(speed_src), scale(scale) {
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
    if (is_paused(paused))
      return;
    AnimationBase::step(canvas);
    // Re-read the live source, keeping the last good speed on a non-finite read:
    // a one-frame NaN/Inf would otherwise poison `mutant` permanently.
    if (speed_src) {
      float s = *speed_src * scale;
      if (std::isfinite(s))
        speed = s;
    }
    mutant.get() += speed;
    if (wrap) {
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
    HS_CHECK(!speed_src, "Driver::set_speed with a live speed_src bound");
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
  bool wrap;                            /**< If true, wraps value to 0-1 range. */
  const bool *paused; /**< Optional pause gate; freezes the driver when set and
                          true. Null = always runs. */
  const float *speed_src =
      nullptr;        /**< Live speed source; null = fixed speed. */
  float scale = 1.0f; /**< Multiplier applied to *speed_src. */
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
        paused(paused) {
    // Reject the perpetual -1 the base permits: step() would clamp t_norm to 0
    // forever, silently freezing the lerp instead of driving it.
    HS_CHECK(duration >= 0, "Lerp duration must be >= 0");
    do_lerp = [](void *subj, const void *s, const void *tgt, float t) {
      static_cast<T *>(subj)->lerp(*static_cast<const T *>(s),
                                   *static_cast<const T *>(tgt), t);
    };
  }

  // Borrow contract: subject/start/target are read across many frames, so they
  // must outlive the Timeline; these deleted overloads reject a temporary.
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
   * @details Freezes while a wired pause flag is set (see Mutation::step). A
   * paused, `.then`-chained lerp never completes, so the whole chain halts.
   */
  void step(Canvas &canvas) override {
    if (is_paused(paused))
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
  const bool *paused = nullptr; /**< Optional pause gate; null = always runs. */
};

/**
 * @brief An animation that smoothly interpolates a GenerativePalette toward a
 * target palette.
 */
class ColorWipe : public AnimationBase<ColorWipe> {
public:
  /**
   * @brief Constructs a ColorWipe animation.
   * @param from_palette The GenerativePalette to animate (snapshotted on the
   * first step(), mirroring Transition).
   * @param to_palette The GenerativePalette to interpolate toward.
   * @param duration The duration in frames.
   * @param easing_fn The easing function.
   */
  ColorWipe(GenerativePalette &from_palette,
            const GenerativePalette &to_palette, int duration,
            EasingFn easing_fn)
      : AnimationBase(duration, false), cur_palette(from_palette),
        to_snap(to_palette.snapshot()), easing_fn(std::move(easing_fn)) {
    // Reject the perpetual -1 the base permits: step() would clamp amount to 0
    // forever, silently freezing the palette instead of driving it.
    HS_CHECK(duration >= 0, "ColorWipe duration must be >= 0");
  }

  /**
   * @brief Steps the animation, blending the palette's colors based on the time
   * factor.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    // Snapshot the start palette once, on the first step (see Transition::step).
    if (!captured) {
      from_snap = cur_palette.get().snapshot();
      captured = true;
    }
    AnimationBase::step(canvas);
    float amount = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    cur_palette.get().lerp(from_snap, to_snap, easing_fn(amount));
  }

private:
  GenerativePalette::Snapshot from_snap{}; /**< Snapshot of starting colors. */
  std::reference_wrapper<GenerativePalette>
      cur_palette;                     /**< The palette being animated. */
  GenerativePalette::Snapshot to_snap; /**< Snapshot of target colors. */
  EasingFn easing_fn;                  /**< Easing curve. */
  bool captured = false; /**< Whether from_snap was taken on the first step. */
};

/**
 * @brief Animates the Mobius parameters for a continuous loxodromic flow.
 */
class MobiusFlow : public AnimationBase<MobiusFlow> {
public:
  /**
   * @brief Constructs a MobiusFlow animation.
   * @param params Reference to the MobiusParams to animate.
   * @param num_rings Reference to the number of rings (scalar).
   * @param num_lines Reference to the number of lines (scalar).
   * @param duration Duration of the flow.
   * @param repeat Whether to repeat.
   */
  MobiusFlow(MobiusParams &params, const float &num_rings,
             const float &num_lines, int duration, bool repeat = true)
      : AnimationBase(duration, repeat), params(params), num_rings(num_rings),
        num_lines(num_lines) {
    HS_CHECK(duration >= 0, "MobiusFlow duration must be >= 0");
  }

  // Borrow contract: num_rings/num_lines are read every frame, so they must
  // outlive the Timeline; these deleted overloads reject a temporary scalar.
  MobiusFlow(MobiusParams &params, const float &&num_rings,
             const float &num_lines, int duration, bool repeat = true) = delete;
  MobiusFlow(MobiusParams &params, const float &num_rings,
             const float &&num_lines, int duration, bool repeat = true) = delete;
  MobiusFlow(MobiusParams &params, const float &&num_rings,
             const float &&num_lines, int duration, bool repeat = true) = delete;

  /**
   * @brief Steps the animation, updating params a and d.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    // Floor rings at 0 so the divisor rings + 1 stays >= 1; a non-finite or -1
    // rings would make logPeriod inf and poison a/d.
    float rings = num_rings;
    if (!std::isfinite(rings) || rings < 0.0f)
      rings = 0.0f;
    float logPeriod = 5.0f / (rings + 1);
    float flowParam = progress * logPeriod;
    float scale = expf(flowParam);
    float s = sqrtf(scale);
    // Clamp lines to >= 1 before dividing (the slider bottoms out at 0 → 2π/0);
    // a non-finite lines falls back to 1.
    float lines = num_lines;
    if (!std::isfinite(lines) || lines < 1.0f)
      lines = 1.0f;
    float angle = progress * (2 * PI_F / lines);

    params.get().a.re = s * cosf(angle);
    params.get().a.im = s * sinf(angle);
    params.get().d.re = (1.0f / s) * cosf(-angle);
    params.get().d.im = (1.0f / s) * sinf(-angle);
  }

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  /** @brief Live ring count driving the log period. */
  std::reference_wrapper<const float> num_rings;
  /** @brief Live line count driving the angular step. */
  std::reference_wrapper<const float> num_lines;
};

/**
 * @brief Animates the Mobius parameters for a warping effect pulling the poles
 * together.
 */
class MobiusWarp : public AnimationBase<MobiusWarp> {
public:
  /**
   * @brief Constructs a MobiusWarp animation.
   * @param params Reference to the MobiusParams to animate.
   * @param scale The magnitude of the warp effect.
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   * @param easing The easing function to use (default: ease_in_out_sin).
   */
  MobiusWarp(MobiusParams &params, float scale, int duration,
             bool repeat = true, EasingFn easing = ease_in_out_sin)
      : AnimationBase(duration, repeat), params(params), scale(scale),
        easing(easing) {
    HS_CHECK(duration >= 0, "MobiusWarp duration must be >= 0");
    HS_CHECK(std::isfinite(scale), "MobiusWarp scale must be finite");
  }

  /**
   * @brief Binds the warp magnitude to a live external float.
   * @param live_scale The external float to read each frame as the magnitude.
   * @details Binding makes step() read the referent every frame, so a wired GUI
   * slider takes effect immediately. The referent must outlive the animation.
   */
  void bind_scale(const float &live_scale) { scale_ref = &live_scale; }
  void bind_scale(const float &&) = delete;

  /**
   * @brief Steps the animation, updating param b.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    float s = scale;
    if (scale_ref) {
      float s2 = *scale_ref;
      if (std::isfinite(s2))
        s = s2;
    }
    params.get().b.re = s * (cosf(angle) - 1.0f);
    params.get().b.im = s * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  float scale;                      /**< Warp magnitude. */
  EasingFn easing;                  /**< Easing curve. */
  const float *scale_ref = nullptr; /**< Optional live magnitude source. */
};

/**
 * @brief Animates the Mobius parameters for a circular warping effect.
 */
class MobiusWarpCircular : public AnimationBase<MobiusWarpCircular> {
public:
  /**
   * @brief Constructs a MobiusWarpCircular animation.
   * @param params Reference to the MobiusParams to animate.
   * @param scale The magnitude of the warp effect.
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   * @param easing The easing function to use (default: ease_in_out_sin).
   */
  MobiusWarpCircular(MobiusParams &params, float scale, int duration,
                     bool repeat = true, EasingFn easing = ease_in_out_sin)
      : AnimationBase(duration, repeat), params(params), scale(scale),
        easing(easing) {
    HS_CHECK(duration >= 0, "MobiusWarpCircular duration must be >= 0");
    HS_CHECK(std::isfinite(scale), "MobiusWarpCircular scale must be finite");
  }

  /**
   * @brief Steps the animation, updating param b.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    params.get().b.re = scale * cosf(angle);
    params.get().b.im = -scale * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  float scale;     /**< Warp magnitude. */
  EasingFn easing; /**< Easing curve. */
};

/**
 * @brief Continuously modulates Mobius parameters to create an evolving warp.
 * @details Uses multiple frequencies for non-repeating chaos.
 *
 * PERPETUAL (duration -1, no repeat): never reaches done(), so a `.then()`
 * callback NEVER fires. Drive follow-on behavior from a finite animation or
 * cancel() it explicitly.
 */
class MobiusWarpEvolving : public AnimationBase<MobiusWarpEvolving> {
public:
  /**
   * @brief Constructs a MobiusWarpEvolving animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   * @note `base` snapshots `params` at construction and is latched at spawn;
   * live edits require a respawn (live `scale`/`speed` go through the setters).
   */
  MobiusWarpEvolving(MobiusParams &params, float scale = 0.5f,
                     float speed = 0.01f)
      : params(params), speed(speed), scale(scale), base(params),
        seed(hs::random()()) {
    HS_CHECK(std::isfinite(scale) && std::isfinite(speed),
             "MobiusWarpEvolving scale and speed must be finite");
  }

  /** @brief Sets the modulation speed (radians of phase per frame unit). */
  void set_speed(float speed) {
    if (std::isfinite(speed))
      this->speed = speed;
  }

  /** @brief Sets the per-channel modulation magnitude. */
  void set_scale(float scale) {
    if (std::isfinite(scale))
      this->scale = scale;
  }

  /**
   * @brief Derives a per-channel phase offset from the seed and channel index.
   * @param i Channel index.
   * @return Phase offset in [0, 100) for that channel.
   */
  float phase(int i) const {
    uint32_t h = seed ^ (static_cast<uint32_t>(i) * 2654435761u);
    return (h & 0xFFFF) * (100.0f / 65536.0f);
  }

  /**
   * @brief Steps the animation, modulating all eight Mobius params.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    // Accepted limit: past t == 2^24 (~77 h at 60 fps, sooner at higher speed)
    // float can't represent consecutive frames and the phase freezes; this
    // animation is perpetual (duration == -1).
    float time = t * speed;
    float s = scale;

    // Use prime-ish number ratios for frequencies to minimize repetition cycle
    params.get().a.re = base.a.re + sinf(time * 1.0f + phase(0)) * s;
    params.get().a.im = base.a.im + cosf(time * 1.13f + phase(1)) * s;

    params.get().b.re = base.b.re + sinf(time * 1.27f + phase(2)) * s;
    params.get().b.im = base.b.im + cosf(time * 1.39f + phase(3)) * s;

    params.get().c.re = base.c.re + sinf(time * 0.71f + phase(4)) * s;
    params.get().c.im = base.c.im + cosf(time * 0.83f + phase(5)) * s;

    params.get().d.re = base.d.re + sinf(time * 0.97f + phase(6)) * s;
    params.get().d.im = base.d.im + cosf(time * 1.09f + phase(7)) * s;
  }

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  float speed; /**< Animation speed (radians of phase per frame unit). */
  float scale; /**< Magnitude of the per-channel modulation. */
  MobiusParams base; /**< Baseline params captured at construction. */
  uint32_t seed;     /**< Seed for the per-channel phase offsets. */
};

/**
 * @brief Parameters for a ripple wave effect.
 */
struct RippleParams {
  Vector center;          /**< Center point of the ripple source. */
  float amplitude = 0.0f; /**< Current height of the wave. */
  float phase = 0.0f;     /**< Current phase offset (time). */
  float decay{5.0};      /**< Spatial decay rate. */
  float thickness{1.0f}; /**< Thickness of the ripple. */

  /** @brief Cached cos(angle) lower fast-reject bound. */
  float cos_threshold_min = 1.0f;
  /** @brief Cached cos(angle) upper fast-reject bound. */
  float cos_threshold_max = -1.0f;

  /**
   * @brief Ricker wavelet half-width, floored so the distance normalization
   * never divides by zero.
   */
  float half_width() const {
    float hw = thickness * 0.5f;
    return hw < 0.001f ? 0.001f : hw;
  }

  /**
   * @brief Recomputes the cos(angle) fast-reject bounds for the wavelet's
   * current phase and thickness, so the renderer can skip points outside the
   * ring.
   */
  void prepare_thresholds() {
    float hw = half_width();
    // Clamp into [0,π]: the active ring lies within the sphere's angular range,
    // so cos(clamped) keeps the fast-reject band engaged past phase=π instead of
    // collapsing both bounds to accept-all. cos(0)=1 and cos(π)=-1 reproduce the
    // out-of-range sentinels at the endpoints.
    float d_min = hs::clamp(phase - hw * 2.0f, 0.0f, PI_F);
    float d_max = hs::clamp(phase + hw * 2.0f, 0.0f, PI_F);
    cos_threshold_min = cosf(d_min);
    cos_threshold_max = cosf(d_max);
  }
};

/**
 * @brief Animates a single ripple event: expanding outward and fading away.
 */
class Ripple : public AnimationBase<Ripple> {
public:
  /**
   * @brief Constructs a Ripple animation.
   * @param params Reference to the params struct to animate. `params.amplitude`
   *        is captured here as the ripple's peak and then reset to 0; set it
   *        before constructing, as later writes are ignored.
   * @param center The center point of the ripple.
   * @param speed How fast the waves travel.
   * @param duration How long the ripple lasts in frames.
   */
  Ripple(RippleParams &params, const Vector &center, float speed = 0.2f,
         int duration = 100)
      : AnimationBase(duration, false), params(params), speed(speed),
        peak_amplitude(params.amplitude) {
    HS_CHECK(duration >= 2, "Ripple duration must be >= 2");
    HS_CHECK(std::isfinite(speed), "Ripple speed must be finite");
    this->params.get().center = center;
    this->params.get().phase = 0.0f;
    // Start at 0 to prevent 1-frame singularity before first step()
    this->params.get().amplitude = 0.0f;
  }

  /**
   * @brief Steps the ripple: advances the wave and updates its envelope.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);

    float progress = static_cast<float>(t) / duration;
    float envelope = 0.0f;

    // Past the duration the envelope is pinned to 0, so skip the wave work; the
    // amplitude is still published below so the last frame renders nothing.
    if (t < static_cast<uint32_t>(duration)) {
      params.get().phase += speed;

      // Attack: fast ramp over ~10% of duration, squared for a parabolic ease-in.
      float attack_dur = 0.1f;
      float attack = std::min(progress / attack_dur, 1.0f);
      attack = attack * attack;

      float decay = 1.0f - progress;

      envelope = attack * decay;

      // Re-prepare the reject thresholds against the phase just advanced, so the
      // render never tests the new phase against thresholds cached at the old one.
      params.get().prepare_thresholds();
    }

    params.get().amplitude = peak_amplitude * envelope;
  }

private:
  std::reference_wrapper<RippleParams> params; /**< Ripple params to animate. */
  float speed;          /**< Wave travel speed (phase increment per frame). */
  float peak_amplitude; /**< Peak wave height captured at construction. */
};

/**
 * @brief Parameters for noise transformation.
 */
struct NoiseParams {
  float amplitude = 0.5f;   /**< Noise output amplitude. */
  float speed = 1.0f;       /**< Temporal evolution speed. */
  float frequency = 0.125f; /**< Spatial frequency of the noise. */
  float time = 0.0f;        /**< Current animation time. */
  float scale = 4.0f;       /**< Spatial scale factor. */
  mutable FastNoiseLite noise; /**< Backing generator; mutable for lazy
                                  init/updates. */

  /**
   * @brief Constructs noise params with an OpenSimplex2 generator.
   */
  NoiseParams() { noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2); }

  /**
   * @brief Syncs the generator's frequency with the `frequency` field.
   */
  void sync() const { noise.SetFrequency(frequency); }

  /**
   * @brief Refreshes live-tunable config from a template snapshot.
   * @param t Template params carrying the current slider values.
   * @details Copies the slider-driven fields but not the `time` axis or the
   * backing generator, so a live edit reaches a spawned entity without resetting
   * its phase. prepare_frame() invokes this before sync().
   */
  void refresh_from(const NoiseParams &t) {
    amplitude = t.amplitude;
    speed = t.speed;
    frequency = t.frequency;
    scale = t.scale;
  }
};

/**
 * @brief Animates noise parameters by updating time.
 */
class Noise : public AnimationBase<Noise> {
public:
  /**
   * @brief Constructs a Noise animation.
   * @param params Reference to the NoiseParams to animate.
   * @param duration Duration in frames (-1 for indefinite). Use -1: this
   *   repeats, so a finite duration rewinds t each cycle and snaps params.time
   *   backward.
   */
  Noise(NoiseParams &params, int duration = -1)
      : AnimationBase(duration, true), params(params) {}

  /**
   * @brief Steps the animation, advancing the noise time field.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    // Accepted limit: past t == 2^24 (~77 h at 60 fps, sooner at higher speed)
    // float can't represent consecutive frames and the noise time axis freezes.
    params.get().time = static_cast<float>(t);
  }

private:
  std::reference_wrapper<NoiseParams> params; /**< Noise params to animate. */
};

/**
 * @brief Parameters for a spherical-cap bump displacement field.
 */
struct BumpParams {
  Vector center;           /**< Bump center direction (unit vector). */
  Vector axis;             /**< Oriented stack axis the drape push acts along. */
  float radius = 0.5f;     /**< Angular radius of the bump footprint (radians). */
  float amplitude = 1.0f;  /**< Drape gain; the weight saturates at full boundary clearance for gains > 1. */
  float envelope = 0.0f;   /**< Footprint scale in [0, 1], animated by BallDrop. */
  float cos_radius = 1.0f; /**< Cached cos(radius) fast-reject bound. */

  /**
   * @brief Recomputes the cos(radius) fast-reject bound after a radius change.
   */
  void prepare_threshold() { cos_radius = cosf(std::min(radius, PI_F)); }

  /**
   * @brief Refreshes live-tunable config from a template snapshot.
   * @param t Template params carrying the current slider values.
   * @details Copies only the drape gain; the footprint and lifecycle fields
   * belong to the spawned fall.
   */
  void refresh_from(const BumpParams &t) { amplitude = t.amplitude; }

  /**
   * @brief Upper bound on |bump_field| for this entity.
   * @return The effective footprint radius (the drape push is capped at the
   * boundary clearance).
   */
  float field_bound() const { return radius * envelope; }

  /**
   * @brief Angular extent of the field's support (the cap radius).
   * @return radius * envelope. Coincides with field_bound() here (the drape push
   * never exceeds the footprint); the pair exists so a generalized prefilter can
   * treat balls and the poi field, whose shift can exceed its cap, uniformly.
   */
  float support_bound() const { return radius * envelope; }
};

/**
 * @brief Animates one bump falling from the world north pole to the south
 * pole along a fixed-azimuth meridian, in screen space.
 */
class BallDrop : public AnimationBase<BallDrop> {
public:
  /**
   * @brief Constructs a BallDrop animation.
   * @param params Bump params to animate. `params.amplitude` and
   *        `params.radius` are read from the spawn-time copy; set them before
   *        constructing.
   * @param orientation Frame whose oriented normal is the displaced stack's
   *        axis (the direction the bump pushes along); retained by pointer, so
   *        it must outlive the animation. The fall path itself is world-fixed.
   * @param normal Un-oriented stack axis.
   * @param azimuth Meridian the bump falls along (radians).
   * @param duration Fall time in frames.
   */
  BallDrop(BumpParams &params, const Orientation<> &orientation,
           const Vector &normal, float azimuth, int duration)
      : AnimationBase(duration, false), params(params),
        orientation(&orientation), normal(normal), azimuth(azimuth) {
    HS_CHECK(duration >= 2, "BallDrop duration must be >= 2");
    HS_CHECK(params.radius > 0.0f, "BallDrop needs a positive bump radius");
    params.prepare_threshold();
    params.envelope = 0.0f;
  }

  /**
   * @brief Steps the fall: advances the bump's colatitude down the world
   * frame, re-derives the displacement axis from the stack's current
   * orientation, and ramps the pole-edge envelope so the bump emerges from
   * and vanishes into the poles smoothly.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float progress = std::min(static_cast<float>(t) / duration, 1.0f);
    float phi = progress * PI_F;
    BumpParams &p = params.get();
    p.center = Vector(sinf(phi) * cosf(azimuth), cosf(phi),
                      sinf(phi) * sinf(azimuth));
    p.axis = make_basis(orientation->get(), normal).v;
    p.envelope = quintic_kernel(progress / EDGE_FRACTION) *
                 quintic_kernel((1.0f - progress) / EDGE_FRACTION);
  }

private:
  static constexpr float EDGE_FRACTION = 0.15f; /**< Fade window at either pole, as a fraction of the fall. */

  std::reference_wrapper<BumpParams> params; /**< Bump params to animate. */
  const Orientation<> *orientation; /**< Stack frame the push axis tracks; not owned. */
  Vector normal; /**< Un-oriented stack axis. */
  float azimuth; /**< Fall meridian (radians). */
};

struct PoiChoreography; // defined in motion.h; PoiDance reads it each frame

/**
 * @brief Parameters for a unidirectional spherical-cap drape field ("poi").
 * @details Like BumpParams but the drape pushes rings one way along the stack
 * axis (never parting them): the sheet slides over each dancer like silk over a
 * hand. The dome center and lifecycle swell are animated by PoiDance.
 */
struct PoiParams {
  Vector center;           /**< Dome center direction (unit vector). */
  Vector axis;             /**< Oriented stack axis the drape push acts along. */
  float radius = 0.3f;     /**< Footprint angular radius (radians), fixed at spawn. */
  float amplitude = 1.0f;  /**< Drape gain; the weight saturates at full clearance for gains > 1. */
  float direction = 1.0f;  /**< Push sign σ = ±1, fixed at phase entry. */
  float envelope = 0.0f;   /**< Lifecycle swell in [0, 1], animated by PoiDance. */
  float cos_radius = 1.0f; /**< Cached cos(radius) fast-reject bound. */

  /** @brief Peak of the drape profile (1−u)·quintic(u) over u∈[0,1], at u≈0.61. */
  static constexpr float PROFILE_MAX = 0.28f;

  /**
   * @brief Recomputes the cos(radius) fast-reject bound after a radius change.
   */
  void prepare_threshold() { cos_radius = cosf(std::min(radius, PI_F)); }

  /**
   * @brief Refreshes live-tunable config from a template snapshot.
   * @param t Template params carrying the current slider values.
   * @details Copies only the drape gain; footprint, direction, and lifecycle
   * belong to the spawned dancer.
   */
  void refresh_from(const PoiParams &t) { amplitude = t.amplitude; }

  /**
   * @brief Upper bound on |poi_field| (the drape shift) for this entity.
   * @return 2·r_eff·min(PROFILE_MAX·amplitude, 1): the unsaturated arm scales the
   * profile-max drape by the gain, the saturated arm caps at the full silhouette
   * shift 2·r_eff. Feeds the clip cull and the scan-band width, so it is kept
   * tight (the naive 2·r_eff·min(amp,1) is ~3.7× loose below saturation).
   */
  float field_bound() const {
    return 2.0f * radius * envelope * std::min(PROFILE_MAX * amplitude, 1.0f);
  }

  /**
   * @brief Angular extent of the field's support (the cap radius).
   * @return radius * envelope. Distinct from field_bound(): the poi shift can
   * exceed the cap radius at saturation, so the reach prefilter needs the cap
   * extent while the band widening needs the shift bound.
   */
  float support_bound() const { return radius * envelope; }
};

/**
 * @brief Animates one poi dancer: a dome tracing epicyclic tangent-plane curves
 * under the ring sheet, swelling up from and deflating back under it.
 * @details The dance center comes from the shared PoiChoreography (read by
 * pointer each frame) combined with this dancer's canon phase offsets; antipodal
 * partners (anchor_index ≥ 6) trace the exact point reflection of their primary.
 * step() is defined out-of-line in motion.h, where PoiChoreography is complete.
 */
class PoiDance : public AnimationBase<PoiDance> {
public:
  /**
   * @brief Constructs a PoiDance animation.
   * @param params Poi params to animate. `params.radius` and `params.direction`
   *        are read from the spawn-time copy; set them before constructing.
   * @param choreography Shared choreography stepped once per frame by the effect;
   *        retained by pointer, so it must outlive the animation.
   * @param orientation Frame whose oriented normal is the displaced stack's axis;
   *        retained by pointer. The dance floor itself is world-fixed.
   * @param normal Un-oriented stack axis.
   * @param anchor_index Dancer index in [0, 12): index % 6 selects the antipodal
   *        pair, index ≥ 6 marks the mirrored partner.
   * @param canon_offsets Per-term constant phase offsets (the canon lag), copied.
   * @param duration Dance lifetime in frames.
   */
  PoiDance(PoiParams &params, const PoiChoreography &choreography,
           const Orientation<> &orientation, const Vector &normal,
           int anchor_index, const float canon_offsets[3], int duration)
      : AnimationBase(duration, false), params(params),
        choreography(&choreography), orientation(&orientation), normal(normal),
        anchor_index(anchor_index) {
    HS_CHECK(duration >= 2, "PoiDance duration must be >= 2");
    HS_CHECK(params.radius > 0.0f, "PoiDance needs a positive poi radius");
    canon[0] = canon_offsets[0];
    canon[1] = canon_offsets[1];
    canon[2] = canon_offsets[2];
    params.prepare_threshold();
    params.envelope = 0.0f;
  }

  /**
   * @brief Steps the dance: ramps the lifecycle envelope, reads the shared
   * choreography to place the dome center, and re-derives the push axis from the
   * stack's current orientation. Defined in motion.h.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override;

private:
  static constexpr float EDGE_FRACTION = 0.10f; /**< Lifecycle swell window at either end of the dance. */

  std::reference_wrapper<PoiParams> params; /**< Poi params to animate. */
  const PoiChoreography *choreography; /**< Shared choreography; not owned. */
  const Orientation<> *orientation;    /**< Stack frame the push axis tracks; not owned. */
  Vector normal;    /**< Un-oriented stack axis. */
  int anchor_index; /**< Dancer index in [0, 12); index % 6 is the antipodal pair. */
  float canon[3];   /**< Per-term constant canon phase offsets. */
};

/**
 * @brief Parameters for a two-octave product noise field.
 */
struct NoiseProductParams {
  float amplitude = 0.0f; /**< Field output amplitude. */
  float scale1 = 1.5f;    /**< Spatial frequency of the envelope octave. */
  float scale2 = 3.0f;    /**< Spatial frequency of the detail octave. */
  float speed = 0.03f;    /**< Time advance per frame. */
  float time = 0.0f;      /**< Current field time, integrated by NoiseProduct. */
  mutable FastNoiseLite noise; /**< Backing generator. */

  /** @brief Spatial offset decorrelating octave 2 from octave 1 at equal scales. */
  static constexpr float OCTAVE2_OFFSET = 50.0f;

  /**
   * @brief Constructs params with an OpenSimplex2 generator at frequency 1.
   */
  NoiseProductParams() {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetFrequency(1.0f);
  }

  /**
   * @brief Refreshes live-tunable config from a template snapshot.
   * @param t Template params carrying the current slider values.
   * @details Copies the slider-driven fields but not the `time` axis or the
   * backing generator, so a live edit reaches a spawned entity without
   * resetting its phase.
   */
  void refresh_from(const NoiseProductParams &t) {
    amplitude = t.amplitude;
    scale1 = t.scale1;
    scale2 = t.scale2;
    speed = t.speed;
  }

  /**
   * @brief Upper bound on |noise_product_field| for this entity.
   * @return |amplitude| (|n1 * n2| <= 1).
   */
  float field_bound() const { return std::fabs(amplitude); }
};

/**
 * @brief Animates a noise-product field by integrating its time axis.
 * @details time += speed per frame keeps the field phase continuous under live
 * speed edits (an absolute time = t * speed would jump the phase).
 */
class NoiseProduct : public AnimationBase<NoiseProduct> {
public:
  /**
   * @brief Constructs a NoiseProduct animation (perpetual).
   * @param params Reference to the NoiseProductParams to animate.
   */
  NoiseProduct(NoiseProductParams &params)
      : AnimationBase(-1, true), params(params) {}

  /**
   * @brief Steps the animation, integrating the field time.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    params.get().time += params.get().speed;
  }

private:
  std::reference_wrapper<NoiseProductParams> params; /**< Noise params to animate. */
};

} // namespace Animation
