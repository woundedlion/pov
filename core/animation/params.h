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
        easing_fn(std::move(easing_fn)), paused(paused) {
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
   * @details By default the magnitude is the `scale` captured at construction.
   * Binding makes step() read the referent every frame instead, so a GUI slider
   * wired to that float takes effect immediately rather than only at the next
   * (re)spawn. The referent must outlive the animation (e.g. an effect's param
   * member). This avoids retaining the animation pointer across frames, which
   * would dangle under timeline compaction.
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
 * PERPETUAL — never completes: the default AnimationBase ctor gives it
 * duration = -1 with no repeat, so it never reaches done() and a `.then()`
 * callback attached to it NEVER fires (same hazard documented on RandomWalk).
 * Do not chain sequencing logic off its `.then()`, and note the Transformer's
 * slot-recycling `.then()` likewise will not fire for a spawned
 * MobiusWarpEvolving — drive follow-on behavior from a finite animation or
 * cancel() it explicitly.
 */
class MobiusWarpEvolving : public AnimationBase<MobiusWarpEvolving> {
public:
  /**
   * @brief Constructs a MobiusWarpEvolving animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   * @note `base` snapshots `params` at construction; MobiusParams has no
   * refresh_from, so Transformer::prepare_frame() never re-reads it from
   * template_params. The baseline is latched at spawn — live edits to it
   * require a respawn (live `scale`/`speed` are mutated via the setters).
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
    // animation is perpetual (duration == -1). A modulo/accumulator fix only
    // trades the slow drift for a sharp artifact.
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
   * @param params Reference to the params struct to animate.
   * @param center The center point of the ripple.
   * @param speed How fast the waves travel.
   * @param duration How long the ripple lasts in frames.
   */
  Ripple(RippleParams &params, const Vector &center, float speed = 0.2f,
         int duration = 100)
      : AnimationBase(duration, false), params(params), speed(speed),
        peak_amplitude(params.amplitude) {
    HS_CHECK(duration >= 0, "Ripple duration must be >= 0");
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
   * @details Copies the slider-driven fields (amplitude, speed, frequency,
   * scale) but not the animation-advanced `time` axis or the backing generator,
   * so a live edit reaches an already-spawned entity without resetting its
   * phase. prepare_frame() invokes this before sync().
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
   * @param duration Duration in frames (-1 for indefinite). Finite durations
   *   are not recommended: this animation repeats, so a finite duration rewinds
   *   t to 0 each cycle and snaps params.time backward — a sawtooth that breaks
   *   the smooth forward flow the time axis exists for. Use -1 (the default).
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
    // A modulo/accumulator fix only trades the slow drift for a sharp artifact.
    params.get().time = static_cast<float>(t);
  }

private:
  std::reference_wrapper<NoiseParams> params; /**< Noise params to animate. */
};

} // namespace Animation
