#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>
#include <array>
#include "3dmath.h"
#include "FastNoiseLite.h"

/**
 * @brief Frames Per Second constant.
 */
static constexpr int FPS = 16;

/**
 * @brief Easing function: Bi-Cubic Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_in_out_bicubic(float t) {
  return t < 0.5 ?
    4 * t * t * t :
    1 - (-2 * t + 2) * (-2 * t + 2) * (-2 * t + 2) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_in_out_sin(float t) {
  return -(cosf(PI_F * t) - 1) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_in_sin(float t) {
  return 1 - cosf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Sine Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_out_sin(float t) {
  return sinf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Cubic Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_in_cubic(float t) {
  return t * t * t;
}

/**
 * @brief Easing function: Circular Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_in_circ(float t) {
  return 1 - sqrtf(1 - t * t);
}

/**
 * @brief Easing function: Linear interpolation (no easing).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_mid(float t) {
  return t;
}

/**
 * @brief Easing function: Exponential Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_out_expo(float t) {
  return t == 1 ? 1 : 1 - powf(2.0f, -10 * t);
}

/**
 * @brief Easing function: Circular Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_out_circ(float t) {
  return sqrtf(1 - (t - 1) * (t - 1));
}

/**
 * @brief Easing function: Cubic easing out.
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_out_cubic(float t) {
  return 1 - powf(1 - t, 3);
}

/**
 * @brief Easing function: Elastic easing out.
 * @param x The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
float ease_out_elastic(float x) {
  const float c4 = (2 * PI_F) / 3;
  return x == 0 ?
    0 : x == 1 ?
    1 : powf(2, -10 * x) * sinf((x * 10 - 0.75f) * c4) + 1;
}

/**
 * @brief Base class for all animations, providing core timing and state management.
 * @tparam Derived The class inheriting from Animation (Curiously Recurring Template Pattern).
 */
template <typename Derived>
class Animation {
public:

  /**
   * @brief Cancels the animation on the next step.
   */
  void cancel() { canceled = true; }
  /**
   * @brief Checks if the animation is finished or canceled.
   * @return True if done.
   */
  virtual bool done() const { return canceled || (duration >= 0 && t >= duration); }
  /**
   * @brief Checks if the animation should repeat after finishing.
   * @return True if repeating.
   */
  virtual bool repeats() const { return repeat; }
  /**
   * @brief Advances the animation state by one frame.
   * @param canvas The canvas buffer (unused by base class, passed to derived classes).
   */
  virtual void step(Canvas& canvas) {
    t++;
  }

  /**
   * @brief Resets the internal timer to zero.
   */
  void rewind() {
    t = 0;
  }

  /**
   * @brief Sets a callback function to be executed when the animation finishes.
   * @param callback The function to execute on completion.
   * @return Reference to the derived animation object.
   */
  Derived& then(std::function<void()> callback) {
    post = callback;
    return static_cast<Derived&>(*this);
  }

  /**
   * @brief Executes the registered post-completion callback.
   */
  void post_callback() const {
    post();
  }

protected:

  /**
   * @brief Constructor for the base animation class.
   * @param duration Total number of frames the animation should run (-1 for indefinite).
   * @param repeat If true, the animation rewinds and restarts when finished.
   */
  Animation(int duration, bool repeat) :
    duration(duration == 0 ? 1 : duration),
    repeat(repeat),
    canceled(false),
    post([]() {})
  {
  }

  int duration; /**< Total length of the animation in frames. */
  bool repeat; /**< Flag indicating if the animation should repeat. */
  int t = 0; /**< Internal frame counter. */

private:

  bool canceled; /**< Flag to signal immediate cancellation. */
  std::function<void()> post; /**< Callback executed when the animation finishes. */
};

/**
 * @brief A placeholder animation that immediately reports being done.
 */
class NullAnimation : public Animation<NullAnimation> {
public:
  /**
   * @brief Constructs a NullAnimation.
   */
  NullAnimation() :
    Animation(0, false)
  {
  }
  /**
   * @brief Performs no action.
   */
  void step(Canvas&) {}
  /**
   * @brief Always reports true.
   */
  bool done() const { return true; }
};

/**
 * @brief An animation that triggers a callback after a random delay.
 */
class RandomTimer : public Animation<RandomTimer> {
public:

  /**
   * @brief Constructs a RandomTimer.
   * @param min Minimum delay in frames.
   * @param max Maximum delay in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  RandomTimer(int min, int max, TimerFn auto f, bool repeat = false) :
    Animation(-1, repeat),
    min(min),
    max(max),
    f(f),
    next(0)
  {
    reset();
  }

  /**
   * @brief Calculates the next random trigger time.
   */
  void reset() {
    next = t + static_cast<int>(std::round(hs::rand_int(min, max)));
  }

  /**
   * @brief Steps the timer, calling the function if the delay has elapsed.
   */
  void step(Canvas& canvas) {
    Animation::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
  }

private:

  int min; /**< Minimum frame delay. */
  int max; /**< Maximum frame delay. */
  std::function<void(Canvas&)> f; /**< The callback function. */
  int next; /**< The target frame count for the next trigger. */
};

/**
 * @brief An animation that triggers a callback at regular intervals.
 */
class PeriodicTimer : public Animation<PeriodicTimer> {
public:
  /**
   * @brief Constructs a PeriodicTimer.
   * @param period The interval between calls, in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  PeriodicTimer(int period, TimerFn auto f, bool repeat = false) :
    Animation(-1, repeat),
    period(period),
    f(f)
  {
    reset();
  }

  /**
   * @brief Calculates the next periodic trigger time.
   */
  void reset() {
    next = t + period;
  }

  /**
   * @brief Steps the timer, calling the function if the period has elapsed.
   */
  void step(Canvas& canvas) {
    Animation::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
  }

private:

  int period; /**< The interval in frames. */
  std::function<void(Canvas&)> f; /**< The callback function. */
  int next; /**< The target frame count for the next trigger. */
};

/**
 * @brief An animation that smoothly transitions a float variable over time.
 */
class Transition : public Animation<Transition> {
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
  Transition(float& mutant, float to, int duration, ScalarFn auto easing_fn, bool quantized = false, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    to(to),
    easing_fn(easing_fn),
    quantized(quantized)
  {
  }

  /**
   * @brief Performs one step of the transition.
   */
  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    Animation::step(canvas);
    auto t = std::min(1.0f, static_cast<float>(this->t) / duration);
    auto n = easing_fn(t) * (to - from) + from;
    if (quantized) {
      n = std::floor(n);
    }
    mutant.get() = n;
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   */
  void rebind_mutant(float& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<float> mutant; /**< Reference to the float variable being animated. */
  float from; /**< Starting value. */
  float to; /**< Target value. */
  std::function<float(float)> easing_fn; /**< Easing curve. */
  bool quantized; /**< Flag to round result to integer. */
};

/**
 * @brief An animation that applies a custom function to a float variable over time.
 */
class Mutation : public Animation<Mutation> {
public:

  /**
   * @brief Constructs a Mutation animation.
   * @param mutant The float variable to modify.
   * @param f The custom mutation function: `float f(float eased_t)`.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to apply to the time factor.
   * @param repeat If true, the mutation repeats indefinitely.
   */
  Mutation(float& mutant, ScalarFn auto f, int duration, ScalarFn auto easing_fn, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    f(f),
    easing_fn(easing_fn)
  {
  }

  /**
   * @brief Performs one step of the mutation.
   */
  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    Animation::step(canvas);
    auto t = std::min(1.0f, static_cast<float>(this->t) / duration);
    mutant.get() = f(easing_fn(t));
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   */
  void rebind_mutant(float& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<float> mutant; /**< Reference to the float variable being modified. */
  float from; /**< Starting value (unused by most MutateFns, but saved for context). */
  std::function<float(float)> f; /**< The custom function to apply. */
  std::function<float(float)> easing_fn; /**< Easing curve. */
};

/**
 * @brief An animation that draws a sprite while managing its fade-in/out effects.
 */
class Sprite : public Animation<Sprite> {
public:

  /**
   * @brief Constructs a Sprite animation.
   * @param draw_fn The function to call each frame for drawing (`void draw_fn(Canvas&, float opacity)`).
   * @param duration The duration the sprite is fully visible (-1 for indefinite).
   * @param fade_in_duration Frames for fading in.
   * @param fade_in_easing_fn Easing for fade-in.
   * @param fade_out_duration Frames for fading out.
   * @param fade_out_easing_fn Easing for fade-out.
   */
  Sprite(SpriteFn auto draw_fn, int duration,
    int fade_in_duration = 0, ScalarFn auto fade_in_easing_fn = ease_mid,
    int fade_out_duration = 0, ScalarFn auto fade_out_easing_fn = ease_mid) :
    Animation(duration, false),
    draw_fn(draw_fn),
    fader(fade_in_duration > 0 ? 0 : 1),
    fade_in_duration(fade_in_duration),
    fade_out_duration(fade_out_duration),
    fade_in(fader, 1, fade_in_duration, fade_in_easing_fn),
    fade_out(fader, 0, fade_out_duration, fade_out_easing_fn)
  {
  }

  /**
   * @brief Move constructor. Rebinds internal references to the new fader variable.
   */
  Sprite(Sprite&& other) noexcept
    : Animation(std::move(other)),
    draw_fn(std::move(other.draw_fn)),
    fader(other.fader),
    fade_in_duration(other.fade_in_duration),
    fade_out_duration(other.fade_out_duration),
    fade_in(std::move(other.fade_in)),
    fade_out(std::move(other.fade_out))
  {
    fade_in.rebind_mutant(this->fader);
    fade_out.rebind_mutant(this->fader);
  }

  /**
   * @brief Updates the drawing function used by the sprite.
   */
  void rebind_draw(SpriteFn auto new_draw_fn) {
    draw_fn = new_draw_fn;
  }

  /**
   * @brief Move assignment operator. Rebinds internal references.
   */
  Sprite& operator=(Sprite&& other) noexcept {
    if (this == &other) return *this;

    Animation::operator=(std::move(other));
    draw_fn = std::move(other.draw_fn);
    fader = other.fader;
    fade_in_duration = other.fade_in_duration;
    fade_out_duration = other.fade_out_duration;
    fade_in = std::move(other.fade_in);
    fade_in.rebind_mutant(this->fader);
    fade_out = std::move(other.fade_out);
    fade_out.rebind_mutant(this->fader);

    return *this;
  }

  /**
   * @brief Steps the animation, updates the fader, and calls the draw function with the current opacity.
   */
  void step(Canvas& canvas) {
    if (t == 0) {
      fade_in.rewind();
      fade_out.rewind();
    }
    Animation::step(canvas);
    if (!fade_in.done()) {
      fade_in.step(canvas);
    }
    else if (duration >= 0 && t >= (duration - fade_out_duration)) {
      fade_out.step(canvas);
    }
    draw_fn(canvas, fader);
  }

private:

  std::function<void(Canvas&, float)> draw_fn; /**< The drawing function functor. */
  float fader; /**< The variable storing the current opacity (0.0 to 1.0). */
  int fade_in_duration; /**< Duration of fade-in phase. */
  int fade_out_duration; /**< Duration of fade-out phase. */
  Transition fade_in; /**< Internal transition managing the fade-in. */
  Transition fade_out; /**< Internal transition managing the fade-out. */
};

/**
 * @brief An animation that moves an Orientation along a defined Path.
 * @details Updated to match JavaScript logic: steps along the path by sub-sampling based on the angle
 * and MAX_ANGLE, ensuring smooth turns even when the path curvature is high.
 * @tparam W The width of the LED display (used for calculating maximum rotation step).
 */
template <int W>
class Motion : public Animation<Motion<W>> {
public:

  /**
   * @brief Constructs a Motion animation.
   * @param orientation The Orientation object to update.
   * @param path The path to follow.
   * @param duration The duration in frames.
   * @param repeat If true, the motion repeats.
   */
  Motion(Orientation& orientation, const Path<W>& path, int duration, bool repeat = false) :
    Animation<Motion<W>>(duration, repeat),
    orientation(orientation),
    path(path)
  {
  }

  /**
   * @brief Steps the animation, calculates intermediate rotation steps along the path,
   * and pushes them to the Orientation.
   */
  void step(Canvas& canvas) {
    Animation<Motion<W>>::step(canvas);
    orientation.get().collapse();
    float t_prev = static_cast<float>(this->t - 1);
    Vector current_v = path.get().get_point(t_prev / this->duration);
    float t_curr = static_cast<float>(this->t);
    Vector target_v = path.get().get_point(t_curr / this->duration);
    float total_angle = angle_between(current_v, target_v);
    int num_steps = static_cast<int>(std::ceil(std::max(1.0f, total_angle / MAX_ANGLE)));
    Quaternion origin = orientation.get().get();
    for (int i = 1; i <= num_steps; ++i) {
      float sub_t = t_prev + (static_cast<float>(i) / num_steps);
      Vector next_v = path.get().get_point(sub_t / this->duration);
      float step_angle = angle_between(current_v, next_v);
      if (step_angle > TOLERANCE) {
        Vector step_axis = cross(current_v, next_v).normalize();
        Quaternion q = make_rotation(step_axis, step_angle);
        origin = (q * origin).normalize();
        orientation.get().push(origin);
      }
      current_v = next_v;
    }
  }

private:

  static constexpr float MAX_ANGLE = 2 * PI_F / W; /**< Maximum rotation angle per step to ensure smoothness. */
  std::reference_wrapper<Orientation> orientation; /**< Reference to the Orientation state. */
  std::reference_wrapper<const Path<W>> path; /**< Reference to the Path object. */
};

/**
 * @brief An animation that applies a fixed, time-eased rotation.
 * @details Updated to use INCREMENTAL rotations (delta) to allow multiple animations to overlap/add.
 * @tparam W The width of the LED display (used for calculating maximum rotation step).
 */
template <int W>
class Rotation : public Animation<Rotation<W>> {
public:

  /**
   * @brief Constructs a Rotation animation.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis (unit vector).
   * @param angle The total rotation angle in radians.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to use.
   * @param repeat If true, the rotation repeats.
   */
  Rotation(Orientation& orientation, const Vector& axis, float angle, int duration, ScalarFn auto easing_fn, bool repeat = false) :
    Animation<Rotation<W>>(duration, repeat),
    orientation(orientation),
    axis(axis),
    total_angle(angle),
    easing_fn(easing_fn),
    last_angle(0)
  {
  }

  /**
   * @brief Steps the animation, calculates the incremental rotation delta, and pushes it to the Orientation.
   */
  void step(Canvas& canvas) {
    if (this->t == 0) {
      last_angle = 0;
    }
    Animation<Rotation<W>>::step(canvas);
    orientation.get().collapse();
    float target_angle = easing_fn(static_cast<float>(this->t) / this->duration) * total_angle;
    float delta = target_angle - last_angle;
    if (std::abs(delta) > TOLERANCE) {
      int num_steps = static_cast<int>(std::ceil(std::abs(delta) / MAX_ANGLE));
      float step_angle = delta / num_steps;
      Quaternion q_step = make_rotation(axis, step_angle);
      for (int i = 0; i < num_steps; ++i) {
        Quaternion current_q = orientation.get().get();
        current_q = (q_step * current_q).normalize();
        orientation.get().push(current_q);
      }
      last_angle = target_angle;
    }
  }

  /**
   * @brief Static helper to apply a rotation immediately (duration = 1 frame).
   * @param canvas The canvas buffer.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis (unit vector).
   * @param angle The rotation angle in radians.
   * @param easing_fn The easing function to use.
   */
  static void animate(Canvas& canvas, Orientation& orientation, const Vector& axis, float_t angle, ScalarFn auto easing_fn) {
    Rotation<W> r(orientation, axis, angle, 1, easing_fn, false);
    r.step(canvas);
  }

private:

  static constexpr float MAX_ANGLE = 2 * PI_F / W; /**< Maximum rotation angle per step to ensure smoothness. */
  std::reference_wrapper<Orientation> orientation; /**< Reference to the Orientation state. */
  Vector axis; /**< The axis of rotation. */
  float total_angle; /**< The total angle to sweep. */
  std::function<float(float)> easing_fn; /**< Easing curve. */
  float last_angle; /**< The angle reached in the previous frame. */
};

/**
 * @brief An animation that simulates a particle/camera performing a random walk across the sphere's surface.
 * @details Uses Perlin noise to create continuous, turbulent pivoting motion.
 * @tparam W The width of the LED display.
 */
template<int W>
class RandomWalk : public Animation<RandomWalk<W>> {
public:
  /**
   * @brief Constructs a RandomWalk animation.
   * @param orientation The Orientation object to update.
   * @param v_start The starting direction vector.
   */
  RandomWalk(Orientation& orientation, const Vector& v_start) :
    Animation<RandomWalk<W>>(-1, false),
    orientation(orientation),
    v(Vector(v_start).normalize())
  {
    Vector u = X_AXIS;
    if (std::abs(dot(v, u)) > 0.99f) {
      u = Y_AXIS;
    }
    direction = cross(v, u).normalize();
    noiseGenerator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noiseGenerator.SetFrequency(NOISE_SCALE);
    noiseGenerator.SetSeed(hs::rand_int(std::numeric_limits<int>::min(), std::numeric_limits<int>::max()));
  }

  /**
   * @brief Steps the walk: pivots direction based on noise, then rotates the view along the calculated axis.
   */
  void step(Canvas& canvas) override {
    Animation<RandomWalk<W>>::step(canvas);
    float pivotAngle = noiseGenerator.GetNoise(static_cast<float>(this->t), 0.0f) * PIVOT_STRENGTH;
    direction = rotate(direction, make_rotation(v, pivotAngle)).normalize();
    Vector walk_axis = cross(v, direction).normalize();
    v = rotate(v, make_rotation(walk_axis, WALK_SPEED)).normalize();
    direction = rotate(direction, make_rotation(walk_axis, WALK_SPEED)).normalize();
    Rotation<W>::animate(canvas, orientation, walk_axis, WALK_SPEED, ease_mid);
  }

private:

  static constexpr float WALK_SPEED = 0.05; /**< Constant linear forward speed. */
  static constexpr float PIVOT_STRENGTH = 0.4; /**< Maximum side-to-side pivot magnitude based on noise. */
  static constexpr float NOISE_SCALE = 0.08; /**< Frequency of the noise used for pivoting. */

  FastNoiseLite noiseGenerator; /**< The noise generator instance. */
  std::reference_wrapper<Orientation> orientation; /**< Reference to the global Orientation state. */
  Vector v; /**< Current forward direction vector. */
  Vector direction; /**< Current pivoting direction (orthogonal to v). */
};

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief An animation that smoothly interpolates a GenerativePalette toward a target palette.
 */
class ColorWipe : public Animation<ColorWipe> {
public:
  /**
   * @brief Constructs a ColorWipe animation.
   * @param from_palette The GenerativePalette to start modifying (will be used for `cur_palette`).
   * @param to_palette The GenerativePalette instance to interpolate toward.
   * @param duration The duration in frames.
   * @param easing_fn The easing function.
   */
  ColorWipe(GenerativePalette& from_palette, const GenerativePalette& to_palette, int duration, ScalarFn auto easing_fn) :
    Animation(duration, false),
    from_palette(from_palette),
    cur_palette(from_palette),
    to_palette(to_palette),
    easing_fn(easing_fn)
  {
  }

  /**
   * @brief Steps the animation, blending the palette's colors based on the time factor.
   */
  void step(Canvas& canvas) {
    if (t == 0) {
      from_palette = cur_palette.get();
    }
    Animation::step(canvas);
    float amount = std::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    cur_palette.get().lerp(from_palette, to_palette, easing_fn(amount));
  }

private:

  GenerativePalette from_palette; /**< A local copy of the starting state of the palette. */
  std::reference_wrapper<GenerativePalette> cur_palette; /**< The actual palette instance being animated. */
  std::reference_wrapper<const GenerativePalette> to_palette; /**< The target final palette. */
  std::function<float(float)> easing_fn; /**< Easing curve. */
};

/**
 * @brief Animates the Mobius parameters for a continuous loxodromic flow.
 */
class MobiusFlow : public Animation<MobiusFlow> {
public:
  /**
   * @brief Constructs a MobiusFlow animation.
   * @param params Reference to the MobiusParams to animate.
   * @param num_rings Reference to the number of rings (scalar).
   * @param num_lines Reference to the number of lines (scalar).
   * @param duration Duration of the flow.
   * @param repeat Whether to repeat.
   */
  MobiusFlow(MobiusParams& params, const float& num_rings, const float& num_lines, int duration, bool repeat = true) :
    Animation(duration, repeat),
    params(params),
    num_rings(num_rings),
    num_lines(num_lines)
  {
  }

  /**
   * @brief Steps the animation, updating params a and d.
   */
  void step(Canvas& canvas) {
    Animation::step(canvas);
    float progress = static_cast<float>(t) / duration;
    float logPeriod = 5.0f / (num_rings + 1);
    float flowParam = progress * logPeriod;
    float scale = expf(flowParam);
    float s = sqrtf(scale);
    float angle = progress * (2 * PI_F / num_lines);

    params.get().aRe = s * cosf(angle);
    params.get().aIm = s * sinf(angle);
    params.get().dRe = (1.0f / s) * cosf(-angle);
    params.get().dIm = (1.0f / s) * sinf(-angle);
  }

private:
  std::reference_wrapper<MobiusParams> params;
  std::reference_wrapper<const float> num_rings;
  std::reference_wrapper<const float> num_lines;
};

/**
 * @brief Animates the Mobius parameters for a warping effect pulling the poles together.
 */
class MobiusWarp : public Animation<MobiusWarp> {
public:
  /**
   * @brief Constructs a MobiusWarp animation.
   * @param params Reference to the MobiusParams to animate.
   * @param num_rings Reference to the number of rings (scalar).
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   */
  MobiusWarp(MobiusParams& params, const float& num_rings, int duration, bool repeat = true) :
    Animation(duration, repeat),
    params(params),
    num_rings(num_rings)
  {
  }

  /**
   * @brief Steps the animation, updating param b.
   */
  void step(Canvas& canvas) {
    Animation::step(canvas);
    float progress = static_cast<float>(t) / duration;
    float angle = progress * 2 * PI_F;
    params.get().bRe = cosf(angle);
    params.get().bIm = sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params;
  std::reference_wrapper<const float> num_rings;
};

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Type alias for a variant that can hold any animation type.
 * @details Used by the Timeline to manage heterogeneous animation objects.
 */
using AnimationVariant = std::variant<
  NullAnimation,
  Sprite,
  Transition,
  Mutation,
  RandomTimer,
  PeriodicTimer,
  Rotation<MAX_W>,
  Motion<MAX_W>,
  RandomWalk<MAX_W>,
  ColorWipe,
  MobiusFlow,
  MobiusWarp
>;

/**
 * @brief Structure linking an animation variant with its starting time.
 */
struct TimelineEvent {
  int start; /**< The global frame count at which the animation should begin. */
  AnimationVariant animation; /**< The actual animation object. */
};

/**
 * @brief Manages all active animations and their execution over time.
 */
class Timeline {
public:

  /**
   * @brief Constructs a Timeline.
   */
  Timeline() :
    num_events(0)
  {
  }

  /**
   * @brief Adds a new animation event to the timeline.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting.
   * @param animation The animation object.
   * @return Reference to the Timeline object.
   */
  template <typename A>
  Timeline& add(float in_frames, A animation) {
    if (num_events >= MAX_EVENTS) {
      Serial.println("Timeline full, failed to add animation!");
      return *this;
    }
    TimelineEvent& e = events[num_events++];
    e.start = t + in_frames;
    e.animation = std::move(animation);
    return *this;
  }

  /**
   * @brief Advances the timeline by one frame, stepping all active or starting animations.
   * @param canvas The current canvas buffer.
   */
  void step(Canvas& canvas) {
    ++t;
    for (int i = 0; i < num_events; ++i) {
      if (t < events[i].start) {
        continue;
      }
      auto& animation = events[i].animation;
      std::visit([&](auto& a) { a.step(canvas); }, animation);
      if (std::visit([&](auto& a) { return a.done(); }, animation)) {
        if (std::visit([&](auto& a) { return a.repeats(); }, animation)) {
          std::visit([&](auto& a) { a.rewind(); }, animation);
          std::visit([&](auto& a) { a.post_callback(); }, animation);
        }
      }
    }

    auto new_logical_end = std::remove_if(
      events.begin(), events.begin() + num_events,
      [this](auto& event) {
        if (t < event.start) {
          return false;
        }
        if (std::visit([](auto& a) { return a.done(); }, event.animation)) {
          std::visit([](auto& a) { a.post_callback(); }, event.animation);
          return true;
        }
        return false;
      }
    );
    num_events = std::distance(events.begin(), new_logical_end);
  }

  int t = 0; /**< The current global frame count. */

private:

  static constexpr int MAX_EVENTS = 32; /**< Maximum number of concurrent animation events. */
  std::array<TimelineEvent, MAX_EVENTS> events; /**< Storage for all animation events. */
  int num_events; /**< Current number of active events. */
};