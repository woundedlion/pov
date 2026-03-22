/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_ANIMATION_H_
#define HOLOSPHERE_CORE_ANIMATION_H_

#include <algorithm>
#include <cmath>
#include <functional>

#include <numeric> // for std::iota
#include <array>
#include "3dmath.h"
#include "platform.h"
#include "FastNoiseLite.h"
#include "generators.h"
#include "geometry.h"
#include "memory.h"
#include "spatial.h"
#include "static_circular_buffer.h"
#include "rotate.h"
namespace Animation {
class Noise;
} // namespace Animation

/**
 * @brief Represents a customizable path.
 * Retains internal buffer for state, but draws to pipeline.
 */
template <int W, int RESOLUTION = 1024> class Path {
public:
  Path() {}

  /**
   * @brief Appends a procedurally generated segment to the path.
   * @param plot The function generating points.
   * @param domain The input domain scale.
   * @param samples The number of samples to take.
   * @param easing The easing function for sample distribution.
   * @return Reference to self for chaining.
   */
  Path &append_segment(PlotFn plot, float domain, float samples,
                       ScalarFn easing) {
    if (!points.is_empty())
      points.pop_back();
    for (float t = 0; t <= samples; t++) {
      points.push_back(plot(easing(t / samples) * domain));
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
    float raw_index = t * (points.size() - 1);
    size_t i = static_cast<size_t>(raw_index);
    float f = raw_index - i;
    if (i >= points.size() - 1)
      return points.back();
    const Vector &p1 = points[i];
    const Vector &p2 = points[i + 1];
    return p1 * (1.0f - f) + p2 * f;
  }

  void collapse() {
    if (points.size() > 1)
      points = {points.back()};
  }

private:
  StaticCircularBuffer<Vector, RESOLUTION> points;
};

/**
 * @brief Represents a path defined by a single procedural function.
 * Matches interface of Path for use in Motion animations.
 */
struct ProceduralPath {
  PlotFn f;

  ProceduralPath() = default;
  ProceduralPath(PlotFn path_fn) : f(std::move(path_fn)) {}

  Vector get_point(float t) const { return f(t); }
};

#include "easing.h"

/// Non-templated interface for all animations.
/// Enables virtual dispatch in Timeline without std::visit.
class IAnimation {
public:
  virtual ~IAnimation() = default;
  virtual void step(Canvas &canvas) = 0;
  virtual bool done() const = 0;
  virtual bool repeats() const = 0;
  virtual void rewind() = 0;
  virtual void post_callback() = 0;
  /// Override in types that own an Orientation to collapse it.
  virtual void collapse_orientation() {}
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
    return canceled || (duration >= 0 && t >= duration);
  }
  /**
   * @brief Checks if the animation should repeat after finishing.
   * @return True if repeating.
   */
  bool repeats() const override { return repeat; }
  /**
   * @brief Advances the animation state by one frame.
   * @param canvas The canvas buffer (unused by base class, passed to derived
   * classes).
   */
  void step(Canvas &canvas) override { t++; }

  /**
   * @brief Resets the internal timer to zero.
   */
  void rewind() override { t = 0; }

  /**
   * @brief Sets a callback function to be executed when the animation finishes.
   * @param callback The function to execute on completion.
   * @return LValue Reference to the derived animation object.
   */
  Derived &then(Fn<void(), 24> callback) & {
    post = std::move(callback);
    return static_cast<Derived &>(*this);
  }

  /**
   * @brief Sets a callback function to be executed when the animation finishes
   * (RValue overload).
   * @param callback The function to execute on completion.
   * @return RValue Reference to the derived animation object.
   */
  Derived &&then(Fn<void(), 24> callback) && {
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
        canceled(false) {}

  AnimationBase() : duration(-1), repeat(false), canceled(false) {}

  int duration; /**< Total length of the animation in frames. */
  bool repeat;  /**< Flag indicating if the animation should repeat. */
  int t = 0;    /**< Internal frame counter. */

private:
  bool canceled;       /**< Flag to signal immediate cancellation. */
  Fn<void(), 24> post; /**< Callback executed when the animation finishes. */
};

/**
 * @brief Manages a history of Orientation states.
 * @tparam CAPACITY The maximum number of snapshots to keep.
 */
template <typename OrientationType, int CAPACITY> class OrientationTrail {
public:
  /**
   * @brief Records a snapshot of the current orientation state.
   * @param source The orientation to copy.
   */
  void record(const OrientationType &source) { snapshots.push_back(source); }

  /**
   * @brief Gets the number of recorded snapshots.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Gets a specific snapshot.
   * @param i Index (0 is newest).
   */
  const OrientationType &get(size_t i) const {
    // JS parity: 0 is oldest
    return snapshots[i];
  }

  /**
   * @brief Gets a specific snapshot (mutable).
   * @param i Index (0 is newest).
   */
  OrientationType &get(size_t i) { return snapshots[i]; }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop(); }

private:
  StaticCircularBuffer<OrientationType, CAPACITY> snapshots;
};

/**
 * @brief Manages a history of world-space Vector positions.
 * @tparam CAPACITY The maximum number of snapshots to keep.
 */
template <int CAPACITY> class VectorTrail {
public:
  void record(const Vector &source) { snapshots.push_back(source); }
  size_t length() const { return snapshots.size(); }
  const Vector &get(size_t i) const { return snapshots[i]; }
  Vector &get(size_t i) { return snapshots[i]; }
  void clear() { snapshots.clear(); }
  void expire() { snapshots.pop(); }

private:
  StaticCircularBuffer<Vector, CAPACITY> snapshots;
};

/**
 * @brief Represents a single particle in a system.
 */
template <int TRAIL_LEN = 8> struct Particle {
  Vector position;         /**< Current 3D position. */
  Vector velocity;         /**< Current velocity vector. */
  uint16_t color_seed = 0; /**< Hue seed for palette offset. */
  uint16_t life;           /**< Remaining life (frames or arbitrary units). */

  VectorTrail<TRAIL_LEN> history; /**< Trail of world-space positions. */

  void init(const Vector &p, const Vector &v, uint16_t seed, float l) {
    position = p;
    velocity = v;
    color_seed = seed;
    life = l;
    history.clear();
  }

  size_t history_length() const { return history.length(); }
};

/**
 * @brief A physics-based particle system with emitters and attractors.
 * @tparam W Width of the display (for Orientation).
 */
template <int W, int CAPACITY, int TRAIL_LEN = 8, int EMITTER_CAP = 8,
          int ATTRACTOR_CAP = 8>
class ParticleSystem
    : public AnimationBase<
          ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP>> {
public:
  ArenaVector<Particle<TRAIL_LEN>> pool;
  uint16_t active_count = 0;

  float friction = 0.85f;
  float gravity = 0.001f;
  uint16_t max_life = 600;

  struct Attractor {
    Vector position;
    float strength;
    float kill_radius;
    float event_horizon;
  };

  using EmitterFn = Fn<void(ParticleSystem &), 32>;

  ArenaVector<Attractor> attractors;
  ArenaVector<EmitterFn> emitters;

  ParticleSystem()
      : AnimationBase<
            ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP>>(
            -1, false) {}

  void init(Arena &arena, float friction = 0.85f, float gravity = 0.001f,
            float max_life = 600.0f) {
    this->friction = friction;
    this->gravity = gravity;
    this->max_life = max_life;
    active_count = 0;
    pool.bind(arena, CAPACITY);
    if (!pool.is_bound()) return;  // OOM — particle system stays inactive
    for (size_t i = 0; i < CAPACITY; ++i) {
      pool.emplace_back();
    }
    attractors.bind(arena, ATTRACTOR_CAP);
    emitters.bind(arena, EMITTER_CAP);
  }

  void add_emitter(EmitterFn fn) { emitters.push_back(fn); }

  void add_attractor(const Vector &pos, float str, float kill, float horizon) {
    attractors.push_back({pos, str, kill, horizon});
  }

  /**
   * @brief Spawns a new particle with a color seed.
   * @param pos Initial position.
   * @param vel Initial velocity.
   * @param seed Color seed for palette offset.
   */
  void spawn(const Vector &pos, const Vector &vel, uint16_t seed) {
    if (!pool.is_bound()) return;
    if (active_count < pool.capacity()) {
      pool[active_count++].init(pos, vel, seed, max_life);
    }
  }

  void step(Canvas &canvas) override {
    AnimationBase<ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP,
                                 ATTRACTOR_CAP>>::step(canvas);

    {
      // Emitters
      for (size_t i = 0; i < emitters.size(); ++i) {
        emitters[i](*this);
      }

      float max_delta = (2 * PI_F) / W;

      // Physics
      for (size_t i = 0; i < active_count; ++i) {
        bool dead = step_particle(pool[i], max_delta);
        if (dead) {
          pool[i] = pool[active_count - 1];
          active_count--;
          i--;
        }
      }
    }
  }

private:
  bool step_particle(Particle<TRAIL_LEN> &p, float max_delta) {
    bool active = p.life > 0;
    if (active) {
      p.life--;
      active = p.life > 0;
    }

    // Physics
    if (active) {
      Vector pos = p.position;

      for (size_t k = 0; k < attractors.size(); ++k) {
        const auto &attr = attractors[k];
        float dist_sq = distance_squared(pos, attr.position);

        if (dist_sq < attr.kill_radius * attr.kill_radius) {
          active = false;
          break; // Killed
        }

        if (dist_sq > 0.0000001f) {
          if (dist_sq < attr.event_horizon * attr.event_horizon) {
            // Steer into center
            Vector torque = (attr.position - pos).normalized();
            float speed = p.velocity.magnitude();
            p.velocity = torque * speed;
          } else {
            // Gravity
            float force = (gravity * attr.strength) / dist_sq;
            Vector torque = cross(pos, attr.position).normalized() * force;
            p.velocity += cross(torque, pos);
          }
        }
      }

      if (active) {
        // Drag
        p.velocity *= friction;

        // Move
        float speed = p.velocity.magnitude();
        if (speed > 0.000001f) {
          Vector axis = cross(pos, p.velocity).normalized();
          Quaternion dq = make_rotation(axis, speed);
          p.position = rotate(p.position, dq);
          p.velocity = rotate(p.velocity, dq);
        }
      }
    }

    // History Management
    if (active) {
      p.history.record(p.position);
    } else {
      if (p.history.length() > 0) {
        p.history.expire();
      }
    }

    return !active && p.history.length() == 0;
  }
};

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
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      } else {
        cancel();
      }
    }
  }

private:
  int min;   /**< Minimum frame delay. */
  int max;   /**< Maximum frame delay. */
  TimerFn f; /**< The callback function. */
  int next;  /**< The target frame count for the next trigger. */
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
      : AnimationBase(-1, repeat), period(period), f(std::move(f)) {
    reset();
  }

  /**
   * @brief Calculates the next periodic trigger time.
   */
  void reset() { next = t + period; }

  /**
   * @brief Steps the timer, calling the function if the period has elapsed.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      } else {
        cancel();
      }
    }
  }

private:
  int period; /**< The interval in frames. */
  TimerFn f;  /**< The callback function. */
  int next;   /**< The target frame count for the next trigger. */
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
        easing_fn(std::move(easing_fn)), quantized(quantized) {}

  /**
   * @brief Performs one step of the transition.
   */
  void step(Canvas &canvas) override {
    if (t == 0) {
      from = mutant;
    }
    AnimationBase::step(canvas);
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
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

private:
  std::reference_wrapper<float>
      mutant;         /**< Reference to the float variable being animated. */
  float from;         /**< Starting value. */
  float to;           /**< Target value. */
  EasingFn easing_fn; /**< Easing curve. */
  bool quantized;     /**< Flag to round result to integer. */
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
   */
  Mutation(float &mutant, ScalarFn f, int duration, EasingFn easing_fn,
           bool repeat = false)
      : AnimationBase(duration, repeat), mutant(mutant), from(mutant),
        f(std::move(f)), easing_fn(std::move(easing_fn)) {}

  /**
   * @brief Performs one step of the mutation.
   */
  void step(Canvas &canvas) override {
    if (t == 0) {
      from = mutant;
    }
    AnimationBase::step(canvas);
    auto t = std::min(1.0f, static_cast<float>(this->t) / duration);
    mutant.get() = f(easing_fn(t));
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   * @param new_mutant The new float variable to modify.
   */
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

private:
  std::reference_wrapper<float>
      mutant; /**< Reference to the float variable being modified. */
  float from; /**< Starting value (unused by most MutateFns, but saved for
                 context). */
  ScalarFn f; /**< The custom function to apply. */
  EasingFn easing_fn; /**< Easing curve. */
};

/**
 * @brief An animation that continuously increments a float variable over time.
 */
class Driver : public AnimationBase<Driver> {
public:
  /**
   * @brief Constructs a Driver animation for continuous progression.
   * @param mutant The float variable to modify.
   * @param speed The amount to add per frame.
   */
  Driver(float &mutant, float speed, bool wrap = true)
      : AnimationBase(1, true), mutant(mutant), speed(speed), wrap_(wrap) {}

  /**
   * @brief Performs one step by adding the speed to the mutant.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    mutant.get() += speed;
    if (wrap_) {
      if (mutant.get() >= 1.0f)
        mutant.get() -= 1.0f;
      else if (mutant.get() < 0.0f)
        mutant.get() += 1.0f;
    }
  }

  /**
   * @brief Rebinds the reference to the float variable being mutated.
   */
  void rebind_mutant(float &new_mutant) { mutant = new_mutant; }

  /**
   * @brief Gets current speed.
   */
  float get_speed() const { return speed; }

  /**
   * @brief Sets current speed.
   */
  void set_speed(float new_speed) { speed = new_speed; }

  const float &get_mutant() const { return mutant.get(); }

private:
  std::reference_wrapper<float> mutant; /**< Reference to the float variable. */
  float speed;                          /**< Amount added per frame. */
  bool wrap_;                           /**< If true, wraps value to 0-1 range. */
};

/**
 * @brief An animation that interpolates between states. The caller owns the
 * start, subject, and target data. Lerp just holds pointers and a type-erased
 * lerp function. Supports any type T that implements lerp(start, target, t).
 */
class Lerp : public AnimationBase<Lerp> {
public:
  template <typename T, typename Easing>
  Lerp(T &subject, const T &start, const T &target, int duration,
       Easing easing_fn)
      : AnimationBase(duration, false), subject_ptr(&subject),
        start_ptr(&start), target_ptr(&target), easing(easing_fn) {
    do_lerp = [](void *subj, const void *s, const void *tgt, float t) {
      static_cast<T *>(subj)->lerp(*static_cast<const T *>(s),
                                   *static_cast<const T *>(tgt), t);
    };
  }

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    do_lerp(subject_ptr, start_ptr, target_ptr, easing(progress));
  }

private:
  void *subject_ptr;
  const void *start_ptr;
  const void *target_ptr;
  EasingFn easing;
  void (*do_lerp)(void *, const void *, const void *, float);
};

/**
 * @brief An animation that draws a sprite while managing its fade-in/out
 * effects. Computes opacity inline rather than embedding Transition objects.
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
   */
  Sprite(SpriteFn draw_fn, int duration, int fade_in_duration = 0,
         EasingFn fade_in_easing_fn = ease_mid, int fade_out_duration = 0,
         EasingFn fade_out_easing_fn = ease_mid)
      : AnimationBase(duration, false), draw_fn(std::move(draw_fn)),
        fade_in_duration(fade_in_duration),
        fade_out_duration(fade_out_duration),
        fade_in_easing(std::move(fade_in_easing_fn)),
        fade_out_easing(std::move(fade_out_easing_fn)) {}

  /**
   * @brief Updates the drawing function used by the sprite.
   */
  void rebind_draw(SpriteFn new_draw_fn) { draw_fn = std::move(new_draw_fn); }

  /**
   * @brief Steps the animation, computes the current opacity inline, and calls
   * the draw function.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float opacity = 1.0f;

    // Fade in
    if (fade_in_duration > 0 && t < fade_in_duration) {
      float progress = static_cast<float>(t) / fade_in_duration;
      opacity = fade_in_easing(hs::clamp(progress, 0.0f, 1.0f));
    }
    // Fade out
    else if (duration >= 0 && fade_out_duration > 0 &&
             t >= (duration - fade_out_duration)) {
      float elapsed = static_cast<float>(t - (duration - fade_out_duration));
      float progress = elapsed / fade_out_duration;
      opacity = 1.0f - fade_out_easing(hs::clamp(progress, 0.0f, 1.0f));
    }

    draw_fn(canvas, opacity);
  }

private:
  SpriteFn draw_fn;         /**< The drawing function functor. */
  int fade_in_duration;     /**< Duration of fade-in phase in frames. */
  int fade_out_duration;    /**< Duration of fade-out phase in frames. */
  EasingFn fade_in_easing;  /**< Easing curve for fade-in. */
  EasingFn fade_out_easing; /**< Easing curve for fade-out. */
};

/**
 * @brief An animation that moves an Orientation along a defined Path.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 */
template <int W, int CAP = 4>
class Motion : public AnimationBase<Motion<W, CAP>> {
public:
  /**
   * @brief Constructs a Motion animation.
   * @param orientation The Orientation object to update.
   * @param path_obj The path to follow (must have get_point(float t)).
   * @param duration The duration in frames.
   * @param repeat If true, the motion repeats.
   */
  template <typename P>
  Motion(Orientation<W, CAP> &orientation, const P &path_obj, int duration,
         bool repeat = false, Space space = Space::World)
      : AnimationBase<Motion<W, CAP>>(duration, repeat),
        orientation(orientation),
        path_fn([&path_obj](float t) { return path_obj.get_point(t); }),
        space(space) {}

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W, CAP> &get_orientation() const { return orientation.get(); }
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Steps the animation, calculates intermediate rotation steps along
   * the path, and pushes them to the Orientation.
   */
  void step(Canvas &canvas) override {
    AnimationBase<Motion<W, CAP>>::step(canvas);
    float t_prev = static_cast<float>(this->t - 1);
    Vector current_v = path_fn(t_prev / this->duration);
    float t_curr = static_cast<float>(this->t);
    Vector target_v = path_fn(t_curr / this->duration);
    float total_angle = angle_between(current_v, target_v);
    int num_steps =
        static_cast<int>(std::ceil(std::max(1.0f, total_angle / MAX_ANGLE)));

    // Ensure sufficient resolution
    orientation.get().upsample(num_steps + 1);
    int len = orientation.get().length();

    // lambda to apply rotation based on space
    auto apply_rotation = [&](Quaternion &target, const Quaternion &source) {
      if (space == Space::Local) {
        target = target * source;
      } else {
        target = source * target;
      }
    };

    Vector prev_v = current_v;
    Quaternion accumulated_q; // Identity by default

    for (int i = 1; i < len; ++i) {
      // Calculate sub-t to sample path
      // t goes from (t-1) to t.
      // i goes from 0 to len-1. i=0 is t-1, i=len-1 is t.
      float sub_t = t_prev + (static_cast<float>(i) / (len - 1));

      Vector next_v = path_fn(sub_t / this->duration);
      float step_angle = angle_between(prev_v, next_v);

      if (step_angle > TOLERANCE) {
        Vector step_axis = cross(prev_v, next_v).normalized();
        Quaternion q_step = make_rotation(step_axis, step_angle);
        apply_rotation(accumulated_q, q_step);
      }

      // Apply the accumulated rotation to the existing orientation frame
      Quaternion &current_q = orientation.get().at(i);
      apply_rotation(current_q, accumulated_q);
      current_q.normalize();

      prev_v = next_v;
    }
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  std::reference_wrapper<Orientation<W, CAP>>
      orientation;               /**< Reference to the Orientation state. */
  Fn<Vector(float), 16> path_fn; /**< Function to retrieve path points. */
  Space space;                   /**< The coordinate space for rotation. */
};

/**
 * @brief An animation that applies a fixed, time-eased rotation.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 */
template <int W, int CAP = 4>
class Rotation : public AnimationBase<Rotation<W, CAP>> {
public:
  /**
   * @brief Default constructor. Creates an inactive/identity rotation.
   */
  Rotation()
      : AnimationBase<Rotation<W, CAP>>(0, false), orientation(nullptr),
        axis(X_AXIS), total_angle(0), easing_fn(ease_mid), last_angle(0),
        space(Space::World) {}

  /**
   * @brief Constructs a Rotation animation.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis (unit vector).
   * @param angle The total rotation angle in radians.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to use.
   * @param repeat If true, the rotation repeats.
   * @param space The coordinate space for rotation ("World" or "Local").
   */
  Rotation(Orientation<W, CAP> &orientation, const Vector &axis, float angle,
           int duration, EasingFn easing_fn, bool repeat = false,
           Space space = Space::World)
      : AnimationBase<Rotation<W, CAP>>(duration, repeat),
        orientation(&orientation), axis(axis), total_angle(angle),
        easing_fn(std::move(easing_fn)), last_angle(0), space(space) {}

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W, CAP> &get_orientation() const { return *orientation; }
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Check if the rotation has a valid orientation bound.
   */
  bool has_orientation() const { return orientation != nullptr; }

  /**
   * @brief Steps the animation, calculates the incremental rotation delta, and
   * pushes it to the Orientation.
   */
  void step(Canvas &canvas) override {
    if (this->t == 0) {
      last_angle = 0;
    }
    AnimationBase<Rotation<W, CAP>>::step(canvas);
    float target_angle =
        easing_fn(static_cast<float>(this->t) / this->duration) * total_angle;
    float delta = target_angle - last_angle;

    if (std::abs(delta) < TOLERANCE) {
      last_angle = target_angle;
      return;
    }

    int num_steps =
        1 + static_cast<int>(std::ceil(std::abs(delta) / MAX_ANGLE));
    orientation->upsample(num_steps + 1);
    int len = orientation->length();

    float step_angle = delta / (len - 1);
    Quaternion accumulated_q;

    // lambda to apply rotation based on space
    auto apply_rotation = [&](Quaternion &target, const Quaternion &source) {
      if (space == Space::Local) {
        target = target * source;
      } else {
        target = source * target;
      }
    };

    for (int i = 1; i < len; ++i) {
      float angle = step_angle * i;
      Quaternion q = make_rotation(axis, angle);

      Quaternion &current_q = orientation->at(i);
      apply_rotation(current_q, q);
      current_q.normalize();
    }
    last_angle = target_angle;
  }

  /**
   * @brief Static helper to apply a rotation immediately (duration = 1 frame).
   * @param canvas The canvas buffer.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis (unit vector).
   * @param angle The rotation angle in radians.
   * @param easing_fn The easing function to use.
   * @param space The coordinate space for rotation.
   */
  static void animate(Canvas &canvas, Orientation<W, CAP> &orientation,
                      const Vector &axis, float_t angle, EasingFn easing_fn,
                      Space space = Space::World) {
    Rotation<W, CAP> r(orientation, axis, angle, 1, easing_fn, false, space);
    r.step(canvas);
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  Orientation<W, CAP> *orientation; /**< Pointer to the Orientation state. */
  Vector axis;                      /**< The axis of rotation. */
  float total_angle;                /**< The total angle to sweep. */
  EasingFn easing_fn;               /**< Easing curve. */
  float last_angle; /**< The angle reached in the previous frame. */
  Space space;      /**< The coordinate space for rotation. */
};

/**
 * @brief An animation that simulates a particle/camera performing a random walk
 * across the sphere's surface.
 * @details Uses Perlin noise to create continuous, turbulent pivoting motion.
 * @tparam W The width of the LED display.
 */
template <int W, int CAP = 4>
class RandomWalk : public AnimationBase<RandomWalk<W, CAP>> {
public:
  struct Options {
    float speed = 0.02f; /**< Movement speed per frame. */
    float pivot_strength =
        0.1f; /**< Strength of the direction change (noise amplitude). */
    float noise_scale = 0.02f;  /**< Frequency of the Perlin noise. */
    float smoothing = 0.85f;    /**< Angular momentum (0 = none, 0.95 = very sluggish). */
    float drift = 0.5f;         /**< Temporal drift speed for spatial noise. */

    static Options Languid() { return {0.02f, 0.1f, 0.02f, 0.85f, 0.5f}; }
    static Options Energetic() { return {0.05f, 0.4f, 0.08f, 0.7f, 1.0f}; }
  };

  /**
   * @brief Constructs a RandomWalk animation.
   * @param orientation The Orientation object to update.
   * @param v_start The starting direction vector.
   * @param noise External noise generator (caller owns lifetime).
   * @param options Configuration options.
   */
  RandomWalk(Orientation<W, CAP> &orientation, const Vector &v_start,
             FastNoiseLite &noise, Options options = Options(),
             int seed = 0)
      : AnimationBase<RandomWalk<W, CAP>>(-1, false), orientation(orientation),
        v(Vector(v_start).normalized()), options(options),
        noiseGenerator(noise) {
    Vector u = X_AXIS;
    if (std::abs(dot(v, u)) > 0.99f) {
      u = Y_AXIS;
    }
    direction = cross(v, u).normalized();
    noiseGenerator.get().SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noiseGenerator.get().SetFrequency(options.noise_scale);
    if (seed == 0) {
      noiseGenerator.get().SetSeed(std::rand());
    } else {
      noiseGenerator.get().SetSeed(seed);
    }
  }

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W, CAP> &get_orientation() const { return orientation.get(); }
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Steps the walk: pivots direction based on noise, then rotates the
   * view along the calculated axis.
   */
  void step(Canvas &canvas) override {
    AnimationBase<RandomWalk<W, CAP>>::step(canvas);
    float target_pivot =
        noiseGenerator.get().GetNoise(
            v.x * options.noise_scale * 100.0f,
            v.y * options.noise_scale * 100.0f,
            v.z * options.noise_scale * 100.0f +
                static_cast<float>(this->t) * options.drift) *
        options.pivot_strength;
    angular_velocity = angular_velocity * options.smoothing +
                       target_pivot * (1.0f - options.smoothing);
    direction = rotate(direction, make_rotation(v, angular_velocity)).normalized();
    Vector walk_axis = cross(v, direction).normalized();
    v = rotate(v, make_rotation(walk_axis, options.speed)).normalized();
    direction =
        rotate(direction, make_rotation(walk_axis, options.speed)).normalized();
    Rotation<W, CAP>::animate(canvas, orientation, walk_axis, options.speed,
                              ease_mid);
  }

private:
  std::reference_wrapper<Orientation<W, CAP>>
      orientation;  /**< Reference to the global Orientation state. */
  Vector v;         /**< Current forward direction vector. */
  Vector direction; /**< Current pivoting direction (orthogonal to v). */
  Options options;  /**< Configuration options. */
  float angular_velocity = 0.0f; /**< Smoothed pivot rate (angular momentum). */
  std::reference_wrapper<FastNoiseLite>
      noiseGenerator; /**< External noise generator. */
};

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief An animation that smoothly interpolates a GenerativePalette toward a
 * target palette.
 */
class ColorWipe : public AnimationBase<ColorWipe> {
public:
  /**
   * @brief Constructs a ColorWipe animation.
   * @param from_palette The GenerativePalette to animate (snapshot taken at
   * t=0).
   * @param to_palette The GenerativePalette to interpolate toward.
   * @param duration The duration in frames.
   * @param easing_fn The easing function.
   */
  ColorWipe(GenerativePalette &from_palette,
            const GenerativePalette &to_palette, int duration,
            EasingFn easing_fn)
      : AnimationBase(duration, false), cur_palette(from_palette),
        to_snap(to_palette.snapshot()), easing_fn(std::move(easing_fn)) {}

  /**
   * @brief Steps the animation, blending the palette's colors based on the time
   * factor.
   */
  void step(Canvas &canvas) override {
    if (t == 0) {
      from_snap = cur_palette.get().snapshot();
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
        num_lines(num_lines) {}

  /**
   * @brief Steps the animation, updating params a and d.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
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
        easing(easing) {}

  /**
   * @brief Steps the animation, updating param b.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    params.get().bRe = scale * (cosf(angle) - 1.0f);
    params.get().bIm = scale * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params;

public:
  float scale;
  EasingFn easing;
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
        easing(easing) {}

  /**
   * @brief Steps the animation, updating param b.
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    params.get().bRe = scale * cosf(angle);
    params.get().bIm = -scale * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params;

public:
  float scale;
  EasingFn easing;
};



/**
 * @brief Animates a vertex-interpolated crossfade between two meshes.
 *
 * Owns its transient state (cloned meshes + SLERP buffers) via a pointer
 * to arena-allocated storage — keeps inline size small for TimelineEvent.
 * The caller provides source/dest MeshState references and an Arena;
 * MeshMorph clones both, builds nearest-vertex correspondence, and
 * interpolates each frame. Transients are released when the animation
 * completes and the TimelineEvent is destroyed; call compact on the
 * arena afterward to reclaim the space.
 */
class MeshMorph : public AnimationBase<MeshMorph> {
public:
  using MorphDrawFn = FunctionRef<void(Canvas &, const MeshState &, float)>;

  /// Two-callback constructor: separate shading for outgoing and incoming.
  MeshMorph(const MeshState &source, const MeshState &dest, Arena &arena,
            MorphDrawFn draw_outgoing, MorphDrawFn draw_incoming, int duration,
            EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(duration, false), easing_fn(easing_fn),
        draw_outgoing(draw_outgoing), draw_incoming(draw_incoming) {
    // Allocate transient storage on the arena
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();

    // Clone both meshes for interpolation
    MeshOps::clone(source, buf_->mesh_A, arena);
    MeshOps::clone(dest, buf_->mesh_B, arena);

    // Allocate SLERP buffers
    buf_->start_pos.bind(arena, dest.vertices.size());
    buf_->end_pos.bind(arena, dest.vertices.size());

    // Symmetry-breaking twist to avoid degenerate nearest-vertex mapping
    Vector twist_axis = Vector(0.0f, 0.0f, 1.0f);
    bool has_poles = false;
    for (const auto &v : source.vertices) {
      if (std::abs(v.z) > 0.99f && std::abs(v.x) < 0.01f)
        has_poles = true;
    }
    if (has_poles) {
      twist_axis = Vector(1.0f, 1.0f, 1.0f).normalized();
    }
    Quaternion twist = make_rotation(twist_axis, 0.05f);

    // Build nearest-vertex correspondence
    for (size_t i = 0; i < dest.vertices.size(); ++i) {
      Vector v_biased = rotate(dest.vertices[i], twist);
      int best_idx = 0;
      float max_dot = -9999.0f;
      for (size_t j = 0; j < source.vertices.size(); ++j) {
        float d = dot(v_biased, source.vertices[j]);
        if (d > max_dot) {
          max_dot = d;
          best_idx = static_cast<int>(j);
        }
      }
      buf_->start_pos.push_back(source.vertices[best_idx]);
      buf_->end_pos.push_back(dest.vertices[i]);
    }
  }

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);

    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    float alpha = easing_fn(progress);

    // Interpolate vertices
    for (size_t i = 0; i < buf_->end_pos.size(); ++i) {
      buf_->mesh_B.vertices[i] =
          slerp(buf_->start_pos[i], buf_->end_pos[i], alpha);
    }

    // Render crossfade
    float op_A = 1.0f - alpha;
    if (op_A > 0.01f)
      draw_outgoing(canvas, buf_->mesh_A, op_A);
    if (alpha > 0.01f)
      draw_incoming(canvas, buf_->mesh_B, alpha);
  }

private:
  /// Arena-allocated transient data — keeps MeshMorph inline size small.
  struct Transients {
    MeshState mesh_A;
    MeshState mesh_B;
    ArenaVector<Vector> start_pos;
    ArenaVector<Vector> end_pos;
  };

  Transients *buf_;
  EasingFn easing_fn;
  MorphDrawFn draw_outgoing;
  MorphDrawFn draw_incoming;
};

/**
 * @brief Continuously modulates Mobius parameters to create an evolving
 * warp. Uses multiple frequencies for non-repeating chaos.
 */
class MobiusWarpEvolving : public AnimationBase<MobiusWarpEvolving> {
public:
  /**
   * @brief Constructs a MobiusGenerate animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   */
  MobiusWarpEvolving(MobiusParams &params, float scale = 0.5f,
                     float speed = 0.01f)
      : speed(speed), scale(scale), params(params), base(params),
        seed(hs::random()()) {}

  /// Derive per-channel phase offset from seed (replaces float phases[8]).
  float phase(int i) const {
    uint32_t h = seed ^ (static_cast<uint32_t>(i) * 2654435761u);
    return (h & 0xFFFF) * (100.0f / 65535.0f);
  }

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float time = t * speed;
    float s = scale;

    // Use prime-ish number ratios for frequencies to minimize repetition cycle
    params.get().aRe = base.aRe + sinf(time * 1.0f + phase(0)) * s;
    params.get().aIm = base.aIm + cosf(time * 1.13f + phase(1)) * s;

    params.get().bRe = base.bRe + sinf(time * 1.27f + phase(2)) * s;
    params.get().bIm = base.bIm + cosf(time * 1.39f + phase(3)) * s;

    params.get().cRe = base.cRe + sinf(time * 0.71f + phase(4)) * s;
    params.get().cIm = base.cIm + cosf(time * 0.83f + phase(5)) * s;

    params.get().dRe = base.dRe + sinf(time * 0.97f + phase(6)) * s;
    params.get().dIm = base.dIm + cosf(time * 1.09f + phase(7)) * s;
  }

  float speed;
  float scale;

private:
  std::reference_wrapper<MobiusParams> params;
  MobiusParams base;
  uint32_t seed;
};

/**
 * @brief Adapts an arbitrary behavior function to the Animation system for
 * Palettes. Uses shared state to allow multiple copies (e.g. in Timeline and
 * AnimatedPalette) to sync.
 */

/**
 * @brief Parameters for a ripple wave effect.
 */
struct RippleParams {
  Vector center;         /**< Center point of the ripple source. */
  float amplitude;       /**< Current height of the wave. */
  float phase;           /**< Current phase offset (time). */
  float frequency{20.0}; /**< Spatial frequency of the wave. */
  float decay{5.0};      /**< Spatial decay rate. */
  float thickness{1.0f}; /**< Thickness of the ripple. */

  // Cached dot-product thresholds for fast rejection
  float cos_threshold_min = 1.0f;
  float cos_threshold_max = -1.0f;

  void prepare_thresholds() {
    float hw = thickness * 0.5f;
    if (hw < 0.001f)
      hw = 0.001f;
    float d_min = phase - hw * 2.0f;
    float d_max = phase + hw * 2.0f;
    cos_threshold_min = (d_min >= 0.0f && d_min <= PI_F) ? cosf(d_min) : 1.0f;
    cos_threshold_max = (d_max >= 0.0f && d_max <= PI_F) ? cosf(d_max) : -1.0f;
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
    this->params.get().center = center;
    this->params.get().phase = 0.0f;
    // Start at 0 to prevent 1-frame singularity before first step()
    this->params.get().amplitude = 0.0f;
  }

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);

    // 1. Move the wave (Emanate outward)
    params.get().phase += speed;

    // 2. Lifecycle Management
    float progress = static_cast<float>(t) / duration;
    float envelope = 0.0f;

    if (t < duration) {
      // Ease In (Attack) - fast ramp up over ~10% of duration
      float attack_dur = 0.1f;
      float attack = std::min(progress / attack_dur, 1.0f);
      // Square the attack for a parabolic ease-in (starts very low)
      attack = attack * attack;

      // Fade Out (Decay)
      float decay = 1.0f - progress;

      envelope = attack * decay;
    }

    params.get().amplitude = peak_amplitude * envelope;
  }

private:
  std::reference_wrapper<RippleParams> params;
  float speed;
  float peak_amplitude;
};

/**
 * @brief Parameters for noise transformation.
 */
struct NoiseParams {
  float amplitude = 0.5f;
  float speed = 1.0f;
  float frequency = 0.125f;
  float time = 0.0f;
  float scale = 4.0f;
  mutable FastNoiseLite noise; // Mutable to allow lazy init/updates

  NoiseParams() { noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2); }

  // Helper to ensure frequency is synced
  void sync() const { noise.SetFrequency(frequency); }
};

/**
 * @brief Animates noise parameters by updating time.
 */
class Noise : public AnimationBase<Noise> {
public:
  Noise(NoiseParams &params, int duration = -1)
      : AnimationBase(duration, true), params(params) {}

  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    params.get().time = static_cast<float>(t);
  }

private:
  std::reference_wrapper<NoiseParams> params;
};

} // namespace Animation

/**
 * @brief Structure linking an animation with its starting time.
 * Stores the animation inline to avoid arena allocation (survives compaction).
 */
struct TimelineEvent {
  static constexpr size_t MAX_ANIM_SIZE = 112;

  int start = 0;
  alignas(std::max_align_t) uint8_t storage[MAX_ANIM_SIZE];

  /// Type-erased operation: dst != nullptr → move src into dst and destroy src.
  ///                        dst == nullptr → just destroy src.
  void (*manager)(TimelineEvent &src, TimelineEvent *dst) = nullptr;

  IAnimation *animation() {
    return manager ? reinterpret_cast<IAnimation *>(storage) : nullptr;
  }

  void move_into(TimelineEvent &dst) {
    dst.start = start;
    dst.manager = manager;
    if (manager) {
      manager(*this, &dst);
      manager = nullptr;
    }
  }

  void destroy() {
    if (manager) {
      manager(*this, nullptr);
      manager = nullptr;
    }
  }
};

/**
 * @brief Global storage for the timeline to prevent template instantiation
 * bloat.
 */
static constexpr int TIMELINE_MAX_EVENTS = 64;
extern DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];

/**
 * @brief Manages all active animations and their execution over time.
 */
template <int W, int CAPACITY = 32> class Timeline {
public:
  /**
   * @brief Constructs a Timeline.
   */
  Timeline() { clear(); }

  /**
   * @brief Cleans up remaining animations, invoked on effect destruction.
   */
  ~Timeline() { clear(); }

  void clear() {
    for (int i = 0; i < num_events; ++i) {
      global_timeline_events[i].destroy();
    }
    num_events = 0;
    t = 0;
  }

  /**
   * @brief Adds a new animation event to the timeline.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting.
   * @param animation The animation object.
   * @return Reference to the Timeline object.
   */
  template <typename A> Timeline &add(float in_frames, A animation) {
    static_assert(sizeof(A) <= TimelineEvent::MAX_ANIM_SIZE,
                  "Animation type exceeds TimelineEvent inline storage");
    if (num_events >= MAX_EVENTS) {
      Serial.println("Timeline full, failed to add animation!");
      return *this;
    }
    auto &e = global_timeline_events[num_events++];
    e.start = t + (int)in_frames;
    new (e.storage) A(std::move(animation));
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      A *obj = reinterpret_cast<A *>(src.storage);
      if (dst) {
        new (dst->storage) A(std::move(*obj));
      }
      obj->~A();
    };
    return *this;
  }

  /**
   * @brief Like add(), but returns the typed pointer to the inline-stored
   * animation. Use when you need to hold a reference for later mutation.
   * @warning Pointer is invalidated if the event is moved during compaction.
   */
  template <typename A> A *add_get(float in_frames, A animation) {
    static_assert(sizeof(A) <= TimelineEvent::MAX_ANIM_SIZE,
                  "Animation type exceeds TimelineEvent inline storage");
    if (num_events >= MAX_EVENTS) {
      Serial.println("Timeline full, failed to add animation!");
      return nullptr;
    }
    auto &e = global_timeline_events[num_events++];
    e.start = t + (int)in_frames;
    auto *ptr = new (e.storage) A(std::move(animation));
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      A *obj = reinterpret_cast<A *>(src.storage);
      if (dst) {
        new (dst->storage) A(std::move(*obj));
      }
      obj->~A();
    };
    return ptr;
  }

  /**
   * @brief Advances the timeline by one frame, stepping all active or starting
   * animations.
   * @param canvas The current canvas buffer.
   */
  void step(Canvas &canvas) {
    ++t;

    int write_idx = 0;
    int active_cnt = num_events; // Snapshot count before callbacks
                                 // potentially add more

    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];

      // 1. Check start time
      if (t < e.start) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
        continue;
      }

      // 2. Collapse Orientation (virtual no-op for most types)
      IAnimation *anim = e.animation();
      if (!anim) {
        write_idx++;
        continue;
      }
      anim->collapse_orientation();

      // 3. Step
      anim->step(canvas);

      // 4. Completion & Cleanup
      bool is_done = anim->done();
      bool keep = true;

      if (is_done) {
        bool does_repeat = anim->repeats();
        if (does_repeat) {
          anim->rewind();
          anim->post_callback();
        } else {
          keep = false;
          // Fire callback for non-repeating event before removing
          anim->post_callback();
        }
      }

      if (keep) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
      } else {
        e.destroy();
      }
    }

    // 5. Move new events (added during callbacks) to fill the gap
    int new_vals_count = num_events - active_cnt;
    if (new_vals_count > 0 && write_idx < active_cnt) {
      for (int i = 0; i < new_vals_count; ++i) {
        global_timeline_events[active_cnt + i].move_into(
            global_timeline_events[write_idx + i]);
      }
    }

    num_events = write_idx + new_vals_count;
  }

  inline static int t = 0;          /**< The current global frame count. */
  inline static int num_events = 0; /**< Current number of active events. */

  static constexpr int MAX_EVENTS =
      TIMELINE_MAX_EVENTS; /**< Must match global_timeline_events array size. */
};

/**
 * @brief Manages double-buffered MeshState transitions on a Timeline.
 *
 * Encapsulates the boilerplate of:
 *   - Two persistent MeshState slots registered with PersistentTracker
 *   - Arena compaction when freeing the old slot
 *   - Scoped scratch generation via a user-provided callback
 *   - Crossfade Sprite scheduling with automatic front/back flip
 *
 * Usage:
 *   MeshCarousel<W> carousel;  // in effect members
 *
 *   // In init():
 *   carousel.load([&](Arena& a, Arena& b) {
 *       return solids[0].generate(a, b);
 *   });
 *
 *   // To transition:
 *   carousel.transition(timeline, generate_fn, draw_fn, duration, fade);
 *
 * @tparam W The display width (passed through to Timeline).
 */
template <int W> class MeshCarousel {
public:
  MeshCarousel() {}

  /**
   * @brief Load the initial shape into the front slot (no crossfade).
   * @param generate_fn Given two scratch arenas, returns a PolyMesh.
   * @param classify If true, runs classify_faces_by_topology after loading.
   */
  void load(Fn<PolyMesh(Arena &, Arena &), 48> generate_fn,
            bool classify = true) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
    {
      ScratchScope _a(scratch_arena_a);
      ScratchScope _b(scratch_arena_b);
      PolyMesh mesh = generate_fn(scratch_arena_a, scratch_arena_b);
      slots_[front_].clear();
      MeshOps::compile(mesh, slots_[front_], persistent_arena);
    }
    if (classify) {
      scratch_arena_a.reset();
      scratch_arena_b.reset();
      MeshOps::classify_faces_by_topology(slots_[front_], scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
  }

  /**
   * @brief Trigger a crossfade transition to a new shape.
   *
   * Internally: flips the double-buffer, compacts the arena, generates the
   * new shape into the back slot, then schedules an overlapping Sprite on
   * the timeline. The front index flips eagerly so the draw lambda captures
   * the correct slot, enabling overlapping sprites.
   *
   * @param timeline The timeline to schedule the crossfade onto.
   * @param generate_fn Produces a PolyMesh into scratch arenas.
   * @param draw_fn Called each frame with (Canvas&, float opacity).
   * @param duration Total sprite duration in frames.
   * @param fade_in Fade-in frames (0 = instant appear).
   * @param fade_out Fade-out frames (0 = instant disappear).
   * @param easing Easing function for fades.
   * @param classify If true, runs classify_faces_by_topology.
   */
  template <typename GenerateFn, typename DrawOutgoingFn,
            typename DrawIncomingFn>
  void transition(Timeline<W> &timeline, GenerateFn generate_fn,
                  DrawOutgoingFn draw_outgoing, DrawIncomingFn draw_incoming,
                  int duration, int fade_in = 0, int fade_out = 0,
                  bool classify = true, EasingFn easing = ease_mid) {
    int back = 1 - front_;

    // Free the old back slot and compact
    {
      slots_[back] = MeshState();
      Persist<MeshState> p(slots_[front_], scratch_arena_b, persistent_arena);
      persistent_arena.reset();
    }

    // Generate new shape into back slot via unified generate()
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh mesh = generate_fn(a, b);
      slots_[back].clear();
      MeshOps::compile(mesh, slots_[back], target);
    });

    if (classify) {
      scratch_arena_a.reset();
      scratch_arena_b.reset();
      MeshOps::classify_faces_by_topology(slots_[back], scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Flip front eagerly so the draw lambda captures the correct slot.
    // This enables overlapping sprites where the next transition()
    // is called before the current sprite finishes.
    front_ = back;

    // Schedule the crossfade sprite using the incoming draw function
    timeline.add(0, Animation::Sprite(std::move(draw_incoming), duration,
                                      fade_in, easing, fade_out, easing));
  }

  /// The currently visible (front) mesh.
  const MeshState &current() const { return slots_[front_]; }
  MeshState &current() { return slots_[front_]; }

  /// The incoming (back) mesh — valid during a transition.
  const MeshState &incoming() const { return slots_[1 - front_]; }
  MeshState &incoming() { return slots_[1 - front_]; }

  /// Direct slot access by index (for effects that need both).
  const MeshState &slot(int i) const { return slots_[i]; }
  MeshState &slot(int i) { return slots_[i]; }

  /// Which index is front (for capture in lambdas).
  int front_index() const { return front_; }

  /// Manually set the front index (for effects that manage transitions
  /// themselves).
  void set_front(int idx) { front_ = idx; }

  /// Compact the persistent arena (evacuates tracked MeshStates, reclaims
  /// fragmented space). Call before allocating new persistent data.
  void compact() {
    Persist<MeshState> p0(slots_[0], scratch_arena_a, persistent_arena);
    Persist<MeshState> p1(slots_[1], scratch_arena_a, persistent_arena);
    persistent_arena.reset();
  }

private:
  MeshState slots_[2];
  int front_ = 0;
};

/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&,
 * float t)`.
 */
template <int W, int CAP>
void tween(const Orientation<W, CAP> &o, TweenFn callback) {
  int len = o.length();
  int start = (len > 1) ? 1 : 0;
  for (int i = start; i < len; ++i) {
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 0.0f;
    callback(o.get(i), t);
  }
}

/**
 * @brief Helper to iterate over a VectorTrail's historical frames.
 * @param trail The trail to iterate.
 * @param callback The function to call for each frame: `void(const Vector&,
 * float t)`.
 */
template <int CAPACITY>
void tween(const Animation::VectorTrail<CAPACITY> &trail,
           VectorTweenFn callback) {
  size_t len = trail.length();
  if (len == 0)
    return;

  for (size_t i = 0; i < len; ++i) {
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 0.0f;
    callback(trail.get(i), t);
  }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or
 * Animation::OrientationTrail).
 * @param o The object to iterate.
 * @param callback The function to call for each step: `void(const T&, float
 * t)`.
 */
void deep_tween(const Tweenable auto &trail, TweenFn callback) {
  size_t trail_len = trail.length();
  if (trail_len == 0)
    return;

  for (size_t i = 0; i < trail_len; ++i) {
    const auto &frame = trail.get(i);
    size_t frame_size = frame.length();
    size_t start_j = (i == 0) ? 0 : 1;

    for (size_t j = start_j; j < frame_size; ++j) {
      const auto &q = frame.get(j);
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      float global_t = (static_cast<float>(i) + sub_t) / trail_len;
      callback(q, global_t);
    }
  }
}
#endif // HOLOSPHERE_CORE_ANIMATION_H_
