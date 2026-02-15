/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <variant>
#include <numeric> // for std::iota
#include <array>
#include "3dmath.h"
#include "FastNoiseLite.h"
#include "geometry.h" 
#include "spatial.h"
#include "static_circular_buffer.h"
#include "rotate.h"
#include "color.h"

/**
 * @brief Represents a customizable path.
 * Retains internal buffer for state, but draws to pipeline.
 */
template <int W>
class Path {
public:
  Path() {}

  static void draw(auto& pipeline, Canvas& canvas, const Path& path, auto color) {
    size_t samples = path.points.size();
    for (size_t i = 0; i < samples; ++i) {
      auto v = path.get_point(static_cast<float>(i) / samples);
      auto c = color(v, i / (samples - 1.0f));
      pipeline.plot(canvas, v, c.color, 0, c.alpha);
    }
  }

  Path& append_segment(PlotFn auto plot, float domain, float samples, ScalarFn auto easing) {
    if (!points.is_empty()) points.pop_back();
    for (float t = 0; t <= samples; t++) {
      points.push_back(plot(easing(t / samples) * domain));
    }
    return *this;
  }

  Vector get_point(float t) const {
    if (points.is_empty()) return Vector(0, 0, 0);
    float raw_index = t * (points.size() - 1);
    size_t i = static_cast<size_t>(raw_index);
    float f = raw_index - i;
    if (i >= points.size() - 1) return points.back();
    const Vector& p1 = points[i];
    const Vector& p2 = points[i + 1];
    return p1 * (1.0f - f) + p2 * f;
  }

  size_t num_points() const { return points.size(); }
  void collapse() { if (points.size() > 1) points = { points.back() }; }
  const Points& get_points() const {
    static Points temp;
    temp.clear();
    for (size_t i = 0; i < points.size(); ++i) temp.push_back(points[i]);
    return temp;
  }

private:
  StaticCircularBuffer<Vector, 4096> points;
};

/**
 * @brief Represents a path defined by a single procedural function.
 * Matches interface of Path for use in Motion animations.
 */
template <PlotFn F>
struct ProceduralPath {
  F f;

  ProceduralPath(F path_fn) : f(path_fn) {}

  Vector get_point(float t) const {
    return f(t);
  }
};

#include "easing.h"


namespace Animation {

/**
 * @brief Defines the coordinate space for rotations.
 */
enum class Space {
  World, /**< Rotations are applied relative to the fixed world axes (Pre-multiply). */
  Local  /**< Rotations are applied relative to the object's current axes (Post-multiply). */
};

/**
 * @brief Base class for all animations, providing core timing and state management.
 * @tparam Derived The class inheriting from Animation (Curiously Recurring Template Pattern).
 */
template <typename Derived>
class Base {
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
  Base(int duration, bool repeat) :
    duration(duration == 0 ? 1 : duration),
    repeat(repeat),
    canceled(false),
    post([]() {})
  {
  }

  Base() :
    duration(-1),
    repeat(false),
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
class NullAnimation : public Base<NullAnimation> {
public:
  /**
   * @brief Constructs a NullAnimation.
   */
  NullAnimation() :
    Base(0, false)
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
 * @brief Manages a history of Orientation states.
 * @tparam CAPACITY The maximum number of snapshots to keep.
 */
template <typename OrientationType, int CAPACITY>
class OrientationTrail {
public:
  /**
   * @brief Records a snapshot of the current orientation state.
   * @param source The orientation to copy.
   */
  void record(const OrientationType& source) {
    snapshots.push_back(source);
  }

  /**
   * @brief Gets the number of recorded snapshots.
   */
  size_t length() const {
    return snapshots.size();
  }

  /**
   * @brief Gets a specific snapshot.
   * @param i Index (0 is newest).
   */
  const OrientationType& get(size_t i) const {
    // JS parity: 0 is oldest
    return snapshots[i];
  }

  /**
   * @brief Gets a specific snapshot (mutable).
   * @param i Index (0 is newest).
   */
  OrientationType& get(size_t i) {
    return snapshots[i];
  }

  /**
   * @brief Clears the history.
   */
  void clear() {
    snapshots.clear();
  }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() {
    snapshots.pop();
  }

private:
  StaticCircularBuffer<OrientationType, CAPACITY> snapshots;
};


template <int W>
struct Particle {
  Vector position;
  Vector velocity;
  PaletteVariant palette;
  float life;
  float max_life;
  Orientation<W> orientation;
  
  OrientationTrail<Orientation<W>, 25> history;

  void init(const Vector& p, const Vector& v, float l, const PaletteVariant& pal) {
    position = p;
    velocity = v;
    palette = pal;
    life = l;
    max_life = l;
    orientation.set(Quaternion());
    history.clear();
  }
  
  size_t history_length() const { return history.length(); }
};

template <int W, int CAPACITY>
class ParticleSystem : public Base<ParticleSystem<W, CAPACITY>> {
public:
  std::array<Particle<W>, CAPACITY> pool;
  int active_count = 0;
  
  int time_scale = 1;
  float friction = 0.85f;
  float gravity = 0.001f;
  float resolution_scale = 1.0f;
  int trail_length = 25; // TODO: this doesn't work as particles have static capacity

  struct Attractor {
    Vector position;
    float strength;
    float kill_radius;
    float event_horizon;
  };

  using EmitterFn = std::function<void(ParticleSystem&)>;
  
  StaticCircularBuffer<Attractor, 32> attractors;
  StaticCircularBuffer<EmitterFn, 32> emitters;

  ParticleSystem(float friction = 0.85f, float gravity = 0.001f, int trail_length = 25) : 
    Base<ParticleSystem<W, CAPACITY>>(-1, false), 
    friction(friction), 
    gravity(gravity),
    trail_length(trail_length)
  {
  }
  
  void reset(float friction = 0.85f, float gravity = 0.001f, int trail_length = 25) {
    this->friction = friction;
    this->gravity = gravity;
    this->trail_length = trail_length;
    
    active_count = 0;
    attractors.clear();
    emitters.clear();
  }

  void add_emitter(EmitterFn fn) {
    emitters.push_back(fn);
  }
  
  void add_attractor(const Vector& pos, float str, float kill, float horizon) {
    attractors.push_back({ pos, str, kill, horizon });
  }

  void spawn(const Vector& pos, const Vector& vel, const Color4& col, float life = 600) {
    if (active_count < CAPACITY) {
      pool[active_count++].init(pos, vel, life, SolidColorPalette(col));
    }
  }

  void spawn(const Vector& pos, const Vector& vel, const PaletteVariant& pal, float life = 600) {
    if (active_count < CAPACITY) {
      pool[active_count++].init(pos, vel, life, pal);
    }
  }

  void step(Canvas& canvas) override {
    Base<ParticleSystem<W, CAPACITY>>::step(canvas);

    for (int k = 0; k < time_scale; ++k) {
        // Emitters
        for (size_t i = 0; i < emitters.size(); ++i) {
          emitters[i](*this);
        }

        float max_delta = (2 * PI_F) / W / resolution_scale;

        // Physics
        for (int i = 0; i < active_count; ++i) {
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

  bool step_particle(Particle<W>& p, float max_delta) {
    p.life--;
    bool active = p.life > 0;

    // Physics
    if (active) {
        Vector pos = p.orientation.orient(p.position);

        for (size_t k = 0; k < attractors.size(); ++k) {
           const auto& attr = attractors[k];
           float dist_sq = distance_squared(pos, attr.position);
           
           if (dist_sq < attr.kill_radius * attr.kill_radius) {
               active = false; 
               break; // Killed
           }
           
           if (dist_sq > 0.0000001f) {
              if (dist_sq < attr.event_horizon * attr.event_horizon) {
                 // Steer into center
                 Vector torque = (attr.position - pos).normalize();
                 float speed = p.velocity.magnitude();
                 p.velocity = torque * speed;
              } else {
                 // Gravity
                 float force = (gravity * attr.strength) / dist_sq;
                 Vector torque = cross(pos, attr.position).normalize() * force;
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
                Vector axis = cross(pos, p.velocity).normalize();
                Quaternion dq = make_rotation(axis, speed);
                // JS Parity: Accumulate rotation (Absolute) instead of pushing Delta
                p.orientation.push(dq * p.orientation.get());
                p.velocity = rotate(p.velocity, dq);
            }
        }
    }
    
    // History Management
    if (active) {
        p.history.record(p.orientation);
    } else {
        if (p.history.length() > 0) {
            p.history.expire();
        }
    }
    
    p.orientation.collapse();
    
    return !active && p.history.length() == 0;
  }
};

/**
 * @brief An animation that triggers a callback after a random delay.
 */
class RandomTimer : public Base<RandomTimer> {
public:

  /**
   * @brief Constructs a RandomTimer.
   * @param min Minimum delay in frames.
   * @param max Maximum delay in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  RandomTimer(int min, int max, TimerFn auto f, bool repeat = false) :
    Base(-1, repeat),
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
    Base::step(canvas);
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
class PeriodicTimer : public Base<PeriodicTimer> {
public:
  /**
   * @brief Constructs a PeriodicTimer.
   * @param period The interval between calls, in frames.
   * @param f The function to call when the timer elapses.
   * @param repeat If true, the timer resets after calling the function.
   */
  PeriodicTimer(int period, TimerFn auto f, bool repeat = false) :
    Base(-1, repeat),
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
    Base::step(canvas);
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
class Transition : public Base<Transition> {
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
    Base(duration, repeat),
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
    Base::step(canvas);
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
class Mutation : public Base<Mutation> {
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
    Base(duration, repeat),
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
    Base::step(canvas);
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
class Sprite : public Base<Sprite> {
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
  Sprite(std::function<void(Canvas&, float)> draw_fn, int duration,
    int fade_in_duration = 0, std::function<float(float)> fade_in_easing_fn = ease_mid,
    int fade_out_duration = 0, std::function<float(float)> fade_out_easing_fn = ease_mid) :
    Base(duration, false),
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
    : Base(std::move(other)),
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

    Base::operator=(std::move(other));
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
    Base::step(canvas);
    if (!fade_in.done()) {
      fade_in.step(canvas);
    }
    else if (duration >= 0 && fade_out_duration > 0 && t >= (duration - fade_out_duration)) {
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
 * @tparam W The width of the LED display (used for calculating maximum rotation step).
 */
template <int W>
class Motion : public Base<Motion<W>> {
public:

  /**
   * @brief Constructs a Motion animation.
   * @param orientation The Orientation object to update.
   * @param path_obj The path to follow (must have get_point(float t)).
   * @param duration The duration in frames.
   * @param repeat If true, the motion repeats.
   */
  template <typename P>
  Motion(Orientation<W>& orientation, const P& path_obj, int duration, bool repeat = false, Space space = Space::World) :
    Base<Motion<W>>(duration, repeat),
    orientation(orientation),
    path_fn([&path_obj](float t) { return path_obj.get_point(t); }),
    space(space)
  {
  }

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W>& get_orientation() const { return orientation.get(); }

  /**
   * @brief Steps the animation, calculates intermediate rotation steps along the path,
   * and pushes them to the Orientation.
   */
  void step(Canvas& canvas) {
    Base<Motion<W>>::step(canvas);
    float t_prev = static_cast<float>(this->t - 1);
    Vector current_v = path_fn(t_prev / this->duration);
    float t_curr = static_cast<float>(this->t);
    Vector target_v = path_fn(t_curr / this->duration);
    float total_angle = angle_between(current_v, target_v);
    int num_steps = static_cast<int>(std::ceil(std::max(1.0f, total_angle / MAX_ANGLE)));

    // Ensure sufficient resolution
    orientation.get().upsample(num_steps + 1);
    int len = orientation.get().length();

    // lambda to apply rotation based on space
    auto apply_rotation = [&](Quaternion& target, const Quaternion& source) {
      if (space == Space::Local) {
        target = target * source;
      }
      else {
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
        Vector step_axis = cross(prev_v, next_v).normalize();
        Quaternion q_step = make_rotation(step_axis, step_angle);
        apply_rotation(accumulated_q, q_step);
      }

      // Apply the accumulated rotation to the existing orientation frame
      Quaternion& current_q = orientation.get().at(i);
      apply_rotation(current_q, accumulated_q);
      current_q.normalize();
      
      prev_v = next_v;
    }
  }

private:

  static constexpr float MAX_ANGLE = 2 * PI_F / W; /**< Maximum rotation angle per step to ensure smoothness. */
  std::reference_wrapper<Orientation<W>> orientation; /**< Reference to the Orientation state. */
  std::function<Vector(float)> path_fn; /**< Function to retrieve path points. */
  Space space; /**< The coordinate space for rotation. */
};

/**
 * @brief An animation that applies a fixed, time-eased rotation.
 * @tparam W The width of the LED display (used for calculating maximum rotation step).
 */
template <int W>
class Rotation : public Base<Rotation<W>> {
public:

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
  Rotation(Orientation<W>& orientation, const Vector& axis, float angle, int duration, ScalarFn auto easing_fn, bool repeat = false, Space space = Space::World) :
    Base<Rotation<W>>(duration, repeat),
    orientation(orientation),
    axis(axis),
    total_angle(angle),
    easing_fn(easing_fn),
    last_angle(0),
    space(space)
  {
  }

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W>& get_orientation() const { return orientation.get(); }

  /**
   * @brief Steps the animation, calculates the incremental rotation delta, and pushes it to the Orientation.
   */
  void step(Canvas& canvas) {
    if (this->t == 0) {
      last_angle = 0;
    }
    Base<Rotation<W>>::step(canvas);
    float target_angle = easing_fn(static_cast<float>(this->t) / this->duration) * total_angle;
    float delta = target_angle - last_angle;
    
    if (std::abs(delta) < TOLERANCE) {
      last_angle = target_angle;
      return;
    }

    int num_steps = 1 + static_cast<int>(std::ceil(std::abs(delta) / MAX_ANGLE));
    orientation.get().upsample(num_steps + 1);
    int len = orientation.get().length();
    
    float step_angle = delta / (len - 1);
    Quaternion accumulated_q;
    
    // lambda to apply rotation based on space
    auto apply_rotation = [&](Quaternion& target, const Quaternion& source) {
      if (space == Space::Local) {
        target = target * source;
      }
      else {
        target = source * target;
      }
    };

    for (int i = 1; i < len; ++i) {
      float angle = step_angle * i;
      Quaternion q = make_rotation(axis, angle);
      
      Quaternion& current_q = orientation.get().at(i);
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
  static void animate(Canvas& canvas, Orientation<W>& orientation, const Vector& axis, float_t angle, ScalarFn auto easing_fn, Space space = Space::World) {
    Rotation<W> r(orientation, axis, angle, 1, easing_fn, false, space);
    r.step(canvas);
  }

private:

  static constexpr float MAX_ANGLE = 2 * PI_F / W; /**< Maximum rotation angle per step to ensure smoothness. */
  std::reference_wrapper<Orientation<W>> orientation; /**< Reference to the Orientation state. */
  Vector axis; /**< The axis of rotation. */
  float total_angle; /**< The total angle to sweep. */
  std::function<float(float)> easing_fn; /**< Easing curve. */
  float last_angle; /**< The angle reached in the previous frame. */
  Space space; /**< The coordinate space for rotation. */
};

/**
 * @brief An animation that simulates a particle/camera performing a random walk across the sphere's surface.
 * @details Uses Perlin noise to create continuous, turbulent pivoting motion.
 * @tparam W The width of the LED display.
 */
template<int W>
class RandomWalk : public Base<RandomWalk<W>> {
public:
  struct Options {
    float speed = 0.02f;          /**< Movement speed per frame. */
    float pivot_strength = 0.1f;  /**< Strength of the direction change (noise amplitude). */
    float noise_scale = 0.02f;    /**< Frequency of the Perlin noise. */
    int seed = 0;                 /**< Random seed (0 for random). */
    Space space = Space::World;   /**< Coordinate space for movement. */

    static Options Languid() { return { 0.02f, 0.1f, 0.02f }; }
    static Options Energetic() { return { 0.05f, 0.4f, 0.08f }; }
  };

  /**
   * @brief Constructs a RandomWalk animation.
   * @param orientation The Orientation object to update.
   * @param v_start The starting direction vector.
   * @param options Configuration options.
   */
  RandomWalk(Orientation<W>& orientation, const Vector& v_start, Options options = Options()) :
    Base<RandomWalk<W>>(-1, false),
    orientation(orientation),
    v(Vector(v_start).normalize()),
    options(options)
  {
    Vector u = X_AXIS;
    if (std::abs(dot(v, u)) > 0.99f) {
      u = Y_AXIS;
    }
    direction = cross(v, u).normalize();
    noiseGenerator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noiseGenerator.SetFrequency(options.noise_scale);
    if (options.seed == 0) {
      noiseGenerator.SetSeed(std::rand());
    } else {
      noiseGenerator.SetSeed(options.seed);
    }
  }

  /**
   * @brief Access the associated Orientation.
   */
  Orientation<W>& get_orientation() const { return orientation.get(); }

  /**
   * @brief Steps the walk: pivots direction based on noise, then rotates the view along the calculated axis.
   */
  void step(Canvas& canvas) override {
    Base<RandomWalk<W>>::step(canvas);
    float pivotAngle = noiseGenerator.GetNoise(static_cast<float>(this->t), 0.0f) * options.pivot_strength;
    direction = rotate(direction, make_rotation(v, pivotAngle)).normalize();
    Vector walk_axis = cross(v, direction).normalize();
    v = rotate(v, make_rotation(walk_axis, options.speed)).normalize();
    direction = rotate(direction, make_rotation(walk_axis, options.speed)).normalize();
    Rotation<W>::animate(canvas, orientation, walk_axis, options.speed, ease_mid, options.space);
  }

private:
  std::reference_wrapper<Orientation<W>> orientation; /**< Reference to the global Orientation state. */
  Vector v; /**< Current forward direction vector. */
  Vector direction; /**< Current pivoting direction (orthogonal to v). */
  Options options; /**< Configuration options. */
  FastNoiseLite noiseGenerator; /**< The noise generator instance. */
};

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief An animation that smoothly interpolates a GenerativePalette toward a target palette.
 */
class ColorWipe : public Base<ColorWipe> {
public:
  /**
   * @brief Constructs a ColorWipe animation.
   * @param from_palette The GenerativePalette to start modifying (will be used for `cur_palette`).
   * @param to_palette The GenerativePalette instance to interpolate toward.
   * @param duration The duration in frames.
   * @param easing_fn The easing function.
   */
  ColorWipe(GenerativePalette& from_palette, const GenerativePalette& to_palette, int duration, ScalarFn auto easing_fn) :
    Base(duration, false),
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
    Base::step(canvas);
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
class MobiusFlow : public Base<MobiusFlow> {
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
    Base(duration, repeat),
    params(params),
    num_rings(num_rings),
    num_lines(num_lines)
  {
  }

  /**
   * @brief Steps the animation, updating params a and d.
   */
  void step(Canvas& canvas) {
    Base::step(canvas);
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
class MobiusWarp : public Base<MobiusWarp> {
public:
  /**
   * @brief Constructs a MobiusWarp animation.
   * @param params Reference to the MobiusParams to animate.
   * @param scale The magnitude of the warp effect.
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   * @param easing The easing function to use (default: ease_in_out_sin).
   */
  MobiusWarp(MobiusParams& params, float scale, int duration, bool repeat = true, std::function<float(float)> easing = ease_in_out_sin) :
    Base(duration, repeat),
    params(params),
    scale(scale),
    easing(easing)
  {
  }

  /**
   * @brief Steps the animation, updating param b.
   */
  void step(Canvas& canvas) {
    Base::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(std::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    params.get().bRe = scale * (cosf(angle) - 1.0f);
    params.get().bIm = scale * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params;
public:
  float scale;
  std::function<float(float)> easing;
};




/*
 * @brief Animates a transition between two meshes using pre-allocated buffers.
 */
// NEW: Morphing Support Structures
struct MorphPath {
    Vector start, end, axis;
    float angle;
};

struct MorphBuffer {
    MeshState source;
    MeshState dest;
    // Paths for Source -> Dest
    std::vector<MorphPath> source_paths;
    // Paths for Dest -> Source
    std::vector<MorphPath> dest_paths;

    void load_dest(const PolyMesh& poly) {
        // Helper to load destination from PolyMesh if needed
        // For now, assuming manual population or via MeshMorph init
    }
};

/*
 * @brief Animates a transition between two meshes using pre-allocated buffers.
 * Implements dual-morph logic: Source deforms to Dest, while Dest deforms from Source.
 */
class MeshMorph : public Base<MeshMorph> {
public:
  using State = MeshState;

  /**
   * @brief Constructs a MeshMorph animation.
   * @param output_source Pointer to the active source mesh state (gets deformed).
   * @param output_dest Pointer to the active dest mesh state (gets deformed).
   * @param buffer Pointer to shared MorphBuffer.
   * @param source The starting geometry of source.
   * @param dest The target geometry of source (and starting geometry of dest).
   * @param duration Duration in frames.
   */
  MeshMorph(MeshState* output_source, MeshState* output_dest, MorphBuffer* buffer, const MeshState& source, const MeshState& dest, int duration, bool repeat = false, ScalarFn auto easing_fn = ease_in_out_sin)
    : Base(duration, repeat), output_source(output_source), output_dest(output_dest), buffer(buffer), easing_fn(easing_fn)
  {
      if (buffer) {
          // Deep copy initial states
          buffer->source = source; 
          buffer->dest = dest;     
          
          // Reset outputs to starting states
          if (output_source) *output_source = buffer->source;
          if (output_dest) *output_dest = buffer->dest;

          init();
      }
  }
  
  // Single output constructor compatibility
  MeshMorph(MeshState* output_source, MorphBuffer* buffer, const MeshState& source, const MeshState& dest, int duration, bool repeat = false, ScalarFn auto easing_fn = ease_in_out_sin)
    : MeshMorph(output_source, nullptr, buffer, source, dest, duration, repeat, easing_fn) {}

  void init() {
      if (!buffer) return;
      
      // 1. Source -> Dest Paths
      size_t src_count = buffer->source.vertices.size();
      buffer->source_paths.resize(src_count);
      for(size_t i=0; i<src_count; ++i) {
         Vector s = buffer->source.vertices[i];
         Vector t = project_to_mesh(s, buffer->dest);
         
         float ang = angle_between(s, t);
         Vector ax = cross(s, t);
         if (dot(ax, ax) > 0.00001f) ax = ax.normalize();
         else ax = Vector(1,0,0);

         buffer->source_paths[i] = { s, t, ax, ang };
      }

      // 2. Dest -> Source Paths (Reverse Morph)
      size_t dest_count = buffer->dest.vertices.size();
      buffer->dest_paths.resize(dest_count);
      for(size_t i=0; i<dest_count; ++i) {
         Vector t = buffer->dest.vertices[i]; 
         Vector s = project_to_mesh(t, buffer->source);
         
         float ang = angle_between(s, t);
         Vector ax = cross(s, t);
         if (dot(ax, ax) > 0.00001f) ax = ax.normalize();
         else ax = Vector(1,0,0);

         buffer->dest_paths[i] = { s, t, ax, ang };
      }
  }
  
  void step(Canvas& canvas) override {
      Base::step(canvas);
      if (!buffer) return;

      float progress = std::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
      float alpha = easing_fn(progress);
      
      // Update Source Output (Deform Source -> Dest)
      if (output_source) {
          size_t count = buffer->source_paths.size();
          if (output_source->vertices.size() != count) output_source->vertices.resize(count);
          
          for(size_t i=0; i<count; ++i) {
             const auto& path = buffer->source_paths[i];
             if (path.angle > 0.00001f) {
                Quaternion R = make_rotation(path.axis, path.angle * alpha);
                output_source->vertices[i] = rotate(path.start, R);
             } else {
                 output_source->vertices[i] = (path.start * (1-alpha) + path.end * alpha).normalize();
             }
          }
      }

      // Update Dest Output (Unwarp Dest from Source)
      if (output_dest) {
          size_t count = buffer->dest_paths.size();
          if (output_dest->vertices.size() != count) output_dest->vertices.resize(count);

          for(size_t i=0; i<count; ++i) {
             const auto& path = buffer->dest_paths[i];
             if (path.angle > 0.00001f) {
                // Dest unwinds from Source: alpha goes 0->1.
                // At alpha=0 (start), using make_rotation(axis, angle * 0) = Identity.
                // rotate(start, I) = start. (start is 's' on Source).
                // This matches Source at t=0. Correct.
                Quaternion R = make_rotation(path.axis, path.angle * alpha);
                output_dest->vertices[i] = rotate(path.start, R);
             } else {
                 output_dest->vertices[i] = (path.start * (1-alpha) + path.end * alpha).normalize();
             }
          }
      }
      
      // Final Swap (Commit Topology to Source Output)
      if (t >= duration && output_source) {
           *output_source = buffer->dest;
      }
  }

private:
  MeshState* output_source;
  MeshState* output_dest;
  MorphBuffer* buffer;
  std::function<float(float)> easing_fn;
};
  



/**
 * @brief Continuously modulates Mobius parameters to create an evolving warp.
 * Uses multiple frequencies for non-repeating chaos.
 */
class MobiusGenerate : public Base<MobiusGenerate> {
public:
  /**
   * @brief Constructs a MobiusGenerate animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   */
  MobiusGenerate(MobiusParams& params, float scale = 0.5f, float speed = 0.01f) :
    Base(-1, true),
    params(params),
    scale(scale),
    speed(speed)
  {
    // Capture initial state as base
    base = params;

    // Random phase offsets
    phases[0] = hs::rand_f() * 100.0f; // aRe
    phases[1] = hs::rand_f() * 100.0f; // aIm
    phases[2] = hs::rand_f() * 100.0f; // bRe
    phases[3] = hs::rand_f() * 100.0f; // bIm
    phases[4] = hs::rand_f() * 100.0f; // cRe
    phases[5] = hs::rand_f() * 100.0f; // cIm
    phases[6] = hs::rand_f() * 100.0f; // dRe
    phases[7] = hs::rand_f() * 100.0f; // dIm
  }

  void step(Canvas& canvas) {
    Base::step(canvas);
    float time = t * speed;
    float s = scale;

    // Use prime-ish number ratios for frequencies to minimize repetition cycle
    params.get().aRe = base.aRe + sinf(time * 1.0f + phases[0]) * s;
    params.get().aIm = base.aIm + cosf(time * 1.13f + phases[1]) * s;

    params.get().bRe = base.bRe + sinf(time * 1.27f + phases[2]) * s;
    params.get().bIm = base.bIm + cosf(time * 1.39f + phases[3]) * s;

    params.get().cRe = base.cRe + sinf(time * 0.71f + phases[4]) * s;
    params.get().cIm = base.cIm + cosf(time * 0.83f + phases[5]) * s;

    params.get().dRe = base.dRe + sinf(time * 0.97f + phases[6]) * s;
    params.get().dIm = base.dIm + cosf(time * 1.09f + phases[7]) * s;
  }

private:
  std::reference_wrapper<MobiusParams> params;
  float scale;
  float speed;
  MobiusParams base;
  std::array<float, 8> phases;
};

/**
 * @brief Adapts an arbitrary behavior function to the Animation system for Palettes.
 * Uses shared state to allow multiple copies (e.g. in Timeline and AnimatedPalette) to sync.
 */


/**
 * @brief Adapts a PaletteModifier to the Animation system.
 * Holds a reference to the modifier, allowing the Timeline to drive it.
 */
class PaletteAnimation : public Base<PaletteAnimation> {
public:
    /**
     * @param modifier Reference to the stateful modifier object.
     * @param duration Duration in frames (-1 for infinite).
     * @param repeat Whether to repeat.
     */
    PaletteAnimation(PaletteModifier& modifier, int duration = -1, bool repeat = true) 
        : Base(duration, repeat), modifier(modifier) 
    {
    }

    /**
     * @brief Steps the modifier.
     */
    void step(Canvas& canvas) {
        Base::step(canvas); // Updates base 't' (optional usage)
        modifier.get().step();   // Updates the modifier's internal state
    }

private:
    std::reference_wrapper<PaletteModifier> modifier;
};

} // namespace Animation

template <int W>
using AnimationVariant = std::variant<
  Animation::NullAnimation,
  Animation::Sprite,
  Animation::Transition,
  Animation::Mutation,
  Animation::RandomTimer,
  Animation::PeriodicTimer,
  Animation::Rotation<W>,
  Animation::Motion<W>,
  Animation::RandomWalk<W>,
  Animation::ColorWipe,
  Animation::MobiusFlow,
  Animation::MobiusGenerate,
  Animation::MobiusWarp,
  Animation::MeshMorph,
  Animation::PaletteAnimation
>;

/**
 * @brief Structure linking an animation variant with its starting time.
 */
template <int W>
struct TimelineEvent {
  int start; /**< The global frame count at which the animation should begin. */
  AnimationVariant<W> animation; /**< The actual animation object. */
};

/**
 * @brief Manages all active animations and their execution over time.
 */
template <int W, int CAPACITY = 32>
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
    TimelineEvent<W>& e = events[num_events++];
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

    // Prep Animations
    std::vector<Orientation<W>*> touched;
    touched.reserve(num_events); // Pre-allocate max size
    for (int i = 0; i < num_events; ++i) {
      if (t < events[i].start) {
        continue;
      }
      std::visit([&](auto& a) { 
        if constexpr (requires { a.get_orientation(); }) {
          auto& o = a.get_orientation();
          bool already_touched = false;
          for (auto* touched_o : touched) {
            if (touched_o == &o) {
              already_touched = true;
              break;
            }
          }
          if (!already_touched) {
            o.collapse();
            touched.push_back(&o);
          }
        }
      }, events[i].animation);
    }

    // Step Animations
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

    // 2. Fire Callbacks for Finished Events
    // Capture original count because callbacks might add new events
    int check_count = num_events;
    for (int i = 0; i < check_count; ++i) {
        if (t < events[i].start) continue;
        
        auto& animation = events[i].animation;
        // If it's done (and wasn't rewound, meaning it's non-repeating)
        if (std::visit([&](auto& a) { return a.done(); }, animation)) {
             std::visit([&](auto& a) { a.post_callback(); }, animation);
        }
    }

    // 3. Remove Finished Events
    // Now we can safely remove any event that is 'done'.
    // New events added by callbacks are not 'done', so they survive.
    auto new_logical_end = std::remove_if(
      events.begin(), events.begin() + num_events,
      [this](auto& event) {
        // Remove if done
        // Note: Repeating events were already rewound in step 1, so they are not done.
        return std::visit([](auto& a) { return a.done(); }, event.animation);
      }
    );
    num_events = std::distance(events.begin(), new_logical_end);
  }

  int t = 0; /**< The current global frame count. */

private:

  static constexpr int MAX_EVENTS = CAPACITY; /**< Maximum number of concurrent animation events. */
  std::array<TimelineEvent<W>, MAX_EVENTS> events; /**< Storage for all animation events. */
  int num_events; /**< Current number of active events. */
};


/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&, float t)`.
 */
template <int W, typename F>
void tween(const Orientation<W>& o, F callback) {
    int len = o.length();
    int start = (len > 1) ? 1 : 0;
    for (int i = start; i < len; ++i) {
        float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 0.0f;
        callback(o.get(i), t);
    }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or Animation::OrientationTrail).
 * @param o The object to iterate.
 * @param callback The function to call for each step: `void(const T&, float t)`.
 */
template <typename T, typename F>
void deep_tween(const T& trail, F callback) {
    size_t trail_len = trail.length();
    if (trail_len == 0) return;

    for (size_t i = 0; i < trail_len; ++i) {
        const auto& frame = trail.get(i);
        size_t frame_size = frame.length();
        size_t start_j = (i == 0) ? 0 : 1;
        
        for (size_t j = start_j; j < frame_size; ++j) {
            const auto& q = frame.get(j);
            float sub_t = (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
            float global_t = (static_cast<float>(i) + sub_t) / trail_len;
            callback(q, global_t);
        }
    }
}
