/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>

#include <numeric> // for std::iota
#include <array>
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
    // points is a fixed-capacity ring buffer: overflowing it silently
    // overwrites the oldest points and corrupts the path. A sizing bug should
    // trap on the bench, not degrade. (Cold path — path construction.)
    // samples >= 1 also keeps the t / samples below from dividing by zero:
    // easing(0/0) = NaN would silently append a garbage point. samples is an
    // integer count, so the final t == samples lands exactly on easing(1.0).
    HS_CHECK(samples >= 1);
    // Capacity check accounts for the pop_back below: a non-empty path drops its
    // last point before appending samples + 1 new ones, so the final size is
    // (size - 1) + samples + 1. Counting the popped point as still-present would
    // reject a final segment that actually fits (over-conservative by one).
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
    // Clamp like every other interpolator here: an out-of-[0,1] t would make
    // raw_index negative, and static_cast<size_t> of a negative float is UB
    // (t > 1 is caught by the i >= size-1 guard below, but t < 0 is not).
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
    // Clear in place and re-push the single survivor rather than assigning a
    // fresh buffer: a StaticCircularBuffer<Vector, RESOLUTION> temporary is
    // ~12.3 KB at the default RESOLUTION=1024, over the WASM build's 8 KB stack.
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
    return canceled || (duration >= 0 && t >= duration);
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
  bool is_finite() const { return duration >= 0; }
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
        canceled(false) {}

  /**
   * @brief Default constructor: an indefinite, non-repeating animation.
   */
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
 * @tparam OrientationType The orientation snapshot type stored in the trail.
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
   * @return Count of live snapshots in the trail.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Gets a specific snapshot.
   * @param i Index into the history: 0 is the OLDEST snapshot, length()-1 the
   *          newest. (record() appends the newest at the end of the underlying
   *          ring buffer, whose operator[](0) is the oldest live element.)
   *          Matches the JS simulator's trail ordering.
   * @return Const reference to the requested snapshot.
   */
  const OrientationType &get(size_t i) const { return snapshots[i]; }

  /**
   * @brief Gets a specific snapshot (mutable). 0 is oldest, length()-1 newest.
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Mutable reference to the requested snapshot.
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
  /**
   * @brief Records a world-space position snapshot.
   * @param source The position to copy into the trail.
   */
  void record(const Vector &source) { snapshots.push_back(source); }

  /**
   * @brief Gets the number of recorded snapshots.
   * @return Count of live positions in the trail.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Gets a specific position snapshot.
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Const reference to the requested position.
   */
  const Vector &get(size_t i) const { return snapshots[i]; }

  /**
   * @brief Gets a specific position snapshot (mutable).
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Mutable reference to the requested position.
   */
  Vector &get(size_t i) { return snapshots[i]; }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop(); }

private:
  StaticCircularBuffer<Vector, CAPACITY> snapshots;
};

/**
 * @brief Represents a single particle in a system.
 * @tparam TRAIL_LEN Maximum number of trail positions retained.
 */
template <int TRAIL_LEN = 8> struct Particle {
  Vector position;         /**< Current 3D position. */
  Vector velocity;         /**< Current velocity vector. */
  uint16_t color_seed = 0; /**< Hue seed for palette offset. */
  uint16_t life = 0;       /**< Remaining life (frames or arbitrary units). */

  VectorTrail<TRAIL_LEN> history; /**< Trail of world-space positions. */

  /**
   * @brief (Re)initializes the particle and clears its trail.
   * @param p Initial world-space position.
   * @param v Initial velocity vector.
   * @param seed Hue seed for palette offset.
   * @param l Life in frames; clamped into uint16_t range.
   */
  void init(const Vector &p, const Vector &v, uint16_t seed, float l) {
    position = p;
    velocity = v;
    color_seed = seed;
    life = static_cast<uint16_t>(std::min(l, 65535.0f));
    history.clear();
  }

  /**
   * @brief Gets the current trail length.
   * @return Number of recorded trail positions.
   */
  size_t history_length() const { return history.length(); }
};

/**
 * @brief A physics-based particle system with emitters and attractors.
 * @tparam W Width of the display (for Orientation).
 * @tparam CAPACITY Maximum number of particles in the pool.
 * @tparam TRAIL_LEN Trail length per particle.
 * @tparam EMITTER_CAP Maximum number of emitters.
 * @tparam ATTRACTOR_CAP Maximum number of attractors.
 */
template <int W, int CAPACITY, int TRAIL_LEN = 8, int EMITTER_CAP = 8,
          int ATTRACTOR_CAP = 8>
class ParticleSystem
    : public AnimationBase<
          ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP>> {
public:
  ArenaVector<Particle<TRAIL_LEN>> pool; /**< Backing pool of particles. */
  uint16_t active_count = 0; /**< Number of live particles in the pool prefix. */

  float friction = 0.85f;  /**< Per-frame velocity damping factor. */
  float gravity = 0.001f;  /**< Base gravitational constant for attractors. */
  uint16_t max_life = 600; /**< Default particle lifetime in frames. */

  /**
   * @brief A point attractor influencing nearby particles.
   */
  struct Attractor {
    Vector position;     /**< World-space center of the attractor. */
    float strength;      /**< Attractive force multiplier. */
    float kill_radius;   /**< Radius within which particles are killed. */
    float event_horizon; /**< Radius within which steering becomes radial. */
  };

  using EmitterFn = Fn<void(ParticleSystem &), 32>; /**< Per-frame emitter functor. */

  ArenaVector<Attractor> attractors; /**< Active attractors. */
  ArenaVector<EmitterFn> emitters;   /**< Active emitters. */

  /**
   * @brief Constructs an indefinite, non-repeating particle system.
   */
  ParticleSystem()
      : AnimationBase<
            ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP>>(
            -1, false) {}

  /**
   * @brief Binds the pool, attractor, and emitter vectors and pre-constructs
   * inactive particles.
   * @param arena Arena providing backing storage for all vectors.
   * @param friction Per-frame velocity damping factor.
   * @param gravity Base gravitational constant for attractors.
   * @param max_life Default particle lifetime in frames; clamped into uint16_t.
   * @details Pre-constructs CAPACITY inactive particles in the pool.
   */
  void init(Arena &arena, float friction = 0.85f, float gravity = 0.001f,
            float max_life = 600.0f) {
    this->friction = friction;
    this->gravity = gravity;
    // Clamp before narrowing: a float outside [0, 65535] is UB on conversion to
    // uint16_t. Particle::init's std::min sits on the far side of this store
    // (its caller already passes the narrowed member), so the guard has to live
    // here. Current callers are in range; this keeps a stray value from spawning
    // particles with a garbage lifetime.
    this->max_life =
        static_cast<uint16_t>(std::min(std::max(max_life, 0.0f), 65535.0f));
    active_count = 0;
    pool.bind(arena, CAPACITY);
    // No is_bound() guard: bind() routes through Arena::allocate, which traps on
    // OOM and never leaves the vector unbound for CAPACITY > 0.
    for (size_t i = 0; i < CAPACITY; ++i) {
      pool.emplace_back();
    }
    attractors.bind(arena, ATTRACTOR_CAP);
    emitters.bind(arena, EMITTER_CAP);
  }

  /**
   * @brief Registers a per-frame emitter functor.
   * @param fn The emitter to add.
   */
  void add_emitter(EmitterFn fn) { emitters.push_back(fn); }

  /**
   * @brief Adds a point attractor to the system.
   * @param pos World-space center of the attractor.
   * @param str Attractive force multiplier.
   * @param kill Kill radius (particles within it die).
   * @param horizon Event-horizon radius (steering becomes radial within it).
   */
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

  /**
   * @brief Advances the particle system by one frame.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details Runs all emitters, then advances each active particle's physics,
   * swap-removing any that died (compacting the live prefix of the pool).
   */
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
  /**
   * @brief Advances one particle on the sphere surface for one frame.
   * @param p The particle to advance (modified in place).
   * @param max_delta Per-frame surface-rotation cap (radians), one display
   * column wide; the move below clamps to it so a fast particle never jumps more
   * than one column per frame (trail/motion-blur aliasing).
   * @return True once the particle is dead AND its trail has fully drained (so
   * the caller can remove it).
   * @details Ages it, applies attractor gravity/steering and drag, rotates
   * position+velocity along the surface, and updates its trail.
   */
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
            // Gravity. pos and the attractor can be ~antipodal, leaving the
            // cross-product axis undefined; guard the normalize.
            float force = (gravity * attr.strength) / dist_sq;
            Vector torque =
                normalized_or(cross(pos, attr.position), Vector(1, 0, 0)) * force;
            p.velocity += cross(torque, pos);
          }
        }
      }

      if (active) {
        // Drag
        p.velocity *= friction;

        // Move. The surface-rotation axis is cross(pos, velocity); it vanishes
        // when the velocity is purely radial (parallel to pos) even though
        // speed is non-zero — e.g. once a non-tangent external force (the flow-
        // field noise) has pushed the velocity off the tangent plane. A radial
        // velocity produces no motion along the sphere, so skip the step rather
        // than normalize a zero-length axis (mirrors the normalized_or guard on
        // the gravity branch above).
        float speed = p.velocity.magnitude();
        Vector axis = cross(pos, p.velocity);
        if (speed > 0.000001f && axis.magnitude() > 0.000001f) {
          // Cap the surface advance at one display column/frame so a high-speed
          // particle can't skip columns (trail/motion-blur aliasing). Velocity
          // keeps its full magnitude; only the per-frame step is clamped.
          speed = std::min(speed, max_delta);
          Quaternion dq = make_rotation(axis.normalized(), speed);
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
    // Cold authoring seam: a negative or inverted range makes the half-open
    // rand_int(min, max + 1) empty/inverted, yielding an implementation-defined
    // garbage delay that would fire at a nondeterministic time. Trap it at
    // construction like the file's other invariants rather than soft-fail later.
    HS_CHECK(min >= 0 && min <= max);
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
        // A repeating timer is constructed with duration=-1 and resets itself
        // here, so it never reaches done() and the Timeline never fires its
        // per-cycle .then(). Fire it directly so an attached callback honors
        // the documented per-cycle contract (the non-repeating branch reaches
        // done() via cancel(), so the Timeline fires its callback for it).
        this->post_callback();
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
        // See RandomTimer::step: a repeating timer never reaches done(), so
        // route its per-cycle .then() through post_callback() here to honor the
        // documented contract.
        this->post_callback();
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
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    // Snapshot the start value once, on the first step — not every time t hits
    // 0. A repeating Transition rewinds t to 0 at the end of each cycle, by
    // which point mutant == to; re-snapshotting would set from = to and every
    // subsequent cycle would write the constant `to` (a silent no-op). Capturing
    // here rather than in the ctor still picks up a mutant changed between
    // construction and the first step.
    if (!captured) {
      from = mutant;
      captured = true;
    }
    AnimationBase::step(canvas);
    // Clamp both ends (not just std::min on the top): a negative duration would
    // make t/duration negative and feed easing a negative arg. The lower bound
    // pins it to 0, matching Lerp/MeshMorph. The base ctor maps duration 0 -> 1,
    // so the divide is never by zero.
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
        easing_fn(std::move(easing_fn)), paused_(paused) {}

  /**
   * @brief Performs one step of the mutation.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details When wired to a pause flag (an effect's `anims_paused_`), the
   * mutation freezes while paused: it neither advances its timer nor writes the
   * mutant, so a GUI slider bound to the same member is the sole writer and the
   * user's value holds. Resuming hands the member back to the curve.
   */
  void step(Canvas &canvas) override {
    if (paused_ && *paused_)
      return;
    AnimationBase::step(canvas);
    // Clamp both ends like Transition: a negative duration must not feed easing a
    // negative arg. Base ctor maps duration 0 -> 1, so the divide is safe.
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
        paused_(paused) {}

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
  Driver(float &mutant, const float *speed_src, float scale, bool wrap = true,
         const bool *paused = nullptr)
      : AnimationBase(1, true), mutant(mutant), speed(0.0f), wrap_(wrap),
        paused_(paused), speed_src_(speed_src), scale_(scale) {
    // Guard the live source before dereferencing it.
    HS_CHECK(speed_src != nullptr, "Driver: live speed_src is null");
    speed = *speed_src * scale;
  }

  /**
   * @brief Performs one step by adding the speed to the mutant.
   * @param canvas The canvas buffer (forwarded to the base step).
   * @details Freezes while a wired pause flag is set (see Mutation::step).
   */
  void step(Canvas &canvas) override {
    if (paused_ && *paused_)
      return;
    AnimationBase::step(canvas);
    // A bound Driver re-reads its live speed source each step (one load +
    // multiply on the once-per-frame timeline path, not the pixel loop).
    if (speed_src_)
      speed = *speed_src_ * scale_;
    mutant.get() += speed;
    if (wrap_) {
      // wrap_t (fract) folds the phase back into [0,1) for any magnitude, so
      // |speed| >= 1 or a resume-jump far outside the range still wraps cleanly.
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
   */
  void set_speed(float new_speed) { speed = new_speed; }

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

  // Borrow contract: subject/start/target are stored as raw pointers and
  // dereferenced in step() across many frames, so they must be caller-owned
  // lvalues that outlive the Timeline. `T &subject` already rejects a temporary
  // subject; these deleted overloads reject a temporary start/target too, so an
  // inline `Lerp(x, Foo{...}, bar, ...)` is a compile error instead of a silent
  // dangling read (mirrors Motion/MeshMorph). The trailing defaulted `paused`
  // mirrors the real ctor so an rvalue start/target is rejected whether or not a
  // pause flag is passed — otherwise a 6-arg paused call would skip these 5-arg
  // overloads and bind the temporary to the real ctor.
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
    if (paused_ && *paused_)
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
        fade_out_easing(std::move(fade_out_easing_fn)), paused_(paused) {}

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
    // Paused: hold the frame — don't advance the timer (so the sprite neither
    // fades nor expires) but keep drawing at the current opacity. Lets a paused,
    // param-driven effect keep rendering its held state instead of going blank
    // when the sprite that renders it would otherwise run out (see HankinSolids,
    // whose render sprite is the same length as the angle Mutation it pauses).
    if (!(paused_ && *paused_))
      AnimationBase::step(canvas);

    // Trapezoid envelope as the MIN of an independent fade-in ramp and fade-out
    // ramp, each 1.0 outside its window. Computing both (rather than an
    // if/else-if where fade-in masks fade-out) keeps opacity continuous when the
    // windows overlap (fade_in_duration + fade_out_duration > duration),
    // degrading smoothly to a triangle instead of jumping at the handoff —
    // duration and fade are independent GUI sliders, so the overlap is
    // user-reachable. In the non-overlapping case only one ramp is ever below
    // 1.0 at a time, so the MIN reduces to a plain trapezoid.
    float fade_in = 1.0f;
    if (fade_in_duration > 0 && t < fade_in_duration) {
      float progress = static_cast<float>(t) / fade_in_duration;
      fade_in = fade_in_easing(hs::clamp(progress, 0.0f, 1.0f));
    }

    // An indefinite sprite (duration < 0) has no end frame to fade toward, so
    // fade_out_duration is intentionally ignored for it — the `duration >= 0`
    // guard is load-bearing, not just overflow safety. Fade a sprite out by
    // giving it a finite duration.
    float fade_out = 1.0f;
    if (duration >= 0 && fade_out_duration > 0 &&
        t >= (duration - fade_out_duration)) {
      float elapsed = static_cast<float>(t - (duration - fade_out_duration));
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

/**
 * @brief Sub-rotation count needed to keep each interpolated step within
 * `max_angle` over a total sweep of `angle` radians.
 * @param angle Total sweep angle in radians.
 * @param max_angle Per-column smoothness threshold in radians.
 * @return Sub-step count, always >= 1 (tight ceil, no extra subdivision).
 * @details The caller upsamples the orientation trail to (result + 1) frames,
 * so the loop's (len-1) sub-intervals each stay <= max_angle. Shared by Motion
 * and Rotation so both size the trail identically for the same angle.
 */
inline int rotation_substeps(float angle, float max_angle) {
  return static_cast<int>(std::ceil(std::max(1.0f, angle / max_angle)));
}

/**
 * @brief An animation that moves an Orientation along a defined Path.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class Motion : public AnimationBase<Motion<W, CAP>> {
  // Motion/Rotation interpolate across (len-1) sub-intervals of the bound
  // Orientation, where len <= CAP. CAP == 1 collapses to len == 1, leaving the
  // sweep with zero sub-intervals (Motion's loop never runs, Rotation divides by
  // zero). At least two sub-frames are required for any motion to be applied.
  static_assert(CAP >= 2, "Motion requires an Orientation capacity >= 2");

public:
  /**
   * @brief Constructs a Motion animation.
   * @tparam P Path type exposing get_point(float t).
   * @param orientation The Orientation object to update.
   * @param path_obj The path to follow (must have get_point(float t)).
   * @param duration The duration in frames.
   * @param repeat If true, the motion repeats.
   * @param space Coordinate space for the applied rotations.
   */
  template <typename P>
  Motion(Orientation<CAP> &orientation, const P &path_obj, int duration,
         bool repeat = false, Space space = Space::World)
      : AnimationBase<Motion<W, CAP>>(duration, repeat),
        orientation(orientation),
        path_fn([&path_obj](float t) { return path_obj.get_point(t); }),
        space(space) {}

  // Borrow contract: path_fn captures &path_obj and dereferences it every
  // frame for the animation's lifetime (which lives in the effect's Timeline).
  // The path must therefore be an effect-owned lvalue that outlives the
  // timeline — reject binding to a temporary at compile time instead of
  // dangling silently. (Non-owning by design: see Fn<>/inline-storage budget.)
  template <typename P,
            typename = std::enable_if_t<!std::is_lvalue_reference_v<P>>>
  Motion(Orientation<CAP> &orientation, P &&path_obj, int duration,
         bool repeat = false, Space space = Space::World) = delete;

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return orientation.get(); }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation.
   */
  const void *orientation_id() const override { return &get_orientation(); }

  /**
   * @brief Live-updates the traversal duration (frames per path loop).
   * @param frames New duration in frames; 0 is clamped to 1.
   * @details Guards against 0 (which would divide-by-zero in step()'s
   * `t / duration`), matching the base ctor's `duration == 0 ? 1` rule that
   * this direct write would bypass.
   */
  void set_duration(int frames) { this->duration = (frames == 0 ? 1 : frames); }

  /**
   * @brief Steps the animation, calculates intermediate rotation steps along
   * the path, and pushes them to the Orientation.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase<Motion<W, CAP>>::step(canvas);

    // Drift-free, composable integration. Motion advances the bound Orientation
    // by RELATIVE deltas only — it never writes the Orientation's absolute state
    // — so it composes with any other animation driving the same Orientation
    // (the multi-animation motion blur the collapse pre-pass in Timeline::step
    // supports). Each delta is the rotation carrying the path's own orientation
    // frame from one (sub)frame to the next, where that frame is a *pure
    // function* of the path parameter: the point path(s) plus the travel-tangent
    // pin a full orthonormal frame (see path_frame()). Drift-freedom comes from
    // that purity: the baseline frame is always a fresh path-derived value, never
    // an accumulating quaternion chain, so the per-frame deltas telescope and no
    // float error builds up — across unbounded cycles each cycle re-derives the
    // same frames. (A two-vector frame is also a global section with no holonomy,
    // so a closed path loops seamlessly.) That is why no per-cycle re-anchor (and
    // so no sole-ownership invariant) is needed: the old design dead-reckoned
    // parallel-transport deltas whose product drifted, and reset the whole
    // Orientation to an absolute anchor each cycle — which clobbered any
    // co-driver. The roll this frame pins is about the head's own axis, invisible
    // for a path-following point head, so the visible motion matches the prior
    // parallel-transport behavior.
    float t_prev = static_cast<float>(this->t - 1);
    Vector current_v = path_fn(t_prev / this->duration);
    float t_curr = static_cast<float>(this->t);
    Vector target_v = path_fn(t_curr / this->duration);
    float total_angle = angle_between(current_v, target_v);
    int num_steps = rotation_substeps(total_angle, MAX_ANGLE);

    // Ensure sufficient resolution. Past CAP this soft-clamps and the smear
    // sub-samples (graceful degradation, never a trap) — see Orientation::upsample.
    orientation.get().upsample(num_steps + 1);
    int len = orientation.get().length();

    // Baseline = the frame Motion last drove the Orientation to (carried across
    // frames, NOT recomputed from t-1). Within a cycle these coincide, but at a
    // repeat seam the phase rewinds to 0 while the Orientation still sits at the
    // end-of-cycle frame; using the carried baseline makes the first delta of
    // the new cycle carry the head from where it actually is (end of path) back
    // to the start, re-seating it for the next traversal. For a closed path the
    // start and end frames coincide and that delta is identity (a seamless
    // loop); for an open arc it is the jump-back the old absolute reset did —
    // but expressed as a *relative* delta, so a co-driver's contribution rides
    // along instead of being clobbered. At i==0 the delta is identity, so the
    // loop starts at i==1 and at(0) is left as the collapsed/co-driven seed.
    if (!have_prev_frame_) {
      prev_frame_ = path_frame(t_prev / this->duration);
      have_prev_frame_ = true;
    }
    Quaternion base_inv = prev_frame_.conjugate();
    Quaternion frame = prev_frame_;

    for (int i = 1; i < len; ++i) {
      // i maps the sub-interval [t-1, t]: i=0 is t-1, i=len-1 is t.
      float sub_t = t_prev + (static_cast<float>(i) / (len - 1));
      frame = path_frame(sub_t / this->duration);

      // The advance from the baseline frame to this substep's frame, applied as
      // a relative delta onto the existing orientation: World pre-multiplies
      // (frame * base_inv), Local post-multiplies (base_inv * frame).
      Quaternion &current_q = orientation.get().at(i);
      if (space == Space::Local) {
        current_q = current_q * (base_inv * frame);
      } else {
        current_q = (frame * base_inv) * current_q;
      }
      current_q.normalize();
    }
    // Carry the final substep's frame (= F at the current phase) as the next
    // frame's baseline. A pure path-derived value, never an accumulating chain,
    // so the per-frame deltas telescope and no float drift builds up.
    prev_frame_ = frame;
  }

  /**
   * @brief The path's own orientation frame at path parameter s.
   * @param s Normalized path parameter (the value passed to path_fn).
   * @return A unit quaternion mapping body +X to the point path(s), body +Y to
   * the travel direction tangent to the sphere, and body +Z to their cross.
   * @details A pure function of s, so consecutive-frame deltas telescope over a
   * cycle (no drift). The tangent is a forward finite difference with a fixed
   * step, so it depends only on s — not on the per-frame substep spacing. A
   * degenerate tangent (a momentarily stationary path point) falls back to a
   * deterministic perpendicular of the point, keeping path_frame a pure function
   * of s rather than reaching for prior state.
   */
  Quaternion path_frame(float s) const {
    Vector point = path_fn(s).normalized();
    Vector ahead = path_fn(s + FRAME_TANGENT_H).normalized();
    Vector travel = ahead - point;
    Vector tangent = travel - dot(travel, point) * point;
    // Fallback when the tangent vanishes: a cross with the body axis least
    // parallel to the point (mirrors make_rotation/RandomWalk's reference pick).
    Vector seed =
        std::abs(dot(point, X_AXIS)) < math::COS_AXIS_PARALLEL ? X_AXIS : Y_AXIS;
    Vector b1 = normalized_or(tangent, cross(seed, point).normalized());
    Vector b2 = cross(point, b1);
    return quaternion_from_basis(point, b1, b2);
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  /** Forward step (in path-parameter space) used to finite-difference the
   * travel tangent in path_frame(). Fixed so the frame is a pure function of s;
   * small enough to track curvature, large enough to clear float cancellation on
   * unit-vector differences. */
  static constexpr float FRAME_TANGENT_H = 1e-3f;
  std::reference_wrapper<Orientation<CAP>>
      orientation;               /**< Reference to the Orientation state. */
  Fn<Vector(float), 16> path_fn; /**< Function to retrieve path points. */
  Space space;                   /**< The coordinate space for rotation. */
  Quaternion prev_frame_;        /**< Path frame Motion last drove to; the
                                      baseline carried across the repeat seam. */
  bool have_prev_frame_ = false; /**< Whether prev_frame_ has been seeded. */
};

/**
 * @brief An animation that applies a fixed, time-eased rotation.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class Rotation : public AnimationBase<Rotation<W, CAP>> {
  // CAP == 1 collapses len to 1 and makes step_angle = delta / (len - 1) a
  // divide-by-zero; at least two sub-frames are required (see Motion).
  static_assert(CAP >= 2, "Rotation requires an Orientation capacity >= 2");

public:
  /**
   * @brief Default constructor. Creates an inactive/identity rotation.
   */
  Rotation()
      : AnimationBase<Rotation<W, CAP>>(0, false), orientation(nullptr),
        axis(X_AXIS), total_angle(0), easing_fn(ease_linear), last_angle(0),
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
  Rotation(Orientation<CAP> &orientation, const Vector &axis, float angle,
           int duration, EasingFn easing_fn, bool repeat = false,
           Space space = Space::World)
      : AnimationBase<Rotation<W, CAP>>(duration, repeat),
        orientation(&orientation), axis(axis), total_angle(angle),
        easing_fn(std::move(easing_fn)), last_angle(0), space(space) {}

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return *orientation; }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation (nullptr if unbound).
   */
  const void *orientation_id() const override { return orientation; }

  /**
   * @brief Checks if the rotation has a valid orientation bound.
   * @return True if an Orientation is bound.
   */
  bool has_orientation() const { return orientation != nullptr; }

  /**
   * @brief Steps the animation, calculates the incremental rotation delta, and
   * pushes it to the Orientation.
   * @param canvas The canvas buffer (forwarded to the base step).
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
      // Leave last_angle untouched so this sub-threshold increment is not
      // dropped: it accumulates into delta on later frames and is applied once
      // the sum crosses TOLERANCE. Advancing last_angle to target_angle here
      // would discard it forever, freezing very slow rotations and
      // systematically undershooting the slow-start portion of eased ones.
      return;
    }

    // Same tight substep count as Motion (shared rotation_substeps).
    int num_steps = rotation_substeps(std::abs(delta), MAX_ANGLE);
    // Past CAP this soft-clamps and the smear sub-samples (graceful
    // degradation, never a trap) — see Orientation::upsample.
    orientation->upsample(num_steps + 1);
    int len = orientation->length();

    float step_angle = delta / (len - 1);

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
  static void animate(Canvas &canvas, Orientation<CAP> &orientation,
                      const Vector &axis, float_t angle, EasingFn easing_fn,
                      Space space = Space::World) {
    Rotation<W, CAP> r(orientation, axis, angle, 1, easing_fn, false, space);
    r.step(canvas);
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  Orientation<CAP> *orientation; /**< Pointer to the Orientation state. */
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
 *
 * PERPETUAL — never completes: it runs with duration = -1 and does NOT repeat
 * or self-reset, so it never reaches done(). A `.then()` callback attached to a
 * RandomWalk therefore NEVER fires (unlike RandomTimer/PeriodicTimer, which are
 * also indefinite but self-trigger their per-cycle hook from step()). Do not
 * chain sequencing logic off a RandomWalk's `.then()`; drive follow-on behavior
 * from a finite animation or cancel() the walk explicitly. (The Transformer's
 * slot-recycling `.then()` likewise will not fire for a spawned RandomWalk — see
 * finding 419 on the repeating-spawn slot leak.)
 * @tparam W The width of the LED display.
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class RandomWalk : public AnimationBase<RandomWalk<W, CAP>> {
public:
  /**
   * @brief Tunable parameters for the random walk.
   */
  struct Options {
    float speed = 0.02f; /**< Movement speed per frame. */
    float pivot_strength =
        0.1f; /**< Strength of the direction change (noise amplitude). */
    float noise_scale = 0.02f;  /**< Frequency of the Perlin noise. */
    float smoothing = 0.85f;    /**< Angular momentum (0 = none, 0.95 = very sluggish). */
    float drift = 0.5f;         /**< Temporal drift speed for spatial noise. */

    /**
     * @brief Preset for slow, gentle motion.
     * @return Options tuned for a languid walk.
     */
    static Options Languid() { return {0.02f, 0.1f, 0.02f, 0.85f, 0.5f}; }

    /**
     * @brief Preset for fast, lively motion.
     * @return Options tuned for an energetic walk.
     */
    static Options Energetic() { return {0.05f, 0.4f, 0.08f, 0.7f, 1.0f}; }
  };

  /**
   * @brief Constructs a RandomWalk animation.
   * @param orientation The Orientation object to update.
   * @param v_start The starting direction vector.
   * @param noise External noise generator (caller owns lifetime).
   * @param options Configuration options.
   * @param seed Noise seed; 0 selects a random seed.
   */
  RandomWalk(Orientation<CAP> &orientation, const Vector &v_start,
             FastNoiseLite &noise, Options options = Options(),
             int seed = 0)
      : AnimationBase<RandomWalk<W, CAP>>(-1, false), orientation(orientation),
        v(Vector(v_start).normalized()), options(options),
        noiseGenerator(noise) {
    Vector u = X_AXIS;
    if (std::abs(dot(v, u)) > math::COS_AXIS_PARALLEL) {
      u = Y_AXIS;
    }
    direction = cross(v, u).normalized();
    noiseGenerator.get().SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noiseGenerator.get().SetFrequency(options.noise_scale);
    if (seed == 0) {
      noiseGenerator.get().SetSeed(static_cast<int>(hs::random()()));
    } else {
      noiseGenerator.get().SetSeed(seed);
    }
  }

  /**
   * @brief Live-updates the per-frame movement speed (e.g. from a GUI slider).
   * @param new_speed The new per-frame movement speed.
   */
  void set_speed(float new_speed) { options.speed = new_speed; }

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return orientation.get(); }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation.
   */
  const void *orientation_id() const override { return &get_orientation(); }

  /**
   * @brief Steps the walk: pivots direction based on noise, then rotates the
   * view along the calculated axis.
   * @param canvas The canvas buffer (forwarded to the base step and rotation).
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
    // If v and direction drift near-parallel the cross collapses to zero; fall
    // back to a deterministic perpendicular of v (the ctor's basis choice)
    // rather than normalizing a zero-length axis into NaN.
    const Vector axis_seed =
        std::abs(dot(v, X_AXIS)) > math::COS_AXIS_PARALLEL ? Y_AXIS : X_AXIS;
    Vector walk_axis =
        normalized_or(cross(v, direction), cross(v, axis_seed).normalized());
    v = rotate(v, make_rotation(walk_axis, options.speed)).normalized();
    direction =
        rotate(direction, make_rotation(walk_axis, options.speed)).normalized();
    Rotation<W, CAP>::animate(canvas, orientation, walk_axis, options.speed,
                              ease_linear);
  }

private:
  std::reference_wrapper<Orientation<CAP>>
      orientation;  /**< Reference to the global Orientation state. */
  Vector v;         /**< Current forward direction vector. */
  Vector direction; /**< Current pivoting direction (orthogonal to v). */
  Options options;  /**< Configuration options. */
  float angular_velocity = 0.0f; /**< Smoothed pivot rate (angular momentum). */
  std::reference_wrapper<FastNoiseLite>
      noiseGenerator; /**< External noise generator. */
};

// ---------------------------------------------------------------------------

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
   * @param canvas The canvas buffer (forwarded to the base step).
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

  // Borrow contract: num_rings/num_lines are stored as reference_wrappers and
  // re-read in step() across many frames (live-tracking), so they must be
  // caller-owned lvalues that outlive the Timeline. `MobiusParams &params`
  // already rejects a temporary; these deleted overloads reject a temporary
  // scalar too, so `MobiusFlow(p, 3.0f, 4.0f, ...)` is a compile error instead
  // of a silent dangling read (mirrors Lerp/Motion/MeshMorph).
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
    // Clamp to [0,1] like every sibling animation (Transition/Mutation/
    // MobiusWarp): repeat=true rewinds t in practice, but a non-repeating or
    // overshooting duration would otherwise feed out-of-range progress to the
    // warp.
    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    // num_rings is live-bound; if it ever reaches -1 the divisor num_rings + 1
    // hits 0 and logPeriod blows up to inf, poisoning a/d with inf/NaN. Floor
    // at 0 so num_rings + 1 stays >= 1, mirroring the num_lines guard below.
    float rings = num_rings;
    if (rings < 0.0f)
      rings = 0.0f;
    float logPeriod = 5.0f / (rings + 1);
    float flowParam = progress * logPeriod;
    float scale = expf(flowParam);
    float s = sqrtf(scale);
    // num_lines is live-bound and its GUI slider bottoms out at 0, so 2π/0
    // would blow angle up to inf. Clamp to >= 1 before dividing, mirroring the
    // defensive floor in set_period/set_duration.
    float lines = num_lines;
    if (lines < 1.0f)
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
        easing(easing) {}

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
  void bind_scale(const float &live_scale) { scale_ref_ = &live_scale; }

  /**
   * @brief Steps the animation, updating param b.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    float s = scale_ref_ ? *scale_ref_ : scale;
    params.get().b.re = s * (cosf(angle) - 1.0f);
    params.get().b.im = s * sinf(angle);
  }

  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  float scale;     /**< Warp magnitude (public: effects read/write directly). */
  EasingFn easing; /**< Easing curve (public: effects read/write directly). */

private:
  const float *scale_ref_ = nullptr; /**< Optional live magnitude source. */
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

  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  float scale;     /**< Warp magnitude (public: effects read/write directly). */
  EasingFn easing; /**< Easing curve (public: effects read/write directly). */
};



/**
 * @brief Animates a vertex-interpolated crossfade between two meshes.
 * @details Owns its transient state (cloned meshes + SLERP buffers) via a
 * pointer to arena-allocated storage — keeps inline size small for
 * TimelineEvent. The caller provides source/dest MeshState references and an
 * Arena; MeshMorph clones both, builds nearest-vertex correspondence, and
 * interpolates each frame. Transients are released when the animation completes
 * and the TimelineEvent is destroyed; call compact on the arena afterward to
 * reclaim the space.
 *
 * Crossfade contract (intentional): only the incoming mesh (mesh_B, carrying
 * dest topology) morphs — each frame its vertices SLERP from their nearest
 * source vertex toward their dest position. The outgoing mesh (mesh_A, a clone
 * of source) holds its geometry and fades out via opacity (op_A = 1 - alpha)
 * while the incoming fades in; the opacities sum to 1 for constant total
 * brightness across the topology swap. The source is cloned (not borrowed) so
 * the animation is self-contained against the caller recycling/mutating source
 * mid-morph, consistent with the borrow contract enforced on the draw callbacks
 * below. The separate draw_outgoing/draw_incoming callbacks exist precisely to
 * shade the two halves independently (see HankinSolids).
 */
class MeshMorph : public AnimationBase<MeshMorph> {
public:
  /**
   * @brief Non-owning per-half draw callback: `void(Canvas&, const MeshState&,
   * float opacity)`. A StoredFunctionRef (not a plain FunctionRef) because it is
   * held as a member and invoked across many frames: the type rejects rvalue
   * temporaries so a dangling inline lambda is a compile error, not a silent
   * use-after-free (see the borrow contract below).
   */
  using MorphDrawFn = StoredFunctionRef<void(Canvas &, const MeshState &, float)>;

  /**
   * @brief Constructs a MeshMorph with separate shading for the two halves.
   * @param source The outgoing mesh (cloned, not borrowed).
   * @param dest The incoming mesh whose topology the morph targets.
   * @param arena Arena providing backing storage for cloned meshes and buffers.
   * @param draw_outgoing Draw callback for the fading-out source clone.
   * @param draw_incoming Draw callback for the fading-in morphing mesh.
   * @param duration The crossfade duration in frames.
   * @param easing_fn The easing function applied to crossfade progress.
   */
  MeshMorph(const MeshState &source, const MeshState &dest, Arena &arena,
            MorphDrawFn draw_outgoing, MorphDrawFn draw_incoming, int duration,
            EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(duration, false), easing_fn(easing_fn),
        draw_outgoing(draw_outgoing), draw_incoming(draw_incoming) {
    // An empty source leaves best_idx pinned at 0 and indexes an empty
    // ArenaVector when building the correspondence below — an OOB read under
    // NDEBUG. The nearest-vertex map is meaningless without source vertices;
    // trap at construction (house style).
    HS_CHECK(!source.vertices.is_empty());
    // Allocate transient storage on the arena
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();

    // Clone both meshes for interpolation
    MeshOps::clone(source, buf_->mesh_A, arena);
    MeshOps::clone(dest, buf_->mesh_B, arena);

    // step() writes interpolated positions into mesh_B.vertices[i] for i in
    // [0, dest.vertices.size()), which is in-bounds only if clone() preserved
    // dest's vertex count one-to-one. If clone ever welds/dedups, trap at
    // construction (cold) rather than writing OOB on the per-frame step.
    HS_CHECK(buf_->mesh_B.vertices.size() == dest.vertices.size());

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

  // Borrow contract: draw_outgoing/draw_incoming are stored as non-owning
  // StoredFunctionRefs (8 bytes each — deliberately not owning, to keep
  // MeshMorph inside the Timeline inline-storage budget) and invoked in step()
  // across many frames. The callables must be effect-owned lvalues that outlive
  // the timeline; StoredFunctionRef rejects a temporary (e.g. an inline lambda)
  // at the MorphDrawFn parameter, turning a dangling bind into a compile error
  // rather than a silent use-after-free. No per-ctor `= delete` overload is
  // needed — the type carries the rule.

  /**
   * @brief Steps the crossfade: interpolates vertices and renders both halves.
   * @param canvas The canvas buffer passed to the draw callbacks.
   */
  void step(Canvas &canvas) override {
    // Increment-first (same ordering as Transition/Mutation) is deliberate: the
    // rendered progress runs 1/duration .. duration/duration, so the final frame
    // lands exactly on the destination mesh (alpha == 1) — the correct terminal
    // state for a morph. The skipped pure-source frame (progress == 0) is
    // immaterial: at progress == 1/duration the easing is still ~0, so the first
    // rendered frame already draws the source at op_A ~= 1.
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
  /**
   * @brief Arena-allocated transient data — keeps MeshMorph inline size small.
   */
  struct Transients {
    MeshState mesh_A;            /**< Outgoing mesh clone. */
    MeshState mesh_B;           /**< Incoming morphing mesh clone. */
    ArenaVector<Vector> start_pos; /**< Per-vertex nearest-source start points. */
    ArenaVector<Vector> end_pos;   /**< Per-vertex dest end points. */
  };

  Transients *buf_;          /**< Pointer to arena-allocated transient state. */
  EasingFn easing_fn;        /**< Easing curve applied to crossfade progress. */
  MorphDrawFn draw_outgoing; /**< Draw callback for the outgoing half. */
  MorphDrawFn draw_incoming; /**< Draw callback for the incoming half. */
};

/**
 * @brief Continuously modulates Mobius parameters to create an evolving warp.
 * @details Uses multiple frequencies for non-repeating chaos.
 */
class MobiusWarpEvolving : public AnimationBase<MobiusWarpEvolving> {
public:
  /**
   * @brief Constructs a MobiusWarpEvolving animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   */
  MobiusWarpEvolving(MobiusParams &params, float scale = 0.5f,
                     float speed = 0.01f)
      : speed(speed), scale(scale), params(params), base(params),
        seed(hs::random()()) {}

  /**
   * @brief Derives a per-channel phase offset from the seed and channel index.
   * @param i Channel index.
   * @return Phase offset in [0, 100) for that channel.
   */
  float phase(int i) const {
    uint32_t h = seed ^ (static_cast<uint32_t>(i) * 2654435761u);
    return (h & 0xFFFF) * (100.0f / 65535.0f);
  }

  /**
   * @brief Steps the animation, modulating all eight Mobius params.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
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

  float speed; /**< Animation speed (radians of phase per frame unit). */
  float scale; /**< Magnitude of the per-channel modulation. */

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
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
  float frequency{20.0}; /**< Spatial frequency of the wave. */
  float decay{5.0};      /**< Spatial decay rate. */
  float thickness{1.0f}; /**< Thickness of the ripple. */

  /** @brief Cached cos(angle) lower fast-reject bound. */
  float cos_threshold_min = 1.0f;
  /** @brief Cached cos(angle) upper fast-reject bound. */
  float cos_threshold_max = -1.0f;

  /**
   * @brief Recomputes the cos(angle) fast-reject bounds for the wavelet's
   * current phase and thickness, so the renderer can skip points outside the
   * ring.
   */
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

  /**
   * @brief Steps the ripple: advances the wave and updates its envelope.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
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

    // Re-prepare the fast-reject thresholds against the phase we just advanced.
    // The transformer's prepare_frame() runs before timeline.step(), so without
    // this the render samples the wavelet at the NEW phase but rejects against
    // thresholds cached at the OLD phase — clipping the leading `speed`-wide
    // slice of every ripple (the whole leading half at extreme Ripp Width / Dur,
    // where speed approaches the 2·half-width window). Recomputing here keeps
    // phase and thresholds consistent in any step ordering. Cold (≤ pool-size
    // ripples per frame, 2 cosf each).
    params.get().prepare_thresholds();
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
};

/**
 * @brief Animates noise parameters by updating time.
 */
class Noise : public AnimationBase<Noise> {
public:
  /**
   * @brief Constructs a Noise animation.
   * @param params Reference to the NoiseParams to animate.
   * @param duration Duration in frames (-1 for indefinite).
   */
  Noise(NoiseParams &params, int duration = -1)
      : AnimationBase(duration, true), params(params) {}

  /**
   * @brief Steps the animation, advancing the noise time field.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    // Accepted divergence: t is an integer frame counter, so once it exceeds
    // 2^24 (~16.7M frames, ~77 h at 60 fps) float can no longer represent
    // consecutive values and the noise time axis quantizes — the field appears
    // to slow then freeze on a permanent install. A modulo-before-cast would
    // hold float precision but OpenSimplex2 is aperiodic, so every wrap injects
    // a visible discontinuity; a float accumulator quantizes the same way once
    // it passes 2^24. Both fixes trade a gradual, barely-perceptible drift for a
    // sharp artifact, so the raw cast is kept and the limit documented.
    params.get().time = static_cast<float>(t);
  }

private:
  std::reference_wrapper<NoiseParams> params; /**< Noise params to animate. */
};

} // namespace Animation

/**
 * @brief Structure linking an animation with its starting time.
 * @details Stores the animation inline to avoid arena allocation (survives
 * compaction).
 */
struct TimelineEvent {
  // Inline storage budget for a type-erased animation. Tuned to 112 B for the
  // device (teensy::inplace_function, 24 B/callable) and the WASM simulator
  // (libc++ std::function, 32 B/callable). The native unit-test build uses
  // MSVC's std::function (64 B/callable), which inflates the same animation
  // types past 112 B; HS_TEST_BUILD widens the budget for that build only so
  // the device's RAM footprint is unchanged. See the `tests` CMake preset
  // (or `just test`).
#ifdef HS_TEST_BUILD
  static constexpr size_t MAX_ANIM_SIZE = 256;
#else
  static constexpr size_t MAX_ANIM_SIZE = 112;
#endif

  int start = 0; /**< Global frame at which the animation begins stepping. */
  /**
   * @brief Whether this event's inline animation pointer was handed out via
   * Timeline::add_get().
   * @details Compaction must never relocate such an event — doing so dangles the
   * caller's cached pointer. (Occupies existing alignment padding before
   * `storage`, so it costs no extra bytes.)
   */
  bool handled = false;
  alignas(std::max_align_t) uint8_t storage[MAX_ANIM_SIZE]; /**< Inline
                                  type-erased animation storage. */

  /**
   * @brief Type-erased move/destroy operation for the stored animation.
   * @details dst != nullptr → move src into dst and destroy src.
   *          dst == nullptr → just destroy src.
   */
  void (*manager)(TimelineEvent &src, TimelineEvent *dst) = nullptr;

  /**
   * @brief The animation's IAnimation* view, captured by static_cast at
   * construction (and re-captured by the manager on every move).
   * @details reinterpret_cast<IAnimation*> on the raw storage would be formally
   * UB: animation types are non-standard-layout (virtual functions), so the
   * standard does not guarantee the IAnimation base subobject sits at offset 0
   * — only a properly-typed upcast performs the (here-zero, but not
   * portably-zero) base adjustment. Trails `manager` so it lands in the existing
   * tail padding (no extra bytes).
   */
  IAnimation *iface = nullptr;

  /**
   * @brief Gets the live animation interface for this slot.
   * @return The IAnimation* view, or nullptr if the slot is empty.
   */
  IAnimation *animation() { return manager ? iface : nullptr; }

  /**
   * @brief Relocates this event's animation into a destination slot.
   * @param dst The destination slot to move into.
   */
  void move_into(TimelineEvent &dst) {
    // Fail-fast on the add_get() dangling-handle hazard: relocating a handled
    // event invalidates the caller's cached animation pointer. Today's callers
    // are safe by an unstated invariant (handled animations are infinite and
    // added before any finite one, so compaction never shifts them), but a
    // future change that violates it would silently corrupt memory. Trap at the
    // move instead — a bench-time crash beats a live use-after-free. (Cold path:
    // once per relocated event per frame, never per pixel.)
    HS_CHECK(!handled);
    dst.start = start;
    dst.handled = handled; // always false past the check; kept for symmetry
    dst.manager = manager;
    if (manager) {
      manager(*this, &dst);
      manager = nullptr;
    }
  }

  /**
   * @brief Destroys the stored animation and empties the slot.
   */
  void destroy() {
    if (manager) {
      manager(*this, nullptr);
      manager = nullptr;
    }
  }
};

// Device inline-storage budget audit (finding 373). The per-type
// static_assert(sizeof(A) <= MAX_ANIM_SIZE) in Timeline::add/add_get only fires
// for types that are actually add()ed in the build being compiled — so if an
// effect stops add()ing a heavy animation, that type can grow past the device
// budget unnoticed. Pin the budget for the largest on-device animation type here
// at namespace scope so it is checked in every build that includes this header,
// independent of any callsite.
//
// MAX_ANIM_SIZE — not a hardcoded 112 — is deliberate. sizeof of an Fn-/vtable-
// bearing animation is wider on the 64-bit native host (8-byte pointers) than on
// the 32-bit device/WASM (4-byte), so a literal `<= 112` would FAIL on the host
// even for a type that fits on device; that pointer-width gap is exactly why the
// native test build widens MAX_ANIM_SIZE to 256 (see TimelineEvent). The
// accurate device check is therefore the 32-bit WASM/device build, where
// MAX_ANIM_SIZE == 112: this assert enforces the real device budget there and is
// a no-op headroom check on the host. The WASM CI build, which compiles every
// registered effect, is the comprehensive gate for all add()ed animation types;
// this line additionally guards MobiusWarpEvolving — the heaviest type (it
// embeds a full MobiusParams snapshot) and the one the finding named — even if
// it ever becomes unreferenced. (Templated animations like Motion/Rotation<W,CAP>
// cannot be sized at namespace scope and stay covered at their add() sites.)
static_assert(sizeof(Animation::MobiusWarpEvolving) <=
                  TimelineEvent::MAX_ANIM_SIZE,
              "MobiusWarpEvolving exceeds the TimelineEvent inline-storage "
              "budget (on the 32-bit WASM/device build MAX_ANIM_SIZE is the "
              "112-byte device budget); shrink the type or raise the budget.");

/**
 * @brief Global storage for the timeline to prevent template instantiation
 * bloat.
 */
static constexpr int TIMELINE_MAX_EVENTS = 64;
extern DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
// True while a Timeline instance is alive. Guards the single-singleton invariant
// (see global_timeline_events): every Timeline shares that one event array, so a
// second live instance would silently stomp the first's events.
extern bool global_timeline_live;
// Bookkeeping cursors into global_timeline_events. Kept as free globals
// alongside the array (not Timeline statics): the storage is one shared
// singleton, so its frame counter and active-event count belong with it.
extern int global_timeline_t;          // current global frame count
extern int global_timeline_num_events; // current number of active events

/**
 * @brief Manages all active animations and their execution over time.
 *
 * Not templated: the entire instance is backed by one global event array
 * (`global_timeline_events`) sized to hold the largest effect, and the
 * live-guard permits exactly one instance at a time, so the capacity is a
 * single process-wide budget (`MAX_EVENTS`), not a per-instance knob.
 */
class Timeline {
public:
  /**
   * @brief Constructs a Timeline.
   *
   * Traps if one is already alive: all Timelines share global_timeline_events,
   * so two live instances would corrupt each other's events. The real app keeps
   * exactly one (destroys the old effect before building the next); this makes
   * the latent footgun a bench-time crash instead of silent stomping.
   */
  Timeline() {
    HS_CHECK(!global_timeline_live);
    global_timeline_live = true;
    clear();
  }

  /**
   * @brief Cleans up remaining animations, invoked on effect destruction.
   */
  ~Timeline() {
    clear();
    global_timeline_live = false;
  }

  // Singleton over global state: copying/moving would either bypass the live
  // guard or leave a moved-from husk whose dtor clears the survivor's events.
  // To reset in place, call clear() rather than reassigning a fresh instance.
  Timeline(const Timeline &) = delete;
  Timeline &operator=(const Timeline &) = delete;
  Timeline(Timeline &&) = delete;
  Timeline &operator=(Timeline &&) = delete;

  /**
   * @brief Destroys all events and resets the global frame counter.
   */
  void clear() {
    for (int i = 0; i < global_timeline_num_events; ++i) {
      global_timeline_events[i].destroy();
    }
    global_timeline_num_events = 0;
    global_timeline_t = 0;
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
    static_assert(alignof(A) <= alignof(std::max_align_t),
                  "Animation type is over-aligned for TimelineEvent inline "
                  "storage (placement-new would be misaligned)");
    // Deliberate soft-drop, not the usual fail-fast trap: pool exhaustion here
    // is a recoverable scheduling condition (an effect over-queues for one
    // frame), not a corrupted invariant, and it is pinned by tests. The chained
    // add() form returns *this regardless; callers that need the handle use
    // add_get(), which returns nullptr here and is null-checked at every site.
    if (global_timeline_num_events >= MAX_EVENTS) {
      hs::log("Timeline full, failed to add animation!");
      return *this;
    }
    auto &e = global_timeline_events[global_timeline_num_events++];
    e.start = global_timeline_t + (int)in_frames;
    e.handled = false; // global slots are reused — clear any stale handled flag
    e.iface = static_cast<IAnimation *>(new (e.storage) A(std::move(animation)));
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      // std::launder: the A was placement-new'd into the raw byte storage, so a
      // bare reinterpret_cast of the array address is not a pointer to the live
      // object. Launder it to recover a usable A* (matches the iface comment's
      // reasoning about typed pointer recovery from this same storage).
      A *obj = std::launder(reinterpret_cast<A *>(src.storage));
      if (dst) {
        dst->iface = static_cast<IAnimation *>(new (dst->storage) A(std::move(*obj)));
      }
      obj->~A();
    };
    return *this;
  }

  /**
   * @brief Like add(), but returns the typed pointer to the inline-stored
   * animation. Use when you need to hold a reference for later mutation.
   * @tparam A The animation type.
   * @param in_frames The number of frames to delay before starting.
   * @param animation The animation object.
   * @param pin If true (default), the caller intends to RETAIN this pointer
   * across frames, so the event is marked handled and step()'s compaction traps
   * (move_into) rather than relocating it out from under the cached pointer.
   * Such a retained handle is only safe when the animation is infinite and added
   * before any finite one (so no earlier event is ever removed to shift it) —
   * the contract the direct callers rely on; the trap enforces it. Pass
   * `pin=false` for a TRANSIENT pointer used only at the call site and not kept
   * across frames (e.g. TransformerPool::spawn, whose finite animations are
   * compacted normally and whose return is typically discarded).
   * @return Typed pointer to the inline-stored animation, or nullptr if full.
   */
  template <typename A> A *add_get(float in_frames, A animation, bool pin = true) {
    static_assert(sizeof(A) <= TimelineEvent::MAX_ANIM_SIZE,
                  "Animation type exceeds TimelineEvent inline storage");
    static_assert(alignof(A) <= alignof(std::max_align_t),
                  "Animation type is over-aligned for TimelineEvent inline "
                  "storage (placement-new would be misaligned)");
    if (global_timeline_num_events >= MAX_EVENTS) {
      hs::log("Timeline full, failed to add animation!");
      return nullptr;
    }
    auto &e = global_timeline_events[global_timeline_num_events++];
    e.start = global_timeline_t + (int)in_frames;
    // A pinned (retained) handle must stay put: step()'s compaction traps in
    // move_into if a later change ever tries to relocate this event.
    e.handled = pin;
    auto *ptr = new (e.storage) A(std::move(animation));
    e.iface = static_cast<IAnimation *>(ptr);
    e.manager = [](TimelineEvent &src, TimelineEvent *dst) {
      // std::launder: the A was placement-new'd into the raw byte storage, so a
      // bare reinterpret_cast of the array address is not a pointer to the live
      // object. Launder it to recover a usable A* (matches the iface comment's
      // reasoning about typed pointer recovery from this same storage).
      A *obj = std::launder(reinterpret_cast<A *>(src.storage));
      if (dst) {
        dst->iface = static_cast<IAnimation *>(new (dst->storage) A(std::move(*obj)));
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
    ++global_timeline_t;

    int write_idx = 0;
    int active_cnt = global_timeline_num_events; // Snapshot count before callbacks
                                 // potentially add more

    // Collapse each distinct Orientation exactly once, before any animation
    // steps it. Collapsing per-animation would discard the sub-frame
    // motion-blur history an earlier animation sharing the same Orientation had
    // already built this frame, breaking the multi-animation motion blur the
    // engine advertises. An Orientation is collapsed only by the first started
    // animation that references it.
    //
    // Track the orientations already collapsed this frame in a fixed-size
    // scratch set (at most MAX_EVENTS distinct ids, so it lives on the stack
    // with no allocation). The earlier form rescanned every prior event for
    // each orientation-owning one — O(active^2), and it paid that cost even
    // though most events own no orientation. Membership is now checked only
    // against the handful of distinct collapsed ids, so the common case (a few
    // orientations shared across many events) is effectively O(active).
    const void *collapsed_ids[MAX_EVENTS];
    int collapsed_n = 0;
    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];
      if (global_timeline_t < e.start)
        continue;
      IAnimation *anim = e.animation();
      if (!anim)
        continue;
      const void *id = anim->orientation_id();
      if (!id)
        continue;
      bool already_collapsed = false;
      for (int j = 0; j < collapsed_n; ++j) {
        if (collapsed_ids[j] == id) {
          already_collapsed = true;
          break;
        }
      }
      if (!already_collapsed) {
        anim->collapse_orientation();
        collapsed_ids[collapsed_n++] = id;
      }
    }

    for (int i = 0; i < active_cnt; ++i) {
      auto &e = global_timeline_events[i];

      // 1. Check start time
      if (global_timeline_t < e.start) {
        if (i != write_idx) {
          e.move_into(global_timeline_events[write_idx]);
        }
        write_idx++;
        continue;
      }

      // 2. Step (Orientation already collapsed once-per-frame above)
      IAnimation *anim = e.animation();
      // A started event always carries a live animation: add() sets
      // manager/iface, and compaction only nulls a slot the loop has already
      // passed (write_idx <= i). A null here is an upstream invariant
      // violation — trap rather than silently counting a moved-out husk as live.
      HS_CHECK(anim);
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
        // Symmetric with move_into's guard, closing the completion side of the
        // pin invariant. The contract is pinned ⇒ infinite, so a pinned event
        // should never reach natural completion; if one does — even as the last
        // event, where no relocation runs and move_into's trap would not fire —
        // destroying it dangles the caller's retained pointer silently. Trap so
        // misuse crashes at the bench instead of live. A deliberate cancel() is
        // exempt: that is the held handle's own sanctioned teardown (like
        // clear()), not a silent dangle. (Cold path: once per completed event,
        // never per pixel.)
        HS_CHECK(!e.handled || anim->is_canceled());
        e.destroy();
      }
    }

    // 5. Move new events (added during callbacks) to fill the gap left by
    //    completed ones. This relocates the new events, so a *pinned* event
    //    spawned inside a callback (add_get/spawn_pinned with pin=true) would
    //    trap here in move_into's HS_CHECK(!handled) — by design: its address
    //    was handed back to the caller and must not move. Callback-spawners use
    //    pin=false today (TransformerPool::spawn), so the gap-fill is safe; the
    //    trap stands as the fail-fast guard if that ever changes.
    int new_vals_count = global_timeline_num_events - active_cnt;
    if (new_vals_count > 0 && write_idx < active_cnt) {
      for (int i = 0; i < new_vals_count; ++i) {
        global_timeline_events[active_cnt + i].move_into(
            global_timeline_events[write_idx + i]);
      }
    }

    global_timeline_num_events = write_idx + new_vals_count;
  }

  // Frame counter (`global_timeline_t`) and active-event count
  // (`global_timeline_num_events`) live as free globals, not statics here, so
  // they stay one shared singleton with global_timeline_events.
  static constexpr int MAX_EVENTS =
      TIMELINE_MAX_EVENTS; /**< Must match global_timeline_events array size. */
};

/**
 * @brief A double-buffered pair of persistent MeshState slots plus the
 *        arena-compaction primitives effects need to swap between them.
 * @details Holds two MeshState slots in `persistent_arena` and a front/back
 * index. Effects own the transition policy themselves — they generate into a
 * slot, schedule whatever animation they want (an opacity-crossfade `Sprite`, a
 * geometric `MeshMorph`, ...), flip the front index, and reclaim the old slot's
 * space. This class only provides the slot storage, index management, and the
 * two compaction helpers (`compact`, `compact_keep_front`) that those policies
 * share; it intentionally does not encode any single transition shape, because
 * the real clients (IslamicStars, HankinSolids, MeshFeedback) diverge on
 * animation type and on the per-slot side-band state they carry.
 *
 * Usage:
 *   MeshCarousel carousel;  // in effect members
 *
 *   // Build the initial shape directly into the front slot:
 *   carousel.current().clear();
 *   MeshOps::compile(mesh, carousel.current(), persistent_arena);
 *
 *   // To transition: generate into the back slot, schedule an animation,
 *   // then flip and compact (see IslamicStars::spawn_shape for the pattern).
 */
class MeshCarousel {
public:
  /**
   * @brief Constructs an empty carousel with front slot index 0.
   */
  MeshCarousel() {}

  /**
   * @brief Gets the currently visible (front) mesh.
   * @return Const reference to the front MeshState slot.
   */
  const MeshState &current() const { return slots_[front_]; }

  /**
   * @brief Gets the currently visible (front) mesh (mutable).
   * @return Mutable reference to the front MeshState slot.
   */
  MeshState &current() { return slots_[front_]; }

  /**
   * @brief Direct slot access by index (for effects that need both).
   * @param i Slot index (0 or 1).
   * @return Const reference to the requested MeshState slot.
   */
  const MeshState &slot(int i) const { return slots_[i]; }

  /**
   * @brief Direct slot access by index (mutable).
   * @param i Slot index (0 or 1).
   * @return Mutable reference to the requested MeshState slot.
   */
  MeshState &slot(int i) { return slots_[i]; }

  /**
   * @brief Gets the front slot index (for capture in lambdas).
   * @return The index (0 or 1) of the front slot.
   */
  int front_index() const { return front_; }

  /**
   * @brief Manually sets the front index (for effects that manage transitions
   * themselves).
   * @param idx The new front slot index (0 or 1).
   */
  void set_front(int idx) { front_ = idx; }

  /**
   * @brief Compacts the persistent arena, reclaiming fragmented space.
   * @details Evacuates tracked MeshStates and resets the arena. Call before
   * allocating new persistent data.
   */
  void compact() {
    // NOTE: Both evacuations share scratch_arena_a. If both slots are
    // populated, scratch_arena_a must have room for both. An OOM here
    // will trigger the Arena::allocate assert.
    Persist<MeshState> p0(slots_[0], scratch_arena_a, persistent_arena);
    Persist<MeshState> p1(slots_[1], scratch_arena_a, persistent_arena);
    persistent_arena.reset();
  }

  /**
   * @brief Frees the back slot and compacts, preserving only the front slot.
   * @tparam AfterReset Callable type invoked as `void(Arena&)`.
   * @param after_reset Callback run immediately after the reset, while the front
   * slot is still evacuated.
   * @details Runs `after_reset(persistent_arena)` immediately after the reset —
   * while the front slot is still evacuated — so the caller can re-bake
   * effect-owned persistent data (e.g. a palette bank) into the fresh arena
   * *before* the front mesh is restored on top of it. Use when only the visible
   * (front) shape must survive a regeneration of the back slot.
   */
  template <typename AfterReset> void compact_keep_front(AfterReset after_reset) {
    int back = 1 - front_;
    slots_[back] = MeshState();
    Persist<MeshState> p(slots_[front_], scratch_arena_b, persistent_arena);
    persistent_arena.reset();
    after_reset(persistent_arena);
  }

private:
  MeshState slots_[2]; /**< Front/back double-buffered mesh slots. */
  int front_ = 0;      /**< Index (0 or 1) of the visible front slot. */
};

/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @tparam CAP Orientation sub-frame capacity.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&,
 * float t)`.
 */
template <int CAP>
void tween(const Orientation<CAP> &o, TweenFn callback) {
  int len = o.length();
  int start = (len > 1) ? 1 : 0;
  for (int i = start; i < len; ++i) {
    // A lone snapshot is the newest sub-position, so t = 1 (age-neutral); the
    // multi-frame sweep already ends at t = 1 for i = len - 1. Emitting t = 0
    // here would push age + 1 on every static orientation (see World::Orient).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(o.get(i), t);
  }
}

/**
 * @brief Helper to iterate over a VectorTrail's historical frames.
 * @tparam CAPACITY Trail capacity.
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
    // A lone sample is the newest (head) position, so t = 1 (age-neutral) — the
    // multi-sample sweep already ends at 1 for i = len - 1. Emitting 0 here
    // would flash a freshly spawned or dying particle at tail brightness for the
    // one frame its history holds a single point (mirrors tween(Orientation)).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(trail.get(i), t);
  }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or
 * Animation::OrientationTrail).
 * @param trail The Tweenable object to iterate.
 * @param callback The function to call for each step: `void(const T&, float
 * t)`.
 */
void deep_tween(const Tweenable auto &trail, TweenFn callback) {
  size_t trail_len = trail.length();
  if (trail_len == 0)
    return;

  // The age ramp must end at 1.0 on the newest *plotted* orientation. A frame
  // recorded with no intra-frame motion collapses to a single sub-frame, whose
  // lone quaternion is the shared boundary the previous frame already plotted —
  // it emits nothing after the sub-frame-0 skip, and plotting it would just
  // double up that junction. So a motionless tail would strand the ramp below
  // 1.0 if we normalized by trail_len. Find the newest frame that actually
  // contributes and normalize against it instead; the backward scan stops at
  // the first contentful frame (O(1) in the common case — not a second pass
  // over the emitted samples), and for an all-moving trail span == trail_len,
  // leaving the per-frame mapping unchanged.
  size_t last = trail_len - 1;
  while (last > 0 && trail.get(last).length() <= 1)
    --last;

  // A fully motionless trail collapses to frame 0's lone sub-position — the
  // newest (and only) plotted orientation. Like the single-snapshot tween()
  // overload, it must read t = 1.0 (age-neutral), not the 0.0 the sub-frame-0
  // fallback would yield; otherwise quintic_kernel(0) = 0 renders a static head
  // invisible. (A length-0 frame 0 emits nothing, matching the main loop.)
  if (last == 0 && trail.get(0).length() == 1) {
    const auto &frame = trail.get(0);
    callback(frame.get(0), 1.0f);
    return;
  }

  float span = static_cast<float>(last + 1);

  for (size_t i = 0; i <= last; ++i) {
    const auto &frame = trail.get(i);
    size_t frame_size = frame.length();
    size_t start_j = (i == 0) ? 0 : 1;

    for (size_t j = start_j; j < frame_size; ++j) {
      const auto &q = frame.get(j);
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      float global_t = (static_cast<float>(i) + sub_t) / span;
      callback(q, global_t);
    }
  }
}
