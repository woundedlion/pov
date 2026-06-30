/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "animation_core.h"

namespace Animation {

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
  static constexpr int kCapacity = CAPACITY; /**< Max retained snapshots. */

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
   * @param l Life in frames; clamped into [0, 65535]. The low clamp matters: a
   *   negative l would otherwise narrow to a huge uint16_t (a near-immortal
   *   particle) rather than a dead one.
   */
  void init(const Vector &p, const Vector &v, uint16_t seed, float l) {
    position = p;
    velocity = v;
    color_seed = seed;
    life = static_cast<uint16_t>(hs::clamp(l, 0.0f, 65535.0f));
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
  static_assert(CAPACITY <= 65535,
                "active_count is uint16_t; CAPACITY must fit in it");

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
    // The arena has no per-allocation free, so re-binding would orphan the first
    // pool/attractor/emitter allocations.
    HS_CHECK(!pool.is_bound(), "ParticleSystem::init called twice");
    this->friction = friction;
    this->gravity = gravity;
    // Clamp before narrowing: a float outside [0, 65535] is UB cast to uint16_t.
    this->max_life =
        static_cast<uint16_t>(std::min(std::max(max_life, 0.0f), 65535.0f));
    active_count = 0;
    pool.bind(arena, CAPACITY);
    for (size_t i = 0; i < CAPACITY; ++i) {
      pool.emplace_back();
    }
    attractors.bind(arena, ATTRACTOR_CAP);
    emitters.bind(arena, EMITTER_CAP);
  }

  /**
   * @brief Registers a per-frame emitter functor.
   * @param fn The emitter to add.
   * @note Traps on overflow: emitters are registered at setup with fixed
   * cardinality, so an overrun is a bug — unlike spawn()'s runtime soft-drop.
   */
  void add_emitter(EmitterFn fn) { emitters.push_back(fn); }

  /**
   * @brief Adds a point attractor to the system.
   * @param pos World-space center of the attractor.
   * @param str Attractive force multiplier.
   * @param kill Kill radius (particles within it die).
   * @param horizon Event-horizon radius (steering becomes radial within it).
   * @note Traps on overflow: attractors are registered at setup with fixed
   * cardinality, so an overrun is a bug — unlike spawn()'s runtime soft-drop.
   */
  void add_attractor(const Vector &pos, float str, float kill, float horizon) {
    attractors.push_back({pos, str, kill, horizon});
  }

  /**
   * @brief Spawns a new particle with a color seed.
   * @param pos Initial position.
   * @param vel Initial velocity.
   * @param seed Color seed for palette offset.
   * @details A full pool is a designed, non-fatal condition (drop the spawn, keep
   * rendering); log it like Timeline::add rather than dropping it silently.
   */
  void spawn(const Vector &pos, const Vector &vel, uint16_t seed) {
    HS_CHECK(pool.is_bound(), "ParticleSystem::spawn before init");
    if (active_count < pool.capacity()) {
      pool[active_count++].init(pos, vel, seed, max_life);
    } else {
      hs::log("ParticleSystem pool full, dropping spawn");
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
      for (size_t i = 0; i < emitters.size(); ++i) {
        emitters[i](*this);
      }

      float max_delta = (2 * PI_F) / W;

      // Swap-remove dead particles, re-testing the same index. The i-- relies on
      // unsigned wraparound (i==0 -> SIZE_MAX -> ++i back to 0): keep i unsigned
      // and do not read i between the decrement and the loop's ++i.
      for (size_t i = 0; i < active_count; ++i) {
        bool dead = step_particle(pool[i], max_delta);
        if (dead) {
          if (i != static_cast<size_t>(active_count - 1)) {
            pool[i] = pool[active_count - 1];
          }
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
   * @details Ages it, drags the carried velocity, applies attractor
   * gravity/steering, rotates position+velocity along the surface, and updates
   * its trail.
   *
   * Integration is forward Euler with an implicit dt = 1 (one step == one
   * frame): forces and the rotation step are applied per-frame with no dt
   * scaling, so the motion is frame-rate dependent — the same emitter/attractor
   * config evolves faster or slower if the step cadence changes. This is
   * intentional for a fixed-cadence display driver; it is NOT a wall-clock
   * physics sim.
   *
   * Ordering: `friction` damps the velocity CARRIED IN from previous frames
   * BEFORE this frame's attractor impulse is added (v <- friction*v + impulse),
   * so a force applied this frame is not also damped the same frame; the drag
   * still precedes the move, which uses the post-impulse velocity. Per-frame
   * impulse magnitude therefore scales with the attractor `strength` alone, no
   * longer with `friction`. Effect presets carry the matching strength (e.g.
   * MindSplatter's Well Str is pre-scaled by its friction) so the shipped look
   * is preserved across this correction.
   */
  bool step_particle(Particle<TRAIL_LEN> &p, float max_delta) {
    bool active = p.life > 0;
    if (active) {
      p.life--;
      active = p.life > 0;
    }

    if (active) {
      Vector pos = p.position;

      // Drag the carried velocity before this frame's impulse so a fresh impulse
      // is not also damped this frame (forward Euler: v <- friction*v + impulse).
      p.velocity *= friction;

      for (size_t k = 0; k < attractors.size(); ++k) {
        const auto &attr = attractors[k];
        float dist_sq = distance_squared(pos, attr.position);

        if (dist_sq < attr.kill_radius * attr.kill_radius) {
          p.life = 0; // Stay dead; a live `life` resurrects the particle next frame.
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
            // Gravity. pos and the attractor can be ~antipodal (undefined cross
            // axis), so guard the normalize.
            float force = (gravity * attr.strength) / dist_sq;
            Vector torque =
                normalized_or(cross(pos, attr.position), Vector(1, 0, 0)) * force;
            p.velocity += cross(torque, pos);
          }
        }
      }

      if (active) {
        // The surface-rotation axis cross(pos, velocity) vanishes for a purely
        // radial velocity (no motion along the sphere), so skip rather than
        // normalize a zero-length axis.
        float speed = p.velocity.magnitude();
        Vector axis = cross(pos, p.velocity);
        if (speed > 0.000001f && axis.magnitude() > 0.000001f) {
          // Cap the per-frame surface advance at one column to avoid aliasing;
          // velocity keeps its full magnitude.
          speed = std::min(speed, max_delta);
          Quaternion dq = make_rotation(axis.normalized(), speed);
          p.position = rotate(p.position, dq);
          p.velocity = rotate(p.velocity, dq);
        }
      }
    }

    if (active) {
      p.history.record(p.position);
    } else {
      if (p.history.length() > 0) {
        p.history.expire();
      }
    }

    // Reclaimable only once dead AND the trail has drained, so a killed particle
    // holds its slot for up to TRAIL_LEN more frames — size CAPACITY accordingly.
    return !active && p.history.length() == 0;
  }
};

} // namespace Animation
