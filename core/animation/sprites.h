/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

namespace Animation {

/** Per-window frame cap that keeps a Sprite's fade-in + fade-out sum in int. */
static constexpr int MAX_FADE_DURATION = 1 << 24;

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
   * @param fade_in_duration Frames for fading in; in [0, MAX_FADE_DURATION].
   * @param fade_in_easing_fn Easing for fade-in.
   * @param fade_out_duration Frames for fading out; in [0, MAX_FADE_DURATION].
   * @param fade_out_easing_fn Easing for fade-out.
   * @param paused Optional pause gate; null = always runs.
   */
  Sprite(SpriteFn draw_fn, int duration, int fade_in_duration = 0,
         EasingFn fade_in_easing_fn = ease_linear, int fade_out_duration = 0,
         EasingFn fade_out_easing_fn = ease_linear,
         const bool *paused = nullptr)
      : AnimationBase(duration, false), draw_fn(std::move(draw_fn)),
        fade_in_duration(fade_in_duration),
        fade_out_duration(fade_out_duration),
        fade_in_easing(std::move(fade_in_easing_fn)),
        fade_out_easing(std::move(fade_out_easing_fn)), paused(paused) {
    HS_CHECK(fade_in_duration >= 0 && fade_in_duration <= MAX_FADE_DURATION,
             "Sprite fade-in duration must be in [0, MAX_FADE_DURATION]");
    HS_CHECK(fade_out_duration >= 0 && fade_out_duration <= MAX_FADE_DURATION,
             "Sprite fade-out duration must be in [0, MAX_FADE_DURATION]");
    // Overlapping windows (durations are independent GUI sliders): scale both
    // fades proportionally to fit the visible duration, so the envelope still
    // peaks at full opacity and stays a continuous triangle (definite sprites
    // only; indefinite ones skip fade-out).
    // The base promotes duration 0 to 1, so read the member rather than the
    // parameter.
    int fade_total = this->fade_in_duration + this->fade_out_duration;
    if (this->duration >= 0 && fade_total > this->duration) {
      const int requested_fade_in = this->fade_in_duration;
      this->fade_in_duration = static_cast<long long>(this->duration) *
                               this->fade_in_duration / fade_total;
      // Keep a requested fade-in visible when the integer scale floors it to 0.
      if (this->fade_in_duration == 0 && requested_fade_in > 0 &&
          this->duration >= 2)
        this->fade_in_duration = 1;
      this->fade_out_duration = this->duration - this->fade_in_duration;
    }
    // A duration-1 sprite's sole frame is t==duration, where the fade-out ramp
    // reaches 0; drop fade-out so that frame isn't drawn fully transparent.
    if (this->duration == 1)
      this->fade_out_duration = 0;
  }

  /**
   * @brief Steps the animation, computes the current opacity inline, and calls
   * the draw function.
   * @param canvas The canvas buffer passed to the draw function.
   */
  void step(Canvas &canvas) override {
    // Paused: hold the frame (don't advance the timer) but keep drawing at the
    // current opacity.
    if (!is_paused(paused))
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
      float elapsed =
          static_cast<float>(t + static_cast<uint32_t>(fade_out_duration) -
                             static_cast<uint32_t>(duration));
      float progress = elapsed / fade_out_duration;
      fade_out = 1.0f - fade_out_easing(hs::clamp(progress, 0.0f, 1.0f));
    }

    // Overshooting easings (ease_out_elastic) exceed 1, so the fade-out ramp
    // goes negative.
    draw_fn(canvas, std::max(std::min(fade_in, fade_out), 0.0f));
  }

private:
  SpriteFn draw_fn;             /**< The drawing function functor. */
  int fade_in_duration;         /**< Duration of fade-in phase in frames. */
  int fade_out_duration;        /**< Duration of fade-out phase in frames. */
  EasingFn fade_in_easing;      /**< Easing curve for fade-in. */
  EasingFn fade_out_easing;     /**< Easing curve for fade-out. */
  const bool *paused = nullptr; /**< Optional pause gate; holds the frame (no
                                    timer advance, still draws) when set. */
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

  /** Trail of world-space positions, snorm16-quantized (unit-sphere domain). */
  QuantizedVectorTrail<TRAIL_LEN> history;

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

static constexpr float ATTRACTOR_MIN_DISTANCE_SQ = 0.0000001f;

/**
 * @brief A point attractor influencing nearby particles.
 */
struct Attractor {
  Vector position;     /**< World-space center of the attractor. */
  float strength;      /**< Attractive force multiplier. */
  float kill_radius;   /**< Radius within which particles are killed. */
  float event_horizon; /**< Radius within which steering becomes radial. */
};

/**
 * @brief Particle/attractor geometry precomputed by the caller.
 */
struct AttractorSample {
  float pos_sq;   /**< dot(pos, pos). */
  float dot_pa;   /**< dot(pos, attractor position). */
  float cross_sq; /**< Squared magnitude of cross(pos, attractor position). */
  float dist_sq;  /**< Squared distance from pos to the attractor. */
};

[[maybe_unused]] HS_COLD static bool
apply_signed_axis_attractor(uint16_t &life, Vector &velocity, const Vector &pos,
                            float max_delta, float gravity,
                            const Attractor &attractor, AttractorSample s) {
  if (s.dist_sq < attractor.kill_radius * attractor.kill_radius) {
    life = 0;
    return false;
  }

  if (s.dist_sq <= ATTRACTOR_MIN_DISTANCE_SQ)
    return true;

  if (s.dist_sq < attractor.event_horizon * attractor.event_horizon) {
    const Vector torque = (attractor.position - pos).normalized();
    const float speed = std::max(velocity.magnitude(), max_delta);
    velocity = torque * speed;
    return true;
  }

  const float force = (gravity * attractor.strength) / s.dist_sq;
  if (s.cross_sq < math::EPS_NORMALIZE_SQ) {
    velocity += cross(Vector(force, 0, 0), pos);
  } else {
    const Vector tangent = attractor.position * s.pos_sq - pos * s.dot_pa;
    velocity += tangent * (force / sqrtf(s.cross_sq));
  }
  return true;
}

/**
 * @brief A physics-based particle system with emitters and attractors.
 *
 * PERPETUAL (duration -1, no repeat): never reaches done(), so a `.then()`
 * callback NEVER fires. Drive follow-on behavior from a finite animation or
 * cancel() it explicitly.
 *
 * @tparam W Denominator of the per-frame surface-advance cap 2*PI/W radians,
 *        i.e. the display width in columns.
 * @tparam CAPACITY Maximum number of particles in the pool.
 * @tparam TRAIL_LEN Trail length per particle.
 * @tparam EMITTER_CAP Maximum number of emitters.
 * @tparam ATTRACTOR_CAP Maximum number of attractors.
 * @tparam SIGNED_AXIS_ATTRACTORS Enable paired unit signed-axis attractor
 *        algebra.
 */
template <int W, int CAPACITY, int TRAIL_LEN = 8, int EMITTER_CAP = 8,
          int ATTRACTOR_CAP = 8, bool SIGNED_AXIS_ATTRACTORS = false>
class ParticleSystem
    : public AnimationBase<
          ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP,
                         SIGNED_AXIS_ATTRACTORS>> {
public:
  static_assert(CAPACITY <= 65535,
                "active_count is uint16_t; CAPACITY must fit in it");
  static_assert(!SIGNED_AXIS_ATTRACTORS || ATTRACTOR_CAP == 6,
                "paired signed-axis physics requires six attractors");

  ArenaVector<Particle<TRAIL_LEN>> pool; /**< Backing pool of particles. */

  /**
   * @brief Number of live particles in the pool prefix.
   * @return Count of currently active particles.
   */
  uint16_t active() const { return active_count; }

  float friction = 0.85f;  /**< Per-frame velocity damping factor. */
  float gravity = 0.001f;  /**< Base gravitational constant for attractors. */
  uint16_t max_life = 600; /**< Default particle lifetime in frames. */

  using Attractor = Animation::Attractor; /**< Point attractor. */

  using EmitterFn =
      Fn<void(ParticleSystem &), 32>; /**< Per-frame emitter functor. */

  ArenaVector<Attractor> attractors; /**< Active attractors. */
  ArenaVector<EmitterFn> emitters;   /**< Active emitters. */
#ifdef HS_TEST_BUILD
  bool reference_signed_axis_physics = false;
#endif

  /**
   * @brief Constructs an indefinite, non-repeating particle system.
   */
  ParticleSystem()
      : AnimationBase<ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP,
                                     ATTRACTOR_CAP, SIGNED_AXIS_ATTRACTORS>>(
            -1, false) {}

  /**
   * @brief Binds the pool, attractor, and emitter vectors and pre-constructs
   * inactive particles.
   * @param arena Arena providing backing storage for all vectors.
   * @param friction Per-frame velocity damping factor.
   * @param gravity Base gravitational constant for attractors.
   * @param max_life Default particle lifetime in frames; must be finite and in
   *        [1, 65535].
   * @details Pre-constructs CAPACITY inactive particles in the pool.
   */
  void init(Arena &arena, float friction = 0.85f, float gravity = 0.001f,
            float max_life = 600.0f) {
    // The arena has no per-allocation free, so re-binding would orphan the first
    // pool/attractor/emitter allocations.
    HS_CHECK(!pool.is_bound(), "ParticleSystem::init called twice");
    HS_CHECK(std::isfinite(max_life) && max_life >= 1.0f &&
                 max_life <= 65535.0f,
             "ParticleSystem max_life must be finite and in [1, 65535]");
    this->friction = friction;
    this->gravity = gravity;
    this->max_life = static_cast<uint16_t>(max_life);
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
    if constexpr (SIGNED_AXIS_ATTRACTORS) {
      static constexpr std::array<Vector, 6> AXES = {X_AXIS,  -X_AXIS, Y_AXIS,
                                                     -Y_AXIS, Z_AXIS,  -Z_AXIS};
      const size_t index = attractors.size();
      HS_CHECK(index < AXES.size() && pos.x == AXES[index].x &&
                   pos.y == AXES[index].y && pos.z == AXES[index].z,
               "ParticleSystem signed-axis attractors must be registered "
               "+X,-X,+Y,-Y,+Z,-Z");
    }
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
    AnimationBase<
        ParticleSystem<W, CAPACITY, TRAIL_LEN, EMITTER_CAP, ATTRACTOR_CAP,
                       SIGNED_AXIS_ATTRACTORS>>::step(canvas);

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
  static constexpr float MOTION_MIN_SPEED = 0.000001f;

  uint16_t active_count =
      0; /**< Number of live particles in the pool prefix. */

  /**
   * @brief Applies every registered attractor to one particle.
   * @param p The particle whose velocity (and life, on a kill) is updated.
   * @param pos The particle's world-space position.
   * @param max_delta Per-frame surface-rotation cap (radians).
   * @return False once an attractor killed the particle.
   */
  __attribute__((always_inline)) bool
  apply_attractors(Particle<TRAIL_LEN> &p, const Vector &pos, float max_delta) {
    for (size_t k = 0; k < attractors.size(); ++k) {
      const Attractor &attr = attractors[k];
      float dist_sq = distance_squared(pos, attr.position);

      if (dist_sq < attr.kill_radius * attr.kill_radius) {
        // Stay dead; a live `life` resurrects the particle next frame.
        p.life = 0;
        return false;
      }

      if (dist_sq > ATTRACTOR_MIN_DISTANCE_SQ) {
        if (dist_sq < attr.event_horizon * attr.event_horizon) {
          // Steer into center. Floor the speed so a friction-drained particle
          // still advances inward to kill_radius instead of stalling.
          Vector torque = (attr.position - pos).normalized();
          float speed = std::max(p.velocity.magnitude(), max_delta);
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
    return true;
  }

  /**
   * @brief Advances one particle on the sphere surface for one frame.
   * @param p The particle to advance (modified in place).
   * @param max_delta Per-frame surface-rotation cap (radians), one display
   * column wide; the move below clamps to it so a fast particle never jumps more
   * than one column per frame (trail/motion-blur aliasing).
   * @return True once the particle is dead AND its trail has fully drained (so
   * the caller can remove it).
   * @details Ages, drags velocity, applies attractor gravity/steering, rotates
   * position+velocity along the surface, and updates the trail.
   *
   * Integration is forward Euler with implicit dt = 1 (one step == one frame),
   * so motion is frame-rate dependent — intentional for a fixed-cadence display
   * driver, not a wall-clock physics sim.
   *
   * Ordering: friction damps the carried-in velocity BEFORE this frame's impulse
   * (v <- friction*v + impulse), so the impulse is not damped the same frame.
   * Impulse magnitude scales with the attractor `strength` alone; presets
   * pre-scale their strength to match (e.g. MindSplatter's Well Str).
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

      if constexpr (SIGNED_AXIS_ATTRACTORS) {
#ifdef HS_TEST_BUILD
        if (!reference_signed_axis_physics) {
#endif
          assert(attractors.size() == 6);
          const float pos_sq = dot(pos, pos);
          const float coordinates[] = {pos.x, pos.y, pos.z};
          for (size_t pair = 0; pair < 3; ++pair) {
            const auto &plus = attractors[pair * 2];
            const auto &minus = attractors[pair * 2 + 1];
            const float q = coordinates[pair];
            float cross_sq;
            if (pair == 0)
              cross_sq = pos.y * pos.y + pos.z * pos.z;
            else if (pair == 1)
              cross_sq = pos.x * pos.x + pos.z * pos.z;
            else
              cross_sq = pos.x * pos.x + pos.y * pos.y;
            const float plus_delta = q - 1.0f;
            const float minus_delta = q + 1.0f;
            const float dist_plus_sq = plus_delta * plus_delta + cross_sq;
            const float dist_minus_sq = minus_delta * minus_delta + cross_sq;
            const float plus_kill_sq = plus.kill_radius * plus.kill_radius;
            const float minus_kill_sq = minus.kill_radius * minus.kill_radius;
            const float plus_horizon_sq =
                plus.event_horizon * plus.event_horizon;
            const float minus_horizon_sq =
                minus.event_horizon * minus.event_horizon;
            const bool common_far = dist_plus_sq > ATTRACTOR_MIN_DISTANCE_SQ &&
                                    dist_minus_sq > ATTRACTOR_MIN_DISTANCE_SQ &&
                                    dist_plus_sq >= plus_horizon_sq &&
                                    dist_minus_sq >= minus_horizon_sq &&
                                    dist_plus_sq >= plus_kill_sq &&
                                    dist_minus_sq >= minus_kill_sq;

            if (!common_far) {
              if (!apply_signed_axis_attractor(
                      p.life, p.velocity, pos, max_delta, gravity, plus,
                      {pos_sq, q, cross_sq, dist_plus_sq}) ||
                  !apply_signed_axis_attractor(
                      p.life, p.velocity, pos, max_delta, gravity, minus,
                      {pos_sq, -q, cross_sq, dist_minus_sq})) {
                active = false;
                break;
              }
              continue;
            }

            const float inv_pair = 1.0f / (dist_plus_sq * dist_minus_sq);
            const float inv_plus = dist_minus_sq * inv_pair;
            const float inv_minus = dist_plus_sq * inv_pair;
            if (cross_sq < math::EPS_NORMALIZE_SQ) {
              const float force = gravity * (plus.strength * inv_plus +
                                             minus.strength * inv_minus);
              p.velocity += cross(Vector(force, 0, 0), pos);
            } else {
              Vector tangent(-pos.x * q, -pos.y * q, -pos.z * q);
              if (pair == 0)
                tangent.x += pos_sq;
              else if (pair == 1)
                tangent.y += pos_sq;
              else
                tangent.z += pos_sq;
              const float inv_cross = 1.0f / sqrtf(cross_sq);
              const float scale =
                  gravity * inv_cross *
                  (plus.strength * inv_plus - minus.strength * inv_minus);
              p.velocity += tangent * scale;
            }
          }
#ifdef HS_TEST_BUILD
        } else {
          active = apply_attractors(p, pos, max_delta);
        }
#endif
      } else {
        active = apply_attractors(p, pos, max_delta);
      }

      if (active) {
        if constexpr (SIGNED_AXIS_ATTRACTORS) {
          const float speed_sq = dot(p.velocity, p.velocity);
          Vector axis = cross(pos, p.velocity);
          const float axis_sq = dot(axis, axis);
          constexpr float MOTION_MIN_SQ = MOTION_MIN_SPEED * MOTION_MIN_SPEED;
          if (speed_sq > MOTION_MIN_SQ && axis_sq > MOTION_MIN_SQ) {
            const float speed = std::min(sqrtf(speed_sq), max_delta);
            axis *= 1.0f / sqrtf(axis_sq);
            const Quaternion dq = make_rotation(axis, speed);
            p.position = rotate(p.position, dq);
            p.velocity = rotate(p.velocity, dq);
          }
        } else {
          // The surface-rotation axis cross(pos, velocity) vanishes for a purely
          // radial velocity (no motion along the sphere), so skip rather than
          // normalize a zero-length axis.
          float speed = p.velocity.magnitude();
          Vector axis = cross(pos, p.velocity);
          if (speed > MOTION_MIN_SPEED && axis.magnitude() > MOTION_MIN_SPEED) {
            // Cap the per-frame surface advance at one column to avoid aliasing;
            // velocity keeps its full magnitude.
            speed = std::min(speed, max_delta);
            Quaternion dq = make_rotation(axis.normalized(), speed);
            p.position = rotate(p.position, dq);
            p.velocity = rotate(p.velocity, dq);
          }
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
