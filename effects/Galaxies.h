/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Galaxies.h
 * @brief Six spiral galaxies, one on each octahedron vertex, built from
 *        particles orbiting central attractors.
 */

#include "core/animation/orientation.h"
#include "core/engine/engine.h"

namespace hs_test {
namespace effects_tests {
struct GalaxiesWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Spiral galaxies on the vertices of an octahedron.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Stars form along rotating logarithmic arms, then orbit a central
 *          attractor. A rotating brightness pattern lights young stars within
 *          the arms; older stars fade into a dim stellar disk.
 */
template <int W, int H> class Galaxies : public Effect {
public:
  static constexpr const char *EFFECT_ID = "Galaxies";

  /**
   * @brief Constructs the effect.
   */
  HS_COLD_MEMBER Galaxies()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})) {}

  /**
   * @brief Registers params, builds the particle system,
   *        and starts the orientation drift.
   */
  HS_COLD_MEMBER void init() override {
    static constexpr size_t SCRATCH_BYTES = 6 * 1024;
    // Scan::Point's row-interval buffers for the core bulges.
    static constexpr size_t SCRATCH_B_BYTES = 2 * 1024;
    ArenaSplit{SCRATCH_BYTES, SCRATCH_B_BYTES}.configure();

    // GLOBAL_ARENA_SIZE is inflated on host; check the device arena.
    static constexpr size_t POOL_BYTES =
        sizeof(Animation::PointParticle) * NUM_PARTICLES;
    static constexpr size_t AUX_RESERVE_BYTES = 6 * 1024;
    static_assert(
        POOL_BYTES + AUX_RESERVE_BYTES <=
            ArenaSplit{SCRATCH_BYTES, SCRATCH_B_BYTES}.device_persistent(),
        "Galaxies particle pool + attractor/emitter storage "
        "overflow the device persistent arena");

    register_param("Friction", &params.friction, 0.999f, 1.0f);
    register_param("Core Mass", &params.core_mass, 0.02f, 0.2f);
    register_param("Black Hole", &params.black_hole);
    register_param("Orbit Spd", &params.orbit_speed, 0.004f, 0.06f);
    register_param("Arm Spin", &params.arm_spin, 0.0f, 0.25f);
    register_param("Arm Pitch", &params.arm_pitch, 0.1f, 0.65f);
    register_param("Arm Sharpness", &params.arm_sharpness, 0.0f, 3.0f);
    register_int_param("Arms", &params.arms, 1, MAX_ARMS);
    register_param("Emission Rate", &params.emission_rate, MIN_EMISSION_RATE,
                   MAX_EMISSION_RATE);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Particles", &params.active_count,
                   ParamSpec<float>{.min = 0.0f,
                                    .max = static_cast<float>(NUM_PARTICLES),
                                    .readonly = true});

    build_particle_system();

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, math::Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /**
   * @brief Steps the timeline and the particle system, then renders.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(gx_timeline_step);
      timeline.step(canvas);
    }

    particle_system.friction = params.friction;
    particle_system.motion_cap =
        std::max((2.0f * math::PI_F) / W, params.orbit_speed * 3.0f);
    for (size_t i = 0; i < particle_system.attractors.size(); ++i)
      particle_system.attractors[i].strength = params.core_mass;
    {
      HS_PROFILE(gx_particle_step);
      particle_system.step(canvas);
    }
    params.active_count = static_cast<float>(particle_system.active());

    // Alpha below one slider LSB: skip rasterizing; the physics still runs.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;
    draw_particles(canvas);
  }

private:
  friend struct ::hs_test::effects_tests::GalaxiesWhiteBox;

  /** @brief Number of galaxies: one per octahedron vertex. */
  static constexpr int NUM_GALAXIES = Solids::Octahedron::NUM_VERTS;
  static constexpr int MAX_ARMS = 8;

  /** @brief Point particles retain no history. */
  static constexpr int TRAIL_LEN = 0;
  /**
   * @brief Fixed particle pool capacity.
   * @details The emitters stop spawning while the pool is full.
   */
  static constexpr int NUM_PARTICLES = 12000;

  using ParticleSystem =
      Animation::ParticleSystem<W, NUM_PARTICLES, TRAIL_LEN, NUM_GALAXIES,
                                NUM_GALAXIES, true, 1,
                                Animation::PointParticle>;

  static constexpr float MIN_EMISSION_RATE = 0.15f;
  static constexpr float MAX_EMISSION_RATE = 4.0f;
  static constexpr float GRAVITY = 0.001f;
  static constexpr float PARTICLE_LIFETIME_FRAMES = 800.0f;
  static constexpr float REFERENCE_ORBIT_SPEED = 0.01285f;
  /** @brief Outer radius of the stellar disk (radians). */
  static constexpr float RING_RADIUS = 0.58f;
  static constexpr float INNER_RADIUS = 0.14f;
  /** @brief Core kill radius (chord) for particles that plunge inward. */
  static constexpr float KILL_RADIUS = 0.008f;
  /** @brief Galaxies do not use the attractor's radial steering zone. */
  static constexpr float EVENT_HORIZON = 0.0f;
  /** @brief Softening radius that limits the core force (chord). */
  static constexpr float SOFTENING_RADIUS = 0.02f;
  /** @brief Radius of the particle fade into each core (radians). */
  static constexpr float HOLE_FADE_RADIUS = 0.035f;
  /** @brief Spawn-angle jitter half-width (radians), thickens the arms. */
  static constexpr float ARM_JITTER = 0.22f;
  /** @brief Orbital speed jitter half-width, as a fraction of the speed. */
  static constexpr float SPEED_JITTER = 0.008f;
  static constexpr float ARM_FADE_START = 300.0f;
  static constexpr float ARM_FADE_END = 650.0f;
  /** @brief Frames over which a new particle fades in. */
  static constexpr float FADE_IN_FRAMES = 2.0f;
  /** @brief Frames over which an expiring particle fades out. */
  static constexpr float FADE_OUT_FRAMES = 20.0f;
  /** @brief Core distance inside which every particle starts to dim (radians). */
  static constexpr float HEAD_FADE_RADIUS = 0.03f;
  /**
   * @brief Radius of the glowing bulge drawn over each core (radians).
   * @details At least 1.5 columns, so it stays visible at low resolution.
   */
  static constexpr float BULGE_RADIUS =
      std::max(0.08f, 1.5f * math::RADIANS_PER_COLUMN<W>);

  /** @brief Spawn geometry and arm state for one galaxy. */
  struct Galaxy {
    math::Vector core;     /**< Unit vector to the attractor. */
    math::Vector u;        /**< First tangent axis at the core. */
    math::Vector w;        /**< Second tangent axis at the core. */
    float phase;           /**< Arm angle (radians, wrapped to [0, 2pi)). */
    float spin;            /**< +1 or -1: orbit direction about the core. */
    float emission_credit; /**< Fractional particles waiting to spawn. */
    uint8_t arm;           /**< Arm the next spawn belongs to. */
  };

  /**
   * @brief Effect parameters and read-only telemetry.
   */
  struct Params {
    float friction = 0.99965f;   /**< Velocity retention per frame. */
    float core_mass = 0.10478f;  /**< Attractor strength. */
    float orbit_speed = 0.0134f; /**< Reference spawn speed (radians/frame). */
    float arm_spin = 0.028f;     /**< Arm rotation (radians/frame). */
    float arm_pitch = 0.2f;      /**< Spiral pitch angle (radians). */
    float arm_sharpness = 0.8f;  /**< Arm lobe narrowing; 1 is a cosine lobe. */
    int arms = 2;                /**< Arms per galaxy. */
    float emission_rate = 2.5f;  /**< Particles per galaxy per frame. */
    float alpha = 1.0f;          /**< Overall opacity. */
    float active_count = 0.0f;   /**< Live particles (engine-written). */
    bool black_hole = false;     /**< Black core. */
  } params;

  /**
   * @brief Places the attractors and emitters on the octahedron vertices.
   * @details Single-shot: ParticleSystem::init traps on a second call.
   */
  HS_COLD_MEMBER void build_particle_system() {
    particle_system.init(persistent_arena, params.friction, GRAVITY,
                         PARTICLE_LIFETIME_FRAMES);
    for (int i = 0; i < NUM_GALAXIES; ++i) {
      Galaxy &g = galaxies[i];
      const math::Basis basis =
          math::make_basis(math::Quaternion(), Solids::Octahedron::vertices[i]);
      // Signed-axis physics requires the exact axis, not a renormalized one.
      g.core = Solids::Octahedron::vertices[i];
      g.u = basis.u;
      g.w = basis.w;
      g.phase = hs::rand_f(0.0f, 2.0f * math::PI_F);
      g.spin = (i & 1) ? -1.0f : 1.0f;
      g.emission_credit = hs::rand_f();
      g.arm = 0;
      particle_system.add_attractor(g.core, params.core_mass, KILL_RADIUS,
                                    EVENT_HORIZON, SOFTENING_RADIUS);
      // EmitterFn's inline capture is too small for a Galaxy; index by i.
      particle_system.add_emitter([this, i](ParticleSystem &) {
        Galaxy &galaxy = galaxies[i];
        galaxy.phase = fmodf(galaxy.phase + params.arm_spin * galaxy.spin,
                             2.0f * math::PI_F);
        if (galaxy.phase < 0.0f)
          galaxy.phase += 2.0f * math::PI_F;
        galaxy.emission_credit += hs::clamp(
            params.emission_rate, MIN_EMISSION_RATE, MAX_EMISSION_RATE);
        while (galaxy.emission_credit >= 1.0f) {
          galaxy.emission_credit -= 1.0f;
          emit(galaxy, i);
        }
      });
    }
  }

  /**
   * @brief Spawns one particle on the next arm of a galaxy.
   * @param g Galaxy to spawn into; its arm advances.
   * @param index Galaxy index, stored in the color seed's low byte.
   */
  HS_FLASH_MEMBER void emit(Galaxy &g, int index) {
    const int arms = hs::clamp(params.arms, 1, MAX_ARMS);
    g.arm = static_cast<uint8_t>((g.arm + 1) % arms);
    // Skip the spawn rather than let spawn() log a dropped particle.
    if (particle_system.active() >= particle_system.pool.capacity())
      return;

    const uint16_t alpha_seed = static_cast<uint16_t>(hs::rand_f() * 255.0f);
    const float spread =
        1.0f - 0.9f * static_cast<float>(alpha_seed) * (1.0f / 255.0f);
    const float ring = hs::rand_f(INNER_RADIUS, RING_RADIUS);
    const float angle =
        g.phase + 2.0f * math::PI_F * g.arm / arms -
        g.spin * logf(ring / RING_RADIUS) / tanf(params.arm_pitch) +
        hs::rand_f(-ARM_JITTER, ARM_JITTER) * spread;
    const math::Vector radial =
        (g.u * cosf(angle) + g.w * sinf(angle)).normalized();
    const math::Vector pos =
        (g.core * cosf(ring) + radial * sinf(ring)).normalized();
    // core and radial are orthonormal, so their cross is the unit tangent.
    const math::Vector tangent = math::cross(g.core, radial);
    const math::Vector outward = math::cross(tangent, pos);
    const float speed = circular_orbit_speed(pos, outward, ring) *
                        (params.orbit_speed / REFERENCE_ORBIT_SPEED) *
                        (1.0f + hs::rand_f(-SPEED_JITTER, SPEED_JITTER));
    const uint16_t color_seed =
        static_cast<uint16_t>(index | (alpha_seed << 8));
    particle_system.spawn(pos, tangent * (speed * g.spin), color_seed);
  }

  /**
   * @brief Brightness pattern that stars move through while orbiting.
   * @details A cosine lobe per arm, narrowed by `Params::arm_sharpness`.
   */
  float arm_density(const math::Vector &position, const Galaxy &galaxy,
                    float cos_distance, float winding) const {
    const float ring = math::fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f));
    const float azimuth = math::fast_atan2(math::dot(position, galaxy.w),
                                           math::dot(position, galaxy.u));
    const float angle =
        azimuth - galaxy.phase +
        galaxy.spin * winding * logf(fmaxf(ring, KILL_RADIUS) / RING_RADIUS);
    const float lobe = angle * hs::clamp(params.arms, 1, MAX_ARMS);
    const float offset =
        lobe - math::TWO_PI_F * floorf(lobe * (1.0f / math::TWO_PI_F) + 0.5f);
    return fmaxf(math::fast_cosf(fminf(fabsf(offset) * params.arm_sharpness,
                                       0.5f * math::PI_F)),
                 0.0f);
  }

  /** @brief Young stars light the arms; old stars form a faint stellar disk. */
  static float particle_alpha(uint16_t color_seed, float age, float arm) {
    const float u = static_cast<float>(color_seed >> 8) * (1.0f / 255.0f);
    const float young = hs::clamp(
        (ARM_FADE_END - age) / (ARM_FADE_END - ARM_FADE_START), 0.0f, 1.0f);
    const float background = 0.006f + 0.008f * u;
    const float broad = arm * arm;
    const float spine = broad * broad * broad * broad;
    const float profile = broad + (spine - broad) * u;
    return background + (0.5f + 0.5f * u * u - background) * young * profile;
  }

  /** @brief Stable white stars with sparse blue, red and pale-yellow accents. */
  static Color4 star_color(uint16_t color_seed) {
    const uint8_t tint = static_cast<uint8_t>((color_seed >> 8) * 73u +
                                              (color_seed & 0xff) * 29u + 41u);
    if (tint < 20)
      return Color4(155, 200, 255);
    if (tint < 40)
      return Color4(255, 150, 140);
    if (tint < 50)
      return Color4(255, 245, 205);
    return Color4(255, 255, 255);
  }

  /** @brief Circular speed from the net inward pull of all six cores. */
  float circular_orbit_speed(const math::Vector &pos,
                             const math::Vector &outward, float ring) const {
    float inward_acceleration = 0.0f;
    for (const auto &attractor : particle_system.attractors) {
      const float cos_distance = math::dot(pos, attractor.position);
      const math::Vector toward = attractor.position - pos * cos_distance;
      const float tangent_length = toward.magnitude();
      const float dist_sq = math::distance_squared(pos, attractor.position);
      if (dist_sq > Animation::ATTRACTOR_MIN_DISTANCE_SQ &&
          tangent_length > math::EPS_NORMALIZE_SQ)
        inward_acceleration -=
            particle_system.gravity * attractor.strength *
            math::dot(toward, outward) /
            ((dist_sq + attractor.softening_sq) * tangent_length);
    }
    return sqrtf(fmaxf(inward_acceleration * tanf(ring), 0.0f));
  }

  /**
   * @brief Renders particles, then glowing or black-hole cores.
   * @param canvas Target canvas.
   */
  void draw_particles(Canvas &canvas) {
    HS_PROFILE(gx_draw_particles);
    const math::RotationMatrix rotation(orientation.get());
    const float cos_hole = math::fast_cosf(HOLE_FADE_RADIUS);
    const float inv_head_span =
        1.0f / (1.0f - math::fast_cosf(HEAD_FADE_RADIUS));
    const float max_life = static_cast<float>(particle_system.max_life);
    const float alpha = params.alpha;
    const float winding = 1.0f / tanf(params.arm_pitch);

    filters.prepare(canvas);
    for (int i = 0; i < particle_system.active(); ++i) {
      const auto &p = particle_system.pool[i];
      const Galaxy &galaxy = galaxies[p.color_seed & 0xff];
      const math::Vector position = p.get_position();
      const float cos_distance = math::dot(position, galaxy.core);
      const float hole =
          cos_distance < cos_hole
              ? 1.0f
              : math::quintic_kernel(
                    math::fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f)) /
                    HOLE_FADE_RADIUS);
      const float head_fade =
          math::quintic_kernel((1.0f - cos_distance) * inv_head_span);
      const float age = max_life - static_cast<float>(p.life);
      const float fade_in = hs::clamp(age / FADE_IN_FRAMES, 0.0f, 1.0f);
      const float fade_out =
          hs::clamp(static_cast<float>(p.life) / FADE_OUT_FRAMES, 0.0f, 1.0f);
      const float arm = age < ARM_FADE_END ? arm_density(position, galaxy,
                                                         cos_distance, winding)
                                           : 0.0f;
      Color4 c = star_color(p.color_seed);
      c.alpha *= hole * head_fade * fade_in * fade_out *
                 particle_alpha(p.color_seed, age, arm) * alpha;
      filters.plot(canvas, rotation.apply(position), c.color, 0.0f, c.alpha);
    }

    // Scan::Point leaves its quintic coverage in v2.
    const Color4 core_color(255, 255, 255);
    auto bulge_shader = [&](const math::Vector &, Fragment &f) {
      Color4 c = core_color;
      c.alpha *= hs::clamp(f.v2, 0.0f, 1.0f) * alpha;
      f.color = c;
    };
    auto hole_shader = [&](const math::Vector &, Fragment &f) {
      f.color = Color4(0, 0, 0, alpha);
    };
    for (const Galaxy &g : galaxies) {
      const math::Vector core = rotation.apply(g.core);
      if (params.black_hole) {
        Scan::Point::draw<W, H>(PipelineRef(filters, canvas), canvas, core,
                                BULGE_RADIUS, hole_shader);
      } else {
        Scan::Point::draw<W, H>(PipelineRef(filters, canvas), canvas, core,
                                BULGE_RADIUS, bulge_shader);
      }
    }
  }

  math::Orientation<> orientation;
  FastNoiseLite noise;
  Filter::Screen::DirectAntiAliasSink<W, H> filters;
  ParticleSystem particle_system;
  std::array<Galaxy, NUM_GALAXIES> galaxies;
  /**
   * @brief Animation timeline for the orientation drift.
   * @details Declared after `orientation` so ~Timeline clears the RandomWalk
   *          that points at it first.
   */
  Timeline timeline;
};
