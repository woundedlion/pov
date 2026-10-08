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
 * @details Each vertex holds an attractor (the galactic core) and an emitter
 *          that spawns particles on a ring around it with an orbital
 *          velocity. A slow inward drift and the emitter's advancing spawn
 *          angle lay successive particles out as rotating arms. Particles are
 *          colored by distance from their own core.
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
   * @brief Registers params, bakes the palette, builds the particle system,
   *        and starts the orientation drift.
   */
  HS_COLD_MEMBER void init() override {
    static constexpr size_t SCRATCH_BYTES = 6 * 1024;
    // Scan::Point's row-interval buffers for the core bulges.
    static constexpr size_t SCRATCH_B_BYTES = 2 * 1024;
    ArenaSplit{SCRATCH_BYTES, SCRATCH_B_BYTES}.configure();

    // GLOBAL_ARENA_SIZE is inflated on host; check the device arena.
    static constexpr size_t POOL_BYTES =
        sizeof(Animation::Particle<TRAIL_LEN>) * NUM_PARTICLES;
    static constexpr size_t AUX_RESERVE_BYTES =
        6 * 1024 + BakedPalette::required_arena_bytes();
    static_assert(
        POOL_BYTES + AUX_RESERVE_BYTES <=
            ArenaSplit{SCRATCH_BYTES, SCRATCH_B_BYTES}.device_persistent(),
        "Galaxies particle pool + attractor/emitter/palette storage "
        "overflow the device persistent arena");

    // One trail at a time stages its fragments, pre-shader positions, and gate
    // arrays in scratch_a. A live tip is appended past the stored history.
    static constexpr size_t TRAIL_POINTS = TRAIL_LEN + 1;
    static_assert(TRAIL_POINTS * (sizeof(Fragment) + sizeof(math::Vector)) +
                          Plot::rasterize_scratch_a_bytes<W>(0, TRAIL_POINTS) <=
                      SCRATCH_BYTES,
                  "Galaxies trail staging exceeds its scratch_a split; retune "
                  "TRAIL_LEN or enlarge the split");

    register_param("Friction", &params.friction, 0.999f, 1.0f);
    register_param("Core Mass", &params.core_mass, 0.02f, 0.2f);
    register_param("Orbit Spd", &params.orbit_speed, 0.004f, 0.06f);
    register_param("Arm Spin", &params.arm_spin, 0.0f, 0.25f);
    register_int_param("Arms", &params.arms, 1, MAX_ARMS);
    register_param("Emission Rate", &params.emission_rate, 0.15f, 0.75f);
    register_int_param("Trail Length", &params.trail_length, 0, TRAIL_LEN - 1);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    // Cosine palette from a warm white core through lavender to blue arms.
    palette.bake(persistent_arena,
                 ProceduralPalette(/*bias*/ {0.62f, 0.6f, 0.9f},
                                   /*amp*/ {0.38f, 0.4f, 0.1f},
                                   /*freq*/ {0.5f, 0.5f, 0.5f},
                                   /*phase*/ {0.0f, 0.05f, 0.4f}));

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
        std::max((2.0f * math::PI_F) / W, params.orbit_speed * 1.4f);
    for (size_t i = 0; i < particle_system.attractors.size(); ++i)
      particle_system.attractors[i].strength = params.core_mass;
    {
      HS_PROFILE(gx_particle_step);
      particle_system.step(canvas);
    }
    trim_histories();

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

  /** @brief Maximum retained trail anchors per particle. */
  static constexpr int TRAIL_LEN = 6;
  /** @brief Frames between stored trail anchors. */
  static constexpr int TRAIL_SAMPLE_STRIDE = 2;
  /**
   * @brief Fixed particle pool capacity.
   * @details The emitters stop spawning while the pool is full.
   */
  static constexpr int NUM_PARTICLES = 3700;

  using ParticleSystem =
      Animation::ParticleSystem<W, NUM_PARTICLES, TRAIL_LEN, NUM_GALAXIES,
                                NUM_GALAXIES, true, TRAIL_SAMPLE_STRIDE>;

  static constexpr float GRAVITY = 0.001f;
  static constexpr float PARTICLE_LIFETIME_FRAMES = 800.0f;
  static constexpr float REFERENCE_ORBIT_SPEED = 0.01285f;
  /** @brief Angular radius of the spawn ring around each core (radians). */
  static constexpr float RING_RADIUS = 0.58f;
  /** @brief Core kill radius (chord) for particles that plunge inward. */
  static constexpr float KILL_RADIUS = 0.008f;
  /** @brief Galaxies do not use the attractor's radial steering zone. */
  static constexpr float EVENT_HORIZON = 0.0f;
  /** @brief Softening radius that limits the core force (chord). */
  static constexpr float SOFTENING_RADIUS = 0.02f;
  /** @brief Radius of the particle fade into each core (radians). */
  static constexpr float HOLE_FADE_RADIUS = 0.035f;
  /** @brief Spawn-angle jitter half-width (radians), thickens the arms. */
  static constexpr float ARM_JITTER = 0.14f;
  /** @brief Spawn-ring radius jitter half-width, as a fraction of the ring. */
  static constexpr float RING_JITTER = 0.015f;
  /** @brief Orbital speed jitter half-width, as a fraction of the speed. */
  static constexpr float SPEED_JITTER = 0.015f;
  /** @brief Opacity at a trail's tail; the head is fully opaque. */
  static constexpr float TRAIL_TAIL_ALPHA = 0.35f;
  static constexpr float MIN_PARTICLE_ALPHA = 0.3f;
  /** @brief Frames over which a new particle fades in. */
  static constexpr float FADE_IN_FRAMES = 16.0f;
  /** @brief Frames over which an expiring particle fades out. */
  static constexpr float FADE_OUT_FRAMES = 20.0f;
  /** @brief Radius where an inward-moving trail starts to dim (radians). */
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
   * @brief User-tunable parameters exposed via register_param.
   */
  struct Params {
    float friction = 0.99965f;   /**< Velocity retention per frame. */
    float core_mass = 0.12f;     /**< Attractor strength. */
    float orbit_speed = 0.0124f; /**< Reference spawn speed (radians/frame). */
    float arm_spin = 0.03f;      /**< Arm rotation (radians/frame). */
    int arms = 2;                /**< Arms per galaxy. */
    float emission_rate = 0.75f; /**< Particles per galaxy per frame. */
    int trail_length = 0;        /**< Visible trail anchors. */
    float alpha = 1.0f;          /**< Overall opacity. */
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
      g.core = basis.v;
      g.u = basis.u;
      g.w = basis.w;
      g.phase = hs::rand_f(0.0f, 2.0f * math::PI_F);
      g.spin = hs::rand_f() < 0.5f ? -1.0f : 1.0f;
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
        galaxy.emission_credit += hs::clamp(params.emission_rate, 0.15f, 0.75f);
        if (galaxy.emission_credit >= 1.0f) {
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
  void emit(Galaxy &g, int index) {
    const int arms = hs::clamp(params.arms, 1, MAX_ARMS);
    g.arm = static_cast<uint8_t>((g.arm + 1) % arms);
    // Skip the spawn rather than let spawn() log a dropped particle.
    if (particle_system.active() >= particle_system.pool.capacity())
      return;

    const float angle = g.phase + 2.0f * math::PI_F * g.arm / arms +
                        hs::rand_f(-ARM_JITTER, ARM_JITTER);
    const float ring =
        RING_RADIUS * (1.0f + hs::rand_f(-RING_JITTER, RING_JITTER));
    // fast_cosf/fast_sinf are approximate; renormalize onto the sphere.
    const math::Vector radial =
        (g.u * math::fast_cosf(angle) + g.w * math::fast_sinf(angle))
            .normalized();
    const math::Vector pos =
        (g.core * math::fast_cosf(ring) + radial * math::fast_sinf(ring))
            .normalized();
    // core and radial are orthonormal, so their cross is the unit tangent.
    const math::Vector tangent = math::cross(g.core, radial);
    const math::Vector outward = math::cross(tangent, pos);
    const float speed = circular_orbit_speed(pos, outward, ring) *
                        (params.orbit_speed / REFERENCE_ORBIT_SPEED) *
                        (1.0f + hs::rand_f(-SPEED_JITTER, SPEED_JITTER));
    const uint16_t alpha_seed = static_cast<uint16_t>(hs::rand_f() * 255.0f);
    const uint16_t color_seed =
        static_cast<uint16_t>(index | (alpha_seed << 8));
    particle_system.spawn(pos, tangent * (speed * g.spin), color_seed);
  }

  /** @brief Opacity stored in the high byte of a particle's color seed. */
  static float particle_alpha(uint16_t color_seed) {
    const float u = static_cast<float>(color_seed >> 8) * (1.0f / 255.0f);
    return MIN_PARTICLE_ALPHA + (1.0f - MIN_PARTICLE_ALPHA) * u * u;
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

  /** @brief Keeps only the number of trail anchors selected by the user. */
  void trim_histories() {
    const int length = hs::clamp(params.trail_length, 0, TRAIL_LEN - 1);
    const size_t limit = static_cast<size_t>(length == 0 ? 0 : length + 1);
    for (int i = 0; i < particle_system.active(); ++i) {
      auto &history = particle_system.pool[i].history;
      while (history.length() > limit)
        history.expire();
    }
  }

  /**
   * @brief Renders particles as points or short trails, then the core bulges.
   * @param canvas Target canvas.
   */
  void draw_particles(Canvas &canvas) {
    HS_PROFILE(gx_draw_particles);
    const math::RotationMatrix rotation(orientation.get());
    const math::Vector *core = nullptr;
    // Maps cos(distance to core) onto the palette: 0 at the core, 1 at the
    // spawn ring.
    const float cos_ring = math::fast_cosf(RING_RADIUS);
    const float ring_span = 1.0f - cos_ring;
    const float inv_ring_span = 1.0f / ring_span;
    const float cos_hole = math::fast_cosf(HOLE_FADE_RADIUS);
    const float inv_head_span =
        1.0f / (1.0f - math::fast_cosf(HEAD_FADE_RADIUS));
    const float max_life = static_cast<float>(particle_system.max_life);
    const float alpha = params.alpha;
    float core_fade = 1.0f;
    float depth_alpha = 1.0f;

    auto vertex_shader = [&](Fragment &f) { f.pos = rotation.apply(f.pos); };

    // v2 holds palette radius; size carries the core fade through rasterization.
    auto radius_shader = [&](FragmentRegisters f,
                             const math::Vector &original_pos) {
      const float cos_distance = math::dot(original_pos, *core);
      f.v2 = (1.0f - cos_distance) * inv_ring_span;
      f.size = cos_distance < cos_hole
                   ? 1.0f
                   : math::quintic_kernel(
                         math::fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f)) /
                         HOLE_FADE_RADIUS);
    };

    auto fragment_shader = [&](const math::Vector &, Fragment &f) {
      // v3 is remaining life over max_life.
      const float fade_in =
          hs::clamp((1.0f - f.v3) * max_life / FADE_IN_FRAMES, 0.0f, 1.0f);
      const float fade_out =
          hs::clamp(f.v3 * max_life / FADE_OUT_FRAMES, 0.0f, 1.0f);
      // The tail keeps TRAIL_TAIL_ALPHA, so the thin inner arms stay visible.
      const float trail = TRAIL_TAIL_ALPHA + (1.0f - TRAIL_TAIL_ALPHA) *
                                                 hs::clamp(f.v0, 0.0f, 1.0f);
      Color4 c = palette.get(hs::clamp(f.v2, 0.0f, 1.0f));
      c.alpha *=
          trail * fade_in * fade_out * f.size * core_fade * depth_alpha * alpha;
      f.color = c;
    };

    // The v2 mapper binds the particle's core before its fragments shade; the
    // deferred shader then overwrites v2.
    auto bind_core = [&](const auto &p, int) {
      core = &galaxies[p.color_seed & 0xff].core;
      core_fade = math::quintic_kernel((1.0f - math::dot(p.position, *core)) *
                                       inv_head_span);
      depth_alpha = particle_alpha(p.color_seed);
      return 0.0f;
    };

    filters.prepare(canvas);
    if (params.trail_length == 0) {
      for (int i = 0; i < particle_system.active(); ++i) {
        const auto &p = particle_system.pool[i];
        const math::Vector &point_core = galaxies[p.color_seed & 0xff].core;
        const float cos_distance = math::dot(p.position, point_core);
        const float radius = (1.0f - cos_distance) * inv_ring_span;
        const float hole =
            cos_distance < cos_hole
                ? 1.0f
                : math::quintic_kernel(
                      math::fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f)) /
                      HOLE_FADE_RADIUS);
        const float head_fade =
            math::quintic_kernel((1.0f - cos_distance) * inv_head_span);
        const float fade_in =
            hs::clamp((max_life - static_cast<float>(p.life)) / FADE_IN_FRAMES,
                      0.0f, 1.0f);
        const float fade_out =
            hs::clamp(static_cast<float>(p.life) / FADE_OUT_FRAMES, 0.0f, 1.0f);
        Color4 c = palette.get(hs::clamp(radius, 0.0f, 1.0f));
        c.alpha *= hole * head_fade * fade_in * fade_out *
                   particle_alpha(p.color_seed) * alpha;
        filters.plot(canvas, rotation.apply(p.position), c.color, 0.0f,
                     c.alpha);
      }
    } else {
      Plot::ParticleSystem::draw_fused_vertex<W, H, true>(
          filters, canvas, particle_system, fragment_shader, vertex_shader,
          radius_shader, bind_core);
    }

    // The bulge covers the final fade around each core.
    // Scan::Point leaves its quintic coverage in v2.
    const Color4 core_color = palette.get(0.0f);
    auto bulge_shader = [&](const math::Vector &, Fragment &f) {
      Color4 c = core_color;
      c.alpha *= hs::clamp(f.v2, 0.0f, 1.0f) * alpha;
      f.color = c;
    };
    for (const Galaxy &g : galaxies)
      Scan::Point::draw<W, H>(PipelineRef(filters, canvas), canvas,
                              rotation.apply(g.core), BULGE_RADIUS,
                              bulge_shader);
  }

  math::Orientation<> orientation;
  FastNoiseLite noise;
  Filter::Screen::DirectAntiAliasSink<W, H> filters;
  BakedPaletteStorage palette;
  ParticleSystem particle_system;
  std::array<Galaxy, NUM_GALAXIES> galaxies;
  /**
   * @brief Animation timeline for the orientation drift.
   * @details Declared after `orientation` so ~Timeline clears the RandomWalk
   *          that points at it first.
   */
  Timeline timeline;
};
