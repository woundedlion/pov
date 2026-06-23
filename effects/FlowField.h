/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Particle flow field over the sphere: particles advect along a 3D
 *        OpenSimplex2 noise field, leaving SLERP-blurred trails.
 * @tparam W Equirectangular canvas width in pixels.
 * @tparam H Equirectangular canvas height in pixels.
 * @details The noise field evolves over time via `t`. The whole field slowly
 *          drifts via a RandomWalk on `orientation`, applied by the
 *          Filter::World::Orient pipeline stage.
 */
template <int W, int H> class FlowField : public Effect {
public:
  /**
   * @brief Constructs the effect, building the palette and filter pipeline.
   * @details Builds the analogous-harmony palette and the Orient + AntiAlias
   *          filter pipeline; noise seeding and emitter setup happen in init().
   */
  FLASHMEM FlowField()
      : Effect(W, H), palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                              BrightnessProfile::ASCENDING),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Registers UI params, seeds noise, and installs the particle emitter.
   * @details Registers UI params, seeds the noise generators, starts the
   *          orientation drift, and installs the emitter that keeps the pool
   *          full and advances particles along the noise force field each step.
   */
  void init() override {
    registerParam("Scale", &params.noise_scale, 0.1f, 10.0f);
    registerParam("Force", &params.force_scale, 0.001f, 0.05f);
    registerParam("Max Spd", &params.max_speed, 0.01f, 0.1f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Time Spd", &params.time_scale, 0.001f, 0.05f);

    noise_generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    // rand_int is half-open [min, max), so use 65536 to make the full 16-bit
    // seed range [0, 65535] reachable rather than stopping one short at 65534.
    noise_generator.SetSeed(hs::rand_int(0, 65536));

    // Drive the whole-sphere drift consumed by the Filter::World::Orient stage.
    // Without an animation the orientation stays at identity and the filter
    // (plus its SLERP motion-blur sweep) is a no-op. `orient_noise` is a
    // dedicated generator so RandomWalk never reconfigures the flow-field noise.
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, random_vector(), orient_noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    // Particle dynamics: light per-frame damping, no gravity (the noise flow
    // field is the only force), and a long lifetime so streamlines read as
    // continuous. Named rather than passed as bare positionals.
    constexpr float kFriction = 0.96f;
    constexpr float kGravity = 0.0f;
    constexpr float kMaxLife = 300.0f;
    particle_system.init(persistent_arena, kFriction, kGravity, kMaxLife);

    // Per-step emitter: keep the pool saturated, then push each particle by the
    // noise force field clamped to max_speed.
    particle_system.add_emitter([this](ParticleSystem &sys) mutable {
      while (sys.active_count < sys.pool.capacity()) {
        sys.spawn(random_vector(), Vector(0, 0, 0), 0);
      }

      for (size_t i = 0; i < sys.active_count; ++i) {
        auto &p = sys.pool[i];

        // Rare random respawn so particles caught in field sinks redistribute.
        // ~0.5% per particle per frame: long enough average dwell that motion
        // reads as continuous, frequent enough to bleed off sink accumulation.
        constexpr float kRespawnProb = 0.005f;
        if (hs::rand_f() < kRespawnProb) {
          p.position = random_vector();
          p.velocity = Vector(0, 0, 0);
          p.history.clear();
        }

        // Sample three noise channels for the force vector's components. They
        // are mutually decorrelated by sampling at large x-offsets (0/100/200),
        // but all three share the same z input (p.z + t), so time evolution
        // advances the field along the z axis: the temporal evolution is
        // coupled to the z spatial structure, not a fully independent 4th axis.
        // FastNoiseLite is only 3D here, so a genuine 4th time dimension isn't
        // available; this shared-z animation is the accepted simplification.
        // The Y/Z channels sample the field shifted along x by these offsets so
        // the three components decorrelate (large enough to clear the noise
        // lattice's correlation length).
        constexpr float kChannelOffsetY = 100.0f;
        constexpr float kChannelOffsetZ = 200.0f;
        float fx =
            noise_generator.GetNoise(p.position.x * params.noise_scale,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        float fy =
            noise_generator.GetNoise(p.position.x * params.noise_scale +
                                         kChannelOffsetY,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        float fz =
            noise_generator.GetNoise(p.position.x * params.noise_scale +
                                         kChannelOffsetZ,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        Vector force(fx, fy, fz);

        p.velocity = p.velocity + force;

        float speed = p.velocity.magnitude();
        if (speed > params.max_speed) {
          p.velocity = (p.velocity / speed) * params.max_speed;
        }
      }
    });
  }

  /**
   * @brief Reports whether a background fill is drawn behind the effect.
   * @return False; trails fully render the frame, so no background fill is
   *         needed behind them.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances the simulation and plots the particle trails for one frame.
   * @details Advances noise time, steps the orientation drift and particles,
   *          then plots the particle trails with a latitude-mapped palette
   *          through the orient/AA filters.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    // Wrap the noise-time accumulator so it never grows large enough for the
    // float ULP to swallow the increment and freeze the field. OpenSimplex2 is
    // not periodic, so the wrap produces a brief flow-field pop once per period;
    // at 2048 the ULP (~2.4e-4) stays well under the slowest Time Spd increment
    // (0.001) so t always advances, and the pop is hours apart at the default
    // rate (and power/show cycles reset it far sooner).
    t = fmodf(t + params.time_scale, TIME_PERIOD);

    timeline.step(canvas);
    particle_system.step(canvas);

    // Rotation is owned by the Filter::World::Orient stage in `filters`, which
    // also sweeps the orientation history for SLERP motion blur. Applying it
    // again in a vertex shader would rotate every point twice, so none is used.
    // Color by latitude (v.y in [-1,1] -> palette t in [0,1]); the trail's
    // per-vertex alpha (f.v0) and the global Alpha param fade the tail.
    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      float alpha = f.v0;
      float palette_t = (v.y + 1.0f) / 2.0f;
      Color4 c = palette.get(palette_t);
      c.alpha = c.alpha * alpha * params.alpha;
      f.color = c;
    };

    Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                     fragment_shader);
  }

private:
  /**
   * @brief Particle pool capacity.
   * @details Sets the per-frame force-loop cost: each active particle samples the
   * noise field 3x per frame (one OpenSimplex2 call per force component), so the
   * loop is a fixed ~k_num_particles*3 (~1800) noise evaluations/frame. Unlike
   * Raymarch's per-pixel marching this cost is independent of W*H -- it does not
   * scale with the pixel/column budget -- so it carries no resolution cost cliff
   * and fits the device frame budget at any supported resolution. Lower this if a
   * future device target needs more margin.
   */
  static constexpr int k_num_particles = 600;
  static constexpr int k_trail_length = 14; /**< Per-particle trail length. */
  /**
   * @brief Wrap period for the noise-time accumulator (see draw_frame).
   * @details Large enough that flow-field pops are far apart, small enough that
   *          the float ULP never swallows the per-step time increment.
   */
  static constexpr float TIME_PERIOD = 2048.0f;

  /** @brief Particle system specialization for this effect's pool and trails. */
  typedef Animation::ParticleSystem<W, k_num_particles, k_trail_length, 1, 0>
      ParticleSystem;

  /** @brief User-tunable parameters exposed through the effect UI. */
  struct Params {
    float noise_scale = 2.0f; /**< Noise sampling frequency over the sphere. */
    float force_scale = 0.005f; /**< Scales noise samples into a force vector. */
    float max_speed = 0.03f; /**< Per-step velocity clamp magnitude. */
    float alpha = 0.8f; /**< Global trail opacity multiplier in [0,1]. */
    float time_scale = 0.005f; /**< Per-step noise-time increment. */
  } params; /**< User-tunable parameter values. */

  float t = 0; /**< Noise-time accumulator, wrapped to TIME_PERIOD. */
  FastNoiseLite noise_generator; /**< Generator for the flow-field force noise. */
  FastNoiseLite orient_noise; /**< Dedicated generator for the orientation drift. */
  GenerativePalette palette; /**< Latitude-mapped trail color palette. */
  ParticleSystem particle_system; /**< Pool of advected particles and trails. */
  Orientation<> orientation; /**< Whole-sphere drift consumed by the Orient stage. */
  Timeline timeline; /**< Drives the orientation RandomWalk animation. */

  /** @brief World-orient then screen anti-alias filter pipeline. */
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(FlowField)
