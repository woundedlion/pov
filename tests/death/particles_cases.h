/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Particles death fixtures and guard cases.

/** @brief Death case: a zero particle lifetime must trap. */
inline void case_particle_lifetime_zero() {
  init_particle_system_with_lifetime(opaque(0.0f));
}

inline void case_particle_friction_nan() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, opaque(std::numeric_limits<float>::quiet_NaN()));
}

inline void case_particle_gravity_nan() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, opaque(std::numeric_limits<float>::quiet_NaN()));
}

/** @brief Death case: a NaN particle lifetime must trap. */
inline void case_particle_lifetime_nan() {
  init_particle_system_with_lifetime(
      opaque(std::numeric_limits<float>::quiet_NaN()));
}

/** @brief Death case: a particle lifetime above uint16_t must trap. */
inline void case_particle_lifetime_over_max() {
  init_particle_system_with_lifetime(opaque(65536.0f));
}

/** @brief Death case: a nonempty render with zero max_life must trap. */
inline void case_particle_render_zero_lifetime() {
  configure_arenas_default();
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, 1.0f);
  ps.spawn(math::Vector(1, 0, 0), math::Vector(), 0);
  ps.max_life = 0;

  DeathEffect fx;
  Canvas canvas(fx);
  DeathPlotPipeline pipeline;
  Plot::ParticleSystem::draw<32, 16>(pipeline, canvas, ps,
                                     [](const math::Vector &, Fragment &) {});
}
