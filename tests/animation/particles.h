/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// ParticleSystem
// ============================================================================

/**
 * @brief Verifies ParticleSystem::spawn adds particles and drops
 * spawns once the fixed pool is at capacity.
 */
inline void test_particle_system_spawn_and_capacity_guard() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps; // CAPACITY = 4
  ps.init(arena);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);

  ps.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);
  ps.spawn(math::Vector(0, 1, 0), math::Vector(0, 0, 0), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 2);

  ps.spawn(math::Vector(0, 0, 1), math::Vector(0, 0, 0), 2);
  ps.spawn(math::Vector(-1, 0, 0), math::Vector(0, 0, 0), 3);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 4);
  HS_EXPECT_EQ(ps.dropped_spawns(), uint32_t{0});
  ps.spawn(math::Vector(0, -1, 0), math::Vector(0, 0, 0),
           4); // capacity is 4 — rejected
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 4);
  HS_EXPECT_EQ(ps.dropped_spawns(), uint32_t{1});
}

/**
 * @brief Verifies both valid ParticleSystem lifetime boundaries are preserved.
 */
inline void test_particle_system_lifetime_boundaries() {
  static uint8_t buf[4096];
  {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 1> ps;
    ps.init(arena, 0.85f, 0.0f, 1.0f);
    HS_EXPECT_EQ(static_cast<int>(ps.max_life), 1);
  }
  {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 1> ps;
    ps.init(arena, 0.85f, 0.0f, 65535.0f);
    HS_EXPECT_EQ(static_cast<int>(ps.max_life), 65535);
  }
}

/**
 * @brief Verifies a particle is reclaimed as soon as its life reaches 0, even
 *        with recorded trail points, so the invisible history holds no slot.
 */
inline void test_particle_system_reclaims_at_life_expiry() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  // Gravity 0 + zero velocity isolates life expiry from physics kills.
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/3.0f);
  ps.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);
  HS_EXPECT_GT(ps.pool[0].history.length(), (size_t)0);
  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);

  ps.spawn(math::Vector(0, 1, 0), math::Vector(0, 0, 0), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);
  HS_EXPECT_EQ(ps.dropped_spawns(), (uint32_t)0);
}

/**
 * @brief Verifies an attractor removes a particle inside its kill_radius via
 * the kill check rather than by life expiry.
 */
inline void test_particle_system_attractor_kills_within_radius() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.001f, /*max_life=*/600.0f);
  ps.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);
  // Attractor co-located (distance 0 < kill_radius): removed by the kill check,
  // not life expiry (life is 600).
  ps.add_attractor(math::Vector(1, 0, 0), /*strength=*/1.0f,
                   /*kill_radius=*/0.5f,
                   /*event_horizon=*/2.0f);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);
}

/**
 * @brief Verifies spawn() initializes a particle's fields and that one step
 * advances it (life decrements, trail records) without dropping it.
 * @details Checks the spawned position, velocity, seed, max_life and empty
 *          trail, then one step: alive, life decremented and one trail point.
 */
inline void test_particle_system_spawn_initializes_and_steps() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  // Zero gravity and no attractor isolate the particle's inertial step.
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/120.0f);

  const math::Vector pos(0.6f, 0.0f, 0.8f);
  const math::Vector vel(0.0f, 0.01f, 0.0f);
  ps.spawn(pos, vel, /*seed=*/42);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  const auto &p = ps.pool[0];
  HS_EXPECT_NEAR(p.position.x, pos.x, 0.0f);
  HS_EXPECT_NEAR(p.position.y, pos.y, 0.0f);
  HS_EXPECT_NEAR(p.position.z, pos.z, 0.0f);
  HS_EXPECT_NEAR(p.velocity.x, vel.x, 0.0f);
  HS_EXPECT_NEAR(p.velocity.y, vel.y, 0.0f);
  HS_EXPECT_NEAR(p.velocity.z, vel.z, 0.0f);
  HS_EXPECT_EQ(static_cast<int>(p.color_seed), 42);
  HS_EXPECT_EQ(static_cast<int>(p.life), 120); // == system max_life
  HS_EXPECT_EQ(static_cast<int>(p.history_length()), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.pool[0].life), 119); // life-- per step
  HS_EXPECT_GT(static_cast<int>(ps.pool[0].history_length()), 0);
}

/** @brief Verifies sparse trails retain anchors at the configured cadence. */
inline void test_particle_system_sparse_trail_sampling() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1, 8, 8, 8, false, 3> ps;
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/30.0f);
  ps.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.pool[0].history_length(), (size_t)1);
  for (int i = 0; i < 2; ++i)
    ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.pool[0].history_length(), (size_t)1);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.pool[0].history_length(), (size_t)2);
  for (int i = 0; i < 18; ++i)
    ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.pool[0].history_length(), (size_t)8);
}

/**
 * @brief Pins the strict attractor kill-radius boundary.
 * @details Probes the start-of-step distance at radius-0.01, radius and
 *          radius+0.01: only the inside particle is killed.
 */
inline void test_particle_system_attractor_kill_radius_boundary() {
  static uint8_t buf[256 * 1024];
  constexpr float kr = 0.5f;

  auto survives = [&](float dist) {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 4> ps;
    ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.001f, /*max_life=*/600.0f);
    ps.spawn(math::Vector(1, 0, 0), math::Vector(0, 0, 0), 0);
    ps.add_attractor(math::Vector(1.0f + dist, 0, 0), /*strength=*/1.0f,
                     /*kill_radius=*/kr, /*event_horizon=*/0.0f);
    ps.step(fake_canvas());
    return ps.active() == 1;
  };

  HS_EXPECT_FALSE(survives(kr - 0.01f)); // inside the radius → killed
  HS_EXPECT_TRUE(survives(kr)); // exactly at the radius → survives (strict <)
  HS_EXPECT_TRUE(survives(kr + 0.01f));
}

template <int CAPACITY, bool SIGNED_AXIS>
using AxisParticleSystem =
    Animation::ParticleSystem<288, CAPACITY, 23, 1, 6, SIGNED_AXIS>;

template <typename PS>
inline void add_signed_axis_attractors(PS &ps, float strength = 0.85f,
                                       float kill_radius = 0.003f,
                                       float event_horizon = 0.2f) {
  const math::Vector axes[] = {math::X_AXIS,  -math::X_AXIS, math::Y_AXIS,
                               -math::Y_AXIS, math::Z_AXIS,  -math::Z_AXIS};
  for (const math::Vector &axis : axes)
    ps.add_attractor(axis, strength, kill_radius, event_horizon);
}

/** @brief Bounds one-step signed-axis physics against the generic path. */
inline void test_particle_system_signed_axis_one_step_equivalence() {
  constexpr int COUNT = 256;
  static uint8_t reference_buf[128 * 1024];
  static uint8_t specialized_buf[128 * 1024];
  Arena reference_arena(reference_buf, sizeof(reference_buf));
  Arena specialized_arena(specialized_buf, sizeof(specialized_buf));
  AxisParticleSystem<COUNT, false> reference;
  AxisParticleSystem<COUNT, true> specialized;
  reference.init(reference_arena, 0.85f, 0.001f, 160.0f);
  specialized.init(specialized_arena, 0.85f, 0.001f, 160.0f);
  add_signed_axis_attractors(reference, 2.55f);
  add_signed_axis_attractors(specialized, 2.55f);
  constexpr float STRENGTHS[] = {0.35f, 1.1f, 2.55f, 0.8f, 3.2f, 1.7f};
  for (size_t i = 0; i < 6; ++i) {
    reference.attractors[i].strength = STRENGTHS[i];
    specialized.attractors[i].strength = STRENGTHS[i];
  }

  const auto saved = hs::random();
  hs::random().seed(0x61786973);
  for (int i = 0; i < COUNT; ++i) {
    const math::Vector pos = rand_unit();
    const float vel_x = hs::rand_f(-0.1f, 0.1f);
    const float vel_y = hs::rand_f(-0.1f, 0.1f);
    const float vel_z = hs::rand_f(-0.1f, 0.1f);
    math::Vector velocity(vel_x, vel_y, vel_z);
    reference.spawn(pos, velocity, static_cast<uint16_t>(i));
    specialized.spawn(pos, velocity, static_cast<uint16_t>(i));
  }
  hs::random() = saved;

  reference.step(fake_canvas());
  specialized.step(fake_canvas());
  HS_EXPECT_EQ(reference.active(), specialized.active());
  float max_position_error = 0.0f;
  float max_velocity_error = 0.0f;
  double max_angle_error = 0.0;
  float max_norm_drift = 0.0f;
  int color_seed_mismatches = 0;
  int life_mismatches = 0;
  for (size_t i = 0; i < reference.active(); ++i) {
    const auto &a = reference.pool[i];
    const auto &b = specialized.pool[i];
    color_seed_mismatches += a.color_seed != b.color_seed;
    life_mismatches += a.life != b.life;
    max_position_error = hs_test::fold_worst(
        max_position_error, max_component_delta(a.position, b.position));
    max_velocity_error = hs_test::fold_worst(
        max_velocity_error, max_component_delta(a.velocity, b.velocity));
    max_angle_error = hs_test::fold_worst(
        max_angle_error, small_angle_between(a.position, b.position));
    max_norm_drift = hs_test::fold_worst(
        max_norm_drift, std::abs(math::dot(a.position, a.position) - 1.0f));
    max_norm_drift = hs_test::fold_worst(
        max_norm_drift, std::abs(math::dot(b.position, b.position) - 1.0f));
  }
  std::printf("axis one-step particles=%u pos=%.9g vel=%.9g angle=%.9g "
              "norm=%.9g\n",
              reference.active(), max_position_error, max_velocity_error,
              max_angle_error, max_norm_drift);
  // Components within COMPONENT_BOUND give a chord of at most
  // sqrt(3)*COMPONENT_BOUND; normalizing each endpoint at most doubles it.
  constexpr float COMPONENT_BOUND = 2e-7f;
  const double angle_bound = 2.0 * std::sqrt(3.0) * COMPONENT_BOUND;
  HS_EXPECT_EQ(color_seed_mismatches, 0);
  HS_EXPECT_EQ(life_mismatches, 0);
  HS_EXPECT_LE(max_position_error, COMPONENT_BOUND);
  HS_EXPECT_LE(max_velocity_error, COMPONENT_BOUND);
  HS_EXPECT_LE(max_angle_error, angle_bound);
  HS_EXPECT_LE(max_norm_drift, 1e-6f);
}

/** @brief Pins signed-axis kill, horizon, and cross-axis fallback boundaries. */
inline void test_particle_system_signed_axis_boundaries() {
  auto compare = [](const math::Vector &position, const math::Vector &velocity,
                    float kill_radius, float event_horizon,
                    int expected_active) {
    uint8_t reference_buf[4096];
    uint8_t specialized_buf[4096];
    Arena reference_arena(reference_buf, sizeof(reference_buf));
    Arena specialized_arena(specialized_buf, sizeof(specialized_buf));
    AxisParticleSystem<1, false> reference;
    AxisParticleSystem<1, true> specialized;
    reference.init(reference_arena, 0.85f, 0.001f, 160.0f);
    specialized.init(specialized_arena, 0.85f, 0.001f, 160.0f);
    add_signed_axis_attractors(reference, 1.0f, kill_radius, event_horizon);
    add_signed_axis_attractors(specialized, 1.0f, kill_radius, event_horizon);
    reference.spawn(position, velocity, 1);
    specialized.spawn(position, velocity, 1);
    reference.step(fake_canvas());
    specialized.step(fake_canvas());
    HS_EXPECT_EQ(reference.active(), expected_active);
    HS_EXPECT_EQ(reference.active(), specialized.active());
    if (reference.active()) {
      HS_EXPECT_EQ(reference.pool[0].life, specialized.pool[0].life);
      HS_EXPECT_LE(max_component_delta(reference.pool[0].position,
                                       specialized.pool[0].position),
                   1e-5f);
      HS_EXPECT_LE(max_component_delta(reference.pool[0].velocity,
                                       specialized.pool[0].velocity),
                   1e-5f);
    }
  };
  auto at_chord_distance = [](float distance) {
    const float x = 1.0f - 0.5f * distance * distance;
    const float DX = x - 1.0f;
    return math::Vector(x, sqrtf(distance * distance - DX * DX), 0.0f);
  };
  auto boundary_positions = [](float radius) {
    return std::array<math::Vector, 3>{
        math::Vector(1.0f, std::nextafter(radius, 0.0f), 0.0f),
        math::Vector(1.0f, radius, 0.0f),
        math::Vector(
            1.0f,
            std::nextafter(radius, std::numeric_limits<float>::infinity()),
            0.0f)};
  };

  constexpr float KILL = 0.003f;
  constexpr float HORIZON = 0.2f;
  for (const math::Vector &axis :
       {math::X_AXIS, -math::X_AXIS, math::Y_AXIS, -math::Y_AXIS, math::Z_AXIS,
        -math::Z_AXIS}) {
    math::Vector tangent = math::cross(axis, math::X_AXIS);
    if (math::dot(tangent, tangent) < 0.5f)
      tangent = math::cross(axis, math::Y_AXIS);
    for (float radius : {KILL, HORIZON}) {
      const auto boundary = boundary_positions(radius);
      for (int side = 0; side < 3; ++side) {
        const auto &local = boundary[side];
        compare(axis * local.x + tangent * local.y, math::Vector(), KILL,
                HORIZON, radius == KILL && side == 0 ? 0 : 1);
      }
      for (float factor : {0.5f, 2.0f}) {
        const auto local = at_chord_distance(radius * factor);
        compare(axis * local.x + tangent * local.y, math::Vector(), KILL,
                HORIZON, radius == KILL && factor < 1.0f ? 0 : 1);
      }
    }
  }

  for (const math::Vector &position : {
           math::X_AXIS,
           math::Vector(1.0f, 1e-7f, 0.0f).normalized(),
           math::Vector(1.0f, 1e-5f, 0.0f).normalized(),
           -math::X_AXIS,
           math::Vector(-1.0f, 1e-7f, 0.0f).normalized(),
           math::Vector(-1.0f, 1e-5f, 0.0f).normalized(),
       })
    compare(position, math::Vector(), 0.0f, 0.0f, 1);
}

/** @brief Bounds deterministic multi-step signed-axis trajectory divergence. */
inline void test_particle_system_signed_axis_trajectory() {
  constexpr int COUNT = 32;
  constexpr int STEPS = 120;
  static uint8_t reference_buf[32 * 1024];
  static uint8_t specialized_buf[32 * 1024];
  Arena reference_arena(reference_buf, sizeof(reference_buf));
  Arena specialized_arena(specialized_buf, sizeof(specialized_buf));
  AxisParticleSystem<COUNT, false> reference;
  AxisParticleSystem<COUNT, true> specialized;
  reference.init(reference_arena, 0.85f, 0.001f, 160.0f);
  specialized.init(specialized_arena, 0.85f, 0.001f, 160.0f);
  add_signed_axis_attractors(reference, 2.55f);
  add_signed_axis_attractors(specialized, 2.55f);

  const math::Vector cube[] = {
      math::Vector(-1, -1, -1).normalized(),
      math::Vector(1, -1, -1).normalized(),
      math::Vector(1, 1, -1).normalized(),
      math::Vector(-1, 1, -1).normalized(),
      math::Vector(-1, -1, 1).normalized(),
      math::Vector(1, -1, 1).normalized(),
      math::Vector(1, 1, 1).normalized(),
      math::Vector(-1, 1, 1).normalized(),
  };
  for (int i = 0; i < COUNT; ++i) {
    const math::Vector pos = cube[i % 8];
    math::Vector velocity =
        math::cross(pos, i % 2 ? math::Y_AXIS : math::Z_AXIS).normalized();
    velocity *= 0.025f + 0.069f * static_cast<float>(i % 7) / 6.0f;
    reference.spawn(pos, velocity, static_cast<uint16_t>(i));
    specialized.spawn(pos, velocity, static_cast<uint16_t>(i));
  }

  float max_position_error = 0.0f;
  float max_velocity_error = 0.0f;
  double max_angle_error = 0.0;
  float max_norm_drift = 0.0f;
  int active_mismatches = 0;
  int color_seed_mismatches = 0;
  for (int step = 0; step < STEPS; ++step) {
    reference.step(fake_canvas());
    specialized.step(fake_canvas());
    active_mismatches += reference.active() != specialized.active();
    for (size_t i = 0; i < reference.active(); ++i) {
      const auto &a = reference.pool[i];
      const auto &b = specialized.pool[i];
      color_seed_mismatches += a.color_seed != b.color_seed;
      max_position_error = hs_test::fold_worst(
          max_position_error, max_component_delta(a.position, b.position));
      max_velocity_error = hs_test::fold_worst(
          max_velocity_error, max_component_delta(a.velocity, b.velocity));
      max_angle_error = hs_test::fold_worst(
          max_angle_error, small_angle_between(a.position, b.position));
      max_norm_drift = hs_test::fold_worst(
          max_norm_drift, std::abs(math::dot(a.position, a.position) - 1.0f));
      max_norm_drift = hs_test::fold_worst(
          max_norm_drift, std::abs(math::dot(b.position, b.position) - 1.0f));
    }
  }
  std::printf("axis trajectory steps=%d particles=%u pos=%.9g vel=%.9g "
              "angle=%.9g norm=%.9g\n",
              STEPS, reference.active(), max_position_error, max_velocity_error,
              max_angle_error, max_norm_drift);
  HS_EXPECT_EQ(active_mismatches, 0);
  HS_EXPECT_EQ(color_seed_mismatches, 0);
  constexpr float COMPONENT_BOUND = 1e-5f;
  const double ANGLE_BOUND = 2.0 * std::sqrt(3.0) * COMPONENT_BOUND;
  HS_EXPECT_LE(max_position_error, COMPONENT_BOUND);
  HS_EXPECT_LE(max_velocity_error, COMPONENT_BOUND);
  HS_EXPECT_LE(max_angle_error, ANGLE_BOUND);
  HS_EXPECT_LE(max_norm_drift, 8e-6f);
}

/** @brief Point storage preserves directions, tangent motion, and lifetime. */
inline void test_point_particle_storage_and_lifetime() {
  HS_EXPECT_EQ(sizeof(Animation::PointParticle), size_t{24});
  for (int x = -1; x <= 1; ++x) {
    for (int y = -1; y <= 1; ++y) {
      for (int z = -1; z <= 1; ++z) {
        if (x == 0 && y == 0 && z == 0)
          continue;
        const math::Vector direction = math::Vector(x, y, z).normalized();
        Animation::PointParticle p;
        p.init(direction, {0.04f, 0.02f, -0.03f}, 1234, 800);
        const math::Vector decoded = p.get_position();
        HS_EXPECT_LE((decoded - direction).magnitude(), 5e-7f);
        HS_EXPECT_NEAR(decoded.magnitude(), 1.0f, 2e-7f);
        HS_EXPECT_NEAR(math::dot(decoded, p.velocity), 0.0f, 1e-8f);
        HS_EXPECT(std::isfinite(p.velocity.magnitude()), "velocity is finite");
        HS_EXPECT_EQ(p.history_length(), size_t{0});
        const auto copy = p;
        HS_EXPECT_EQ(copy.color_seed, 1234);
        HS_EXPECT_EQ(copy.life, 800);
        HS_EXPECT_NEAR((copy.get_position() - decoded).magnitude(), 0.0f, 0.0f);
      }
    }
  }
  for (const math::Vector source :
       {math::Vector(1.0f, 0.3f, 1e-7f), math::Vector(1.0f, 0.3f, -1e-7f),
        math::Vector(-0.3f, -1.0f, 1e-7f), math::Vector(-0.3f, -1.0f, -1e-7f),
        math::Vector(1e-7f, -1e-7f, -1.0f),
        math::Vector(-1e-7f, 1e-7f, -1.0f)}) {
    const math::Vector direction = source.normalized();
    Animation::PointParticle p;
    p.init(direction, {0.04f, 0.02f, -0.03f}, 0, 800);
    const math::Vector decoded = p.get_position();
    HS_EXPECT_LE((decoded - direction).magnitude(), 5e-7f);
    HS_EXPECT_NEAR(decoded.magnitude(), 1.0f, 2e-7f);
    HS_EXPECT_NEAR(math::dot(decoded, p.velocity), 0.0f, 1e-8f);
  }
  Animation::PointParticle p;
  p.init(math::X_AXIS, {}, 0, std::numeric_limits<float>::infinity());
  HS_EXPECT_EQ(p.life, 0);
  p.init(math::Y_AXIS, {}, 0, 70000);
  HS_EXPECT_EQ(p.life, 65535);

  static uint8_t storage[4096];
  Arena arena(storage, sizeof(storage));
  Animation::ParticleSystem<288, 2, 0, 1, 1, false, 1, Animation::PointParticle>
      ps;
  ps.init(arena, 1.0f, 0.0f, 3.0f);
  ps.spawn(math::X_AXIS, {0.0f, 0.01f, 0.0f}, 17);
  ps.spawn(math::Y_AXIS, {0.0f, 0.0f, 0.01f}, 29);
  ps.pool[0].life = 1;
  ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.active(), 1);
  HS_EXPECT_EQ(ps.pool[0].color_seed, 29);
  HS_EXPECT_EQ(ps.pool[0].life, 2);
  HS_EXPECT_GT(ps.pool[0].get_position().z, 0.005f);
  HS_EXPECT_NEAR(ps.pool[0].get_position().magnitude(), 1.0f, 2e-7f);
  HS_EXPECT_NEAR(math::dot(ps.pool[0].get_position(), ps.pool[0].velocity),
                 0.0f, 1e-8f);
  ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.active(), 1);
  ps.step(fake_canvas());
  HS_EXPECT_EQ(ps.active(), 0);
}
