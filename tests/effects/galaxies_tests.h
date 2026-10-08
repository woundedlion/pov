/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Galaxies: octahedral cores, spawn geometry and galaxy containment.
// ---------------------------------------------------------------------------

/** @brief White-box accessor for Galaxies' spawn state and particle pool. */
struct GalaxiesWhiteBox {
  template <int W, int H> using FX = Galaxies<W, H>;

  static constexpr int NUM_GALAXIES = FX<1, 1>::NUM_GALAXIES;
  static constexpr float RING_RADIUS = FX<1, 1>::RING_RADIUS;
  static constexpr float SPEED_JITTER = FX<1, 1>::SPEED_JITTER;

  template <int W, int H> static void emit(Galaxies<W, H> &fx, int galaxy) {
    fx.emit(fx.galaxies[galaxy], galaxy);
  }
  template <int W, int H>
  static const math::Vector &core(const Galaxies<W, H> &fx, int galaxy) {
    return fx.galaxies[galaxy].core;
  }
  template <int W, int H> static int arm(const Galaxies<W, H> &fx, int galaxy) {
    return fx.galaxies[galaxy].arm;
  }
  template <int W, int H> static auto &system(Galaxies<W, H> &fx) {
    return fx.particle_system;
  }
  template <int W, int H>
  static float orbit_speed(const Galaxies<W, H> &fx, const math::Vector &pos,
                           const math::Vector &outward, float ring) {
    return fx.circular_orbit_speed(pos, outward, ring) *
           (fx.params.orbit_speed / fx.REFERENCE_ORBIT_SPEED);
  }
  template <int W, int H> static void set_arms(Galaxies<W, H> &fx, int arms) {
    fx.params.arms = arms;
  }
  template <int W, int H>
  static void set_emission_rate(Galaxies<W, H> &fx, float rate) {
    fx.params.emission_rate = rate;
  }
  template <int W, int H> static void hide(Galaxies<W, H> &fx) {
    fx.params.alpha = 0.0f;
  }
  template <int W, int H> static auto &galaxy(Galaxies<W, H> &fx, int index) {
    return fx.galaxies[index];
  }
  template <int W, int H>
  static void set_pitch(Galaxies<W, H> &fx, float pitch) {
    fx.params.arm_pitch = pitch;
  }
  static float alpha(uint16_t seed, float age, float arm) {
    return FX<1, 1>::particle_alpha(seed, age, arm);
  }
  template <int W, int H>
  static float density(const Galaxies<W, H> &fx, int index,
                       const math::Vector &position) {
    const auto &galaxy = fx.galaxies[index];
    return fx.arm_density(position, galaxy, math::dot(position, galaxy.core),
                          1.0f / tanf(fx.params.arm_pitch));
  }
};

/** @brief Stars span logarithmic arms with tangential orbital velocities. */
inline void test_galaxies_spawn_along_spiral_with_orbital_velocity() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<SMALL_W, SMALL_H> fx;
  fx.init();
  WB::set_arms(fx, 3);
  WB::set_pitch(fx, 0.5f);

  const auto &ps = WB::system(fx);
  HS_EXPECT_EQ(WB::NUM_GALAXIES, Solids::Octahedron::NUM_VERTS);
  int radial_bins[4] = {};
  for (int sample = 0; sample < 600; ++sample) {
    const int g = sample % WB::NUM_GALAXIES;
    const int before = WB::arm(fx, g);
    const uint16_t index = ps.active();
    WB::emit(fx, g);
    HS_EXPECT_EQ(static_cast<int>(ps.active()), index + 1);
    HS_EXPECT_EQ(WB::arm(fx, g), (before + 1) % 3);

    const auto &p = ps.pool[index];
    const math::Vector &core = WB::core(fx, g);
    HS_EXPECT_NEAR(math::dot(core, Solids::Octahedron::vertices[g]), 1.0f,
                   1e-5f);
    HS_EXPECT_EQ(static_cast<int>(p.color_seed & 0xff), g);
    HS_EXPECT_NEAR(p.position.magnitude(), 1.0f, 1e-4f);

    const float ring =
        acosf(hs::clamp(math::dot(p.position, core), -1.0f, 1.0f));
    HS_EXPECT_GE(ring, 0.139f);
    HS_EXPECT_LE(ring, 0.581f);
    ++radial_bins[hs::clamp(static_cast<int>((ring - 0.14f) / 0.11f), 0, 3)];
    const auto &galaxy = WB::galaxy(fx, g);
    const float azimuth = atan2f(math::dot(p.position, galaxy.w),
                                 math::dot(p.position, galaxy.u));
    const float expected = galaxy.phase + 2.0f * math::PI_F * galaxy.arm / 3 -
                           galaxy.spin * logf(ring / 0.58f) / tanf(0.5f);
    const float error =
        atan2f(sinf(azimuth - expected), cosf(azimuth - expected));
    const float spread = 1.0f - 0.9f * (p.color_seed >> 8) / 255.0f;
    HS_EXPECT_LE(fabsf(error), 0.22f * spread + 0.005f);

    // Tangent to the sphere and perpendicular to the core direction: a pure
    // orbit about the core, with no radial component.
    HS_EXPECT_NEAR(math::dot(p.velocity, p.position), 0.0f, 1e-4f);
    HS_EXPECT_NEAR(math::dot(p.velocity, core), 0.0f, 1e-4f);
    const math::Vector outward =
        (p.position * math::dot(p.position, core) - core).normalized();
    const float target_speed = WB::orbit_speed(fx, p.position, outward, ring);
    const float speed = p.velocity.magnitude();
    HS_EXPECT_GE(speed, target_speed * (1.0f - WB::SPEED_JITTER) - 1e-5f);
    HS_EXPECT_LE(speed, target_speed * (1.0f + WB::SPEED_JITTER) + 1e-5f);
  }
  for (const int count : radial_bins)
    HS_EXPECT_GE(count, 100);
}

/** @brief A rotating spiral lights young stars; old stars stay dim. */
inline void test_galaxies_arm_contrast_and_age_fade() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<DEFAULT_W, DEFAULT_H> fx;
  fx.init();
  WB::set_pitch(fx, 0.5f);
  for (int arms : {1, 2, 4, 8}) {
    WB::set_arms(fx, arms);
    for (int g = 0; g < WB::NUM_GALAXIES; ++g) {
      auto &galaxy = WB::galaxy(fx, g);
      for (float ring : {0.18f, 0.35f, 0.54f}) {
        const float angle =
            galaxy.phase - galaxy.spin * logf(ring / 0.58f) / tanf(0.5f);
        const math::Vector position =
            galaxy.core * cosf(ring) +
            (galaxy.u * cosf(angle) + galaxy.w * sinf(angle)) * sinf(ring);
        const float on_arm = WB::density(fx, g, position);
        HS_EXPECT_GE(on_arm, 0.99f);
        HS_EXPECT_GE(WB::alpha(0xff00, 100, on_arm), 0.92f);
        HS_EXPECT_NEAR(WB::alpha(0xff00, 475, on_arm),
                       (WB::alpha(0xff00, 100, on_arm) + 0.014f) * 0.5f, 1e-5f);
        HS_EXPECT_NEAR(WB::alpha(0xff00, 700, on_arm), 0.014f, 1e-5f);
        galaxy.phase += math::PI_F / arms;
        const float between_arms = WB::density(fx, g, position);
        HS_EXPECT_NEAR(between_arms, 0.0f, 1e-5f);
        HS_EXPECT_NEAR(WB::alpha(0xff00, 100, between_arms), 0.014f, 1e-5f);
        HS_EXPECT_NEAR(WB::alpha(0, 100, between_arms), 0.006f, 1e-5f);
        galaxy.phase -= math::PI_F / arms;
      }
    }
  }
  HS_EXPECT_LE(WB::alpha(0xff00, 100, 0.7f) / WB::alpha(0xff00, 100, 1.0f),
               WB::alpha(0, 100, 0.7f) / WB::alpha(0, 100, 1.0f));
}

/** @brief Emission consumes whole credits and discards births while full. */
inline void test_galaxies_emission_rate_and_capacity() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<SMALL_W, SMALL_H> fx;
  fx.init();
  bool found = false;
  for (const auto &param : fx.getParameters()) {
    if (std::string_view(param.name) == "Emission Rate") {
      found = true;
      HS_EXPECT_EQ(param.min, 0.15f);
      HS_EXPECT_EQ(param.max, 4.0f);
      HS_EXPECT_EQ(param.get(), 0.75f);
    }
  }
  HS_EXPECT(found, "Emission Rate is registered");

  auto &ps = WB::system(fx);
  auto &galaxy = WB::galaxy(fx, 0);
  galaxy.emission_credit = 0.0f;
  galaxy.phase = 0.0f;
  galaxy.spin = 1.0f;
  WB::set_emission_rate(fx, 2.5f);
  ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), 2);
  HS_EXPECT_EQ(galaxy.emission_credit, 0.5f);
  HS_EXPECT_NEAR(galaxy.phase, 0.028f, 1e-6f);
  ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), 5);
  HS_EXPECT_EQ(galaxy.emission_credit, 0.0f);
  HS_EXPECT_NEAR(galaxy.phase, 0.056f, 1e-6f);

  WB::set_emission_rate(fx, 4.0f);
  ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), 9);
  WB::set_emission_rate(fx, 10.0f);
  ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), 13);
  HS_EXPECT_EQ(galaxy.emission_credit, 0.0f);

  const int capacity = static_cast<int>(ps.pool.capacity());
  HS_EXPECT_EQ(capacity, 6000);
  for (int i = ps.active(); i < capacity; ++i)
    WB::emit(fx, 0);
  HS_EXPECT_EQ(ps.active(), capacity);
  galaxy.emission_credit = 0.5f;
  WB::set_emission_rate(fx, 4.0f);
  for (int i = 0; i < 8; ++i)
    ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), capacity);
  HS_EXPECT_EQ(galaxy.emission_credit, 0.5f);

  for (int i = 0; i < capacity; ++i)
    ps.pool[i].life = 0;
  Canvas canvas(fx);
  ps.step(canvas);
  HS_EXPECT_EQ(ps.active(), 0);
  WB::set_emission_rate(fx, 0.0f);
  ps.emitters[0](ps);
  HS_EXPECT_EQ(ps.active(), 0);
  HS_EXPECT_NEAR(galaxy.emission_credit, 0.65f, 1e-6f);
}

/**
 * @brief Runs the default settings and checks that particles stay near their
 *        own core.
 * @details Checks both 96 and 288 columns because their per-frame motion caps
 *          differ.
 */
template <int W, int H> inline void check_galaxies_stay_contained() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<W, H> fx;
  fx.init();
  // Physics only; the containment check needs no pixels.
  WB::hide(fx);
  for (int f = 0; f < 300; ++f) {
    pin_frame_clock(f);
    fx.draw_frame();
    fx.advance_display();
  }

  const auto &ps = WB::system(fx);
  const int live = ps.active();
  HS_EXPECT_GE(live, 200);
  int strays = 0;
  int outside = 0;
  for (int i = 0; i < live; ++i) {
    const auto &p = ps.pool[i];
    const float own = math::dot(p.position, WB::core(fx, p.color_seed & 0xff));
    for (int g = 0; g < WB::NUM_GALAXIES; ++g) {
      if (math::dot(p.position, WB::core(fx, g)) > own) {
        ++strays;
        break;
      }
    }
    if (acosf(hs::clamp(own, -1.0f, 1.0f)) > 1.5f * WB::RING_RADIUS)
      ++outside;
  }
  HS_EXPECT_LE(strays * 100, live);
  HS_EXPECT_LE(outside * 100, live);
}

inline void test_galaxies_stay_contained_small() {
  check_galaxies_stay_contained<SMALL_W, SMALL_H>();
}

inline void test_galaxies_stay_contained_default() {
  check_galaxies_stay_contained<DEFAULT_W, DEFAULT_H>();
}
