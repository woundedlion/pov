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
  static constexpr float RING_JITTER = FX<1, 1>::RING_JITTER;
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
  template <int W, int H> static const auto &system(const Galaxies<W, H> &fx) {
    return fx.particle_system;
  }
  template <int W, int H> static float orbit_speed(const Galaxies<W, H> &fx) {
    return fx.params.orbit_speed;
  }
  template <int W, int H> static void set_arms(Galaxies<W, H> &fx, int arms) {
    fx.params.arms = arms;
  }
  template <int W, int H> static void hide(Galaxies<W, H> &fx) {
    fx.params.alpha = 0.0f;
  }
};

/**
 * @brief Cores follow octahedron vertices; spawns land on their galaxy's ring
 *        with an orbital velocity, and arms cycle.
 */
inline void test_galaxies_spawn_on_ring_with_orbital_velocity() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<SMALL_W, SMALL_H> fx;
  fx.init();
  WB::set_arms(fx, 3);

  const auto &ps = WB::system(fx);
  HS_EXPECT_EQ(WB::NUM_GALAXIES, Solids::Octahedron::NUM_VERTS);
  for (int g = 0; g < WB::NUM_GALAXIES; ++g) {
    const int before = WB::arm(fx, g);
    const uint16_t index = ps.active();
    WB::emit(fx, g);
    HS_EXPECT_EQ(static_cast<int>(ps.active()), index + 1);
    HS_EXPECT_EQ(WB::arm(fx, g), (before + 1) % 3);

    const auto &p = ps.pool[index];
    const math::Vector &core = WB::core(fx, g);
    HS_EXPECT_NEAR(math::dot(core, Solids::Octahedron::vertices[g]), 1.0f,
                   1e-5f);
    HS_EXPECT_EQ(static_cast<int>(p.color_seed), g);
    HS_EXPECT_NEAR(p.position.magnitude(), 1.0f, 1e-4f);

    const float ring =
        acosf(hs::clamp(math::dot(p.position, core), -1.0f, 1.0f));
    HS_EXPECT_GE(ring, WB::RING_RADIUS * (1.0f - WB::RING_JITTER) - 1e-3f);
    HS_EXPECT_LE(ring, WB::RING_RADIUS * (1.0f + WB::RING_JITTER) + 1e-3f);

    // Tangent to the sphere and perpendicular to the core direction: a pure
    // orbit about the core, with no radial component.
    HS_EXPECT_NEAR(math::dot(p.velocity, p.position), 0.0f, 1e-4f);
    HS_EXPECT_NEAR(math::dot(p.velocity, core), 0.0f, 1e-4f);
    const float speed = p.velocity.magnitude();
    HS_EXPECT_GE(speed,
                 WB::orbit_speed(fx) * (1.0f - WB::SPEED_JITTER) - 1e-5f);
    HS_EXPECT_LE(speed,
                 WB::orbit_speed(fx) * (1.0f + WB::SPEED_JITTER) + 1e-5f);
  }
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
    const float own = math::dot(p.position, WB::core(fx, p.color_seed));
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
