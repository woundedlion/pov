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
  template <int W, int H> static void render(Galaxies<W, H> &fx) {
    Canvas canvas(fx);
    fx.draw_particles(canvas);
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
  static Color4 color(uint16_t seed) { return FX<1, 1>::star_color(seed); }
  template <int W, int H>
  static float density(const Galaxies<W, H> &fx, int index,
                       const math::Vector &position) {
    const auto &galaxy = fx.galaxies[index];
    return fx.arm_density(position, galaxy, math::dot(position, galaxy.core),
                          1.0f / tanf(fx.params.arm_pitch));
  }
};

/** @brief Every galaxy has mostly white stars and sparse blue/red accents. */
inline void test_galaxies_star_color_mix() {
  for (int galaxy = 0; galaxy < GalaxiesWhiteBox::NUM_GALAXIES; ++galaxy) {
    int white = 0;
    int blue = 0;
    int red = 0;
    for (int seed = 0; seed < 256; ++seed) {
      const auto color =
          GalaxiesWhiteBox::color(static_cast<uint16_t>((seed << 8) | galaxy));
      HS_EXPECT_EQ(color.alpha, 1.0f);
      if (color.color.r == color.color.g && color.color.g == color.color.b)
        ++white;
      else if (color.color.b > color.color.g && color.color.g > color.color.r)
        ++blue;
      else if (color.color.r > color.color.g && color.color.g > color.color.b)
        ++red;
      else
        HS_EXPECT(false, "star is white, blue or red");
    }
    HS_EXPECT_EQ(white, 216);
    HS_EXPECT_EQ(blue, 20);
    HS_EXPECT_EQ(red, 20);
  }
}

/** @brief A white star remains neutral while fading with age. */
inline void test_galaxies_white_star_stays_white() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<DEFAULT_W, DEFAULT_H> fx;
  fx.init();
  fx.updateParameter("Black Hole", 1.0f);
  auto &ps = WB::system(fx);
  const auto &galaxy = WB::galaxy(fx, 0);
  ps.spawn(galaxy.core * cosf(0.35f) + galaxy.u * sinf(0.35f), math::Vector(),
           0xc000);
  std::array<std::array<uint64_t, 3>, 3> colors{};
  const int AGES[] = {80, 520, 640};
  for (int stage = 0; stage < 3; ++stage) {
    ps.pool[0].life = ps.max_life - AGES[stage];
    WB::render(fx);
    fx.advance_display();
    for (int y = 0; y < DEFAULT_H; ++y)
      for (int x = 0; x < DEFAULT_W; ++x) {
        const Pixel &pixel = fx.get_pixel(x, y);
        colors[stage][0] += pixel.r;
        colors[stage][1] += pixel.g;
        colors[stage][2] += pixel.b;
      }
    HS_EXPECT(colors[stage][0] > 0, "star is visible");
  }
  for (const auto &color : colors) {
    HS_EXPECT_EQ(color[0], color[1]);
    HS_EXPECT_EQ(color[1], color[2]);
  }
}

/** @brief The core toggle removes the core glow and restores it. */
inline void test_galaxies_black_hole_core_toggle() {
  reset_effect_globals();
  Galaxies<DEFAULT_W, DEFAULT_H> fx;
  fx.init();
  const auto *toggle = fx.getParameters().find("Black Hole");
  HS_EXPECT(toggle != nullptr, "Black Hole is registered");
  if (!toggle)
    return;
  HS_EXPECT(toggle->is_bool(), "Black Hole is a toggle");
  HS_EXPECT_EQ(toggle->get(), 0.0f);

  uint64_t original_energy = 0;
  for (int mode : {0, 1, 0}) {
    fx.updateParameter("Black Hole", static_cast<float>(mode));
    GalaxiesWhiteBox::render(fx);
    fx.advance_display();
    uint64_t center_energy = 0;
    uint64_t rim_energy = 0;
    for (int y = 0; y < DEFAULT_H; ++y) {
      for (int x = 0; x < DEFAULT_W; ++x) {
        const math::Vector position =
            math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(x, y);
        const float distance = acosf(hs::clamp(position.x, -1.0f, 1.0f));
        const Pixel &pixel = fx.get_pixel(x, y);
        const uint64_t energy = pixel.r + pixel.g + pixel.b;
        if (distance < 0.025f)
          center_energy += energy;
        if (distance > 0.045f && distance < 0.075f)
          rim_energy += energy;
      }
    }
    if (mode == 1) {
      HS_EXPECT_EQ(center_energy, 0u);
      HS_EXPECT_EQ(rim_energy, 0u);
    } else {
      HS_EXPECT(center_energy > 0, "glowing core lights the center");
      if (original_energy != 0)
        HS_EXPECT_EQ(center_energy, original_energy);
      original_energy = center_energy;
    }
  }
}

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
    HS_EXPECT_NEAR(p.get_position().magnitude(), 1.0f, 1e-4f);

    const float ring =
        acosf(hs::clamp(math::dot(p.get_position(), core), -1.0f, 1.0f));
    HS_EXPECT_GE(ring, 0.139f);
    HS_EXPECT_LE(ring, 0.581f);
    ++radial_bins[hs::clamp(static_cast<int>((ring - 0.14f) / 0.11f), 0, 3)];
    const auto &galaxy = WB::galaxy(fx, g);
    const float azimuth = atan2f(math::dot(p.get_position(), galaxy.w),
                                 math::dot(p.get_position(), galaxy.u));
    const float expected = galaxy.phase + 2.0f * math::PI_F * galaxy.arm / 3 -
                           galaxy.spin * logf(ring / 0.58f) / tanf(0.5f);
    const float error =
        atan2f(sinf(azimuth - expected), cosf(azimuth - expected));
    const float spread = 1.0f - 0.9f * (p.color_seed >> 8) / 255.0f;
    HS_EXPECT_LE(fabsf(error), 0.22f * spread + 0.005f);

    // Tangent to the sphere and perpendicular to the core direction: a pure
    // orbit about the core, with no radial component.
    HS_EXPECT_NEAR(math::dot(p.velocity, p.get_position()), 0.0f, 1e-4f);
    HS_EXPECT_NEAR(math::dot(p.velocity, core), 0.0f, 1e-4f);
    const math::Vector outward =
        (p.get_position() * math::dot(p.get_position(), core) - core)
            .normalized();
    const float target_speed =
        WB::orbit_speed(fx, p.get_position(), outward, ring);
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
      HS_EXPECT_EQ(param.get(), 2.5f);
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
  HS_EXPECT_EQ(capacity, 12000);
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

/** @brief Packed stellar orbits stay within an eighth-column of float orbits. */
inline void test_galaxies_packed_orbits_match_float() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<DEFAULT_W, DEFAULT_H> fx;
  fx.init();
  auto &packed = WB::system(fx);
  packed.emitters.clear();
  packed.motion_cap = 0.0372f;
  packed.max_life = 801;
  static uint8_t storage[16384];
  Arena arena(storage, sizeof(storage));
  Animation::ParticleSystem<DEFAULT_W, 96, 1, 1, 6, true> full;
  full.init(arena, packed.friction, packed.gravity, 801);
  full.motion_cap = packed.motion_cap;
  for (const auto &a : packed.attractors)
    full.add_attractor(a.position, a.strength, a.kill_radius, a.event_horizon,
                       sqrtf(a.softening_sq));
  for (int i = 0; i < 96; ++i) {
    WB::emit(fx, i % WB::NUM_GALAXIES);
    const auto &p = packed.pool[i];
    full.spawn(p.get_position(), p.velocity, p.color_seed);
  }
  Canvas canvas(fx);
  float max_angle_error = 0.0f;
  for (int frame = 0; frame < 800; ++frame) {
    full.step(canvas);
    packed.step(canvas);
    HS_EXPECT_EQ(full.active(), 96);
    HS_EXPECT_EQ(packed.active(), 96);
    for (int i = 0; i < packed.active(); ++i) {
      const auto &p = packed.pool[i];
      const math::Vector a = p.get_position();
      const math::Vector b = full.pool[i].get_position();
      const float error =
          atan2f(math::cross(a, b).magnitude(), math::dot(a, b));
      HS_EXPECT(std::isfinite(error) && std::isfinite(p.velocity.magnitude()),
                "orbit remains finite");
      HS_EXPECT_NEAR(a.magnitude(), 1.0f, 3e-7f);
      HS_EXPECT_NEAR(math::dot(a, p.velocity), 0.0f, 2e-9f);
      max_angle_error = std::max(max_angle_error, error);
    }
  }
  std::printf("packed orbit maximum angular error: %.9g rad\n",
              max_angle_error);
  HS_EXPECT_LE(max_angle_error, math::RADIANS_PER_COLUMN<DEFAULT_W> / 8.0f);
}

/**
 * @brief Runs the default settings and checks that particles stay near their
 *        own core.
 * @details Checks both 96 and 288 columns because their per-frame motion caps
 *          differ.
 */
template <int W, int H>
inline void check_galaxies_stay_contained(float pitch = -1.0f) {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<W, H> fx;
  fx.init();
  if (pitch > 0.0f)
    WB::set_pitch(fx, pitch);
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
    const math::Vector position = p.get_position();
    const int owner = p.color_seed & 0xff;
    const float own = math::dot(position, WB::core(fx, owner));
    for (int g = 0; g < WB::NUM_GALAXIES; ++g) {
      if (g != owner && math::dot(position, WB::core(fx, g)) > own) {
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

/** @brief The narrowest arm pitch keeps stars within their own galaxies. */
inline void test_galaxies_min_pitch_stays_contained() {
  check_galaxies_stay_contained<SMALL_W, SMALL_H>(0.1f);
  check_galaxies_stay_contained<DEFAULT_W, DEFAULT_H>(0.1f);
}

/** @brief Live particle telemetry follows births and deaths with alpha zero. */
inline void test_galaxies_live_particle_count() {
  using WB = GalaxiesWhiteBox;
  reset_effect_globals();
  Galaxies<SMALL_W, SMALL_H> fx;
  fx.init();
  const auto *count = fx.getParameters().find("Particles");
  HS_EXPECT(count != nullptr, "Particles telemetry is registered");
  if (!count)
    return;
  HS_EXPECT(count->readonly, "Particles telemetry is read-only");
  HS_EXPECT_EQ(count->min, 0.0f);
  HS_EXPECT_EQ(count->max, 12000.0f);
  HS_EXPECT_EQ(count->get(), 0.0f);
  HS_EXPECT(fx.updateParameter("Particles", 100.0f) == ParamSetResult::READONLY,
            "client writes cannot replace live particle count");
  WB::hide(fx);
  WB::set_emission_rate(fx, 4.0f);
  for (int g = 0; g < WB::NUM_GALAXIES; ++g)
    WB::galaxy(fx, g).emission_credit = 0.0f;
  for (int frame = 0; frame < 2; ++frame) {
    pin_frame_clock(frame);
    fx.draw_frame();
    HS_EXPECT_EQ(count->get(), 24.0f * (frame + 1));
    fx.advance_display();
  }
  auto &ps = WB::system(fx);
  for (int i = 0; i < ps.active(); ++i)
    ps.pool[i].life = 1;
  pin_frame_clock(2);
  fx.draw_frame();
  HS_EXPECT_EQ(count->get(), 24.0f);
  fx.advance_display();
}
