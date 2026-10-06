/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// World/Screen::Trails: quantization, capacity replacement and ttl lifecycle
// ============================================================================

/**
 * @brief Verifies plot() passes the frame through and stores it, and flush()
 *        decodes the int16 entry and re-emits it (then ages) within the
 *        quantization error bound (< 1/32767).
 */
inline void test_world_trails_int16_quantization_roundtrip() {
  constexpr int Cap = 8;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  const math::Vector v0 = math::Vector(0.3f, -0.6f, 0.74f).normalized();

  // plot() passes the frame through and (age=0 -> ttl=10>0) stores it.
  int passthru = 0;
  trails.plot(
      v0, Pixel(100, 100, 100), 0.0f, 1.0f,
      [&](const math::Vector &, const Pixel &, float, float) { ++passthru; });
  HS_EXPECT_EQ(passthru, 1);
  HS_EXPECT_EQ(trails.size(), (size_t)1);

  // flush() decodes the int16 entry and re-emits it, then ages (10->9, alive).
  math::Vector decoded(0, 0, 0);
  int emitted = 0;
  auto trail = [](const math::Vector &, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(WorldTrailFn(trail), 1.0f,
               [&](const math::Vector &v, const Pixel &, float, float) {
                 decoded = v;
                 ++emitted;
               });
  HS_EXPECT_EQ(emitted, 1);
  // Quantization is v*32767 truncated to int16 then *(1/32767): |err| < 1/32767.
  HS_EXPECT_NEAR(decoded.x, v0.x, 1.0f / 32767.0f + 1e-7f);
  HS_EXPECT_NEAR(decoded.y, v0.y, 1.0f / 32767.0f + 1e-7f);
  HS_EXPECT_NEAR(decoded.z, v0.z, 1.0f / 32767.0f + 1e-7f);
}

/**
 * @brief Verifies a component pushed past the unit cube saturates rather than
 *        overflowing int16 and wrapping to a garbage point on the far side of
 *        the sphere.
 * @details encode() clamps to [-1, 1] before quantizing.
 */
inline void test_world_trails_clamps_out_of_range() {
  constexpr int Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  // 1.8*32767 and -1.5*32767 overflow int16; both must saturate to +/-1.
  const math::Vector v = math::Vector(1.8f, 0.5f, -1.5f);
  trails.plot(v, Pixel(1, 1, 1), 0.0f, 1.0f,
              [](const math::Vector &, const Pixel &, float, float) {});

  math::Vector decoded(0, 0, 0);
  auto trail = [](const math::Vector &, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(
      WorldTrailFn(trail), 1.0f,
      [&](const math::Vector &d, const Pixel &, float, float) { decoded = d; });
  HS_EXPECT_NEAR(decoded.x, 1.0f, 1e-3f);  // saturated, not wrapped negative
  HS_EXPECT_NEAR(decoded.y, 0.5f, 1e-3f);  // in range: untouched
  HS_EXPECT_NEAR(decoded.z, -1.0f, 1e-3f); // saturated
}

/** @brief Capacity eviction retains a bounded live set and accepts new points. */
inline void test_world_trails_capacity_evicts_one_slot() {
  constexpr int Cap = 4;
  constexpr int Overflow = 3;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/100);
  trails.init_storage(arena);

  auto noop = [](const math::Vector &, const Pixel &, float, float) {};
  math::Vector pushed[Cap + Overflow];
  for (int i = 0; i < Cap + Overflow; ++i) {
    pushed[i] =
        math::Vector(static_cast<float>(i + 1), 1.0f, 0.5f).normalized();
    trails.plot(pushed[i], Pixel(1, 1, 1), 0.0f, 1.0f, noop);
  }
  HS_EXPECT_EQ(trails.size(), (size_t)Cap);

  auto trail = [](const math::Vector &, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  std::vector<math::Vector> decoded;
  trails.flush(WorldTrailFn(trail), 1.0f,
               [&](const math::Vector &v, const Pixel &, float, float) {
                 decoded.push_back(v);
               });
  HS_EXPECT_SIZE_OR_RETURN(decoded, Cap);
  constexpr float TOLERANCE = 1.0f / 32767.0f + 1e-7f;
  auto matches = [&](const math::Vector &want) {
    int n = 0;
    for (const math::Vector &v : decoded)
      if (std::abs(v.x - want.x) <= TOLERANCE &&
          std::abs(v.y - want.y) <= TOLERANCE &&
          std::abs(v.z - want.z) <= TOLERANCE)
        ++n;
    return n;
  };
  // Each overflow push replaces the newest slot; the oldest Cap - 1 survive.
  for (int i = 0; i < Cap + Overflow; ++i) {
    const bool kept = i < Cap - 1 || i == Cap + Overflow - 1;
    HS_EXPECT_EQ(matches(pushed[i]), (kept ? 1 : 0));
  }
}

/**
 * @brief Verifies each flush decrements an entry's ttl and the entry is popped
 *        once ttl reaches 0.
 */
inline void test_world_trails_ttl_expiry() {
  constexpr int Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/2);
  trails.init_storage(arena);

  auto noop = [](const math::Vector &, const Pixel &, float, float) {};
  trails.plot(math::Vector(0, 1, 0), Pixel(1, 1, 1), 0.0f, 1.0f,
              noop); // ttl = 2
  HS_EXPECT_EQ(trails.size(), (size_t)1);

  auto trail = [](const math::Vector &, float) {
    return Color4(Pixel(1, 1, 1), 1.0f);
  };
  auto sink = [](const math::Vector &, const Pixel &, float, float) {};
  trails.flush(WorldTrailFn(trail), 1.0f, sink); // ttl 2 -> 1, alive
  HS_EXPECT_EQ(trails.size(), (size_t)1);
  trails.flush(WorldTrailFn(trail), 1.0f, sink); // ttl 1 -> 0, popped
  HS_EXPECT_EQ(trails.size(), (size_t)0);
}

/**
 * @brief Verifies set_lifetime() caps buffered World::Trails ttl and restarts
 *        fade progress at zero for points above the new lifetime.
 */
inline void test_world_trails_set_lifetime_caps_ttl() {
  constexpr int Cap = 8;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  const math::Vector v0 = math::Vector(0.3f, -0.6f, 0.74f).normalized();
  trails.plot(
      v0, Pixel(100, 100, 100), 0.0f, 1.0f,
      [](const math::Vector &, const Pixel &, float, float) {}); // ttl = 10

  trails.set_lifetime(2);

  float captured_t = -999.0f;
  auto trail = [&](const math::Vector &, float t) {
    captured_t = t;
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(WorldTrailFn(trail), 1.0f,
               [](const math::Vector &, const Pixel &, float, float) {});

  HS_EXPECT_EQ(captured_t, 0.0f);
  trails.flush(WorldTrailFn(trail), 1.0f,
               [](const math::Vector &, const Pixel &, float, float) {});
  HS_EXPECT_NEAR(captured_t, 0.5f, 1e-6f);
  captured_t = -1.0f;
  trails.flush(WorldTrailFn(trail), 1.0f,
               [](const math::Vector &, const Pixel &, float, float) {});
  HS_EXPECT_EQ(captured_t, -1.0f);
}

/**
 * @brief Verifies flush() reclaims a dead mid-buffer item (heterogeneous TTLs),
 *        not just dead items at the head.
 * @details A short-lived point buffered behind a long-lived older one dies in
 *          the middle of the array; flush()'s swap-remove cull frees its slot, so
 *          live capacity is preserved and the next plot() fills the freed slot
 *          without evicting any live point.
 */
inline void test_world_trails_midbuffer_expiry_reclaims_slot() {
  constexpr int Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<Cap> trails(/*lifetime=*/100);
  trails.init_storage(arena);

  // Orthogonal/antipodal unit vectors so int16-quantized decodes stay trivially
  // identifiable by dot product.
  const math::Vector p0(1, 0, 0), p1(0, 1, 0), p2(0, 0, 1), p3(-1, 0, 0),
      p4(0, -1, 0);
  auto noop = [](const math::Vector &, const Pixel &, float, float) {};

  trails.plot(p0, Pixel(1, 1, 1), 0.0f, 1.0f,
              noop); // ttl 100 — oldest, long-lived
  trails.plot(p1, Pixel(1, 1, 1), 99.0f, 1.0f,
              noop); // ttl 1 — dies on next flush
  trails.plot(p2, Pixel(1, 1, 1), 0.0f, 1.0f, noop); // ttl 100
  trails.plot(p3, Pixel(1, 1, 1), 0.0f, 1.0f, noop); // ttl 100
  HS_EXPECT_EQ(trails.size(), (size_t)Cap);          // [p0, p1, p2, p3] — full

  int live_drawn = 0;
  bool saw_p0 = false;
  auto trail = [&](const math::Vector &v, float) {
    ++live_drawn;
    if (math::dot(v, p0) > 0.9f)
      saw_p0 = true;
    return Color4(Pixel(1, 1, 1), 1.0f);
  };
  // p1 is drawn this frame (ttl 1, its final render) and then ages to 0. The
  // front (p0) is still alive; the swap-remove cull reclaims p1 mid-buffer.
  trails.flush(WorldTrailFn(trail), 1.0f, noop);
  HS_EXPECT_EQ(live_drawn, 4); // all 4 still live at emit time
  HS_EXPECT_TRUE(saw_p0);      // the oldest live point still draws
  HS_EXPECT_EQ(trails.size(), (size_t)(Cap - 1)); // dead p1's slot reclaimed

  // The reclaimed slot admits p4 without evicting any live point.
  trails.plot(p4, Pixel(1, 1, 1), 0.0f, 1.0f, noop);

  bool saw_p0_after = false, saw_p2_after = false, saw_p3_after = false,
       saw_p4_after = false;
  int drawn_after = 0;
  auto trail2 = [&](const math::Vector &v, float) {
    ++drawn_after;
    if (math::dot(v, p0) > 0.9f)
      saw_p0_after = true;
    if (math::dot(v, p2) > 0.9f)
      saw_p2_after = true;
    if (math::dot(v, p3) > 0.9f)
      saw_p3_after = true;
    if (math::dot(v, p4) > 0.9f)
      saw_p4_after = true;
    return Color4(Pixel(1, 1, 1), 1.0f);
  };
  trails.flush(WorldTrailFn(trail2), 1.0f, noop);
  HS_EXPECT_EQ(drawn_after, 4);
  HS_EXPECT_TRUE(saw_p0_after);
  HS_EXPECT_TRUE(saw_p2_after);
  HS_EXPECT_TRUE(saw_p3_after);
  HS_EXPECT_TRUE(saw_p4_after);
}

/**
 * @brief Verifies the Screen::Trails store / emit / decay lifecycle.
 * @details Screen::Trails stores float DecayPixels with no int16 quantization
 *          (that path is World::Trails-specific).
 */
inline void test_screen_trails_store_emit_decay() {
  constexpr int W = 32, MAXP = 16;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<MAXP> trails(/*lifetime=*/3);
  trails.init_storage(arena);

  hs_test::StubEffect fx(W,
                         8); // flush takes a Canvas& (unused by Screen::Trails)
  Canvas c(fx);

  // age=1 (0<age<lifetime): forwarded this frame AND stored.
  int passthru = 0;
  float fwd_age = -1.0f;
  trails.plot(10.0f, 4.0f, Pixel(5, 6, 7), 1.0f, 1.0f,
              [&](float, float, const Pixel &, float a, float) {
                fwd_age = a;
                ++passthru;
              });
  HS_EXPECT_EQ(passthru, 1);
  HS_EXPECT_NEAR(fwd_age, 1.0f, 1e-6f); // forwarded verbatim

  auto trail = [](float, float, float) { return Color4(Pixel(9, 9, 9), 1.0f); };
  // Stored ttl = lifetime - age = 2. Each flush emits then decays (--ttl).
  int emitted = 0;
  auto counting_sink = [&](float, float, const Pixel &, float, float) {
    ++emitted;
  };
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // emit, ttl 2->1
  HS_EXPECT_EQ(emitted, 1);
  emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // emit, ttl 1->0
  HS_EXPECT_EQ(emitted, 1);
  emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // decayed out
  HS_EXPECT_EQ(emitted, 0);
}

/**
 * @brief Verifies Screen::Trails clamps fade progress for a negative age.
 * @details Negative age seeds ttl above lifetime, giving unclamped t = -0.8.
 */
inline void test_screen_trails_negative_age_clamps_t() {
  constexpr int W = 32, MAXP = 8;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<MAXP> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  hs_test::StubEffect fx(W, 8);
  Canvas c(fx);
  trails.plot(4.0f, 2.0f, Pixel(100, 100, 100), -8.0f, 1.0f,
              [](float, float, const Pixel &, float, float) {}); // ttl = 18

  float captured_t = -999.0f;
  auto trail = [&](float, float, float t) {
    captured_t = t;
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [](float, float, const Pixel &, float, float) {});

  HS_EXPECT_EQ(captured_t, 0.0f);
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [](float, float, const Pixel &, float, float) {});
  HS_EXPECT_EQ(captured_t, 0.0f);
  for (int i = 0; i < 16; ++i)
    trails.flush(c, ScreenTrailFn(trail), 1.0f,
                 [](float, float, const Pixel &, float, float) {});
  captured_t = -1.0f;
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [](float, float, const Pixel &, float, float) {});
  HS_EXPECT_EQ(captured_t, -1.0f);
}

/**
 * @brief Verifies Screen::Trails forwards an already-aged emission, matching
 *        World::Trails.
 * @details A point with 0 < age < lifetime is both passed through the current
 *          frame and seeded into storage. A point at/past lifetime is still
 *          forwarded, but ttl<=0 keeps it out of storage.
 */
inline void test_screen_trails_forwards_aged_emission() {
  constexpr int W = 32, MAXP = 16;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<MAXP> trails(/*lifetime=*/5);
  trails.init_storage(arena);

  hs_test::StubEffect fx(W, 8);
  Canvas c(fx);
  auto trail = [](float, float, float) { return Color4(Pixel(9, 9, 9), 1.0f); };

  // 0 < age < lifetime: forwarded this frame and seeded.
  int fwd = 0;
  trails.plot(3.0f, 4.0f, Pixel(1, 2, 3), 2.0f, 1.0f,
              [&](float, float, const Pixel &, float, float) { ++fwd; });
  HS_EXPECT_EQ(fwd, 1);

  // age == lifetime: dead point is still forwarded, but ttl<=0 so not seeded.
  fwd = 0;
  trails.plot(7.0f, 4.0f, Pixel(1, 2, 3), 5.0f, 1.0f,
              [&](float, float, const Pixel &, float, float) { ++fwd; });
  HS_EXPECT_EQ(fwd, 1);

  // Only the live (age=2) point was seeded -> exactly one flush emission.
  int emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [&](float, float, const Pixel &, float, float) { ++emitted; });
  HS_EXPECT_EQ(emitted, 1);
}

/**
 * @brief Capacity overflow replaces the last slot and preserves the others.
 */
inline void test_screen_trails_at_capacity_replaces_last_slot() {
  constexpr int W = 32, MAXP = 4;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<MAXP> trails(/*lifetime=*/100);
  trails.init_storage(arena);

  hs_test::StubEffect fx(W, 8);
  Canvas c(fx);
  auto noop = [](float, float, const Pixel &, float, float) {};
  for (int i = 0; i < MAXP + 1; ++i)
    trails.plot(static_cast<float>(i + 1), 2.0f, Pixel(1, 1, 1), 0.0f, 1.0f,
                noop);

  auto trail = [](float, float, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  std::vector<float> emitted;
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [&](float x, float, const Pixel &, float, float) {
                 emitted.push_back(x);
               });

  HS_EXPECT_SIZE_OR_RETURN(emitted, MAXP);
  std::sort(emitted.begin(), emitted.end());
  for (int i = 0; i < MAXP - 1; ++i)
    HS_EXPECT_EQ(emitted[i], static_cast<float>(i + 1));
  HS_EXPECT_EQ(emitted[MAXP - 1], static_cast<float>(MAXP + 1));
}

/** @brief Shortening screen trails caps remaining lifetime and fade progress. */
inline void test_screen_trails_set_lifetime_caps_ttl() {
  uint8_t buf[Filter::Screen::Trails<4>::STORAGE_BYTES];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<4> trails(10);
  trails.init_storage(arena);
  hs_test::StubEffect fx(32, 8);
  Canvas canvas(fx);
  auto pass = [](float, float, const Pixel &, float, float) {};
  trails.plot(3.0f, 4.0f, Pixel(1, 2, 3), 0.0f, 1.0f, pass);
  trails.set_lifetime(2);
  std::vector<float> ages;
  auto trail = [&](float, float, float t) {
    ages.push_back(t);
    return Color4(Pixel(1, 2, 3), 1.0f);
  };
  for (int frame = 0; frame < 3; ++frame)
    trails.flush(canvas, ScreenTrailFn(trail), 1.0f, pass);
  HS_EXPECT_SIZE_OR_RETURN(ages, 2);
  HS_EXPECT_EQ(ages[0], 0.0f);
  HS_EXPECT_EQ(ages[1], 0.5f);
}

/** @brief Pins screen seeding and both domains' emission alpha floors. */
inline void test_trails_alpha_gates() {
  uint8_t screen_buf[Filter::Screen::Trails<4>::STORAGE_BYTES];
  uint8_t world_buf[Filter::World::Trails<4>::STORAGE_BYTES];
  Arena screen_arena(screen_buf, sizeof(screen_buf));
  Arena world_arena(world_buf, sizeof(world_buf));
  Filter::Screen::Trails<4> screen(10);
  Filter::World::Trails<4> world(10);
  screen.init_storage(screen_arena);
  world.init_storage(world_arena);
  hs_test::StubEffect fx(32, 8);
  Canvas canvas(fx);
  int forwards = 0, emitted = 0;
  auto pass = [&](float, float, const Pixel &, float, float) { ++forwards; };
  screen.plot(1, 2, Pixel(1, 2, 3), 0, 0, pass);
  screen.plot(1, 2, Pixel(1, 2, 3), 0, Filter::TRAIL_EMIT_ALPHA_FLOOR, pass);
  HS_EXPECT_EQ(forwards, 2);
  auto visible = [](float, float, float) {
    return Color4(Pixel(1, 2, 3), 1.0f);
  };
  auto screen_emit = [&](float, float, const Pixel &, float, float) {
    ++emitted;
  };
  screen.flush(canvas, ScreenTrailFn(visible), 1, screen_emit);
  HS_EXPECT_EQ(emitted, 0);
  screen.plot(1, 2, Pixel(1, 2, 3), 0, 1, pass);
  world.plot(math::X_AXIS, Pixel(1, 2, 3), 0, 1,
             [](const math::Vector &, const Pixel &, float, float) {});
  for (float alpha : {0.0f, Filter::TRAIL_EMIT_ALPHA_FLOOR, 1.0f}) {
    emitted = 0;
    auto screen_trail = [=](float, float, float) {
      return Color4(Pixel(1, 2, 3), alpha);
    };
    auto world_trail = [=](const math::Vector &, float) {
      return Color4(Pixel(1, 2, 3), alpha);
    };
    screen.flush(canvas, ScreenTrailFn(screen_trail), 1, screen_emit);
    world.flush(
        WorldTrailFn(world_trail), 1,
        [&](const math::Vector &, const Pixel &, float, float) { ++emitted; });
    HS_EXPECT_EQ(emitted, alpha == 1.0f ? 2 : 0);
  }
}

/** @brief Drains a pipeline carrying history in both domains, world first. */
inline void test_mixed_domain_flush_drains_both_buffers() {
  constexpr int W = 32, H = 16, CAP = 8, MAXP = 64, LIFETIME = 3;
  static uint8_t buf[CAP * 16 + MAXP * 32];
  Arena arena(buf, sizeof(buf));
  using WorldTrails = Filter::World::Trails<CAP>;
  using ScreenTrails = Filter::Screen::Trails<MAXP>;

  Pipeline<W, H, WorldTrails, ScreenTrails> pipe{WorldTrails(LIFETIME),
                                                 ScreenTrails(LIFETIME)};
  pipe.get<WorldTrails>().init_storage(arena);
  pipe.get<ScreenTrails>().init_storage(arena);

  hs_test::StubEffect fx(W, H);
  Canvas c(fx);

  pipe.plot(c, math::Vector(0.3f, -0.6f, 0.74f).normalized(), Pixel(1, 2, 3),
            0.0f, 1.0f);
  HS_EXPECT_EQ(pipe.get<WorldTrails>().size(), (size_t)1);

  int world_emits = 0, screen_emits = 0;
  auto world_trail = [&](const math::Vector &, float) {
    ++world_emits;
    return Color4(Pixel(9, 9, 9), 1.0f);
  };
  auto screen_trail = [&](float, float, float) {
    ++screen_emits;
    return Color4(Pixel(9, 9, 9), 1.0f);
  };

  // The original point plus one re-emission per frame, seeded at ttl 3, 2, 1.
  constexpr int EXPECTED_SCREEN[LIFETIME] = {2, 3, 4};
  for (int frame = 0; frame < LIFETIME; ++frame) {
    world_emits = screen_emits = 0;
    pipe.flush(c, WorldTrailFn(world_trail), ScreenTrailFn(screen_trail), 1.0f);
    HS_EXPECT_EQ(world_emits, 1);
    HS_EXPECT_EQ(screen_emits, EXPECTED_SCREEN[frame]);
  }

  // ttl reached 0 on the last pass: the world buffer aged out.
  HS_EXPECT_EQ(pipe.get<WorldTrails>().size(), (size_t)0);

  world_emits = screen_emits = 0;
  pipe.flush(c, WorldTrailFn(world_trail), ScreenTrailFn(screen_trail), 1.0f);
  HS_EXPECT_EQ(world_emits, 0);
  HS_EXPECT_EQ(screen_emits, 0);
}
