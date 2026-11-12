/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Pool slots count only committed rings and survive alternate lookup. */
inline void test_ring_stack_pool_lifecycle_and_bounds() {
  constexpr int N = 64;
  using Pool = Scan::RingStackPool<3, N + 1>;
  alignas(std::max_align_t) uint8_t storage[Pool::BYTES];
  alignas(std::max_align_t) uint8_t value_storage[2048], table_storage[2048];
  Arena arena(storage, sizeof(storage));
  Arena values(value_storage, sizeof(value_storage)),
      tables(table_storage, sizeof(table_storage));
  typename Pool::Alternates alternates;
  alternates.allocate(values, tables, 3, N);
  Pool pool;
  pool.init_storage(arena);
  pool.begin_frame(&alternates);
  auto pending = pool.next();
  HS_EXPECT_EQ(pool.next().shift_row, pending.shift_row);
  HS_EXPECT_EQ(pool.size(), 0);
  alternates.begin_ring(1);
  for (int k = 0; k < N; ++k)
    pending.shift_row[k] = 0.01f;
  alternates.retain(1, 0.01f, 0.07f, 0.2f, pending.shift_row);
  alternates.retain(33, -0.08f, 0.01f, 0.8f, pending.shift_row);
  alternates.finish_ring(1);
  pending.shift_row[N] = pending.shift_row[0];
  const auto basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  pool.commit(1, 0.7f, N, basis, 0.5f, 0.03f, pending.shift_row, N, 0.0f,
              nullptr);
  HS_EXPECT_EQ(pool.size(), 1);
  HS_EXPECT_EQ(pool.slot_of(0), Pool::CULLED);
  HS_EXPECT_EQ(pool.slot_of(1), 0);
  float lo, hi;
  pool.shapes().knot_range(0, 1, lo, hi);
  HS_EXPECT_NEAR(lo, 0.01f, 1e-7f);
  HS_EXPECT_NEAR(hi, 0.07f, 1e-7f);
  pool.shapes().knot_range(0, 33, lo, hi);
  HS_EXPECT_NEAR(lo, -0.08f, 1e-7f);
  HS_EXPECT_NEAR(hi, 0.01f, 1e-7f);
  lo = hi = 0.3f;
  alternates.widen(1, 32, lo, hi);
  HS_EXPECT_EQ(lo, 0.3f);
  HS_EXPECT_EQ(hi, 0.3f);
  HS_EXPECT_FALSE(pool.shapes().table_ready());
  pool.shapes().table_built();
  HS_EXPECT_TRUE(pool.shapes().table_ready());
  pool.release();
  HS_EXPECT_EQ(pool.size(), 0);
  HS_EXPECT_EQ(pool.destroyed, 1);
  pool.begin_frame();
  HS_EXPECT_EQ(pool.slot_of(1), Pool::CULLED);
  HS_EXPECT_FALSE(pool.shapes().table_ready());
  HS_EXPECT_EQ(alternates.count, 2);
  HS_EXPECT_TRUE(alternates.fits(N - 2));
  HS_EXPECT_FALSE(alternates.fits(N - 1));
  alternates.reset();
  HS_EXPECT_EQ(alternates.count, 0);
  HS_EXPECT_FALSE(alternates.table_ready);
}

/** @brief Two-window scan matches independent draws with manually installed rows. */
inline void test_swept_ring_stack_matches_independent_windows() {
  constexpr int W = 32, H = 16, N = 16;
  auto render = [=](bool swept) {
    using Pool = Scan::RingStackPool<2, W + 1>;
    alignas(std::max_align_t) uint8_t storage[Pool::BYTES];
    alignas(std::max_align_t) uint8_t value_storage[1024], table_storage[1024];
    Arena arena(storage, sizeof(storage));
    Arena values(value_storage, sizeof(value_storage)),
        tables(table_storage, sizeof(table_storage));
    typename Pool::Alternates alternates;
    alternates.allocate(values, tables, 2, W);
    Pool pool;
    pool.init_storage(arena);
    pool.begin_frame(&alternates);
    const auto basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    for (int i = 0; i < 2; ++i) {
      auto pending = pool.next();
      alternates.begin_ring(i);
      for (int k = 0; k < N; ++k) {
        const float early = 0.03f * sinf(math::TWO_PI_F * k / N);
        const float late = -early;
        alternates.retain(k, early, late, 0.25f, pending.shift_row);
        pending.hue_row[k] = Pixel(30000 + i * 10000, 20000, 40000);
      }
      pending.shift_row[N] = pending.shift_row[0];
      pending.hue_row[N] = pending.hue_row[0];
      alternates.finish_ring(i);
      pool.commit(i, 0.8f, N, basis, 0.6f + i * 0.4f, 0.1f, pending.shift_row,
                  N, 0.0f, nullptr);
    }
    Scan::DistortedRingStack::CandidateTable<W, H> table;
    Scan::HueRowShader<Pool> shader{&pool};
    int reshaded = 0;
    int finished = 0;
    struct Reshade {
      Pool &pool;
      int &reshaded;
      int &finished;
      int slot = 0;
      int knots = 0;

      void begin(int ring, int s) {
        HS_EXPECT_EQ(s, ring);
        slot = s;
        knots = 0;
        ++reshaded;
      }
      void operator()(int k, float shift) {
        HS_EXPECT_EQ(k, knots);
        HS_EXPECT_NEAR(shift, -0.03f * sinf(math::TWO_PI_F * k / N), 1e-7f);
        ++knots;
      }
      void finish() {
        HS_EXPECT_EQ(knots, N);
        HS_EXPECT_EQ(pool.shift_row(slot)[N], pool.shift_row(slot)[0]);
        ++finished;
      }
    };
    const auto result =
        render_stub_frame<W, H>([&](auto &pipeline, Canvas &canvas) {
          if (swept) {
            Reshade reshade{pool, reshaded, finished};
            Scan::SweptRingStack::scan_two_windows<W, H>(
                pipeline, canvas, 2, pool, alternates, table, shader, reshade);
            HS_EXPECT_EQ(reshade.knots, 0);
          } else {
            for (int window_index = 0; window_index < 2; ++window_index) {
              if (window_index)
                for (int i = 0; i < 2; ++i) {
                  float *row = pool.shift_row(i);
                  for (int k = 0; k < N; ++k)
                    row[k] = -0.03f * sinf(math::TWO_PI_F * k / N);
                  row[N] = row[0];
                }
              TimeWindowScope window(canvas, window_index * W / 4,
                                     (window_index + 1) * W / 4);
              Scan::DistortedRingStack::draw<W, H>(
                  pipeline, canvas, 2, pool.shapes(), pool.slot_map(),
                  pool.size(), table, shader);
            }
            pool.release();
          }
        });
    if (swept) {
      HS_EXPECT_EQ(reshaded, 2);
      HS_EXPECT_EQ(finished, 2);
      HS_EXPECT_EQ(alternates.count, 0);
      HS_EXPECT_EQ(pool.size(), 0);
      HS_EXPECT_EQ(pool.destroyed, 2);
      HS_EXPECT_EQ(alternates.stack_flushes, 1);
    }
    return result;
  };
  const auto expected = render(false), actual = render(true);
  HS_EXPECT_EQ(actual.size(), expected.size());
  int lit = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    lit += actual[i].r != 0 || actual[i].g != 0 || actual[i].b != 0;
    HS_EXPECT_EQ(actual[i].r, expected[i].r);
    HS_EXPECT_EQ(actual[i].g, expected[i].g);
    HS_EXPECT_EQ(actual[i].b, expected[i].b);
  }
  HS_EXPECT_GT(lit, 0);
}

/** @brief Boundary masks cover sampled displaced geometry and its supported azimuths. */
inline void test_swept_ring_stack_boundary_is_conservative() {
  constexpr int W = 96, MAX_KNOTS = 96;
  constexpr double PI_EXACT = 3.14159265358979323846;
  constexpr int POLAR_SAMPLES = 25, AZIMUTH_SAMPLES = 33;
  uint32_t seed = 0x5e27a19bu;
  auto unit = [&]() {
    seed = 1664525u * seed + 1013904223u;
    return static_cast<float>(seed >> 8) / 16777216.0f;
  };
  int supported_hits = 0, unmarked_knots = 0;
  for (int trial = 0; trial < 64; ++trial) {
    const math::Vector axis =
        math::Vector(2.0f * unit() - 1.0f, 0.2f + unit(), 2.0f * unit() - 1.0f)
            .normalized();
    const auto basis = math::make_basis(
        math::make_rotation(axis, math::TWO_PI_F * unit()), math::Y_AXIS);
    const int n = MAX_KNOTS / (1 << (trial % 4));
    const float poles[] = {
        0.0f, 0.01f, 0.4f, math::PI_F / 2.0f, math::PI_F - 0.01f, math::PI_F};
    const float theta = trial < 6 ? poles[trial] : math::PI_F * unit();
    const float displacement = trial % 7 == 0 ? 0.0f : 1.2f * unit();
    const float thickness = 0.0005f + 0.09f * unit();
    const float motion = trial % 3 == 0 ? 0.0f : 0.3f * unit();
    const int support = trial % 5;
    const float cell = math::TWO_PI_F / n;
    const SDF::KnotStep step{cosf(theta), sinf(theta), cosf(cell), sinf(cell)};
    uint8_t marked[MAX_KNOTS];
    Scan::SweptRingStack::mark_arm_boundary<W>(basis, step, theta, n,
                                               displacement, thickness, motion,
                                               support, marked);

    const double low = std::max(0.0, double(theta) - displacement - thickness);
    const double high =
        std::min(PI_EXACT, double(theta) + displacement + thickness);
    const double min_sine =
        std::max(std::sqrt(double(SDF::DistortedRing::POLE_SIN2_FLOOR)),
                 std::min(std::sin(low), std::sin(high)));
    const int search_limit =
        std::min(SDF::DistortedRing::MAX_SEARCH_CELLS, n / 2 + 1) + 1;
    const double search = std::min(
        double(search_limit),
        std::ceil(thickness / (min_sine * (2.0 * PI_EXACT / n))) + 1.0);
    const double reach =
        std::min(PI_EXACT, (search + support) * (2.0 * PI_EXACT / n));
    const double epsilon =
        std::sin(std::min(PI_EXACT / 2.0, motion + 2.0 * PI_EXACT / W));
    double polar_cos[POLAR_SAMPLES], polar_sin[POLAR_SAMPLES];
    for (int p = 0; p < POLAR_SAMPLES; ++p) {
      const double polar = low + (high - low) * p / (POLAR_SAMPLES - 1);
      polar_cos[p] = std::cos(polar);
      polar_sin[p] = std::sin(polar);
    }
    for (int k = 0; k < n; ++k) {
      bool touches = false;
      for (int a = 0; a < AZIMUTH_SAMPLES && !touches; ++a) {
        const double azimuth = 2.0 * PI_EXACT * k / n - reach +
                               2.0 * reach * a / (AZIMUTH_SAMPLES - 1);
        const double radial =
            basis.u.z * std::cos(azimuth) + basis.w.z * std::sin(azimuth);
        for (int p = 0; p < POLAR_SAMPLES; ++p) {
          const double z = basis.v.z * polar_cos[p] + radial * polar_sin[p];
          if (std::abs(z) <= epsilon) {
            touches = true;
            break;
          }
        }
      }
      if (touches) {
        ++supported_hits;
        HS_EXPECT_NE(marked[k], 0);
      }
      unmarked_knots += marked[k] == 0;
    }
  }
  HS_EXPECT_GT(supported_hits, 0);
  HS_EXPECT_GT(unmarked_knots, 0);
}

/** @brief Sweep times follow half-turn columns within each clamped window. */
inline void test_sweep_fraction_matches_half_turn_windows() {
  constexpr double PI_EXACT = 3.14159265358979323846;
  for (int k = 0; k < 1001; ++k) {
    const double angle = 2.0 * PI_EXACT * (k + 0.173) / 1001;
    const math::Vector p(static_cast<float>(std::cos(angle)), 0.0f,
                         static_cast<float>(std::sin(angle)));
    for (float middle : {0.25f, 0.75f}) {
      double u =
          std::atan2(static_cast<double>(p.z), static_cast<double>(p.x)) /
          PI_EXACT;
      u -= std::floor(u);
      if (std::fabs(std::fabs(u - middle) - 0.5) < 1.3e-3)
        continue;
      if (u - middle >= 0.5)
        u -= 1.0;
      else if (middle - u > 0.5)
        u += 1.0;
      const float expected = static_cast<float>(std::clamp(u, 0.0, 1.0));
      const float actual = Scan::SweptRingStack::sweep_fraction(p, middle);
      HS_EXPECT_NEAR(actual, expected, 1.3e-3f);
      HS_EXPECT_GE(actual, 0.0f);
      HS_EXPECT_LE(actual, 1.0f);
    }
  }
}

/** @brief Alternate storage preserves caller allocation order and aligned footprints. */
inline void test_knot_alternates_interleaved_storage_footprint() {
  constexpr int KNOTS = 65, SLOTS = 3, CAPACITY = 128;
  using Alternates = Scan::KnotAlternates<KNOTS>;
  alignas(std::max_align_t) uint8_t value_storage[2048], table_storage[2048];
  Arena values(value_storage, sizeof(value_storage)),
      tables(table_storage, sizeof(table_storage));
  values.allocate_n<uint8_t>(3);
  tables.allocate_n<uint8_t>(3);
  Alternates alternates;
  uint64_t *caller = nullptr;
  alternates.allocate(values, tables, SLOTS, CAPACITY, [&](Arena &arena) {
    HS_EXPECT_EQ(arena.get_high_water_mark(),
                 size_t{4 + SLOTS * Alternates::WORDS * 6});
    caller = arena.allocate_n<uint64_t>(3);
  });
  HS_EXPECT_TRUE(reinterpret_cast<uint8_t *>(alternates.starts) <
                 reinterpret_cast<uint8_t *>(caller));
  HS_EXPECT_EQ(alternates.boundary, reinterpret_cast<uint8_t *>(caller + 3));
  constexpr auto FOOTPRINTS = [] {
    ArenaFootprint packed, table;
    packed.add<uint8_t>(3);
    table.add<uint8_t>(3);
    Alternates::add_footprint(
        packed, table, SLOTS, CAPACITY,
        [](ArenaFootprint &footprint) { footprint.add<uint64_t>(3); });
    return std::pair{packed.bytes, table.bytes};
  }();
  HS_EXPECT_EQ(values.get_high_water_mark(), FOOTPRINTS.first);
  HS_EXPECT_EQ(tables.get_high_water_mark(), FOOTPRINTS.second);
  ArenaFootprint plain_values, plain_table;
  Alternates::add_footprint(plain_values, plain_table, SLOTS, CAPACITY);
  HS_EXPECT_EQ(plain_values.bytes, size_t{CAPACITY * sizeof(float)});
  HS_EXPECT_EQ(plain_table.bytes,
               size_t{SLOTS * Alternates::WORDS * 6 + KNOTS + 1});
}
