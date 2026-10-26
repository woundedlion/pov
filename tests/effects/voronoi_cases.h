/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Voronoi site sampling, union candidates and segmented render parity.
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for Voronoi's seeded sites and adaptive block floor.
 */
struct VoronoiWhiteBox {
  using VO = Voronoi<DEFAULT_W, DEFAULT_H>;
  static constexpr int MAX_SITES = VO::MAX_SITES;

  /** @brief Adaptive block floor at render resolution W x H. */
  template <int W, int H> static constexpr int coherence_block_min() {
    return Voronoi<W, H>::COHERENCE_BLOCK_MIN;
  }

  /** @brief Number of currently seeded sites. */
  template <int W, int H> static size_t site_count(const Voronoi<W, H> &v) {
    return v.sites_buffer.size();
  }
  /** @brief Spin axis of seeded site @p i. */
  template <int W, int H>
  static math::Vector site_axis(const Voronoi<W, H> &v, size_t i) {
    return v.sites_buffer[i].axis;
  }

  /** @brief Position of seeded site @p i. */
  template <int W, int H>
  static math::Vector site_position(const Voronoi<W, H> &v, size_t i) {
    return v.sites_buffer[i].pos;
  }

  /** @brief Installs deterministic sites and index-coded colors. */
  template <int W, int H>
  static void set_sites(Voronoi<W, H> &v, std::span<const math::Vector> sites) {
    v.sites_buffer.clear();
    for (size_t i = 0; i < sites.size(); ++i) {
      v.sites_buffer.push_back(
          {sites[i], math::Vector(1, 0, 0),
           Color4(Pixel(static_cast<uint16_t>(i + 1), 0, 0))});
    }
    v.current_num_sites = static_cast<int>(sites.size());
    v.params.num_sites = static_cast<float>(sites.size());
    v.params.speed = 0.0f;
    v.params.sharpness = 0.0f;
    v.params.border_thickness = 0.0f;
  }

  /** @brief Installs sites spinning about their own axes, index-coded. */
  template <int W, int H>
  static void
  set_moving_sites(Voronoi<W, H> &v, std::span<const math::Vector> sites,
                   std::span<const math::Vector> axes, float speed) {
    set_sites(v, sites);
    for (size_t i = 0; i < sites.size(); ++i)
      v.sites_buffer[i].axis = axes[i];
    v.params.speed = speed;
  }

  /** @brief Forces both endpoint pairs to be searched on the next draw. */
  template <int W, int H> static void forget_corners(Voronoi<W, H> &v) {
    v.clear_corner_history();
  }

  /** @brief Nearest pairs at both endpoints of every pixel, in scan order,
   *  as the column-time draw classifies its block corners.
   *  @return Per pixel: previous pair, then current pair. */
  template <int W, int H>
  static std::vector<std::array<uint16_t, 4>>
  endpoint_pairs(const Voronoi<W, H> &v) {
    const size_t n = v.sites_buffer.size();
    std::vector<math::Vector> current(n);
    for (size_t i = 0; i < n; ++i)
      current[i] = v.sites_buffer[i].pos;
    ScratchScope scope(scratch_arena_a);
    const KDTree now(scratch_arena_a, std::span<const math::Vector>(current));
    auto near = v.first_nearby();
    std::vector<std::array<uint16_t, 4>> out;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out.push_back(
            v.endpoint_pairs(now, math::pixel_to_vector<W, H>(x, y), near));
    return out;
  }

  /** @brief Spin of one frame at a speed, by the effect's mapping.
   *  @param speed Speed slider value.
   *  @return Radians per frame. */
  static float spin(float speed) {
    return logf(speed + 1.0f) * VO::SITE_SPIN_RADIANS;
  }
};

/**
 * @brief Verifies Voronoi seeds every spin axis through random_vector().
 */
inline void test_voronoi_axes_use_uniform_sampler() {
  using WB = VoronoiWhiteBox;
  reset_effect_globals();
  Voronoi<SMALL_W, SMALL_H> effect;
  effect.init();

  const size_t sites = WB::site_count(effect);
  HS_EXPECT_GT(sites, 0u);
  hs::random().seed(1337u);
  for (size_t i = 0; i < sites; ++i)
    HS_EXPECT_VEC(WB::site_axis(effect, i), math::random_vector(), 0.0f);
}

/**
 * @brief Renders production Voronoi and compares each pixel with exact nearest.
 * @tparam W,H Render resolution.
 * @param sites Site positions on the unit sphere.
 * @param max_deficit Out: worst dot(p, true nearest) - dot(p, union nearest)
 *        over the mismatched pixels (0 when every pixel matches).
 * @return Fraction of pixels whose union-of-corner-pairs nearest matches the
 *         true nearest site.
 */
template <int W, int H>
inline double voronoi_render_nearest_match(std::span<const math::Vector> sites,
                                           float &max_deficit) {
  using WB = VoronoiWhiteBox;
  reset_effect_globals();
  Voronoi<W, H> effect;
  effect.init();
  WB::set_sites(effect, sites);
  effect.draw_frame();
  effect.advance_display();

  long matched = 0;
  max_deficit = 0.0f;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      float exact = -2.0f;
      for (size_t i = 0; i < WB::site_count(effect); ++i)
        exact = std::max(exact, math::dot(p, WB::site_position(effect, i)));

      const Pixel rendered = effect.get_pixel(x, y);
      const size_t rendered_site = rendered.r - 1u;
      const bool encoded = rendered.r > 0 && rendered.g == 0 &&
                           rendered.b == 0 &&
                           rendered_site < WB::site_count(effect);
      const float deficit =
          encoded
              ? exact - math::dot(p, WB::site_position(effect, rendered_site))
              : 4.0f;
      if (deficit <= 0.0f)
        ++matched;
      else
        max_deficit = std::max(max_deficit, deficit);
    }
  }
  return static_cast<double>(matched) / (static_cast<double>(W) * H);
}

/**
 * @brief Pins Voronoi's block-candidate-union coverage across the adaptive
 *        block regime at both render resolutions: the union of a block's four
 *        corner pairs contains the true nearest site at every pixel in the
 *        low-density octahedral case, and at >= 99.9% of pixels (with a
 *        sub-visibility dot deficit on the rest) at the MAX_SITES Fibonacci
 *        spread.
 * @details The short-canvas resolution, where the adaptive block floors at 1,
 *          is checked for exact coverage.
 */
inline void test_voronoi_union_candidates_cover_nearest() {
  static_assert(VoronoiWhiteBox::coherence_block_min<SMALL_W, SMALL_H>() == 1);
  static_assert(VoronoiWhiteBox::coherence_block_min<DEFAULT_W, DEFAULT_H>() ==
                4);
  float deficit = 0.0f;

  const math::Vector octahedral[] = {
      math::Vector(1, 0, 0),  math::Vector(-1, 0, 0), math::Vector(0, 1, 0),
      math::Vector(0, -1, 0), math::Vector(0, 0, 1),  math::Vector(0, 0, -1),
  };
  const size_t octa_count = sizeof(octahedral) / sizeof(octahedral[0]);
  const std::span<const math::Vector> sparse(octahedral, octa_count);
  const double octa_match =
      voronoi_render_nearest_match<DEFAULT_W, DEFAULT_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match, 1.0);
  const double octa_match_dev =
      voronoi_render_nearest_match<SMALL_W, SMALL_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match_dev, 1.0);

  // Dense regime: seed MAX_SITES on a Fibonacci sphere exactly as
  // Voronoi::seed_sites places them, so the adaptive block floors at
  // COHERENCE_BLOCK_MIN.
  constexpr int N = VoronoiWhiteBox::MAX_SITES;
  static math::Vector fib[N];
  for (int i = 0; i < N; ++i)
    fib[i] = math::fib_spiral(N, /*eps=*/0.5f, i);
  const std::span<const math::Vector> dense(fib, N);
  const double fib_match =
      voronoi_render_nearest_match<DEFAULT_W, DEFAULT_H>(dense, deficit);
  HS_EXPECT_GE(fib_match, 0.999);
  HS_EXPECT_LE(deficit, 0.005f);

  const double fib_match_dev =
      voronoi_render_nearest_match<SMALL_W, SMALL_H>(dense, deficit);
  HS_EXPECT_EQ(fib_match_dev, 1.0);
  HS_EXPECT_EQ(deficit, 0.0f);
}

/**
 * @brief Requires a Voronoi segment band to shade every pixel exactly as the
 *        full-canvas render does, on column-time frames.
 * @details The coarse-coherence grid decides per block which sites reach a
 *          pixel's candidate union, so its phase must not follow the clip
 *          origin. The clipped edges cut through blocks, and bands cross both
 *          arm boundaries.
 */
inline void test_voronoi_segment_render_matches_full_frame() {
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;

  auto render = [](int x0, int x1, int y0, int y1) {
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    effect.set_margin(3);
    effect.set_clip(y0, y1, x0, x1);
    for (int frame = 0; frame < 3; ++frame) {
      effect.draw_frame();
      effect.advance_display();
    }
    std::vector<Pixel> band;
    band.reserve(static_cast<size_t>(x1 - x0) * (y1 - y0));
    for (int y = y0; y < y1; ++y)
      for (int x = x0; x < x1; ++x)
        band.push_back(effect.get_pixel(x, y));
    return band;
  };

  const std::vector<Pixel> full = render(0, W, 0, H);

  struct Band {
    int x0, x1, y0, y1;
  };
  const Band bands[] = {{0, 100, 0, H},   {100, W, 0, H},
                        {0, W, 50, H},    {37, 205, 11, 93},
                        {140, 150, 0, H}, {W / 2 - 1, W / 2 + 1, 20, 30}};

  for (const Band &b : bands) {
    HS_CONTEXT("band", b.x0, b.y0);
    size_t lit = 0;
    const std::vector<Pixel> banded = render(b.x0, b.x1, b.y0, b.y1);
    HS_EXPECT_EQ(banded.size(),
                 static_cast<size_t>(b.x1 - b.x0) * (b.y1 - b.y0));
    size_t different = 0;
    size_t i = 0;
    for (int y = b.y0; y < b.y1; ++y)
      for (int x = b.x0; x < b.x1; ++x, ++i) {
        const Pixel &reference = full[static_cast<size_t>(y) * W + x];
        if (banded[i] != reference)
          ++different;
        if (!is_black(reference))
          ++lit;
      }
    if (different)
      std::printf("  VORONOI SEAM band x[%d,%d) y[%d,%d): %zu of %zu pixels "
                  "differ from the full-canvas render\n",
                  b.x0, b.x1, b.y0, b.y1, different, banded.size());
    HS_EXPECT_EQ(different, static_cast<size_t>(0));
    HS_EXPECT_GT(lit, size_t{0});
  }
}

/**
 * @brief Renders three Voronoi frames with column-time sampling disabled.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param sites Site count.
 * @param border Border thickness.
 * @return Pin folded over the frames.
 */
template <int W, int H>
inline PixelGolden voronoi_instantaneous_render(int sites, float border) {
  reset_effect_globals();
  Voronoi<W, H> effect;
  effect.init();
  HS_EXPECT_EQ(effect.updateParameter("Num Sites", static_cast<float>(sites)),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("Border Thick", border),
               ParamSetResult::APPLIED);
  Spacetime::set_flat_surface(true);
  PixelGolden pin;
  for (int frame = 0; frame < 3; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    fold_pixel_golden<W, H>(pin, effect);
  }
  Spacetime::set_flat_surface(false);
  return pin;
}

/** @brief Voronoi renders as recorded with column-time sampling off. */
inline void test_voronoi_instantaneous_render_pins() {
  const PixelGolden pins[] = {
      voronoi_instantaneous_render<DEFAULT_W, DEFAULT_H>(200, 0.0f),
      voronoi_instantaneous_render<DEFAULT_W, DEFAULT_H>(400, 0.02f),
      voronoi_instantaneous_render<SMALL_W, SMALL_H>(200, 0.02f),
  };
  static constexpr RecordedGolden RECORDED[] = {
      {{0x48c42b8111d71c22ull,
        {{3224501316ull, 1023501533ull, 145258570ull},
         158823194083ull,
         630753351407ull}},
       {0x3e54c2a6ca7db594ull,
        {{3224501786ull, 1023501922ull, 145258571ull},
         158823193975ull,
         630753433049ull}}},
      {{0x3569d4d4a598f681ull,
        {{2496742876ull, 790651062ull, 113499247ull},
         122790894831ull,
         486978999414ull}},
       {0x0ad1754efb5d0dffull,
        {{2496756535ull, 790656005ull, 113495457ull},
         122791077284ull,
         486980241206ull}}},
      {{0xe2cc79e93c41179aull,
        {{130302290ull, 42439639ull, 5904153ull}, 799440922ull, 8436662516ull}},
       {0x5cffd0c32c547be7ull,
        {{130302275ull, 42439633ull, 5904155ull},
         799440797ull,
         8436661516ull}}},
  };
  for (size_t i = 0; i < std::size(RECORDED); ++i) {
    HS_CONTEXT("scenario", static_cast<long long>(i));
    expect_pixel_golden("Voronoi", pins[i], RECORDED[i]);
  }
}

/** @brief Brute-force agreement of a render with the sites at column time. */
struct VoronoiTimeMatch {
  long matched = 0;         ///< Pixels drawn with a nearest site.
  long total = 0;           ///< Pixels compared.
  long moved = 0;           ///< Pixels whose nearest site at column time
                            ///< differs from that at the frame's end.
  float max_deficit = 0.0f; ///< Worst dot shortfall of a mismatched pixel.
};

/**
 * @brief Compares rendered sites with an exhaustive search at each pixel's
 *        column time.
 * @tparam W,H Render resolution.
 * @param effect Effect after draw_frame() and advance_display().
 * @param from Site positions at the frame's start.
 * @param axes Site spin axes.
 * @param spin The frame's spin, radians.
 * @param x0,x1,y0,y1 Compared pixel window.
 * @return Match statistics; a pixel matches within 1e-5 of the best dot.
 */
template <int W, int H>
inline VoronoiTimeMatch voronoi_time_match(const Voronoi<W, H> &effect,
                                           std::span<const math::Vector> from,
                                           std::span<const math::Vector> axes,
                                           float spin, int x0, int x1, int y0,
                                           int y1) {
  VoronoiTimeMatch out;
  std::vector<math::Vector> at(from.size()), end(from.size());
  for (size_t i = 0; i < from.size(); ++i)
    end[i] = math::rotate(from[i], math::make_rotation(axes[i], spin));
  for (int x = x0; x < x1; ++x) {
    const float FRACTION = static_cast<float>(x % (W / 2)) / (W / 2);
    for (size_t i = 0; i < from.size(); ++i)
      at[i] =
          math::rotate(from[i], math::make_rotation(axes[i], FRACTION * spin));
    for (int y = y0; y < y1; ++y) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      size_t best = 0;
      size_t best_end = 0;
      for (size_t i = 1; i < from.size(); ++i) {
        if (math::dot(p, at[i]) > math::dot(p, at[best]))
          best = i;
        if (math::dot(p, end[i]) > math::dot(p, end[best_end]))
          best_end = i;
      }
      const Pixel rendered = effect.get_pixel(x, y);
      const size_t site = rendered.r - 1u;
      const bool ENCODED = rendered.r > 0 && rendered.g == 0 &&
                           rendered.b == 0 && site < from.size();
      const float DEFICIT =
          ENCODED ? math::dot(p, at[best]) - math::dot(p, at[site]) : 4.0f;
      ++out.total;
      out.moved += best != best_end;
      if (DEFICIT <= 1e-5f)
        ++out.matched;
      else
        out.max_deficit = std::max(out.max_deficit, DEFICIT);
    }
  }
  return out;
}

/** @brief Site positions of an effect. */
template <int W, int H>
inline std::vector<math::Vector> voronoi_positions(const Voronoi<W, H> &v) {
  std::vector<math::Vector> out(VoronoiWhiteBox::site_count(v));
  for (size_t i = 0; i < out.size(); ++i)
    out[i] = VoronoiWhiteBox::site_position(v, i);
  return out;
}

/** @brief Unit spin axes drawn from a fixed generator. */
inline std::vector<math::Vector> voronoi_test_axes(size_t count) {
  hs::Pcg32 rng(20261018u);
  std::vector<math::Vector> axes(count);
  for (auto &axis : axes) {
    const float x = rand_uniform(rng, -1.0f, 1.0f);
    const float y = rand_uniform(rng, -1.0f, 1.0f);
    const float z = rand_uniform(rng, -1.0f, 1.0f);
    axis = math::Vector(x, y, z).normalized();
  }
  return axes;
}

/**
 * @brief Renders moving sites and checks each column against an exhaustive
 *        search at its time, across the seam into the next frame.
 * @tparam W,H Render resolution.
 * @param sites Initial site positions.
 * @param minimum Least fraction of pixels that must match.
 * @param max_deficit Largest dot shortfall allowed on a mismatch.
 */
template <int W, int H>
inline void check_voronoi_column_time(std::span<const math::Vector> sites,
                                      double minimum, float max_deficit) {
  using WB = VoronoiWhiteBox;
  constexpr float SPEED = 100.0f;
  reset_effect_globals();
  Voronoi<W, H> effect;
  effect.init();
  const auto axes = voronoi_test_axes(sites.size());
  WB::set_moving_sites(effect, sites, axes, SPEED);
  effect.draw_frame();
  effect.advance_display();
  const auto first = voronoi_positions(effect);
  const auto primed =
      voronoi_time_match<W, H>(effect, first, axes, 0.0f, 0, W, 0, H);
  HS_EXPECT_GE(static_cast<double>(primed.matched) / primed.total, minimum);

  auto from = first;
  for (int frame = 0; frame < 2; ++frame) {
    HS_CONTEXT("frame", frame);
    effect.draw_frame();
    effect.advance_display();
    const auto match = voronoi_time_match<W, H>(effect, from, axes,
                                                WB::spin(SPEED), 0, W, 0, H);
    HS_EXPECT_GT(match.moved, 0);
    HS_EXPECT_GE(static_cast<double>(match.matched) / match.total, minimum);
    HS_EXPECT_LE(match.max_deficit, max_deficit);
    const auto end = voronoi_positions(effect);
    for (size_t i = 0; i < from.size(); ++i)
      HS_EXPECT_VEC(
          math::rotate(from[i], math::make_rotation(axes[i], WB::spin(SPEED))),
          end[i], 2e-6f);
    from = end;
  }
}

/** @brief Voronoi shades each column from its sites at the column's time. */
inline void test_voronoi_column_time_matches_brute_force() {
  const math::Vector octahedral[] = {
      math::Vector(1, 0, 0),  math::Vector(-1, 0, 0), math::Vector(0, 1, 0),
      math::Vector(0, -1, 0), math::Vector(0, 0, 1),  math::Vector(0, 0, -1),
  };
  check_voronoi_column_time<DEFAULT_W, DEFAULT_H>(octahedral, 1.0, 0.0f);
  check_voronoi_column_time<SMALL_W, SMALL_H>(octahedral, 1.0, 0.0f);
  constexpr int N = VoronoiWhiteBox::MAX_SITES;
  static math::Vector fib[N];
  for (int i = 0; i < N; ++i)
    fib[i] = math::fib_spiral(N, /*eps=*/0.5f, i);
  check_voronoi_column_time<DEFAULT_W, DEFAULT_H>(fib, 0.999, 0.005f);
  check_voronoi_column_time<SMALL_W, SMALL_H>(fib, 1.0, 0.0f);
}

/**
 * @brief Corner classification finds each endpoint's nearest pair, as an
 *        exhaustive search ranks them.
 */
inline void test_voronoi_endpoint_pairs_match_brute_force() {
  using WB = VoronoiWhiteBox;
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  for (const int count : {1, 2, 3, 400}) {
    HS_CONTEXT("sites", count);
    std::vector<math::Vector> sites(count);
    for (int i = 0; i < count; ++i)
      sites[i] = math::fib_spiral(count, /*eps=*/0.5f, i);
    const auto axes = voronoi_test_axes(sites.size());
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    WB::set_moving_sites(effect, std::span<const math::Vector>(sites), axes,
                         100.0f);
    effect.draw_frame();
    effect.advance_display();
    const auto from = voronoi_positions(effect);
    effect.draw_frame();
    effect.advance_display();
    const auto to = voronoi_positions(effect);
    const auto pairs = WB::endpoint_pairs(effect);
    const auto ranked = [&](const std::vector<math::Vector> &at,
                            const math::Vector &p) {
      std::vector<uint16_t> order(at.size());
      for (size_t i = 0; i < order.size(); ++i)
        order[i] = static_cast<uint16_t>(i);
      std::sort(order.begin(), order.end(), [&](uint16_t a, uint16_t b) {
        const float DA = math::distance_squared(at[a], p);
        const float DB = math::distance_squared(at[b], p);
        return DA != DB ? DA < DB : a < b;
      });
      order.resize(2, order[0]);
      return order;
    };
    long changed = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const math::Vector p = math::pixel_to_vector<W, H>(x, y);
        const auto before = ranked(from, p);
        const auto after = ranked(to, p);
        const auto &got = pairs[y * W + x];
        changed += before != after;
        HS_EXPECT_EQ(got[0], before[0]);
        HS_EXPECT_EQ(got[1], before[1]);
        HS_EXPECT_EQ(got[2], after[0]);
        HS_EXPECT_EQ(got[3], after[1]);
      }
    if (count > 2)
      HS_EXPECT_GT(changed, 0);
  }
}

/** @brief Disabled sampling draws the current sites while history advances. */
inline void test_voronoi_temporal_toggle() {
  using WB = VoronoiWhiteBox;
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;
  constexpr float SPEED = 100.0f;
  const math::Vector octahedral[] = {
      math::Vector(1, 0, 0),  math::Vector(-1, 0, 0), math::Vector(0, 1, 0),
      math::Vector(0, -1, 0), math::Vector(0, 0, 1),  math::Vector(0, 0, -1),
  };
  const auto axes = voronoi_test_axes(std::size(octahedral));
  reset_effect_globals();
  Voronoi<W, H> effect;
  effect.init();
  WB::set_moving_sites(effect, std::span<const math::Vector>(octahedral), axes,
                       SPEED);
  effect.draw_frame();
  effect.advance_display();

  Spacetime::set_flat_surface(true);
  effect.draw_frame();
  effect.advance_display();
  const auto held = voronoi_positions(effect);
  const auto instantaneous =
      voronoi_time_match<W, H>(effect, held, axes, 0.0f, 0, W, 0, H);
  HS_EXPECT_EQ(instantaneous.matched, instantaneous.total);

  Spacetime::set_flat_surface(false);
  effect.draw_frame();
  effect.advance_display();
  const auto resumed =
      voronoi_time_match<W, H>(effect, held, axes, WB::spin(SPEED), 0, W, 0, H);
  HS_EXPECT_EQ(resumed.matched, resumed.total);
  HS_EXPECT_GT(resumed.moved, 0);
}

/**
 * @brief Renders the frame that reseeds Voronoi.
 * @param enabled Column-time sampling during that frame.
 * @return The displayed frame.
 */
inline std::vector<Pixel> voronoi_reseed_frame(bool enabled) {
  reset_effect_globals();
  Voronoi<DEFAULT_W, DEFAULT_H> effect;
  effect.init();
  for (int frame = 0; frame < 2; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.updateParameter("Num Sites", 12.0f),
               ParamSetResult::APPLIED);
  Spacetime::set_flat_surface(!enabled);
  effect.draw_frame();
  effect.advance_display();
  Spacetime::set_flat_surface(false);
  std::vector<Pixel> out;
  capture_frame<DEFAULT_W, DEFAULT_H>(effect, out);
  return out;
}

/** @brief A reseed draws its first frame at the new sites without motion. */
inline void test_voronoi_reseed_resets_history() {
  const auto temporal = voronoi_reseed_frame(true);
  const auto instantaneous = voronoi_reseed_frame(false);
  HS_EXPECT_TRUE(compare_buffers(instantaneous, temporal).exact());
  HS_EXPECT_GT(count_lit(temporal), size_t{0});
}

/** @brief Retained corner pairs preserve pixels across clip and history changes. */
inline void test_voronoi_retained_corners_match_search() {
  constexpr int W = DEFAULT_W, H = DEFAULT_H;
  const auto render = [](bool retain) {
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    std::vector<Pixel> frames;
    for (int frame = 0; frame < 20; ++frame) {
      if (frame == 4)
        effect.set_clip(17, 103, 40, 151);
      if (frame == 6)
        Spacetime::set_flat_surface(true);
      if (frame == 7)
        Spacetime::set_flat_surface(false);
      if (frame == 9)
        HS_EXPECT_EQ(effect.updateParameter("Num Sites", 83.0f),
                     ParamSetResult::APPLIED);
      if (frame >= 12) {
        const int X0 = (frame % 2) * (W / 2);
        effect.set_clip(17, 103, X0, X0 + W / 2);
      }
      if (!retain)
        VoronoiWhiteBox::forget_corners(effect);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.push_back(effect.get_pixel(x, y));
    }
    return frames;
  };
  const auto searched = render(false);
  const auto retained = render(true);
  HS_EXPECT_TRUE(compare_buffers(searched, retained).exact());
}
