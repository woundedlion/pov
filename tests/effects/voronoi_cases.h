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
 * @details The dense cases floor the adaptive block at COHERENCE_BLOCK_MIN. A
 *          block edge that outruns the cell pixel extent straddles whole cells
 *          and collapses the match fraction, so the short-canvas resolution —
 *          where MAX_SITES cells are sub-pixel and the floor drops to 1 — is
 *          checked for exact coverage.
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
 *        full-canvas render does.
 * @details The coarse-coherence grid decides per block which sites reach a
 *          pixel's candidate union, so its phase must not follow the clip
 *          origin. The clipped edges cut through blocks.
 */
inline void test_voronoi_segment_render_matches_full_frame() {
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;

  auto render = [](int x0, int x1, int y0, int y1,
                   const std::vector<Pixel> *reference = nullptr) {
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    effect.set_margin(3);
    effect.set_clip(y0, y1, x0, x1);
    effect.draw_frame();
    effect.advance_display();
    if (reference)
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          if (effect.clip().contains_y(y) && effect.clip().contains_x(x))
            HS_EXPECT_EQ(effect.get_pixel(x, y),
                         (*reference)[static_cast<size_t>(y) * W + x]);
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
  const Band bands[] = {
      {0, 100, 0, H}, {100, W, 0, H}, {0, W, 50, H}, {37, 205, 11, 93}};

  for (const Band &b : bands) {
    HS_CONTEXT("band", b.x0, b.y0);
    size_t lit = 0;
    const std::vector<Pixel> banded = render(b.x0, b.x1, b.y0, b.y1, &full);
    HS_EXPECT_EQ(banded.size(),
                 static_cast<size_t>(b.x1 - b.x0) * (b.y1 - b.y0));
    size_t different = 0;
    size_t i = 0;
    for (int y = b.y0; y < b.y1; ++y)
      for (int x = b.x0; x < b.x1; ++x, ++i) {
        const Pixel &reference = full[static_cast<size_t>(y) * W + x];
        if (banded[i] != reference)
          ++different;
        if (reference.r | reference.g | reference.b)
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
