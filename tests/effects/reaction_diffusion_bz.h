/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Belousov-Zhabotinsky reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for BZReactionDiffusion's private fixed-point and
 *        physics internals.
 * @details The lattice is independent of <W,H>.
 */
struct BZWhiteBox {
  using BZ = BZReactionDiffusion<SMALL_W, SMALL_H>;
  using Grid = Scan::Shader::SsaaGrid<SMALL_W, SMALL_H>;
  static constexpr int N = BZ::RD_N;

  static Pixel palette_color(const BZ &bz, float t) {
    const auto &color = t == 0.0f   ? bz.color_a
                        : t == 0.5f ? bz.color_b
                                    : bz.color_c;
    return Pixel(static_cast<uint16_t>(color.r), static_cast<uint16_t>(color.g),
                 static_cast<uint16_t>(color.b));
  }
  static uint16_t to_q16(float v) { return BZ::to_q16(v); }
  static float from_q16(uint16_t v) { return BZ::from_q16(v); }
  static void set_params(BZ &bz, float alpha, float D, float dt) {
    bz.params.alpha = alpha;
    bz.params.D = D;
    bz.params.dt = dt;
  }
  static uint16_t advance_species(const BZ &bz, float conc, float predator,
                                  float laplacian) {
    return bz.advance_species(conc, predator, laplacian);
  }
  static void perturb(BZ &bz, uint16_t *nA, uint16_t *nB, uint16_t *nC) {
    bz.perturb_state(nA, nB, nC);
  }
  static int num_perturbations() { return BZ::NUM_PERTURBATIONS; }
  static constexpr int perturb_amount() { return BZ::PERTURB_AMOUNT; }
  static void step(BZ &bz, uint16_t *sA, uint16_t *sB, uint16_t *sC) {
    std::vector<float> fA(N), fB(N), fC(N);
    bz.step_physics(sA, sB, sC, fA.data(), fB.data(), fC.data());
  }
  static void set_state(BZ &bz, uint16_t *a, uint16_t *b, uint16_t *c) {
    bz.state.A = a;
    bz.state.B = b;
    bz.state.C = c;
  }
  static int nearest_node(const math::Vector &v, const math::Vector *nodes) {
    int nearest = 0;
    float nearest_d2 = BZ::dist2(v, nodes[0]);
    for (int i = 1; i < N; ++i) {
      float d2 = BZ::dist2(v, nodes[i]);
      if (d2 < nearest_d2) {
        nearest = i;
        nearest_d2 = d2;
      }
    }
    return nearest;
  }
  struct CenterError {
    int pixels = 0;
    int refined_seeds = 0;
    int center_mismatches = 0;
  };

  /**
   * @brief Sweeps a <W,H> frame of cubemap-LUT seeds and compares the certified
   *        render-center early-out against the unconditional argmin walk.
   * @details Seeds come from seed_face_lut, so a quantized seed that is not
   *          the nearest node reaches the walk; refined_seeds counts those.
   */
  template <int W, int H> static CenterError render_center_error(BZ &bz) {
    ScratchScope guard(scratch_arena_a);
    auto lattice = bz.orient_lattice();
    math::Vector *world_nodes = lattice.get();
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    CenterError error;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        Fragment frag;
        frag.pos = math::pixel_to_vector<W, H>(x, y);
        bz.seed_face_lut(frag);
        int seed = static_cast<int>(frag.v0);
        int reference = BZ::refine_center(frag.pos, world_nodes, seed);
        if (reference != seed)
          ++error.refined_seeds;
        if (BZ::refine_render_center(frag.pos, world_nodes, seed) != reference)
          ++error.center_mismatches;
        ++error.pixels;
      }
    }
    return error;
  }
  static Pixel shade(const BZ &bz, int seed, const math::Vector &center,
                     const math::Vector *nodes, const Grid &grid, int x,
                     const Color4 &ca, const Color4 &cb, const Color4 &cc) {
    BZ::FloatRgb fa(ca.color), fb(cb.color), fc(cc.color);
    return bz.shade_pixel(seed, center, nodes, grid, x, fa, fb, fc);
  }
  static Pixel reference_shade(const BZ &bz, int seed,
                               const math::Vector &center_rv,
                               const math::Vector *world_nodes,
                               const Grid &grid, int x, const Color4 &ca,
                               const Color4 &cb, const Color4 &cc) {
    int center = BZ::refine_center(center_rv, world_nodes, seed);
    math::Vector spos[BZ::RD_K + 1];
    uint16_t sa[BZ::RD_K + 1], sb[BZ::RD_K + 1], sc[BZ::RD_K + 1];
    spos[0] = world_nodes[center];
    sa[0] = bz.state.A[center];
    sb[0] = bz.state.B[center];
    sc[0] = bz.state.C[center];
    int k = 1;
    BZ::for_each_neighbor(center, [&](int ni) {
      spos[k] = world_nodes[ni];
      sa[k] = bz.state.A[ni];
      sb[k] = bz.state.B[ni];
      sc[k] = bz.state.C[ni];
      ++k;
    });

    constexpr float INV_SAMPLES = 1.0f / 4.0f;
    Pixel accum(0, 0, 0);
    for (int i = 0; i < 4; ++i) {
      math::Vector v = grid.at(x, i);
      float tw = 0, wa = 0, wb = 0, wc = 0;
      for (int j = 0; j < BZ::RD_K + 1; ++j)
        BZ::with_biweight_weight(BZ::dist2(v, spos[j]), [&](float w) {
          wa += sa[j] * w;
          wb += sb[j] * w;
          wc += sc[j] * w;
          tw += w;
        });
      if (tw <= BZ::KERNEL_MIN_TOTAL_WEIGHT)
        continue;
      float inv = BZ::Q16_INV / tw;
      float a = wa * inv, b = wb * inv, c = wc * inv;
      float total = a + b + c;
      if (total < BZ::SPECIES_EMPTY_EPS)
        continue;
      float color_inv = 1.0f / total;
      Pixel color(
          static_cast<uint16_t>(
              (ca.color.r * a + cb.color.r * b + cc.color.r * c) * color_inv +
              0.5f),
          static_cast<uint16_t>(
              (ca.color.g * a + cb.color.g * b + cc.color.g * c) * color_inv +
              0.5f),
          static_cast<uint16_t>(
              (ca.color.b * a + cb.color.b * b + cc.color.b * c) * color_inv +
              0.5f));
      accum += color * (hs::clamp(total, 0.0f, 1.0f) * INV_SAMPLES);
    }
    return accum;
  }
};

/** @brief Pins the authored BZ species colors. */
inline void test_bz_authored_palette() {
  reset_effect_globals();
  BZWhiteBox::BZ bz;
  bz.init();
  const auto expect_color = [&](float position, const Pixel &expected) {
    const Pixel actual = BZWhiteBox::palette_color(bz, position);
    HS_EXPECT_NEAR(static_cast<float>(actual.r), static_cast<float>(expected.r),
                   4.0f);
    HS_EXPECT_NEAR(static_cast<float>(actual.g), static_cast<float>(expected.g),
                   4.0f);
    HS_EXPECT_NEAR(static_cast<float>(actual.b), static_cast<float>(expected.b),
                   4.0f);
  };
  expect_color(0.0f, Pixel(36844, 10770, 3));
  expect_color(0.5f, Pixel(0, 8112, 5753));
  expect_color(1.0f, Pixel(2059, 0, 9668));
}

/**
 * @brief Verifies the Q16 fixed-point round-trip and the +0.5 rounding/clamp
 *        boundaries.
 * @details to_q16(from_q16(v)) must be the identity over every representable
 *          value, and to_q16 must clamp out-of-range floats and round to nearest
 *          (so 1.0 tops out at 65535 with no overflow).
 */
inline void test_bz_q16_roundtrip() {
  HS_EXPECT_EQ(BZWhiteBox::to_q16(0.0f), (uint16_t)0);
  HS_EXPECT_EQ(BZWhiteBox::to_q16(1.0f), (uint16_t)65535);
  HS_EXPECT_EQ(BZWhiteBox::to_q16(2.0f), (uint16_t)65535); // clamp high
  HS_EXPECT_EQ(BZWhiteBox::to_q16(-0.5f), (uint16_t)0);    // clamp low
  HS_EXPECT_NEAR(BZWhiteBox::from_q16(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(BZWhiteBox::from_q16(65535), 1.0f, 1e-9f);
  int bad = 0;
  for (int v = 0; v <= 65535; ++v)
    if (BZWhiteBox::to_q16(BZWhiteBox::from_q16((uint16_t)v)) != (uint16_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies the state resolution carries Diff's minimum step at default Speed.
 * @details A cold node beside one saturated neighbour at Diff's 0.001 floor
 *          and Speed 0.35 receives D·lap·dt = 3.5e-4 of full scale; the stored
 *          sample must move.
 */
inline void test_bz_min_diffusion_step_survives_quantization() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.001f, /*dt*/ 0.35f);
  // Cold node, one saturated neighbour out of RD_K: lap = 1, reaction term 0.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, /*lap*/ 1.0f),
               0);
  // The same gradient one LSB above the floor still resolves.
  const float lsb = BZWhiteBox::from_q16(1);
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, lsb, 0.0f, /*lap*/ 1.0f),
               (int)BZWhiteBox::to_q16(lsb));
}

/**
 * @brief Verifies advance_species has the right reaction/diffusion signs and
 *        that the Q16 clamp backstop holds even past the Euler stability bound.
 * @details advance_species is the single-species core of the BZ update:
 *          conc + (D·laplacian + conc·(1 − conc − α·predator))·dt, mapped through
 *          to_q16. Extreme over- and under-shoots clamp to the [0, 65535]
 *          rails rather than wrapping.
 */
inline void test_bz_advance_species_signs_and_clamp() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 0.35f);

  // Empty rest cell: no diffusion, no reaction -> stays 0.
  HS_EXPECT_EQ((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, 0.0f), 0);

  // Diffusion from higher neighbors lifts an empty cell above 0.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, /*lap*/ 6.0f),
               0);

  // Logistic growth: a half-filled, predator-free cell grows above its start.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.5f, 0.0f, 0.0f),
               (int)BZWhiteBox::to_q16(0.5f));

  // Predation drives a saturated cell negative; it must clamp to 0, not wrap.
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1.0f, 0.0f), 0);

  // Backstop past the Euler bound: a huge positive laplacian clamps to 65535, a
  // huge predation clamps to 0 — to_q16 keeps every written state in range.
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, 0.0f, /*lap*/ 1000.0f), 65535);
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1000.0f, 0.0f),
      0);
}

/**
 * @brief Verifies perturb_state nudges nodes by a fixed Q16 amount and saturates
 *        at the 65535 rail without wrapping.
 * @details At dt = 1 the nudge is the full PERTURB_AMOUNT. A saturated field
 *          stays saturated; on a zero field every touched entry is a multiple
 *          of the nudge step, and at least one entry is touched.
 */
inline void test_bz_perturb_state_saturates_and_nudges() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 1.0f);
  // Saturation / no-wrap: all rails stay at the rail.
  {
    std::vector<uint16_t> a(BZWhiteBox::N, 65535), b(BZWhiteBox::N, 65535),
        c(BZWhiteBox::N, 65535);
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int wrapped = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      if (a[i] != 65535 || b[i] != 65535 || c[i] != 65535)
        ++wrapped;
    HS_EXPECT_EQ(wrapped, 0);
  }
  // Zero field: touched entries are small positive multiples of the step.
  {
    const int step = BZWhiteBox::perturb_amount();
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int touched = 0, malformed = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      for (uint16_t v : {a[i], b[i], c[i]}) {
        if (v == 0)
          continue;
        ++touched;
        if (v % step != 0) // accumulations stay multiples of the nudge step
          ++malformed;
      }
    HS_EXPECT_GT(touched, 0);
    HS_EXPECT_EQ(malformed, 0);
  }
}

/**
 * @brief Pins perturb_state's per-frame draw count on the shared RNG stream.
 * @details perturb_state advances hs::random() by exactly 2*NUM_PERTURBATIONS
 *          draws (idx + species per nudge) at both ends of the Speed slider;
 *          downstream stream positions depend on that count.
 */
inline void test_bz_perturb_state_draw_count_pinned() {
  const int expected_draws = 2 * BZWhiteBox::num_perturbations();

  constexpr uint64_t SEED = 1337u;
  for (float dt : {1.0f, 0.0f}) {
    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, dt);
    hs::random().seed(SEED);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());

    // A private generator from the same seed, advanced by the contracted count,
    // must now be at the same stream position as the global generator.
    hs::Pcg32 ref(SEED);
    for (int i = 0; i < expected_draws; ++i)
      (void)ref();
    HS_EXPECT_EQ(hs::random()(), ref());
    // Off-by-one in either direction would have landed at a different output.
    HS_EXPECT_EQ(hs::random()(), ref());
  }
  reset_effect_globals();
}

/**
 * @brief Verifies the stochastic nudge scales with the Speed slider and reaches
 *        zero where the integrator freezes.
 * @details At dt = 0 nothing may move; at a mid-slider dt the touched entries
 * must be multiples of the scaled step, between zero and the full-rate step.
 */
inline void test_bz_perturb_scales_with_timestep() {
  constexpr int PASSES = 64;
  {
    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 0.0f);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    for (int p = 0; p < PASSES; ++p)
      BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int moved = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      if (a[i] != 0 || b[i] != 0 || c[i] != 0)
        ++moved;
    HS_EXPECT_EQ(moved, 0);
  }
  {
    constexpr float DT = 0.35f;
    constexpr int full = BZWhiteBox::perturb_amount();
    constexpr int step = static_cast<int>(full * DT);
    static_assert(step > 0 && full > step);

    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, DT);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    for (int p = 0; p < PASSES; ++p)
      BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int touched = 0, malformed = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      for (uint16_t v : {a[i], b[i], c[i]}) {
        if (v == 0)
          continue;
        ++touched;
        if (v % step != 0)
          ++malformed;
      }
    HS_EXPECT_GT(touched, 0);
    HS_EXPECT_EQ(malformed, 0);
  }
}

/**
 * @brief Verifies one fused physics substep diffuses a seeded species into its
 *        neighborhood with the right sign, in place.
 * @details A single step from one fully-seeded lattice node: A must diffuse
 *          into at least one empty neighbor and the seed must stay lit. The step
 *          writes the new generation over the state it read.
 */
inline void test_bz_substep_diffuses() {
  std::vector<uint16_t> sA(BZWhiteBox::N, 0), sB(BZWhiteBox::N, 0),
      sC(BZWhiteBox::N, 0);
  const int seed = 4000;
  sA[seed] = 65535;

  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, 3.0f, 0.05f, 0.35f);
  BZWhiteBox::step(bz, sA.data(), sB.data(), sC.data());

  HS_EXPECT_GT((int)sA[seed], 0); // the seed decays but does not vanish/wrap
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    HS_EXPECT_GE(sA[nb], BZWhiteBox::advance_species(bz, 0, 0, 1.0f));
    spread += sA[nb] > 0;
  }
  HS_EXPECT_GT(spread, 0); // A diffused into at least one empty neighbor
}

/**
 * @brief Pins the optimized BZ raster against its scalar sampling contract.
 * @details reference_shade normalizes each SSAA sample to concentrations,
 *          blends the palette into a
 *          uint16 Pixel, premultiplies by that sample's coverage, and adds the
 *          result into a uint16 accumulator. shade_pixel fuses the same algebra
 *          into float species coefficients and quantizes once at the end, so
 *          the two agree only up to the reference's intermediate rounding:
 *          four palette-blend roundings weighted by coverages that sum to at
 *          most 1 (<= 0.5 LSB), four premultiply roundings (<= 2.0 LSB), and
 *          the single final rounding both paths pay (<= 0.5 LSB). Any channel
 *          past that 3 LSB envelope is a formula difference, not round-off.
 */
inline void test_bz_raster_matches_reference() {
  using WhiteBox = BZWhiteBox;
  constexpr int MAX_ROUNDING_DRIFT = 3;
  auto drift = [](uint16_t p, uint16_t q) { return p > q ? p - q : q - p; };
  std::vector<math::Vector> nodes(WhiteBox::N);
  std::vector<uint16_t> a(WhiteBox::N), b(WhiteBox::N), c(WhiteBox::N);
  for (int i = 0; i < WhiteBox::N; ++i) {
    nodes[i] = ReactionGraph::node(i);
    a[i] = static_cast<uint16_t>((i * 4051u + 123u) & 0xffffu);
    b[i] = static_cast<uint16_t>((i * 7919u + 4567u) & 0xffffu);
    c[i] = static_cast<uint16_t>((i * 10429u + 8901u) & 0xffffu);
  }

  WhiteBox::BZ bz;
  WhiteBox::set_state(bz, a.data(), b.data(), c.data());
  if (!math::TrigLUT<SMALL_W, SMALL_H>::initialized)
    math::TrigLUT<SMALL_W, SMALL_H>::init();
  WhiteBox::Grid grid;
  const Color4 ca(Pixel(61123, 913, 17771), 1.0f);
  const Color4 cb(Pixel(2819, 59731, 1207), 1.0f);
  const Color4 cc(Pixel(4513, 7781, 62927), 1.0f);

  auto compare = [&] {
    int mismatches = 0;
    for (int y = 0; y < SMALL_H; y += 2) {
      grid.set_row(y);
      for (int x = 0; x < SMALL_W; x += 3) {
        math::Vector center = math::pixel_to_vector<SMALL_W, SMALL_H>(x, y);
        int seed = WhiteBox::nearest_node(center, nodes.data());
        Pixel expected = WhiteBox::reference_shade(
            bz, seed, center, nodes.data(), grid, x, ca, cb, cc);
        Pixel actual = WhiteBox::shade(bz, seed, center, nodes.data(), grid, x,
                                       ca, cb, cc);
        if (drift(actual.r, expected.r) > MAX_ROUNDING_DRIFT ||
            drift(actual.g, expected.g) > MAX_ROUNDING_DRIFT ||
            drift(actual.b, expected.b) > MAX_ROUNDING_DRIFT)
          ++mismatches;
      }
    }
    return mismatches;
  };

  int mismatches = compare();
  for (int i = 0; i < WhiteBox::N; ++i) {
    a[i] >>= 6;
    b[i] >>= 6;
    c[i] >>= 6;
  }
  mismatches += compare();
  HS_EXPECT_EQ(mismatches, 0);
}

/**
 * @brief Requires BZ's certified render-center early-out to agree with the
 *        unconditional argmin on production-resolution, LUT-seeded pixels.
 * @details refine_render_center returns the seed unwalked inside
 *          BULK_CERTIFIED_D2 / POLE_CERTIFIED_D2. Only a cubemap-LUT seed can be
 *          a non-nearest node. A few frames of the orientation random walk move
 *          the lattice off its initial alignment first; refined_seeds must be
 *          nonzero.
 */
inline void test_bz_render_center_matches_reference() {
  hs_test::reset_globals();
  {
    BZWhiteBox::BZ bz;
    bz.init();
    for (int frame = 0; frame < 4; ++frame) {
      bz.draw_frame();
      bz.advance_display();
    }
    auto error = BZWhiteBox::render_center_error<DEFAULT_W, DEFAULT_H>(bz);
    std::printf("BZ render centers: pixels=%d refined_seeds=%d mismatches=%d\n",
                error.pixels, error.refined_seeds, error.center_mismatches);
    HS_EXPECT_EQ(error.pixels, DEFAULT_W * DEFAULT_H);
    HS_EXPECT_GT(error.refined_seeds, 0);
    HS_EXPECT_EQ(error.center_mismatches, 0);
  }
  reset_effect_globals();
}
