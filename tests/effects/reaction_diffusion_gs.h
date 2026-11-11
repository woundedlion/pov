/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Gray-Scott reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/** @brief Verifies source palette recipes are opaque before the RGB bake. */
inline void test_gs_palette_is_opaque() {
  hs_test::reset_globals();
  for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed) {
    auto palette = GSWhiteBox::source_palette();
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
      HS_EXPECT_EQ(
          palette.get(static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1))
              .alpha,
          1.0f);
  }
}

/** @brief Verifies each reseeded reaction receives a freshly generated palette. */
inline void test_gs_reseed_generates_palette() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  std::array<Pixel, BakedPalette::LUT_SIZE> before;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
    before[i] = GSWhiteBox::palette_sample(gs, static_cast<float>(i) /
                                                   (BakedPalette::LUT_SIZE - 1))
                    .color;

  const size_t arena_offset = persistent_arena.get_offset();
  GSWhiteBox::start_reaction(gs);

  int changed = 0;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const Color4 after = GSWhiteBox::palette_sample(
        gs, static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1));
    changed += after.color != before[i];
  }
  HS_EXPECT_GT(changed, 0);
  HS_EXPECT_EQ(persistent_arena.get_offset(), arena_offset);
}

/** @brief Shares cube projection without changing noise on faces, edges or corners. */
inline void test_gs_shared_cube_projection() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (float phase : {0.0f, 0.123f, 0.9999f}) {
    GSWhiteBox::set_noise_phase(gs, phase);
    GSWhiteBox::advance_color_noise(gs);
    for (int face = 0; face < Pattern::SphereLutView::FACE_COUNT; ++face)
      for (float u : {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f})
        for (float v : {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) {
          const auto direction = Pattern::sphere_face_direction(face, u, v);
          const auto projection = ReactionGraph::CubemapLUT::project(direction);
          if (fabsf(u) < 1.0f && fabsf(v) < 1.0f) {
            HS_EXPECT_EQ(projection.face, face);
            HS_EXPECT_NEAR(projection.u, face < 2 ? -u : u, 2e-7f);
            HS_EXPECT_NEAR(projection.v, face >= 2 && face < 4 ? -v : v, 2e-7f);
          }
          HS_EXPECT_EQ(GSWhiteBox::projected_color_noise(gs, direction),
                       GSWhiteBox::color_noise(gs, direction));
        }
    for (const auto &direction : ReactionGraph::node_positions)
      HS_EXPECT_EQ(GSWhiteBox::projected_color_noise(gs, direction),
                   GSWhiteBox::color_noise(gs, direction));
  }
}

/** @brief
 * Seed palettes share the same spatial noise modifiers for hue and shimmer.
 */
inline void test_gs_shared_noise_palette_modifiers() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  const size_t offset = persistent_arena.get_offset();
  int rotated = 0;
  for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed) {
    for (float t : {0.0f, 0.37f, 1.0f})
      for (float shift : {-0.3f, 0.0f, 0.25f})
        for (float lift : {0.0f, 0.5f, 1.0f})
          HS_EXPECT_EQ(GSWhiteBox::modified_palette(gs, seed, t, shift, lift),
                       GSWhiteBox::wrapped_palette(gs, seed, t, shift, lift));
    Pixel base = GSWhiteBox::palette_sample(gs, 0.5f, seed).color;
    HS_EXPECT_EQ(GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.0f, 0.0f),
                 base);
    Pixel hue = GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.25f, 0.0f);
    Pixel expected = hue_rotate_lut_gamut(Color4(base, 1.0f), 0.25f).color;
    HS_EXPECT_EQ(hue, expected);
    rotated += hue != base;
    Pixel shimmer = GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.0f, 0.5f);
    const LinRGB before = pixel_to_linrgb(base);
    const LinRGB after = pixel_to_linrgb(shimmer);
    HS_EXPECT_GT(linear_rgb_to_oklab(after.r, after.g, after.b).L,
                 linear_rgb_to_oklab(before.r, before.g, before.b).L);
  }
  HS_EXPECT_EQ(rotated, GSWhiteBox::SEEDS);

  GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  const float still = GSWhiteBox::color_noise(gs, math::X_AXIS);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_EQ(GSWhiteBox::noise_phase(gs), 0.0f);
  HS_EXPECT_EQ(GSWhiteBox::color_noise(gs, math::X_AXIS), still);
  GSWhiteBox::set_color_params(gs, 0.001f, 2.0f, 0.35f, 0.25f);
  for (int i = 0; i < 40; ++i)
    GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_GT(fabsf(GSWhiteBox::color_noise(gs, math::X_AXIS) - still), 1e-3f);
  const float moving = GSWhiteBox::color_noise(gs, math::X_AXIS);
  GSWhiteBox::set_color_params(gs, 0.0f, 4.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_GT(fabsf(GSWhiteBox::color_noise(gs, math::X_AXIS) - moving),
               1e-3f);
  GSWhiteBox::set_noise_phase(gs, 0.0005f);
  GSWhiteBox::set_color_params(gs, -0.001f, 4.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_NEAR(GSWhiteBox::noise_phase(gs), 0.9995f, 1e-6f);
  HS_EXPECT(!GSWhiteBox::reaction_edited(gs),
            "color controls must not restart the reaction");
  HS_EXPECT_EQ(persistent_arena.get_offset(), offset);
}

/** @brief The baked color noise samples match the direct reference. */
inline void test_gs_noise_bake_matches_reference() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (float scale : {1.0f / 64.0f, 4.941984f, 8.0f}) {
    for (float phase : {0.0f, 1.0f / 65536.0f, 1.0f / 6.0f, 0.5f, 0.999f}) {
      HS_EXPECT_TRUE(
          GSWhiteBox::noise_bake_matches_reference(gs, scale, phase));
    }
  }
}

/** @brief
 * Pigment travels between seeded regions while preserving distinct seed
 * palettes.
 */
inline void test_gs_seed_palettes_and_pigment_transport() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int seed = 1; seed < GSWhiteBox::SEEDS; ++seed) {
    bool distinct = false;
    for (int j = 0; j < 9; ++j)
      distinct |= GSWhiteBox::palette_sample(gs, j / 8.0f, seed).color !=
                  GSWhiteBox::palette_sample(gs, j / 8.0f, seed - 1).color;
    HS_EXPECT(distinct, "each seed must receive an independent palette");
  }

  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f);
  constexpr int NODE = GSWhiteBox::N / 2;
  const int neighbor = ReactionGraph::neighbors[NODE][0];
  b[NODE] = 0.5f;
  GSWhiteBox::set_pigment(gs, NODE, 7);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, neighbor, 7), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 7), 1.0f, 1e-6f);

  b[neighbor] = 0.5f;
  GSWhiteBox::set_pigment(gs, neighbor, 12);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  float mix = GSWhiteBox::pigment_weight(gs, NODE, 12);
  HS_EXPECT_GT(mix, 0.0f);
  HS_EXPECT_LT(mix, 1.0f);
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 7) + mix, 1.0f, 1e-6f);
  Pixel first = GSWhiteBox::palette_sample(gs, 0.5f, 7).color;
  Pixel second = GSWhiteBox::palette_sample(gs, 0.5f, 12).color;
  Pixel blended = GSWhiteBox::pigment_color(gs, NODE, 0.5f);
  HS_EXPECT_NEAR(blended.r, first.r * (1.0f - mix) + second.r * mix, 2.0f);
  HS_EXPECT_NEAR(blended.g, first.g * (1.0f - mix) + second.g * mix, 2.0f);
  HS_EXPECT_NEAR(blended.b, first.b * (1.0f - mix) + second.b * mix, 2.0f);

  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.0f, 2.5f);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 12), mix, 1e-6f);

  const size_t offset = persistent_arena.get_offset();
  GSWhiteBox::start_reaction(gs);
  HS_EXPECT_EQ(persistent_arena.get_offset(), offset);
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    if (GSWhiteBox::b_field(gs)[i] != 0)
      continue;
    HS_EXPECT_EQ(GSWhiteBox::pigment_weight(gs, i, 0), 1.0f);
  }
}

/** @brief The unit-vector dot kernel agrees with squared-distance evaluation. */
inline void test_gs_dot_kernel_matches_squared_distance() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  HS_EXPECT_GT(GSWhiteBox::check_render_kernel(gs), 0);
  hs_test::run_frames(gs, 40);
  HS_EXPECT_GT(GSWhiteBox::check_render_kernel(gs), 0);
}

/** @brief The direct Gray-Scott draw matches the grid reference. */
inline void test_gs_direct_draw_matches_grid() {
  hs_test::reset_globals();
  GSReactionDiffusion<DEFAULT_W, DEFAULT_H> gs;
  gs.init();
  HS_EXPECT_GT(GSWhiteBox::check_direct_draw(gs), size_t{0});
  hs_test::run_frames(gs, 40);
  HS_EXPECT_GT(GSWhiteBox::check_direct_draw(gs), size_t{0});
}

/** @brief Sparse pigment transport agrees with the dense reference. */
inline void test_gs_sparse_pigment_matches_dense() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
  std::vector<uint16_t> expected(GSWhiteBox::N);
  const float PARAMS[][5] = {{0.04f, 0.06f, 0.02f, 0.01f, 2.5f},
                             {0.1f, 0.1f, 0.05f, 0.05f, 3.0f},
                             {0.04f, 0.06f, 0.02f, 0.0f, 0.1f},
                             {0.04f, 0.06f, 0.02f, -0.0f, 0.1f}};
  for (const auto &params : PARAMS) {
    GSWhiteBox::set_params(gs, params[0], params[1], params[2], params[3],
                           params[4]);
    for (int mode = 0; mode < 5; ++mode) {
      for (int i = 0; i < GSWhiteBox::N; ++i) {
        uint32_t h = (static_cast<uint32_t>(i) + 173u) * 2654435761u;
        a[i] = static_cast<float>(h & 65535u) / 65535.0f;
        b[i] = mode == 0 ? 0.0f
                         : static_cast<float>((h >> 16) & 65535u) / 65535.0f;
        int first = mode == 1 ? 17 : static_cast<int>(h % GSWhiteBox::SEEDS);
        int second =
            mode == 1 ? 0 : static_cast<int>((h >> 8) % GSWhiteBox::SEEDS);
        int mix = mode <= 2   ? GSWhiteBox::Pigment::FULL_MASS
                  : mode == 3 ? 32
                              : static_cast<int>(h & 63u);
        GSWhiteBox::pigments(gs)[i] =
            GSWhiteBox::Pigment::pack(first, second, mix);
      }
      for (int step = 0; step < GSWhiteBox::STEPS_PER_FRAME; ++step) {
        GSWhiteBox::pigment_reference(gs, a.data(), b.data(), expected.data());
        GSWhiteBox::step_pigment(gs, a.data(), b.data());
        for (int i = 0; i < GSWhiteBox::N; ++i)
          HS_EXPECT_EQ(GSWhiteBox::pigments(gs)[i], expected[i]);
        GSWhiteBox::pigment_reference<GSWhiteBox::STEPS_PER_FRAME>(
            gs, a.data(), b.data(), expected.data());
        GSWhiteBox::step_pigment<GSWhiteBox::STEPS_PER_FRAME>(gs, a.data(),
                                                              b.data());
        for (int i = 0; i < GSWhiteBox::N; ++i)
          HS_EXPECT_EQ(GSWhiteBox::pigments(gs)[i], expected[i]);
      }
    }
  }
}

/** @brief Verifies the pigment codec's bit layout and mass quantization. */
inline void test_gs_pigment_codec_layout() {
  using Pigment = GSWhiteBox::Pigment;
  HS_EXPECT_EQ(Pigment::pack(7, 12, 40), 7u | (12u << 5) | (40u << 10));
  const auto UNPACKED = Pigment::unpack(0xffffu);
  HS_EXPECT_EQ(UNPACKED.first, 31);
  HS_EXPECT_EQ(UNPACKED.second, 31);
  HS_EXPECT_EQ(UNPACKED.mass, 63);
  HS_EXPECT_EQ(Pigment::unpack(Pigment::pack(3, 29, 17)).second, 29);
  HS_EXPECT_EQ(Pigment::pack(5, 9, 21) & Pigment::SOLID_KEY_MASK,
               Pigment::pack(5, 0, 21));
  HS_EXPECT_EQ(Pigment::first_mass(Pigment::FULL_MASS), 1.0f);
  HS_EXPECT_EQ(Pigment::first_mass(0), 0.0f);
  HS_EXPECT_EQ(Pigment::quantize_mass(1.0f, 1.0f), 63);
  HS_EXPECT_EQ(Pigment::quantize_mass(1.0f, 2.0f), 32);
  HS_EXPECT_EQ(Pigment::quantize_mass(0.0f, 2.0f), 0);
}

/** @brief Compares both cull rings, as byte flags and as the slab shader's
 *         bit sets, to the original directed adjacency table. */
inline void test_gs_hot_flags_match_directed_graph() {
  constexpr int N = GSWhiteBox::N;
  constexpr uint8_t GUARD = 0xa5;
  std::vector<uint16_t> b(N);
  std::vector<uint8_t> hot1(N + 2), hot2(N + 2), ref1(N), ref2(N);
  std::vector<uint32_t> bits1(GSWhiteBox::HOT_WORDS),
      bits2(GSWhiteBox::HOT_WORDS);
  const uint16_t thresholds[] = {
      0, 1, GSWhiteBox::to_q16(GSWhiteBox::B_CULL_THRESHOLD), 65535};
  for (uint16_t threshold : thresholds) {
    const uint16_t below = threshold == 0 ? 0 : threshold - 1;
    const uint16_t above = threshold == 65535 ? 65535 : threshold + 1;
    for (int pattern = 0; pattern < 9; ++pattern) {
      for (int i = 0; i < N; ++i) {
        const uint32_t hash = static_cast<uint32_t>(i) * 2654435761u;
        switch (pattern) {
        case 0:
          b[i] = 0;
          break;
        case 1:
          b[i] = 65535;
          break;
        case 2:
          b[i] = below;
          break;
        case 3:
          b[i] = threshold;
          break;
        case 4:
          b[i] = hash % 257 == 0 ? above : below;
          break;
        case 5:
          b[i] = hash % 7 != 0 ? above : below;
          break;
        case 6:
          b[i] = hash % 3 == 0 ? below : hash % 3 == 1 ? threshold : above;
          break;
        case 7:
          b[i] = i == 0 ? above : below;
          break;
        case 8:
          b[i] = i == N - 1 ? above : below;
          break;
        }
      }
      for (int i = 0; i < N; ++i) {
        ref1[i] = b[i] >= threshold;
        for (int k = 0; k < ReactionGraph::RD_K; ++k)
          ref1[i] |= b[ReactionGraph::neighbors[i][k]] >= threshold;
      }
      for (int i = 0; i < N; ++i) {
        ref2[i] = ref1[i];
        for (int k = 0; k < ReactionGraph::RD_K; ++k)
          ref2[i] |= ref1[ReactionGraph::neighbors[i][k]];
      }
      std::fill(hot1.begin(), hot1.end(), GUARD);
      std::fill(hot2.begin(), hot2.end(), GUARD);
      GSWhiteBox::fill_hot_flags(b.data(), hot1.data() + 1, hot2.data() + 1, N,
                                 threshold);
      HS_EXPECT(std::equal(ref1.begin(), ref1.end(), hot1.begin() + 1),
                "GS one-ring flags match directed graph");
      HS_EXPECT(std::equal(ref2.begin(), ref2.end(), hot2.begin() + 1),
                "GS two-ring flags match directed graph");
      HS_EXPECT_EQ(hot1.front(), GUARD);
      HS_EXPECT_EQ(hot1.back(), GUARD);
      HS_EXPECT_EQ(hot2.front(), GUARD);
      HS_EXPECT_EQ(hot2.back(), GUARD);

      std::fill(bits1.begin(), bits1.end(), 0xa5a5a5a5u);
      std::fill(bits2.begin(), bits2.end(), 0xa5a5a5a5u);
      GSWhiteBox::fill_hot_bits(b.data(), threshold, bits1.data(),
                                bits2.data());
      int mismatched = 0;
      for (int i = 0; i < N; ++i)
        mismatched += GSWhiteBox::hot_bit(bits1.data(), i) != (ref1[i] != 0) ||
                      GSWhiteBox::hot_bit(bits2.data(), i) != (ref2[i] != 0);
      HS_EXPECT_EQ(mismatched, 0);
      for (int w = 0; w < GSWhiteBox::HOT_GUARD_WORDS; ++w)
        for (const auto *bits : {&bits1, &bits2}) {
          HS_EXPECT_EQ((*bits)[w], 0u);
          HS_EXPECT_EQ((*bits)[bits->size() - 1 - w], 0u);
        }
    }
  }
  int reach = 0;
  for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r)
    for (int k = 0; k < ReactionGraph::RD_K; ++k)
      reach = std::max(reach, std::abs(static_cast<int>(
                                  ReactionGraph::neighbor_runs[r].delta[k])));
  HS_EXPECT_LE(reach + 32, 32 * GSWhiteBox::HOT_GUARD_WORDS);
}

/**
 * @brief Re-measures the lattice and requires refine_render_center's early-out
 *        certificates to bound the true minimum node spacing.
 * @details The early-out returns the seed unwalked whenever the query lies
 *          within sqrt(safe_d2) of it, which is the nearest node only while
 *          safe_d2 stays at or under a quarter of the seed's true
 *          nearest-neighbor distance squared. node()'s y is strictly
 *          decreasing in the index, so the outward scan stops once |dy|
 *          reaches the running best.
 */
inline void test_gs_render_certificates_bound_lattice() {
  const int N = GSWhiteBox::N;
  std::vector<math::Vector> lattice(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i)
    lattice[static_cast<size_t>(i)] = ReactionGraph::node(i);

  const float BULK_FLOOR =
      GSWhiteBox::BULK_MIN_SPACING_FRAC * ReactionGraph::D_AVG;
  float bulk_min_d2 = 4.0f; // antipodal chord squared: the loosest possible
  float pole_min_d2 = 4.0f;
  int tight_band = 0; // band width the sub-bulk-spacing nodes actually need
  for (int i = 0; i < N; ++i) {
    const math::Vector &p = lattice[static_cast<size_t>(i)];
    float best = 4.0f;
    for (int step = -1; step <= 1; step += 2) {
      for (int j = i + step; j >= 0 && j < N; j += step) {
        const math::Vector &q = lattice[static_cast<size_t>(j)];
        float dy = p.y - q.y;
        if (dy * dy >= best)
          break;
        float dx = p.x - q.x, dz = p.z - q.z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best)
          best = d2;
      }
    }
    if (i < GSWhiteBox::POLE_BAND || i >= N - GSWhiteBox::POLE_BAND) {
      pole_min_d2 = std::min(pole_min_d2, best);
    } else {
      bulk_min_d2 = std::min(bulk_min_d2, best);
    }
    if (best < BULK_FLOOR * BULK_FLOOR)
      tight_band = std::max(tight_band, std::min(i, N - 1 - i) + 1);
  }

  const float bulk_frac = std::sqrt(bulk_min_d2) / ReactionGraph::D_AVG;
  std::printf("  [info] GS lattice: bulk min spacing %.6f*D_AVG, bulk d2/4 "
              "%.9f, pole d2/4 %.9f, tight band %d\n",
              bulk_frac, 0.25f * bulk_min_d2, 0.25f * pole_min_d2, tight_band);

  HS_EXPECT_LE(GSWhiteBox::BULK_MIN_SPACING_FRAC, bulk_frac);
  HS_EXPECT_LE(GSWhiteBox::BULK_CERTIFIED_D2, 0.25f * bulk_min_d2);
  HS_EXPECT_LE(GSWhiteBox::POLE_CERTIFIED_D2, 0.25f * pole_min_d2);
  // Every node packed tighter than the bulk minimum must fall inside the band
  // that gets the smaller certificate.
  HS_EXPECT_LE(tight_band, GSWhiteBox::POLE_BAND);
}

/**
 * @brief Bounds the shared stencil against independent SSAA refinement.
 * @details Compares production-resolution pixels after 4, 16, and 40 rendered
 * frames. Coverage changes stay within 2% of lit pixels, errors past 6.25%
 * within 3%, and changes beyond two RGB16 quantization units within half.
 * Max channel error stays within 23000 and mean within 3/512 full scale.
 */
inline void test_gs_shared_stencil_error_is_bounded() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr int PROBE_FRAMES[] = {4, 16, 40};
  int next_probe = 0;
  for (int frame = 1; frame <= PROBE_FRAMES[2]; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != PROBE_FRAMES[next_probe])
      continue;
    auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, {.center_pigment = true, .full_kernel = true});
    std::printf("GS shared stencil frame=%d: different=%d above_rounding=%d "
                "lit=%d coverage=%d "
                "hard=%d center_changes=%d/%d center_mismatches=%d max=%d "
                "total=%llu\n",
                frame, error.different, error.above_rounding, error.lit,
                error.coverage, error.hard, error.stencil_changes,
                error.stencil_samples, error.center_mismatches,
                error.max_channel,
                static_cast<unsigned long long>(error.total_channel));
    HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_EQ(error.center_mismatches, 0);
    HS_EXPECT(error.above_rounding * 2 <= error.lit,
              "shared stencil changed over half of lit pixels beyond rounding");
    HS_EXPECT(error.coverage * 50 <= error.lit,
              "shared stencil moved over 2% of lit coverage");
    HS_EXPECT(error.hard * 100 <= error.lit * 3,
              "shared stencil put over 3% of lit pixels past 6.25% error");
    HS_EXPECT_LE(error.max_channel, 23000);
    HS_EXPECT(error.total_channel <=
                  static_cast<uint64_t>(error.lit) * 3u * 384u,
              "shared stencil mean channel error exceeds 3/512 full scale");
    if (next_probe < 2)
      ++next_probe;
  }
}

/** @brief
 * Partial color bakes update the requested rows and retain other rows.
 */
inline void test_gs_partial_color_palette_rows() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr uint16_t ALL_ROWS = (1u << 15) - 1;
  // A non-empty run of adjacent rows.
  const auto contiguous = [](uint16_t mask) {
    if (mask == 0)
      return false;
    while (!(mask & 1u))
      mask >>= 1;
    return (mask & (mask + 1u)) == 0;
  };
  HS_EXPECT_EQ(GSWhiteBox::color_palette_rows(gs), ALL_ROWS);
  const size_t offset = persistent_arena.get_offset();
  for (int generation = 0; generation < 3; ++generation) {
    if (generation == 1)
      GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, -0.2f, 0.3f);
    else
      GSWhiteBox::start_reaction(gs);
    GSWhiteBox::refresh_color_palettes(gs);
    uint16_t rows = GSWhiteBox::color_palette_rows(gs);
    HS_EXPECT_TRUE(rows != ALL_ROWS);
    HS_EXPECT_TRUE((rows >> (GSWhiteBox::COLOR_NOISE_STEPS / 2)) & 1u);
    HS_EXPECT_TRUE(contiguous(rows));
    HS_EXPECT_TRUE(GSWhiteBox::exact_color_sample(gs, -1.0f));
    HS_EXPECT_TRUE(GSWhiteBox::exact_color_sample(gs, 1.0f));
    float hue = generation == 0 ? 0.2f : -0.2f;
    float shimmer = generation == 0 ? 0.4f : 0.3f;
    for (float noise : {-1.0f + 12.0f / 14.0f, 0.0f}) {
      HS_EXPECT_FALSE(GSWhiteBox::exact_color_sample(gs, noise));
      for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed)
        for (float value : {0.0f, 1.0f})
          HS_EXPECT_EQ(
              GSWhiteBox::staged_palette(gs, seed, value, noise),
              GSWhiteBox::modified_palette(gs, seed, value, noise * hue,
                                           std::max(noise, 0.0f) * shimmer));
    }
    for (int frame = 0; frame < 7; ++frame) {
      for (int n = 0; n < GSWhiteBox::COLOR_NOISE_STEPS; ++n) {
        const float noise =
            -1.0f + 2.0f * n / (GSWhiteBox::COLOR_NOISE_STEPS - 1);
        const bool exact = GSWhiteBox::exact_color_sample(gs, noise);
        const float KNOT = 5.0f / (GSWhiteBox::COLOR_VALUE_STEPS - 1);
        for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed)
          for (float value : {0.0f, exact ? 0.37f : KNOT, 1.0f})
            HS_EXPECT_EQ(
                GSWhiteBox::staged_palette(gs, seed, value, noise),
                GSWhiteBox::modified_palette(gs, seed, value, noise * hue,
                                             std::max(noise, 0.0f) * shimmer));
      }
      GSWhiteBox::refresh_color_palettes(gs);
      const uint16_t grown = GSWhiteBox::color_palette_rows(gs);
      HS_EXPECT_EQ(grown & rows, rows);
      HS_EXPECT_TRUE(contiguous(grown));
      rows = grown;
    }
    HS_EXPECT_EQ(rows, ALL_ROWS);
    HS_EXPECT_EQ(persistent_arena.get_offset(), offset);
  }
}

/**
 * @brief Verifies a reseed at unchanged hue and shimmer rebakes the shimmer
 *        cache from the new palettes.
 */
inline void test_gs_reseed_invalidates_shimmer_cache() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  const auto controls = GSWhiteBox::color_params(gs);
  const float hue = controls[2], shimmer = controls[3];
  HS_EXPECT_FALSE(GSWhiteBox::exact_color_sample(gs, 0.0f));
  const Pixel before = GSWhiteBox::staged_palette(gs, 0, 1.0f, 0.0f);
  GSWhiteBox::start_reaction(gs);
  HS_EXPECT(GSWhiteBox::modified_palette(gs, 0, 1.0f, 0.0f, 0.0f) != before,
            "reseed kept the first seed palette");
  GSWhiteBox::refresh_color_palettes(gs, true);
  int stale = 0;
  for (int n = 0; n < GSWhiteBox::COLOR_NOISE_STEPS; ++n) {
    const float noise = -1.0f + 2.0f * n / (GSWhiteBox::COLOR_NOISE_STEPS - 1);
    HS_EXPECT_FALSE(GSWhiteBox::exact_color_sample(gs, noise));
    for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed)
      for (float value : {0.0f, 1.0f})
        stale += GSWhiteBox::staged_palette(gs, seed, value, noise) !=
                 GSWhiteBox::modified_palette(gs, seed, value, noise * hue,
                                              std::max(noise, 0.0f) * shimmer);
  }
  HS_EXPECT_EQ(stale, 0);
}

/** @brief Checks certified support against the fixed-profile SSAA offsets. */
inline void test_gs_support_certificate_matches_display_geometry() {
  HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<96, 20>()),
               GSWhiteBox::RENDER_MIN_WEIGHT);
  HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<288, 16>()),
               GSWhiteBox::RENDER_MIN_WEIGHT);
  HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<288, 144>()),
               GSWhiteBox::RENDER_MIN_WEIGHT);
  HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<4096, 16>()),
               GSWhiteBox::RENDER_MIN_WEIGHT);
}

/** @brief Preserves four-sample coverage while bounding reciprocal estimation error. */
inline void test_gs_signed_coverage_and_concentration() {
  const uint16_t floor = GSWhiteBox::to_q16(GSWhiteBox::B_CULL_THRESHOLD);
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int frame = 1; frame <= 650; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 40 && frame != 200 && frame != 400 &&
        frame != 450 && frame != 470 && frame != 600 && frame != 650)
      continue;
    const auto error =
        GSWhiteBox::concentration_shader_error<DEFAULT_W, DEFAULT_H>(gs);
    HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_EQ(error.coverage, 0);
    HS_EXPECT_LE(error.max_channel, 6000);
    HS_EXPECT_LE(error.total_channel,
                 static_cast<uint64_t>(error.lit) * 3u * 125u);
  }
  int pattern_lit = 0, fallback_lit = 0;
  for (int pattern = 0; pattern < 8; ++pattern) {
    for (int node = 0; node < GSWhiteBox::N; ++node) {
      const uint32_t hash = static_cast<uint32_t>(node) * 2654435761u;
      const uint16_t b = pattern < 4
                             ? static_cast<uint16_t>(floor - 2 + pattern)
                             : static_cast<uint16_t>(
                                   floor - 2 + ((hash >> (pattern - 4)) & 3));
      GSWhiteBox::set_node(gs, node, 65535, b);
    }
    const auto error =
        GSWhiteBox::concentration_shader_error<DEFAULT_W, DEFAULT_H>(gs);
    pattern_lit += error.lit;
    HS_EXPECT_EQ(error.coverage, 0);
    const auto fallback = GSWhiteBox::concentration_shader_error<16, 8>(gs);
    fallback_lit += fallback.lit;
    HS_EXPECT_EQ(fallback.max_channel, 0);
    HS_EXPECT_EQ(fallback.total_channel, 0u);
  }
  HS_EXPECT_GT(pattern_lit, 0);
  HS_EXPECT_GT(fallback_lit, 0);
}

/** @brief Bounds nearest-pigment color fidelity against procedural and aggregate references. */
inline void test_gs_nearest_pigment_shader_fidelity() {
  // Float kernel rearrangement can move one pixel across the render floor.
  constexpr int MAX_COVERAGE_DIFFERENCES = 1;
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int frame = 1; frame <= 600; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 40 && frame != 200 && frame != 400 &&
        frame != 500 && frame != 600)
      continue;
    const auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, {.fixed_stencil = true,
             .direct_color = true,
             .aggregate_color = false});
    HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_LE(error.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(error.max_channel, 22000);
    HS_EXPECT_LE(error.total_channel,
                 static_cast<uint64_t>(error.lit) * 3u * 850u);
    HS_EXPECT_LE(error.hard * 10, error.lit);
    const auto rounding = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs,
        {.fixed_stencil = true, .center_pigment = true, .full_kernel = true});
    HS_EXPECT_GT(rounding.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_LE(rounding.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(rounding.max_channel, 2048);
    HS_EXPECT_LE(rounding.total_channel,
                 static_cast<uint64_t>(rounding.lit) * 6u);
    const auto aggregate =
        GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
            gs, {.fixed_stencil = true});
    HS_EXPECT_GT(aggregate.lit, DEFAULT_W * DEFAULT_H / 100);
    std::printf(
        "GS nearest pigment frame=%d procedural_mae=%.2f procedural_max=%d "
        "aggregate_mae=%.2f aggregate_max=%d rounding_max=%d\n",
        frame, static_cast<double>(error.total_channel) / (error.lit * 3),
        error.max_channel,
        static_cast<double>(aggregate.total_channel) / (aggregate.lit * 3),
        aggregate.max_channel, rounding.max_channel);
    HS_EXPECT_LE(aggregate.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(aggregate.max_channel, 22000);
    HS_EXPECT_LE(aggregate.total_channel,
                 static_cast<uint64_t>(aggregate.lit) * 3u * 400u);
    if (frame != 200 && frame != 400 && frame != 600)
      continue;
    const auto controls = GSWhiteBox::color_params(gs);
    for (float hue : {-4.0f, -0.25f, 0.0f, 0.25f, 4.0f}) {
      for (float shimmer : {0.0f, 0.4f, 1.0f}) {
        GSWhiteBox::set_color_params(gs, controls[0], controls[1], hue,
                                     shimmer);
        GSWhiteBox::refresh_color_palettes(gs);
        const auto extreme =
            GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
                gs, {.fixed_stencil = true,
                     .direct_color = true,
                     .aggregate_color = false});
        if (extreme.max_channel > 30000)
          std::printf(
              "GS nearest pigment extreme frame=%d hue=%.2f shimmer=%.2f max=%d\n",
              frame, static_cast<double>(hue), static_cast<double>(shimmer),
              extreme.max_channel);
        HS_EXPECT_LE(extreme.coverage, MAX_COVERAGE_DIFFERENCES);
        HS_EXPECT_LE(extreme.max_channel, 35000);
        HS_EXPECT_LE(extreme.total_channel,
                     static_cast<uint64_t>(extreme.lit) * 3u * 1000u);
      }
    }
    GSWhiteBox::set_color_params(gs, controls[0], controls[1], controls[2],
                                 controls[3]);
  }
}

/** @brief Exact palette paths enforce hue and shimmer boundaries. */
inline void test_gs_exact_color_path_boundaries() {
  reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (float hue : {-4.0f, -0.3f, 0.0f, 0.3f, 4.0f}) {
    const float shimmer = hue == 0.0f ? 0.0f : 0.4f;
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, hue, shimmer);
    for (float value : {0.0f, 0.37f, 1.0f})
      for (float noise : {-1.0f, -0.25f, 0.0f, 0.8f, 1.0f})
        HS_EXPECT_EQ(
            GSWhiteBox::cached_palette(gs, 0, value, noise),
            GSWhiteBox::modified_palette(gs, 0, value, hue * noise,
                                         std::max(noise, 0.0f) * shimmer));
  }
  constexpr float SHIMMER_LIMIT = GSWhiteBox::CACHED_SHIMMER_LIMIT;
  constexpr float HUE_LIMIT = GSWhiteBox::CACHED_HUE_LIMIT;
  const float INSIDE_HUE = 0.8f * HUE_LIMIT;
  for (float shimmer :
       {0.0f, SHIMMER_LIMIT, std::nextafter(SHIMMER_LIMIT, 2.0f), 1.0f}) {
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, INSIDE_HUE, shimmer);
    HS_EXPECT_EQ(GSWhiteBox::exact_color_path(gs), shimmer > SHIMMER_LIMIT);
    if (shimmer > SHIMMER_LIMIT)
      for (float noise : {-1.0f, 0.0f, 0.31f, 1.0f})
        HS_EXPECT_EQ(
            GSWhiteBox::cached_palette(gs, 0, 0.37f, noise),
            GSWhiteBox::modified_palette(gs, 0, 0.37f, INSIDE_HUE * noise,
                                         std::max(noise, 0.0f) * shimmer));
  }
  const float past_hue = std::nextafter(HUE_LIMIT, 2.0f);
  for (float hue : {-past_hue, -HUE_LIMIT, HUE_LIMIT, past_hue}) {
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, hue, SHIMMER_LIMIT);
    HS_EXPECT_EQ(GSWhiteBox::exact_color_path(gs), fabsf(hue) > HUE_LIMIT);
  }
}

/** @brief Checks nearest-pigment shading against independent scalar SSAA samples. */
inline void test_gs_nearest_pigment_shader_matches_scalar_reference() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, 0.0f, 0.0f);
  for (int frame = 1; frame <= 40; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 16 && frame != 40)
      continue;
    auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs,
        {.fixed_stencil = true, .center_pigment = true, .full_kernel = true});
    HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_EQ(error.coverage, 0);
    HS_EXPECT_LE(error.max_channel, 2048);
    HS_EXPECT_LE(error.total_channel, static_cast<uint64_t>(error.lit) * 3u);
  }
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    uint16_t b = i < GSWhiteBox::N / 2
                     ? GSWhiteBox::to_q16(GSWhiteBox::B_CULL_THRESHOLD)
                     : 32768;
    GSWhiteBox::set_node(gs, i, 65535, b);
  }
  auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
      gs, {.fixed_stencil = true, .center_pigment = true, .full_kernel = true});
  HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
  HS_EXPECT_EQ(error.coverage, 0);
  HS_EXPECT_LE(error.max_channel, 2048);
}

/** @brief The dissolve frontier fades before clearing its underlying state. */
inline void test_gs_dissolve_frontier_fades_before_clear() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr float PHASE = 0.5f;
  constexpr uint32_t SEED = 0x7b31d2a5u;
  int cleared = -1, fading = -1;
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    float h = GSWhiteBox::dissolve_hash(i, SEED);
    if (cleared < 0 && h < PHASE)
      cleared = i;
    float fade = (h - PHASE) / GSWhiteBox::DISSOLVE_FADE_FRACTION;
    if (fading < 0 && fade > 0.25f && fade < 0.75f)
      fading = i;
  }
  HS_EXPECT(cleared >= 0 && fading >= 0,
            "dissolve probe did not find clear and fade nodes");
  if (cleared < 0 || fading < 0)
    return;
  GSWhiteBox::set_node(gs, cleared, 0, 65535);
  GSWhiteBox::set_node(gs, fading, 0, 65535);
  GSWhiteBox::convert(gs, PHASE, SEED);
  HS_EXPECT_EQ(GSWhiteBox::a_field(gs)[cleared], (uint16_t)65535);
  HS_EXPECT_EQ(GSWhiteBox::b_field(gs)[cleared], (uint16_t)0);
  HS_EXPECT(GSWhiteBox::a_field(gs)[fading] > 0 &&
                GSWhiteBox::a_field(gs)[fading] < 65535,
            "dissolve fade did not ease A toward rest");
  HS_EXPECT(GSWhiteBox::b_field(gs)[fading] > 0 &&
                GSWhiteBox::b_field(gs)[fading] < 65535,
            "dissolve fade cleared B without an intermediate value");
}

/**
 * @brief The two-frame staged reseed matches start_reaction() fields, palettes
 *        and RNG position.
 */
inline void test_gs_staged_reseed_matches_synchronous() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::check_staged_reseed(gs);
}

/**
 * @brief Verifies dissolve coverage falls below one eighth and reseeds.
 * @details On the last dissolve frame, coverage is below one eighth of the
 *          grown pattern. The closing frame restores the fresh seed pattern.
 */
inline void test_gs_dissolve_clears_and_reseeds() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  auto covered = [&]() {
    const uint16_t *b = GSWhiteBox::b_field(gs);
    int n = 0;
    for (int i = 0; i < GSWhiteBox::N; ++i)
      if (b[i] >= GSWhiteBox::to_q16(GSWhiteBox::B_CULL_THRESHOLD))
        ++n;
    return n;
  };

  // Run until the stalled field trips a dissolve; the default reaction settles
  // well inside this budget.
  const int BUDGET =
      GSWhiteBox::MIN_GROW_FRAMES + GSWhiteBox::STABLE_HOLD_FRAMES + 400;
  int f = 0;
  for (; f < BUDGET && GSWhiteBox::dissolve_frame(gs) < 0; ++f) {
    gs.draw_frame();
    gs.advance_display();
  }
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "stalled field did not start a dissolve");
  const int grown = covered();
  HS_EXPECT(grown > GSWhiteBox::N / 20,
            "default reaction grew no pattern to dissolve");

  // Sample the last frame before the dissolve window closes.
  hs_test::run_frames(gs, GSWhiteBox::DISSOLVE_FRAMES - 1);
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "dissolve ended before its window closed");
  HS_EXPECT(covered() < grown / 8,
            "dissolve left the sphere covered: cleared nodes healed");

  // The closing frame converts the remainder and reseeds.
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  HS_EXPECT(covered() > 0, "reseed left no nucleation sites");
}

/**
 * @brief Verifies editing the reaction constants starts a dissolve.
 * @details Feed/Kill/dA/dB define the reaction, so an edit clears and restarts
 *          it; Speed only sets the integration rate and must not.
 */
inline void test_gs_reaction_edit_starts_dissolve() {
  hs_test::reset_globals();
  {
    GSWhiteBox::GS gs;
    gs.init();
    gs.draw_frame();
    gs.advance_display();
    HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);

    for (float speed : {0.1f, 1.0f}) {
      GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, speed);
      gs.draw_frame();
      gs.advance_display();
      HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
    }
  }

  for (int parameter = 0; parameter < 4; ++parameter) {
    HS_CONTEXT("reaction parameter", parameter);
    hs_test::reset_globals();
    GSWhiteBox::GS edited;
    edited.init();
    edited.draw_frame();
    edited.advance_display();
    HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(edited), -1);
    float values[] = {0.04f, 0.06f, 0.02f, 0.01f};
    values[parameter] *= 0.75f;
    GSWhiteBox::set_params(edited, values[0], values[1], values[2], values[3],
                           1.0f);
    edited.draw_frame();
    edited.advance_display();
    HS_EXPECT_GE(GSWhiteBox::dissolve_frame(edited), 0);
  }
}

/** @brief Frames of synthetic mean |dB| until a dissolve starts, or -1. */
template <typename MeanDb>
inline int gs_frames_to_dissolve(float dt, MeanDb &&mean_db, int limit) {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, dt);
  for (int frame = 1; frame <= limit; ++frame) {
    GSWhiteBox::advance_transition(gs, mean_db(frame));
    if (GSWhiteBox::dissolve_frame(gs) >= 0)
      return frame;
  }
  return -1;
}

/**
 * @brief Verifies a dt-scaled sub-floor mean |dB| held for STABLE_HOLD_FRAMES
 *        after MIN_GROW_FRAMES starts a dissolve, and an above-floor frame
 *        restarts the hold.
 */
inline void test_gs_stall_hold_starts_dissolve() {
  constexpr int GROW = GSWhiteBox::MIN_GROW_FRAMES;
  constexpr int HOLD = GSWhiteBox::STABLE_HOLD_FRAMES;
  constexpr int LIMIT = GROW + 3 * HOLD;
  constexpr int SPIKE = GROW + HOLD / 2;
  constexpr float DT = 2.5f;
  const float floor = GSWhiteBox::stable_floor(DT);
  HS_EXPECT_EQ(gs_frames_to_dissolve(
                   DT, [](int) { return 0.0f; }, LIMIT),
               GROW + HOLD - 1);
  HS_EXPECT_EQ(gs_frames_to_dissolve(
                   DT, [&](int) { return 0.75f * floor; }, LIMIT),
               GROW + HOLD - 1);
  HS_EXPECT_EQ(
      gs_frames_to_dissolve(
          DT, [&](int frame) { return frame == SPIKE ? 2.0f * floor : 0.0f; },
          LIMIT),
      SPIKE + HOLD);
  HS_EXPECT_EQ(gs_frames_to_dissolve(
                   DT / 2.0f, [&](int) { return 0.75f * floor; }, LIMIT),
               -1);
}

/** @brief The preset pause flag does not hold reaction growth or reseeding. */
inline void test_gs_preset_pause_does_not_hold_lifecycle() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr int TO_DISSOLVE =
      GSWhiteBox::MIN_GROW_FRAMES + GSWhiteBox::STABLE_HOLD_FRAMES - 1;
  gs.setAnimationsPaused(true);
  for (int frame = 0; frame < TO_DISSOLVE; ++frame)
    GSWhiteBox::advance_transition(gs, 0.0f);
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), 0);
  constexpr uint32_t DISSOLVE_SEED = 12345;
  GSWhiteBox::convert(gs, 0.0f, DISSOLVE_SEED);
  for (int i = 0; i < GSWhiteBox::N; ++i)
    GSWhiteBox::set_node(gs, i, 32768, 16384);
  for (int frame = 1; frame <= 4; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), frame);
    const float phase = static_cast<float>(frame) / GSWhiteBox::DISSOLVE_FRAMES;
    int cleared = 0;
    for (int i = 0; i < GSWhiteBox::N; ++i)
      if (GSWhiteBox::dissolve_hash(i, DISSOLVE_SEED) < phase) {
        HS_EXPECT_EQ(GSWhiteBox::a_field(gs)[i], 65535);
        HS_EXPECT_EQ(GSWhiteBox::b_field(gs)[i], 0);
        ++cleared;
      }
    HS_EXPECT_GT(cleared, 0);
  }
  for (int frame = 4; frame < GSWhiteBox::DISSOLVE_FRAMES; ++frame)
    GSWhiteBox::advance_transition(gs, 0.0f);
  HS_EXPECT_TRUE(GSWhiteBox::reseed_pending(gs));
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_FALSE(GSWhiteBox::reseed_pending(gs));
  HS_EXPECT_EQ(GSWhiteBox::grow_frames(gs), 0);
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::grow_frames(gs), 1);

  gs.setAnimationsPaused(true);
  GSWhiteBox::set_params(gs, 0.035f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::advance_transition(gs, 1.0f);
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), 0);
  for (int frame = 0; frame < GSWhiteBox::DISSOLVE_FRAMES; ++frame)
    GSWhiteBox::advance_transition(gs, 1.0f);
  HS_EXPECT_TRUE(GSWhiteBox::reseed_pending(gs));
}

/**
 * @brief Verifies a dissolve ends in a staged reseed after DISSOLVE_FRAMES
 *        frames and edits made before the reseed completes start no dissolve.
 */
inline void test_gs_dissolve_absorbs_edits_and_stages_reseed() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::set_params(gs, 0.035f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::advance_transition(gs, 1.0f);
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), 0);
  GSWhiteBox::set_params(gs, 0.03f, 0.06f, 0.02f, 0.01f, 2.5f);
  for (int frame = 1; frame < GSWhiteBox::DISSOLVE_FRAMES; ++frame) {
    GSWhiteBox::advance_transition(gs, 1.0f);
    HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), frame);
    HS_EXPECT_FALSE(GSWhiteBox::reseed_pending(gs));
  }
  GSWhiteBox::advance_transition(gs, 1.0f);
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  HS_EXPECT_TRUE(GSWhiteBox::reseed_pending(gs));
  HS_EXPECT_EQ(GSWhiteBox::grow_frames(gs), 0);
  GSWhiteBox::set_params(gs, 0.045f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::finish_reseed(gs);
  HS_EXPECT_FALSE(GSWhiteBox::reseed_pending(gs));
  GSWhiteBox::advance_transition(gs, 1.0f);
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  HS_EXPECT_EQ(GSWhiteBox::grow_frames(gs), 1);
}

/**
 * @brief Verifies the homogeneous rest state (A=1, B=0) is a Gray-Scott fixed
 *        point: one substep leaves it exactly unchanged.
 * @details With B=0 the reaction term A·B² and both Laplacians vanish and the
 *          feed term feed·(1-A) is zero at A=1, so the state must not move. A
 *          stray nonzero term in the update assembly perturbs it here.
 */
inline void test_gs_rest_state_is_fixed_point() {
  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f),
      na(GSWhiteBox::N), nb(GSWhiteBox::N);
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f); // defaults
  GSWhiteBox::step_physics(gs, a.data(), b.data(), na.data(), nb.data());
  int moved = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (na[i] != 1.0f || nb[i] != 0.0f)
      ++moved;
  HS_EXPECT_EQ(moved, 0);
}

/**
 * @brief Pins the optimized float substep to the scalar equation bit-for-bit.
 */
inline void test_gs_substep_matches_scalar_reference() {
  std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
  std::vector<float> gotA(GSWhiteBox::N), gotB(GSWhiteBox::N);
  std::vector<float> refA(GSWhiteBox::N), refB(GSWhiteBox::N);
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    a[i] = static_cast<float>((i * 73) & 65535) * (1.0f / 65535.0f);
    b[i] = static_cast<float>((i * 151 + 17) & 65535) * (1.0f / 65535.0f);
  }
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.037f, 0.061f, 0.019f, 0.013f, 2.3f);
  GSWhiteBox::step_physics(gs, a.data(), b.data(), gotA.data(), gotB.data());
  GSWhiteBox::step_float_reference(gs, a.data(), b.data(), refA.data(),
                                   refB.data());
  HS_EXPECT(
      std::memcmp(gotA.data(), refA.data(), gotA.size() * sizeof(float)) == 0,
      "optimized A substep differs from scalar reference");
  HS_EXPECT(
      std::memcmp(gotB.data(), refB.data(), gotB.size() * sizeof(float)) == 0,
      "optimized B substep differs from scalar reference");
}

/** @brief Pins delayed writes to full float Jacobi across a complete frame. */
inline void test_gs_inplace_frame_matches_jacobi() {
  GSWhiteBox::validate_physics_neighbors();
  for (int i = 0; i < GSWhiteBox::N; ++i)
    for (int neighbor : ReactionGraph::neighbors[i])
      HS_EXPECT_LE(i - neighbor, GSWhiteBox::PHYSICS_NEIGHBOR_REACH);

  const float PARAMS[][5] = {{0.04f, 0.06f, 0.02f, 0.01f, 2.5f},
                             {0.1f, 0.1f, 0.05f, 0.05f, 3.0f},
                             {0.0f, 0.0f, 0.0f, 0.0f, 0.1f}};
  GSWhiteBox::GS gs;
  for (const auto &params : PARAMS) {
    GSWhiteBox::set_params(gs, params[0], params[1], params[2], params[3],
                           params[4]);
    std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
    for (int i = 0; i < GSWhiteBox::N; ++i) {
      a[i] = static_cast<float>((i * 73) & 65535) * (1.0f / 65535.0f);
      b[i] = static_cast<float>((i * 151 + 17) & 65535) * (1.0f / 65535.0f);
    }
    std::vector<float> ref_a = a, ref_b = b;
    std::vector<float> next_a(GSWhiteBox::N), next_b(GSWhiteBox::N);
    for (int step = 0; step < GSWhiteBox::STEPS_PER_FRAME; ++step) {
      GSWhiteBox::step_physics(gs, ref_a.data(), ref_b.data(), next_a.data(),
                               next_b.data());
      ref_a.swap(next_a);
      ref_b.swap(next_b);
      GSWhiteBox::step_float_inplace(gs, a.data(), b.data());
      HS_EXPECT(std::memcmp(a.data(), ref_a.data(), a.size() * sizeof(float)) ==
                    0,
                "delayed A writes differ from full Jacobi");
      HS_EXPECT(std::memcmp(b.data(), ref_b.data(), b.size() * sizeof(float)) ==
                    0,
                "delayed B writes differ from full Jacobi");
    }
  }
}

/**
 * @brief Verifies substep reaction/diffusion signs and the float [0, 1] clamp.
 */
inline void test_gs_substep_signs_and_clamp() {
  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f),
      next_a(GSWhiteBox::N), next_b(GSWhiteBox::N);
  constexpr int SEED = 4000;
  b[SEED] = 1.0f;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::step_physics(gs, a.data(), b.data(), next_a.data(),
                           next_b.data());

  HS_EXPECT_EQ(next_a[SEED], 0.0f);
  HS_EXPECT_EQ(next_b[SEED], 1.0f);
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[SEED][k];
    if (next_b[nb] > 0.0f)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0);
}

/** @brief Runs the boundedness fixture through 256 physics substeps. */
inline std::pair<std::vector<float>, std::vector<float>>
bounded_gs_fields(float feed, float kill, float diffusion_a, float diffusion_b,
                  float timestep) {
  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f),
      next_a(GSWhiteBox::N), next_b(GSWhiteBox::N);
  for (int seed : {500, 2500, 4500, 6500})
    b[seed] = 1.0f;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, feed, kill, diffusion_a, diffusion_b, timestep);
  for (int step = 0; step < 256; ++step) {
    GSWhiteBox::step_physics(gs, a.data(), b.data(), next_a.data(),
                             next_b.data());
    a.swap(next_a);
    b.swap(next_b);
  }
  return {std::move(a), std::move(b)};
}

/**
 * @brief Verifies the explicit-Euler integrator does not diverge over many
 *        substeps at a high-diffusion stable setting.
 * @details The stability product dt·D·|λ|max stays within the Euler bound, so
 *          fewer than 5% of nodes may sit at the upper rail. Whether B persists or
 *          decays is regime-dependent and not asserted.
 */
inline void test_gs_evolution_stays_bounded() {
  const auto [a, b] = bounded_gs_fields(0.04f, 0.06f, 0.03f, 0.03f, 3.0f);
  int saturated = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (b[i] >= 1.0f)
      ++saturated;
  HS_EXPECT_LT(saturated, GSWhiteBox::N / 20);
}

/**
 * @brief Verifies the substep stays finite and inside [0, 1] at the joint
 *        feed/k corner of the slider box, past the Euler stability bound.
 * @details At top Speed/diffusion both the diffusion and reaction terms exceed
 * the Euler bound; float fields must stay finite and in [0, 1] before Q16
 * output conversion.
 */
inline void test_gs_reaction_corner_stays_bounded() {
  const auto [a, b] = bounded_gs_fields(0.1f, 0.1f, 0.05f, 0.05f, 3.0f);
  int escaped = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (!std::isfinite(a[i]) || !std::isfinite(b[i]) || a[i] < 0.0f ||
        a[i] > 1.0f || b[i] < 0.0f || b[i] > 1.0f)
      ++escaped;
  HS_EXPECT_EQ(escaped, 0);
}

/**
 * @brief With temporal sampling disabled, GS renders lit, and a
 *        banded clip matches the full render.
 */
inline void test_gs_untimed_band_matches_full_render() {
  expect_untimed_rd_render_properties<GSReactionDiffusion>("GS");
}

/**
 * @brief Each column of a GS frame shaded through the lattice basis matches the
 *        lattice oriented at that column's sweep time and shaded at the fixed
 *        meridian, within a fixed mean channel error.
 * @details The oriented render and the render whose lattice is rotated by an
 * explicit matrix take the same rotation path and agree exactly.
 */
inline void test_gs_lattice_basis_matches_rotated_lattice() {
  using WhiteBox = GSWhiteBox;
  constexpr int W = SMALL_W, H = SMALL_H;
  reset_globals();
  WhiteBox::GS gs;
  gs.init();
  for (int f = 0; f < 120; ++f) {
    gs.draw_frame();
    gs.advance_display();
  }
  const math::Quaternion previous =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 0.7f);
  const math::Quaternion current =
      (math::make_rotation(math::Vector(-0.6f, 0.2f, 0.7f).normalized(),
                           0.35f) *
       previous)
          .normalized();
  auto capture = [&](auto &&shade) {
    shade();
    std::vector<Pixel> frame;
    capture_frame<W, H>(gs, frame);
    return frame;
  };

  WhiteBox::set_orientation(gs, current, current);
  const std::vector<Pixel> oriented =
      capture([&] { WhiteBox::shade_field(gs, false); });
  const FrameErrorStats rounding = compare_buffers(
      oriented, capture([&] { WhiteBox::shade_matrix_oriented_field(gs); }));
  const FrameErrorStats still = compare_buffers(
      oriented, capture([&] { WhiteBox::shade_field(gs, true); }));

  WhiteBox::set_orientation(gs, previous, current);
  const std::vector<Pixel> actual =
      capture([&] { WhiteBox::shade_field(gs, true); });
  const auto expected = column_time_oracle<W, H>([&](float u) {
    const math::Quaternion at = math::slerp(previous, current, u).normalized();
    WhiteBox::set_orientation(gs, at, at);
    return capture([&] { WhiteBox::shade_field(gs, false); });
  });
  const FrameErrorStats moving = compare_buffers(expected, actual);
  std::printf("  [info] GS lattice basis: rounding mean=%.5f max=%u; still "
              "mean=%.5f max=%u; moving lit=%zu mean=%.5f max=%u\n",
              rounding.mean_absolute_error(), rounding.max_absolute_error,
              still.mean_absolute_error(), still.max_absolute_error,
              moving.lit_reference_pixels, moving.mean_absolute_error(),
              moving.max_absolute_error);
  HS_EXPECT_EQ(rounding.max_absolute_error, 0u);
  HS_EXPECT_GT(moving.lit_reference_pixels, static_cast<size_t>(W * H / 8));
  for (const FrameErrorStats *stats : {&still, &moving}) {
    HS_EXPECT_LE(stats->max_absolute_error, 16u);
    HS_EXPECT_LE(stats->mean_absolute_error(), 0.15);
  }
}

/**
 * @brief Each sweep slab shows the field after the substep that ends its share
 *        of the frame, through the dissolve mask at that substep, with the
 *        slab's colour-noise blend and the pigment of the nearer frame
 *        endpoint, on both arms; a frame that ends a dissolve stays black, a
 *        reseed frame shows the seeded field whole, and toggling temporal
 *        sampling leaves the simulation unchanged.
 */
inline void test_gs_slabs_shade_the_field_after_their_substep() {
  using WhiteBox = GSWhiteBox;
  constexpr int W = SMALL_W, H = SMALL_H, STEPS = WhiteBox::STEPS_PER_FRAME;
  constexpr int N = WhiteBox::N, GROW = 40;
  constexpr bool TOGGLED[] = {true, false, true, true, true, true, true, true};
  // Frames 0 and 1 convert halfway through a dissolve; frame 2 starts the last
  // three frames of one: two converting frames, the black frame that ends it,
  // and the reseed frame.
  constexpr int DISSOLVE_AT = 2;
  int earlier = 0, earlier_pigment = 0, earlier_dissolve = 0, masked = 0,
      black = 0, reseeded = 0;
  auto run = [&](bool check) {
    reset_globals();
    WhiteBox::GS gs;
    gs.init();
    Spacetime::set_flat_surface(!check);
    for (int f = 0; f < GROW; ++f) {
      gs.draw_frame();
      gs.advance_display();
    }
    for (int f = 0; f < static_cast<int>(std::size(TOGGLED)); ++f) {
      const bool enabled = check && TOGGLED[f];
      HS_CONTEXT("frame, temporal sampling", f, enabled);
      Spacetime::set_flat_surface(!enabled);
      if (f == 0)
        WhiteBox::dissolve_within(gs, WhiteBox::DISSOLVE_FRAMES / 2);
      if (f == DISSOLVE_AT)
        WhiteBox::dissolve_within(gs, 3);
      const WhiteBox::Snapshot before = WhiteBox::snapshot(gs);
      gs.draw_frame();
      gs.advance_display();
      if (!enabled)
        continue;
      std::vector<Pixel> frame;
      capture_frame<W, H>(gs, frame);
      std::vector<std::vector<Pixel>> fields(STEPS);
      const int dissolve = before.transition.dissolve_frames;
      if (before.transition.reseed_pending()) {
        ++reseeded;
        const WhiteBox::Snapshot after = WhiteBox::snapshot(gs);
        for (auto &field : fields)
          field = WhiteBox::shade_slab_field(gs, after.b, -1.0f, 1.0f,
                                             after.pigment);
      } else if (dissolve >= 0 && dissolve + 1 >= WhiteBox::DISSOLVE_FRAMES) {
        ++black;
        for (auto &field : fields)
          field.assign(static_cast<size_t>(W) * H, Pixel(0, 0, 0));
      } else {
        masked += dissolve >= 0;
        const WhiteBox::Snapshot after = WhiteBox::snapshot(gs);
        std::vector<float> a(N), b(N), next_a(N), next_b(N);
        for (int i = 0; i < N; ++i) {
          a[i] = WhiteBox::from_q16(before.a[i]);
          b[i] = WhiteBox::from_q16(before.b[i]);
        }
        std::vector<uint16_t> q(N);
        for (int k = 0; k < STEPS; ++k) {
          WhiteBox::step_physics(gs, a.data(), b.data(), next_a.data(),
                                 next_b.data());
          a.swap(next_a);
          b.swap(next_b);
          for (int i = 0; i < N; ++i)
            q[i] = WhiteBox::to_q16(b[i]);
          const float ends = static_cast<float>(k + 1) / STEPS;
          const float phase = dissolve >= 0
                                  ? (static_cast<float>(dissolve) + ends) /
                                        WhiteBox::DISSOLVE_FRAMES
                                  : -1.0f;
          fields[k] = WhiteBox::shade_slab_field(
              gs, q, phase, ends,
              2 * (k + 1) <= STEPS ? before.pigment : after.pigment);
          if (2 * (k + 1) <= STEPS)
            earlier_pigment +=
                WhiteBox::shade_slab_field(gs, q, phase, ends, after.pigment) !=
                fields[k];
          if (k == 0 && dissolve >= 0)
            earlier_dissolve +=
                WhiteBox::shade_slab_field(gs, q,
                                           static_cast<float>(dissolve + 1) /
                                               WhiteBox::DISSOLVE_FRAMES,
                                           ends, before.pigment) != fields[k];
        }
      }
      for (int x = 0; x < W; ++x) {
        const int slab =
            static_cast<int>(Spacetime::SamplingSurface<W>::time(x) * STEPS);
        for (int y = 0; y < H; ++y) {
          HS_CONTEXT("px", x, y);
          const size_t i = static_cast<size_t>(y) * W + x;
          HS_EXPECT_TRUE(frame[i] == fields[slab][i]);
          earlier += fields[slab][i] != fields[STEPS - 1][i];
        }
      }
    }
    Spacetime::set_flat_surface(false);
    return WhiteBox::snapshot(gs);
  };
  const WhiteBox::Snapshot toggled = run(true);
  HS_EXPECT_GT(earlier, 0);
  HS_EXPECT_GT(earlier_pigment, 0);
  HS_EXPECT_GT(earlier_dissolve, 0);
  HS_EXPECT_EQ(masked, 3);
  HS_EXPECT_EQ(black, 1);
  HS_EXPECT_EQ(reseeded, 1);
  HS_EXPECT(run(false) == toggled,
            "temporal sampling must not change the simulation");
}

/**
 * @brief Under temporal sampling the previous frame's colour-noise table is
 *        the bake at the previous phase, and a noise scale change drops it.
 */
inline void test_gs_slab_noise_keeps_the_previous_bake() {
  using WhiteBox = GSWhiteBox;
  reset_globals();
  WhiteBox::GS gs;
  gs.init();
  for (int f = 0; f < 3; ++f) {
    const float previous = WhiteBox::noise_phase(gs);
    gs.draw_frame();
    gs.advance_display();
    HS_CONTEXT("frame", f);
    HS_EXPECT_TRUE(WhiteBox::noise_history(gs));
    HS_EXPECT_TRUE(previous != WhiteBox::noise_phase(gs));
    HS_EXPECT(WhiteBox::previous_noise(gs) ==
                  WhiteBox::bake_noise(gs, previous),
              "the retained table must be the previous phase's bake");
    HS_EXPECT(WhiteBox::current_noise(gs) ==
                  WhiteBox::bake_noise(gs, WhiteBox::noise_phase(gs)),
              "the current table must be the current phase's bake");
  }
  WhiteBox::set_noise_scale(gs, 2.0f);
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_FALSE(WhiteBox::noise_history(gs));
}

/** @brief GS slab bands match full-canvas frames. */
inline void test_gs_segmented_slabs_match_full_frame() {
  expect_slab_bands_match_full_frame<GSWhiteBox::GS>(40);
}
