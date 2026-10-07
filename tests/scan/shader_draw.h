/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan::Shader::draw — full-sphere per-pixel shader
// ============================================================================

/** @brief A standalone bounding sphere prepares its row lookup table. */
inline void test_bounding_sphere_initializes_trig() {
  constexpr int W = 36, H = 19;
  math::TrigLUT<W, H>::initialized = false;
  Scan::BoundingSphere<W, H> bounds(math::Vector(1, 0, 0), 0.2f);
  HS_EXPECT_TRUE((math::TrigLUT<W, H>::initialized));
  int intervals = 0;
  bounds.get_intervals(H / 2, [&](float start, float end) {
    ++intervals;
    HS_EXPECT_GT(end - start, 0.0f);
    HS_EXPECT_LT(end - start, W * 0.5f);
  });
  HS_EXPECT_EQ(intervals, 1);
}

/** @brief Alpha at the cutoff is discarded; the next representable value draws. */
inline void test_min_alpha_boundary() {
  constexpr int W = 32, H = 16;
  const math::Basis BASIS = math::make_basis(math::Quaternion(), math::UP);
  SDF::PlanarPolygon shape(BASIS, 0.5f, 5, 0.0f);
  for (float alpha : {std::nextafter(Scan::MIN_ALPHA, 0.0f), Scan::MIN_ALPHA,
                      std::nextafter(Scan::MIN_ALPHA, 1.0f)}) {
    for (bool constant : {false, true}) {
      hs_test::StubEffect fx(W, H);
      {
        Canvas canvas(fx);
        Pipeline<W, H> pipe;
        const Color4 COLOR(Pixel(65535, 65535, 65535), alpha);
        if (constant)
          Scan::rasterize_solid<W, H>(pipe, canvas, shape, COLOR);
        else
          Scan::rasterize<W, H, false>(
              pipe, canvas, shape,
              [&](const math::Vector &, Fragment &f) { f.color = COLOR; });
      }
      fx.advance_display();
      const size_t lit = count_lit_region<W, H>(fx);
      HS_EXPECT_EQ(lit > 0, alpha > Scan::MIN_ALPHA);
    }
  }
}

/**
 * @brief Verifies a constant-color shader fills every pixel of the full sphere.
 */
inline void test_shader_constant_fills_canvas() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(c, [](const math::Vector &) {
      return Color4(Pixel(40000, 20000, 10000), 1.0f);
    });
  }
  fx.advance_display();

  // At alpha=1 the blend is the identity, so channels survive verbatim.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      HS_EXPECT_EQ((int)p.r, 40000);
      HS_EXPECT_EQ((int)p.g, 20000);
      HS_EXPECT_EQ((int)p.b, 10000);
    }
  }
}

/**
 * @brief Verifies SAMPLES==4 SSAA premultiplies each sub-sample before averaging.
 */
inline void test_shader_ssaa_premultiplies_partial_coverage() {
  constexpr int W = 16, H = 8;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    // Opacity keys on sub-sample position: the 2x2 grid's +/-0.25 px x-offsets
    // land at theta-grid phase 0.25 vs 0.75, so two of four samples are opaque.
    Scan::Shader::draw<W, H, 4>(c, [](const math::Vector &v) -> Color4 {
      float theta = std::atan2(v.z, v.x);
      if (theta < 0.0f)
        theta += 2.0f * math::PI_F;
      float g = theta * W / (2.0f * math::PI_F);
      float frac = g - std::floor(g);
      bool opaque = frac < 0.5f;
      return opaque ? Color4(Pixel(60000, 0, 0), 1.0f)
                    : Color4(Pixel(0, 0, 0), 0.0f);
    });
  }
  fx.advance_display();

  // Premultiplied: (60000*1 + 60000*1 + 0 + 0) / 4 = 30000.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      HS_EXPECT_NEAR((int)p.r, 30000, 2);
      HS_EXPECT_EQ((int)p.g, 0);
      HS_EXPECT_EQ((int)p.b, 0);
    }
  }
}

/**
 * @brief Verifies the split vertex/fragment draw averages SAMPLES==4 sub-samples
 *        over one per-pixel vertex seed.
 * @details Each of the four sub-fragments inherits the registers the vertex
 * shader seeded at the pixel center.
 */
inline void test_shader_split_ssaa_averages_subsamples() {
  constexpr int W = 16, H = 8;
  int vertex_calls = 0;
  int fragment_calls = 0;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 4>(
        c,
        [&](const math::Vector &v, Fragment &f) {
          ++fragment_calls;
          float theta = std::atan2(v.z, v.x);
          if (theta < 0.0f)
            theta += 2.0f * math::PI_F;
          float g = theta * W / (2.0f * math::PI_F);
          bool opaque = (g - std::floor(g)) < 0.5f;
          // Green carries the vertex seed: a sub-fragment that missed it reads 0.
          f.color = opaque
                        ? Color4(Pixel(0, static_cast<uint16_t>(f.v0), 0), 1.0f)
                        : Color4(Pixel(0, 0, 0), 0.0f);
        },
        [&](Fragment &f) {
          ++vertex_calls;
          f.v0 = 60000.0f;
        });
  }
  fx.advance_display();

  HS_EXPECT_EQ(vertex_calls, W * H);
  HS_EXPECT_EQ(fragment_calls, 4 * W * H);

  // Premultiplied: (60000*1 + 60000*1 + 0 + 0) / 4 = 30000.
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      HS_CONTEXT("px", x, y);
      const Pixel &p = fx.get_pixel(x, y);
      HS_EXPECT_NEAR((int)p.g, 30000, 2);
      HS_EXPECT_EQ((int)p.r, 0);
      HS_EXPECT_EQ((int)p.b, 0);
    }
  }
}

/**
 * @brief Verifies a position-reading shader maps the sphere's latitude.
 * @details The +Y pole (top row) renders brighter than the -Y pole (bottom row).
 */
inline void test_shader_positional_maps_latitude() {
  constexpr int W = 32, H = 32;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    // Green encodes latitude: north pole (v.y≈+1) bright, south (v.y≈-1) dark.
    Scan::Shader::draw<W, H, 1>(c, [](const math::Vector &v) {
      uint16_t g = (uint16_t)((v.y * 0.5f + 0.5f) * 60000.0f);
      return Color4(Pixel(0, g, 0), 1.0f);
    });
  }
  fx.advance_display();

  // Top row (y=0) is the north pole, brighter than the bottom row.
  HS_EXPECT_GT((int)fx.get_pixel(0, 0).g, (int)fx.get_pixel(0, H - 1).g);
  HS_EXPECT_GT((int)fx.get_pixel(0, 0).g, 40000);
  HS_EXPECT_LT((int)fx.get_pixel(0, H - 1).g, 20000);
}

/**
 * @brief Verifies the shader writes only inside the active clip band.
 */
inline void test_shader_respects_clip_band() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  fx.set_clip(5, 10, 0, W); // rows [5,10)
  fx.set_margin(0);         // no render-margin expansion

  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(c, [](const math::Vector &) {
      return Color4(Pixel(0, 0, 50000), 1.0f);
    });
  }
  fx.advance_display();

  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      HS_EXPECT_EQ(!is_black(fx.get_pixel(x, y)), y >= 5 && y < 10);
}

/** @brief SSAA grids sample at the expected subpixel positions. */
inline void test_ssaa_grid_sample_positions() {
  constexpr int W = 32, H = 16;
  math::TrigLUT<W, H>::init();
  Scan::Shader::SsaaGrid<W, H> grid;
  for (int y = 1; y < H - 1; ++y) {
    grid.set_row(y);
    for (int x = 0; x < W; ++x)
      for (int i = 0; i < 4; ++i) {
        const math::Vector expected = math::pixel_to_vector<W, H>(
            x + ((i & 1) ? -0.25f : 0.25f), y + ((i & 2) ? -0.25f : 0.25f));
        const math::Vector actual = grid.at(x, i);
        HS_EXPECT_NEAR(actual.x, expected.x, 1e-6f);
        HS_EXPECT_NEAR(actual.y, expected.y, 1e-6f);
        HS_EXPECT_NEAR(actual.z, expected.z, 1e-6f);
      }
  }
}

/**
 * @brief Verifies every Shader entry point paints exactly the columns the clip
 *        arc admits, with the values an unclipped draw produced.
 * @details Compared against the unclipped render filtered by XClip::clipped,
 * under a plain arc and one whose margin wraps the band across the seam.
 */
inline void test_shader_clip_arc_matches_predicate() {
  constexpr int W = 32, H = 16;

  auto positional = [](const math::Vector &v) {
    return Color4(Pixel(static_cast<uint16_t>((v.x * 0.5f + 0.5f) * 60000.0f),
                        static_cast<uint16_t>((v.z * 0.5f + 0.5f) * 60000.0f),
                        30000),
                  1.0f);
  };

  auto draw_variant = [&](Canvas &c, int variant) {
    switch (variant) {
    case 0:
      Scan::Shader::draw<W, H, 1>(c, positional);
      break;
    case 1:
      Scan::Shader::draw<W, H, 4>(c, positional);
      break;
    case 2:
      Scan::Shader::draw<W, H, 1>(
          c,
          [&](const math::Vector &v, Fragment &f) { f.color = positional(v); },
          [](Fragment &f) { f.v0 = 1.0f; });
      break;
    case 3:
      Scan::Shader::draw<W, H, 4>(
          c,
          [&](const math::Vector &v, Fragment &f) { f.color = positional(v); },
          [](Fragment &f) { f.v0 = 1.0f; });
      break;
    case 5:
      Scan::Shader::draw_cached<W, H>(
          c,
          [&](const math::Vector &v, int x, int y) {
            const auto EXPECTED = math::pixel_to_vector<W, H>(x, y);
            HS_EXPECT_NEAR(v.x, EXPECTED.x, 1e-6f);
            HS_EXPECT_NEAR(v.y, EXPECTED.y, 1e-6f);
            HS_EXPECT_NEAR(v.z, EXPECTED.z, 1e-6f);
            return positional(v);
          },
          [](int) {});
      break;
    case 6: {
      int row = -1;
      Scan::Shader::draw_cached<W, H>(
          c,
          [&](const math::Vector &v, int, int y) {
            HS_EXPECT_EQ(row, y);
            const Color4 COLOR = positional(v);
            return COLOR.color * COLOR.alpha;
          },
          [&](int y) {
            HS_EXPECT_EQ(y, row + 1);
            row = y;
          });
      HS_EXPECT_EQ(row, H - 1);
      break;
    }
    case 7:
      Scan::Shader::walk_grid<W, H>(
          c, [&](const math::Vector &center, const auto &grid, int x) {
            const Color4 COLOR = positional(grid.at(x, 0));
            HS_EXPECT_NEAR(center.magnitude(), 1.0f, 1e-6f);
            return COLOR.color * COLOR.alpha;
          });
      break;
    default:
      Scan::Shader::draw_grid<W, H>(
          c, [](Fragment &) {},
          [&](Fragment &, const Scan::Shader::SsaaGrid<W, H> &grid, int x) {
            Color4 s = positional(grid.at(x, 0));
            return s.color * s.alpha;
          });
      break;
    }
  };

  auto render = [&](int variant, int x0, int x1, int margin,
                    std::vector<Pixel> &out) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(0, H, x0, x1);
    fx.set_margin(margin);
    {
      Canvas c(fx);
      draw_variant(c, variant);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, out);
  };

  struct Band {
    int x0, x1, margin;
  };
  // A plain arc, and one whose margin underflows column 0 into a wrapped band.
  const Band bands[] = {{8, 20, 2}, {0, 10, 3}};

  for (int variant = 0; variant < 8; ++variant) {
    HS_CONTEXT("variant", variant, 0);
    // One live Effect at a time: read the unclipped render back before the
    // clipped fixture exists.
    std::vector<Pixel> reference;
    render(variant, 0, W, 0, reference);
    size_t painted = 0;
    for (const Pixel &p : reference)
      painted += p.b == 30000;
    HS_EXPECT_EQ(painted, static_cast<size_t>(W * H));

    for (const Band &b : bands) {
      ClipRegion cr;
      cr.w = W;
      cr.h = H;
      cr.y_end = H;
      cr.x_start = b.x0;
      cr.x_end = b.x1;
      cr.margin = b.margin;
      const ClipRegion::XClip xc = cr.x_clip();
      HS_EXPECT_TRUE(xc.active);

      std::vector<Pixel> got;
      render(variant, b.x0, b.x1, b.margin, got);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          HS_CONTEXT("px", x, y);
          // Outside the arc the canvas keeps the frame clear.
          const Pixel expected =
              xc.clipped(x) ? Pixel(0, 0, 0) : reference[y * W + x];
          const Pixel &g = got[y * W + x];
          HS_EXPECT_EQ((int)g.r, (int)expected.r);
          HS_EXPECT_EQ((int)g.g, (int)expected.g);
          HS_EXPECT_EQ((int)g.b, (int)expected.b);
        }
    }
  }
}
