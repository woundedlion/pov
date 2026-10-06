/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan::rasterize_face vs scan_region, seams, pole LOD, region clip arc
// ============================================================================

/**
 * @brief Verifies rasterize_face walks the pixels scan_region walks.
 * @details The same SDF::Face drawn through Scan::rasterize and
 * Scan::rasterize_face lights the same pixels with the same coverage. Scoped to
 * Render::pole_lod_aggressiveness 0.
 */
inline void test_face_rasterize_matches_scan_region() {
  constexpr int W = 96, H = 64;
  constexpr int HV = H + hs::H_OFFSET;
  const ScopedPoleLod lod(0.0f);

  // x0 == x1 leaves the frame unclipped; a margin that underflows x0 gives the
  // seam-wrapping band.
  auto run_case = [&](const math::Vector &axis, float rho, int sides,
                      float phase, int x0, int x1, int margin) {
    HS_CONTEXT("face", sides, x0);
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    math::Vector verts[8];
    uint16_t idx[8];
    for (int i = 0; i < sides; ++i) {
      float a = (math::TWO_PI_F * i) / sides + phase;
      verts[i] = (basis.v * cosf(rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                     .normalized();
      idx[i] = static_cast<uint16_t>(i);
    }
    std::span<const math::Vector> vspan(verts, sides);
    std::span<const uint16_t> ispan(idx, sides);
    const Color4 color(Pixel(60000, 20000, 40000), 0.75f);

    auto draw = [&](hs_test::StubEffect &fx, bool fused) {
      if (x0 != x1) {
        fx.set_clip(4, H - 3, x0, x1);
        fx.set_margin(margin);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(fx);
        SDF::FaceScratchBuffer scratch;
        SDF::Face face(vspan, ispan, scratch, HV, H, &canvas.clip());
        auto shader = [&](const math::Vector &, Fragment &f) {
          f.color = color;
        };
        if (fused)
          Scan::rasterize_face<W, H, Pipeline<W, H>>(pipeline, canvas, face,
                                                     shader);
        else
          Scan::rasterize<W, H, false>(pipeline, canvas, face, shader);
      }
      fx.advance_display();
    };

    // One live Effect at a time: read the generic path back before the fused
    // fixture exists.
    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect generic(W, H);
      draw(generic, false);
      capture_frame<W, H>(generic, expected);
    }

    hs_test::StubEffect fused(W, H);
    draw(fused, true);

    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = fused.get_pixel(x, y);
        HS_CONTEXT("px", x, y);
        HS_EXPECT_EQ(static_cast<int>(a.r), static_cast<int>(b.r));
        HS_EXPECT_EQ(static_cast<int>(a.g), static_cast<int>(b.g));
        HS_EXPECT_EQ(static_cast<int>(a.b), static_cast<int>(b.b));
        if (a.r || a.g || a.b)
          ++lit;
      }
    }
    // Guard against both paths drawing nothing.
    HS_EXPECT_GT(lit, (size_t)40);
  };

  // Mid-latitude hexagon, whose azimuth wedge sits around column 80.
  const math::Vector mid = math::Vector(0.3f, 0.8f, -0.5f).normalized();
  run_case(mid, 0.45f, 6, 0.37f, 0, 0, 0);
  run_case(mid, 0.45f, 6, 0.37f, 62, 94, 0);
  // Centred on theta=0, so the wedge straddles the seam. The band [90, 96) u
  // [0, 46) that the margin wraps around the seam covers it.
  const math::Vector seam = math::Vector(1.0f, 0.15f, 0.0f).normalized();
  run_case(seam, 0.5f, 5, 0.0f, 0, 0, 0);
  run_case(seam, 0.5f, 5, 0.0f, 0, 40, 6);
  // Face over the pole: runs are rebuilt per row and the rows it cannot bound
  // fall back to a full-width scan.
  const math::Vector pole = math::Vector(0.02f, 1.0f, 0.0f).normalized();
  run_case(pole, 0.35f, 4, 0.11f, 0, 0, 0);
  run_case(pole, 0.35f, 4, 0.11f, 12, 62, 0);
  run_case(pole, 0.35f, 4, 0.11f, 0, 40, 6);
}

/**
 * @brief Verifies the scan_region seam coalescer avoids double-plotting.
 * @details A sorted two-span row (a low span plus a seam-crosser) plots the
 * wrapped columns both spans share exactly once.
 */
inline void test_scan_region_seam_no_double_plot() {
  constexpr int W = 96, H = 20;
  int counts[W];
  for (int i = 0; i < W; ++i)
    counts[i] = 0;
  const int y = 10;

  Scan::scan_region<W, H>(
      y, y,
      [](int, auto &&out) {
        out(1.0f, 5.0f);                     // low span      -> 1,2,3,4
        out((float)(W - 2), (float)(W + 2)); // seam-crosser  -> W-2,W-1,0,1
        return true;
      },
      [&](int wx, int, const math::Vector &, int run) {
        for (int i = 0; i < run; ++i)
          if (wx + i >= 0 && wx + i < W)
            counts[wx + i]++;
        return run;
      });

  // The shared column x=1 is plotted once.
  for (int x = 0; x < W; ++x)
    HS_EXPECT_LE(counts[x], 1);

  // Coverage is exactly {0,1,2,3,4, W-2, W-1}.
  const int covered[] = {0, 1, 2, 3, 4, W - 2, W - 1};
  for (int x : covered)
    HS_EXPECT_EQ(counts[x], 1);
  HS_EXPECT_EQ(counts[5], 0);
  HS_EXPECT_EQ(counts[W - 3], 0);
}

/**
 * @brief Verifies the scan_region forward coalescer handles fractional bounds.
 * @details Two abutting spans whose shared boundary falls fractionally inside
 * one pixel column plot that column once.
 */
inline void test_scan_region_fractional_boundary_no_double_plot() {
  constexpr int W = 96, H = 20;
  int counts[W];
  for (int i = 0; i < W; ++i)
    counts[i] = 0;
  const int y = 10;

  Scan::scan_region<W, H>(
      y, y,
      [](int, auto &&out) {
        out(2.0f, 5.4f); // -> 2,3,4,5
        out(5.6f, 8.0f); // floor(5.6)=5 would re-plot x=5 without the clamp
        return true;
      },
      [&](int wx, int, const math::Vector &, int run) {
        for (int i = 0; i < run; ++i)
          if (wx + i >= 0 && wx + i < W)
            counts[wx + i]++;
        return run;
      });

  for (int x = 0; x < W; ++x)
    HS_EXPECT_LE(counts[x], 1);

  // Coverage is exactly {2,3,4,5,6,7}; x=5 covered once, not twice.
  const int covered[] = {2, 3, 4, 5, 6, 7};
  for (int x : covered)
    HS_EXPECT_EQ(counts[x], 1);
  HS_EXPECT_EQ(counts[1], 0);
  HS_EXPECT_EQ(counts[8], 0);
}

/**
 * @brief Verifies near-pole runs are anchored to canvas columns, not to the
 *        span or the clip arc.
 * @details A column is settled from its own block's anchor while the whole
 * block fits inside the walked run, and walked at its own column otherwise.
 * Records the column each shade was probed at.
 */
inline void test_pole_lod_runs_are_canvas_anchored() {
  constexpr int W = 288;
  constexpr int H = 144;
  const int y = 2; // near the north pole, so sin(phi) is small and stride > 1

  math::TrigLUT<W, H>::init();

  const ScopedPoleLod lod(1.0f);
  const int lod_stride = Scan::pole_lod_run(math::TrigLUT<W, H>::sin_phi[y]);
  HS_EXPECT_GT(lod_stride, 1);

  // probe[x] = the column whose shade settled x, -1 if unvisited.
  auto scan_probes = [&](ClipRegion::XClip xc, std::array<int, W> &probe) {
    probe.fill(-1);
    Scan::scan_region<W, H>(
        y, y,
        [](int, auto &&out) {
          out(0.0f, (float)W);
          return true;
        },
        [&](int wx, int, const math::Vector &, int run) {
          for (int i = 0; i < run; ++i) {
            const int x = wx + i;
            if (x >= 0 && x < W)
              probe[static_cast<size_t>(x)] = wx;
          }
          // Every offer lies inside one canvas-anchored block, and a
          // multi-column offer starts on its block boundary.
          HS_EXPECT_EQ(wx / lod_stride, (wx + run - 1) / lod_stride);
          if (run > 1)
            HS_EXPECT_EQ(wx % lod_stride, 0);
          return run;
        },
        xc);
  };

  auto expected_probe = [&](int x, int rs, int re) {
    if (x < rs || x >= re)
      return -1;
    const int anchor = x - x % lod_stride;
    return (anchor >= rs && anchor + lod_stride <= re) ? anchor : x;
  };

  auto check = [&](ClipRegion::XClip xc, int rs, int re,
                   std::array<int, W> &probe) {
    HS_CONTEXT("arc", rs, re);
    scan_probes(xc, probe);
    for (int x = 0; x < W; ++x) {
      HS_CONTEXT("col", x);
      HS_EXPECT_EQ(probe[static_cast<size_t>(x)], expected_probe(x, rs, re));
    }
  };

  std::array<int, W> full{};
  check(ClipRegion::XClip{}, 0, W, full);

  for (int q = 0; q < 4; ++q) {
    ClipRegion clip{};
    clip.w = W;
    clip.h = H;
    clip.y_end = H;
    clip.x_start = q * (W / 4);
    clip.x_end = clip.x_start + (W / 4);
    clip.margin = 0;
    std::array<int, W> part{};
    check(clip.x_clip(), clip.x_start, clip.x_end, part);
  }

  // Arc ends inside a block: the columns the truncated blocks keep are walked
  // per column, so their shade comes from a different source than unclipped.
  ClipRegion::XClip cut{};
  cut.active = true;
  cut.wrap = false;
  cut.rs = lod_stride / 2;
  cut.re = W / 2 + lod_stride / 2 + 1;
  std::array<int, W> part{};
  check(cut, cut.rs, cut.re, part);
  size_t moved = 0;
  for (int x = cut.rs; x < cut.re; ++x)
    if (part[static_cast<size_t>(x)] != full[static_cast<size_t>(x)])
      ++moved;
  HS_EXPECT_GT(moved, (size_t)0);

  // Aggressiveness 0 is exactly one column per run.
  Render::pole_lod_aggressiveness = 0.0f;
  HS_EXPECT_EQ(Scan::pole_lod_run(math::TrigLUT<W, H>::sin_phi[y]), 1);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f), 1);
}

/**
 * @brief Pins pole_lod_run on both sides of the Render::POLE_LOD_MAX_RUN clamp.
 * @details The run is aggressiveness over |sin(phi)|, truncated: 1 wherever
 * that quotient is under 2, the quotient itself below the clamp, and
 * Render::POLE_LOD_MAX_RUN past it and at the pole, where sin(phi) is zero. At
 * aggressiveness 0 every row is one column, the pole included.
 */
inline void test_pole_lod_run_clamps_to_max_run() {
  static_assert(16 < Render::POLE_LOD_MAX_RUN);
  const ScopedPoleLod lod(1.0f);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f), 1);
  HS_EXPECT_EQ(Scan::pole_lod_run(0.75f), 1);
  HS_EXPECT_EQ(Scan::pole_lod_run(0.25f), 4);
  HS_EXPECT_EQ(Scan::pole_lod_run(-0.25f), 4);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f / 16.0f), 16);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f / Render::POLE_LOD_MAX_RUN),
               Render::POLE_LOD_MAX_RUN);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f / (2 * Render::POLE_LOD_MAX_RUN)),
               Render::POLE_LOD_MAX_RUN);
  HS_EXPECT_EQ(Scan::pole_lod_run(0.0f), Render::POLE_LOD_MAX_RUN);
  Render::pole_lod_aggressiveness = 4.0f;
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f), 4);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f / 16.0f), Render::POLE_LOD_MAX_RUN);
  Render::pole_lod_aggressiveness = 0.0f;
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f / (2 * Render::POLE_LOD_MAX_RUN)), 1);
  HS_EXPECT_EQ(Scan::pole_lod_run(0.0f), 1);
}

/**
 * @brief Verifies near-pole decimation shades the same pixels as an
 *        undecimated walk.
 * @details A block is settled from one probe only where that probe bounds the
 * whole block, so a constant-color draw lands the same framebuffer at
 * Render::pole_lod_aggressiveness 1.0 as at 0.
 */
inline void test_pole_lod_shading_matches_undecimated() {
  constexpr int W = 96, H = 64;
  constexpr int HV = H + hs::H_OFFSET;
  const ScopedPoleLod scoped_lod(0.0f);
  const Color4 color(Pixel(60000, 45000, 30000), 1.0f);

  auto readback = [](hs_test::StubEffect &fx) {
    std::vector<Pixel> out(W * H);
    capture_frame<W, H>(fx, out);
    return out;
  };

  auto compare = [&](const std::vector<Pixel> &plain,
                     const std::vector<Pixel> &decimated) {
    size_t lit = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const Pixel &a = plain[y * W + x];
        const Pixel &b = decimated[y * W + x];
        HS_CONTEXT("px", x, y);
        HS_EXPECT_EQ(static_cast<int>(a.r), static_cast<int>(b.r));
        HS_EXPECT_EQ(static_cast<int>(a.g), static_cast<int>(b.g));
        HS_EXPECT_EQ(static_cast<int>(a.b), static_cast<int>(b.b));
        if (a.r || a.g || a.b)
          ++lit;
      }
    HS_EXPECT_GT(lit, (size_t)40);
  };

  // Ring axis tilted just off the canvas pole, so the stroke band's two arcs
  // cross the rows whose stride exceeds 1.
  auto draw_ring = [&](float lod, const math::Vector &axis, float radius) {
    Render::pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      math::Basis basis = math::make_basis(math::Quaternion(), axis);
      Scan::Ring::draw<W, H, false>(
          pipe, c, basis, radius, /*thickness=*/0.05f,
          [&](const math::Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  auto draw_face = [&](float lod, float tilt, float rho, int sides, float phase,
                       bool fused) {
    Render::pole_lod_aggressiveness = lod;
    const math::Vector axis = math::Vector(tilt, 1.0f, 0.0f).normalized();
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    math::Vector verts[8];
    uint16_t idx[8];
    for (int i = 0; i < sides; ++i) {
      float a = (math::TWO_PI_F * i) / sides + phase;
      verts[i] = (basis.v * cosf(rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                     .normalized();
      idx[i] = static_cast<uint16_t>(i);
    }
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas canvas(fx);
      SDF::FaceScratchBuffer scratch;
      SDF::Face face(std::span<const math::Vector>(verts, sides),
                     std::span<const uint16_t>(idx, sides), scratch, HV, H,
                     &canvas.clip());
      auto shader = [&](const math::Vector &, Fragment &f) { f.color = color; };
      if (fused)
        Scan::rasterize_face<W, H, Pipeline<W, H>>(pipe, canvas, face, shader);
      else
        Scan::rasterize<W, H, false>(pipe, canvas, face, shader);
    }
    fx.advance_display();
    return readback(fx);
  };

  auto draw_folded = [&](auto tag, float lod, const math::Vector &axis,
                         float radius, int sides, bool typed) {
    using ShapeT = typename decltype(tag)::type;
    Render::pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      math::Basis basis = math::make_basis(math::Quaternion(), axis);
      if (typed)
        ShapeT::template draw_solid<W, H>(pipe, c, basis, radius, sides, color);
      else
        ShapeT::template draw<W, H, false>(
            pipe, c, basis, radius, sides,
            [&](const math::Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  // A ring carved out of a polygon. Past its stroke the ring reports the
  // sentinel, which loses Subtract's max, so a probe one column from a carved
  // column reports the polygon's own depth and an ungated splat paints the
  // carve shut.
  auto draw_carved = [&](float lod, bool typed) {
    Render::pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      const math::Basis poly_basis =
          math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
      const math::Vector ring_axis =
          math::Vector(0.12f, 1.0f, 0.0f).normalized();
      const math::Basis ring_basis =
          math::make_basis(math::Quaternion(), ring_axis);
      SDF::PlanarPolygon poly(poly_basis,
                              /*radius=*/0.35f / (math::PI_F / 2.0f),
                              /*sides=*/5, 0.0f);
      // Centerline through the canvas pole: the stroke crosses the decimated
      // rows the polygon's interior covers.
      SDF::Ring ring(ring_basis, /*radius=*/0.12f / (math::PI_F / 2.0f),
                     /*thickness=*/0.05f);
      SDF::Subtract<SDF::PlanarPolygon, SDF::Ring> carved(poly, ring);
      if (typed)
        Scan::rasterize_solid<W, H>(pipe, c, carved, color);
      else
        Scan::rasterize<W, H, false>(
            pipe, c, carved,
            [&](const math::Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  // A polygon repeated about the canvas pole with its centre on a sector
  // boundary. Across a boundary distance() reports the folded copy rather than
  // the nearest one. Sector boundaries converge at the pole, so a block on a
  // decimated row spans them.
  auto draw_repeat = [&](float lod, bool typed, int reps) {
    Render::pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      // Azimuth in the fold's own perpendicular frame, where the sector
      // boundaries sit at +/- PI / reps. rho holds the copies over the
      // decimated rows, close enough to the pole that the circumradius spans
      // most of a sector in azimuth.
      const float rho = 0.30f;
      const float theta = math::PI_F / static_cast<float>(reps);
      const math::Vector centre(sinf(rho) * cosf(theta), cosf(rho),
                                -sinf(rho) * sinf(theta));
      const math::Basis child_basis =
          math::make_basis(math::Quaternion(), centre);
      SDF::PlanarPolygon child(child_basis, /*circumradius=*/0.15f, /*sides=*/5,
                               0.0f);
      SDF::AngularRepeat<SDF::PlanarPolygon> rep(child, reps,
                                                 math::Vector(0, 1, 0));
      if (typed)
        Scan::rasterize_solid<W, H>(pipe, c, rep, color);
      else
        Scan::rasterize<W, H, false>(
            pipe, c, rep,
            [&](const math::Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  // A sliver, posed so the decimated rows are the ones its projection stretches
  // most: the circumradius reaches 0.7 gnomonic-plane units while the inradius
  // keeps the face on the plane-distance path, whose report runs up to 1 + r^2
  // times an angular step.
  auto draw_sliver_face = [&](float lod, bool fused, float tilt, float spin) {
    Render::pole_lod_aggressiveness = lod;
    const math::Vector axis = math::Vector(sinf(tilt), cosf(tilt), 0.0f);
    const math::Basis basis = math::make_basis(math::Quaternion(), axis);
    const float px[3] = {0.7f, -0.35f, -0.35f};
    const float py[3] = {0.0f, 0.06f, -0.06f};
    math::Vector verts[3];
    uint16_t idx[3];
    for (int i = 0; i < 3; ++i) {
      const float rx = px[i] * cosf(spin) - py[i] * sinf(spin);
      const float ry = px[i] * sinf(spin) + py[i] * cosf(spin);
      verts[i] = (basis.v + basis.u * rx + basis.w * ry).normalized();
      idx[i] = static_cast<uint16_t>(i);
    }
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas canvas(fx);
      SDF::FaceScratchBuffer scratch;
      SDF::Face face(std::span<const math::Vector>(verts, 3),
                     std::span<const uint16_t>(idx, 3), scratch, HV, H,
                     &canvas.clip());
      HS_EXPECT_TRUE(face.linear_dist);
      HS_EXPECT_GT(face.max_dist_sq, 0.25f);
      auto shader = [&](const math::Vector &, Fragment &f) { f.color = color; };
      if (fused)
        Scan::rasterize_face<W, H, Pipeline<W, H>>(pipe, canvas, face, shader);
      else
        Scan::rasterize<W, H, false>(pipe, canvas, face, shader);
    }
    fx.advance_display();
    return readback(fx);
  };

  // The rows the draws reach must actually be decimated, or the comparison is
  // vacuous.
  Render::pole_lod_aggressiveness = 1.0f;
  math::TrigLUT<W, H>::init();
  HS_EXPECT_GT(Scan::pole_lod_run(math::TrigLUT<W, H>::sin_phi[2]), 1);

  // Poses whose decimated rows carry a per-arc report rate the shared
  // change-per-arc factor alone undercuts.
  struct SliverCase {
    float tilt;
    int spin_step;
  };
  for (SliverCase sc : {SliverCase{0.32f, 9}, SliverCase{0.36f, 8},
                        SliverCase{0.46f, 7}, SliverCase{0.56f, 7}})
    for (float lod : {1.0f, 4.0f})
      for (bool fused : {false, true}) {
        const float spin = (math::TWO_PI_F * sc.spin_step) / 48.0f;
        HS_CONTEXT("sliver", static_cast<int>(sc.tilt * 100.0f),
                   static_cast<int>(lod) * 2 + fused);
        compare(draw_sliver_face(0.0f, fused, sc.tilt, spin),
                draw_sliver_face(lod, fused, sc.tilt, spin));
      }

  for (float lod : {1.0f, 4.0f})
    for (bool typed : {false, true}) {
      HS_CONTEXT("carved", static_cast<int>(lod), typed);
      compare(draw_carved(0.0f, typed), draw_carved(lod, typed));
    }

  for (int reps : {4, 5, 6})
    for (float lod : {1.0f, 4.0f})
      for (bool typed : {false, true}) {
        HS_CONTEXT("repeat", reps, static_cast<int>(lod) * 2 + typed);
        compare(draw_repeat(0.0f, typed, reps), draw_repeat(lod, typed, reps));
      }

  for (float tilt : {0.12f, 0.3f}) {
    HS_CONTEXT("ring", static_cast<int>(tilt * 100.0f));
    const math::Vector axis = math::Vector(tilt, 1.0f, 0.0f).normalized();
    const float radius = tilt / (math::PI_F / 2.0f);
    compare(draw_ring(0.0f, axis, radius), draw_ring(1.0f, axis, radius));
  }

  // Triangles: the widest gap between inradius and circumradius. Both take the
  // linear_dist path.
  struct FaceCase {
    float tilt, rho, phase;
  };
  for (FaceCase fc :
       {FaceCase{0.10f, 0.20f, 0.0f}, FaceCase{0.05f, 0.12f, 0.40f}})
    for (bool fused : {false, true}) {
      HS_CONTEXT("face", static_cast<int>(fc.tilt * 100.0f), fused);
      compare(draw_face(0.0f, fc.tilt, fc.rho, 3, fc.phase, fused),
              draw_face(1.0f, fc.tilt, fc.rho, 3, fc.phase, fused));
    }

  // Cap axis tilted off the canvas pole at a near-hemispherical radius, so the
  // rim -- where the fold's polar/sin(polar) term peaks -- crosses the
  // decimated rows.
  const math::Vector rim_axis = math::Vector(1.0f, 0.5f, 0.0f).normalized();
  // Flower's petals meet at the antipode of basis.v, so that axis sits by the
  // canvas pole: across it the sign alternates over a vanishing arc, which no
  // finite slack bounds.
  const math::Vector fold_axis = math::Vector(0.10f, -1.0f, 0.0f).normalized();
  // 4.0 lengthens every run.
  for (float lod : {1.0f, 4.0f})
    for (bool typed : {false, true}) {
      HS_CONTEXT("folded", static_cast<int>(lod), typed);
      compare(draw_folded(std::type_identity<Scan::PlanarPolygon>{}, 0.0f,
                          rim_axis, 0.99f, 3, typed),
              draw_folded(std::type_identity<Scan::PlanarPolygon>{}, lod,
                          rim_axis, 0.99f, 3, typed));
      compare(draw_folded(std::type_identity<Scan::Star>{}, 0.0f, rim_axis,
                          0.99f, 4, typed),
              draw_folded(std::type_identity<Scan::Star>{}, lod, rim_axis,
                          0.99f, 4, typed));
      for (int petals : {3, 6})
        compare(draw_folded(std::type_identity<Scan::Flower>{}, 0.0f, fold_axis,
                            0.6f, petals, typed),
                draw_folded(std::type_identity<Scan::Flower>{}, lod, fold_axis,
                            0.6f, petals, typed));
    }
}

/** @brief Concave near-pole faces retain their coverage under block settling. */
inline void test_pole_lod_concave_face_matches_undecimated() {
  constexpr int W = 288, H = 144;
  const ScopedPoleLod scoped_lod(0.0f);
  struct FaceCase {
    float tilt, rho, rho_inner, phase;
  };
  for (const FaceCase fc : {FaceCase{0.1f, 0.12f, 0.036f, math::PI_F / 3.0f},
                            FaceCase{0.45f, 0.65f, 0.30f, math::PI_F / 4.0f}}) {
    const math::Basis basis = math::make_basis(
        math::Quaternion(), math::Vector(fc.tilt, 1.0f, 0.0f).normalized());
    math::Vector vertices[8];
    uint16_t indices[8];
    for (int i = 0; i < 8; ++i) {
      const float angle = math::TWO_PI_F * i / 8.0f + fc.phase;
      const float rho = i % 2 ? fc.rho_inner : fc.rho;
      vertices[i] =
          (basis.v * cosf(rho) +
           (basis.u * cosf(angle) + basis.w * sinf(angle)) * sinf(rho))
              .normalized();
      indices[i] = static_cast<uint16_t>(i);
    }
    auto draw = [&](float lod) {
      Render::pole_lod_aggressiveness = lod;
      hs_test::StubEffect fx(W, H);
      Pipeline<W, H> pipe;
      {
        Canvas canvas(fx);
        SDF::FaceScratchBuffer scratch;
        SDF::Face face(std::span<const math::Vector>(vertices, 8),
                       std::span<const uint16_t>(indices, 8), scratch,
                       H + hs::H_OFFSET, H, &canvas.clip());
        HS_EXPECT_FALSE(face.convex);
        HS_EXPECT_EQ(face.linear_dist, fc.rho < 0.2f);
        auto shade = [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 45000, 30000), 1.0f);
        };
        Scan::rasterize_face<W, H>(pipe, canvas, face, shade);
      }
      fx.advance_display();
      std::vector<Pixel> pixels(W * H);
      capture_frame<W, H>(fx, pixels);
      return pixels;
    };
    const auto plain = draw(0.0f);
    size_t lit = 0;
    for (const Pixel &pixel : plain)
      lit += (pixel.r | pixel.g | pixel.b) != 0;
    HS_EXPECT_GT(lit, size_t(40));
    for (float lod : {1.0f, 4.0f}) {
      const auto decimated = draw(lod);
      for (size_t i = 0; i < plain.size(); ++i) {
        HS_CONTEXT("concave face", static_cast<int>(fc.rho * 100.0f), i);
        HS_EXPECT_EQ(plain[i], decimated[i]);
      }
    }
  }
}

/**
 * @brief Verifies scan_region's clip arc matches per-column XClip::clipped.
 * @details Runs the interval path (spans straddling the seam and both wrap
 * pieces) and the full-row path under a plain arc, a seam-wrapping arc, and no
 * clip; each column's coverage must equal the unclipped coverage filtered by
 * XClip::clipped, with no column walked twice.
 */
inline void test_scan_region_clip_arc_matches_predicate() {
  constexpr int W = 96, H = 20;
  const int y = 10;

  auto run = [&](ClipRegion::XClip xc, bool handled, int (&counts)[W]) {
    for (int i = 0; i < W; ++i)
      counts[i] = 0;
    Scan::scan_region<W, H>(
        y, y,
        [&](int, auto &&out) {
          if (!handled)
            return false;
          out(1.0f, 5.0f);
          out((float)(W - 4), (float)(W + 2)); // seam-crosser
          return true;
        },
        [&](int wx, int, const math::Vector &, int run) {
          for (int i = 0; i < run; ++i)
            if (wx + i >= 0 && wx + i < W)
              counts[wx + i]++;
          return run;
        },
        xc);
  };

  const ClipRegion::XClip no_clip{};
  const ClipRegion::XClip arc{3, 8, true, false};
  const ClipRegion::XClip wrap_arc{W - 2, 2, true, true};

  for (bool handled : {true, false}) {
    int reference[W];
    run(no_clip, handled, reference);
    for (int x = 0; x < W; ++x)
      HS_EXPECT_EQ(reference[x], !handled || x <= 4 || x >= W - 4 ? 1 : 0);
    for (const auto &xc : {arc, wrap_arc}) {
      int counts[W];
      run(xc, handled, counts);
      int kept = 0;
      for (int x = 0; x < W; ++x) {
        HS_EXPECT_LE(counts[x], 1);
        HS_EXPECT_EQ(counts[x], xc.clipped(x) ? 0 : reference[x]);
        kept += counts[x];
      }
      HS_EXPECT_GT(kept, 0);
    }
  }
}
