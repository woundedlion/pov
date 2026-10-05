/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// ClipRegion::arcs_overlap + the col-span helpers — column-arc clip cull
// ============================================================================

/**
 * @brief Pins the shared cylindrical arc-overlap helper across the wrap
 *        topologies: disjoint, overlapping, seam-crossing, containment,
 *        full-width, and empty arcs.
 */
inline void test_clip_arcs_overlap() {
  constexpr int W = 96;
  // Disjoint / touching / overlapping, no wrap.
  HS_EXPECT_FALSE(ClipRegion::arcs_overlap(10, 20, 40, 20, W));
  HS_EXPECT_FALSE(ClipRegion::arcs_overlap(10, 20, 30, 10, W)); // touch at 30
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(10, 21, 30, 10, W));  // share col 30
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(10, 20, 5, 10, W));
  // One arc contains the other.
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(10, 50, 20, 5, W));
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(20, 5, 10, 50, W));
  // Seam-crossing arc [90, 96) U [0, 10).
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(90, 16, 0, 5, W));
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(90, 16, 92, 2, W));
  HS_EXPECT_FALSE(ClipRegion::arcs_overlap(90, 16, 10, 20, W));
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(0, 5, 90, 16, W)); // symmetric
  // Full-width and empty arcs.
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(0, W, 50, 1, W));
  HS_EXPECT_TRUE(ClipRegion::arcs_overlap(30, W + 5, 0, 1, W));
  HS_EXPECT_FALSE(ClipRegion::arcs_overlap(10, 0, 10, 5, W));
  HS_EXPECT_FALSE(ClipRegion::arcs_overlap(10, 5, 10, 0, W));
}

/**
 * @brief geodesic_col_span declines an arc pole within AXIS_Y_EPS of horizontal.
 */
inline void test_col_span_rejects_ill_conditioned_pole() {
  const math::Vector a(1.0f, 0.0f, 0.0f);
  const math::Vector b = math::Vector(-1.0f, 0.0002f, 0.00000001f).normalized();
  const auto span = Plot::make_geodesic_edge_span(a, b);
  HS_EXPECT_FALSE(span.antipodal);
  HS_EXPECT_LT(std::abs(span.axis.y), Plot::AXIS_Y_EPS);
  int start = 0, length = 0;
  HS_EXPECT_FALSE(Plot::geodesic_col_span<288>(a, b, span, start, length));
}

/**
 * @brief Verifies the col-span helpers conservatively cover the rendered arc's
 *        screen-column sweep and stay within the half-width sweep bound.
 * @details Densely samples the renderer's own circle (same axis selection as
 *          rasterize_geodesic_strategy) and asserts modular containment of
 *          every sample column, plus the antipodal-symmetry bound: a geodesic
 *          arc sweeps at most half the canvas in longitude, so every span must
 *          fit W/2 plus padding. Near-half sweeps (pole-grazing near-meridian
 *          arcs) stress the direction logic: covering the wrong side of an
 *          ambiguous near-half separation fails on the mid-arc samples.
 *          Non-vacuity counters require many genuinely cullable spans and many
 *          near-half sweeps.
 *          The ~7.2M dense containment samples aggregate into escape counters;
 *          sample counters keep a sweep that stopped generating samples from
 *          passing vacuously.
 */
inline void test_col_span_covers_arc() {
  constexpr int TW = 288;
  auto col_of = [](const math::Vector &v) {
    return math::vector_to_theta<TW>(v.normalized());
  };
  auto contains = [](int s, int len, float c) {
    int ci = static_cast<int>(floorf(c));
    int d = ((ci - s) % TW + TW) % TW;
    return d < len;
  };

  hs::random().seed(20260714);
  int cullable = 0;  // spans narrower than half the canvas
  int near_half = 0; // sweeps close to the half-width bound
  int fallbacks = 0; // near-meridian edges that decline to bound
  int geodesic_escapes = 0, geodesic_samples = 0;
  int over_bound = 0;

  for (int trial = 0; trial < 4000; ++trial) {
    math::Vector a, b;
    if (trial % 3 == 2) {
      // Near-meridian circle: an axis close to the equator plane produces
      // pole-grazing arcs whose longitude sweeps far past the endpoints.
      const float adx = hs::rand_f(-1, 1);
      const float ady = hs::rand_f(-0.05f, 0.05f);
      const float adz = hs::rand_f(-1, 1);
      math::Vector ad(adx, ady, adz);
      if (ad.length() < 0.1f)
        continue;
      math::Basis cb = basis_from_normal(ad.normalized());
      float a0 = hs::rand_f(0, 2 * math::PI_F);
      float a1 = a0 + hs::rand_f(0.5f, 3.0f);
      a = (cb.u * cosf(a0) + cb.w * sinf(a0)).normalized();
      b = (cb.u * cosf(a1) + cb.w * sinf(a1)).normalized();
    } else {
      a = rand_unit();
      b = rand_unit();
    }
    float ang = math::angle_between(a, b);
    if (ang < 0.05f)
      continue;

    int s, len;
    if (!Plot::geodesic_col_span<TW>(a, b, Plot::make_geodesic_edge_span(a, b),
                                     s, len)) {
      fallbacks++; // meridian fallback skips the cull; nothing to verify
      continue;
    }

    // Dense ground truth along the renderer's own circle.
    math::Vector axis = Plot::make_geodesic_edge_span(a, b).axis;
    math::Vector vperp = math::cross(axis, a);
    constexpr int N = 1000;
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector p = a * cosf(ang * t) + vperp * sinf(ang * t);
      ++geodesic_samples;
      geodesic_escapes += !contains(s, len, col_of(p));
    }

    // Antipodal symmetry caps a geodesic arc's longitude sweep at half the
    // canvas; the span may only exceed it by its own padding.
    over_bound += len > TW / 2 + 6;
    if (len < TW / 2)
      cullable++;
    if (len > TW / 2 - 10)
      near_half++;
  }
  HS_EXPECT_EQ(geodesic_escapes, 0);
  HS_EXPECT_GT(geodesic_samples, 1000000);
  HS_EXPECT_EQ(over_bound, 0);
  HS_EXPECT_GT(cullable, 500);
  HS_EXPECT_GT(near_half, 50);
  HS_EXPECT_LT(fallbacks, 400); // the guard must stay a rare escape hatch

  // Exact-antipodal edges: the span must cover the semicircle the renderer
  // bulges about stable_perpendicular_axis, not just the endpoint columns.
  int antipodal_escapes = 0, antipodal_samples = 0;
  for (int trial = 0; trial < 500; ++trial) {
    math::Vector a = rand_unit();
    math::Vector b = a * -1.0f;
    int s, len;
    if (!Plot::geodesic_col_span<TW>(a, b, Plot::make_geodesic_edge_span(a, b),
                                     s, len))
      continue;
    math::Vector axis = Plot::stable_perpendicular_axis(a);
    math::Vector vperp = math::cross(axis, a);
    constexpr int N = 1000;
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector p = a * cosf(math::PI_F * t) + vperp * sinf(math::PI_F * t);
      ++antipodal_samples;
      antipodal_escapes += !contains(s, len, col_of(p));
    }
  }
  HS_EXPECT_EQ(antipodal_escapes, 0);
  HS_EXPECT_GT(antipodal_samples, 400000);

  // Planar (azimuthal-equidistant) edges: ground truth is the chart line the
  // renderer walks, unprojected densely. Charts centered near a pole force the
  // near-pole fallback; the rest must produce bounded, containing spans.
  int planar_bounded = 0, planar_fallbacks = 0;
  int planar_escapes = 0, planar_samples = 0;
  for (int trial = 0; trial < 3000; ++trial) {
    math::Vector center = rand_unit();
    math::Basis basis = basis_from_normal(center);
    math::Vector a, b;
    random_disk_edge(basis, a, b);
    // Antipode-seam segments render geodesic (use_planar is false there).
    if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;

    int s, len;
    if (!Plot::planar_col_span<TW>(
            a, basis, Plot::make_planar_edge_span(a, b, basis), s, len)) {
      planar_fallbacks++;
      continue;
    }
    planar_bounded++;

    auto p1 = Plot::azimuthal_project(a, basis);
    auto p2 = Plot::azimuthal_project(b, basis);
    constexpr int N = 1000;
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector p = Plot::azimuthal_unproject(
          p1.first + (p2.first - p1.first) * t,
          p1.second + (p2.second - p1.second) * t, basis);
      ++planar_samples;
      planar_escapes += !contains(s, len, col_of(p));
    }
  }
  HS_EXPECT_EQ(planar_escapes, 0);
  HS_EXPECT_GT(planar_samples, 1500000);
  // Both outcomes must be exercised: bounded spans (the cull works) and the
  // near-pole/short-way fallbacks (the escape hatch fires when it must).
  HS_EXPECT_GT(planar_bounded, 1500);
  HS_EXPECT_GT(planar_fallbacks, 20);
}

// has_world_cull gates ParticleSystem::draw's hoisted per-point gate path: it
// must be false for screen-only pipelines and true whenever any stage re-emits
// clip-cull edges (cull_edge).
static_assert(!Pipeline<96, 48>::has_world_cull);
static_assert(
    !Pipeline<96, 48, Filter::Screen::AntiAlias<96, 48>>::has_world_cull);
static_assert(Pipeline<96, 48, Filter::World::Orient>::has_world_cull);
static_assert(Pipeline<96, 48, Filter::World::Orient,
                       Filter::Screen::AntiAlias<96, 48>>::has_world_cull);

// has_world_stage gates rasterize's precomputed screen-coordinate shortcut: it
// must be true for every world-space stage, including the ones that define no
// cull_edge.
static_assert(!Pipeline<96, 48>::has_world_stage);
static_assert(
    !Pipeline<96, 48, Filter::Screen::AntiAlias<96, 48>>::has_world_stage);
static_assert(!Filter::Screen::DirectAntiAliasSink<96, 48>::has_world_stage);
static_assert(Pipeline<96, 48, Filter::World::Mobius>::has_world_stage);
static_assert(Pipeline<96, 48, Filter::World::Hole,
                       Filter::Screen::AntiAlias<96, 48>>::has_world_stage);
static_assert(Pipeline<96, 48, Filter::World::Orient>::has_world_stage);

/**
 * @brief Pins edge_visible_in_clip's decision to the composed row-span /
 *        col-span cull: y-reject first, then the column arc, with either span's
 *        no-bound fallback reading as visible.
 * @details Clip bands cover the device quadrant shapes: seam-wrapping
 *          (margin pushes rs past the seam), interior non-wrapping, full-width
 *          x (XClip inactive), and the full canvas. The geodesic corpus
 *          includes antipodal, near-collapsed, and near-meridian edges; the
 *          planar corpus draws chart-line edges on random disks, including
 *          near-pole charts that force the col-span fallback.
 */
inline void test_edge_visible_in_clip_matches_span_composition() {
  constexpr int TW = 288, TH = 144;
  Pipeline<TW, TH> sink;

  hs::random().seed(20260716);
  const int bands[][4] = {
      {0, 72, 0, 144},     // segment quadrant; margin wraps rs past the seam
      {36, 108, 144, 288}, // opposite half; also seam-wrapping via margin
      {0, 144, 60, 200},   // interior wedge, non-wrapping
      {100, 144, 0, 288},  // full-width x: XClip inactive, y-only cull
      {0, 144, 0, 288},    // full canvas
  };
  int visible = 0, culled = 0;
  for (const auto &bd : bands) {
    ClipRegion cr;
    cr.w = TW;
    cr.h = TH;
    cr.y_start = bd[0];
    cr.y_end = bd[1];
    cr.x_start = bd[2];
    cr.x_end = bd[3];
    const auto xc = cr.x_clip();
    const int band_len = xc.length(TW);

    for (int trial = 0; trial < 2000; ++trial) {
      math::Vector a = rand_unit();
      math::Vector b;
      switch (trial % 7) {
      case 0:
        b = a * -1.0f; // antipodal
        break;
      case 1: // near-collapsed
        b = (a + math::Vector(1e-4f, 0.0f, 0.0f)).normalized();
        break;
      case 2: { // near-meridian arc (pole-grazing, axis.y ~ 0)
        const float adx = hs::rand_f(-1, 1);
        const float ady = hs::rand_f(-0.05f, 0.05f);
        const float adz = hs::rand_f(-1, 1);
        math::Vector ad(adx, ady, adz);
        if (ad.length() < 0.1f) {
          b = rand_unit();
          break;
        }
        math::Basis cb = basis_from_normal(ad.normalized());
        float a0 = hs::rand_f(0, 2 * math::PI_F);
        a = (cb.u * cosf(a0) + cb.w * sinf(a0)).normalized();
        b = (cb.u * cosf(a0 + 1.5f) + cb.w * sinf(a0 + 1.5f)).normalized();
        break;
      }
      default:
        b = rand_unit();
      }

      const bool got =
          Plot::edge_visible_in_clip<TW, TH>(sink, cr, xc, a, b, nullptr);
      const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a, b);
      float row_lo, row_hi;
      Plot::geodesic_row_span<TH>(a, b, es, row_lo, row_hi);
      bool want;
      if (!cr.could_intersect_y(row_lo, row_hi + Plot::GEODESIC_ROW_AA_PAD)) {
        want = false;
      } else if (!xc.active) {
        want = true;
      } else {
        int col_s, col_len;
        want = !Plot::geodesic_col_span<TW>(a, b, es, col_s, col_len) ||
               ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, TW);
      }
      HS_EXPECT_TRUE(got == want);
      (got ? visible : culled)++;
    }

    for (int trial = 0; trial < 2000; ++trial) {
      math::Vector center = rand_unit();
      math::Basis basis = basis_from_normal(center);
      math::Vector a, b;
      random_disk_edge(basis, a, b);
      // Antipode-seam segments render geodesic (use_planar is false there).
      if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
          math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
        continue;

      const bool got =
          Plot::edge_visible_in_clip<TW, TH>(sink, cr, xc, a, b, &basis);
      const Plot::PlanarEdgeSpan ps = Plot::make_planar_edge_span(a, b, basis);
      float row_lo, row_hi;
      Plot::planar_row_span<TH>(a, b, ps, row_lo, row_hi);
      bool want;
      if (!cr.could_intersect_y(row_lo, row_hi)) {
        want = false;
      } else if (!xc.active) {
        want = true;
      } else {
        int col_s, col_len;
        want = !Plot::planar_col_span<TW>(a, basis, ps, col_s, col_len) ||
               ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, TW);
      }
      HS_EXPECT_TRUE(got == want);
      (got ? visible : culled)++;
    }
  }
  // Both verdicts must be exercised across the band topologies.
  HS_EXPECT_GT(visible, 2000);
  HS_EXPECT_GT(culled, 2000);
}

/** @brief Pixel counts from a clipped render-band/reference comparison. */
struct RenderBandDiff {
  int lit = 0;            /**< Lit reference pixels in the render band. */
  int margin_lit = 0;     /**< Lit reference pixels outside the display band. */
  int row_margin_lit = 0; /**< Lit reference pixels above or below display. */
  int diff = 0;           /**< Pixels differing from the reference. */
  int margin_diff = 0;    /**< Differing pixels outside the display band. */
  int first_x = -1;       /**< Column of the first difference, or -1. */
  int first_y = -1;       /**< Row of the first difference, or -1. */
};

/** @brief Initializes every pixel in the pending frame to black. */
template <int W, int H> inline void initialize_parity_frame(Canvas &canvas) {
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      canvas(x, y) = Pixel{};
}

/** @brief Compares a frame with its full-canvas reference over the render band. */
template <int W>
inline RenderBandDiff render_band_diff(const hs_test::StubEffect &fx,
                                       const std::vector<Pixel> &reference) {
  RenderBandDiff out;
  const ClipRegion &clip = fx.clip();
  for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
    for (int x = 0; x < W; ++x) {
      if (!clip.contains_x(x))
        continue;
      const Pixel &actual = fx.get_pixel(x, y);
      const Pixel &expected = reference[static_cast<size_t>(y) * W + x];
      const bool in_display_row = y >= clip.y_start && y < clip.y_end;
      const bool in_display =
          in_display_row && x >= clip.x_start && x < clip.x_end;
      const bool lit = (expected.r | expected.g | expected.b) != 0;
      out.lit += lit;
      out.margin_lit += !in_display && lit;
      out.row_margin_lit += !in_display_row && lit;
      if (actual == expected)
        continue;
      if (out.first_x < 0) {
        out.first_x = x;
        out.first_y = y;
      }
      ++out.diff;
      out.margin_diff += !in_display;
    }
  }
  return out;
}

/** @brief Reports one render-band parity result with its first mismatch. */
inline void expect_render_band_parity(const char *label,
                                      const RenderBandDiff &diff) {
  if (diff.diff != 0)
    std::printf("  [%s] diff=%d margin_diff=%d first=(%d,%d)\n", label,
                diff.diff, diff.margin_diff, diff.first_x, diff.first_y);
  HS_CONTEXT(label, diff.first_x, diff.first_y);
  HS_EXPECT_EQ(diff.diff, 0);
}

/**
 * @brief End-to-end conservativeness of the wireframe edge gate: a
 *        quadrant/wedge-clipped Plot::Mesh::draw reproduces the full render's
 *        strokes inside the render band, including its margin ring.
 * @details Covers the whole-edge reject and the per-piece bits that replace
 *          rasterize's own cull. Random orientations put edges across both
 *          poles and the seam; an over-cull drops a whole stroke, which shows
 *          as a long run of unlit pixels along a row.
 *
 *          Cut pieces collapse into one windowed segment, preserving the full
 *          edge's step schedule and sample positions inside the render band.
 */
inline void test_mesh_edge_gate_pixel_parity() {
  constexpr int W = 96, H = 48;
  configure_arenas_default();

  alignas(32) static uint8_t seed_a[24 * 1024];
  alignas(32) static uint8_t seed_b[24 * 1024];
  alignas(32) static uint8_t geom[16 * 1024];
  Arena sa(seed_a, sizeof(seed_a));
  Arena sb(seed_b, sizeof(seed_b));
  Arena ga(geom, sizeof(geom));

  MeshState mesh;
  build_icosahedron_meshstate(sa, sb, ga, mesh);

  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(ga, mesh.faces.size());
  Plot::Mesh::extract_edges(mesh, edges);

  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.9f);
  };

  hs::random().seed(0x5EED);
  const int clips[4][4] = {
      {0, H / 2, 0, W / 2}, {H / 2, H, W / 2, W}, {0, H, 10, 34}, {0, H, 0, 8}};
  int lit_total = 0, margin_lit_total = 0;

  MeshState posed;
  posed.vertices.bind(ga, mesh.vertices.size());
  for (size_t i = 0; i < mesh.vertices.size(); ++i)
    posed.vertices.push_back(mesh.vertices[i]);

  for (int trial = 0; trial < 12; ++trial) {
    // Re-orient the shell so edges sweep the poles and the wrap seam.
    const float axis_x = hs::rand_f(-1, 1);
    const float axis_y = hs::rand_f(-1, 1);
    const float axis_z = hs::rand_f(-1, 1);
    math::Vector axis = math::Vector(axis_x, axis_y, axis_z);
    if (axis.length() < 0.1f)
      axis = math::Y_AXIS;
    math::Quaternion q =
        math::make_rotation(axis.normalized(), hs::rand_f(0, 2 * math::PI_F));
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
      posed.vertices[i] = math::rotate(mesh.vertices[i], q);

    auto render = [&](hs_test::StubEffect &fx) {
      Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters{
          Filter::Screen::AntiAlias<W, H>()};
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::Mesh::draw<W, H>(filters, c, posed, edges, shade);
    };

    std::vector<Pixel> ref(static_cast<size_t>(W) * H);
    {
      hs_test::StubEffect fx(W, H);
      render(fx);
      fx.advance_display();
      hs_test::capture_frame<W, H>(fx, ref);
    }

    for (auto &cl : clips) {
      hs_test::StubEffect fx(W, H);
      fx.set_clip(cl[0], cl[1], cl[2], cl[3]);
      render(fx);
      fx.advance_display();
      const ClipRegion &clip = fx.clip();
      for (const auto &edge : edges) {
        const math::PixelCoords endpoint =
            math::vector_to_pixel<W, H>(posed.vertices[edge.v]);
        const int x = static_cast<int>(endpoint.x);
        const int y = static_cast<int>(endpoint.y);
        if (!clip.contains_x(x) || !clip.contains_y(y))
          continue;
        const Pixel &expected = ref[static_cast<size_t>(y) * W + x];
        if ((expected.r | expected.g | expected.b) == 0)
          continue;
        const Pixel &actual = fx.get_pixel(x, y);
        HS_CONTEXT("mesh endpoint", edge.u, edge.v);
        HS_EXPECT_TRUE((actual.r | actual.g | actual.b) != 0);
      }
      const RenderBandDiff diff = render_band_diff<W>(fx, ref);
      expect_render_band_parity("mesh edge gate", diff);
      lit_total += diff.lit;
      margin_lit_total += diff.margin_lit;
    }
  }
  HS_EXPECT_GT(lit_total, 200);
  HS_EXPECT_GT(margin_lit_total, 20);
}

/** @brief Short geodesic spans report their true, monotone arc length. */
inline void test_short_geodesic_arc_lengths() {
  float previous = 0.0f;
  for (float angle : {0.0001f, 0.0002f, 0.00025f, 0.0004f, 0.0005f}) {
    const math::Vector a(1.0f, 0.0f, 0.0f);
    const math::Vector b(cosf(angle), 0.0f, sinf(angle));
    const auto span = Plot::make_geodesic_edge_span(a, b);
    HS_EXPECT_NEAR(span.total, angle, 1e-9f);
    HS_EXPECT_GT(span.total, previous);
    previous = span.total;
  }
}

/** @brief Endpoint shortcuts obey the same arc window as the adaptive walk. */
inline void test_rasterize_short_edge_windows() {
  constexpr int W = 96, H = 48;
  const math::Basis BASIS =
      math::make_basis(math::Quaternion(), math::Vector(1, 0, 0));
  for (bool planar : {false, true}) {
    for (float distance : {0.0f, 0.0001f}) {
      ScratchScope scope(plot_arena());
      Fragments points;
      points.bind(plot_arena(), 2);
      Fragment start, end;
      start.pos = math::Vector(1, 0, 0);
      end.pos = math::Vector(cosf(distance), 0, sinf(distance));
      points.push_back(start);
      points.push_back(end);
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      CapturePipeline sink;
      int shaded = 0;
      auto shader = [&](const math::Vector &, Fragment &) { ++shaded; };
      Plot::RasterOptions options;
      if (planar)
        options.projection = Plot::RasterProjection::planar(BASIS);
      options.plot_t_start = 0.25f;
      options.plot_t_end = 0.75f;
      Plot::rasterize<W, H>(sink, canvas, points, shader, options);
      HS_EXPECT_EQ(shaded, 0);
      HS_EXPECT_TRUE(sink.plotted.empty());
      options.plot_t_start = 0;
      options.plot_t_end = 1;
      Plot::rasterize<W, H>(sink, canvas, points, shader, options);
      HS_EXPECT_GT(shaded, 0);
      HS_EXPECT_FALSE(sink.plotted.empty());
    }
  }
}

/** @brief An open upper window retains the whole edge's terminal sample. */
inline void test_rasterize_window_preserves_terminal_sample() {
  constexpr int W = 96, H = 48;
  const std::pair<math::Vector, math::Vector> EDGES[] = {
      {math::Vector(-0x1.6a5c08p-1f, -0x1.38938ap-2f, 0x1.463602p-1f),
       math::Vector(0x1.d7291p-6f, 0x1.835d7cp-2f, 0x1.d9b94cp-1f)},
      {math::Vector(-0x1.55a9acp-2f, 0x1.84c37cp-1f, -0x1.1e0c5p-1f),
       math::Vector(-0x1.8aa97ep-1f, 0x1.f1a1bap-2f, -0x1.a5ca0ep-2f)},
  };
  for (const auto &[start, end] : EDGES) {
    ScratchScope sc(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 2);
    Fragment a, b;
    a.pos = start;
    b.pos = end;
    a.v0 = 0.0f;
    b.v0 = 1.0f;
    points.push_back(a);
    points.push_back(b);

    auto check = [&]<bool SinglePass>() {
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      CapturePipeline full, clipped;
      float terminal_t = 0.0f;
      auto shade = [&](const math::Vector &, Fragment &f) {
        terminal_t = f.v0;
      };
      Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = SinglePass}>(
          full, canvas, points, shade);
      const float FULL_TERMINAL_T = terminal_t;
      Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = SinglePass}>(
          clipped, canvas, points, shade, {.plot_t_start = 0.5f});
      HS_EXPECT_GT(clipped.plotted.size(), size_t{2});
      if (clipped.plotted.size() <= 2 || full.plotted.empty())
        return;
      HS_EXPECT_LT(clipped.plotted.size(), full.plotted.size());
      HS_EXPECT_EQ(terminal_t, FULL_TERMINAL_T);
      HS_EXPECT_EQ(clipped.plotted.back().x, full.plotted.back().x);
      HS_EXPECT_EQ(clipped.plotted.back().y, full.plotted.back().y);
      HS_EXPECT_EQ(clipped.plotted.back().z, full.plotted.back().z);
      CapturePipeline bounded;
      Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = SinglePass}>(
          bounded, canvas, points, shade,
          {.plot_t_start = 0.5f, .plot_t_end = 0.75f});
      HS_EXPECT_GT(bounded.plotted.size(), size_t{1});
      HS_EXPECT_LT(bounded.plotted.size(), clipped.plotted.size());
      HS_EXPECT_LE(terminal_t, 0.75f);
      CapturePipeline near_end;
      terminal_t = -1.f;
      const float BEFORE_END = std::nextafter(1.0f, 0.0f);
      Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = SinglePass}>(
          near_end, canvas, points, shade, {.plot_t_end = BEFORE_END});
      HS_EXPECT_GT(near_end.plotted.size(), size_t{1});
      HS_EXPECT_GT(terminal_t, .75f);
      HS_EXPECT_LE(terminal_t, BEFORE_END);
    };
    check.template operator()<false>();
    check.template operator()<true>();
  }
}

/**
 * @brief The complementary masks of a Segue::Dissolve partition a wireframe's
 *        edges exactly: every edge is drawn by one sprite and skipped by the
 *        other, at every phase.
 * @details Both draws use the same edge list. The shader records each drawn
 *          edge's index (register v2), so the check is on the drawn set itself.
 */
inline void test_mesh_dissolve_masks_partition_edges() {
  constexpr int W = 96, H = 48;
  configure_arenas_default();

  alignas(32) static uint8_t seed_a[24 * 1024];
  alignas(32) static uint8_t seed_b[24 * 1024];
  alignas(32) static uint8_t geom[16 * 1024];
  Arena sa(seed_a, sizeof(seed_a));
  Arena sb(seed_b, sizeof(seed_b));
  Arena ga(geom, sizeof(geom));

  MeshState mesh;
  build_icosahedron_meshstate(sa, sb, ga, mesh);

  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(ga, mesh.faces.size());
  Plot::Mesh::extract_edges(mesh, edges);
  const size_t num_edges = edges.size();

  Segue::Dissolve dissolve;
  hs::random().seed(0xD155);
  dissolve.retarget(math::Y_AXIS);

  auto drawn_set = [&](const DissolveMask &mask) {
    std::vector<bool> seen(num_edges, false);
    auto shade = [&](const math::Vector &, Fragment &f) {
      const int ei = static_cast<int>(f.v2);
      HS_EXPECT_TRUE(ei >= 0 && static_cast<size_t>(ei) < num_edges);
      if (ei < 0 || static_cast<size_t>(ei) >= num_edges)
        return;
      seen[static_cast<size_t>(ei)] = true;
      f.color = Color4(Pixel(65535, 65535, 65535), 0.9f);
    };
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters{
        Filter::Screen::AntiAlias<W, H>()};
    Canvas c(fx);
    Plot::Mesh::draw<W, H>(filters, c, mesh, edges, shade, {}, &mask);
    return seen;
  };

  const float phases[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (uint32_t frame = 0; frame < 3; ++frame) {
    for (float p : phases) {
      const auto masks = dissolve.mask_pair(p, frame);
      auto in_set = drawn_set(masks.incoming);
      auto out_set = drawn_set(masks.outgoing);
      size_t in_count = 0;
      for (size_t e = 0; e < num_edges; ++e) {
        HS_EXPECT_TRUE(in_set[e] != out_set[e]);
        in_count += in_set[e] ? 1 : 0;
      }
      // The endpoints are exact: nothing incoming at phase 0, everything at 1.
      if (p == 0.0f)
        HS_EXPECT_EQ(in_count, size_t{0});
      if (p == 1.0f)
        HS_EXPECT_EQ(in_count, num_edges);
      // Mid-transition both halves must be non-empty, or the split is vacuous.
      if (p == 0.5f) {
        HS_EXPECT_GT(in_count, size_t{0});
        HS_EXPECT_GT(num_edges - in_count, size_t{0});
      }
    }
  }
}

/**
 * @brief End-to-end conservativeness of the rasterizer's column cull: a
 *        quadrant/wedge-clipped render is pixel-identical to the full render
 *        inside the render band, including its margin ring.
 * @details Random trail-like geodesic polylines (the MindSplatter stack) and
 *          planar disk polylines (the ShapeShifter planar-shape stack) through the AntiAlias
 *          pipeline. Clips cover both device quadrants, a narrow interior
 *          wedge, and a seam-adjacent wedge whose margin expansion wraps
 *          (rs > re). A cull false-negative drops in-band pixels and breaks
 *          the comparison.
 */
inline void test_rasterize_column_cull_pixel_parity() {
  constexpr int W = 96, H = 48;
  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.9f);
  };

  hs::random().seed(0xC01C);
  const int clips[4][4] = {
      {0, H / 2, 0, W / 2}, {H / 2, H, W / 2, W}, {0, H, 10, 34}, {0, H, 0, 8}};
  int lit_total = 0, margin_lit_total = 0, row_margin_lit_total = 0;

  for (int trial = 0; trial < 50; ++trial) {
    const bool planar = (trial & 1);
    constexpr size_t WALK = 6;
    math::Vector walk[WALK];
    math::Basis chart;

    if (planar) {
      // Planar polyline: points on a chart disk, as ShapeShifter planar shapes emit them.
      chart = basis_from_normal(rand_unit());
      float radius = hs::rand_f(0.3f, 1.3f);
      float a0 = hs::rand_f(0, 2 * math::PI_F);
      for (size_t i = 0; i < WALK; ++i) {
        float ang = a0 + hs::rand_f(0.2f, 1.0f) * static_cast<float>(i);
        math::Vector dir = chart.u * cosf(ang) + chart.w * sinf(ang);
        walk[i] = (chart.v * cosf(radius) + dir * sinf(radius)).normalized();
      }
    } else {
      // Trail-like random walk: successive short geodesic hops.
      walk[0] = rand_unit();
      for (size_t i = 1; i < WALK; ++i) {
        math::Vector step = math::cross(walk[i - 1], rand_unit());
        if (step.length() < 0.05f) {
          walk[i] = walk[i - 1];
          continue;
        }
        float hop = hs::rand_f(0.15f, 0.7f);
        walk[i] = (walk[i - 1] * cosf(hop) + step.normalized() * sinf(hop))
                      .normalized();
      }
    }

    auto render = [&](hs_test::StubEffect &fx) {
      Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters{
          Filter::Screen::AntiAlias<W, H>()};
      ScratchScope sc(plot_arena());
      Fragments pts;
      pts.bind(plot_arena(), WALK);
      for (const math::Vector &v : walk) {
        Fragment f;
        f.pos = v;
        pts.push_back(f);
      }
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::rasterize<W, H>(
          filters, c, pts, shade,
          {.projection = planar ? Plot::RasterProjection::planar(chart)
                                : Plot::RasterProjection{}});
    };

    std::vector<Pixel> ref(static_cast<size_t>(W) * H);
    {
      hs_test::StubEffect fx(W, H);
      render(fx);
      fx.advance_display();
      hs_test::capture_frame<W, H>(fx, ref);
    }

    for (auto &cl : clips) {
      hs_test::StubEffect fx(W, H);
      fx.set_clip(cl[0], cl[1], cl[2], cl[3]);
      render(fx);
      fx.advance_display();
      const RenderBandDiff diff = render_band_diff<W>(fx, ref);
      lit_total += diff.lit;
      margin_lit_total += diff.margin_lit;
      row_margin_lit_total += diff.row_margin_lit;
      expect_render_band_parity("raster column cull", diff);
    }
  }
  HS_EXPECT_GT(lit_total, 200); // the sweep must actually exercise the bands
  HS_EXPECT_GT(margin_lit_total, 20);
  HS_EXPECT_GT(row_margin_lit_total, 0);
}

/**
 * @brief Pins the whole-trail column cull to the per-edge column bound.
 * @details An edge whose great-circle axis is near-horizontal has no bounded
 *          azimuth span, so the per-edge tier reports it visible. The
 *          whole-trail walk must not cull it from the endpoint columns, which
 *          do not bound it either. The geometry below has |axis.y| ~ 1e-4,
 *          just under geodesic_col_span_cols' threshold, and sits far enough
 *          from the poles to clear the cull's own sin(phi) guard.
 */
inline void test_gate_trail_column_cull_honors_unbounded_edge() {
  constexpr int TW = 288, TH = 144;
  Pipeline<TW, TH, Filter::Screen::AntiAlias<TW, TH>> pipeline{
      Filter::Screen::AntiAlias<TW, TH>()};
  const math::Vector a(-0.616987944f, -0.142148912f, 0.774028182f);
  const math::Vector b(-0.623163402f, 0.021294117f, 0.78180176f);

  HS_EXPECT_FALSE(Plot::make_geodesic_edge_span(a, b).azimuth_bounded);
  size_t excluding_bands = 0;
  ScratchScope sc(plot_arena());
  Fragments trail;
  trail.bind(plot_arena(), 2);
  for (const math::Vector &v : {a, b}) {
    Fragment f;
    f.pos = v;
    trail.push_back(f);
  }

  // Bands across the sphere, so at least one excludes the endpoint columns.
  for (int x0 = 0; x0 < TW; x0 += 24) {
    ClipRegion cr;
    cr.w = TW;
    cr.h = TH;
    cr.y_start = 0;
    cr.y_end = TH;
    cr.x_start = x0;
    cr.x_end = std::min(x0 + 96, TW);
    const auto xc = cr.x_clip();

    const auto A = math::vector_to_pixel<TW, TH>(a);
    const auto B = math::vector_to_pixel<TW, TH>(b);
    excluding_bands += !cr.contains_x(static_cast<int>(A.x)) &&
                       !cr.contains_x(static_cast<int>(B.x));
    uint8_t bits[2];
    const bool any =
        Plot::gate_trail_edges<TW, TH>(pipeline, cr, xc, trail, bits);
    const bool want =
        Plot::edge_visible_in_clip<TW, TH>(pipeline, cr, xc, a, b, nullptr);
    HS_EXPECT_EQ(any, want);
    HS_EXPECT_EQ(bits[0] != 0, want);
  }
  HS_EXPECT_GT(excluding_bands, size_t{0});
}

/**
 * @brief Compares the raw-cross long-edge gate with the normalized span gate.
 */
inline void test_raw_geodesic_edge_gate_parity() {
  constexpr int W = 288, H = 144;
  auto exact_gate = [](const ClipRegion &cr, const ClipRegion::XClip &xc,
                       float ra, float rb, float ca, float cb,
                       const math::Vector &a, const math::Vector &b) {
    const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a, b);
    float row_lo, row_hi;
    Plot::geodesic_row_span_rows<H>(ra, rb, a, b, es, row_lo, row_hi);
    if (!cr.could_intersect_y(row_lo, row_hi + Plot::GEODESIC_ROW_AA_PAD))
      return false;
    if (!xc.active)
      return true;
    int col_s, col_len;
    return !Plot::geodesic_col_span_cols<W>(ca, cb, a, es, col_s, col_len) ||
           ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W);
  };
  auto run = [&](const ClipRegion &cr, const math::Vector &a,
                 const math::Vector &b) {
    const auto xc = cr.x_clip();
    const float ra = Plot::y_to_screen_row<H>(a.y);
    const float rb = Plot::y_to_screen_row<H>(b.y);
    const float ca = math::vector_to_theta<W>(a);
    const float cb = math::vector_to_theta<W>(b);
    const bool exact = exact_gate(cr, xc, ra, rb, ca, cb, a, b);
    const auto raw =
        Plot::raw_geodesic_edge_gate<W, H>(cr, xc, ra, rb, ca, cb, a, b);
    if (raw != Plot::RawGeodesicGateResult::EXACT_FALLBACK)
      HS_EXPECT_EQ(raw == Plot::RawGeodesicGateResult::VISIBLE, exact);
    return raw;
  };

  ClipRegion cr;
  cr.w = W;
  cr.h = H;
  cr.y_start = 0;
  cr.y_end = H / 2;
  cr.x_start = 0;
  cr.x_end = W / 2;

  const math::Vector a = math::X_AXIS;
  auto arc = [&](float angle, const math::Vector &tangent) {
    return (a * cosf(angle) + tangent * sinf(angle)).normalized();
  };
  HS_EXPECT_EQ(run(cr, a, arc(0.0005f, math::Z_AXIS)),
               Plot::RawGeodesicGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(run(cr, a, arc(math::PI_F - 0.0005f, math::Z_AXIS)),
               Plot::RawGeodesicGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(run(cr, a, arc(0.4f, math::Y_AXIS)),
               Plot::RawGeodesicGateResult::EXACT_FALLBACK);

  hs::random().seed(0xC2055);
  const int clips[][4] = {
      {0, H / 2, 0, W / 2}, {0, H / 2, W / 2, W}, {H / 2, H, 0, W / 2},
      {H / 2, H, W / 2, W}, {24, 120, 250, W},    {0, H, 40, 112},
  };
  int raw_count = 0;
  for (const auto &bounds : clips) {
    cr.y_start = bounds[0];
    cr.y_end = bounds[1];
    cr.x_start = bounds[2];
    cr.x_end = bounds[3];
    for (int trial = 0; trial < 5000; ++trial) {
      math::Vector p;
      do {
        const float px = hs::rand_f(-1.0f, 1.0f);
        const float py = hs::rand_f(-1.0f, 1.0f);
        const float pz = hs::rand_f(-1.0f, 1.0f);
        p = math::Vector(px, py, pz);
      } while (p.length() < 0.1f);
      p = p.normalized();
      math::Vector tangent;
      do {
        const float rx = hs::rand_f(-1.0f, 1.0f);
        const float ry = hs::rand_f(-1.0f, 1.0f);
        const float rz = hs::rand_f(-1.0f, 1.0f);
        math::Vector r(rx, ry, rz);
        tangent = r - p * math::dot(r, p);
      } while (tangent.length() < 0.05f);
      tangent = tangent.normalized();
      const float angle = hs::rand_f(0.03f, math::PI_F - 0.03f);
      const math::Vector q =
          (p * cosf(angle) + tangent * sinf(angle)).normalized();
      const auto result = run(cr, p, q);
      if (result != Plot::RawGeodesicGateResult::EXACT_FALLBACK)
        ++raw_count;
    }

    for (int pole = 0; pole < 200; ++pole) {
      const float epsilon = hs::rand_f(1.0e-4f, 0.05f);
      const float azimuth = hs::rand_f(0.0f, 2.0f * math::PI_F);
      const float sp = sinf(epsilon);
      const math::Vector p(sp * cosf(azimuth),
                           (pole & 1) ? cosf(epsilon) : -cosf(epsilon),
                           sp * sinf(azimuth));
      math::Vector tangent = math::cross(p, math::Y_AXIS);
      if (tangent.length() < 1.0e-4f)
        tangent = math::cross(p, math::X_AXIS);
      tangent = tangent.normalized();
      const math::Vector q =
          (p * cosf(0.08f) + tangent * sinf(0.08f)).normalized();
      const auto result = run(cr, p, q);
      if (result != Plot::RawGeodesicGateResult::EXACT_FALLBACK)
        ++raw_count;
    }

    const float guarded_angles[] = {1.0e-6f, 5.0e-4f, 1.5e-3f, 2.5e-3f};
    for (float angle : guarded_angles) {
      for (float end_angle : {angle, math::PI_F - angle}) {
        const auto result = run(cr, a, arc(end_angle, math::Z_AXIS));
        if (result != Plot::RawGeodesicGateResult::EXACT_FALLBACK)
          ++raw_count;
      }
    }
  }
  HS_EXPECT_GT(raw_count, 25000);
}

/**
 * @brief Pins the raw gate's bounded column normalization to the shared span.
 */
inline void test_finish_col_span_one_period() {
  constexpr int W = 288;
  for (int s = 0; s < 4 * W; ++s) {
    const float start = static_cast<float>(s) * 0.25f;
    for (int l = 0; l <= W; l += 3) {
      for (float length :
           {static_cast<float>(l), static_cast<float>(l) + 0.125f}) {
        int raw_s, raw_len, exact_s, exact_len;
        Plot::finish_col_span_one_period<W>(start, length, raw_s, raw_len);
        Plot::finish_col_span<W>(start, length, exact_s, exact_len);
        HS_EXPECT_GE(raw_s, 0);
        HS_EXPECT_LT(raw_s, W);
        HS_EXPECT_EQ(raw_s, exact_s);
        HS_EXPECT_EQ(raw_len, exact_len);
      }
    }
  }
}

/**
 * @brief The conditional-add column wrap is bit-exact with wrap() over (-W, W).
 * @details Probes aggregate into a divergence counter.
 */
inline void test_wrap_one_period_matches_modulo() {
  constexpr int W = 288;
  constexpr float FW = static_cast<float>(W);
  int divergent = 0;
  int probed = 0;
  float first_bad = 2.0f * FW; // outside the domain, so no probe can match it
  auto check = [&](float d) {
    ++probed;
    const float fast = Plot::wrap_one_period<W>(d);
    const float exact = math::wrap(d, FW);
    if (std::bit_cast<uint32_t>(fast) != std::bit_cast<uint32_t>(exact)) {
      if (divergent == 0)
        first_bad = d;
      ++divergent;
    }
  };

  // Deltas as geodesic_col_span_cols forms them: two vector_to_theta results in
  // [0, W) subtracted.
  for (int i = 0; i < 4 * W; ++i) {
    const float a = static_cast<float>(i) * 0.25f;
    for (int j = 0; j < 4 * W; j += 7)
      check(static_cast<float>(j) * 0.25f - a);
  }

  // The seam: magnitudes where d + W rounds up to exactly W, plus both zeros
  // and the extremes of the domain.
  check(0.0f);
  check(-0.0f);
  for (uint32_t bits = 1; bits <= (1u << 20); bits <<= 1) {
    const float tiny = std::bit_cast<float>(bits);
    check(tiny);
    check(-tiny);
  }
  const float top = std::bit_cast<float>(std::bit_cast<uint32_t>(FW) - 1u);
  for (uint32_t back = 0; back < 64; ++back) {
    const float d = std::bit_cast<float>(std::bit_cast<uint32_t>(top) - back);
    check(d);
    check(-d);
  }

  HS_EXPECT_EQ(divergent, 0);
  HS_EXPECT_EQ(first_bad, 2.0f * FW);
  HS_EXPECT_GT(probed, 100000);
}

/**
 * @brief Exercises the Cartesian quadrant gate's rejecting and fallback cases.
 */
inline void test_cartesian_quadrant_gate_classification() {
  constexpr int W = 288, H = 144;
  auto clip = [](int y0, int y1, int x0, int x1) {
    ClipRegion cr;
    cr.w = W;
    cr.h = H;
    cr.y_start = y0;
    cr.y_end = y1;
    cr.x_start = x0;
    cr.x_end = x1;
    return cr;
  };
  auto classify = [](const ClipRegion &cr,
                     std::initializer_list<math::Vector> points) {
    ScratchScope sc(plot_arena());
    Fragments trail;
    trail.bind(plot_arena(), points.size());
    for (const math::Vector &p : points) {
      Fragment f;
      f.pos = p.normalized();
      trail.push_back(f);
    }
    return Plot::cartesian_quadrant_trail_gate(
        Plot::make_cartesian_quadrant_clip<W, H>(cr), trail);
  };

  const ClipRegion north_left = clip(0, H / 2, 0, W / 2);
  const ClipRegion south_right = clip(H / 2, H, W / 2, W);
  HS_EXPECT_EQ(classify(north_left, {math::Vector(0.5f, -0.8f, 0.3f),
                                     math::Vector(0.51f, -0.79f, 0.31f)}),
               Plot::CartesianTrailGateResult::LATITUDE_REJECT);
  HS_EXPECT_EQ(classify(north_left, {math::Vector(0.3f, 0.5f, -0.8f),
                                     math::Vector(0.31f, 0.51f, -0.79f)}),
               Plot::CartesianTrailGateResult::MERIDIAN_REJECT);
  HS_EXPECT_EQ(classify(south_right, {math::Vector(0.5f, 0.8f, 0.3f),
                                      math::Vector(0.51f, 0.79f, 0.31f)}),
               Plot::CartesianTrailGateResult::LATITUDE_REJECT);
  HS_EXPECT_EQ(classify(south_right, {math::Vector(0.3f, -0.5f, 0.8f),
                                      math::Vector(0.31f, -0.51f, 0.79f)}),
               Plot::CartesianTrailGateResult::MERIDIAN_REJECT);

  // Poles and both quadrant boundaries remain exact-fallback cases.
  HS_EXPECT_EQ(classify(north_left, {math::Y_AXIS, math::Y_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(classify(south_right, {-math::Y_AXIS, -math::Y_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(classify(north_left, {math::X_AXIS, math::X_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(classify(north_left, {math::Z_AXIS, math::Z_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);

  // Large and antipodal arcs retain enough slack to fall back; a tiny trail
  // well outside still takes the cheap rejection.
  HS_EXPECT_EQ(classify(north_left, {math::Vector(0.6f, -0.8f, 0.0f),
                                     math::Vector(-0.6f, -0.8f, 0.0f)}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(classify(north_left, {math::X_AXIS, -math::X_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
  HS_EXPECT_EQ(classify(north_left, {math::Vector(0.0f, -1.0f, 0.001f),
                                     math::Vector(0.0f, -1.0f, 0.00101f)}),
               Plot::CartesianTrailGateResult::LATITUDE_REJECT);

  ClipRegion wedge = clip(0, H / 2, 10, 100);
  HS_EXPECT_EQ(classify(wedge, {math::Y_AXIS, math::Y_AXIS}),
               Plot::CartesianTrailGateResult::EXACT_FALLBACK);
}

/**
 * @brief Checks Cartesian rejections against per-edge bounds and dense arc taps.
 * @details Random tiny, ordinary, large, polar, seam, and antipodal edges are
 *          swept over all four hardware quadrants. A Cartesian rejection must
 *          contain no bilinear tap in the render region, including unbounded-pole cases.
 */
inline void test_cartesian_quadrant_gate_is_conservative() {
  constexpr int W = 288, H = 144;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline{
      Filter::Screen::AntiAlias<W, H>()};
  hs::random().seed(0xCA47);

  const int clips[][4] = {
      {0, H / 2, 0, W / 2},
      {0, H / 2, W / 2, W},
      {H / 2, H, 0, W / 2},
      {H / 2, H, W / 2, W},
  };
  int latitude_rejects = 0, meridian_rejects = 0;
  for (const auto &bounds : clips) {
    ClipRegion cr;
    cr.w = W;
    cr.h = H;
    cr.y_start = bounds[0];
    cr.y_end = bounds[1];
    cr.x_start = bounds[2];
    cr.x_end = bounds[3];
    const auto xc = cr.x_clip();
    const auto cartesian = Plot::make_cartesian_quadrant_clip<W, H>(cr);

    for (int trial = 0; trial < 1000; ++trial) {
      ScratchScope sc(plot_arena());
      constexpr size_t N = 6;
      Fragments trail;
      trail.bind(plot_arena(), N);
      math::Vector p =
          trial % 13 == 0
              ? math::Vector(0.001f, trial % 26 == 0 ? 1.0f : -1.0f, 0.001f)
                    .normalized()
              : rand_unit();
      for (size_t k = 0; k < N; ++k) {
        Fragment f;
        f.pos = p;
        trail.push_back(f);
        math::Vector next = rand_unit();
        if (k == 2 && trial % 17 == 0)
          p = -p;
        else {
          const float hop = trial % 5 == 0   ? 1e-5f
                            : trial % 7 == 0 ? 2.7f
                                             : hs::rand_f(0.01f, 0.7f);
          math::Vector tangent = next - p * math::dot(next, p);
          if (tangent.length() > 1e-4f)
            p = (p * cosf(hop) + tangent.normalized() * sinf(hop)).normalized();
        }
      }

      const auto result = Plot::cartesian_quadrant_trail_gate(cartesian, trail);
      if (result == Plot::CartesianTrailGateResult::EXACT_FALLBACK)
        continue;
      if (result == Plot::CartesianTrailGateResult::LATITUDE_REJECT)
        ++latitude_rejects;
      else
        ++meridian_rejects;
      for (size_t e = 0; e + 1 < trail.size(); ++e) {
        const bool visible = Plot::edge_visible_in_clip<W, H>(
            pipeline, cr, xc, trail[e].pos, trail[e + 1].pos, nullptr);
        if (visible) {
          const auto &a = trail[e].pos;
          const auto &b = trail[e + 1].pos;
          const auto span = Plot::make_geodesic_edge_span(a, b);
          const auto tangent = math::cross(span.axis, a);
          for (int sample = 0; sample <= 4096; ++sample) {
            const float angle = span.total * (sample / 4096.0f);
            const auto p = a * cosf(angle) + tangent * sinf(angle);
            const auto screen = math::vector_to_pixel<W, H>(p);
            const int x = static_cast<int>(floorf(screen.x));
            const int y = static_cast<int>(floorf(screen.y));
            for (int dy = 0; dy <= 1; ++dy)
              for (int dx = 0; dx <= 1; ++dx) {
                const int wx = (x + dx + W) % W;
                HS_EXPECT_FALSE(cr.contains_y(y + dy) && cr.contains_x(wx));
              }
          }
        }
      }
    }
  }
  HS_EXPECT_GT(latitude_rejects, 100);
  HS_EXPECT_GT(meridian_rejects, 100);
}

/**
 * @brief Pins gate_trail_edges to the per-edge edge_visible_in_clip verdicts.
 * @details Random geodesic step-walk trails over the device band shapes. A
 *          false return must leave every byte zero AND every edge individually
 *          invisible (the whole-trail bound is conservative); a true return's
 *          bytes must equal the per-edge predicate exactly (rasterize consumes
 *          them as its cull).
 */
inline void test_gate_trail_edges_matches_edge_visible() {
  constexpr int TW = 288, TH = 144;
  Pipeline<TW, TH, Filter::Screen::AntiAlias<TW, TH>> pipeline{
      Filter::Screen::AntiAlias<TW, TH>()};
  hs::random().seed(0x60FE);

  const int bands[][4] = {
      {0, 72, 0, 144},
      {36, 108, 144, 288},
      {0, 144, 60, 200},
      {100, 144, 0, 288},
  };
  int rejects = 0, visible = 0, culled = 0;
  for (const auto &bd : bands) {
    ClipRegion cr;
    cr.w = TW;
    cr.h = TH;
    cr.y_start = bd[0];
    cr.y_end = bd[1];
    cr.x_start = bd[2];
    cr.x_end = bd[3];
    const auto xc = cr.x_clip();

    for (int trial = 0; trial < 500; ++trial) {
      ScratchScope sc(plot_arena());
      const size_t n = 2 + static_cast<size_t>(hs::rand_f(0.0f, 38.0f));
      Fragments trail;
      trail.bind(plot_arena(), n);
      math::Vector p = rand_unit();
      for (size_t k = 0; k < n; ++k) {
        Fragment f;
        f.pos = p;
        trail.push_back(f);
        math::Vector step = math::cross(p, rand_unit());
        if (step.length() < 0.05f)
          continue;
        const float hop = hs::rand_f(0.005f, 0.4f);
        p = (p * cosf(hop) + step.normalized() * sinf(hop)).normalized();
      }

      uint8_t bits[40];
      const bool any =
          Plot::gate_trail_edges<TW, TH>(pipeline, cr, xc, trail, bits);
      for (size_t e = 0; e + 1 < n; ++e) {
        const bool want = Plot::edge_visible_in_clip<TW, TH>(
            pipeline, cr, xc, trail[e].pos, trail[e + 1].pos, nullptr);
        if (any) {
          HS_EXPECT_TRUE((bits[e] != 0) == want);
        } else {
          HS_EXPECT_EQ(bits[e], 0);
          HS_EXPECT_FALSE(want);
        }
        (want ? visible : culled)++;
      }
      if (!any)
        ++rejects;
    }
  }
  // All three outcomes must be exercised: whole-trail rejects, per-edge
  // culls, and visible edges.
  HS_EXPECT_GT(rejects, 20);
  HS_EXPECT_GT(visible, 1000);
  HS_EXPECT_GT(culled, 1000);
}

/**
 * @brief Verifies visible arc samples belong to pieces the clip gate keeps.
 * @details Conservativeness proof for the clip cut: sweeps the rendered great circle of
 *          random edges against bands covering both seam topologies: a sample
 *          whose plotted pixel falls in the render region must belong to a kept
 *          piece. Kept and culled piece counts are floored too, so a cut that
 *          stops separating the band (or stops cutting at all) cannot pass by
 *          keeping everything.
 */
inline void test_mesh_clip_cut_separates_band() {
  constexpr int TW = 288, TH = 144;
  constexpr int SWEEP = 128;
  Pipeline<TW, TH, Filter::Screen::AntiAlias<TW, TH>> pipeline{
      Filter::Screen::AntiAlias<TW, TH>()};
  hs::random().seed(0xC07);

  const int bands[][4] = {
      {0, 72, 0, 144},
      {36, 108, 144, 288},
      {0, 144, 60, 200},
      {40, 100, 2, 30},
  };
  int cuts = 0, kept = 0, culled = 0;
  long shown_arc = 0, drawn_arc = 0;
  for (const auto &bd : bands) {
    ClipRegion cr;
    cr.w = TW;
    cr.h = TH;
    cr.y_start = bd[0];
    cr.y_end = bd[1];
    cr.x_start = bd[2];
    cr.x_end = bd[3];
    const auto xc = cr.x_clip();
    const Plot::ClipCutBounds cb = Plot::make_clip_cut_bounds<TW, TH>(cr, xc);

    for (int trial = 0; trial < 300; ++trial) {
      ScratchScope sc(plot_arena());
      Fragment fa, fb;
      fa.pos = rand_unit();
      fb.pos = rand_unit();
      const Plot::GeodesicEdgeSpan es =
          Plot::make_geodesic_edge_span(fa.pos, fb.pos);
      if (!es.have_axis)
        continue;

      float ts[Plot::GEODESIC_CLIP_MAX_SPLITS];
      const int n = Plot::geodesic_clip_splits(fa.pos, fb.pos, es, cb, ts);
      cuts += n;

      const math::Vector perp = math::cross(es.axis, fa.pos);
      Fragments points;
      points.bind(plot_arena(), Plot::Mesh::EDGE_MAX_POINTS);
      points.push_back(Plot::Line::sample_point(fa, fb, es, perp, 0.0f));
      for (int i = 0; i < n; ++i)
        points.push_back(Plot::Line::sample_point(fa, fb, es, perp, ts[i]));
      points.push_back(Plot::Line::sample_point(fa, fb, es, perp, 1.0f));

      uint8_t bits[Plot::Mesh::EDGE_MAX_POINTS - 1];
      Plot::gate_trail_edges<TW, TH>(pipeline, cr, xc, points, bits);
      for (int i = 0; i <= n; ++i)
        (bits[i] != 0 ? kept : culled)++;

      for (int s = 0; s <= SWEEP; ++s) {
        const float t = static_cast<float>(s) / SWEEP;
        const math::Vector p =
            Plot::Line::sample_point(fa, fb, es, perp, t).pos;
        const math::PixelCoords px = math::vector_to_pixel<TW, TH>(p);
        const bool shown = cr.contains_x(static_cast<int>(px.x)) &&
                           cr.contains_y(static_cast<int>(px.y));
        int piece = 0;
        while (piece < n && t > ts[piece])
          ++piece;
        const bool drawn = bits[piece] != 0;
        if (shown)
          HS_EXPECT_TRUE(drawn);
        shown_arc += shown ? 1 : 0;
        drawn_arc += drawn ? 1 : 0;
      }
    }
  }
  HS_EXPECT_GT(cuts, 400);
  HS_EXPECT_GT(kept, 400);
  HS_EXPECT_GT(culled, 400);
  // Tightness, the half the gate cannot supply on its own: it keeps any piece
  // straddling the band edge, so a cut that stops separating shows only as
  // drawn arc past what the region shows. Measured 1.08x; dropping the cut
  // entirely gives 1.96x and cutting inside the band 1.76x.
  HS_EXPECT_LT(drawn_arc * 2, shown_arc * 3);
}

/**
 * @brief Pins rasterize's precomputed-bits path to its inline gate: rendering
 *        with gate_trail_edges bytes must be pixel-identical to rendering the
 *        same polyline with the per-edge cull evaluated in place.
 */
inline void test_rasterize_gate_bits_pixel_parity() {
  constexpr int W = 96, H = 48;
  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.9f);
  };
  hs::random().seed(0x617E);

  const int clips[3][4] = {
      {0, H / 2, 0, W / 2}, {H / 2, H, W / 2, W}, {0, H, 10, 34}};
  int lit_total = 0;

  for (int trial = 0; trial < 40; ++trial) {
    constexpr size_t WALK = 8;
    math::Vector walk[WALK];
    walk[0] = rand_unit();
    for (size_t i = 1; i < WALK; ++i) {
      math::Vector step = math::cross(walk[i - 1], rand_unit());
      if (step.length() < 0.05f) {
        walk[i] = walk[i - 1];
        continue;
      }
      const float hop = hs::rand_f(0.05f, 0.6f);
      walk[i] = (walk[i - 1] * cosf(hop) + step.normalized() * sinf(hop))
                    .normalized();
    }

    for (auto &cl : clips) {
      auto render = [&](hs_test::StubEffect &fx, bool use_bits) {
        Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters{
            Filter::Screen::AntiAlias<W, H>()};
        ScratchScope sc(plot_arena());
        Fragments pts;
        pts.bind(plot_arena(), WALK);
        for (const math::Vector &v : walk) {
          Fragment f;
          f.pos = v;
          pts.push_back(f);
        }
        Canvas c(fx);
        uint8_t bits[WALK - 1];
        std::span<const uint8_t> vis;
        if (use_bits) {
          const ClipRegion &cr = fx.clip();
          const auto xc = cr.x_clip();
          // A whole-trail reject renders nothing; the inline-gate reference
          // paints nothing for it too (per-edge conservativeness), so the
          // buffers still compare equal.
          if (!Plot::gate_trail_edges<W, H>(filters, cr, xc, pts, bits))
            return;
          vis = {bits, pts.size() - 1};
        }
        Plot::rasterize<W, H>(
            filters, c, pts, shade,
            {.projection = Plot::RasterProjection::geodesic(vis)});
      };

      std::vector<Pixel> ref(static_cast<size_t>(W) * H);
      {
        hs_test::StubEffect fx(W, H);
        fx.set_clip(cl[0], cl[1], cl[2], cl[3]);
        render(fx, false);
        fx.advance_display();
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            const Pixel p = fx.get_pixel(x, y);
            ref[static_cast<size_t>(y) * W + x] = p;
            if (p.r | p.g | p.b)
              ++lit_total;
          }
      }
      {
        hs_test::StubEffect fx(W, H);
        fx.set_clip(cl[0], cl[1], cl[2], cl[3]);
        render(fx, true);
        fx.advance_display();
        int diff = 0;
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            const Pixel p = fx.get_pixel(x, y);
            const Pixel &r = ref[static_cast<size_t>(y) * W + x];
            if (p.r != r.r || p.g != r.g || p.b != r.b)
              ++diff;
          }
        HS_EXPECT_EQ(diff, 0);
      }
    }
  }
  HS_EXPECT_GT(lit_total, 200); // the sweep must actually light the bands
}
