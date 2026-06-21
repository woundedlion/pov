/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the core/scan.h rasterizer. Exercises two entry points end to
 * end into a live Canvas:
 *   - Scan::Shader::draw  : full-sphere per-pixel shader (constant, positional,
 *                           and clip-respecting).
 *   - Scan::Ring::draw    : the SDF rasterize() path (bounding-box scan ->
 *                           interval generation -> process_pixel -> Pipeline
 *                           sink plot), verified to produce bounded output.
 */
#pragma once

#include "core/scan.h"
#include "core/plot.h"
#include "core/filter.h"
#include "core/canvas.h"
#include "core/geometry.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace scan_tests {

/**
 * @brief Minimal Effect backing a Canvas for these tests.
 * @details Performs no per-frame drawing and shows no background, so the canvas
 * starts black and displays only what the test explicitly plots.
 */
struct ScanFx : public Effect {
  /**
   * @brief Constructs the test effect at the given canvas resolution.
   * @param W Canvas width in pixels.
   * @param H Canvas height in pixels.
   */
  ScanFx(int W, int H) : Effect(W, H) {}
  /**
   * @brief Per-frame draw hook; intentionally does nothing for these tests.
   */
  void draw_frame() override {}
  /**
   * @brief Reports whether the effect renders a background.
   * @return Always false, so the canvas starts black.
   */
  bool show_bg() const override { return false; }
};

/**
 * @brief Tests whether a pixel is fully unwritten (cleared-frame black).
 * @param p Pixel to inspect.
 * @return True when all RGB channels are zero.
 */
inline bool is_black(const Pixel &p) { return p.r == 0 && p.g == 0 && p.b == 0; }

// ============================================================================
// Scan::Shader::draw — full-sphere per-pixel shader
// ============================================================================

/**
 * @brief Verifies a constant-color shader fills every pixel of the full sphere.
 */
inline void test_shader_constant_fills_canvas() {
  constexpr int W = 32, H = 16;
  ScanFx fx(W, H);
  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(
        c, [](const Vector &) { return Color4(Pixel(40000, 20000, 10000), 1.0f); });
  }
  fx.advance_display();

  // Every pixel in the full canvas got the constant color.
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
 * @details On a partial-coverage pixel whose sub-samples vary in BOTH color and
 * alpha (here: two opaque red, two transparent black per pixel), correct
 * premultiplied SSAA writes (sum of color*alpha) / N. The old straight-alpha
 * model — average color and alpha separately, then re-multiply — would apply
 * coverage twice and darken the result (red/4 instead of red/2). This pins the
 * premultiplied result.
 */
inline void test_shader_ssaa_premultiplies_partial_coverage() {
  constexpr int W = 16, H = 8;
  ScanFx fx(W, H);
  {
    Canvas c(fx);
    // The 2x2 grid calls the shader 4× per pixel; alternating opaque/transparent
    // by call index yields exactly two opaque red + two transparent samples per
    // pixel (4 is even, so the phase realigns at every pixel).
    int call = 0;
    Scan::Shader::draw<W, H, 4>(c, [&call](const Vector &) -> Color4 {
      bool opaque = (call++ % 2) == 0;
      return opaque ? Color4(Pixel(60000, 0, 0), 1.0f)
                    : Color4(Pixel(0, 0, 0), 0.0f);
    });
  }
  fx.advance_display();

  // Premultiplied: (60000*1 + 60000*1 + 0 + 0) / 4 = 30000.
  // Straight-alpha (old, wrong) would have been (60000/2) * 0.5 = 15000.
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
 * @brief Verifies a position-reading shader maps the sphere's latitude.
 * @details The +Y pole (top row) must render brighter than the -Y pole (bottom
 * row), confirming the shader receives the correct surface position.
 */
inline void test_shader_positional_maps_latitude() {
  constexpr int W = 32, H = 32;
  ScanFx fx(W, H);
  {
    Canvas c(fx);
    // Green encodes latitude: north pole (v.y≈+1) bright, south (v.y≈-1) dark.
    Scan::Shader::draw<W, H, 1>(c, [](const Vector &v) {
      uint16_t g = (uint16_t)((v.y * 0.5f + 0.5f) * 60000.0f);
      return Color4(Pixel(0, g, 0), 1.0f);
    });
  }
  fx.advance_display();

  // Top row (y=0) is the north pole and must be brighter than the bottom row.
  HS_EXPECT_GT((int)fx.get_pixel(0, 0).g, (int)fx.get_pixel(0, H - 1).g);
  HS_EXPECT_GT((int)fx.get_pixel(0, 0).g, 40000);
  HS_EXPECT_LT((int)fx.get_pixel(0, H - 1).g, 20000);
}

/**
 * @brief Verifies the shader writes only inside the active clip band.
 * @details Rows outside the clip band must stay black (untouched by the clear).
 */
inline void test_shader_respects_clip_band() {
  constexpr int W = 32, H = 16;
  ScanFx fx(W, H);
  fx.set_clip(5, 10, 0, W); // rows [5,10)
  fx.clip.margin = 0;       // no render-margin expansion

  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(
        c, [](const Vector &) { return Color4(Pixel(0, 0, 50000), 1.0f); });
  }
  fx.advance_display();

  // Inside the band: written. Outside: untouched (black from the frame clear).
  HS_EXPECT_FALSE(is_black(fx.get_pixel(0, 7)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(0, 0)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(0, 12)));
}

// ============================================================================
// Scan::Ring::draw — SDF rasterize() path through a Pipeline sink
// ============================================================================

/**
 * @brief Verifies the SDF rasterize() path plots a bounded ring.
 * @details The ring must cover a nonempty subset of the canvas, but never the
 * whole canvas.
 */
inline void test_ring_rasterize_produces_bounded_output() {
  constexpr int W = 64, H = 48;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe; // bare 2D sink (no filters)

  size_t plotted = 0;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Ring::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.5f, /*thickness=*/0.4f,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++plotted;

  // The rasterizer plotted a ring: some pixels, but not the whole canvas.
  HS_EXPECT_GT(plotted, (size_t)0);
  HS_EXPECT_LT(plotted, (size_t)(W * H));
}

/**
 * @brief Geometric oracle: every lit pixel of a rasterized ring lies on the
 *        ring's latitude band — placement, not just "something was drawn".
 * @details The ring axis is +Y, so a pixel's world direction v sits at polar
 * angle acos(v.y) from the pole. A correct stroke lights only pixels whose polar
 * angle is within `thickness` (plus the AA fringe and one row of pixel
 * quantization) of the ring centre `target = radius·π/2`. A pixel far off the
 * band — near the pole (a fill that spilled inward) or in the wrong hemisphere
 * (a projection sign error) — fails this, where the bounded-count check above
 * would not. The interior cap (polar < target − band) staying dark proves the
 * ring is a hollow stroke, not a filled disk.
 */
inline void test_ring_rasterize_lit_pixels_on_band() {
  constexpr int W = 96, H = 64;
  constexpr float radius = 0.5f, thickness = 0.2f;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Ring::draw<W, H, false>(
        pipe, c, basis, radius, thickness,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  const float target = radius * (PI_F / 2.0f); // ring-centre polar angle
  const float row = PI_F / (H - 1);            // angular height of one row
  const float band = thickness + 2.0f * row;   // stroke + AA + quantization
  size_t lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      if (is_black(fx.get_pixel(x, y)))
        continue;
      ++lit;
      Vector v = pixel_to_vector<W, H>(x, y);
      float polar = acosf(hs::clamp(v.y, -1.0f, 1.0f));
      HS_EXPECT_LE(fabsf(polar - target), band);
    }
  HS_EXPECT_GT(lit, (size_t)0); // the ring is actually drawn
}

/**
 * @brief Verifies a degenerate clip band makes rasterize plot nothing.
 * @details With y_start == y_end the rasterizer must early-out and leave the
 * canvas black.
 */
inline void test_ring_rasterize_empty_clip_draws_nothing() {
  constexpr int W = 64, H = 48;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;

  // Degenerate clip (y_start > y_end) → rasterize must early-out, plot nothing.
  fx.set_clip(30, 30, 0, W);
  fx.clip.margin = 0;

  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Ring::draw<W, H, false>(
        pipe, c, basis, 0.5f, 0.4f, [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  size_t plotted = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++plotted;
  HS_EXPECT_EQ(plotted, (size_t)0);
}

/**
 * @brief Verifies the scan_region seam coalescer avoids double-plotting.
 * @details A span crossing x=0 must not double-plot the wrapped overlap shared
 * with another span. Drives scan_region with a sorted two-span row (a low span
 * plus a seam-crosser, as a shape's sorted output emits) so the wrapped columns
 * shared by both spans are plotted exactly once.
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
      [&](int wx, int, const Vector &) {
        if (wx >= 0 && wx < W)
          counts[wx]++;
      });

  // No pixel plotted more than once (the wrapped overlap at x=0,1 is not doubled).
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
 * one pixel column must not both plot that column. Float span merging alone
 * would let prev end 5.4 (ceil paints x=5) and next start 5.6 (floor gives 5)
 * each touch x=5, doubling process_pixel / alpha; integer-space clamping
 * (last_x2) keeps x=5 single.
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
      [&](int wx, int, const Vector &) {
        if (wx >= 0 && wx < W)
          counts[wx]++;
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
 * @brief Verifies a geodesic line through the north pole plots the pole row.
 * @details map_geodesic/map_planar build interpolated points with
 * fast_sinf/fast_cosf, which are ~0.04% non-unit; vector_to_pixel takes
 * phi = acos(v.y) directly, and acos's infinite slope at y=1 amplifies that
 * tiny error into a multi-row shift unless interpolated positions are
 * re-normalized before mapping. The drawing phase re-normalizes, so the pole
 * lands on row 0.
 */
inline void test_plot_line_over_pole_reaches_row0() {
  constexpr int W = 288, H = 144;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe; // bare sink (no AA) so we see raw sample placement

  // Geodesic from 0.4 rad down the +Z side of the N pole to 0.4 rad down the
  // -Z side; its midpoint is the pole (row 0).
  Fragment f1, f2;
  f1.pos = Vector(0.0f, cosf(0.4f), sinf(0.4f));
  f2.pos = Vector(0.0f, cosf(0.4f), -sinf(0.4f));
  {
    Canvas c(fx);
    Plot::Line::draw<W, H>(pipe, c, f1, f2, [](const Vector &, Fragment &f) {
      f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
    });
  }
  fx.advance_display();

  size_t row0 = 0;
  for (int x = 0; x < W; ++x)
    if (!is_black(fx.get_pixel(x, 0)))
      ++row0;
  // The line reaches the pole row.
  HS_EXPECT_GT(row0, (size_t)0);
}

/**
 * @brief Capturing plot sink that records the AA alpha process_pixel forwards.
 * @details Bypasses the canvas-blend round trip. The recorded alpha is
 * frag.alpha (=1) times the AA alpha; count tracks whether the pixel was drawn
 * at all.
 */
struct AlphaSink {
  float last_alpha = -1.0f; /**< AA alpha from the most recent plot, -1 if none. */
  int count = 0;            /**< Number of times plot() was invoked. */
  /**
   * @brief Records the forwarded AA alpha and increments the plot count.
   * @param a Anti-aliasing alpha forwarded by process_pixel.
   */
  void plot(Canvas &, int, int, const Pixel &, float, float a) {
    last_alpha = a;
    ++count;
  }
};

/**
 * @brief Runs process_pixel for a shape at a surface point and returns its AA alpha.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param shape SDF shape to rasterize at the sample point.
 * @param p Surface point on the unit sphere to evaluate.
 * @param c Canvas providing clip and projection context.
 * @param count Optional out-param; receives how many times the pixel was plotted (0 or 1).
 * @return The forwarded AA alpha, or -1 if the pixel was not drawn.
 */
template <int W, int H>
inline float scan_alpha_at(const auto &shape, const Vector &p, Canvas &c,
                           int *count = nullptr) {
  AlphaSink sink;
  SDF::DistanceResult res;
  Fragment frag;
  Scan::process_pixel<W, H, false>(
      0, 0, p, sink, c, shape,
      [](const Vector &, Fragment &f) {
        f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
      },
      /*debug_bb=*/false, res, frag);
  if (count)
    *count = sink.count;
  return sink.last_alpha;
}

/**
 * @brief Verifies a CSG composite renders each stroke with its own thickness.
 * @details Each stroke must render with its OWN thickness, not as a hard solid
 * band and not scaled by a sibling's thickness. The formula-independent
 * invariant: a Union<thin, thick> evaluated at a point inside the thin stroke
 * yields the SAME AA alpha as the bare thin Line, and a DIFFERENT alpha from
 * the bare thick Line.
 */
inline void test_csg_stroke_aa_uses_winning_child_thickness() {
  constexpr int W = 288, H = 144;
  ScanFx fx(W, H);
  Canvas c(fx);

  const float thin = 0.05f, thick = 0.30f;
  // Thin line: equatorial arc +X -> +Z (great circle in the y=0 plane).
  SDF::Line thin_line(Vector(1, 0, 0), Vector(0, 0, 1), thin);
  // Thick line: opposite equatorial quadrant, far from the test point so the
  // Union always selects the thin line; it exists only to push the wrapper's
  // max-thickness up to `thick`.
  SDF::Line thick_line(Vector(-1, 0, 0), Vector(0, 0, -1), thick);
  SDF::Union<SDF::Line, SDF::Line> u(thin_line, thick_line);
  // Same geometry as `thick_line` but standalone, for the contrast check.
  SDF::Line thick_solo(Vector(1, 0, 0), Vector(0, 0, 1), thick);

  // A composite of two strokes is itself a stroke.
  HS_EXPECT_FALSE((SDF::Union<SDF::Line, SDF::Line>::is_solid));

  // Point at geodesic distance 0.025 (half the thin thickness) north of the
  // thin arc, projecting to azimuth 45 deg (well inside the arc).
  const float dist = 0.025f;
  Vector p(cosf(dist) * cosf(PI_F / 4), sinf(dist),
           cosf(dist) * sinf(PI_F / 4));

  int n_union = 0, n_bare = 0;
  float a_union = scan_alpha_at<W, H>(u, p, c, &n_union);
  float a_thin = scan_alpha_at<W, H>(thin_line, p, c, &n_bare);
  float a_thick = scan_alpha_at<W, H>(thick_solo, p, c);

  HS_EXPECT_EQ(n_union, 1);
  HS_EXPECT_EQ(n_bare, 1);
  // The composite reproduces the winning child's own AA exactly...
  HS_EXPECT_NEAR(a_union, a_thin, 1e-4f);
  // ...and the thin/thick falloffs are clearly separable, so reproducing thin
  // rather than the wrapper-max thick is a meaningful distinction.
  HS_EXPECT_GT(fabsf(a_thin - a_thick), 0.1f);
}

/**
 * @brief Verifies the rasterized ring lights the analytically predicted row.
 * @details A ring of normalized radius r is centered on the basis axis at polar
 * angle target = r*(PI/2), lighting a single latitude band whose center row is
 * phi_to_y(target). Asserts the rasterizer output lands there (an analytic
 * position check, not just "some pixels, not all") and that rows well away from
 * the band stay dark.
 */
inline void test_ring_rasterize_lights_expected_row() {
  constexpr int W = 96, H = 48;

  auto centroid_and_band = [](float radius) {
    ScanFx fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Basis basis = make_basis(Quaternion(), Y_AXIS); // axis = north pole (+Y)
      Scan::Ring::draw<W, H, false>(
          pipe, c, basis, radius, /*thickness=*/0.05f,
          [](const Vector &, Fragment &f) {
            f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
          });
    }
    fx.advance_display();

    int lit[H] = {0};
    long total = 0, weighted = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        if (!is_black(fx.get_pixel(x, y))) {
          lit[y]++;
          total++;
          weighted += y;
        }
    /** @brief Per-radius result: lit-pixel centroid row, total lit count, and per-row lit counts. */
    struct R { float centroid; int total; int lit[H]; };
    R r;
    r.centroid = total ? static_cast<float>(weighted) / total : -1.0f;
    r.total = static_cast<int>(total);
    for (int y = 0; y < H; ++y) r.lit[y] = lit[y];
    return r;
  };

  for (float radius : {0.5f, 1.0f}) {
    float target = radius * (PI_F / 2.0f);
    float expected_y = phi_to_y<H>(target);
    auto r = centroid_and_band(radius);

    // The ring is a full circle of latitude: most columns light up its row.
    HS_EXPECT_GT(r.total, W / 2);
    // Lit-pixel centroid lands on the analytically predicted row.
    HS_EXPECT_NEAR(r.centroid, expected_y, 1.0f);
    // Rows far from the band (> 4 px away) are dark.
    int ey = static_cast<int>(expected_y + 0.5f);
    for (int y = 0; y < H; ++y)
      if (std::abs(y - ey) > 4)
        HS_EXPECT_EQ(r.lit[y], 0);
  }
}

/**
 * @brief Verifies the stroke anti-aliasing alpha is a monotone ramp.
 * @details The alpha falls from ~1 at the ring centerline to 0 at its outer
 * surface, not a hard binary edge. Samples the AA alpha process_pixel produces
 * while marching a point radially outward across the band and asserts it
 * decreases monotonically through intermediate values.
 */
inline void test_stroke_aa_is_monotone_ramp() {
  constexpr int W = 288, H = 144;
  ScanFx fx(W, H);
  Canvas c(fx);

  const float radius = 0.5f;            // centerline at polar PI/4
  const float thickness = 0.10f;        // band half-width in radians
  const float target = radius * (PI_F / 2.0f);
  Basis basis = make_basis(Quaternion(), Y_AXIS);
  SDF::Ring ring(basis, radius, thickness);

  // March outward from the centerline along the az=0 meridian. raw distance to
  // the centerline grows with the offset, so the signed distance crosses 0 at
  // the surface and the alpha should fall 1 -> 0.
  const int N = 12;
  float prev = 2.0f;
  bool saw_one = false, saw_zero = false, saw_mid = false;
  for (int i = 0; i < N; ++i) {
    float delta = (thickness * 1.4f) * i / (N - 1); // 0 .. 1.4*thickness
    float ph = target + delta;
    Vector p(sinf(ph), cosf(ph), 0.0f); // az=0, polar angle ph
    int count = 0;
    float a = scan_alpha_at<W, H>(ring, p, c, &count);
    if (count == 0) a = 0.0f; // outside the stroke -> not drawn -> alpha 0

    // Monotone non-increasing as we move outward.
    HS_EXPECT_LE(a, prev + 1e-4f);
    prev = a;

    if (a > 0.95f) saw_one = true;
    if (a < 0.05f) saw_zero = true;
    if (a > 0.1f && a < 0.9f) saw_mid = true;
  }

  // Centerline ~opaque, far edge ~transparent, with a genuine ramp between
  // (a hard edge would jump 1 -> 0 with no intermediate sample).
  HS_EXPECT_TRUE(saw_one);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_TRUE(saw_mid);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every scan test under the "scan" module scope.
 * @return The number of test failures recorded by the module.
 */
inline int run_scan_tests() {
  auto scope = hs_test::begin_module("scan");

  test_shader_constant_fills_canvas();
  test_shader_ssaa_premultiplies_partial_coverage();
  test_shader_positional_maps_latitude();
  test_shader_respects_clip_band();
  test_ring_rasterize_produces_bounded_output();
  test_ring_rasterize_lit_pixels_on_band();
  test_ring_rasterize_lights_expected_row();
  test_stroke_aa_is_monotone_ramp();
  test_ring_rasterize_empty_clip_draws_nothing();
  test_scan_region_seam_no_double_plot();
  test_scan_region_fractional_boundary_no_double_plot();
  test_plot_line_over_pole_reaches_row0();
  test_csg_stroke_aa_uses_winning_child_thickness();

  return hs_test::end_module(scope);
}

} // namespace scan_tests
} // namespace hs_test

