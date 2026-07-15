/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the core/render/scan.h rasterizer. Exercises two entry points end to
 * end into a live Canvas:
 *   - Scan::Shader::draw  : full-sphere per-pixel shader (constant, positional,
 *                           and clip-respecting).
 *   - Scan::Ring::draw    : the SDF rasterize() path (bounding-box scan ->
 *                           interval generation -> process_pixel -> Pipeline
 *                           sink plot), verified to produce bounded output.
 */
#pragma once

#include "core/render/scan.h"
#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/canvas.h"
#include "core/math/geometry.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <vector>

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

  // Bit-exact readback: at alpha=1 the blend is the identity (src*1 + dst*0), so
  // channels survive verbatim.
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
    // Opacity keys on the sub-sample's position, not call order: the 2x2 grid's
    // +0.25/-0.25 px x-offsets land at fractional theta-grid phase 0.25 vs 0.75,
    // so two of the four samples per pixel are opaque regardless of iteration.
    Scan::Shader::draw<W, H, 4>(c, [](const Vector &v) -> Color4 {
      float theta = std::atan2(v.z, v.x);
      if (theta < 0.0f) theta += 2.0f * PI_F;
      float g = theta * W / (2.0f * PI_F);
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

  // Top row (y=0) is the north pole, brighter than the bottom row.
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
  fx.set_margin(0);       // no render-margin expansion

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
  HS_EXPECT_GT(lit, (size_t)0);
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
  fx.set_margin(0);

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
 * @brief Verifies flat-ring rasterization matches zero knots within one code.
 */
inline void test_distorted_ring_flat_matches_zero_knot_raster() {
  constexpr int W = 96, H = 64, LUT_N = 16;
  float knots[LUT_N + 1] = {};
  Basis basis =
      make_basis(Quaternion(), Vector(0.3f, 0.8f, -0.5f).normalized());

  auto check = [&](bool partial_clip) {
    auto shader = [](const Vector &, Fragment &f) {
      f.color = Color4(Pixel(51000, 27000, 9000), 0.8f * f.v2);
    };
    std::vector<Pixel> expected(W * H);
    {
      ScanFx legacy(W, H);
      if (partial_clip) {
        legacy.set_clip(9, 53, 17, 81);
        legacy.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(legacy);
        Scan::DistortedRing::draw<W, H>(pipeline, canvas, basis, 1.7f, 0.12f,
                                        knots, LUT_N, 0.0f, shader);
      }
      legacy.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          expected[y * W + x] = legacy.get_pixel(x, y);
    }

    ScanFx flat(W, H);
    if (partial_clip) {
      flat.set_clip(9, 53, 17, 81);
      flat.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(flat);
      Scan::DistortedRing::draw_flat<W, H, false>(pipeline, canvas, basis, 1.7f,
                                                  0.12f, shader);
    }
    flat.advance_display();

    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = flat.get_pixel(x, y);
        HS_EXPECT_NEAR(static_cast<int>(a.r), static_cast<int>(b.r), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.g), static_cast<int>(b.g), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.b), static_cast<int>(b.b), 1);
      }
    }
  };

  check(false);
  check(true);
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
        [&](int wx, int, const Vector &) {
          if (wx >= 0 && wx < W)
            counts[wx]++;
        },
        xc);
  };

  const ClipRegion::XClip no_clip{};
  const ClipRegion::XClip arc{3, 8, true, false};
  const ClipRegion::XClip wrap_arc{W - 2, 2, true, true};

  for (bool handled : {true, false}) {
    int reference[W];
    run(no_clip, handled, reference);
    for (const auto &xc : {arc, wrap_arc}) {
      int counts[W];
      run(xc, handled, counts);
      for (int x = 0; x < W; ++x) {
        HS_EXPECT_LE(counts[x], 1);
        HS_EXPECT_EQ(counts[x], xc.clipped(x) ? 0 : reference[x]);
      }
    }
  }
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
 * @brief Verifies Scan's v2 contract for stroke and solid shaders.
 * @details A stroke shader that writes v2 as alpha produces coverage squared
 * at the sink; solid shaders receive zero.
 */
inline void test_scan_shader_v2_contract() {
  constexpr int W = 288, H = 144;
  ScanFx fx(W, H);
  Canvas canvas(fx);
  Basis basis = make_basis(Quaternion(), Y_AXIS);
  constexpr float RADIUS = 0.5f;
  constexpr float THICKNESS = 0.1f;
  const float phi = RADIUS * (PI_F / 2.0f) + THICKNESS * 0.5f;
  const Vector p(sinf(phi), cosf(phi), 0.0f);
  SDF::Ring ring(basis, RADIUS, THICKNESS);

  SDF::DistanceResult expected_result;
  ring.distance<false>(p, expected_result);
  const float expected_coverage =
      quintic_kernel(-expected_result.dist / expected_result.size);
  const float legacy_factor = quintic_kernel(
      1.0f - hs::clamp(expected_result.raw_dist / expected_result.size, 0.0f,
                       1.0f));

  AlphaSink sink;
  SDF::DistanceResult result;
  Fragment frag;
  float shader_coverage = -1.0f;
  Scan::process_pixel<W, H, false>(
      0, 0, p, sink, canvas, ring,
      [&](const Vector &, Fragment &f) {
        shader_coverage = f.v2;
        f.color = Color4(Pixel(60000, 60000, 60000), f.v2);
      },
      false, result, frag);

  HS_EXPECT_EQ(sink.count, 1);
  HS_EXPECT_GT(expected_coverage, 0.1f);
  HS_EXPECT_LT(expected_coverage, 0.9f);
  HS_EXPECT_NEAR(shader_coverage, expected_coverage, 1e-6f);
  HS_EXPECT_NEAR(sink.last_alpha, expected_coverage * expected_coverage, 1e-6f);
  HS_EXPECT_NEAR(sink.last_alpha, expected_coverage * legacy_factor, 2e-6f);

  SDF::PlanarPolygon solid(basis, 0.5f, 6, 0.0f);
  AlphaSink solid_sink;
  float solid_v2 = -1.0f;
  Scan::process_pixel<W, H, false>(
      0, 0, basis.v, solid_sink, canvas, solid,
      [&](const Vector &, Fragment &f) {
        solid_v2 = f.v2;
        f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
      },
      false, result, frag);
  HS_EXPECT_EQ(solid_sink.count, 1);
  HS_EXPECT_EQ(solid_v2, 0.0f);
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
  // The composite reproduces the winning child's own AA, and the thin/thick
  // falloffs are clearly separable so reproducing thin (not the wrapper-max
  // thick) is a meaningful distinction.
  HS_EXPECT_NEAR(a_union, a_thin, 1e-4f);
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

    HS_EXPECT_LE(a, prev + 1e-4f); // monotone non-increasing outward
    prev = a;

    if (a > 0.95f) saw_one = true;
    if (a < 0.05f) saw_zero = true;
    if (a > 0.1f && a < 0.9f) saw_mid = true;
  }

  // Centerline ~opaque, far edge ~transparent, with a genuine ramp between.
  HS_EXPECT_TRUE(saw_one);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_TRUE(saw_mid);
}

// ============================================================================
// Scan::Star / PlanarPolygon / Flower — filled-shape placement oracle
// ============================================================================

/**
 * @brief Reports whether any lit pixel lies in the given row.
 * @tparam W Canvas width in pixels.
 * @param fx Effect whose canvas is sampled.
 * @param y Row to scan.
 * @return True if a non-black pixel is found in row y.
 */
template <int W>
inline bool row_has_lit(const ScanFx &fx, int y) {
  for (int x = 0; x < W; ++x)
    if (!is_black(fx.get_pixel(x, y)))
      return true;
  return false;
}

/**
 * @brief Placement oracle for a filled cap: the pole the shape is centred on is
 *        lit and the opposite pole stays dark.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param fx Effect whose canvas was drawn into.
 * @param cap_north True if the shape caps the +Y pole (row 0); false if it caps
 *        the -Y pole (row H-1). Star/PlanarPolygon centre on basis.v (north);
 *        Flower centres on its antipode -basis.v (south).
 * @details Drawn to actual pixels — the placement, not just distance()/cull
 *          coverage. With an angular radius well under PI/2 the far pole sits
 *          outside the fill, so a projection sign flip (wrong hemisphere) or a
 *          collapsed fill is caught here where a bare draw-without-asserting is
 *          not.
 */
template <int W, int H>
inline void expect_filled_cap(const ScanFx &fx, bool cap_north) {
  const int near_row = cap_north ? 0 : H - 1;
  const int far_row = cap_north ? H - 1 : 0;
  HS_EXPECT_TRUE((row_has_lit<W>(fx, near_row)));
  HS_EXPECT_FALSE((row_has_lit<W>(fx, far_row)));
}

/** @brief Verifies a filled Star caps the basis.v (+Y) pole and not the other. */
inline void test_star_pixel_placement() {
  constexpr int W = 96, H = 64;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Star::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/5,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  expect_filled_cap<W, H>(fx, /*cap_north=*/true);
}

/** @brief Verifies a filled PlanarPolygon caps the basis.v (+Y) pole, not the other. */
inline void test_planar_polygon_pixel_placement() {
  constexpr int W = 96, H = 64;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  expect_filled_cap<W, H>(fx, /*cap_north=*/true);
}

/**
 * @brief Verifies a filled Flower caps the antipode (-Y) pole, not basis.v.
 * @details Flower's SDF scans from antipode = -basis.v, so the fill lands on the
 *          opposite pole from Star/PlanarPolygon — a placement detail only a
 *          pixel oracle pins.
 */
inline void test_flower_pixel_placement() {
  constexpr int W = 96, H = 64;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Flower::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  expect_filled_cap<W, H>(fx, /*cap_north=*/false);
}

/**
 * @brief Verifies overlapping strokes composite via the over operator at the
 *        shared pixel.
 * @details The 2D sink blends dst.lerp16(src, alpha) = dst*(1-a) + src*a. Two
 *          filled polygons both capping the +Y pole are drawn in order over a
 *          black frame, each at frag.alpha 0.5; their interiors are fully
 *          covered (AA alpha 1), so the pole pixel sees plot alpha exactly 0.5.
 *          First (red over black) yields red*0.5; second (green over that)
 *          yields red*0.25 + green*0.5. A draw-order or coverage bug (e.g.
 *          replacing instead of blending, or doubling coverage) moves the pole
 *          pixel off this composite.
 */
inline void test_overlapping_strokes_composite_blend() {
  constexpr int W = 96, H = 64;
  ScanFx fx(W, H);
  Pipeline<W, H> pipe;

  constexpr uint16_t red = 60000, green = 50000;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    // First fill: red at half alpha over black -> red*0.5.
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(red, 0, 0), 0.5f);
        });
    // Second fill: green at half alpha over the red -> red*0.25 + green*0.5.
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(0, green, 0), 0.5f);
        });
  }
  fx.advance_display();

  // Sample an interior column at the pole row, where both fills are fully
  // covered (AA alpha 1) so plot alpha is exactly frag.alpha 0.5.
  const Pixel &p = fx.get_pixel(W / 2, 0);
  // dst.lerp16 rounds to nearest; allow a few LSB of quantization slack.
  HS_EXPECT_NEAR((int)p.r, (int)(red * 0.25f), 64);
  HS_EXPECT_NEAR((int)p.g, (int)(green * 0.5f), 64);
  HS_EXPECT_EQ((int)p.b, 0);
}

// ============================================================================
// Scan::Volume / TransformedVolume — orthographic ray-march
// ============================================================================

/**
 * @brief Minimal analytic SDF: a sphere of `radius` centred at the local origin.
 * @details Supplies the distance() half of the Volume shape concept; wrapped in a
 * TransformedVolume to gain ray_to_local()/origin_to_local().
 */
struct SphereSDF {
  float radius; /**< Sphere radius in local units. */
  /**
   * @brief Signed distance to the sphere surface.
   * @param p Query point in local space.
   * @return |p| - radius (negative inside, zero on the surface).
   */
  float distance(const Vector &p) const { return p.length() - radius; }
};

/**
 * @brief Analytic SDF of two spheres unioned by min(): a small foreground sphere
 *        floating in front of a larger background sphere along the view axis.
 * @details The foreground silhouette has empty space immediately behind its edge
 * and the background surface a short march deeper, so a grazing ray at that edge
 * self-occludes the background — the case Volume::draw's occluder probe handles.
 */
struct TwoSphereSDF {
  Vector fg_center; /**< Foreground sphere centre (nearer the camera). */
  float fg_radius;  /**< Foreground sphere radius. */
  Vector bg_center; /**< Background sphere centre (deeper along the ray). */
  float bg_radius;  /**< Background sphere radius. */
  /**
   * @brief Signed distance to the union of the two spheres.
   * @param p Query point in local space.
   * @return The nearer of the two surface distances.
   */
  float distance(const Vector &p) const {
    return std::min((p - fg_center).length() - fg_radius,
                    (p - bg_center).length() - bg_radius);
  }
};

/**
 * @brief Capturing volume sink: records every plotted world position and its
 *        composited alpha.
 * @details Provides BOTH plot() overloads so it can back the type-erased
 * PipelineRef that Volume::draw takes; the 2D overload is never reached by the
 * volume path (it plots 3D world points) and is a no-op.
 */
struct VolumeSink {
  std::vector<Vector> plotted; /**< World positions handed to the 3D plot(). */
  std::vector<float> alpha;    /**< Composited alpha per plotted pixel. */
  void plot(Canvas &, const Vector &v, const Pixel &, float, float a) {
    plotted.push_back(v);
    alpha.push_back(a);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/**
 * @brief Verifies TransformedVolume's world<->local contract: the round trip is
 *        the identity, the mapped ray direction stays unit length, and distance()
 *        delegates to the wrapped SDF.
 * @details Volume::draw's per-step bounding-sphere cull is only correct when
 * ray_to_local is a rigid (length-preserving) map and bounds_center lands at the
 * local origin — both HS_CHECKed once per draw. This pins the transform math
 * those asserts rely on, independent of the rasterizer.
 */
inline void test_transformed_volume_world_local_roundtrip() {
  SphereSDF sphere{0.3f};
  const Vector center(0.2f, -0.5f, 0.8f);
  const Quaternion q =
      make_rotation(Vector(0.3f, 1.0f, -0.2f).normalized(), 0.7f);
  Scan::TransformedVolume vol(sphere, center, q);

  // bounds_center maps to the local origin (the cull precondition).
  Vector local_bc = vol.origin_to_local(center);
  HS_EXPECT_NEAR(local_bc.length(), 0.0f, 1e-5f);

  // A local point pushed out to world and back is recovered exactly.
  const Vector lp(0.1f, 0.2f, -0.25f);
  Vector world = center + rotate(lp, q);
  Vector back = vol.origin_to_local(world);
  HS_EXPECT_NEAR(back.x, lp.x, 1e-5f);
  HS_EXPECT_NEAR(back.y, lp.y, 1e-5f);
  HS_EXPECT_NEAR(back.z, lp.z, 1e-5f);

  // ray_to_local maps the origin like origin_to_local and keeps a unit direction
  // unit (rigid map, no scale) — the |local_vd| == 1 precondition.
  const Vector vd(0.0f, 0.0f, -1.0f);
  auto [lro, lvd] = vol.ray_to_local(world, vd);
  HS_EXPECT_NEAR(lro.x, lp.x, 1e-5f);
  HS_EXPECT_NEAR(lro.y, lp.y, 1e-5f);
  HS_EXPECT_NEAR(lro.z, lp.z, 1e-5f);
  HS_EXPECT_NEAR(lvd.length(), 1.0f, 1e-5f);

  // distance() forwards straight to the wrapped SDF.
  HS_EXPECT_NEAR(vol.distance(lp), sphere.distance(lp), 1e-6f);
}

/**
 * @brief Verifies Volume::draw ray-marches a sphere SDF into a bounded silhouette
 *        whose every shaded fragment's hit registers land on the surface, on the
 *        camera-facing cap.
 * @details Pins what the smoke loop never checks: (1) the rendered silhouette is
 * non-empty and strictly smaller than the canvas (a real hit set, not a full
 * clear or an empty frame), and the plotted set is a subset of the shaded set;
 * (2) each hit's frag.pos (closest_local) sits within the AA band of the sphere
 * surface (|pos| ≈ radius) and frag.size (closest_d) is inside the AA width — the
 * registers handed to the shader are genuine surface hits; (3) the hit centroid
 * lies on the +Z hemisphere, i.e. the visible cap faces the camera (rays travel
 * along -Z), so a projection/back-face sign flip would be caught.
 */
inline void test_volume_raymarch_silhouette_and_registers() {
  constexpr int W = 96, H = 64;
  const Vector center(0.0f, 0.0f, 1.0f);    // bounds centre in LED space
  const Vector view_dir(0.0f, 0.0f, -1.0f); // camera -> scene
  const float bounds_radius = 0.35f;
  const float sphere_r = 0.28f; // < bounds so the SDF fits the cull sphere
  const float aa_width = 0.01f;

  SphereSDF sphere{sphere_r};
  Scan::TransformedVolume vol(sphere, center, Quaternion());

  ScanFx fx(W, H);
  VolumeSink sink;

  int hits = 0;
  float max_surf_err = 0.0f; // worst |‖pos‖ - radius| over all hits
  float max_reg_d = 0.0f;    // worst |frag.size| (closest_d) over all hits
  Vector centroid_sum(0.0f, 0.0f, 0.0f);
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, center, bounds_radius, view_dir, vol,
        [&](const Vector &loc, Fragment &frag) {
          ++hits;
          max_surf_err =
              std::max(max_surf_err, std::fabs(loc.length() - sphere_r));
          max_reg_d = std::max(max_reg_d, std::fabs(frag.size));
          centroid_sum = centroid_sum + loc;
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/24, aa_width);
  }

  // (1) Real silhouette: some fragments, never the whole canvas; plotted ⊆ shaded.
  HS_EXPECT_GT(hits, 0);
  HS_EXPECT_LT((size_t)hits, (size_t)(W * H));
  HS_EXPECT_GT(sink.plotted.size(), (size_t)0);
  HS_EXPECT_LE(sink.plotted.size(), (size_t)hits);

  // (2) Every shaded fragment is a genuine surface hit inside the AA band.
  HS_EXPECT_LE(max_surf_err, aa_width + 1e-3f);
  HS_EXPECT_LE(max_reg_d, aa_width);

  // (3) The hit centroid is on the camera-facing (+Z) cap.
  Vector centroid = centroid_sum * (1.0f / static_cast<float>(hits));
  HS_EXPECT_GT(centroid.z, 0.1f);
}

/**
 * @brief Verifies Volume::draw antialiases a self-occlusion edge over the surface
 *        behind it rather than fading the edge to black.
 * @details A small foreground sphere floats in front of a larger background sphere.
 * Around the foreground silhouette the grazing ray lands in the AA band while the
 * background sits a short march deeper, so probe_occluder reports a solid surface
 * and Volume::draw takes the occluded-edge branch: it lays the shaded background
 * down, then blends the foreground over it by the edge coverage — emitting TWO
 * plots at the same world position (background first, foreground second). This
 * pins that branch (the most intricate, previously untested logic in scan.h): the
 * duplicate-position signature must appear, the background must be laid down
 * opaque, and the foreground must be a partial (0<α<1) blend over it, not a fade
 * to black.
 */
inline void test_volume_draw_occluded_edge_blends_over_background() {
  constexpr int W = 96, H = 64;
  const Vector center(0.0f, 0.0f, 1.0f);
  const Vector view_dir(0.0f, 0.0f, -1.0f); // rays travel along -Z
  const float bounds_radius = 0.50f;
  const float aa_width = 0.01f;

  // Small foreground sphere nearer the camera (+Z local); larger background
  // sphere deeper and wider so it sits behind the whole foreground silhouette.
  // The modest step budget lets the grazing ray stall on the foreground edge
  // (landing in the AA band) rather than reaching the background solidly, so the
  // occluder probe — not the main trace — is what discovers the surface behind.
  TwoSphereSDF shape{Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     Vector(0.0f, 0.0f, -0.20f), 0.30f};
  Scan::TransformedVolume vol(shape, center, Quaternion());

  ScanFx fx(W, H);
  VolumeSink sink;
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, center, bounds_radius, view_dir, vol,
        [&](const Vector &, Fragment &frag) {
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/12, aa_width);
  }

  // The occluder branch is the only path that plots the same world position
  // twice back-to-back (background then foreground); every other path plots each
  // pixel once, and distinct pixels never share a world position.
  int occ_pairs = 0;
  float bg_alpha = 0.0f, fg_alpha = 1.0f;
  for (size_t i = 1; i < sink.plotted.size(); ++i) {
    if ((sink.plotted[i] - sink.plotted[i - 1]).length() < 1e-5f) {
      ++occ_pairs;
      if (occ_pairs == 1) {
        bg_alpha = sink.alpha[i - 1];
        fg_alpha = sink.alpha[i];
      }
    }
  }

  HS_EXPECT_GT(occ_pairs, 0);
  // Background laid down opaque; foreground blended over it as a partial edge, so
  // the edge reads over the surface instead of fading to black.
  HS_EXPECT_GT(bg_alpha, 0.5f);
  HS_EXPECT_GT(fg_alpha, 0.0f);
  HS_EXPECT_LT(fg_alpha, bg_alpha);
}

/**
 * @brief Verifies trace_closest stops at the first silhouette graze instead of
 *        letting an occluded surface behind it steal the closest approach.
 * @details A ray grazes the foreground sphere's edge inside the AA band and then
 * passes solidly through the background sphere. The returned distance must be
 * the foreground graze and the returned point must sit on the foreground
 * silhouette, independent of the step budget — with a generous budget an
 * unguarded march reaches the background, collapsing the edge alpha to opaque
 * and moving the shading point to the wrong surface.
 */
inline void test_volume_trace_closest_stops_at_first_graze() {
  const float aa_width = 0.01f;
  TwoSphereSDF shape{Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     Vector(0.0f, 0.0f, -0.20f), 0.30f};

  // Grazing ray: passes 0.005 outside the foreground silhouette, then through
  // the background sphere's interior.
  const Vector ro(0.185f, 0.0f, 1.0f);
  const Vector vd(0.0f, 0.0f, -1.0f);

  for (int max_steps : {14, 40}) {
    Vector closest_local;
    float closest_d = Scan::Volume::trace_closest(shape, ro, vd, 0.5f,
                                                  max_steps, aa_width,
                                                  closest_local);
    HS_EXPECT_GT(closest_d, 0.003f);
    HS_EXPECT_LT(closest_d, 0.0075f);
    HS_EXPECT_GT(closest_local.z, 0.1f);
  }
}

/**
 * @brief Verifies probe_occluder reports the grazed background edge's own
 *        closest-approach point, so the corner fill is shaded on the background
 *        surface rather than reusing the foreground fragment.
 * @details The ray passes through the corner where the two spheres' silhouettes
 * cross, inside both AA bands but hitting neither: the trace grazes the
 * foreground sphere, and the probe must come back non-solid with the graze
 * coverage and a `behind` point at the background sphere's edge (negative z),
 * not the foreground graze it was seeded with.
 */
inline void test_volume_probe_occluder_reports_background_graze_point() {
  const float aa_width = 0.01f;
  const float hit_threshold = aa_width * 0.1f;
  TwoSphereSDF shape{Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     Vector(0.28f, 0.0f, -0.20f), 0.30f};

  // Silhouette-crossing corner, offset outward of both circles by ~0.004.
  const Vector ro(0.0357f, 0.1797f, 1.0f);
  const Vector vd(0.0f, 0.0f, -1.0f);

  Vector closest_local;
  float closest_d = Scan::Volume::trace_closest(shape, ro, vd, 0.6f, 40,
                                                aa_width, closest_local);
  HS_EXPECT_GT(closest_d, hit_threshold);
  HS_EXPECT_LT(closest_d, aa_width);
  HS_EXPECT_GT(closest_local.z, 0.1f);

  auto occ = Scan::Volume::probe_occluder(shape, closest_local, vd, 0.6f,
                                          hit_threshold, aa_width);
  HS_EXPECT_FALSE(occ.solid);
  // The ray passes 0.0033 outside the background sphere, so the analytic
  // coverage is quintic(1 - (0.0033 - 0.001)/0.009) ~= 0.89; the parabolic
  // refinement must land near it (the coarse stride alone reads ~0.65).
  HS_EXPECT_GT(occ.soft, 0.8f);
  HS_EXPECT_LT(occ.soft, 0.95f);
  HS_EXPECT_LT(occ.behind.z, -0.05f);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every scan test under the "scan" module scope.
 * @return The number of test failures recorded by the module.
 */
inline int run_scan_tests() {
  hs_test::ModuleFixture fixture("scan");

  test_shader_constant_fills_canvas();
  test_shader_ssaa_premultiplies_partial_coverage();
  test_shader_positional_maps_latitude();
  test_shader_respects_clip_band();
  test_ring_rasterize_produces_bounded_output();
  test_ring_rasterize_lit_pixels_on_band();
  test_ring_rasterize_lights_expected_row();
  test_stroke_aa_is_monotone_ramp();
  test_ring_rasterize_empty_clip_draws_nothing();
  test_distorted_ring_flat_matches_zero_knot_raster();
  test_scan_shader_v2_contract();
  test_scan_region_seam_no_double_plot();
  test_scan_region_fractional_boundary_no_double_plot();
  test_scan_region_clip_arc_matches_predicate();
  test_plot_line_over_pole_reaches_row0();
  test_csg_stroke_aa_uses_winning_child_thickness();

  test_star_pixel_placement();
  test_planar_polygon_pixel_placement();
  test_flower_pixel_placement();
  test_overlapping_strokes_composite_blend();

  test_transformed_volume_world_local_roundtrip();
  test_volume_raymarch_silhouette_and_registers();
  test_volume_draw_occluded_edge_blends_over_background();
  test_volume_trace_closest_stops_at_first_graze();
  test_volume_probe_occluder_reports_background_graze_point();

  return fixture.result();
}

} // namespace scan_tests
} // namespace hs_test

