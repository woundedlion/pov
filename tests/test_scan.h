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
#include "core/filter.h"
#include "core/canvas.h"
#include "core/geometry.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace scan_tests {

struct ScanFx : public Effect {
  ScanFx(int W, int H) : Effect(W, H) {}
  void draw_frame() override {}
  bool show_bg() const override { return false; }
};

inline bool is_black(const Pixel &p) { return p.r == 0 && p.g == 0 && p.b == 0; }

// ============================================================================
// Scan::Shader::draw — full-sphere per-pixel shader
// ============================================================================

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

// scan_region seam coalescer: a span crossing x=0 must not double-plot the
// wrapped overlap with another span. Drives scan_region directly with a sorted
// two-span row (a low span + a seam-crosser, as a shape's sorted output emits).
// The old coalescer worked in unwrapped space and double-plotted the overlap.
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

  // No pixel plotted more than once (the wrapped overlap at x=0,1 is not
  // doubled — without the fix counts[1] == 2).
  for (int x = 0; x < W; ++x)
    HS_EXPECT_LE(counts[x], 1);

  // Coverage is exactly {0,1,2,3,4, W-2, W-1}.
  const int covered[] = {0, 1, 2, 3, 4, W - 2, W - 1};
  for (int x : covered)
    HS_EXPECT_EQ(counts[x], 1);
  HS_EXPECT_EQ(counts[5], 0);
  HS_EXPECT_EQ(counts[W - 3], 0);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_scan_tests() {
  auto scope = hs_test::begin_module("scan");

  test_shader_constant_fills_canvas();
  test_shader_positional_maps_latitude();
  test_shader_respects_clip_band();
  test_ring_rasterize_produces_bounded_output();
  test_ring_rasterize_empty_clip_draws_nothing();
  test_scan_region_seam_no_double_plot();

  return hs_test::end_module(scope);
}

} // namespace scan_tests
} // namespace hs_test

