/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/filter.h.
 *
 * Focus: compile-time trait machinery + pure per-call helper logic that runs
 * without a live Canvas/Effect render.
 *
 * Coverage:
 *   - Filter traits Is2D / Is3D / Is2DWithHistory / Is3DWithHistory
 *     (is_2d / has_history members)
 *   - Trait inheritance on representative filters (AntiAlias, Blur,
 *     ChromaticShift, World::Replicate, World::Trails, Screen::Trails,
 *     Pixel::Feedback)
 *   - Pipeline<W,H>::is_2d sink flag and Pipeline::get<T>() type-correct lookup
 *   - Screen::AntiAlias::plot — bilinear weight partition (sums to alpha),
 *     pole snap behaviour
 *   - Screen::Blur::plot — kernel passthrough (factor=0 → identity center,
 *     factor=1 → 3x3 weights sum to alpha)
 *   - Pixel::ChromaticShift::plot — channel-split fan-out
 *   - Pixel::Feedback — Style binding accessor, set_enabled
 *
 * End-to-end (live Canvas via the test_canvas/test_scan advance_display pattern,
 * which dissolves the buffer_free() ctor spin):
 *   - Pipeline 2D sink plot — int + float overloads, alpha blend, x-wrap, clip
 *   - Pipeline 3D sink plot — vector_to_pixel routing to the Canvas
 *   - World filter routing — World::Replicate fans out through the 3D->2D sink
 *   - 2D->3D mismatch — a 2D coord into a 3D-headed pipeline round-trips
 *     pixel_to_vector -> vector_to_pixel
 *   - Screen filter routing — Screen::AntiAlias forwards to the sink
 *   - Pixel::Feedback::flush — warp-field flush blends the (faded) prev frame
 *
 * STILL SKIPPED (covered elsewhere or need separate scaffolding):
 *   - World filters that call tween(Orientation, ...) — need a populated
 *     Orientation and exercise geometry, covered by test_geometry.
 *   - World::Trails encode/decode + ring buffer (private; needs init_storage
 *     with an Arena and plot/flush driving).
 */
#pragma once

#include "core/filter.h"
#include "core/canvas.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace filter_tests {

// ============================================================================
// Trait structs — direct member values
// ============================================================================

inline void test_trait_member_values() {
  HS_EXPECT_TRUE(Is2D::is_2d);
  HS_EXPECT_FALSE(Is2D::has_history);

  HS_EXPECT_FALSE(Is3D::is_2d);
  HS_EXPECT_FALSE(Is3D::has_history);

  HS_EXPECT_TRUE(Is2DWithHistory::is_2d);
  HS_EXPECT_TRUE(Is2DWithHistory::has_history);

  HS_EXPECT_FALSE(Is3DWithHistory::is_2d);
  HS_EXPECT_TRUE(Is3DWithHistory::has_history);
}

// ============================================================================
// Trait inheritance on representative filters
// ============================================================================

inline void test_filter_trait_inheritance() {
  constexpr int W = 32, H = 32;

  // 2D screen-space, stateless.
  HS_EXPECT_TRUE((Filter::Screen::AntiAlias<W, H>::is_2d));
  HS_EXPECT_FALSE((Filter::Screen::AntiAlias<W, H>::has_history));
  HS_EXPECT_TRUE((Filter::Screen::Blur<W, H>::is_2d));
  HS_EXPECT_FALSE((Filter::Screen::Blur<W, H>::has_history));

  // Pixel-space, stateless (tagged Is2D).
  HS_EXPECT_TRUE((Filter::Pixel::ChromaticShift<W>::is_2d));
  HS_EXPECT_FALSE((Filter::Pixel::ChromaticShift<W>::has_history));

  // Pixel-space feedback is 2D *with* history.
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::is_2d));
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::has_history));

  // 3D world-space, stateless.
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::is_2d));
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::has_history));
  HS_EXPECT_FALSE((Filter::World::Hole<W>::is_2d));
  HS_EXPECT_FALSE((Filter::World::Hole<W>::has_history));

  // History-bearing trail filters.
  HS_EXPECT_FALSE((Filter::World::Trails<W, 16>::is_2d));
  HS_EXPECT_TRUE((Filter::World::Trails<W, 16>::has_history));
  HS_EXPECT_TRUE((Filter::Screen::Trails<W>::is_2d));
  HS_EXPECT_TRUE((Filter::Screen::Trails<W>::has_history));

  // Static-assert form (compile-time): mirrors the runtime checks above.
  static_assert(Filter::Screen::AntiAlias<W, H>::is_2d, "AntiAlias is 2D");
  static_assert(!Filter::World::Replicate<W>::is_2d, "Replicate is 3D");
  static_assert(Filter::Pixel::Feedback<W, H>::has_history,
                "Feedback keeps history");
}

// ============================================================================
// Pipeline sink + get<T>()
// ============================================================================

inline void test_pipeline_sink_is_2d() {
  HS_EXPECT_TRUE((Pipeline<32, 32>::is_2d));
}

inline void test_pipeline_get_returns_correct_filter() {
  constexpr int W = 32, H = 32;
  using AA = Filter::Screen::AntiAlias<W, H>;
  using CS = Filter::Pixel::ChromaticShift<W>;
  using Blur = Filter::Screen::Blur<W, H>;

  // Compose a 3-filter pipeline (all constructible without a Canvas).
  Pipeline<W, H, AA, Blur, CS> pipe(AA{}, Blur{1.0f}, CS{});

  // get<T>() must return a reference whose static type is T and which is the
  // base subobject of the pipeline node.
  AA &aa = pipe.get<AA>();
  Blur &bl = pipe.get<Blur>();
  CS &cs = pipe.get<CS>();

  static_assert(std::is_same_v<decltype(pipe.get<AA>()), AA &>,
                "get<AA>() returns AA&");
  static_assert(std::is_same_v<decltype(pipe.get<Blur>()), Blur &>,
                "get<Blur>() returns Blur&");
  static_assert(std::is_same_v<decltype(pipe.get<CS>()), CS &>,
                "get<CS>() returns CS&");

  // The head filter is also the pipeline itself upcast to AA.
  HS_EXPECT_TRUE(static_cast<AA *>(&pipe) == &aa);
  // Tail nodes resolve to distinct subobjects.
  HS_EXPECT_TRUE(static_cast<const void *>(&bl) !=
                 static_cast<const void *>(&cs));

  // const overload of get<T>() is usable.
  const Pipeline<W, H, AA, Blur, CS> &cpipe = pipe;
  const Blur &cbl = cpipe.get<Blur>();
  HS_EXPECT_TRUE(&cbl == &bl);
}

// ============================================================================
// Screen::AntiAlias::plot — pure bilinear weight partition
// ============================================================================

inline void test_antialias_weights_partition() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // Interior, non-pole sample with a clear fractional offset.
  // The 4 distributed alphas must sum back to the input alpha (partition of
  // unity), because v00+v10+v01+v11 == 1 by construction.
  float sum = 0.0f;
  int count = 0;
  const float in_alpha = 0.8f;
  aa.plot(10.3f, 20.6f, Pixel(1, 2, 3), 0.0f, in_alpha,
          [&](float, float, const Pixel &, float, float a) {
            sum += a;
            ++count;
          });
  HS_EXPECT_GE(count, 1);
  HS_EXPECT_LE(count, 4);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

inline void test_antialias_integer_coord_single_tap() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // Exact integer coordinates → fractional parts are 0 → only the (x0,y0) tap
  // carries weight (v00 == 1, the rest are 0 and fall below 1e-8 threshold).
  int count = 0;
  float kept_alpha = 0.0f;
  float kept_x = -1.0f, kept_y = -1.0f;
  aa.plot(12.0f, 24.0f, Pixel(0, 0, 0), 0.0f, 0.5f,
          [&](float x, float y, const Pixel &, float, float a) {
            ++count;
            kept_alpha = a;
            kept_x = x;
            kept_y = y;
          });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kept_alpha, 0.5f, 1e-5f);
  HS_EXPECT_NEAR(kept_x, 12.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_y, 24.0f, 1e-5f);
}

inline void test_antialias_seam_wraps_left_column() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // A sub-pixel coordinate just left of the theta=0 seam (x = -0.3) must
  // distribute across the wrapped columns W-1 and 0, not collapse onto a single
  // unwrapped column. y at the equator (y_m = 0) isolates the two X taps.
  // Regression: the old truncating modf put a negative fraction onto column 0
  // only and left the left neighbor unwrapped.
  const float in_alpha = 1.0f;
  float sum = 0.0f;
  bool all_in_range = true, saw_w_minus_1 = false, saw_zero = false;
  aa.plot(-0.3f, 32.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
          [&](float x, float y, const Pixel &, float, float a) {
            (void)y;
            sum += a;
            if (x < 0.0f || x >= static_cast<float>(W)) all_in_range = false;
            int xi = static_cast<int>(x);
            if (xi == W - 1) saw_w_minus_1 = true;
            if (xi == 0) saw_zero = true;
          });
  HS_EXPECT_TRUE(all_in_range);
  HS_EXPECT_TRUE(saw_w_minus_1);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

// ============================================================================
// Screen::Blur::plot — kernel passthrough
// ============================================================================

inline void test_blur_factor_zero_is_identity() {
  constexpr int W = 32, H = 32;
  // factor=0 → center weight c = 1.0, edge/corner weights 0 → single tap.
  Filter::Screen::Blur<W, H> blur(0.0f);

  int count = 0;
  float kept_alpha = 0.0f;
  float kept_x = -1.0f, kept_y = -1.0f;
  blur.plot(8.0f, 16.0f, Pixel(4, 5, 6), 2.0f, 1.0f,
            [&](float x, float y, const Pixel &, float, float a) {
              ++count;
              kept_alpha = a;
              kept_x = x;
              kept_y = y;
            });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kept_alpha, 1.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_x, 8.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_y, 16.0f, 1e-5f);
}

inline void test_blur_full_kernel_sums_to_alpha() {
  constexpr int W = 32, H = 32;
  // factor=1 → full Gaussian (center 0.25, edge 0.125, corner 0.0625).
  // For an interior pixel all 9 taps are in-bounds; weights sum to 1, so the
  // distributed alpha sums to the input alpha.
  Filter::Screen::Blur<W, H> blur(1.0f);

  float sum = 0.0f;
  int count = 0;
  const float in_alpha = 0.6f;
  blur.plot(15.0f, 16.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
            [&](float, float, const Pixel &, float, float a) {
              sum += a;
              ++count;
            });
  HS_EXPECT_EQ(count, 9);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

inline void test_blur_update_changes_kernel() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(1.0f);

  // After update(0) the kernel collapses to identity (single center tap).
  blur.update(0.0f);
  int count = 0;
  blur.plot(15.0f, 16.0f, Pixel(1, 1, 1), 0.0f, 1.0f,
            [&](float, float, const Pixel &, float, float) { ++count; });
  HS_EXPECT_EQ(count, 1);
}

// ============================================================================
// Pixel::ChromaticShift::plot — channel-split fan-out
// ============================================================================

inline void test_chromatic_shift_fanout() {
  constexpr int W = 64;
  Filter::Pixel::ChromaticShift<W> cs;

  struct Tap {
    float x, y;
    Pixel c;
    float alpha;
  };
  Tap taps[8];
  int count = 0;
  Pixel src(100, 150, 200);
  cs.plot(10.0f, 5.0f, src, 0.0f, 1.0f,
          [&](float x, float y, const Pixel &c, float, float a) {
            if (count < 8) taps[count] = {x, y, c, a};
            ++count;
          });

  // Original + 3 channel-shifted copies.
  HS_EXPECT_EQ(count, 4);

  // First tap is the unmodified colour at the original coordinate.
  HS_EXPECT_NEAR(taps[0].x, 10.0f, 1e-5f);
  HS_EXPECT_NEAR(taps[0].y, 5.0f, 1e-5f);
  HS_EXPECT_EQ(taps[0].c.r, src.r);
  HS_EXPECT_EQ(taps[0].c.g, src.g);
  HS_EXPECT_EQ(taps[0].c.b, src.b);

  // Red-only copy at x+1.
  HS_EXPECT_NEAR(taps[1].x, 11.0f, 1e-5f);
  HS_EXPECT_EQ(taps[1].c.r, src.r);
  HS_EXPECT_EQ(taps[1].c.g, 0);
  HS_EXPECT_EQ(taps[1].c.b, 0);

  // Green-only copy at x+2.
  HS_EXPECT_NEAR(taps[2].x, 12.0f, 1e-5f);
  HS_EXPECT_EQ(taps[2].c.r, 0);
  HS_EXPECT_EQ(taps[2].c.g, src.g);
  HS_EXPECT_EQ(taps[2].c.b, 0);

  // Blue-only copy at x+3.
  HS_EXPECT_NEAR(taps[3].x, 13.0f, 1e-5f);
  HS_EXPECT_EQ(taps[3].c.r, 0);
  HS_EXPECT_EQ(taps[3].c.g, 0);
  HS_EXPECT_EQ(taps[3].c.b, src.b);
}

// ============================================================================
// Pixel::Feedback — Style binding + enable flag (no flush / no Canvas)
// ============================================================================

inline void test_feedback_style_binding() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  // The accessor returns the same Style object that was bound.
  HS_EXPECT_TRUE(&fb.style() == &style);

  const Filter::Pixel::Feedback<W, H> &cfb = fb;
  HS_EXPECT_TRUE(&cfb.style() == &style);

  // Smoke preset values are reflected through the binding.
  HS_EXPECT_NEAR(fb.style().fade, 0.9f, 1e-6f);
  HS_EXPECT_EQ(fb.style().downsample, 4);

  // Mutating through the accessor is visible on the original Style.
  fb.style().fade = 0.5f;
  HS_EXPECT_NEAR(style.fade, 0.5f, 1e-6f);

  // set_enabled compiles and is callable.
  fb.set_enabled(false);
  fb.set_enabled(true);
}

inline void test_feedback_plot_is_passthrough() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  // plot() is a pure pass-through (no Canvas touched): one tap, unchanged.
  int count = 0;
  float kx = -1, ky = -1, ka = -1;
  Pixel src(7, 8, 9);
  fb.plot(3.0f, 4.0f, src, 1.5f, 0.75f,
          [&](float x, float y, const Pixel &, float, float a) {
            ++count;
            kx = x;
            ky = y;
            ka = a;
          });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kx, 3.0f, 1e-6f);
  HS_EXPECT_NEAR(ky, 4.0f, 1e-6f);
  HS_EXPECT_NEAR(ka, 0.75f, 1e-6f);
}

// ============================================================================
// End-to-end pipeline routing through a live Canvas
//
// The test_canvas/test_scan pattern: a Canvas spins in its ctor while
// !buffer_free(), so every drawn frame MUST advance_display() before the next
// Canvas is constructed. With that, the full Pipeline::plot/flush routing
// (which the per-call helper tests above cannot reach) is exercisable.
// ============================================================================

// Minimal concrete Effect for binding a live Canvas.
struct PipeFx : public Effect {
  PipeFx(int W, int H) : Effect(W, H) {}
  void draw_frame() override {}
  bool show_bg() const override { return false; }
};

inline bool px_black(const Pixel &p) {
  return p.r == 0 && p.g == 0 && p.b == 0;
}

inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

inline size_t count_lit(const PipeFx &fx) {
  size_t n = 0;
  for (int y = 0; y < fx.height(); ++y)
    for (int x = 0; x < fx.width(); ++x)
      if (!px_black(fx.get_pixel(x, y))) ++n;
  return n;
}

// Bare 2D sink: int + float overloads, exact (alpha=1) write, x-wrap, clip.
inline void test_pipeline_sink_2d_plot_blends_wraps_clips() {
  constexpr int W = 16, H = 8;
  PipeFx fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    // Integer overload, full alpha over a black buffer -> exact source colour.
    pipe.plot(c, 4, 3, Pixel(100, 200, 300), 0.0f, 1.0f);
    // x = W+2 is in the sink's [-W, 2W) contract window and wraps to column 2.
    pipe.plot(c, W + 2, 5, Pixel(10, 20, 30), 0.0f, 1.0f);
    // Float overload rounds to the nearest pixel.
    pipe.plot(c, 7.4f, 1.6f, Pixel(1, 2, 3), 0.0f, 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(4, 3), 100, 200, 300));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 5), 10, 20, 30)); // wrapped
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(7, 2), 1, 2, 3));    // rounded

  // A second frame: a plot outside the clip band is dropped.
  PipeFx fx2(W, H);
  fx2.set_clip(2, 5, 0, W); // rows [2,5)
  fx2.clip.margin = 0;
  {
    Canvas c(fx2);
    pipe.plot(c, 8, 6, Pixel(500, 0, 0), 0.0f, 1.0f); // row 6 outside band
    pipe.plot(c, 8, 3, Pixel(0, 500, 0), 0.0f, 1.0f); // row 3 inside band
  }
  fx2.advance_display();
  HS_EXPECT_TRUE(px_black(fx2.get_pixel(8, 6)));
  HS_EXPECT_TRUE(pix_eq(fx2.get_pixel(8, 3), 0, 500, 0));
}

// 3D sink overload: a unit vector is routed via vector_to_pixel and written.
inline void test_pipeline_sink_3d_plot_routes_to_canvas() {
  constexpr int W = 32, H = 16;
  PipeFx fx(W, H);
  Pipeline<W, H> pipe;
  Vector v = Vector(0.6f, 0.4f, 0.69f).normalized();
  {
    Canvas c(fx);
    pipe.plot(c, v, Pixel(40000, 20000, 10000), 0.0f, 1.0f);
  }
  fx.advance_display();

  // Exactly one pixel lit, carrying the source colour, at the coordinate
  // vector_to_pixel predicts.
  HS_EXPECT_EQ(count_lit(fx), (size_t)1);
  PixelCoords pc = vector_to_pixel<W, H>(v);
  int ex = fast_wrap(static_cast<int>(std::round(pc.x)), W);
  int ey = static_cast<int>(std::round(pc.y));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(ex, ey), 40000, 20000, 10000));
}

// World filter (3D) head fans out through the 3D->2D sink conversion.
inline void test_pipeline_world_replicate_fans_out() {
  constexpr int W = 32, H = 16;
  PipeFx fx(W, H);
  // Replicate(2): original + one copy rotated 180 deg about Y (same latitude,
  // longitude + W/2) -> two distinct columns -> two lit pixels.
  Pipeline<W, H, Filter::World::Replicate<W>> pipe(Filter::World::Replicate<W>(2));
  Vector v = Vector(0.6f, 0.4f, 0.69f).normalized();
  {
    Canvas c(fx);
    pipe.plot(c, v, Pixel(60000, 60000, 60000), 0.0f, 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_EQ(count_lit(fx), (size_t)2);
}

// 2D coordinate into a 3D-headed pipeline: the float-plot mismatch branch maps
// pixel_to_vector, the identity filter passes it, and the sink maps it back.
inline void test_pipeline_2d_into_3d_head_roundtrips() {
  constexpr int W = 32, H = 16;
  PipeFx fx(W, H);
  // Replicate(1) clamps to a single emission -> identity pass-through.
  Pipeline<W, H, Filter::World::Replicate<W>> pipe(Filter::World::Replicate<W>(1));
  const int px = 16, py = 8; // interior, away from poles/seam
  {
    Canvas c(fx);
    pipe.plot(c, static_cast<float>(px), static_cast<float>(py),
              Pixel(0, 60000, 0), 0.0f, 1.0f);
  }
  fx.advance_display();

  // One pixel lit, within a pixel of the input after the trig round-trip.
  HS_EXPECT_EQ(count_lit(fx), (size_t)1);
  bool near_input = false;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!px_black(fx.get_pixel(x, y)))
        near_input = (std::abs(x - px) <= 1 && std::abs(y - py) <= 1);
  HS_EXPECT_TRUE(near_input);
}

// Screen filter (2D) head forwards to the sink. An integer coordinate collapses
// AntiAlias to a single tap, so the routed result is one exact pixel.
inline void test_pipeline_screen_antialias_routes_to_sink() {
  constexpr int W = 32, H = 16;
  PipeFx fx(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipe;
  {
    Canvas c(fx);
    pipe.plot(c, 10.0f, 5.0f, Pixel(50000, 0, 25000), 0.0f, 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_EQ(count_lit(fx), (size_t)1);
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(10, 5), 50000, 0, 25000));
}

// Pixel::Feedback::flush warp-field path. A Smoke style with no bound
// NoiseParams makes noise_warp an identity map, so the flush blends the previous
// frame back at the same location, faded by style.fade — deterministic.
inline void test_feedback_flush_blends_prev_frame() {
  constexpr int W = 32, H = 16; // both divisible by Smoke's downsample (4)
  PipeFx fx(W, H);
  ::Feedback::Style style = ::Feedback::Style::Smoke(); // noise stays nullptr
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  auto trail = [](float, float, float) { return Color4(Pixel(0, 0, 0), 0.0f); };

  // Frame 1: a bright horizontal band becomes the "previous" frame.
  {
    Canvas c(fx);
    for (int y = 6; y <= 9; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(40000, 40000, 40000);
  }
  fx.advance_display();

  // Frame 2: empty buffer; flush warps+blends the prev band into it.
  {
    Canvas c(fx);
    pipe.flush(c, ScreenTrailFn(trail), 1.0f);
  }
  fx.advance_display();

  // Interior band row carried over and was faded (fade=0.9 -> ~36000), not the
  // original brightness. A row far from the band stays black (identity warp,
  // no spill that far).
  const Pixel &band = fx.get_pixel(W / 2, 8);
  HS_EXPECT_FALSE(px_black(band));
  HS_EXPECT_GT((int)band.r, 30000);
  HS_EXPECT_LT((int)band.r, 40000);
  HS_EXPECT_TRUE(px_black(fx.get_pixel(W / 2, 0)));

  // Disabled feedback short-circuits flush: a fresh frame stays black.
  pipe.get<Filter::Pixel::Feedback<W, H>>().set_enabled(false);
  {
    Canvas c(fx);
    pipe.flush(c, ScreenTrailFn(trail), 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_TRUE(px_black(fx.get_pixel(W / 2, 8)));
}

// ============================================================================
// Runner
// ============================================================================

inline int run_filter_tests() {
  auto scope = hs_test::begin_module("filter");

  test_trait_member_values();
  test_filter_trait_inheritance();

  test_pipeline_sink_is_2d();
  test_pipeline_get_returns_correct_filter();

  test_antialias_weights_partition();
  test_antialias_integer_coord_single_tap();
  test_antialias_seam_wraps_left_column();

  test_blur_factor_zero_is_identity();
  test_blur_full_kernel_sums_to_alpha();
  test_blur_update_changes_kernel();

  test_chromatic_shift_fanout();

  test_feedback_style_binding();
  test_feedback_plot_is_passthrough();

  test_pipeline_sink_2d_plot_blends_wraps_clips();
  test_pipeline_sink_3d_plot_routes_to_canvas();
  test_pipeline_world_replicate_fans_out();
  test_pipeline_2d_into_3d_head_roundtrips();
  test_pipeline_screen_antialias_routes_to_sink();
  test_feedback_flush_blends_prev_frame();

  return hs_test::end_module(scope);
}

} // namespace filter_tests
} // namespace hs_test

