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
 *     factor=1 → 3x3 weights sum to alpha), pole-row clip renormalization
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
 *   - World::Orient tween over a populated Orientation history
 *     (test_world_orient_motion_blur_sweep_ages)
 *   - World::Trails int16 encode/decode + ring buffer / ttl lifecycle
 *     (test_world_trails_* — init_storage(arena) + plot/flush driving)
 */
#pragma once

#include <array>
#include <span>

#include "core/filter.h"
#include "core/canvas.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace filter_tests {

// ============================================================================
// Trait structs — direct member values
// ============================================================================

/**
 * @brief Verifies each trait tag exposes the expected is_2d / has_history member
 *        constants.
 */
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

/**
 * @brief Verifies representative filters inherit the right is_2d / has_history
 *        traits, at both runtime and compile time.
 */
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

/**
 * @brief Verifies the bare (filter-free) pipeline sink is 2D.
 */
inline void test_pipeline_sink_is_2d() {
  HS_EXPECT_TRUE((Pipeline<32, 32>::is_2d));
}

/**
 * @brief Verifies get<T>() resolves each composed filter to the correctly-typed
 *        base subobject, for both the head node and tail nodes, in const and
 *        non-const pipelines.
 */
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

  // Absent-filter case: get<T>() on a pipeline lacking T must fail to compile
  // cleanly (the dependent-false static_assert "Filter type T not found"). It
  // cannot be exercised at runtime — instantiating the base get<T>() is a hard
  // error, not a SFINAE-detectable one — so it lives behind a default-off macro.
  // Verify by hand: build with -DHS_TEST_PIPELINE_GET_ABSENT and confirm the
  // sole diagnostic is the intended message.
#ifdef HS_TEST_PIPELINE_GET_ABSENT
  using Absent = Filter::Pixel::Feedback<W, H>; // not in the pipeline above
  (void)pipe.get<Absent>();
#endif
}

// ============================================================================
// Screen::AntiAlias::plot — pure bilinear weight partition
// ============================================================================

/**
 * @brief Verifies a fractional interior sample splits into 1..4 taps whose
 *        alphas sum back to the input alpha.
 * @details Partition of unity: the four bilinear weights sum to 1, so the per-tap
 *          alphas sum to the input alpha.
 */
inline void test_antialias_weights_partition() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // Interior, non-pole sample with a clear fractional offset; the four bilinear
  // weights v00+v10+v01+v11 == 1 by construction.
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

/**
 * @brief Verifies integer coordinates collapse AntiAlias to a single full-weight
 *        tap at the input pixel.
 */
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

/**
 * @brief Verifies a sub-pixel sample across the theta=0 seam wraps its X taps
 *        onto columns W-1 and 0 rather than collapsing onto one unwrapped column.
 */
inline void test_antialias_seam_wraps_left_column() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // x = -0.3 sits just left of the seam; y at the equator (y_m = 0) isolates
  // the two X taps so only the wrap behaviour is under test.
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

/**
 * @brief Verifies AntiAlias CLIPS a virtual sub-pole sample (y >= H) instead of
 *        clamping it onto the last physical row — the device-only H_OFFSET path.
 * @details On the device, H_OFFSET appends virtual rows below the LED ring so a
 *          sample landing at y >= H must be dropped (the LEDs stop short of the
 *          south pole) rather than stretched onto row H-1 at full weight — the
 *          y0_ok/y1_ok clip in AntiAlias::plot (filter.h). The host/sim build
 *          sets H_OFFSET == 0, so projected geometry never yields y >= H and that
 *          branch is otherwise unreachable under the native asserts (finding
 *          375). Rather than make the filter's compile-time offset injectable,
 *          drive the branch directly with an out-of-range row coordinate: the
 *          clip decision (y0 < H) is identical regardless of how y was produced.
 *          A sub-pole y emits no tap; the bottom physical row keeps its y0 taps
 *          and never emits onto a row >= H (its y1 == H neighbor is clipped, not
 *          clamped back onto H-1).
 */
inline void test_antialias_clips_virtual_subpole_row() {
  constexpr int W = 32, H = 16;
  Filter::Screen::AntiAlias<W, H> aa;

  // Virtual sub-pole sample: y0 = H and y1 = H+1 are both >= H, so every tap is
  // clipped (clamping would have splatted it onto row H-1 instead).
  int subpole_taps = 0;
  aa.plot(10.5f, static_cast<float>(H) + 0.3f, Pixel(40000, 0, 0), 0.0f, 1.0f,
          [&](float, float, const Pixel &, float, float) { ++subpole_taps; });
  HS_EXPECT_EQ(subpole_taps, 0);

  // Bottom physical row (y0 = H-1): its taps survive, but none lands on a row
  // >= H — the y1 = H neighbor is clipped, confirming clip-not-clamp.
  int row_taps = 0;
  int max_row = -1;
  aa.plot(10.5f, static_cast<float>(H - 1) + 0.3f, Pixel(40000, 0, 0), 0.0f,
          1.0f, [&](float, float y, const Pixel &, float, float) {
            ++row_taps;
            int yi = static_cast<int>(y);
            if (yi > max_row)
              max_row = yi;
          });
  HS_EXPECT_TRUE(row_taps > 0);
  HS_EXPECT_EQ(max_row, H - 1);
}

// ============================================================================
// Screen::Blur::plot — kernel passthrough
// ============================================================================

/**
 * @brief Verifies factor=0 makes Blur an identity: a single full-weight tap at
 *        the input pixel.
 */
inline void test_blur_factor_zero_is_identity() {
  constexpr int W = 32, H = 32;
  // center weight c = 1.0, edge/corner weights 0 → single tap.
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

/**
 * @brief Verifies factor=1 gives the full 3x3 Gaussian: all 9 taps land
 *        in-bounds for an interior pixel and their alphas sum to the input alpha.
 */
inline void test_blur_full_kernel_sums_to_alpha() {
  constexpr int W = 32, H = 32;
  // Gaussian weights: center 0.25, edge 0.125, corner 0.0625 (sum to 1).
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

/**
 * @brief Verifies update() rebuilds the kernel: update(0) collapses a full blur
 *        back to identity.
 */
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

/**
 * @brief Verifies the pole-clip renormalization: a full-blur tap on a pole row
 *        drops the off-pole neighbor row (poles don't wrap) yet still deposits
 *        the full input alpha, with no tap landing outside [0, H).
 * @details Mirrors AntiAlias's Y-clip renorm. Without it the surviving 6 taps
 *          sum to 1 - (2*corner + edge) < alpha and the pole rows darken.
 */
inline void test_blur_pole_row_renormalizes() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(1.0f);

  const float in_alpha = 0.6f;
  for (float py : {0.0f, static_cast<float>(H - 1)}) {
    float sum = 0.0f;
    int count = 0;
    bool all_in_bounds = true;
    blur.plot(15.0f, py, Pixel(1, 1, 1), 0.0f, in_alpha,
              [&](float, float yy, const Pixel &, float, float a) {
                sum += a;
                ++count;
                if (yy < 0.0f || yy >= static_cast<float>(H))
                  all_in_bounds = false;
              });
    // Top/bottom row clipped → only the 2 surviving rows (6 taps) emit.
    HS_EXPECT_EQ(count, 6);
    HS_EXPECT_TRUE(all_in_bounds);
    HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
  }
}

// ============================================================================
// Pixel::ChromaticShift::plot — channel-split fan-out
// ============================================================================

/**
 * @brief Verifies ChromaticShift emits the original colour plus three
 *        single-channel copies shifted to x+1/x+2/x+3 (R, G, B respectively).
 */
inline void test_chromatic_shift_fanout() {
  constexpr int W = 64;
  Filter::Pixel::ChromaticShift<W> cs;

  /**
   * @brief One recorded fan-out tap: position, colour, and alpha.
   */
  struct Tap {
    float x, y;       /**< Emitted pixel coordinate. */
    Pixel c;          /**< Emitted colour. */
    float alpha;      /**< Emitted alpha. */
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

/**
 * @brief Verifies style() returns a live reference to the bound Style (reads,
 *        mutation, and the const overload all alias the original); set_enabled
 *        is callable.
 */
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

/**
 * @brief Verifies Feedback::plot() forwards its input unchanged: one tap with
 *        the original coord and alpha, no Canvas touched.
 * @details The feedback effect lives in flush, not plot.
 */
inline void test_feedback_plot_is_passthrough() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

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
// World filters — direct plot() coverage (no Canvas; capture the PassFn3D taps)
// ============================================================================

/**
 * @brief A small recorder for the (vector, color, age, alpha) tuples a World
 *        filter emits.
 */
struct Tap3D {
  Vector v;            /**< Emitted world-space position. */
  Pixel c;             /**< Emitted colour. */
  float age, alpha;    /**< Emitted age and alpha. */
};

/**
 * @brief Verifies Hole masks a spherical cap: outside the radius the point
 *        passes through untouched; inside, the colour is scaled by
 *        quintic_kernel(d/r) so the very center is fully extinguished.
 */
inline void test_world_hole_masks_cap() {
  constexpr int W = 32;
  Filter::World::Hole<W> hole(Vector(0, 1, 0), 0.5f); // cap at +Y, radius 0.5 rad

  // Far point (south pole) is well outside -> verbatim passthrough.
  int n = 0;
  Tap3D got{};
  hole.plot(Vector(0, -1, 0), Pixel(10000, 20000, 30000), 4.0f, 0.8f,
            [&](const Vector &v, const Pixel &c, float age, float a) {
              got = {v, c, age, a};
              ++n;
            });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_EQ((int)got.c.r, 10000);
  HS_EXPECT_EQ((int)got.c.g, 20000);
  HS_EXPECT_NEAR(got.age, 4.0f, 1e-6f);
  HS_EXPECT_NEAR(got.alpha, 0.8f, 1e-6f);
  HS_EXPECT_NEAR(got.v.y, -1.0f, 1e-6f);

  // Exact center: d=0 -> quintic_kernel(0)=0 -> color fully masked to black.
  hole.plot(Vector(0, 1, 0), Pixel(10000, 20000, 30000), 0.0f, 1.0f,
            [&](const Vector &, const Pixel &c, float, float) { got.c = c; });
  HS_EXPECT_EQ((int)got.c.r, 0);
  HS_EXPECT_EQ((int)got.c.g, 0);
  HS_EXPECT_EQ((int)got.c.b, 0);

  // Half-radius (d = 0.25 rad): scaled by quintic_kernel(0.5) = 0.5.
  Vector half(sinf(0.25f), cosf(0.25f), 0.0f); // 0.25 rad from +Y
  hole.plot(half, Pixel(10000, 20000, 30000), 0.0f, 1.0f,
            [&](const Vector &, const Pixel &c, float, float) { got.c = c; });
  HS_EXPECT_NEAR((float)got.c.r, 5000.0f, 50.0f);
  HS_EXPECT_NEAR((float)got.c.g, 10000.0f, 50.0f);
}

/**
 * @brief Verifies Orient rotates by the bound Orientation, and a single-frame
 *        (stationary) orientation emits one tap with age left untouched.
 * @details A lone snapshot is the newest sub-position (t = 1), so the (1 - t)
 *          age offset is zero. A static orientation must not drift the temporal
 *          channel frame over frame.
 */
inline void test_world_orient_rotates_and_offsets_age() {
  constexpr int W = 32;
  Quaternion q = make_rotation(Y_AXIS, PI_F / 2); // 90 deg about +Y
  Orientation<> ori(q);
  Filter::World::Orient<W> orient(ori);

  int n = 0;
  Tap3D got{};
  orient.plot(X_AXIS, Pixel(1, 2, 3), 5.0f, 1.0f,
              [&](const Vector &v, const Pixel &c, float age, float a) {
                got = {v, c, age, a};
                ++n;
              });
  HS_EXPECT_EQ(n, 1);
  Vector expected = rotate(X_AXIS, q);
  HS_EXPECT_NEAR(got.v.x, expected.x, 1e-4f);
  HS_EXPECT_NEAR(got.v.y, expected.y, 1e-4f);
  HS_EXPECT_NEAR(got.v.z, expected.z, 1e-4f);
  HS_EXPECT_NEAR(got.age, 5.0f, 1e-4f); // age + (1 - t), t = 1 (age-neutral)
}

/**
 * @brief Verifies that with a 3-frame history the tween sweeps the trailing
 *        sub-positions (i = 1..n-1) and spreads age across one frame.
 * @details Age offsets are (1 - t) for t in {0.5, 1.0}.
 */
inline void test_world_orient_motion_blur_sweep_ages() {
  constexpr int W = 32;
  Orientation<> ori; // identity, 1 frame
  ori.push(make_rotation(Y_AXIS, PI_F / 4));
  ori.push(make_rotation(Y_AXIS, PI_F / 2)); // now 3 frames
  Filter::World::Orient<W> orient(ori);

  int n = 0;
  float ages[4] = {0};
  orient.plot(X_AXIS, Pixel(1, 1, 1), 10.0f, 1.0f,
              [&](const Vector &, const Pixel &, float age, float) {
                if (n < 4) ages[n] = age;
                ++n;
              });
  // tween emits len-1 = 2 taps (i=1,2), t in {0.5, 1.0} -> age + {0.5, 0.0}.
  HS_EXPECT_EQ(n, 2);
  HS_EXPECT_NEAR(ages[0], 10.5f, 1e-4f);
  HS_EXPECT_NEAR(ages[1], 10.0f, 1e-4f);
}

/**
 * @brief Verifies OrientSlice picks an orientation from a list by the point's
 *        projection onto an axis, then rotates by it.
 * @details With two distinct orientations a point near +axis selects the last,
 *          near -axis selects the first; disabled is a passthrough.
 */
inline void test_world_orient_slice_selects_by_projection() {
  constexpr int W = 32;
  Orientation<> oris[2];
  oris[0].set(make_rotation(X_AXIS, PI_F / 2)); // index 0
  oris[1].set(make_rotation(Z_AXIS, PI_F / 2)); // index 1
  std::span<const Orientation<>> span(oris, 2);
  Filter::World::OrientSlice<W> slice(span, Y_AXIS);

  auto first_tap = [&](const Vector &probe) {
    Vector out{};
    slice.plot(probe, Pixel(1, 1, 1), 0.0f, 1.0f,
               [&](const Vector &v, const Pixel &, float, float) { out = v; });
    return out;
  };

  // Probe near +Y (projection ~ +1 -> t ~ 1 -> last index 1 = Z rotation).
  Vector near_pos = Vector(0.15f, 0.98f, 0.0f).normalized();
  Vector exp_pos = rotate(near_pos, oris[1].get());
  Vector got_pos = first_tap(near_pos);
  HS_EXPECT_NEAR(got_pos.x, exp_pos.x, 1e-3f);
  HS_EXPECT_NEAR(got_pos.y, exp_pos.y, 1e-3f);
  HS_EXPECT_NEAR(got_pos.z, exp_pos.z, 1e-3f);

  // Probe near -Y (projection ~ -1 -> t ~ 0 -> first index 0 = X rotation).
  Vector near_neg = Vector(0.15f, -0.98f, 0.0f).normalized();
  Vector exp_neg = rotate(near_neg, oris[0].get());
  Vector got_neg = first_tap(near_neg);
  HS_EXPECT_NEAR(got_neg.x, exp_neg.x, 1e-3f);
  HS_EXPECT_NEAR(got_neg.y, exp_neg.y, 1e-3f);

  // Disabled -> verbatim passthrough (no rotation).
  slice.enabled = false;
  Vector pass = first_tap(near_pos);
  HS_EXPECT_NEAR(pass.x, near_pos.x, 1e-6f);
  HS_EXPECT_NEAR(pass.y, near_pos.y, 1e-6f);
  HS_EXPECT_NEAR(pass.z, near_pos.z, 1e-6f);
}

/**
 * @brief Verifies VertexReplicate fans a point onto N copies via rotations from
 *        vertices[0] to each vertex, every copy carrying the same source age.
 * @details Replication is spatial, not temporal. Probing with vertices[0] maps
 *          copy i back to vertices[i] exactly.
 */
inline void test_world_vertex_replicate_fanout_and_age() {
  constexpr int W = 32, N = 3;
  std::array<Vector, N> verts = {X_AXIS, Y_AXIS, Z_AXIS};
  Filter::World::VertexReplicate<W, N> vr(verts);

  Tap3D taps[N];
  int n = 0;
  vr.plot(X_AXIS, Pixel(1, 1, 1), 7.0f, 1.0f,
          [&](const Vector &v, const Pixel &c, float age, float a) {
            if (n < N) taps[n] = {v, c, age, a};
            ++n;
          });
  HS_EXPECT_EQ(n, N);
  // rotate(vertices[0], make_rotation(v0, v_i)) == v_i.
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_NEAR(taps[i].v.x, verts[i].x, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.y, verts[i].y, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.z, verts[i].z, 1e-4f);
    // Age is unchanged for every copy (no per-copy +i offset).
    HS_EXPECT_NEAR(taps[i].age, 7.0f, 1e-6f);
  }
}

/**
 * @brief Verifies Mobius warps via stereographic -> Mobius -> inverse
 *        stereographic, with the default identity map and a non-identity map.
 * @details The default MobiusParams is the identity (a=1,b=0,c=0,d=1), so an
 *          interior point round-trips back to itself; a non-identity map
 *          actually moves it.
 */
inline void test_world_mobius_identity_and_transform() {
  constexpr int W = 32;
  MobiusParams identity; // a=1,b=0,c=0,d=1
  Filter::World::Mobius<W> mob(identity);

  const Vector v = Vector(0.4f, 0.3f, 0.86f).normalized();
  Vector out{};
  int n = 0;
  mob.plot(v, Pixel(1, 2, 3), 2.0f, 0.5f,
           [&](const Vector &o, const Pixel &, float, float) {
             out = o;
             ++n;
           });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(out.x, v.x, 1e-3f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-3f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-3f);

  // A translation map f(z) = z + 1 moves the point and keeps it on the sphere.
  MobiusParams shift(1, 0, 1, 0, 0, 0, 1, 0); // a=1, b=1, c=0, d=1
  Filter::World::Mobius<W> mob2(shift);
  Vector out2{};
  mob2.plot(v, Pixel(1, 1, 1), 0.0f, 1.0f,
            [&](const Vector &o, const Pixel &, float, float) { out2 = o; });
  HS_EXPECT_NEAR(out2.length(), 1.0f, 1e-3f);
  HS_EXPECT_GT(distance_between(out2, v), 0.05f); // actually transformed
}

// ============================================================================
// End-to-end pipeline routing through a live Canvas
//
// The test_canvas/test_scan pattern: a Canvas spins in its ctor while
// !buffer_free(), so every drawn frame MUST advance_display() before the next
// Canvas is constructed. With that, the full Pipeline::plot/flush routing
// (which the per-call helper tests above cannot reach) is exercisable.
// ============================================================================

/**
 * @brief Minimal concrete Effect for binding a live Canvas.
 */
struct PipeFx : public Effect {
  /**
   * @brief Constructs the effect at the given framebuffer dimensions.
   * @param W Framebuffer width in pixels.
   * @param H Framebuffer height in pixels.
   */
  PipeFx(int W, int H) : Effect(W, H) {}
  /**
   * @brief Draws one frame (no-op; the tests populate the Canvas directly).
   */
  void draw_frame() override {}
  /**
   * @brief Reports whether the effect paints a background.
   * @return Always false (tests rely on a black backdrop).
   */
  bool show_bg() const override { return false; }
};

/**
 * @brief Tests whether a pixel is exactly black (all channels zero).
 * @param p Pixel to test.
 * @return True when r, g, and b are all zero.
 */
inline bool px_black(const Pixel &p) {
  return p.r == 0 && p.g == 0 && p.b == 0;
}

/**
 * @brief Tests whether a pixel matches the given RGB triple exactly.
 * @param p Pixel to test.
 * @param r Expected red channel value.
 * @param g Expected green channel value.
 * @param b Expected blue channel value.
 * @return True when every channel equals its expected value.
 */
inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

/**
 * @brief Counts the non-black pixels in the effect's framebuffer.
 * @param fx Effect whose framebuffer is scanned.
 * @return Number of lit (non-black) pixels.
 */
inline size_t count_lit(const PipeFx &fx) {
  size_t n = 0;
  for (int y = 0; y < fx.height(); ++y)
    for (int x = 0; x < fx.width(); ++x)
      if (!px_black(fx.get_pixel(x, y))) ++n;
  return n;
}

/**
 * @brief Verifies the bare 2D sink: int + float overloads, exact (alpha=1)
 *        write, x-wrap, and clip.
 */
inline void test_pipeline_sink_2d_plot_blends_wraps_clips() {
  constexpr int W = 16, H = 8;
  Pipeline<W, H> pipe;
  // fx and fx2 alias the same static double buffer, so only one may be live at
  // a time (Effect's single-live guard) — scope the first frame's effect closed
  // before constructing the second.
  {
    PipeFx fx(W, H);
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
  }

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

/**
 * @brief Verifies the 3D sink overload routes a unit vector via vector_to_pixel
 *        and writes it.
 */
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

/**
 * @brief Verifies a World filter (3D) head fans out through the 3D->2D sink
 *        conversion.
 */
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

/**
 * @brief Verifies a 2D coordinate into a 3D-headed pipeline round-trips through
 *        the float-plot mismatch branch.
 * @details The mismatch branch maps pixel_to_vector, the identity filter passes
 *          it, and the sink maps it back.
 */
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
        // OR-accumulate: the intent is "at least one lit pixel is within a pixel
        // of the input". Plain assignment would record only the last lit pixel
        // in scan order, so the check would stay green even if a regression lit
        // a far-away pixel (caught only by the separate count_lit == 1 above).
        near_input |= (std::abs(x - px) <= 1 && std::abs(y - py) <= 1);
  HS_EXPECT_TRUE(near_input);
}

/**
 * @brief Verifies a Screen filter (2D) head forwards to the sink.
 * @details An integer coordinate collapses AntiAlias to a single tap, so the
 *          routed result is one exact pixel.
 */
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

/**
 * @brief Verifies the Pixel::Feedback::flush warp-field path.
 * @details A Smoke style with no bound NoiseParams makes noise_warp an identity
 *          map, so the flush blends the previous frame back at the same location,
 *          faded by style.fade — deterministic.
 */
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

/**
 * @brief Verifies flush() honors the segment clip like every other rasterizer.
 * @details On segmented hardware each board owns a Y-band, and a feedback flush
 *          that iterated the full canvas would composite the whole sphere into
 *          every board's buffer (wrong output + wasted work). Here the prev frame
 *          is bright at every row but the clip restricts rendering to a sub-band;
 *          rows outside the margin-expanded render band must stay untouched.
 */
inline void test_feedback_flush_respects_clip() {
  constexpr int W = 32, H = 16; // both divisible by Smoke's downsample (4)
  PipeFx fx(W, H);
  ::Feedback::Style style = ::Feedback::Style::Smoke(); // identity warp
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  auto trail = [](float, float, float) { return Color4(Pixel(0, 0, 0), 0.0f); };

  // Frame 1: bright across the WHOLE canvas becomes the "previous" frame.
  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(40000, 40000, 40000);
  }
  fx.advance_display();

  // Restrict this board to the y-band [8,12). With the default margin of 1 the
  // render band is [7,13); rows outside it must not be written.
  fx.set_clip(8, 12, 0, W);

  // Frame 2: empty buffer; clipped flush warps+blends only within the band.
  {
    Canvas c(fx);
    pipe.flush(c, ScreenTrailFn(trail), 1.0f);
  }
  fx.advance_display();

  // Inside the band: carried over and faded.
  HS_EXPECT_FALSE(px_black(fx.get_pixel(W / 2, 9)));
  // Outside the render band: untouched despite the prev frame being lit there.
  HS_EXPECT_TRUE(px_black(fx.get_pixel(W / 2, 0)));
  HS_EXPECT_TRUE(px_black(fx.get_pixel(W / 2, 5)));
  HS_EXPECT_TRUE(px_black(fx.get_pixel(W / 2, 15)));
}

// ============================================================================
// World::Trails — int16 quantization round-trip + ring buffer / ttl lifecycle
//
// Trails are history-bearing filters: plot() encodes the world position into a
// quantized int16 ring entry, flush() ages, evicts, decodes and re-emits. These
// tests cover the quantization round-trip and the ring mechanics directly.
// ============================================================================

/**
 * @brief Verifies plot() passes the frame through and stores it, and flush()
 *        ages, decodes the int16 entry and re-emits it within the quantization
 *        error bound (< 1/32767).
 */
inline void test_world_trails_int16_quantization_roundtrip() {
  constexpr int W = 32, Cap = 8;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<W, Cap> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  const Vector v0 = Vector(0.3f, -0.6f, 0.74f).normalized();

  // plot() passes the current frame through and (age=0 -> ttl=10>0) stores it.
  int passthru = 0;
  trails.plot(v0, Pixel(100, 100, 100), 0.0f, 1.0f,
              [&](const Vector &, const Pixel &, float, float) { ++passthru; });
  HS_EXPECT_EQ(passthru, 1);
  HS_EXPECT_EQ(trails.size(), (size_t)1);

  // flush() ages (10->9, still alive), decodes the int16 entry and re-emits it.
  Vector decoded(0, 0, 0);
  int emitted = 0;
  auto trail = [](const Vector &, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(WorldTrailFn(trail), 1.0f,
               [&](const Vector &v, const Pixel &, float, float) {
                 decoded = v;
                 ++emitted;
               });
  HS_EXPECT_EQ(emitted, 1);
  // Quantization is v*32767 truncated to int16 then *(1/32767): |err| < 1/32767.
  HS_EXPECT_NEAR(decoded.x, v0.x, 1e-4f);
  HS_EXPECT_NEAR(decoded.y, v0.y, 1e-4f);
  HS_EXPECT_NEAR(decoded.z, v0.z, 1e-4f);
}

/**
 * @brief Verifies a component pushed past the unit cube saturates rather than
 *        overflowing int16 and wrapping to a garbage point on the far side of
 *        the sphere.
 * @details encode() clamps to [-1, 1] before quantizing.
 */
inline void test_world_trails_clamps_out_of_range() {
  constexpr int W = 32, Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<W, Cap> trails(/*lifetime=*/10);
  trails.init_storage(arena);

  // x = 1.8 > 1: 1.8*32767 = 58980 would overflow int16 and wrap to ~-0.2 on
  // decode without the clamp. z = -1.5 likewise. Both must saturate to +/-1.
  const Vector v = Vector(1.8f, 0.5f, -1.5f);
  trails.plot(v, Pixel(1, 1, 1), 0.0f, 1.0f,
              [](const Vector &, const Pixel &, float, float) {});

  Vector decoded(0, 0, 0);
  auto trail = [](const Vector &, float) {
    return Color4(Pixel(60000, 60000, 60000), 1.0f);
  };
  trails.flush(WorldTrailFn(trail), 1.0f,
               [&](const Vector &d, const Pixel &, float, float) {
                 decoded = d;
               });
  HS_EXPECT_NEAR(decoded.x, 1.0f, 1e-3f);   // saturated, not wrapped negative
  HS_EXPECT_NEAR(decoded.y, 0.5f, 1e-3f);   // in range: untouched
  HS_EXPECT_NEAR(decoded.z, -1.0f, 1e-3f);  // saturated
}

/**
 * @brief Verifies the ring is a hard-bounded buffer: pushing past Cap evicts the
 *        oldest entries so size() saturates at Cap.
 */
inline void test_world_trails_ring_evicts_oldest() {
  constexpr int W = 32, Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<W, Cap> trails(/*lifetime=*/100);
  trails.init_storage(arena);

  auto noop = [](const Vector &, const Pixel &, float, float) {};
  for (int i = 0; i < Cap + 3; ++i) {
    Vector v = Vector(static_cast<float>(i + 1), 1.0f, 0.5f).normalized();
    trails.plot(v, Pixel(1, 1, 1), 0.0f, 1.0f, noop);
  }
  // Capacity is a hard ring bound; the oldest entries were dropped on overflow.
  HS_EXPECT_EQ(trails.size(), (size_t)Cap);
}

/**
 * @brief Verifies each flush decrements an entry's ttl and the entry is popped
 *        once ttl reaches 0.
 */
inline void test_world_trails_ttl_expiry() {
  constexpr int W = 32, Cap = 4;
  static uint8_t buf[Cap * 16];
  Arena arena(buf, sizeof(buf));
  Filter::World::Trails<W, Cap> trails(/*lifetime=*/2);
  trails.init_storage(arena);

  auto noop = [](const Vector &, const Pixel &, float, float) {};
  trails.plot(Vector(0, 1, 0), Pixel(1, 1, 1), 0.0f, 1.0f, noop); // ttl = 2
  HS_EXPECT_EQ(trails.size(), (size_t)1);

  auto trail = [](const Vector &, float) {
    return Color4(Pixel(1, 1, 1), 1.0f);
  };
  auto sink = [](const Vector &, const Pixel &, float, float) {};
  trails.flush(WorldTrailFn(trail), 1.0f, sink); // ttl 2 -> 1, alive
  HS_EXPECT_EQ(trails.size(), (size_t)1);
  trails.flush(WorldTrailFn(trail), 1.0f, sink); // ttl 1 -> 0, popped
  HS_EXPECT_EQ(trails.size(), (size_t)0);
}

/**
 * @brief Verifies the Screen::Trails store / emit / decay lifecycle.
 * @details Screen::Trails stores float DecayPixels with no int16 quantization
 *          (that path is World::Trails-specific).
 */
inline void test_screen_trails_store_emit_decay() {
  constexpr int W = 32, MAXP = 16;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<W, MAXP> trails(/*lifetime=*/3);
  trails.init_storage(arena);

  PipeFx fx(W, 8); // flush takes a Canvas& (unused by Screen::Trails)
  Canvas c(fx);

  // age=1 (0<age<lifetime): forwarded this frame AND stored. The forward
  // mirrors World::Trails — every live point paints the current frame.
  int passthru = 0;
  float fwd_age = -1.0f;
  trails.plot(10.0f, 4.0f, Pixel(5, 6, 7), 1.0f, 1.0f,
              [&](float, float, const Pixel &, float a, float) {
                fwd_age = a;
                ++passthru;
              });
  HS_EXPECT_EQ(passthru, 1);
  HS_EXPECT_NEAR(fwd_age, 1.0f, 1e-6f); // forwarded verbatim

  auto trail = [](float, float, float) {
    return Color4(Pixel(9, 9, 9), 1.0f);
  };
  // Stored ttl = lifetime - age = 2. Each flush emits then decays (--ttl).
  int emitted = 0;
  auto counting_sink = [&](float, float, const Pixel &, float, float) {
    ++emitted;
  };
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // emit, ttl 2->1
  HS_EXPECT_EQ(emitted, 1);
  emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // emit, ttl 1->0
  HS_EXPECT_EQ(emitted, 1);
  emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f, counting_sink); // decayed out
  HS_EXPECT_EQ(emitted, 0);
}

/**
 * @brief Verifies Screen::Trails forwards an already-aged emission, matching
 *        World::Trails.
 * @details A point with 0 < age < lifetime is both passed through the current
 *          frame and seeded into storage. A point at/past lifetime is still
 *          forwarded, but ttl<=0 keeps it out of storage.
 */
inline void test_screen_trails_forwards_aged_emission() {
  constexpr int W = 32, MAXP = 16;
  static uint8_t buf[MAXP * 32];
  Arena arena(buf, sizeof(buf));
  Filter::Screen::Trails<W, MAXP> trails(/*lifetime=*/5);
  trails.init_storage(arena);

  PipeFx fx(W, 8);
  Canvas c(fx);
  auto trail = [](float, float, float) { return Color4(Pixel(9, 9, 9), 1.0f); };

  // 0 < age < lifetime: forwarded this frame and seeded.
  int fwd = 0;
  trails.plot(3.0f, 4.0f, Pixel(1, 2, 3), 2.0f, 1.0f,
              [&](float, float, const Pixel &, float, float) { ++fwd; });
  HS_EXPECT_EQ(fwd, 1);

  // age == lifetime: dead point is still forwarded, but ttl<=0 so not seeded.
  fwd = 0;
  trails.plot(7.0f, 4.0f, Pixel(1, 2, 3), 5.0f, 1.0f,
              [&](float, float, const Pixel &, float, float) { ++fwd; });
  HS_EXPECT_EQ(fwd, 1);

  // Only the live (age=2) point was seeded -> exactly one flush emission.
  int emitted = 0;
  trails.flush(c, ScreenTrailFn(trail), 1.0f,
               [&](float, float, const Pixel &, float, float) { ++emitted; });
  HS_EXPECT_EQ(emitted, 1);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every filter test case under the "filter" module scope.
 * @return The module's failure count.
 */
inline int run_filter_tests() {
  auto scope = hs_test::begin_module("filter");

  test_trait_member_values();
  test_filter_trait_inheritance();

  test_pipeline_sink_is_2d();
  test_pipeline_get_returns_correct_filter();

  test_antialias_weights_partition();
  test_antialias_integer_coord_single_tap();
  test_antialias_seam_wraps_left_column();
  test_antialias_clips_virtual_subpole_row();

  test_blur_factor_zero_is_identity();
  test_blur_full_kernel_sums_to_alpha();
  test_blur_update_changes_kernel();
  test_blur_pole_row_renormalizes();

  test_chromatic_shift_fanout();

  test_feedback_style_binding();
  test_feedback_plot_is_passthrough();

  test_world_hole_masks_cap();
  test_world_orient_rotates_and_offsets_age();
  test_world_orient_motion_blur_sweep_ages();
  test_world_orient_slice_selects_by_projection();
  test_world_vertex_replicate_fanout_and_age();
  test_world_mobius_identity_and_transform();

  test_pipeline_sink_2d_plot_blends_wraps_clips();
  test_pipeline_sink_3d_plot_routes_to_canvas();
  test_pipeline_world_replicate_fans_out();
  test_pipeline_2d_into_3d_head_roundtrips();
  test_pipeline_screen_antialias_routes_to_sink();
  test_feedback_flush_blends_prev_frame();
  test_feedback_flush_respects_clip();

  test_world_trails_int16_quantization_roundtrip();
  test_world_trails_clamps_out_of_range();
  test_world_trails_ring_evicts_oldest();
  test_world_trails_ttl_expiry();
  test_screen_trails_store_emit_decay();
  test_screen_trails_forwards_aged_emission();

  return hs_test::end_module(scope);
}

} // namespace filter_tests
} // namespace hs_test

