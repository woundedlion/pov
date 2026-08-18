/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
#include "core/render/sdf/volume.h"
#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/canvas.h"
#include "core/math/geometry.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cfloat>
#include <vector>

namespace hs_test {
namespace scan_tests {

/**
 * @brief Pins pole_lod_aggressiveness for a scope and restores it on exit.
 * @details The knob is process-global, so an early return between a hand-rolled
 * save and restore would hand every later case in the module a different LOD
 * regime.
 */
struct ScopedPoleLod {
  float saved; /**< Value in force at construction. */
  /**
   * @brief Pins the knob.
   * @param value Aggressiveness held for the scope.
   */
  explicit ScopedPoleLod(float value) : saved(pole_lod_aggressiveness) {
    pole_lod_aggressiveness = value;
  }
  ~ScopedPoleLod() { pole_lod_aggressiveness = saved; }
  ScopedPoleLod(const ScopedPoleLod &) = delete;
  ScopedPoleLod &operator=(const ScopedPoleLod &) = delete;
};

// ============================================================================
// Scan::Shader::draw — full-sphere per-pixel shader
// ============================================================================

/**
 * @brief Verifies a constant-color shader fills every pixel of the full sphere.
 */
inline void test_shader_constant_fills_canvas() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    Scan::Shader::draw<W, H, 1>(c, [](const Vector &) {
      return Color4(Pixel(40000, 20000, 10000), 1.0f);
    });
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
  hs_test::StubEffect fx(W, H);
  {
    Canvas c(fx);
    // Opacity keys on the sub-sample's position, not call order: the 2x2 grid's
    // +0.25/-0.25 px x-offsets land at fractional theta-grid phase 0.25 vs 0.75,
    // so two of the four samples per pixel are opaque regardless of iteration.
    Scan::Shader::draw<W, H, 4>(c, [](const Vector &v) -> Color4 {
      float theta = std::atan2(v.z, v.x);
      if (theta < 0.0f)
        theta += 2.0f * PI_F;
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
  hs_test::StubEffect fx(W, H);
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
  hs_test::StubEffect fx(W, H);
  fx.set_clip(5, 10, 0, W); // rows [5,10)
  fx.set_margin(0);         // no render-margin expansion

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

/**
 * @brief Verifies every Shader entry point paints exactly the columns the clip
 *        arc admits, with the values an unclipped draw produced.
 * @details The draws walk the arc's pieces instead of testing every column, so
 * pin the painted output against the unclipped render filtered by
 * XClip::clipped — under a plain arc and under one whose margin wraps the band
 * across the seam.
 */
inline void test_shader_clip_arc_matches_predicate() {
  constexpr int W = 32, H = 16;

  auto positional = [](const Vector &v) {
    return Color4(Pixel(static_cast<uint16_t>((v.x * 0.5f + 0.5f) * 60000.0f),
                        static_cast<uint16_t>((v.z * 0.5f + 0.5f) * 60000.0f),
                        30000),
                  1.0f);
  };

  // 0/1: single-callback draw at 1x and 4x. 2: split vertex/fragment draw.
  // 3: draw_grid.
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
          c, [&](const Vector &v, Fragment &f) { f.color = positional(v); },
          [](Fragment &f) { f.v0 = 1.0f; });
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
    out.resize(W * H);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out[y * W + x] = fx.get_pixel(x, y);
  };

  struct Band {
    int x0, x1, margin;
  };
  // A plain arc, and one whose margin underflows column 0 into a wrapped band.
  const Band bands[] = {{8, 20, 2}, {0, 10, 3}};

  for (int variant = 0; variant < 4; ++variant) {
    HS_CONTEXT("variant", variant, 0);
    // One live Effect at a time: read the unclipped render back before the
    // clipped fixture exists.
    std::vector<Pixel> reference;
    render(variant, 0, W, 0, reference);

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
  hs_test::StubEffect fx(W, H);
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
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Ring::draw<W, H, false>(
        pipe, c, basis, radius, thickness, [](const Vector &, Fragment &f) {
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
  hs_test::StubEffect fx(W, H);
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
      hs_test::StubEffect legacy(W, H);
      if (partial_clip) {
        legacy.set_clip(9, 53, 17, 81);
        legacy.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(legacy);
        Scan::DistortedRing::draw<W, H>(pipeline, canvas, basis, 1.7f, 0.12f,
                                        knots, LUT_N, shader);
      }
      legacy.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          expected[y * W + x] = legacy.get_pixel(x, y);
    }

    hs_test::StubEffect flat(W, H);
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

    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = flat.get_pixel(x, y);
        if (!is_black(a))
          ++lit;
        HS_EXPECT_NEAR(static_cast<int>(a.r), static_cast<int>(b.r), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.g), static_cast<int>(b.g), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.b), static_cast<int>(b.b), 1);
      }
    }
    // Both paths drawing nothing would satisfy every comparison above.
    HS_EXPECT_GT(lit, (size_t)0);
  };

  check(false);
  check(true);
}

// AA-dust envelope of the fused ring-group scan. Both paths currently agree
// bit-exactly on all three cases (the per-case line reports the live figure),
// so these bound the divergence the doc below permits rather than one that is
// spent: a handful of stroke-edge pixels at 1/128 of full scale.
constexpr int GROUP_MAX_CHANNEL_DELTA = 512;
constexpr int GROUP_MAX_DIFF_PIXELS = 16;

/**
 * @brief Verifies RingGroup::draw matches per-ring sequential rasterizes
 *        up to interval-clip AA dust.
 * @details Four near-coincident rings with distinct colors, alphas, and
 * thicknesses, drawn sequentially vs as one fused group. Per-pixel blend
 * order is slot order in both paths; the only permitted divergence is the
 * handful of low-alpha stroke-edge pixels that a ring's own fast-path
 * interval clip drops but the group's covering-ring scan paints (the group
 * is the truer SDF coverage), so mismatches must be both rare and small.
 * Covered: full frame, a partial clip with an x band, and a near-pole axis
 * that forces the group's full-row-scan fallback. The claim is scoped to
 * pole_lod_aggressiveness 0, which the test pins: the fused scan shades every
 * column while the per-ring path decimates near-pole rows.
 */
inline void test_ring_group_matches_sequential() {
  constexpr int W = 96, H = 64;
  constexpr int N = 4;
  const ScopedPoleLod lod(0.0f);

  auto run_case = [&](const Vector &normal, bool partial_clip) {
    Basis bases[N];
    const float ths[N] = {0.08f, 0.04f, 0.04f, 0.08f};
    const Color4 colors[N] = {Color4(Pixel(60000, 10000, 5000), 0.9f),
                              Color4(Pixel(5000, 60000, 10000), 0.6f),
                              Color4(Pixel(10000, 5000, 60000), 0.4f),
                              Color4(Pixel(30000, 30000, 30000), 0.7f)};
    alignas(SDF::Ring) unsigned char mem[N * sizeof(SDF::Ring)];
    auto *shapes = reinterpret_cast<SDF::Ring *>(mem);
    for (int s = 0; s < N; ++s) {
      Quaternion q =
          make_rotation(Vector(0.2f, 0.5f, 0.8f).normalized(), 0.02f * s);
      bases[s] = make_basis(q, normal);
      new (&shapes[s]) SDF::Ring(bases[s], 1.0f, ths[s]);
    }

    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect seq(W, H);
      if (partial_clip) {
        seq.set_clip(9, 53, 17, 81);
        seq.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(seq);
        for (int s = 0; s < N; ++s) {
          auto shader = [&](const Vector &, Fragment &f) {
            f.color = colors[s];
          };
          Scan::rasterize<W, H, false>(pipeline, canvas, shapes[s], shader);
        }
      }
      seq.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          expected[y * W + x] = seq.get_pixel(x, y);
    }

    hs_test::StubEffect fused(W, H);
    if (partial_clip) {
      fused.set_clip(9, 53, 17, 81);
      fused.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fused);
      Scan::RingGroup::draw<W, H>(
          pipeline, canvas, shapes, N,
          [&](int s, const Vector &, Fragment &f) { f.color = colors[s]; });
    }
    fused.advance_display();

    int diff_px = 0;
    int worst_delta = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = fused.get_pixel(x, y);
        if (a.r == b.r && a.g == b.g && a.b == b.b)
          continue;
        ++diff_px;
        const int deltas[3] = {std::abs(static_cast<int>(a.r) - b.r),
                               std::abs(static_cast<int>(a.g) - b.g),
                               std::abs(static_cast<int>(a.b) - b.b)};
        for (int d : deltas) {
          worst_delta = std::max(worst_delta, d);
          HS_EXPECT_LE(d, GROUP_MAX_CHANNEL_DELTA);
        }
      }
    }
    std::printf("  [ring-group] diff_px=%d worst_delta=%d\n", diff_px,
                worst_delta);
    HS_EXPECT_LE(diff_px, GROUP_MAX_DIFF_PIXELS);
  };

  run_case(Vector(0.3f, 0.8f, -0.5f).normalized(), false);
  run_case(Vector(0.3f, 0.8f, -0.5f).normalized(), true);
  // Near-pole axis: r_val under the horizontal-projection floor forces the
  // group's full-row-scan fallback.
  run_case(Vector(0.005f, 1.0f, 0.0f).normalized(), false);
}

/**
 * @brief Verifies DistortedRingStack::draw matches rasterizing the stack's
 *        rings one by one.
 * @details Five evenly spaced same-axis knot rings with distinct thicknesses,
 * colors, and alphas, drawn sequentially through Scan::DistortedRing::draw with
 * suppress_pole_fill (the per-ring path the stack's doc claims to mirror) vs as
 * one fused stack scan. The shader keys its green channel on the azimuth v0 and
 * its alpha on the coverage v2, so a divergence in either register shows up in
 * the pixels. Covered: full frame, a partial clip with an x band, a culled
 * middle ring (slot_by_ring -1), and a near-pole axis that forces the
 * full-row-scan fallback on both paths. The claim is scoped to
 * pole_lod_aggressiveness 0, which the test pins: the fused scan shades every
 * column while the per-ring path decimates near-pole rows.
 *
 * Which pixels light is exact. Channel values carry a tolerance: the shipping
 * builds are -ffast-math, which inlines the shared frame and coverage
 * arithmetic into two different loops and reassociates each copy on its own, so
 * the two blend weights agree to float rounding rather than bit-exactly. The
 * tolerance is derived, not fitted: blend_alpha quantizes the weight to 16 bits
 * and the reassociation divergence stays far inside one 1/65535 quantum, so
 * each blend's weight shifts by at most one step; lerp16's slope in both the
 * weight and the destination is at most 1 at 16-bit scale, so one blend moves a
 * channel by at most 1 and passes an incoming difference through undamped. At
 * most N_RINGS blends land on a pixel, so N_RINGS bounds the channel delta.
 */
inline void test_distorted_ring_stack_matches_sequential() {
  constexpr int W = 96, H = 64;
  constexpr int N_RINGS = 5, LUT_N = 32;
  const ScopedPoleLod lod(0.0f);

  const float ths[N_RINGS] = {0.07f, 0.05f, 0.09f, 0.05f, 0.07f};
  const Color4 colors[N_RINGS] = {Color4(Pixel(60000, 0, 5000), 0.9f),
                                  Color4(Pixel(5000, 0, 10000), 0.6f),
                                  Color4(Pixel(10000, 0, 60000), 0.4f),
                                  Color4(Pixel(30000, 0, 30000), 0.7f),
                                  Color4(Pixel(45000, 0, 20000), 1.0f)};

  auto run_case = [&](const Vector &normal, bool partial_clip, int culled) {
    Basis basis = make_basis(
        make_rotation(Vector(0.2f, 0.5f, 0.8f).normalized(), 0.3f), normal);

    float knots[N_RINGS][LUT_N + 1];
    for (int i = 0; i < N_RINGS; ++i) {
      for (int k = 0; k < LUT_N; ++k) {
        float t = 2.0f * PI_F * k / LUT_N;
        knots[i][k] = 0.06f * sinf((i + 2) * t) + 0.03f * cosf(3.0f * t + i);
      }
      knots[i][LUT_N] = knots[i][0];
    }
    // Ring i's centerline colatitude must be PI * (i + 1) / (N_RINGS + 1);
    // target_angle is radius * PI/2.
    auto ring_radius = [](int i) { return 2.0f * (i + 1) / (N_RINGS + 1); };

    alignas(SDF::DistortedRing) unsigned char
        mem[N_RINGS * sizeof(SDF::DistortedRing)];
    auto *shapes = reinterpret_cast<SDF::DistortedRing *>(mem);
    SDF::KnotPrefilter prefilters[N_RINGS];
    int8_t slot_by_ring[N_RINGS];
    Color4 slot_color[N_RINGS];
    int n_slots = 0;
    for (int i = 0; i < N_RINGS; ++i) {
      if (i == culled) {
        slot_by_ring[i] = -1;
        continue;
      }
      new (&shapes[n_slots])
          SDF::DistortedRing(basis, ring_radius(i), ths[i], knots[i], LUT_N,
                             0.0f, prefilters[n_slots]);
      slot_color[n_slots] = colors[i];
      slot_by_ring[i] = static_cast<int8_t>(n_slots);
      ++n_slots;
    }

    auto shade = [&](int s, Fragment &f) {
      const uint16_t g = static_cast<uint16_t>(wrap_t(f.v0) * 60000.0f);
      f.color = Color4(Pixel(slot_color[s].color.r, g, slot_color[s].color.b),
                       slot_color[s].alpha * f.v2);
    };

    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect seq(W, H);
      if (partial_clip) {
        seq.set_clip(9, 53, 17, 81);
        seq.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(seq);
        for (int i = 0; i < N_RINGS; ++i) {
          const int s = slot_by_ring[i];
          if (s < 0)
            continue;
          auto shader = [&](const Vector &, Fragment &f) { shade(s, f); };
          Scan::DistortedRing::draw<W, H>(pipeline, canvas, basis,
                                          ring_radius(i), ths[i], knots[i],
                                          LUT_N, shader, 0.0f, false,
                                          /*suppress_pole_fill=*/true);
        }
      }
      seq.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          expected[y * W + x] = seq.get_pixel(x, y);
    }

    hs_test::StubEffect fused(W, H);
    if (partial_clip) {
      fused.set_clip(9, 53, 17, 81);
      fused.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fused);
      Scan::DistortedRingStack::draw<W, H>(
          pipeline, canvas, N_RINGS, shapes, slot_by_ring, n_slots,
          [&](int s, const Vector &, Fragment &f) { shade(s, f); });
    }
    fused.advance_display();

    // ScalarFn's inplace_function member is not trivially destructible.
    for (int s = 0; s < n_slots; ++s)
      shapes[s].~DistortedRing();

    constexpr int CHANNEL_TOL = N_RINGS; // one 16-bit step per composited blend
    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = fused.get_pixel(x, y);
        if (!is_black(a))
          ++lit;
        HS_EXPECT_EQ(is_black(a), is_black(b));
        HS_EXPECT_NEAR(static_cast<int>(a.r), static_cast<int>(b.r),
                       CHANNEL_TOL);
        HS_EXPECT_NEAR(static_cast<int>(a.g), static_cast<int>(b.g),
                       CHANNEL_TOL);
        HS_EXPECT_NEAR(static_cast<int>(a.b), static_cast<int>(b.b),
                       CHANNEL_TOL);
      }
    }
    // Guard against both paths drawing nothing.
    HS_EXPECT_GT(lit, (size_t)200);
  };

  const Vector axis = Vector(0.3f, 0.8f, -0.5f).normalized();
  run_case(axis, false, -1);
  run_case(axis, true, -1);
  run_case(axis, false, 2);
  // Near-pole axis: r_val under the horizontal-projection floor forces the
  // full-row-scan fallback on both paths.
  run_case(Vector(0.005f, 1.0f, 0.0f).normalized(), false, -1);
}

/**
 * @brief Verifies rasterize_face walks the pixels scan_region walks.
 * @details rasterize_face carries its own copy of scan_region's
 * wrap/coalesce/clip run builder. SDF::Face satisfies ScanShape, so the same
 * face drawn through Scan::rasterize and through Scan::rasterize_face must
 * light the same pixels with the same coverage; a constant-color shader makes
 * any divergence in either the run set or the AA band a framebuffer difference.
 * Covered: a mid-latitude face, one whose azimuth wedge straddles theta=0, a
 * pole-touching face (per-row runs plus the full-row fallback), and both a
 * plain and a seam-wrapping x clip. Scoped to pole_lod_aggressiveness 0, which
 * the test pins: only rasterize_face scales its block slack by the face's plane
 * stretch, so the two decimate near-pole rows differently.
 */
inline void test_face_rasterize_matches_scan_region() {
  constexpr int W = 96, H = 64;
  constexpr int HV = H + hs::H_OFFSET;
  const ScopedPoleLod lod(0.0f);

  // x0 == x1 leaves the frame unclipped; a margin that underflows x0 gives the
  // seam-wrapping band.
  auto run_case = [&](const Vector &axis, float rho, int sides, float phase,
                      int x0, int x1, int margin) {
    HS_CONTEXT("face", sides, x0);
    Basis basis = make_basis(Quaternion(), axis);
    Vector verts[8];
    uint16_t idx[8];
    for (int i = 0; i < sides; ++i) {
      float a = (TWO_PI_F * i) / sides + phase;
      verts[i] = (basis.v * cosf(rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                     .normalized();
      idx[i] = static_cast<uint16_t>(i);
    }
    std::span<const Vector> vspan(verts, sides);
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
        auto shader = [&](const Vector &, Fragment &f) { f.color = color; };
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
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          expected[y * W + x] = generic.get_pixel(x, y);
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
  const Vector mid = Vector(0.3f, 0.8f, -0.5f).normalized();
  run_case(mid, 0.45f, 6, 0.37f, 0, 0, 0);
  run_case(mid, 0.45f, 6, 0.37f, 62, 94, 0);
  // Centred on theta=0, so the wedge straddles the seam. The band [90, 96) u
  // [0, 46) that the margin wraps around the seam covers it.
  const Vector seam = Vector(1.0f, 0.15f, 0.0f).normalized();
  run_case(seam, 0.5f, 5, 0.0f, 0, 0, 0);
  run_case(seam, 0.5f, 5, 0.0f, 0, 40, 6);
  // Face over the pole: runs are rebuilt per row and the rows it cannot bound
  // fall back to a full-width scan.
  const Vector pole = Vector(0.02f, 1.0f, 0.0f).normalized();
  run_case(pole, 0.35f, 4, 0.11f, 0, 0, 0);
  run_case(pole, 0.35f, 4, 0.11f, 12, 62, 0);
  run_case(pole, 0.35f, 4, 0.11f, 0, 40, 6);
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
      [&](int wx, int, const Vector &, int run) {
        for (int i = 0; i < run; ++i)
          if (wx + i >= 0 && wx + i < W)
            counts[wx + i]++;
        return run;
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
      [&](int wx, int, const Vector &, int run) {
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
 * @details A run that straddled a stride block would shade from a column the
 * block does not own. Records the column each shade was probed at, not just the
 * block it fell in, and pins it against the offer rule: a column is settled from
 * its own block's anchor while the whole block fits inside the walked run, and
 * walked at its own column otherwise. A clip arc that cuts a block therefore
 * moves the shade source for the columns it leaves behind, so the probe column
 * is what a seam comparison has to read. Drives a near-pole row unclipped, under
 * each quarter clip, and under an arc whose ends fall inside blocks.
 */
inline void test_pole_lod_runs_are_canvas_anchored() {
  constexpr int W = 288;
  constexpr int H = 144;
  const int y = 2; // near the north pole, so sin(phi) is small and stride > 1

  TrigLUT<W, H>::init();

  const ScopedPoleLod lod(1.0f);
  const int lod_stride = Scan::pole_lod_run(TrigLUT<W, H>::sin_phi[y]);
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
        [&](int wx, int, const Vector &, int run) {
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
    ClipRegion::XClip xc{};
    xc.active = true;
    xc.wrap = false;
    xc.rs = q * (W / 4);
    xc.re = xc.rs + (W / 4);
    std::array<int, W> part{};
    check(xc, xc.rs, xc.re, part);
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
  pole_lod_aggressiveness = 0.0f;
  HS_EXPECT_EQ(Scan::pole_lod_run(TrigLUT<W, H>::sin_phi[y]), 1);
  HS_EXPECT_EQ(Scan::pole_lod_run(1.0f), 1);
}

/**
 * @brief Verifies near-pole decimation shades the same pixels as an
 *        undecimated walk.
 * @details A block is settled from one probe only where that probe bounds the
 * whole block, so a constant-color draw lands the same framebuffer at
 * pole_lod_aggressiveness 1.0 as at 0: a cleared block paints nothing either
 * way, and a splatted interior block is at full coverage in every column. A
 * shape that reports FAR_SENTINEL past a zero-margin reject band bounds
 * nothing, so an ungated block test drops whole runs of opaque columns here.
 * Drives the stroke path (Scan::Ring, whose bounding annulus is the stroke band
 * itself), both face rasterizers, and the folded-sector solids down both the
 * shaded and the constant-color scan, across rows whose stride exceeds 1.
 */
inline void test_pole_lod_shading_matches_undecimated() {
  constexpr int W = 96, H = 64;
  constexpr int HV = H + hs::H_OFFSET;
  const ScopedPoleLod scoped_lod(0.0f);
  const Color4 color(Pixel(60000, 45000, 30000), 1.0f);

  auto readback = [](hs_test::StubEffect &fx) {
    std::vector<Pixel> out(W * H);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out[y * W + x] = fx.get_pixel(x, y);
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
  auto draw_ring = [&](float lod, const Vector &axis, float radius) {
    pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Basis basis = make_basis(Quaternion(), axis);
      Scan::Ring::draw<W, H, false>(
          pipe, c, basis, radius, /*thickness=*/0.05f,
          [&](const Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  auto draw_face = [&](float lod, float tilt, float rho, int sides, float phase,
                       bool fused) {
    pole_lod_aggressiveness = lod;
    const Vector axis = Vector(tilt, 1.0f, 0.0f).normalized();
    Basis basis = make_basis(Quaternion(), axis);
    Vector verts[8];
    uint16_t idx[8];
    for (int i = 0; i < sides; ++i) {
      float a = (TWO_PI_F * i) / sides + phase;
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
      SDF::Face face(std::span<const Vector>(verts, sides),
                     std::span<const uint16_t>(idx, sides), scratch, HV, H,
                     &canvas.clip());
      auto shader = [&](const Vector &, Fragment &f) { f.color = color; };
      if (fused)
        Scan::rasterize_face<W, H, Pipeline<W, H>>(pipe, canvas, face, shader);
      else
        Scan::rasterize<W, H, false>(pipe, canvas, face, shader);
    }
    fx.advance_display();
    return readback(fx);
  };

  auto draw_folded = [&](auto tag, float lod, const Vector &axis, float radius,
                         int sides, bool typed) {
    using ShapeT = typename decltype(tag)::type;
    pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Basis basis = make_basis(Quaternion(), axis);
      if (typed)
        ShapeT::template draw_solid<W, H>(pipe, c, basis, radius, sides, color);
      else
        ShapeT::template draw<W, H, false>(
            pipe, c, basis, radius, sides,
            [&](const Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  // A ring carved out of a polygon. Past its stroke the ring reports the
  // sentinel, which loses Subtract's max, so a probe one column from a carved
  // column reports the polygon's own depth and an ungated splat paints the
  // carve shut.
  auto draw_carved = [&](float lod, bool typed) {
    pole_lod_aggressiveness = lod;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      const Basis poly_basis = make_basis(Quaternion(), Vector(0, 1, 0));
      const Vector ring_axis = Vector(0.12f, 1.0f, 0.0f).normalized();
      const Basis ring_basis = make_basis(Quaternion(), ring_axis);
      SDF::PlanarPolygon poly(poly_basis, /*circumradius=*/0.35f, /*sides=*/5,
                              0.0f);
      // Centerline through the canvas pole: the stroke crosses the decimated
      // rows the polygon's interior covers.
      SDF::Ring ring(ring_basis, /*radius=*/0.12f / (PI_F / 2.0f),
                     /*thickness=*/0.05f);
      SDF::Subtract<SDF::PlanarPolygon, SDF::Ring> carved(poly, ring);
      if (typed)
        Scan::rasterize_solid<W, H>(pipe, c, carved, color);
      else
        Scan::rasterize<W, H, false>(
            pipe, c, carved,
            [&](const Vector &, Fragment &f) { f.color = color; });
    }
    fx.advance_display();
    return readback(fx);
  };

  // A sliver, posed so the decimated rows are the ones its projection stretches
  // most: the circumradius reaches 0.7 gnomonic-plane units while the inradius
  // keeps the face on the plane-distance path, whose report runs up to 1 + r^2
  // times an angular step.
  auto draw_sliver_face = [&](float lod, bool fused, float tilt, float spin) {
    pole_lod_aggressiveness = lod;
    const Vector axis = Vector(sinf(tilt), cosf(tilt), 0.0f);
    const Basis basis = make_basis(Quaternion(), axis);
    const float px[3] = {0.7f, -0.35f, -0.35f};
    const float py[3] = {0.0f, 0.06f, -0.06f};
    Vector verts[3];
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
      SDF::Face face(std::span<const Vector>(verts, 3),
                     std::span<const uint16_t>(idx, 3), scratch, HV, H,
                     &canvas.clip());
      HS_EXPECT_TRUE(face.linear_dist);
      HS_EXPECT_GT(face.max_dist_sq, 0.25f);
      auto shader = [&](const Vector &, Fragment &f) { f.color = color; };
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
  pole_lod_aggressiveness = 1.0f;
  HS_EXPECT_GT(Scan::pole_lod_run(TrigLUT<W, H>::sin_phi[2]), 1);

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
        const float spin = (TWO_PI_F * sc.spin_step) / 48.0f;
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

  for (float tilt : {0.12f, 0.3f}) {
    HS_CONTEXT("ring", static_cast<int>(tilt * 100.0f));
    const Vector axis = Vector(tilt, 1.0f, 0.0f).normalized();
    const float radius = tilt / (PI_F / 2.0f);
    compare(draw_ring(0.0f, axis, radius), draw_ring(1.0f, axis, radius));
  }

  // Triangles: the widest gap between inradius and circumradius, so the cull
  // disk clears a vertex by least and a block probe just outside it can sit
  // within the AA fringe's reach of the polygon. Both are small enough to take
  // the linear_dist path, on which distance() reports the plane distance the
  // cull margin is measured in.
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
  const Vector rim_axis = Vector(1.0f, 0.5f, 0.0f).normalized();
  // Flower's petals meet at the antipode of basis.v, so that axis sits by the
  // canvas pole: across it the sign alternates over a vanishing arc, which no
  // finite slack bounds.
  const Vector fold_axis = Vector(0.10f, -1.0f, 0.0f).normalized();
  // 4.0 is inside the live WASM range and lengthens every run, so a slack that
  // undercuts the shape's own stretch mis-shades a fringe column there.
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
        [&](int wx, int, const Vector &, int run) {
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
  hs_test::StubEffect fx(W, H);
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
  float last_alpha =
      -1.0f;     /**< AA alpha from the most recent plot, -1 if none. */
  int count = 0; /**< Number of times plot() was invoked. */
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
  hs_test::StubEffect fx(W, H);
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
      1.0f -
      hs::clamp(expected_result.raw_dist / expected_result.size, 0.0f, 1.0f));

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
  hs_test::StubEffect fx(W, H);
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

  static_assert(!SDF::Union<SDF::Line, SDF::Line>::is_solid,
                "a Union of strokes is not solid");

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
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Basis basis = make_basis(Quaternion(), Y_AXIS); // axis = north pole (+Y)
      Scan::Ring::draw<W, H, false>(pipe, c, basis, radius, /*thickness=*/0.05f,
                                    [](const Vector &, Fragment &f) {
                                      f.color = Color4(
                                          Pixel(60000, 60000, 60000), 1.0f);
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
    struct R {
      float centroid;
      int total;
      int lit[H];
    };
    R r;
    r.centroid = total ? static_cast<float>(weighted) / total : -1.0f;
    r.total = static_cast<int>(total);
    for (int y = 0; y < H; ++y)
      r.lit[y] = lit[y];
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
  hs_test::StubEffect fx(W, H);
  Canvas c(fx);

  const float radius = 0.5f;     // centerline at polar PI/4
  const float thickness = 0.10f; // band half-width in radians
  const float target = radius * (PI_F / 2.0f);
  Basis basis = make_basis(Quaternion(), Y_AXIS);
  SDF::Ring ring(basis, radius, thickness);

  // March outward from the centerline along the az=0 meridian. raw distance to
  // the centerline grows with the offset, so the signed distance crosses 0 at
  // the surface and the alpha should fall 1 -> 0.
  const int N = 12;
  float prev = 1.0f;
  bool saw_one = false, saw_zero = false, saw_mid = false;
  for (int i = 0; i < N; ++i) {
    float delta = (thickness * 1.4f) * i / (N - 1); // 0 .. 1.4*thickness
    float ph = target + delta;
    Vector p(sinf(ph), cosf(ph), 0.0f); // az=0, polar angle ph
    int count = 0;
    float a = scan_alpha_at<W, H>(ring, p, c, &count);
    if (count == 0)
      a = 0.0f; // outside the stroke -> not drawn -> alpha 0

    HS_EXPECT_LE(a, prev + 1e-4f); // monotone non-increasing outward
    prev = a;

    if (a > 0.95f)
      saw_one = true;
    if (a < 0.05f)
      saw_zero = true;
    if (a > 0.1f && a < 0.9f)
      saw_mid = true;
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
template <int W> inline bool row_has_lit(const hs_test::StubEffect &fx, int y) {
  for (int x = 0; x < W; ++x)
    if (!is_black(fx.get_pixel(x, y)))
      return true;
  return false;
}

/**
 * @brief Placement oracle for a filled cap: the pole the shape is centred on is
 *        lit, the opposite pole stays dark, and the lit area matches the
 *        shape's angular extent.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param fx Effect whose canvas was drawn into.
 * @param cap_north True if the shape caps the +Y pole (row 0); false if it caps
 *        the -Y pole (row H-1). Star/PlanarPolygon centre on basis.v (north);
 *        Flower centres on its antipode -basis.v (south).
 * @param r_min Angular radius of the largest cap inscribed in the shape (the
 *        shape's smallest boundary radius, radians).
 * @param r_max Angular radius of the cap circumscribing the shape (radians).
 * @details Drawn to actual pixels — the placement, not just distance()/cull
 *          coverage. With an angular radius well under PI/2 the far pole sits
 *          outside the fill, so a projection sign flip (wrong hemisphere) is
 *          caught here where a bare draw-without-asserting is not.
 *
 *          Rows are uniform in polar angle, so a pole cap of angular radius r
 *          lights W*H*r/PI pixels and the shape's own count falls between its
 *          inscribed and circumscribed caps. The slack is 3 rows either way:
 *          half a row of quantization plus the antialiased fringe, a ~1 pixel
 *          band along a perimeter of at most a couple of canvas widths. A fill
 *          collapsed to a handful of pixels falls under the floor; a flooded
 *          hemisphere (W*H/2) sits over the ceiling for every radius here.
 */
template <int W, int H>
inline void expect_filled_cap(const hs_test::StubEffect &fx, bool cap_north,
                              float r_min, float r_max) {
  const int near_row = cap_north ? 0 : H - 1;
  const int far_row = cap_north ? H - 1 : 0;
  HS_EXPECT_TRUE((row_has_lit<W>(fx, near_row)));
  HS_EXPECT_FALSE((row_has_lit<W>(fx, far_row)));

  const float rows_per_radian = static_cast<float>(H) / PI_F;
  const float slack_rows = 3.0f;
  const int lit = static_cast<int>(count_lit_region<W, H>(fx));
  HS_EXPECT_GE(lit,
               static_cast<int>(W * (r_min * rows_per_radian - slack_rows)));
  HS_EXPECT_LE(lit,
               static_cast<int>(W * (r_max * rows_per_radian + slack_rows)));
}

/** @brief Verifies a filled Star caps the basis.v (+Y) pole and not the other. */
inline void test_star_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Star::draw<W, H, false>(
        pipe, c, basis, RADIUS, /*sides=*/5, [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  // Tips at the circumradius, notches at STAR_INNER_RATIO of it.
  const float r_max = RADIUS * (PI_F / 2.0f);
  expect_filled_cap<W, H>(fx, /*cap_north=*/true, r_max * STAR_INNER_RATIO,
                          r_max);
}

/** @brief Verifies a filled PlanarPolygon caps the basis.v (+Y) pole, not the other. */
inline void test_planar_polygon_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  constexpr int SIDES = 6;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, RADIUS, SIDES, [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  // Vertices at the circumradius, edge midpoints at the apothem.
  const float r_max = RADIUS * (PI_F / 2.0f);
  expect_filled_cap<W, H>(fx, /*cap_north=*/true, r_max * cosf(PI_F / SIDES),
                          r_max);
}

/**
 * @brief Verifies a filled Flower caps the antipode (-Y) pole, not basis.v.
 * @details Flower's SDF scans from antipode = -basis.v, so the fill lands on the
 *          opposite pole from Star/PlanarPolygon — a placement detail only a
 *          pixel oracle pins.
 */
inline void test_flower_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  constexpr int SIDES = 6;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Basis basis = make_basis(Quaternion(), Y_AXIS);
    Scan::Flower::draw<W, H, false>(
        pipe, c, basis, RADIUS, SIDES, [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  // Petal boundary (PI - r) * cos(local) == PI - r_max: the tip at local 0,
  // narrowest at the sector edge.
  const float r_max = RADIUS * (PI_F / 2.0f);
  const float r_min = PI_F - (PI_F - r_max) / cosf(PI_F / SIDES);
  expect_filled_cap<W, H>(fx, /*cap_north=*/false, r_min, r_max);
}

enum class SolidShape { PLANAR_POLYGON, SPHERICAL_POLYGON, FLOWER, STAR };

/**
 * @brief Verifies the typed solid-color path matches generic Scan output.
 */
inline void test_solid_color_path_matches_generic() {
  constexpr int W = 64, H = 48;
  constexpr float RADIUS = 0.72f;
  constexpr int SIDES = 5;
  constexpr float PHASE = 0.31f;
  const Color4 color(Pixel(51000, 23000, 9000), 0.37f);
  const SolidShape shapes[] = {SolidShape::PLANAR_POLYGON,
                               SolidShape::SPHERICAL_POLYGON,
                               SolidShape::FLOWER, SolidShape::STAR};

  for (SolidShape shape : shapes) {
    for (bool debug_bb : {false, true}) {
      std::vector<Pixel> generic_pixels;
      generic_pixels.reserve(W * H);

      {
        hs_test::StubEffect generic_fx(W, H);
        Pipeline<W, H> generic_pipeline;
        {
          Canvas canvas(generic_fx);
          Basis basis = make_basis(
              make_rotation(Vector(0.3f, 0.7f, 0.2f).normalized(), 0.4f),
              Y_AXIS);
          auto shader = [&](const Vector &, Fragment &f) { f.color = color; };
          switch (shape) {
          case SolidShape::PLANAR_POLYGON:
            Scan::PlanarPolygon::draw<W, H, false>(generic_pipeline, canvas,
                                                   basis, RADIUS, SIDES, shader,
                                                   PHASE, debug_bb);
            break;
          case SolidShape::SPHERICAL_POLYGON:
            Scan::SphericalPolygon::draw<W, H, false>(generic_pipeline, canvas,
                                                      basis, RADIUS, SIDES,
                                                      shader, PHASE, debug_bb);
            break;
          case SolidShape::FLOWER:
            Scan::Flower::draw<W, H, false>(generic_pipeline, canvas, basis,
                                            RADIUS, SIDES, shader, PHASE,
                                            debug_bb);
            break;
          case SolidShape::STAR:
            Scan::Star::draw<W, H, false>(generic_pipeline, canvas, basis,
                                          RADIUS, SIDES, shader, PHASE,
                                          debug_bb);
            break;
          }
        }
        generic_fx.advance_display();
        // A pair of all-black frames would satisfy the parity check below.
        const size_t generic_lit = count_lit_region<W, H>(generic_fx);
        HS_EXPECT_GT(generic_lit, (size_t)0);
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x)
            generic_pixels.push_back(generic_fx.get_pixel(x, y));
      }

      {
        hs_test::StubEffect solid_fx(W, H);
        Pipeline<W, H> solid_pipeline;
        {
          Canvas canvas(solid_fx);
          Basis basis = make_basis(
              make_rotation(Vector(0.3f, 0.7f, 0.2f).normalized(), 0.4f),
              Y_AXIS);
          switch (shape) {
          case SolidShape::PLANAR_POLYGON:
            Scan::PlanarPolygon::draw_solid<W, H>(solid_pipeline, canvas, basis,
                                                  RADIUS, SIDES, color, PHASE,
                                                  debug_bb);
            break;
          case SolidShape::SPHERICAL_POLYGON:
            Scan::SphericalPolygon::draw_solid<W, H>(solid_pipeline, canvas,
                                                     basis, RADIUS, SIDES,
                                                     color, PHASE, debug_bb);
            break;
          case SolidShape::FLOWER:
            Scan::Flower::draw_solid<W, H>(solid_pipeline, canvas, basis,
                                           RADIUS, SIDES, color, PHASE,
                                           debug_bb);
            break;
          case SolidShape::STAR:
            Scan::Star::draw_solid<W, H>(solid_pipeline, canvas, basis, RADIUS,
                                         SIDES, color, PHASE, debug_bb);
            break;
          }
        }
        solid_fx.advance_display();
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            const Pixel &generic = generic_pixels[y * W + x];
            const Pixel &solid = solid_fx.get_pixel(x, y);
            HS_EXPECT_EQ(generic.r, solid.r);
            HS_EXPECT_EQ(generic.g, solid.g);
            HS_EXPECT_EQ(generic.b, solid.b);
          }
      }
    }
  }
}

/**
 * @brief Bounds spherical sine-distance framebuffer error at device resolution.
 * @details The two paths' distance gap is fast_acos' ~1.3e-4 rad wherever the
 *   circumscribed-disc clamp wins (see the sine-domain note on
 *   SphericalPolygon::sine_distance); the coverage ramp scales that by the
 *   quintic kernel's slope over 2*pixel_width, so a channel may swing a few
 *   hundred ppm of full scale on the handful of pixels straddling a vertex.
 */
inline void test_spherical_sine_distance_framebuffer_error() {
  constexpr int W = 288;
  constexpr int H = 144;
  const Color4 color(Pixel(61000, 43000, 17000), 0.73f);
  Basis basis = make_basis(
      make_rotation(Vector(-0.2f, 0.9f, 0.4f).normalized(), 0.63f), Y_AXIS);
  struct Case {
    float radius;
    int sides;
    float phase;
  };
  const Case cases[] = {{0.24f, 3, -2.2f},
                        {0.74f, 5, 0.37f},
                        {0.98f, 12, 4.8f},
                        {1.38f, 7, -0.61f}};

  auto render = [&]<bool SineDistance>(const Case &c) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      Scan::SphericalPolygon::draw_solid<W, H, SineDistance>(
          pipeline, canvas, basis, c.radius, c.sides, color, c.phase);
    }
    fx.advance_display();
    // A pair of all-black frames would satisfy the error bounds below.
    const size_t lit = count_lit_region<W, H>(fx);
    HS_EXPECT_GT(lit, (size_t)0);
    std::vector<Pixel> pixels;
    pixels.reserve(W * H);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        pixels.push_back(fx.get_pixel(x, y));
    return pixels;
  };

  size_t different_pixels = 0;
  int max_channel_error = 0;
  for (const Case &c : cases) {
    const std::vector<Pixel> exact = render.template operator()<false>(c);
    const std::vector<Pixel> sine = render.template operator()<true>(c);
    for (size_t i = 0; i < exact.size(); ++i) {
      int dr = std::abs(static_cast<int>(exact[i].r) - sine[i].r);
      int dg = std::abs(static_cast<int>(exact[i].g) - sine[i].g);
      int db = std::abs(static_cast<int>(exact[i].b) - sine[i].b);
      int pixel_error = std::max(dr, std::max(dg, db));
      if (pixel_error != 0)
        ++different_pixels;
      max_channel_error = std::max(max_channel_error, pixel_error);
    }
  }

  std::printf("spherical sine framebuffer samples=%d different=%zu max=%d\n",
              W * H * static_cast<int>(std::size(cases)), different_pixels,
              max_channel_error);
  HS_EXPECT_LE(different_pixels, static_cast<size_t>(512));
  HS_EXPECT_LE(max_channel_error, 128);
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
  hs_test::StubEffect fx(W, H);
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
  /** Pixel coordinates handed to the integer plot(). */
  std::vector<std::pair<int, int>> plotted;
  std::vector<float> alpha; /**< Composited alpha per plotted pixel. */
  void plot(Canvas &, int x, int y, const Pixel &, float, float a) {
    plotted.push_back({x, y});
    alpha.push_back(a);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
  void plot(Canvas &, const Vector &, const Pixel &, float, float) {}
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

  hs_test::StubEffect fx(W, H);
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
 * plots at the same pixel (background first, foreground second). This pins that
 * branch: the duplicate-position signature must appear, the background must be
 * laid down opaque, and the foreground must be a partial (0<α<1) blend over it,
 * not a fade to black.
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

  hs_test::StubEffect fx(W, H);
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

  // The occluder branch is the only path that plots the same pixel twice
  // back-to-back (background then foreground); every other path plots each pixel
  // once.
  int occ_pairs = 0;
  float bg_alpha = 0.0f, fg_alpha = 1.0f;
  for (size_t i = 1; i < sink.plotted.size(); ++i) {
    if (sink.plotted[i] == sink.plotted[i - 1]) {
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
 * @brief Verifies Volume::draw's pole-LOD block skip drops no covered column.
 * @details The ray marcher settles a whole block from one clear trace, the only
 * pole-LOD site that skips without splatting, so an over-wide skip erases
 * columns outright. The cap sits on the canvas pole, where the walk offers the
 * widest blocks, and the sphere is offset off the view axis so its silhouette
 * runs across each decimated row rather than around it — on-axis every column
 * of a row carries the same clearance and the comparison is vacuous. The
 * framebuffer must match the undecimated walk at every aggressiveness.
 */
inline void test_volume_pole_lod_matches_undecimated() {
  constexpr int W = 96, H = 64;
  const ScopedPoleLod scoped_lod(0.0f);
  // Radial view of a cap centred on the canvas pole.
  const Vector center(0.0f, 1.0f, 0.0f);
  const Vector view_dir(0.0f, -1.0f, 0.0f);
  const float bounds_radius = 0.45f;

  // Counts the traces the walk actually runs, so the block skip is shown to
  // fire rather than assumed.
  struct CountingSphere {
    Vector center;
    float radius;
    long long *calls;
    float distance(const Vector &p) const {
      ++*calls;
      return (p - center).length() - radius;
    }
  };

  long long traced = 0;
  auto render = [&](float lod) {
    pole_lod_aggressiveness = lod;
    traced = 0;
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Pipeline<W, H> pipe;
      CountingSphere sphere{Vector(0.18f, 0.0f, 0.0f), 0.18f, &traced};
      Scan::TransformedVolume vol(sphere, center, Quaternion());
      Scan::Volume::draw<W, H>(
          pipe, c, center, bounds_radius, view_dir, vol,
          [](const Vector &, Fragment &f) {
            f.color = Color4(Pixel(60000, 45000, 30000), 1.0f);
          },
          /*max_steps=*/24, /*aa_width=*/0.01f);
    }
    fx.advance_display();
    std::vector<Pixel> out(W * H);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out[y * W + x] = fx.get_pixel(x, y);
    return out;
  };

  const std::vector<Pixel> plain = render(0.0f);
  const long long plain_traced = traced;
  bool blocks_settled = false;
  for (float lod : {1.0f, 4.0f}) {
    HS_CONTEXT("lod", static_cast<int>(lod * 100.0f));
    // The rows the cap reaches must actually be decimated, or this is vacuous.
    pole_lod_aggressiveness = lod;
    HS_EXPECT_GT(Scan::pole_lod_run(TrigLUT<W, H>::sin_phi[2]), 1);

    const std::vector<Pixel> decimated = render(lod);
    // A wide run demands slack the clearance cannot clear, so only the modest
    // aggressiveness actually settles blocks here.
    blocks_settled = blocks_settled || traced < plain_traced;
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
  }
  // Blocks were settled from one trace, so the parity above is not vacuous.
  HS_EXPECT_TRUE(blocks_settled);
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
    float closest_d = Scan::Volume::trace_closest(
        shape, ro, vd, 0.5f, max_steps, aa_width, closest_local);
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

/**
 * @brief Verifies overrelaxed sphere tracing never steps over a surface.
 * @details Sweeps rays across a twisted torus whose Lipschitz-divided distance
 * badly underestimates the true one — the case overrelaxation exploits — and
 * compares each trace against a dense fixed-step scan of the same ray. A ray
 * whose true closest approach lies well inside the AA band must be reported as a
 * hit: an unsafe step scheme strides over the ring and reports a miss instead.
 * No ray may report a hit the dense scan cannot corroborate.
 */
inline void test_volume_trace_closest_overrelax_never_skips_surface() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{{0.45f, 0.14f},
                                                        {2, 0.35f, 0.45f}};
  const float bounds_radius = 0.72f;
  const float aa_width = 0.07f;
  const Vector vd(0.0f, 0.0f, -1.0f);

  int covered = 0, hits = 0;
  for (int iy = -14; iy <= 14; ++iy) {
    for (int ix = -14; ix <= 14; ++ix) {
      Vector ro(ix * 0.045f, iy * 0.045f, bounds_radius);

      Vector closest_local;
      float closest_d = Scan::Volume::trace_closest(
          torus, ro, vd, bounds_radius, 18, aa_width, closest_local);

      // Dense reference: true closest approach along the same ray segment.
      float ref_min = FLT_MAX;
      for (int s = 0; s <= 4000; ++s) {
        Vector p(ro.x, ro.y, ro.z - s * (2.0f * bounds_radius / 4000.0f));
        float d = torus.distance(p);
        if (d < ref_min)
          ref_min = d;
      }

      if (closest_d < aa_width) {
        ++hits;
        // A reported hit must correspond to a real approach on the ray.
        HS_EXPECT_LT(ref_min, aa_width);
      }
      // Rays grazing the band edge may legitimately land either side of it;
      // anything this far inside is unambiguous coverage the march owes.
      if (ref_min < aa_width * 0.9f) {
        ++covered;
        HS_EXPECT_LT(closest_d, aa_width);
      }
    }
  }
  HS_EXPECT_GT(covered, 400);
  HS_EXPECT_GT(hits, covered);
}

// ============================================================================
// Scan::Circle::draw and Scan::Point::draw — the radius-0 ring regime
//
// Both wrappers build an SDF::Ring with radius 0, where target_angle is 0, the
// linearized-distance shortcut is disabled, and the axis' horizontal projection
// can fall under MIN_HORIZONTAL_PROJ so every row takes the full-row scan. The
// shape that regime draws is a spherical cap of angular radius `thickness`
// centred on the basis axis, with quintic coverage from 1 at the centre to 0 at
// the rim — the analytic oracle these cases compare against.
// ============================================================================

/**
 * @brief Analytic coverage of a radius-0 ring at one direction.
 * @param v Unit direction of the pixel centre.
 * @param axis Cap centre (the ring's basis axis).
 * @param thickness Cap angular radius in radians.
 * @return Stroke coverage in [0, 1]; 0 at or beyond the rim.
 * @details Mirrors SDF::Ring::distance at radius 0 op for op — the same
 * fast_acos and the same quintic — so the comparison is exact, not approximate.
 */
inline float cap_coverage(const Vector &v, const Vector &axis,
                          float thickness) {
  const float d = dot(v, axis);
  const float angle = fast_acos(hs::clamp(d, -1.0f, 1.0f));
  const float sd = angle - thickness;
  if (sd >= 0.0f || thickness <= 0.0f)
    return 0.0f;
  return quintic_kernel(-sd / thickness);
}

/**
 * @brief Verifies Point::draw paints exactly the analytic spherical cap.
 * @details Sweeps a pole-centred axis (where the ring's horizontal projection
 * collapses and every row full-row scans) and an equatorial axis (where the row
 * interval math runs), and requires both to reproduce cap_coverage per pixel:
 * every covered direction lit, every uncovered one black, and the plotted
 * channel proportional to the coverage.
 */
inline void test_point_draws_the_analytic_cap() {
  constexpr int W = 96, H = 64;
  constexpr uint16_t LEVEL = 60000;
  constexpr float THICKNESS = 0.35f;
  TrigLUT<W, H>::init();

  const Vector axes[] = {Y_AXIS, Vector(1.0f, 0.0f, 0.0f),
                         Vector(0.3f, -0.6f, 0.74f).normalized()};
  for (const Vector &axis : axes) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Scan::Point::draw<W, H>(
          pipe, c, axis, THICKNESS, [](const Vector &, Fragment &f) {
            f.color = Color4(Pixel(LEVEL, LEVEL, LEVEL), 1.0f);
          });
    }
    fx.advance_display();

    // The cap centre must be the basis axis the wrapper built, not the raw
    // vector: make_basis can reorient, so read it back rather than assuming.
    const Basis basis = make_basis(Quaternion(), axis);
    size_t lit = 0, covered = 0;
    float worst_value_error = 0.0f;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const Vector v = pixel_to_vector<W, H>(x, y);
        const float coverage = cap_coverage(v, basis.v, THICKNESS);
        const bool is_lit = !is_black(fx.get_pixel(x, y));
        if (is_lit)
          ++lit;
        if (coverage > 0.0f)
          ++covered;
        // Coverage 0 is a hard cut in the shape, so an outside pixel is black.
        if (coverage == 0.0f)
          HS_EXPECT_FALSE(is_lit);
        // Inside the cut, the plotted channel is the level scaled by coverage.
        if (coverage > 0.02f) {
          HS_EXPECT_TRUE(is_lit);
          const float want = static_cast<float>(LEVEL) * coverage;
          worst_value_error = std::max(
              worst_value_error,
              std::fabs(static_cast<float>(fx.get_pixel(x, y).r) - want));
        }
      }
    HS_EXPECT_GT(covered, (size_t)0);
    HS_EXPECT_GT(lit, (size_t)0);
    // Blend rounding only; the coverage itself is reproduced exactly.
    HS_EXPECT_LT(worst_value_error, 2.0f);
  }
}

/**
 * @brief Verifies a pole-centred cap really exercises the full-row-scan regime.
 * @details Point::draw at +Y builds a ring whose axis has no horizontal
 * projection, so get_horizontal_intervals refuses every row it is asked for and
 * the rasterizer falls back to scanning whole rows. Without this the cap case
 * above would be indistinguishable from the ordinary interval path.
 */
inline void test_pole_centred_cap_takes_the_full_row_scan() {
  constexpr int W = 96, H = 64;
  TrigLUT<W, H>::init();
  const Basis pole = make_basis(Quaternion(), Y_AXIS);
  const SDF::Ring degenerate(pole, 0.0f, 0.35f);

  const auto rows = degenerate.get_vertical_bounds<H>();
  HS_EXPECT_EQ(rows.y_min, 0);
  HS_EXPECT_GT(rows.y_max, 0);
  int full_scan_rows = 0;
  for (int y = rows.y_min; y <= rows.y_max; ++y)
    if (degenerate.needs_full_row_scan(TrigLUT<W, H>::sin_phi[y]))
      ++full_scan_rows;
  HS_EXPECT_EQ(full_scan_rows, rows.y_max - rows.y_min + 1);

  // An equatorial cap of the same size answers with intervals on most rows, so
  // the two axes above cover both sides of the branch.
  const Basis equator = make_basis(Quaternion(), Vector(1.0f, 0.0f, 0.0f));
  const SDF::Ring ordinary(equator, 0.0f, 0.35f);
  const auto eq_rows = ordinary.get_vertical_bounds<H>();
  int interval_rows = 0;
  for (int y = eq_rows.y_min; y <= eq_rows.y_max; ++y)
    if (!ordinary.needs_full_row_scan(TrigLUT<W, H>::sin_phi[y]))
      ++interval_rows;
  HS_EXPECT_GT(interval_rows, 0);
}

/**
 * @brief Verifies the Circle and Point wrappers are their documented rings.
 * @details Circle is a radius-0 ring whose stroke half-width is radius * pi/2;
 * Point is a radius-0 ring of the given thickness. Both must rasterize
 * bit-identically to the ring they claim to be, so the wrapper cannot drift
 * from the shape its callers (ShapeShifter, Comets) were tuned against.
 */
inline void test_circle_and_point_match_their_rings() {
  constexpr int W = 96, H = 64;
  Pipeline<W, H> pipe;
  auto shader = [](const Vector &p, Fragment &f) {
    // Position-dependent so a mismatched basis or phase shows up as color, not
    // only as coverage.
    f.color =
        Color4(Pixel(static_cast<uint16_t>(30000.0f + 20000.0f * p.x),
                     static_cast<uint16_t>(30000.0f + 20000.0f * p.y), 50000),
               1.0f);
  };

  const Basis basis = make_basis(Quaternion(), Vector(0.2f, 0.8f, -0.5f));
  std::vector<Pixel> circle_frame, ring_frame;
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Circle::draw<W, H>(pipe, c, basis, 0.4f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, circle_frame);
  }
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Ring::draw<W, H>(pipe, c, basis, 0.0f, 0.4f * (PI_F / 2.0f),
                             shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, ring_frame);
  }
  size_t circle_lit = 0, circle_diff = 0;
  for (size_t i = 0; i < circle_frame.size(); ++i) {
    if (!is_black(circle_frame[i]))
      ++circle_lit;
    if (!(circle_frame[i] == ring_frame[i]))
      ++circle_diff;
  }
  HS_EXPECT_GT(circle_lit, (size_t)0);
  HS_EXPECT_EQ(circle_diff, (size_t)0);

  const Vector center(-0.4f, 0.5f, 0.766f);
  std::vector<Pixel> point_frame, point_ring_frame;
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Point::draw<W, H>(pipe, c, center, 0.3f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, point_frame);
  }
  {
    hs_test::StubEffect fx(W, H);
    const Basis point_basis = make_basis(Quaternion(), center);
    {
      Canvas c(fx);
      Scan::Ring::draw<W, H>(pipe, c, point_basis, 0.0f, 0.3f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, point_ring_frame);
  }
  size_t point_lit = 0, point_diff = 0;
  for (size_t i = 0; i < point_frame.size(); ++i) {
    if (!is_black(point_frame[i]))
      ++point_lit;
    if (!(point_frame[i] == point_ring_frame[i]))
      ++point_diff;
  }
  HS_EXPECT_GT(point_lit, (size_t)0);
  HS_EXPECT_EQ(point_diff, (size_t)0);
}

/**
 * @brief Verifies a Circle's painted extent tracks its radius argument.
 * @details radius maps to a cap of angular half-width radius * pi/2, so the
 * lit area must grow with radius and the farthest lit direction must sit at
 * that angle. A wrapper that dropped the pi/2 or passed radius as the ring's
 * radius rather than its thickness fails both.
 */
inline void test_circle_extent_follows_its_radius() {
  constexpr int W = 96, H = 64;
  const Basis basis = make_basis(Quaternion(), Vector(0.0f, 0.0f, 1.0f));
  size_t previous_lit = 0;
  for (float radius : {0.15f, 0.3f, 0.6f}) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Scan::Circle::draw<W, H>(
          pipe, c, basis, radius, [](const Vector &, Fragment &f) {
            f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
          });
    }
    fx.advance_display();

    const float rim = radius * (PI_F / 2.0f);
    size_t lit = 0;
    float widest = 0.0f;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        if (!is_black(fx.get_pixel(x, y))) {
          ++lit;
          const Vector v = pixel_to_vector<W, H>(x, y);
          widest = std::max(widest,
                            fast_acos(hs::clamp(dot(v, basis.v), -1.0f, 1.0f)));
        }
    HS_EXPECT_GT(lit, previous_lit);
    // No lit pixel past the rim, and the cap is sampled close enough to it
    // that the widest lit direction is within a row of the rim.
    HS_EXPECT_LT(widest, rim);
    HS_EXPECT_GT(widest, rim - 2.0f * (PI_F / (H - 1)));
    previous_lit = lit;
  }
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
  test_shader_clip_arc_matches_predicate();
  test_ring_rasterize_produces_bounded_output();
  test_ring_rasterize_lit_pixels_on_band();
  test_ring_rasterize_lights_expected_row();
  test_stroke_aa_is_monotone_ramp();
  test_ring_rasterize_empty_clip_draws_nothing();
  test_distorted_ring_flat_matches_zero_knot_raster();
  test_ring_group_matches_sequential();
  test_distorted_ring_stack_matches_sequential();
  test_face_rasterize_matches_scan_region();
  test_scan_shader_v2_contract();
  test_scan_region_seam_no_double_plot();
  test_scan_region_fractional_boundary_no_double_plot();
  test_scan_region_clip_arc_matches_predicate();
  test_pole_lod_runs_are_canvas_anchored();
  test_pole_lod_shading_matches_undecimated();
  test_plot_line_over_pole_reaches_row0();
  test_csg_stroke_aa_uses_winning_child_thickness();

  test_star_pixel_placement();
  test_planar_polygon_pixel_placement();
  test_flower_pixel_placement();
  test_solid_color_path_matches_generic();
  test_spherical_sine_distance_framebuffer_error();
  test_overlapping_strokes_composite_blend();

  test_point_draws_the_analytic_cap();
  test_pole_centred_cap_takes_the_full_row_scan();
  test_circle_and_point_match_their_rings();
  test_circle_extent_follows_its_radius();

  test_transformed_volume_world_local_roundtrip();
  test_volume_raymarch_silhouette_and_registers();
  test_volume_draw_occluded_edge_blends_over_background();
  test_volume_pole_lod_matches_undecimated();
  test_volume_trace_closest_stops_at_first_graze();
  test_volume_probe_occluder_reports_background_graze_point();
  test_volume_trace_closest_overrelax_never_skips_surface();

  return fixture.result();
}

} // namespace scan_tests
} // namespace hs_test
