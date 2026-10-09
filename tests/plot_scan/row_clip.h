/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// ClipRegion::could_intersect_y  (clip.h — pure clip culling)
// ============================================================================

/**
 * @brief Verifies ClipRegion::could_intersect_y culls by screen-row range:
 *        segments overlapping the clip band pass, fully-above/below segments are
 *        rejected, and the test is order-independent in its two y arguments.
 */
inline void test_clip_could_intersect_y() {
  ClipRegion cr;
  cr.y_start = 40;
  cr.y_end = 80;
  cr.x_start = 0;
  cr.x_end = MAX_W;
  cr.margin = 0;
  cr.h = MAX_H;
  cr.w = MAX_W;

  HS_EXPECT_FALSE(cr.is_full());

  HS_EXPECT_EQ(cr.render_y_start(), 40);
  HS_EXPECT_EQ(cr.render_y_end(), 80);

  HS_EXPECT_TRUE(cr.could_intersect_y(50.0f, 60.0f));
  HS_EXPECT_TRUE(cr.could_intersect_y(30.0f, 45.0f));

  HS_EXPECT_FALSE(cr.could_intersect_y(0.0f, 39.0f));
  HS_EXPECT_FALSE(cr.could_intersect_y(80.0f, 120.0f));

  // Order independence: swapped arguments give the same answer.
  HS_EXPECT_TRUE(cr.could_intersect_y(60.0f, 50.0f));
  HS_EXPECT_FALSE(cr.could_intersect_y(120.0f, 80.0f));
}

/**
 * @brief Asserts contains_x() and !XClip::clipped() agree on every column.
 */
inline void expect_xclip_parity(const ClipRegion &cr) {
  const ClipRegion::XClip xc = cr.x_clip();
  for (int x = 0; x < cr.w; ++x) {
    HS_EXPECT_EQ(cr.contains_x(x), !xc.clipped(x));
  }
}

/**
 * @brief Exercises the cylindrical x-clip predicates (render_x_*, contains_x,
 *        x_clip/XClip) across every documented band topology.
 */
inline void test_clip_x_band_topologies() {
  constexpr int W = 96;

  auto make = [](int x0, int x1, int margin) {
    ClipRegion cr;
    cr.y_start = 0;
    cr.y_end = MAX_H;
    cr.x_start = x0;
    cr.x_end = x1;
    cr.margin = margin;
    cr.w = W;
    cr.h = MAX_H;
    return cr;
  };

  // Non-wrapping sub-band: [20,60) expanded by 2 -> render band [18,62).
  {
    ClipRegion cr = make(20, 60, 2);
    HS_EXPECT_EQ(cr.render_x_start(), 18);
    HS_EXPECT_EQ(cr.render_x_end(), 62);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_TRUE(xc.active);
    HS_EXPECT_FALSE(xc.wrap);
    HS_EXPECT_EQ(xc.rs, 18);
    HS_EXPECT_EQ(xc.re, 62);
    HS_EXPECT_FALSE(cr.contains_x(17));
    HS_EXPECT_TRUE(cr.contains_x(18));
    HS_EXPECT_TRUE(cr.contains_x(61));
    HS_EXPECT_FALSE(cr.contains_x(62)); // re exclusive
    expect_xclip_parity(cr);
  }

  // Partial seam-crossing band: [2,90) expanded by 3 -> [95, w) U [0, 93);
  // 88 + 2*3 = 94 < w leaves the gap {93,94} clipped.
  {
    ClipRegion cr = make(2, 90, 3);
    HS_EXPECT_EQ(cr.render_x_start(), 95);
    HS_EXPECT_EQ(cr.render_x_end(), 93);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_TRUE(xc.active);
    HS_EXPECT_TRUE(xc.wrap);
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(92));  // last column before the gap
    HS_EXPECT_FALSE(cr.contains_x(93)); // re exclusive
    HS_EXPECT_FALSE(cr.contains_x(94)); // interior of the clipped gap
    HS_EXPECT_TRUE(cr.contains_x(95));  // rs inclusive
    expect_xclip_parity(cr);
  }

  // Exact wrap to full width: 40 + 2*28 = w, so rs == re == 78 yet the band is
  // full, not empty.
  {
    ClipRegion cr = make(10, 50, 28);
    HS_EXPECT_EQ(cr.render_x_start(), 78);
    HS_EXPECT_EQ(cr.render_x_end(), 78);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active);
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(78));
    HS_EXPECT_TRUE(cr.contains_x(95));
    expect_xclip_parity(cr);
  }

  // Explicit full-width band.
  {
    ClipRegion cr = make(0, W, 0);
    HS_EXPECT_TRUE(cr.is_full());
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active);
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(W - 1));
    expect_xclip_parity(cr);
  }

  // Over-wrap to full coverage: 88 + 2*8 = 104 >= w, though rs != re.
  {
    ClipRegion cr = make(2, 90, 8);
    HS_EXPECT_EQ(cr.render_x_start(), 90);
    HS_EXPECT_EQ(cr.render_x_end(), 2);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active);
    HS_EXPECT_TRUE(cr.contains_x(2));  // display left edge
    HS_EXPECT_TRUE(cr.contains_x(50)); // display interior
    HS_EXPECT_TRUE(cr.contains_x(89)); // display right edge
    HS_EXPECT_TRUE(cr.contains_x(90)); // gap, but reached from both margins
    expect_xclip_parity(cr);
  }

  // Empty band: zero width, no margin.
  {
    ClipRegion cr = make(30, 30, 0);
    HS_EXPECT_EQ(cr.render_x_start(), 30);
    HS_EXPECT_EQ(cr.render_x_end(), 30);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_TRUE(xc.active);
    HS_EXPECT_FALSE(xc.wrap);
    HS_EXPECT_FALSE(cr.contains_x(0));
    HS_EXPECT_FALSE(cr.contains_x(29));
    HS_EXPECT_FALSE(cr.contains_x(30));
    HS_EXPECT_FALSE(cr.contains_x(W - 1));
    expect_xclip_parity(cr);
  }
}

/**
 * @brief Pins render_x_start/render_x_end (and contains_x) to modular
 *        arithmetic over the whole domain Canvas::set_clip/set_margin allow.
 */
inline void test_clip_x_wrap_matches_modulo() {
  auto ref_start = [](int x_start, int margin, int w) {
    return (x_start - margin + w) % w;
  };
  auto ref_end = [](int x_end, int margin, int w) {
    return (x_end + margin) % w;
  };

  // Edge accessors: full (x, margin) sweep.
  for (int w : {1, 2, 3, 7, 96, 288}) {
    for (int x = 0; x <= w; ++x) {
      for (int margin = 0; margin < w; ++margin) {
        ClipRegion cr;
        cr.x_start = x;
        cr.x_end = x;
        cr.margin = margin;
        cr.w = w;
        HS_EXPECT_EQ(cr.render_x_start(), ref_start(x, margin, w));
        HS_EXPECT_EQ(cr.render_x_end(), ref_end(x, margin, w));
      }
    }
  }

  // Whole predicate: every band and margin at a small width, every column.
  constexpr int MAX_SWEEP_W = 13;
  for (int w : {1, 3, MAX_SWEEP_W}) {
    for (int x0 = 0; x0 <= w; ++x0) {
      for (int x1 = x0; x1 <= w; ++x1) {
        for (int margin = 0; margin < w; ++margin) {
          ClipRegion cr;
          cr.x_start = x0;
          cr.x_end = x1;
          cr.margin = margin;
          cr.w = w;
          cr.h = MAX_H;
          cr.y_end = MAX_H;
          bool covered[MAX_SWEEP_W] = {};
          for (int k = x0 - margin; k < x1 + margin; ++k)
            covered[((k % w) + w) % w] = true;
          for (int x = 0; x < w; ++x) {
            HS_EXPECT_EQ(cr.contains_x(x), covered[x]);
          }
          expect_xclip_parity(cr);
        }
      }
    }
  }
}
