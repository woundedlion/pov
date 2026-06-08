/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/canvas.h — the Effect double-buffer state machine, the
 * parameter system (registerParam / updateParameter / ParamList), the clip
 * setters, and the Canvas scoped drawing context.
 *
 * Frame protocol note: a Canvas spins in its ctor while !buffer_free(), so every
 * test that draws a frame MUST advance_display() before constructing the next
 * Canvas (otherwise the host build would spin forever on a still-pending frame).
 */
#pragma once

#include "core/canvas.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace canvas_tests {

// Minimal concrete Effect for exercising the base-class machinery. Exposes the
// protected hooks (registerParam / persist_pixels) needed by the tests.
struct TestEffect : public Effect {
  float speed = 1.5f;
  bool flag = false;

  TestEffect(int W, int H) : Effect(W, H) {}
  void draw_frame() override {}
  bool show_bg() const override { return false; }

  void set_persist(bool b) { persist_pixels = b; }
  void add_float(const char *n, float *p, float mn, float mx) {
    registerParam(n, p, mn, mx);
  }
  void add_bool(const char *n, bool *p, bool d) { registerParam(n, p, d); }
};

inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

// ============================================================================
// Construction
// ============================================================================

inline void test_construction_dims_and_clear() {
  TestEffect fx(96, 20);
  HS_EXPECT_EQ(fx.width(), 96);
  HS_EXPECT_EQ(fx.height(), 20);
  // Both buffers start black; display is ready (no frame queued).
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(10, 5), 0, 0, 0));
  HS_EXPECT_TRUE(fx.buffer_free());
  // Clip initialized to the full canvas.
  HS_EXPECT_TRUE(fx.clip.is_full());
  HS_EXPECT_EQ(fx.clip.y_end, 20);
  HS_EXPECT_EQ(fx.clip.x_end, 96);
}

// ============================================================================
// Double-buffer state machine
// ============================================================================

inline void test_frame_visible_only_after_advance_display() {
  TestEffect fx(8, 8);

  {
    Canvas c(fx);
    c(3, 3) = Pixel(100, 200, 300);
    // Mid-frame: the display still shows the OLD (black) buffer.
    HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 0, 0, 0));
  } // ~Canvas queues the frame

  // Queued but not yet displayed: buffer is "busy" and pixel still old.
  HS_EXPECT_FALSE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 0, 0, 0));

  fx.advance_display();
  // Now the drawn frame is live, and the buffer is free for the next frame.
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 100, 200, 300));
}

inline void test_consecutive_frames_alternate_buffers() {
  TestEffect fx(8, 8);

  { Canvas c(fx); c(1, 1) = Pixel(11, 0, 0); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 11, 0, 0));

  // Second frame draws a different pixel; non-persist clears the new buffer, so
  // the first pixel is gone.
  { Canvas c(fx); c(2, 2) = Pixel(0, 22, 0); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 0, 22, 0));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 0, 0, 0)); // cleared
}

inline void test_persist_pixels_copies_previous_frame() {
  TestEffect fx(8, 8);
  fx.set_persist(true);

  { Canvas c(fx); c(2, 2) = Pixel(10, 20, 30); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 10, 20, 30));

  // Next frame draws nothing; with persist the previous content carries over.
  { Canvas c(fx); /* draw nothing */ }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 10, 20, 30));
}

// ============================================================================
// Canvas access
// ============================================================================

inline void test_canvas_2d_and_1d_access_and_prev() {
  TestEffect fx(8, 8);
  {
    Canvas c(fx);
    HS_EXPECT_EQ(c.width(), 8);
    HS_EXPECT_EQ(c.height(), 8);
    c(1, 1) = Pixel(5, 6, 7);
    c(2 * 8 + 3) = Pixel(8, 9, 10); // 1D index → (x=3, y=2)
    // prev() reads the previously displayed (black) buffer.
    HS_EXPECT_TRUE(pix_eq(c.prev(1, 1), 0, 0, 0));
  }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 5, 6, 7));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 2), 8, 9, 10));
}

// ============================================================================
// Parameter system
// ============================================================================

inline void test_register_float_and_bool_params() {
  TestEffect fx(4, 4);
  fx.add_float("Speed", &fx.speed, 0.0f, 10.0f);
  fx.add_bool("Flag", &fx.flag, true);

  const auto &params = fx.getParameters();
  HS_EXPECT_EQ(params.size(), (size_t)2);

  const auto *sp = params.find("Speed");
  HS_EXPECT_TRUE(sp != nullptr);
  HS_EXPECT_FALSE(sp->is_bool());
  HS_EXPECT_NEAR(sp->get(), 1.5f, 1e-6f); // captured current value as default
  HS_EXPECT_NEAR(sp->min, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(sp->max, 10.0f, 1e-6f);

  const auto *fl = params.find("Flag");
  HS_EXPECT_TRUE(fl != nullptr);
  HS_EXPECT_TRUE(fl->is_bool());
  HS_EXPECT_NEAR(fl->get(), 1.0f, 1e-6f); // default true → bool reads as 1.0
  HS_EXPECT_TRUE(fx.flag);                // registerParam(bool) writes *ptr

  HS_EXPECT_TRUE(params.find("Missing") == nullptr);
}

inline void test_update_parameter_by_name() {
  TestEffect fx(4, 4);
  fx.add_float("Speed", &fx.speed, 0.0f, 10.0f);
  fx.add_bool("Flag", &fx.flag, false);

  HS_EXPECT_TRUE(fx.updateParameter("Speed", 7.25f));
  HS_EXPECT_NEAR(fx.speed, 7.25f, 1e-6f);

  // Bool target uses a 0.5 threshold.
  HS_EXPECT_TRUE(fx.updateParameter("Flag", 0.2f));
  HS_EXPECT_FALSE(fx.flag);
  HS_EXPECT_TRUE(fx.updateParameter("Flag", 0.8f));
  HS_EXPECT_TRUE(fx.flag);

  // Unknown name is a safe no-op and now reports false to the caller.
  HS_EXPECT_FALSE(fx.updateParameter("Nope", 99.0f));
  HS_EXPECT_NEAR(fx.speed, 7.25f, 1e-6f);

  // Non-finite values are rejected and also report false.
  HS_EXPECT_FALSE(fx.updateParameter("Speed", std::numeric_limits<float>::quiet_NaN()));
  HS_EXPECT_NEAR(fx.speed, 7.25f, 1e-6f);
}

inline void test_paramlist_fills_to_capacity() {
  TestEffect fx(4, 4);
  static float vals[32];
  // Registering exactly ParamList's capacity (std::array<ParamDef, 32>) is
  // valid. Registering a 33rd now TRAPS (an effect-authoring bug, fail-fast),
  // so the overflow path can't be exercised here without death-test infra.
  for (int i = 0; i < 32; ++i) {
    vals[i] = static_cast<float>(i);
    fx.add_float("p", &vals[i], 0.0f, 100.0f);
  }
  HS_EXPECT_EQ(fx.getParameters().size(), (size_t)32);
}

// ============================================================================
// Clip setters
// ============================================================================

inline void test_clip_setters() {
  TestEffect fx(96, 20);
  fx.set_clip(2, 10, 5, 40);
  HS_EXPECT_EQ(fx.clip.y_start, 2);
  HS_EXPECT_EQ(fx.clip.y_end, 10);
  HS_EXPECT_EQ(fx.clip.x_start, 5);
  HS_EXPECT_EQ(fx.clip.x_end, 40);
  HS_EXPECT_FALSE(fx.clip.is_full());

  fx.set_clip_x(7, 50);
  HS_EXPECT_EQ(fx.clip.x_start, 7);
  HS_EXPECT_EQ(fx.clip.x_end, 50);
  HS_EXPECT_EQ(fx.clip.y_start, 2); // unchanged

  fx.set_margin(3);
  HS_EXPECT_EQ(fx.clip.margin, 3);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_canvas_tests() {
  auto scope = hs_test::begin_module("canvas");

  test_construction_dims_and_clear();
  test_frame_visible_only_after_advance_display();
  test_consecutive_frames_alternate_buffers();
  test_persist_pixels_copies_previous_frame();
  test_canvas_2d_and_1d_access_and_prev();
  test_register_float_and_bool_params();
  test_update_parameter_by_name();
  test_paramlist_fills_to_capacity();
  test_clip_setters();

  return hs_test::end_module(scope);
}

} // namespace canvas_tests
} // namespace hs_test

