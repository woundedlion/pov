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

#include <atomic>
#include <chrono>
#include <thread>

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
  void add_bool(const char *n, bool *p, bool d) {
    *p = d; // registerParam(bool) captures *ptr as the default; set it first
    registerParam(n, p);
  }
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

// The relaxed-atomic double-buffer hand-off (canvas.h ctor reasoning): the
// writer (main loop, cur_) must never claim the buffer the display side (ISR,
// prev_) is reading, and a queued-but-not-displayed frame must not disturb the
// live frame. True ISR concurrency isn't deterministically unit-testable, but
// the single-threaded state machine that the relaxed atomics implement is —
// drive many cycles and assert the non-aliasing / no-torn-read guarantee
// observably.
inline void test_double_buffer_handoff_no_aliasing() {
  TestEffect fx(8, 4);
  const int N = 8 * 4;

  const Pixel *seen[2] = {nullptr, nullptr};
  int distinct = 0;
  auto record_ptr = [&](const Pixel *p) {
    for (int i = 0; i < distinct; ++i)
      if (seen[i] == p) return;
    if (distinct < 2) seen[distinct++] = p;
  };

  const Pixel colors[6] = {Pixel(10, 0, 0), Pixel(0, 20, 0), Pixel(0, 0, 30),
                           Pixel(40, 0, 0), Pixel(0, 50, 0), Pixel(0, 0, 60)};

  // Frame 0 establishes a displayed baseline.
  { Canvas c(fx); for (int i = 0; i < N; ++i) c(i) = colors[0]; }
  fx.advance_display();
  record_ptr(fx.display_buffer());
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), 10, 0, 0));

  for (int f = 1; f < 6; ++f) {
    const Pixel *display_before = fx.display_buffer();

    // Draw the next frame but do NOT advance the display.
    { Canvas c(fx); for (int i = 0; i < N; ++i) c(i) = colors[f]; }

    // Frame queued: the writer claimed the OTHER buffer, so the displayed
    // pointer is unchanged and still shows the PREVIOUS frame — the newly
    // written frame is invisible until the display side picks it up.
    HS_EXPECT_FALSE(fx.buffer_free());
    HS_EXPECT_TRUE(fx.display_buffer() == display_before);
    HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), colors[f - 1].r, colors[f - 1].g,
                          colors[f - 1].b));

    fx.advance_display();

    // Display side advanced to the queued frame; the writer's next buffer is a
    // DIFFERENT physical buffer than the one just displayed.
    HS_EXPECT_TRUE(fx.buffer_free());
    HS_EXPECT_TRUE(
        pix_eq(fx.get_pixel(0, 0), colors[f].r, colors[f].g, colors[f].b));
    HS_EXPECT_TRUE(fx.display_buffer() != display_before);
    record_ptr(fx.display_buffer());
  }

  // Only two physical buffers were ever displayed across all cycles.
  HS_EXPECT_EQ(distinct, 2);
}

// Directly exercises the Canvas ctor's buffer_free() spin-wait — the one
// synchronization gate the rest of the suite deliberately steps around by
// calling advance_display() before every ctor. On real hardware the display
// ISR consumes frames asynchronously; here a helper thread plays that ISR.
//
// With a frame queued-but-not-displayed the buffer is busy, so the ctor MUST
// block. The release is deterministic, not timing-based: the ctor can only
// return once prev_ == next_, and the sole writer of prev_ is advance_display()
// — which only the helper runs. So the helper's "ctor has not returned yet"
// assertion holds by construction (it checks before advancing), and an inverted
// gate (spin while buffer_free()) would let the ctor return early and fail it.
// The ctor's 2 s watchdog bounds the spin, so a logic break traps loudly here
// rather than hanging the suite.
inline void test_ctor_spin_waits_for_buffer_free() {
  hs::clear_mock_time(); // use the real wall clock so the spin/watchdog are live
  TestEffect fx(8, 8);

  // Queue frame 0 WITHOUT advancing the display: the buffer is now busy.
  {
    Canvas c(fx);
    c(0, 0) = Pixel(1, 2, 3);
  }
  HS_EXPECT_FALSE(fx.buffer_free());

  std::atomic<bool> ctor_returned{false};
  std::atomic<bool> released{false};

  std::thread display_isr([&] {
    // Give the main thread time to actually enter the spin, then confirm the
    // ctor is still blocked (it cannot have returned while the buffer is busy)
    // before consuming the frame.
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    HS_EXPECT_FALSE(ctor_returned.load(std::memory_order_acquire));
    released.store(true, std::memory_order_relaxed);
    fx.advance_display(); // the "ISR" frees the buffer
  });

  // Blocks in the ctor spin-wait until the helper advances the display.
  { Canvas c(fx); }
  ctor_returned.store(true, std::memory_order_release);

  display_isr.join();
  HS_EXPECT_TRUE(released.load(std::memory_order_relaxed));
  // The helper's advance_display() promoted frame 0 to the displayed buffer, so
  // its pixel is now live. (buffer_free() is back to false here — the second
  // Canvas's dtor already queued its own frame — so we assert visibility, not
  // the gate.)
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), 1, 2, 3));
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
  HS_EXPECT_NEAR(fl->get(), 1.0f, 1e-6f); // captured *ptr (true) → reads as 1.0
  HS_EXPECT_TRUE(fx.flag);                // registerParam(bool) leaves *ptr as-is

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
  test_double_buffer_handoff_no_aliasing();
  test_ctor_spin_waits_for_buffer_free();
  test_canvas_2d_and_1d_access_and_prev();
  test_register_float_and_bool_params();
  test_update_parameter_by_name();
  test_paramlist_fills_to_capacity();
  test_clip_setters();

  return hs_test::end_module(scope);
}

} // namespace canvas_tests
} // namespace hs_test

