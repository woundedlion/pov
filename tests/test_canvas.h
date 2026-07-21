/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/render/canvas.h — the Effect double-buffer state machine, the
 * parameter system (register_param / updateParameter / ParamList), the clip
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

#include "core/render/canvas.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace canvas_tests {

/**
 * @brief Minimal concrete Effect for exercising the base-class machinery.
 * @details Exposes the protected hooks (register_param / persist_pixels) needed
 * by the tests.
 */
struct TestEffect : public Effect {
  float speed = 1.5f; /**< Sample float param backing store. */
  bool flag = false;  /**< Sample bool param backing store. */

  /**
   * @brief Constructs the test effect with the given canvas dimensions.
   * @param W Canvas width in pixels.
   * @param H Canvas height in pixels.
   */
  TestEffect(int W, int H, EffectConfig cfg = {}) : Effect(W, H, cfg) {}
  /**
   * @brief Per-frame draw hook; intentionally a no-op for these tests.
   */
  void draw_frame() override {}

  /**
   * @brief Toggles frame-to-frame pixel persistence.
   * @param b True to copy the previous frame into each new buffer.
   */
  void set_persist(bool b) { persist_pixels = b; }
  /**
   * @brief Registers a float parameter with the base Effect.
   * @param n Parameter name.
   * @param p Pointer to the float backing store (current value is the default).
   * @param mn Minimum allowed value.
   * @param mx Maximum allowed value.
   */
  void add_float(const char *n, float *p, float mn, float mx) {
    register_param(n, p, mn, mx);
  }
  /**
   * @brief Registers a bool parameter, seeding its default first.
   * @param n Parameter name.
   * @param p Pointer to the bool backing store.
   * @param d Default value written into *p before registration.
   */
  void add_bool(const char *n, bool *p, bool d) {
    *p = d; // register_param(bool) captures *ptr as the default; set it first
    register_param(n, p);
  }
  /**
   * @brief Marks a registered parameter as readonly (engine-written).
   * @param n Parameter name to mark.
   */
  void mark_readonly(const char *n) { markReadonly(n); }
  /**
   * @brief Registers an enumerated parameter with the base Effect.
   * @param n Parameter name.
   * @param p Pointer to the float backing store holding the option index.
   * @param options Option label array (static storage).
   * @param count Number of labels.
   */
  void add_enum(const char *n, float *p, const char *const *options,
                int count) {
    register_param(n, p, options, count);
  }
};

/**
 * @brief Tests whether a pixel has exactly the given channel values.
 * @param p Pixel to inspect.
 * @param r Expected red channel value.
 * @param g Expected green channel value.
 * @param b Expected blue channel value.
 * @return True if all three channels match exactly.
 */
inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

// ============================================================================
// Construction
// ============================================================================

/**
 * @brief Verifies a freshly constructed Effect reports its dimensions, starts
 * fully black with a free display buffer, and has its clip set to the whole
 * canvas.
 */
inline void test_construction_dims_and_clear() {
  TestEffect fx(96, 20);
  HS_EXPECT_EQ(fx.width(), 96);
  HS_EXPECT_EQ(fx.height(), 20);
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(10, 5), 0, 0, 0));
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(fx.clip().is_full());
  HS_EXPECT_EQ(fx.clip().y_end, 20);
  HS_EXPECT_EQ(fx.clip().x_end, 96);
}

/**
 * @brief Verifies the smallest and largest supported dimensions are accepted.
 */
inline void test_construction_dimension_boundaries() {
  {
    TestEffect fx(1, 1);
    HS_EXPECT_EQ(fx.width(), 1);
    HS_EXPECT_EQ(fx.height(), 1);
    HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), 0, 0, 0));
  }
  {
    TestEffect fx(MAX_W, MAX_H);
    HS_EXPECT_EQ(fx.width(), MAX_W);
    HS_EXPECT_EQ(fx.height(), MAX_H);
    HS_EXPECT_TRUE(
        pix_eq(fx.get_pixel(MAX_W - 1, MAX_H - 1), 0, 0, 0));
  }
}

// ============================================================================
// Double-buffer state machine
// ============================================================================

/**
 * @brief Verifies a drawn frame stays invisible until advance_display()
 * promotes it.
 * @details The display shows the old buffer while drawing and while the frame
 * is merely queued.
 */
inline void test_frame_visible_only_after_advance_display() {
  TestEffect fx(8, 8);

  {
    Canvas c(fx);
    c(3, 3) = Pixel(100, 200, 300);
    // Mid-frame: the display still shows the OLD (black) buffer.
    HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 0, 0, 0));
  } // ~Canvas queues the frame

  // Queued but not displayed: buffer busy, pixel still old.
  HS_EXPECT_FALSE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 0, 0, 0));

  fx.advance_display();
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(3, 3), 100, 200, 300));
}

/**
 * @brief Verifies that without persist, each frame draws into a freshly cleared
 * buffer, so pixels from the previous frame do not survive into the next.
 */
inline void test_consecutive_frames_alternate_buffers() {
  TestEffect fx(8, 8);

  { Canvas c(fx); c(1, 1) = Pixel(11, 0, 0); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 11, 0, 0));

  // Non-persist clears the new buffer, so the first frame's pixel is gone.
  { Canvas c(fx); c(2, 2) = Pixel(0, 22, 0); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 0, 22, 0));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 0, 0, 0));
}

/** @brief Clip clearing touches the display rectangle but not its margin. */
inline void test_clip_clear_exact_rectangle() {
  TestEffect fx(8, 4);

  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(11, 22, 33); }
  fx.advance_display();
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(44, 55, 66); }
  fx.advance_display();

  fx.set_clip(1, 3, 2, 6);
  fx.set_margin(2);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int y = 0; y < 4; ++y) {
      for (int x = 0; x < 8; ++x) {
        const bool inside = y >= 1 && y < 3 && x >= 2 && x < 6;
        HS_EXPECT_TRUE(inside ? pix_eq(c(x, y), 0, 0, 0)
                              : pix_eq(c(x, y), 11, 22, 33));
      }
    }
  }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 1), 0, 0, 0));
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(1, 1), 11, 22, 33));
}

/** @brief Alternating POV clips clear the recycled buffer before publication. */
inline void test_clip_clear_alternates_clips_and_buffers() {
  TestEffect fx(8, 4);

  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(10, 0, 0); }
  fx.advance_display();
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(0, 20, 0); }
  fx.advance_display();

  fx.set_clip(0, 4, 0, 4);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int y = 0; y < 4; ++y)
      for (int x = 0; x < 8; ++x)
        HS_EXPECT_TRUE(x < 4 ? pix_eq(c(x, y), 0, 0, 0)
                             : pix_eq(c(x, y), 10, 0, 0));
    for (int y = 0; y < 4; ++y)
      for (int x = 0; x < 4; ++x)
        c(x, y) = Pixel(0, 0, 30);
  }
  fx.advance_display();

  fx.set_clip(0, 4, 4, 8);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int y = 0; y < 4; ++y)
      for (int x = 0; x < 8; ++x)
        HS_EXPECT_TRUE(x >= 4 ? pix_eq(c(x, y), 0, 0, 0)
                              : pix_eq(c(x, y), 0, 20, 0));
  }
  fx.advance_display();

  fx.set_clip(0, 4, 0, 4);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int y = 0; y < 4; ++y)
      for (int x = 0; x < 8; ++x)
        HS_EXPECT_TRUE(x < 4 ? pix_eq(c(x, y), 0, 0, 0)
                             : pix_eq(c(x, y), 10, 0, 0));
  }
}

/** @brief A full display clip is equivalent to the generic full clear. */
inline void test_clip_clear_full_clip_matches_full_clear() {
  TestEffect fx(8, 4);
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(10, 20, 30); }
  fx.advance_display();
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(40, 50, 60); }
  fx.advance_display();

  Canvas c(fx, Canvas::ClearDisplayClipTag{});
  for (int i = 0; i < 32; ++i)
    HS_EXPECT_TRUE(pix_eq(c(i), 0, 0, 0));
}

/** @brief A full-width band clears contiguously at the correct row offset. */
inline void test_clip_clear_full_width_band() {
  TestEffect fx(8, 4);
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(10, 20, 30); }
  fx.advance_display();
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(40, 50, 60); }
  fx.advance_display();

  fx.set_clip(1, 3, 0, 8);
  Canvas c(fx, Canvas::ClearDisplayClipTag{});
  for (int y = 0; y < 4; ++y)
    for (int x = 0; x < 8; ++x)
      HS_EXPECT_TRUE(y >= 1 && y < 3 ? pix_eq(c(x, y), 0, 0, 0)
                                     : pix_eq(c(x, y), 10, 20, 30));
}

/** @brief Effect construction still fully zeroes both shared buffers. */
inline void test_clip_clear_constructor_zeroes_both_buffers() {
  {
    TestEffect dirty(8, 4);
    { Canvas c(dirty); for (int i = 0; i < 32; ++i) c(i) = Pixel(1, 2, 3); }
    dirty.advance_display();
    { Canvas c(dirty); for (int i = 0; i < 32; ++i) c(i) = Pixel(4, 5, 6); }
    dirty.advance_display();
  }

  TestEffect fx(8, 4);
  fx.set_clip(1, 3, 2, 6);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int i = 0; i < 32; ++i)
      HS_EXPECT_TRUE(pix_eq(c(i), 0, 0, 0));
  }
  fx.advance_display();
  fx.set_clip(0, 1, 0, 1);
  {
    Canvas c(fx, Canvas::ClearDisplayClipTag{});
    for (int i = 0; i < 32; ++i)
      HS_EXPECT_TRUE(pix_eq(c(i), 0, 0, 0));
  }
}

/** @brief Adding the tagged path leaves persistent-frame copying unchanged. */
inline void test_clip_clear_does_not_change_persistence() {
  TestEffect fx(8, 4, {.persist = true});
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(7, 8, 9); }
  fx.advance_display();

  fx.set_clip(1, 3, 2, 6);
  Canvas c(fx);
  for (int i = 0; i < 32; ++i)
    HS_EXPECT_TRUE(pix_eq(c(i), 7, 8, 9));
}

/** @brief Adding the tagged path leaves generic full-frame clearing unchanged. */
inline void test_clip_clear_does_not_change_full_frame_clear() {
  TestEffect fx(8, 4, {.full_frame = true});
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(1, 2, 3); }
  fx.advance_display();
  { Canvas c(fx); for (int i = 0; i < 32; ++i) c(i) = Pixel(4, 5, 6); }
  fx.advance_display();

  fx.set_clip(1, 3, 2, 6);
  Canvas c(fx);
  for (int i = 0; i < 32; ++i)
    HS_EXPECT_TRUE(pix_eq(c(i), 0, 0, 0));
}

/**
 * @brief Verifies that with persist_pixels set, each new frame inherits the
 * previous frame's contents, so undrawn pixels retain their prior color.
 */
inline void test_persist_pixels_copies_previous_frame() {
  TestEffect fx(8, 8);
  fx.set_persist(true);

  { Canvas c(fx); c(2, 2) = Pixel(10, 20, 30); }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 10, 20, 30));

  // With persist, the undrawn next frame carries the previous content over.
  { Canvas c(fx); /* draw nothing */ }
  fx.advance_display();
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(2, 2), 10, 20, 30));
}

/**
 * @brief Verifies the relaxed-atomic double-buffer hand-off via its
 * single-threaded state machine.
 * @details The writer (main loop, cur_) must never claim the buffer the display
 * side (ISR, prev_) is reading, and a queued-but-not-displayed frame must not
 * disturb the live frame. True ISR concurrency isn't deterministically
 * unit-testable, but the single-threaded state machine that the relaxed atomics
 * implement is — drive many cycles and assert the non-aliasing / no-torn-read
 * guarantee observably.
 */
inline void test_double_buffer_handoff_no_aliasing() {
  TestEffect fx(8, 4);
  const int N = 8 * 4;

  // Capacity 4 (> the 2 expected) so a third/fourth buffer would be recorded.
  const Pixel *seen[4] = {nullptr, nullptr, nullptr, nullptr};
  int distinct = 0;
  auto record_ptr = [&](const Pixel *p) {
    for (int i = 0; i < distinct; ++i)
      if (seen[i] == p) return;
    HS_EXPECT_TRUE(distinct < 4);
    if (distinct < 4) seen[distinct++] = p;
  };

  const Pixel colors[6] = {Pixel(10, 0, 0), Pixel(0, 20, 0), Pixel(0, 0, 30),
                           Pixel(40, 0, 0), Pixel(0, 50, 0), Pixel(0, 0, 60)};

  { Canvas c(fx); for (int i = 0; i < N; ++i) c(i) = colors[0]; }
  fx.advance_display();
  record_ptr(fx.display_buffer());
  HS_EXPECT_TRUE(fx.buffer_free());
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), 10, 0, 0));

  for (int f = 1; f < 6; ++f) {
    const Pixel *display_before = fx.display_buffer();

    // Draw the next frame but do NOT advance the display.
    { Canvas c(fx); for (int i = 0; i < N; ++i) c(i) = colors[f]; }

    // Queued frame is invisible: the displayed pointer is unchanged and still
    // shows the previous frame.
    HS_EXPECT_FALSE(fx.buffer_free());
    HS_EXPECT_TRUE(fx.display_buffer() == display_before);
    HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), colors[f - 1].r, colors[f - 1].g,
                          colors[f - 1].b));

    fx.advance_display();

    HS_EXPECT_TRUE(fx.buffer_free());
    HS_EXPECT_TRUE(
        pix_eq(fx.get_pixel(0, 0), colors[f].r, colors[f].g, colors[f].b));
    HS_EXPECT_TRUE(fx.display_buffer() != display_before);
    record_ptr(fx.display_buffer());
  }

  HS_EXPECT_EQ(distinct, 2);
}

/**
 * @brief Hammers the double-buffer hand-off under real producer/consumer
 * contention, asserting no torn read and no out-of-order frame.
 * @details test_double_buffer_handoff_no_aliasing drives the same state machine
 * single-threaded; this runs the two roles on separate threads so the relaxed
 * atomics are exercised under genuine concurrency. A PRODUCER thread plays the
 * main loop: it constructs a Canvas (whose ctor blocks on buffer_free(), so the
 * consumer rate-limits it), fills every pixel with a sentinel encoding the frame
 * index, and lets ~Canvas queue it. A CONSUMER thread plays the display ISR:
 * when a frame is queued (!buffer_free()) it advance_display()s and reads the
 * whole displayed buffer, asserting all pixels carry one sentinel (a torn read
 * mixing two frames trips torn_read) and that the sentinel strictly advances (a
 * doubly-displayed or regressed frame trips out_of_order). The buffer_free()
 * gate forces the producer one frame ahead at most, so every frame is displayed
 * exactly once — distinct_displayed must reach the full frame count, proving the
 * hand-off really happened rather than frames coalescing. The harness counters
 * are single-threaded, so both threads record into atomics and the main thread
 * runs all HS_EXPECT_* after join(). The Canvas ctor's 2 s watchdog (live under
 * the real wall clock) bounds the producer, and the consumer terminates once it
 * has promoted every frame, so a logic break traps loudly instead of hanging.
 */
inline void test_double_buffer_handoff_concurrent() {
  hs::clear_mock_time(); // real wall clock keeps the ctor's spin watchdog live
  TestEffect fx(8, 4);
  const int N = 8 * 4;
  const int FRAMES = 3000;
  // r holds the frame index (FRAMES < 65536, so r alone is a unique sentinel).
  auto sentinel = [](int f) { return Pixel((uint16_t)f, 0, 0); };

  std::atomic<bool> producer_done{false};
  std::atomic<bool> torn_read{false};
  std::atomic<bool> out_of_order{false};
  std::atomic<int> distinct_displayed{0};

  std::thread producer([&] {
    for (int f = 0; f < FRAMES; ++f) {
      Canvas c(fx); // blocks until the consumer frees the buffer
      for (int i = 0; i < N; ++i) c(i) = sentinel(f);
    } // ~Canvas queues the frame
    producer_done.store(true, std::memory_order_release);
  });

  std::thread consumer([&] {
    int last_sentinel = -1;
    int promoted = 0;
    while (promoted < FRAMES) {
      if (fx.buffer_free()) {
        // Nothing queued; if the producer has finished there is no more work.
        if (producer_done.load(std::memory_order_acquire) && fx.buffer_free())
          break;
        std::this_thread::yield();
        continue;
      }
      fx.advance_display();
      ++promoted;
      const Pixel *buf = fx.display_buffer();
      int s = buf[0].r;
      bool uniform = true;
      for (int i = 1; i < N; ++i)
        if (buf[i].r != s) { uniform = false; break; }
      if (!uniform) torn_read.store(true, std::memory_order_relaxed);
      if (s <= last_sentinel) out_of_order.store(true, std::memory_order_relaxed);
      last_sentinel = s;
      distinct_displayed.store(promoted, std::memory_order_relaxed);
    }
  });

  producer.join();
  consumer.join();

  HS_EXPECT_FALSE(torn_read.load(std::memory_order_relaxed));
  HS_EXPECT_FALSE(out_of_order.load(std::memory_order_relaxed));
  // The gate admits no coalescing, so every queued frame is displayed once.
  HS_EXPECT_EQ(distinct_displayed.load(std::memory_order_relaxed), FRAMES);
}

/**
 * @brief Verifies the Canvas ctor's buffer_free() spin-wait blocks until the
 * display side frees the buffer.
 * @details Directly exercises the one synchronization gate the rest of the
 * suite deliberately steps around by calling advance_display() before every
 * ctor. On real hardware the display ISR consumes frames asynchronously; here a
 * helper thread plays that ISR. With a frame queued-but-not-displayed the
 * buffer is busy, so the ctor MUST block. The release is deterministic, not
 * timing-based: the ctor can only return once prev_ == next_, and the sole
 * writer of prev_ is advance_display() — which only the helper runs. So the
 * helper's "ctor has not returned yet" assertion holds by construction (it
 * checks before advancing), and an inverted gate (spin while buffer_free())
 * would let the ctor return early and fail it. The ctor's 2 s watchdog bounds
 * the spin, so a logic break traps loudly here rather than hanging the suite.
 */
inline void test_ctor_spin_waits_for_buffer_free() {
  hs::clear_mock_time(); // use the real wall clock so the spin/watchdog are live
  TestEffect fx(8, 8);

  // Queue frame 0 WITHOUT advancing: the buffer is now busy.
  {
    Canvas c(fx);
    c(0, 0) = Pixel(1, 2, 3);
  }
  HS_EXPECT_FALSE(fx.buffer_free());

  std::atomic<bool> ctor_returned{false};
  std::atomic<bool> released{false};
  // The harness counters are single-threaded, so the helper records here and
  // the assertions run on the main thread after join().
  std::atomic<bool> ctor_blocked_when_checked{false};

  const unsigned long spins_before = Canvas::buffer_free_spin_count();
  std::thread display_isr([&] {
    // Wait until the ctor is observably spinning on buffer_free() before freeing
    // the buffer, so the release is gated on real progress rather than a fixed
    // sleep a loaded CI host could outrun. Observing a spin iteration proves the
    // ctor is in the wait loop and has not yet returned. The ctor_returned exit
    // keeps a broken/inverted gate (ctor returns without spinning) from hanging
    // here: it falls through and the assertions below fail loudly instead.
    while (Canvas::buffer_free_spin_count() == spins_before &&
           !ctor_returned.load(std::memory_order_acquire))
      std::this_thread::yield();
    ctor_blocked_when_checked.store(
        !ctor_returned.load(std::memory_order_acquire),
        std::memory_order_relaxed);
    released.store(true, std::memory_order_relaxed);
    fx.advance_display(); // the "ISR" frees the buffer
  });

  // Blocks in the ctor spin-wait until the helper advances the display.
  { Canvas c(fx); }
  ctor_returned.store(true, std::memory_order_release);

  display_isr.join();
  HS_EXPECT_TRUE(ctor_blocked_when_checked.load(std::memory_order_relaxed));
  HS_EXPECT_TRUE(released.load(std::memory_order_relaxed));
  // Frame 0 is now the displayed buffer, so its pixel is live.
  HS_EXPECT_TRUE(pix_eq(fx.get_pixel(0, 0), 1, 2, 3));
}

// ============================================================================
// Canvas access
// ============================================================================
//
// OOB pixel access is not a death-harness case: Canvas index accessors guard
// with a debug-only `assert` (SIGABRT), not HS_CHECK's __builtin_trap (SIGILL)
// that the death harness keys on. On device (NDEBUG) the access is unchecked.

/**
 * @brief Verifies Canvas exposes 2D (x,y) and 1D (row-major index) writes that
 * target the same buffer, and that prev() reads the previously displayed frame.
 */
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

/**
 * @brief Verifies register_param captures each param's current pointee as its
 * default and exposes type, value, and min/max through getParameters().
 * @details Also checks that find() of an unregistered name returns null.
 */
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
  HS_EXPECT_TRUE(fx.flag);                // register_param(bool) leaves *ptr as-is

  HS_EXPECT_TRUE(params.find("Missing") == nullptr);
}

/**
 * @brief Verifies updateParameter writes a float param by name, thresholds at
 * 0.5 for bool params, and returns false (leaving values untouched) for unknown
 * names or non-finite inputs.
 */
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

  // Unknown name is a no-op returning false.
  HS_EXPECT_FALSE(fx.updateParameter("Nope", 99.0f));
  HS_EXPECT_NEAR(fx.speed, 7.25f, 1e-6f);

  // Non-finite values are rejected, returning false.
  HS_EXPECT_FALSE(fx.updateParameter("Speed", std::numeric_limits<float>::quiet_NaN()));
  HS_EXPECT_NEAR(fx.speed, 7.25f, 1e-6f);
}

/**
 * @brief Verifies the untrusted JS boundary cannot write readonly
 * (engine-written telemetry) params, even though the GUI already disables their
 * editing.
 */
inline void test_update_parameter_rejects_readonly() {
  TestEffect fx(4, 4);
  float telemetry = 1.0f;
  fx.add_float("Telemetry", &telemetry, 0.0f, 10.0f);
  fx.mark_readonly("Telemetry");

  HS_EXPECT_FALSE(fx.updateParameter("Telemetry", 5.0f));
  HS_EXPECT_NEAR(telemetry, 1.0f, 1e-6f);

  // An editable param in the same effect still writes.
  fx.add_float("Speed", &fx.speed, 0.0f, 10.0f);
  HS_EXPECT_TRUE(fx.updateParameter("Speed", 4.0f));
  HS_EXPECT_NEAR(fx.speed, 4.0f, 1e-6f);
}

/**
 * @brief Verifies enum registration derives its range from the option count and
 * exposes the labels, and that updateParameter snaps writes to a valid index.
 * @details Fractional writes round to the nearest option; out-of-range writes
 *          clamp to the ends, so the target always holds a real option index.
 */
inline void test_register_and_update_enum_param() {
  TestEffect fx(4, 4);
  static const char *const MODES[] = {"None", "Warp", "Sparkle"};
  float mode = 1.0f;
  fx.add_enum("Mode", &mode, MODES, 3);

  const auto *def = fx.getParameters().find("Mode");
  HS_EXPECT_TRUE(def != nullptr);
  HS_EXPECT_TRUE(def->is_enum());
  HS_EXPECT_FALSE(def->is_bool());
  HS_EXPECT_EQ(def->option_count, 3);
  HS_EXPECT_TRUE(def->options == MODES);
  HS_EXPECT_NEAR(def->min, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(def->max, 2.0f, 1e-6f);
  HS_EXPECT_NEAR(def->get(), 1.0f, 1e-6f); // captured current value as default

  // A fractional write snaps to the nearest option index.
  HS_EXPECT_TRUE(fx.updateParameter("Mode", 1.7f));
  HS_EXPECT_NEAR(mode, 2.0f, 1e-6f);

  // Out-of-range writes clamp to the option range.
  HS_EXPECT_TRUE(fx.updateParameter("Mode", 9.0f));
  HS_EXPECT_NEAR(mode, 2.0f, 1e-6f);
  HS_EXPECT_TRUE(fx.updateParameter("Mode", -3.0f));
  HS_EXPECT_NEAR(mode, 0.0f, 1e-6f);

  // Plain params stay non-enum.
  fx.add_float("Speed", &fx.speed, 0.0f, 10.0f);
  HS_EXPECT_FALSE(fx.getParameters().find("Speed")->is_enum());
}

/**
 * @brief Verifies ParamList holds its full capacity of registered params.
 */
inline void test_paramlist_fills_to_capacity() {
  TestEffect fx(4, 4);
  static float vals[32];
  // Unique names "pNN": register_param traps on duplicates, and ParamDef stores
  // the char* by pointer, so the names need static storage.
  static char names[32][4];
  for (int i = 0; i < 32; ++i) {
    vals[i] = static_cast<float>(i);
    names[i][0] = 'p';
    names[i][1] = static_cast<char>('0' + i / 10);
    names[i][2] = static_cast<char>('0' + i % 10);
    names[i][3] = '\0';
    fx.add_float(names[i], &vals[i], 0.0f, 100.0f);
  }
  HS_EXPECT_EQ(fx.getParameters().size(), (size_t)32);
}

// ============================================================================
// Clip setters
// ============================================================================

/**
 * @brief Verifies the clip setters write valid bounds through the canvas edge.
 * @details set_clip sets all four bounds (and clears is_full), set_clip_x
 * touches only x bounds, and set_margin sets margin.
 */
inline void test_clip_setters() {
  TestEffect fx(96, 20);
  fx.set_clip(2, 10, 5, 96);
  HS_EXPECT_EQ(fx.clip().y_start, 2);
  HS_EXPECT_EQ(fx.clip().y_end, 10);
  HS_EXPECT_EQ(fx.clip().x_start, 5);
  HS_EXPECT_EQ(fx.clip().x_end, 96);
  HS_EXPECT_FALSE(fx.clip().is_full());

  fx.set_clip_x(7, 96);
  HS_EXPECT_EQ(fx.clip().x_start, 7);
  HS_EXPECT_EQ(fx.clip().x_end, 96);
  HS_EXPECT_EQ(fx.clip().y_start, 2);

  fx.set_margin(3);
  HS_EXPECT_EQ(fx.clip().margin, 3);
}

struct StubPipe {
  int integer_hits = 0;
  int float_hits = 0;
  int integer_x = 0;
  int integer_y = 0;
  float float_x = 0.0f;
  float float_y = 0.0f;

  void plot(Canvas &, int x, int y, const Pixel &, float, float) {
    ++integer_hits;
    integer_x = x;
    integer_y = y;
  }
  void plot(Canvas &, float x, float y, const Pixel &, float, float) {
    ++float_hits;
    float_x = x;
    float_y = y;
  }
  void plot(Canvas &, const Vector &, const Pixel &, float, float) {}
};

/**
 * @brief Verifies integer and floating-point coordinates retain their overload.
 */
inline void test_pipeline_ref_routes_screen_coordinate_overloads() {
  StubPipe pipe;
  PipelineRef pipeline(pipe);
  TestEffect fx(8, 4);
  Canvas canvas(fx);

  pipeline.plot(canvas, -3, 7, Pixel(0, 0, 0), 0.0f, 0.0f);
  HS_EXPECT_EQ(pipe.integer_hits, 1);
  HS_EXPECT_EQ(pipe.float_hits, 0);
  HS_EXPECT_EQ(pipe.integer_x, -3);
  HS_EXPECT_EQ(pipe.integer_y, 7);

  pipeline.plot(canvas, 1.25f, -2.5f, Pixel(0, 0, 0), 0.0f, 0.0f);
  HS_EXPECT_EQ(pipe.integer_hits, 1);
  HS_EXPECT_EQ(pipe.float_hits, 1);
  HS_EXPECT_NEAR(pipe.float_x, 1.25f, 0.0f);
  HS_EXPECT_NEAR(pipe.float_y, -2.5f, 0.0f);
}

/**
 * @brief Verifies a PipelineRef copy remains bound to its original target.
 */
inline void test_pipeline_ref_copy_is_independent_of_source() {
  StubPipe s1, s2;
  PipelineRef a(s1);
  PipelineRef b(a);
  a = PipelineRef(s2);
  {
    TestEffect fx(8, 4);
    Canvas c(fx);
    b.plot(c, 2, 3, Pixel(0, 0, 0), 0.0f, 0.0f);
  }
  HS_EXPECT_EQ(s1.integer_hits, 1);
  HS_EXPECT_EQ(s2.integer_hits, 0);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every canvas test in order.
 * @return The module's failure count reported by end_module().
 */
inline int run_canvas_tests() {
  hs_test::ModuleFixture fixture("canvas");

  test_construction_dims_and_clear();
  test_construction_dimension_boundaries();
  test_frame_visible_only_after_advance_display();
  test_consecutive_frames_alternate_buffers();
  test_clip_clear_exact_rectangle();
  test_clip_clear_alternates_clips_and_buffers();
  test_clip_clear_full_clip_matches_full_clear();
  test_clip_clear_full_width_band();
  test_clip_clear_constructor_zeroes_both_buffers();
  test_clip_clear_does_not_change_persistence();
  test_clip_clear_does_not_change_full_frame_clear();
  test_persist_pixels_copies_previous_frame();
  test_double_buffer_handoff_no_aliasing();
  test_double_buffer_handoff_concurrent();
  test_ctor_spin_waits_for_buffer_free();
  test_canvas_2d_and_1d_access_and_prev();
  test_register_float_and_bool_params();
  test_update_parameter_by_name();
  test_update_parameter_rejects_readonly();
  test_register_and_update_enum_param();
  test_paramlist_fills_to_capacity();
  test_clip_setters();
  test_pipeline_ref_routes_screen_coordinate_overloads();
  test_pipeline_ref_copy_is_independent_of_source();

  return fixture.result();
}

} // namespace canvas_tests
} // namespace hs_test

