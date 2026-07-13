/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the segmented-POV index math (hardware/pov_segment_map.h,
 * the pure arithmetic the Arduino-only pov_segmented.h driver derives its ISR
 * mapping from). This is the one place index arithmetic reaches the physical
 * LEDs, so an off-by-one silently mis-paints the sphere. Covers per-segment
 * derivation (arm side and north/south band direction), active-low ID decode,
 * the arm-B half-image x offset,
 * and that each segment's (i -> x_col, y) writes tile the canvas exactly once.
 */
#pragma once

#include "hardware/pov_handoff.h"
#include "hardware/pov_segment_map.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <vector>

namespace hs_test {
namespace pov_segmented_tests {

using pov::EffectHandoff;
using pov::decode_segment_id;
using pov::segment_clip;
using pov::segment_id_strap_count;
using pov::segment_map;
using pov::SegmentClip;
using pov::SegmentMap;
using pov::segment_x_col;
using pov::segment_y;

// Compile-time proof the mapping is genuinely constexpr (the driver relies on
// it folding at boot, off the ISR hot path).
static_assert(segment_map(1, 288, 4).y_base == 143);
static_assert(segment_map(1, 288, 4).y_step == -1);
static_assert(segment_map(1, 288, 8).y_base == 36);
static_assert(segment_map(2, 288, 8).y_base == 143);
static_assert(segment_map(3, 288, 8).y_base == 107);
static_assert(segment_id_strap_count(8) == 3);
static_assert(decode_segment_id(0b010, 8) == 5);
static_assert(segment_x_col(true, 0, 288) == 144);
static_assert(segment_clip(segment_map(0, 288, 4), true, 288, 4, 288).x0 == 0);
static_assert(segment_clip(segment_map(0, 288, 4), true, 288, 4, 288).x1 == 144);
static_assert(segment_clip(segment_map(2, 288, 4), true, 288, 4, 288).x0 == 144);
static_assert(segment_clip(segment_map(1, 288, 4), true, 288, 4, 288).y0 == 72);

/**
 * @brief Assert the N segments tile both arm columns exactly once at column @p x.
 * @param S Total LED count across all segments.
 * @param N Number of segments.
 * @param w Canvas width in columns (rotation resolution).
 * @param x Rotation column sampled by arm A (arm B samples (x+w/2)%w).
 * @details At a fixed rotation column, every segment's per-LED write must land
 * on a distinct canvas pixel and the N segments together must tile exactly the
 * two columns the arms sample (x for arm A, (x+w/2)%w for arm B), each row
 * covered once. The S LED writes must map onto 2*ROWS canvas pixels.
 */
inline void check_tiling(int S, int N, int w, int x) {
  const int PPS = S / N;
  const int ROWS = S / 2;
  const int col_a = x;
  const int col_b = (x + w / 2) % w;

  std::vector<int> cover(static_cast<size_t>(w) * ROWS, 0);
  int writes = 0;

  for (int seg = 0; seg < N; ++seg) {
    const SegmentMap m = segment_map(seg, S, N);
    const int x_col = segment_x_col(m.arm_b, x, w);
    HS_EXPECT_TRUE(x_col >= 0 && x_col < w);
    if (x_col < 0 || x_col >= w)
      continue; // keep an out-of-range column out of the cover grid

    int prev_y = -1;
    for (int i = 0; i < PPS; ++i) {
      const int y = segment_y(m, i);
      HS_EXPECT_TRUE(y >= 0 && y < ROWS);
      // Interior ordering: rows advance by the segment's unit y_step, so a
      // scrambled-but-bijective interior fails here, not just covered-once.
      if (i > 0)
        HS_EXPECT_EQ(y - prev_y, m.y_step);
      prev_y = y;
      cover[static_cast<size_t>(x_col) * ROWS + y]++;
      ++writes;
    }
  }

  HS_EXPECT_EQ(writes, S);

  int covered = 0;
  bool double_paint = false;
  bool stray_column = false;
  for (int c = 0; c < w; ++c) {
    for (int y = 0; y < ROWS; ++y) {
      const int hits = cover[static_cast<size_t>(c) * ROWS + y];
      if (hits > 1)
        double_paint = true;
      if (hits == 1)
        ++covered;
      if (hits != 0 && c != col_a && c != col_b)
        stray_column = true;
    }
  }
  HS_EXPECT_FALSE(double_paint);
  HS_EXPECT_FALSE(stray_column);
  HS_EXPECT_EQ(covered, 2 * ROWS);

  for (int y = 0; y < ROWS; ++y) {
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_a) * ROWS + y], 1);
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_b) * ROWS + y], 1);
  }
}

/**
 * @brief Verify per-segment SegmentMap fields for the 4-segment layout.
 * @details Checks arm side (A/B) and the y_base/y_step that distinguish top
 * strips (ascending from the pole) from bottom strips (descending, reversed),
 * plus that the strip endpoints map to the expected pole/junction rows.
 */
inline void test_segment_derivation() {
  // Phantasm config: N=4, S=288 -> ROWS=144, PPS=72.
  const int S = 288, N = 4;

  const SegmentMap s0 = segment_map(0, S, N); // arm A, top
  HS_EXPECT_FALSE(s0.arm_b);
  HS_EXPECT_EQ(s0.y_base, 0);
  HS_EXPECT_EQ(s0.y_step, 1);

  const SegmentMap s1 = segment_map(1, S, N); // arm A, bottom (reversed)
  HS_EXPECT_FALSE(s1.arm_b);
  HS_EXPECT_EQ(s1.y_base, 143);
  HS_EXPECT_EQ(s1.y_step, -1);

  const SegmentMap s2 = segment_map(2, S, N); // arm B, top
  HS_EXPECT_TRUE(s2.arm_b);
  HS_EXPECT_EQ(s2.y_base, 0);
  HS_EXPECT_EQ(s2.y_step, 1);

  const SegmentMap s3 = segment_map(3, S, N); // arm B, bottom (reversed)
  HS_EXPECT_TRUE(s3.arm_b);
  HS_EXPECT_EQ(s3.y_base, 143);
  HS_EXPECT_EQ(s3.y_step, -1);

  // Strip endpoints: the top strip runs N pole (y=0) -> junction (y=PPS-1);
  // the bottom strip runs S pole (y=ROWS-1) -> junction (y=ROWS-PPS), reversed.
  HS_EXPECT_EQ(segment_y(s0, 0), 0);
  HS_EXPECT_EQ(segment_y(s0, 71), 71);
  HS_EXPECT_EQ(segment_y(s1, 0), 143);
  HS_EXPECT_EQ(segment_y(s1, 71), 72);
}

/**
 * @brief Verify every row band and strip direction in the 8-segment layout.
 */
inline void test_eight_segment_derivation() {
  const int S = 288;
  const int N = 8;
  const int EXPECTED_BASE[N] = {0, 36, 143, 107, 0, 36, 143, 107};
  const int EXPECTED_STEP[N] = {1, 1, -1, -1, 1, 1, -1, -1};

  for (int seg = 0; seg < N; ++seg) {
    const SegmentMap m = segment_map(seg, S, N);
    HS_EXPECT_EQ(m.arm_b, seg >= N / 2);
    HS_EXPECT_EQ(m.y_base, EXPECTED_BASE[seg]);
    HS_EXPECT_EQ(m.y_step, EXPECTED_STEP[seg]);
  }

  HS_EXPECT_EQ(segment_y(segment_map(0, S, N), 35), 35);
  HS_EXPECT_EQ(segment_y(segment_map(1, S, N), 35), 71);
  HS_EXPECT_EQ(segment_y(segment_map(2, S, N), 35), 108);
  HS_EXPECT_EQ(segment_y(segment_map(3, S, N), 35), 72);
}

/**
 * @brief Verify strap counts and active-low decoding for every supported mode.
 */
inline void test_segment_id_straps() {
  HS_EXPECT_EQ(segment_id_strap_count(2), 1);
  HS_EXPECT_EQ(segment_id_strap_count(4), 2);
  HS_EXPECT_EQ(segment_id_strap_count(8), 3);

  for (int id = 0; id < 8; ++id) {
    const int RAW = (~id) & 7;
    HS_EXPECT_EQ(decode_segment_id(RAW, 8), id);
  }
  for (int id = 0; id < 4; ++id) {
    const int RAW = (~id) & 3;
    HS_EXPECT_EQ(decode_segment_id(RAW, 4), id);
  }
  for (int id = 0; id < 2; ++id) {
    const int RAW = (~id) & 1;
    HS_EXPECT_EQ(decode_segment_id(RAW, 2), id);
  }
}

/**
 * @brief Verify arm-B column offset and wrap behavior.
 * @details Arm-B columns sit half an image (w/2) ahead of the rotation column
 * and wrap modulo w; arm A samples the column unshifted.
 */
inline void test_arm_b_offset() {
  const int w = 288;
  HS_EXPECT_EQ(segment_x_col(false, 0, w), 0);
  HS_EXPECT_EQ(segment_x_col(false, 17, w), 17);
  HS_EXPECT_EQ(segment_x_col(true, 0, w), 144);
  HS_EXPECT_EQ(segment_x_col(true, 144, w), 0);   // wraps to the seam
  HS_EXPECT_EQ(segment_x_col(true, 287, w), 143); // (287 + 144) % 288
}

/**
 * @brief Verify per-window clip rectangles tile the canvas for one layout.
 * @param S Total LED count.
 * @param N Segment count.
 * @param w Canvas width.
 */
inline void check_segment_clips(int S, int N, int w) {
  const int ROWS = S / 2;
  const int PPS = S / N;

  for (bool arm_a_left : {true, false}) {
    std::vector<int> cover(static_cast<size_t>(w) * ROWS, 0);
    for (int seg = 0; seg < N; ++seg) {
      const SegmentMap m = segment_map(seg, S, N);
      const SegmentClip c = segment_clip(m, arm_a_left, S, N, w);

      // Half the width, one segment's row band.
      HS_EXPECT_EQ(c.x1 - c.x0, w / 2);
      HS_EXPECT_EQ(c.y1 - c.y0, PPS);
      const int EXPECTED_Y0 =
          m.y_step > 0 ? m.y_base : m.y_base - (PPS - 1);
      HS_EXPECT_EQ(c.y0, EXPECTED_Y0);
      // Arm B paints the opposite column half from arm A in the same window.
      HS_EXPECT_EQ(c.x0 == 0, m.arm_b ? !arm_a_left : arm_a_left);

      for (int y = c.y0; y < c.y1; ++y)
        for (int x = c.x0; x < c.x1; ++x)
          cover[static_cast<size_t>(x) * ROWS + y]++;
    }
    for (size_t i = 0; i < cover.size(); ++i)
      HS_EXPECT_EQ(cover[i], 1);
  }
}

/**
 * @brief Verify clips for every supported segment count.
 */
inline void test_segment_clip() {
  check_segment_clips(/*S=*/288, /*N=*/2, /*w=*/288);
  check_segment_clips(/*S=*/288, /*N=*/4, /*w=*/288);
  check_segment_clips(/*S=*/288, /*N=*/8, /*w=*/288);
}

// Stand-in for Effect: the handoff never dereferences the pointee, only tracks
// ownership by address, so an empty tag type exercises every code path.
struct FakeEffect {};

/**
 * @brief Teardown handshake: request, then ISR service acks and drops the live
 * pointer.
 * @details The foreground bumps the request counter and spins on
 * release_complete(); a single ISR service_release() must ack and null the live
 * pointer — the use-after-free guard, since the foreground frees the instance
 * only once the ISR can no longer reach it.
 */
inline void test_release_handshake() {
  EffectHandoff<FakeEffect> h;
  FakeEffect e;
  h.adopt(&e, 7);
  HS_EXPECT_EQ(h.live(), &e);

  h.request_release();
  HS_EXPECT_FALSE(h.release_complete()); // req bumped, not yet acked
  HS_EXPECT_EQ(h.live(), &e);            // ISR still holds it until serviced

  h.service_release();
  HS_EXPECT_TRUE(h.release_complete());
  HS_EXPECT_EQ(h.live(), nullptr); // dropped before the foreground frees it

  // A second service with no outstanding request is a no-op.
  h.service_release();
  HS_EXPECT_TRUE(h.release_complete());
  HS_EXPECT_EQ(h.live(), nullptr);
}

/**
 * @brief Two release requests before a single service still reconcile.
 * @details ISR wakes are advisory; several foreground requests may accumulate
 * before one service runs. The ack copies the whole request count, so one
 * service clears the backlog.
 */
inline void test_release_backlog_reconciles() {
  EffectHandoff<FakeEffect> h;
  FakeEffect e;
  h.adopt(&e, 1);
  h.request_release();
  h.request_release();
  HS_EXPECT_FALSE(h.release_complete());
  h.service_release();
  HS_EXPECT_TRUE(h.release_complete());
  HS_EXPECT_EQ(h.live(), nullptr);
}

/**
 * @brief Commit path: publish then adopt with a matching generation.
 * @details committable() is the commit-time use-after-free guard — true only
 * when the pending slot holds an effect whose generation matches the wire's
 * advertised build.
 */
inline void test_commit_adopt_matching_gen() {
  EffectHandoff<FakeEffect> h;
  FakeEffect e;
  h.publish(&e, 42);

  const auto p = h.pending_acquire();
  HS_EXPECT_EQ(p.effect, &e);
  HS_EXPECT_EQ(p.gen, 42u);
  HS_EXPECT_TRUE(h.committable(p, 42));
  HS_EXPECT_FALSE(h.committable(p, 43)); // stale wire generation

  h.adopt(p.effect, p.gen);
  HS_EXPECT_EQ(h.live(), &e);
  HS_EXPECT_TRUE(h.consumed(42));
  HS_EXPECT_FALSE(h.consumed(41));
}

/**
 * @brief Commit guard trips on an empty or stale pending slot.
 * @details After clear_pending() the slot is null, so committable() is false —
 * the commit-time HS_CHECK that would otherwise trap on a use-after-free.
 */
inline void test_commit_guard_rejects_empty_and_stale() {
  EffectHandoff<FakeEffect> h;
  const auto empty = h.pending_acquire();
  HS_EXPECT_EQ(empty.effect, nullptr);
  HS_EXPECT_FALSE(h.committable(empty, 0));

  FakeEffect e;
  h.publish(&e, 5);
  h.clear_pending();
  const auto cleared = h.pending_acquire();
  HS_EXPECT_EQ(cleared.effect, nullptr);
  HS_EXPECT_FALSE(h.committable(cleared, 5));
}

/**
 * @brief Join path: adopt only a present, unconsumed, wire-matching generation.
 * @details A late joiner takes the pending effect live, but only when its
 * generation still matches the wire and has not already been consumed; a
 * mismatch simply waits for the next join grid step.
 */
inline void test_join_adopt_and_gen_gating() {
  EffectHandoff<FakeEffect> h;
  FakeEffect e;
  h.publish(&e, 9);

  const auto p = h.pending_acquire();
  HS_EXPECT_FALSE(h.joinable(p, 8)); // wire advertises a different generation
  HS_EXPECT_TRUE(h.joinable(p, 9));

  h.adopt(p.effect, p.gen);
  HS_EXPECT_EQ(h.live(), &e);

  // Already consumed: the same generation must not re-adopt.
  const auto again = h.pending_acquire();
  HS_EXPECT_FALSE(h.joinable(again, 9));
}

/**
 * @brief Full teardown→publish→commit cycle across two generations.
 * @details Mirrors the run_show/flywheel_isr loop single-threaded: adopt gen1,
 * tear it down via the handshake, clear, publish gen2, commit gen2. Exercises
 * the ordering the device relies on without a live effect ever being adopted
 * while a release is outstanding.
 */
inline void test_full_handoff_cycle() {
  EffectHandoff<FakeEffect> h;
  FakeEffect e1, e2;

  // Generation 1 goes live.
  h.publish(&e1, 1);
  auto p1 = h.pending_acquire();
  HS_EXPECT_TRUE(h.committable(p1, 1));
  h.adopt(p1.effect, p1.gen);
  HS_EXPECT_EQ(h.live(), &e1);
  HS_EXPECT_TRUE(h.consumed(1));

  // Foreground tears e1 down before freeing it.
  h.request_release();
  h.service_release();
  HS_EXPECT_TRUE(h.release_complete());
  HS_EXPECT_EQ(h.live(), nullptr);
  h.clear_pending();
  HS_EXPECT_EQ(h.pending_acquire().effect, nullptr);

  // Generation 2 is built, published, and committed.
  h.publish(&e2, 2);
  auto p2 = h.pending_acquire();
  HS_EXPECT_EQ(p2.effect, &e2);
  HS_EXPECT_TRUE(h.committable(p2, 2));
  h.adopt(p2.effect, p2.gen);
  HS_EXPECT_EQ(h.live(), &e2);
  HS_EXPECT_TRUE(h.consumed(2));
}

/**
 * @brief Display-window (clip) alternation publishes the swept half.
 * @details The window defaults to arm-A-left (1); a ZERO-crossing flip keeps it
 * 1, a HALF flip clears it to 0 — the value the foreground reads to clip the
 * next frame to the opposite quadrant.
 */
inline void test_window_alternation() {
  EffectHandoff<FakeEffect> h;
  HS_EXPECT_EQ(h.window_left(), 1u); // default: arm-A-left
  h.set_window_left(false);          // HALF flip
  HS_EXPECT_EQ(h.window_left(), 0u);
  h.set_window_left(true); // ZERO flip
  HS_EXPECT_EQ(h.window_left(), 1u);
}

/**
 * @brief Module entry point: run the segmented-POV cases.
 * @return Number of failures recorded by the module.
 */
inline int run_pov_segmented_tests() {
  hs_test::ModuleFixture fixture("pov_segmented");

  test_segment_derivation();
  test_eight_segment_derivation();
  test_segment_id_straps();
  test_arm_b_offset();
  test_segment_clip();

  test_release_handshake();
  test_release_backlog_reconciles();
  test_commit_adopt_matching_gen();
  test_commit_guard_rejects_empty_and_stale();
  test_join_adopt_and_gen_gating();
  test_full_handoff_cycle();
  test_window_alternation();

  // Full-canvas tiling across representative rotation columns, including the
  // x=0 and x=w/2 frame boundaries and the wrap seam.
  for (int x : {0, 1, 143, 144, 287})
    check_tiling(/*S=*/288, /*N=*/4, /*w=*/288, x);

  for (int x : {0, 1, 143, 144, 287})
    check_tiling(/*S=*/288, /*N=*/8, /*w=*/288, x);

  // N=2: a single (non-reversed) segment owns each arm's full column.
  for (int x : {0, 1, 47, 48, 95})
    check_tiling(/*S=*/288, /*N=*/2, /*w=*/96, x);

  // Small config swept over every rotation column: S=8, N=4 -> ROWS=4, PPS=2.
  for (int x = 0; x < 8; ++x)
    check_tiling(/*S=*/8, /*N=*/4, /*w=*/8, x);

  for (int x = 0; x < 16; ++x)
    check_tiling(/*S=*/16, /*N=*/8, /*w=*/16, x);

  return fixture.result();
}

} // namespace pov_segmented_tests
} // namespace hs_test
