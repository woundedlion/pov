/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host arena high-water-mark probe across every effect.
 *
 * Reads Arena::get_lifetime_high_water_mark() for the three global arenas after
 * running each effect's init + hs_test::smoke_frames() frames, to size the
 * device GLOBAL_ARENA_SIZE against what effects actually touch. The lifetime
 * peak, not the windowed one: an effect that re-splits the arenas mid-window
 * (IslamicStars, per shape spawn) discards the windowed mark on every re-split,
 * leaving the single post-window sample reporting only the peak since the last
 * one. The host build uses an
 * 8 MiB host global arena so nothing OOMs mid-measure. The
 * window is HS_SMOKE_FRAMES (default 8); CI drives the 120-frame window the
 * effects sweep uses, which is where late-lifecycle allocation (slot reuse,
 * FIFO expiry, arena compaction) reaches its high-water mark. Frames advance
 * under the injected fixed-cadence clock (hs_test::pin_frame_clock), so a
 * time-driven effect allocates the same way whatever the runner's speed.
 *
 * CI gate: fails (non-zero exit) if the worst single-effect total (persistent +
 * both scratch arenas, the three partitions of the one device pool) exceeds
 * DEVICE_GLOBAL_ARENA_SIZE. The host is 64-bit, so pointer-bearing pooled
 * headers inflate this figure above the 32-bit device (memory.h) — the gate is
 * therefore a conservative upper bound: a host pass guarantees the device fits.
 */
#include <cstdint>
#include <cstdio>

#include "engine/effects.h"
#include "engine/memory.h"
#include "tests/test_fixture.h"

namespace {
constexpr int W = 288, H = 144;
const int FRAMES = hs_test::smoke_frames();

size_t g_max_p = 0, g_max_a = 0, g_max_b = 0, g_worst_total = 0;
const char *g_worst_name = "";
int g_measured = 0;

template <typename Effect> void measure(const char *name) {
  ++g_measured;
  hs_test::reset_globals();
  hs_test::pin_frame_clock(0);
  persistent_arena.reset_peak_tracking();
  scratch_arena_a.reset_peak_tracking();
  scratch_arena_b.reset_peak_tracking();

  Effect effect;
  effect.init();
  for (int f = 0; f < FRAMES; ++f) {
    hs_test::pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
  }
  volatile auto px = effect.get_pixel(0, 0); // keep the render live
  (void)px;

  size_t p = persistent_arena.get_lifetime_high_water_mark();
  size_t a = scratch_arena_a.get_lifetime_high_water_mark();
  size_t b = scratch_arena_b.get_lifetime_high_water_mark();
  size_t tot = p + a + b;
  std::printf(
      "  %-22s persist=%7zu  scratchA=%6zu  scratchB=%6zu  total=%7zu\n", name,
      p, a, b, tot);
  if (p > g_max_p)
    g_max_p = p;
  if (a > g_max_a)
    g_max_a = a;
  if (b > g_max_b)
    g_max_b = b;
  if (tot > g_worst_total) {
    g_worst_total = tot;
    g_worst_name = name;
  }
}
} // namespace

int main() {
  std::printf("=== arena high-water per effect (host -Os, %dx%d, %d frames) "
              "===\n",
              W, H, FRAMES);
#define HS_ARENA_ONE(name) measure<name<W, H>>(#name);
  HS_EFFECT_LIST(HS_ARENA_ONE)
#undef HS_ARENA_ONE
  if (g_measured != HS_EFFECT_COUNT) {
    std::printf(
        "measured %d effects but HS_EFFECT_COUNT = %d — roster empty or "
        "measure() calls dropped\n",
        g_measured, HS_EFFECT_COUNT);
    return 1;
  }
  std::printf("\nMAX persist=%zu  scratchA=%zu  scratchB=%zu\n", g_max_p,
              g_max_a, g_max_b);
  std::printf("WORST single-effect total = %zu B (%s)\n", g_worst_total,
              g_worst_name);
  std::printf("sum of per-arena maxima  = %zu B   (device GLOBAL = %zu B)\n",
              g_max_p + g_max_a + g_max_b, DEVICE_GLOBAL_ARENA_SIZE);
  const bool over = g_worst_total > DEVICE_GLOBAL_ARENA_SIZE;
  std::printf("WORST total %zu B   budget %zu B   [%s]\n", g_worst_total,
              DEVICE_GLOBAL_ARENA_SIZE, over ? "FAIL" : "PASS");
  if (over)
    std::printf("  arena budget exceeded — %s outgrew the device global arena; "
                "see tests/arena_measure.cpp header.\n",
                g_worst_name);
  return over ? 1 : 0;
}
