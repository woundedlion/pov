/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host arena high-water-mark probe across every effect.
 *
 * Reads Arena::get_high_water_mark() for the three global arenas after running
 * each effect's init + a few frames, to size the device GLOBAL_ARENA_SIZE
 * (343,040 B = 335 KiB) against what effects actually touch. The host build uses
 * an 8 MB global arena (memory.h, HS_TEST_BUILD) so nothing OOMs mid-measure.
 *
 * Host caveat vs device: arena payloads are mostly POD geometry (Vector =
 * 3 floats, uint16 indices — same size both), so the figure tracks device
 * closely; the few pointer-bearing structures make host a slight over-estimate
 * (8-byte vs 4-byte pointers) — conservative for sizing. Not a CTest.
 */
#include <cstdint>
#include <cstdio>

#include "effects.h"
#include "memory.h"

namespace {
constexpr int W = 288, H = 144;
constexpr int kFrames = 8;

size_t g_max_p = 0, g_max_a = 0, g_max_b = 0, g_worst_total = 0;
const char *g_worst_name = "";

template <typename Effect> void measure(const char *name) {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;
  // Zero the peaks after the (re)configure so we capture only this effect's
  // init + frame allocations, including any per-effect configure_arenas() split.
  persistent_arena.reset_high_water_mark();
  scratch_arena_a.reset_high_water_mark();
  scratch_arena_b.reset_high_water_mark();

  Effect effect;
  effect.init();
  for (int f = 0; f < kFrames; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  volatile auto px = effect.get_pixel(0, 0); // keep the render live
  (void)px;

  size_t p = persistent_arena.get_high_water_mark();
  size_t a = scratch_arena_a.get_high_water_mark();
  size_t b = scratch_arena_b.get_high_water_mark();
  size_t tot = p + a + b;
  std::printf("  %-22s persist=%7zu  scratchA=%6zu  scratchB=%6zu  total=%7zu\n",
              name, p, a, b, tot);
  if (p > g_max_p) g_max_p = p;
  if (a > g_max_a) g_max_a = a;
  if (b > g_max_b) g_max_b = b;
  if (tot > g_worst_total) { g_worst_total = tot; g_worst_name = name; }
}
} // namespace

int main() {
  std::printf("=== arena high-water per effect (host -Os, %dx%d) ===\n", W, H);
#define HS_ARENA_ONE(name) measure<name<W, H>>(#name);
  HS_EFFECT_LIST(HS_ARENA_ONE)
#undef HS_ARENA_ONE
  std::printf("\nMAX persist=%zu  scratchA=%zu  scratchB=%zu\n", g_max_p,
              g_max_a, g_max_b);
  std::printf("WORST single-effect total = %zu B (%s)\n", g_worst_total,
              g_worst_name);
  std::printf("sum of per-arena maxima  = %zu B   (device GLOBAL = %zu B)\n",
              g_max_p + g_max_a + g_max_b, DEVICE_GLOBAL_ARENA_SIZE);
  return 0;
}
