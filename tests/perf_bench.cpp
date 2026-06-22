/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host render micro-benchmark: per-frame wall time for every effect.
 *
 * Purpose: bracket the cost of the size build's per-pixel codegen. Build the two
 * targets (perf_bench_os at -Os = the Phantasm size config, perf_bench_o3 at
 * -O3) and diff: the gap is what optimization buys the render (filter-plot
 * inlining + per-pixel math vectorization; the FunctionRef pipeline indirection
 * is un-inlinable in BOTH, so it is NOT in the gap). Host x86-64 wall time is a
 * relative proxy for the Teensy device, not an absolute — use it to compare
 * configs and to catch per-effect render regressions, not as a device frame time.
 *
 * Not a CTest: wall time is environment-dependent and would flake a CI gate.
 */
#include <chrono>
#include <cstdint>
#include <cstdio>

#include "effects.h"
#include "memory.h"

namespace {
constexpr int W = 288, H = 144;
constexpr int kFrames = 60, kWarm = 10;

template <typename Effect> double bench(const char *name) {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;
  Effect effect;
  effect.init();
  for (int i = 0; i < kWarm; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < kFrames; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  auto t1 = std::chrono::steady_clock::now();
  double us =
      std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() /
      1000.0 / kFrames;
  volatile auto px = effect.get_pixel(0, 0); // keep the render live
  (void)px;
  std::printf("  %-22s %9.1f us/frame\n", name, us);
  return us;
}
} // namespace

int main() {
#ifdef __OPTIMIZE_SIZE__
  const char *opt = "-Os";
#else
  const char *opt = "-O3";
#endif
  std::printf("=== perf_bench [%s] us/frame, %dx%d (%d frames) ===\n", opt, W, H,
              kFrames);
  double total = 0;
#define HS_BENCH_ONE(name) total += bench<name<W, H>>(#name);
  HS_EFFECT_LIST(HS_BENCH_ONE)
#undef HS_BENCH_ONE
  std::printf("  %-22s %9.1f us/frame (sum)\n", "TOTAL", total);
  return 0;
}
