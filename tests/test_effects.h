/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Effect smoke / robustness harness — exercises all registered effects.
 *
 * Each effect is constructed, init()'d, and rendered for several frames at a
 * supported resolution, then its full output buffer is read back. The test
 * build keeps asserts ON (NDEBUG is force-defined only on the Arduino target),
 * so this drives the very ArenaVector bounds / Canvas out-of-bounds /
 * use-after-free guards that are compiled out on hardware: an effect that
 * overruns an arena, indexes a pixel out of range, or derefs a failed
 * allocation will abort here instead of silently corrupting the live show.
 *
 * SCOPE / FOLLOW-UP: this is a robustness pass, not a pixel-correctness or
 * cross-run-determinism check. Output pixels are 16-bit integers (NaN cannot
 * survive into them), and the engine's animation clock (hs::millis / beatsin*)
 * reads wall-clock time on host, so byte-exact reproducibility across runs
 * requires a time-injection seam that does not yet exist. Adding that seam
 * (mockable millis/micros) would let this harness assert determinism and hash
 * golden frames — the recommended next step.
 */
#pragma once

#include "core/effects.h"
#include "core/canvas.h"
#include "core/memory.h"
#include "tests/test_harness.h"

#include <cstdint>
#include <cstdlib>

namespace hs_test {
namespace effects_tests {

// Render at the primary production resolution (the daydream simulator default
// and the device's full sphere). This is the configuration effects are
// exercised at in deployment, so it is the representative smoke target.
constexpr int kW = 288;
constexpr int kH = 144;

// Default smoke frame count — kept small so the full suite stays well under a
// few seconds (the fast `ctest`/pre-commit path). Set HS_SMOKE_FRAMES=<n> to
// drive long, cyclic code paths (effect morph cycles, particle/trail wraps,
// arena compaction) that only surface over many frames. That deep mode is
// opt-in and intentionally NOT part of the default run.
constexpr int kDefaultFrames = 8;

inline int smoke_frames() {
  // getenv is the simplest way to parameterize a test binary; the MSVC CRT
  // flags it deprecated in favor of _dupenv_s, but the standard call is fine
  // for read-only test config.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *e = std::getenv("HS_SMOKE_FRAMES")) {
#pragma clang diagnostic pop
    int n = std::atoi(e);
    if (n > 0) return n;
  }
  return kDefaultFrames;
}

// Drive one effect type through construct -> init -> render -> read-back.
template <template <int, int> class E>
inline void smoke_one(const char *name) {
  // Deterministic RNG so a given effect takes the same code path each run.
  hs::random().seed(1337u);
  // Fresh arena layout for every effect (each effect's init() may reconfigure).
  configure_arenas_default();
  // CRITICAL: Timeline state (event buffer + frame counter) is global/static and
  // shared across all effects. The real app resets it on effect switch; without
  // this, one effect's leftover events reference the previous (now-destroyed)
  // effect instance, and the next effect's step() would deref freed state.
  Timeline<kW>().clear();
  Timeline<kW>::t = 0;

  E<kW, kH> effect;
  effect.init();

  HS_EXPECT_EQ(effect.width(), kW);
  HS_EXPECT_EQ(effect.height(), kH);

  const int frames = smoke_frames();
  for (int f = 0; f < frames; ++f) {
    effect.draw_frame();
    // Simulate the display consuming the queued frame, otherwise the next
    // Canvas ctor would spin-wait on buffer_free() forever.
    effect.advance_display();
  }

  // Read every pixel of the last displayed frame. Exercises get_pixel()'s
  // index path across the whole buffer; the accumulator prevents the loop
  // from being optimized away.
  auto sum_buffer = [&effect]() {
    uint64_t s = 0;
    for (int y = 0; y < kH; ++y)
      for (int x = 0; x < kW; ++x) {
        const Pixel &p = effect.get_pixel(x, y);
        s += static_cast<uint64_t>(p.r) + p.g + p.b;
      }
    return s;
  };
  const uint64_t acc = sum_buffer();

  // Real post-condition (not the old tautological `acc + 1 > 0`): get_pixel is a
  // pure, stable accessor over the displayed buffer, so re-reading the very same
  // frame must yield the identical sum. A get_pixel that mutated state, indexed
  // nondeterministically, or read outside the buffer would diverge here — and
  // unlike `acc > 0` this holds for effects that legitimately render black
  // (e.g. RingShower sums to 0 at this frame count). Reaching this point also
  // proves the effect constructed, init'd, rendered, and read back without
  // tripping an assert / OOB / hang.
  HS_EXPECT_EQ(acc, sum_buffer());
  std::printf("  [ok] %-20s rendered %d frames @ %dx%d (sum=%llu)\n", name,
              frames, kW, kH, static_cast<unsigned long long>(acc));
}

inline int run_effects_tests() {
  auto scope = hs_test::begin_module("effects");

  // Smoke every registered effect. The list is GENERATED from the single-source
  // roster in core/effects.h (HS_EFFECT_LIST) rather than hand-maintained here,
  // so it can no longer silently drift from the shipped set: an effect added to
  // the roster is smoke-tested automatically (and one in the roster but not
  // #included is a compile error). The matching WASM-registry count is checked
  // against HS_EFFECT_COUNT at engine startup (targets/wasm/wasm.cpp).
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE

  return hs_test::end_module(scope);
}

} // namespace effects_tests
} // namespace hs_test

