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
 * COVERAGE: two passes over the full effect roster.
 *   1. smoke_one  — construct/init/render/read-back robustness (above).
 *   2. determinism_one — renders each effect twice under an injected, fixed
 *      per-frame clock (hs::set_mock_time, the seam in core/platform.h that
 *      hs::millis/micros and beatsin* honor) and asserts the two final frames
 *      are byte-identical. This proves *what* is rendered is reproducible — the
 *      gap the in-run get_pixel-stability check cannot cover — and directly
 *      guards the time/state-dependent defect class (drifting trail clocks,
 *      stale globals, uninitialized reads) without any stored reference frame.
 *
 * Deliberately NOT golden-frame hashing: a stored hash over the actively-tuned
 * generative effects inverts the signal (every intentional look change is a red
 * test), carries zero diagnostic value, and is not even bit-reproducible across
 * the native + WASM targets (fast_* trig / float rounding diverge). Self-
 * referential cross-run comparison catches the real bug classes with none of
 * that maintenance tax. Pixel-aesthetic correctness remains a human judgment.
 */
#pragma once

#include "core/effects.h"
#include "core/canvas.h"
#include "core/memory.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace hs_test {
namespace effects_tests {

/**
 * @brief Primary production render width in pixels.
 * @details The daydream simulator default and the device's full sphere; the
 * configuration effects are exercised at in deployment, so it is the
 * representative smoke target.
 */
constexpr int kW = 288;
/**
 * @brief Primary production render height in pixels.
 * @details Paired with kW for the full-sphere production resolution.
 */
constexpr int kH = 144;

/**
 * @brief Holosphere device render width in pixels.
 * @details The hardware target instantiates every effect at <96,20>. The native
 * suite is the only place that specialization runs under asserts (the device
 * forces NDEBUG; the CI WASM smoke runs assert-free), so a second roster pass at
 * this resolution exercises height-20-specific paths (PhiLUT<20> indexing,
 * small-aspect arena sizing, H_OFFSET interactions) that bypass every
 * assert-enabled layer otherwise. Both are in HS_WASM_RESOLUTIONS (wasm.cpp).
 */
constexpr int kDeviceW = 96;
/**
 * @brief Holosphere device render height in pixels.
 * @details Paired with kDeviceW for the <96,20> device specialization.
 */
constexpr int kDeviceH = 20;

/**
 * @brief Default per-effect smoke frame count.
 * @details Kept small so the effects smoke pass itself stays quick. (The full
 * pre-commit suite still runs ~25 s, dominated by the death-test subprocess
 * spawns and the multi-board sync simulator, not these smoke frames.) Set
 * HS_SMOKE_FRAMES=<n> to drive long, cyclic code paths (effect morph cycles,
 * particle/trail wraps, arena compaction, and the effect-lifecycle transitions
 * — RingShower slot reuse, Thrusters fire/FIFO expiry, ShapeShifter's 48-frame
 * cut) that only surface over many frames. 8 frames never reaches those
 * windows, so CI sets HS_SMOKE_FRAMES=120 (.github/workflows/ci.yml) to
 * exercise them on every push/PR while local commits keep the fast 8-frame
 * path.
 */
constexpr int kDefaultFrames = 8;

/**
 * @brief Resolves the per-effect frame count from the environment.
 * @return HS_SMOKE_FRAMES if set to a positive int, else kDefaultFrames.
 */
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

/**
 * @brief Forward declaration of the registered-but-unread param lint.
 * @param effect Effect instance whose editable params are probed.
 * @param name Effect name used in diagnostic output.
 */
inline void lint_dead_sliders(Effect &effect, const char *name);

/**
 * @brief Drives one effect type through construct -> init -> render -> read-back.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to kW).
 * @tparam H Render height in pixels (defaults to kH).
 * @param name Effect name used in the [ok] / diagnostic output.
 * @details Verifies the effect constructs, init's, renders smoke_frames()
 * frames, and reads back every pixel without tripping an assert/OOB/hang, and
 * that get_pixel is a stable pure accessor. Runs the dead-slider lint once on
 * the primary <kW,kH> pass.
 */
template <template <int, int> class E, int W = kW, int H = kH>
inline void smoke_one(const char *name) {
  // Deterministic RNG so a given effect takes the same code path each run.
  hs::random().seed(1337u);
  // Fresh arena layout for every effect (each effect's init() may reconfigure).
  configure_arenas_default();
  // CRITICAL: Timeline state (event buffer + frame counter) is global/static and
  // shared across all effects. The real app resets it on effect switch; without
  // this, one effect's leftover events reference the previous (now-destroyed)
  // effect instance, and the next effect's step() would deref freed state.
  Timeline().clear();
  global_timeline_t = 0;

  E<W, H> effect;
  effect.init();

  HS_EXPECT_EQ(effect.width(), W);
  HS_EXPECT_EQ(effect.height(), H);

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
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const Pixel &p = effect.get_pixel(x, y);
        s += static_cast<uint64_t>(p.r) + p.g + p.b;
      }
    return s;
  };
  const uint64_t acc = sum_buffer();

  // Post-condition: get_pixel is a pure, stable accessor over the displayed
  // buffer, so re-reading the very same frame must yield the identical sum. A
  // get_pixel that mutated state, indexed nondeterministically, or read outside
  // the buffer would diverge here. Unlike `acc > 0` this also holds for effects
  // that legitimately render black (e.g. RingShower sums to 0 at this frame
  // count). Reaching this point also proves the effect constructed, init'd,
  // rendered, and read back without tripping an assert / OOB / hang.
  HS_EXPECT_EQ(acc, sum_buffer());
  std::printf("  [ok] %-20s rendered %d frames @ %dx%d (sum=%llu)\n", name,
              frames, W, H, static_cast<unsigned long long>(acc));

  // The dead-slider lint is resolution-independent (param persistence), so run
  // it once on the primary pass rather than redundantly at every resolution.
  if constexpr (W == kW && H == kH)
    lint_dead_sliders(effect, name);
}

/**
 * @brief Build-time "registered-but-unread" lint for the live-art param system.
 * @param effect Effect instance whose editable params are probed.
 * @param name Effect name used in DEAD SLIDER diagnostic output.
 * @details Contract: a registered, editable param — one NOT flagged
 * markAnimated() and NOT flagged markReadonly() — must be genuinely
 * user-controllable, i.e. a value written through updateParameter() must persist
 * across frames. If the engine overwrites it every frame (a Mutation/Driver/Lerp
 * bound to the same member, or output-only telemetry), the slider is dead: the
 * author must drive a private member and markAnimated() it, or markReadonly()
 * pure telemetry. This is the build-time gate for the theme-4 dead-slider class
 * (it catches the per-frame overwrite mechanism behind
 * MobiusGrid/ShapeShifter/the preset-lerp group; a value that merely has no
 * rendered effect can't be detected without flaky golden-image diffing and is
 * out of scope).
 */
inline void lint_dead_sliders(Effect &effect, const char *name) {
  for (const auto &def : effect.getParameters()) {
    if (def.is_bool() || def.animated || def.readonly)
      continue;
    const float range = def.max - def.min;
    if (range <= 0.0f)
      continue;
    const float cur = def.get();
    // An in-range target well clear of the current value, so a revert is visible.
    const float target = (cur - def.min) > (def.max - cur)
                             ? def.min + 0.25f * range
                             : def.min + 0.75f * range;
    effect.updateParameter(def.name, target);
    for (int f = 0; f < 3; ++f) {
      effect.draw_frame();
      effect.advance_display();
    }
    const float eps = fmaxf(1e-3f, 1e-3f * range);
    const bool persisted = fabsf(def.get() - target) <= eps;
    if (!persisted)
      std::printf("  DEAD SLIDER %s::%s — wrote %.4f, engine reverted to %.4f "
                  "(markAnimated / markReadonly / drive a private member)\n",
                  name, def.name, static_cast<double>(target),
                  static_cast<double>(def.get()));
    HS_EXPECT(persisted, "editable param must persist across frames");
  }
}

/**
 * @brief Per-frame clock advance in milliseconds for the determinism pass (~30fps).
 * @details Identical across both runs, so any frame-to-frame animation driven by
 * hs::millis/micros or beatsin* is reproduced exactly rather than tracking the
 * wall clock.
 */
constexpr unsigned long kFrameMs = 33;
/**
 * @brief Per-frame clock advance in microseconds, paired with kFrameMs.
 */
constexpr unsigned long kFrameUs = 33000;

/**
 * @brief Renders one effect under the injected clock and copies the final buffer out.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to kW).
 * @tparam H Render height in pixels (defaults to kH).
 * @param out Receives the final displayed frame, sized W*H pixels (row-major).
 * @param frames Number of frames to render before capture.
 * @details Resets every shared global the smoke path does (RNG seed, arenas,
 * Timeline) plus the generative-hue cursor and mock clock, so two calls start
 * from an identical state.
 */
template <template <int, int> class E, int W = kW, int H = kH>
inline void render_capture(std::vector<Pixel> &out, int frames) {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;
  // The global generative-hue cursor drifts across palette constructions by
  // design (see GenerativePalette); pin it so both runs start from identical
  // global state — without this, the second run's palettes are hue-rotated.
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * kFrameMs,
                      static_cast<unsigned long>(f) * kFrameUs);
    effect.draw_frame();
    effect.advance_display();
  }

  out.resize(static_cast<size_t>(W) * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      out[static_cast<size_t>(y) * W + x] = effect.get_pixel(x, y);
}

/**
 * @brief Cross-run determinism test: renders an effect twice and requires byte-identical frames.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to kW).
 * @tparam H Render height in pixels (defaults to kH).
 * @param name Effect name used in the NONDETERMINISTIC diagnostic output.
 * @details The clock seam neutralizes wall-time, so a divergence here is real
 * nondeterminism (uninitialized read, stale global, address-dependent path) —
 * the defect class smoke coverage cannot see.
 */
template <template <int, int> class E, int W = kW, int H = kH>
inline void determinism_one(const char *name) {
  const int frames = smoke_frames();
  std::vector<Pixel> a, b;
  render_capture<E, W, H>(a, frames);
  render_capture<E, W, H>(b, frames);
  hs::clear_mock_time();

  int first_diff = -1;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i].r != b[i].r || a[i].g != b[i].g || a[i].b != b[i].b) {
      first_diff = static_cast<int>(i);
      break;
    }

  if (first_diff >= 0) {
    const int x = first_diff % W, y = first_diff / W;
    std::printf("  NONDETERMINISTIC %-20s pixel (%d,%d): runA (%d,%d,%d) != "
                "runB (%d,%d,%d) over %d frames\n",
                name, x, y, static_cast<int>(a[first_diff].r),
                static_cast<int>(a[first_diff].g),
                static_cast<int>(a[first_diff].b),
                static_cast<int>(b[first_diff].r),
                static_cast<int>(b[first_diff].g),
                static_cast<int>(b[first_diff].b), frames);
  }
  HS_EXPECT(first_diff < 0,
            "effect must render identically across runs under a fixed clock");
}

/**
 * @brief Verifies SHMath::decode_lm yields a valid spherical-harmonic order for every flat index.
 * @details Requires l = floor(sqrt(idx)) and m in [-l, l], with idx == l*l + l +
 * m. The risk is float sqrtf rounding at perfect squares truncating l low and
 * pushing m out of band; sweeps a wide range (well past the effect's idx<=23) so
 * every perfect square and its neighbours — including the large ones where sqrtf
 * actually misrounds — are exercised.
 */
inline void test_sh_decode_lm_valid_order() {
  auto check = [](int idx) {
    auto [l, m] = SHMath::decode_lm(idx);
    HS_EXPECT_EQ(l * l + l + m, idx);   // round-trips the index
    HS_EXPECT_TRUE(m >= -l && m <= l);  // order within the level's band
    HS_EXPECT_TRUE(l * l <= idx && idx < (l + 1) * (l + 1)); // l is the floor
  };

  // Dense sweep over the range the effect actually uses and well beyond.
  for (int idx = 0; idx <= 4096; ++idx)
    check(idx);

  // Targeted large indices where float sqrtf genuinely misrounds: once k^2
  // exceeds 2^24 it is no longer exactly representable, and for k beyond ~2900 a
  // just-below-square index sqrt-rounds *up* to k. Both would truncate l one
  // level low and push m out of [-l, l] if decode_lm did not compensate. (k kept
  // under ~46340 so (l+1)^2 stays within int.)
  for (long long k : {3000LL, 5000LL, 12345LL, 40000LL}) {
    long long sq = k * k;
    for (long long idx : {sq - 1, sq, sq + 1, sq + k, sq + 2 * k})
      check(static_cast<int>(idx));
  }
}

/**
 * @brief Module entry point for the effects test suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Runs the SH-decode check, then both smoke and determinism passes over
 * the full effect roster at the production and device resolutions.
 */
inline int run_effects_tests() {
  auto scope = hs_test::begin_module("effects");

  test_sh_decode_lm_valid_order();

  // Smoke every registered effect. The list is GENERATED from the single-source
  // roster in core/effects.h (HS_EFFECT_LIST), so it cannot drift from the
  // shipped set: an effect added to the roster is smoke-tested automatically
  // (and one in the roster but not #included is a compile error). The matching
  // WASM-registry count is checked against HS_EFFECT_COUNT at engine startup
  // (targets/wasm/wasm.cpp).
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE

  // Second pass: cross-run determinism under the injected clock. Same roster,
  // so a new effect is automatically held to byte-exact reproducibility too.
#define HS_DET_ONE(name) determinism_one<name>(#name);
  HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE

  // Repeat both passes at the Holosphere device resolution <96,20>, the only
  // place that specialization runs under native asserts (see kDeviceW/kDeviceH).
  std::printf("  -- device resolution %dx%d --\n", kDeviceW, kDeviceH);
#define HS_SMOKE_ONE_DEV(name) smoke_one<name, kDeviceW, kDeviceH>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE_DEV)
#undef HS_SMOKE_ONE_DEV
#define HS_DET_ONE_DEV(name) determinism_one<name, kDeviceW, kDeviceH>(#name);
  HS_EFFECT_LIST(HS_DET_ONE_DEV)
#undef HS_DET_ONE_DEV

  return hs_test::end_module(scope);
}

} // namespace effects_tests
} // namespace hs_test

