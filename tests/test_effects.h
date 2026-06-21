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

// ---------------------------------------------------------------------------
// Gray-Scott reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for GSReactionDiffusion's private fixed-point and
 *        physics internals (befriended in effects/GSReactionDiffusion.h).
 * @details The generic smoke/determinism harness only proves the effect doesn't
 *          crash and renders reproducibly — it cannot see a Q16 round-trip error
 *          or a numerically-unstable-but-deterministic blow-up. This seam pins
 *          the conversions and one Gray-Scott substep directly. The lattice is a
 *          fixed 7680-node graph, independent of <W,H>, so the device resolution
 *          is used arbitrarily.
 */
struct GSWhiteBox {
  using GS = GSReactionDiffusion<kDeviceW, kDeviceH>;
  static constexpr int N = GS::RD_N;

  static uint16_t to_q16(float v) { return GS::to_q16(v); }
  static float from_q16(uint16_t v) { return GS::from_q16(v); }
  static void set_params(GS &gs, float feed, float k, float dA, float dB,
                         float dt) {
    gs.params.feed = feed;
    gs.params.k = k;
    gs.params.dA = dA;
    gs.params.dB = dB;
    gs.params.dt = dt;
  }
  static void step(GS &gs, const uint16_t *cA, const uint16_t *cB, uint16_t *nA,
                   uint16_t *nB) {
    gs.step_physics(cA, cB, nA, nB);
  }
};

/**
 * @brief Verifies the Q16 fixed-point round-trip and the documented +0.5
 *        rounding/clamp boundaries.
 * @details to_q16(from_q16(v)) must be the identity over every representable
 *          value, and to_q16 must clamp out-of-range floats and round to nearest
 *          (so 1.0 tops out at 65535 with no overflow). A truncating or
 *          unclamped regression — which would bias the RD dynamics — fails here.
 */
inline void test_gs_q16_roundtrip() {
  HS_EXPECT_EQ(GSWhiteBox::to_q16(0.0f), (uint16_t)0);
  HS_EXPECT_EQ(GSWhiteBox::to_q16(1.0f), (uint16_t)65535);
  HS_EXPECT_EQ(GSWhiteBox::to_q16(2.0f), (uint16_t)65535);  // clamp high
  HS_EXPECT_EQ(GSWhiteBox::to_q16(-0.5f), (uint16_t)0);     // clamp low
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(65535), 1.0f, 1e-9f);
  // Round-trip is exact for every Q16 value (the +0.5 absorbs the float error).
  int bad = 0;
  for (int v = 0; v <= 65535; ++v)
    if (GSWhiteBox::to_q16(GSWhiteBox::from_q16((uint16_t)v)) != (uint16_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies the homogeneous rest state (A=1, B=0) is a Gray-Scott fixed
 *        point: one substep leaves it exactly unchanged.
 * @details With B=0 the reaction term A·B² and both Laplacians vanish and the
 *          feed term feed·(1-A) is zero at A=1, so the state must not move. A
 *          sign error or stray term in the update assembly perturbs it here.
 */
inline void test_gs_rest_state_is_fixed_point() {
  std::vector<uint16_t> cA(GSWhiteBox::N, 65535), cB(GSWhiteBox::N, 0),
      nA(GSWhiteBox::N), nB(GSWhiteBox::N);
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f); // defaults
  GSWhiteBox::step(gs, cA.data(), cB.data(), nA.data(), nB.data());
  int moved = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (nA[i] != 65535 || nB[i] != 0)
      ++moved;
  HS_EXPECT_EQ(moved, 0);
}

/**
 * @brief Verifies one substep has the right reaction/diffusion signs and that
 *        the Q16 clamp is actually applied.
 * @details Seed a single saturated-B nucleus on the otherwise-rest field. After
 *          one default-parameter step: A at the seed is fully consumed (the
 *          1 - dt update underflows and must clamp to 0, not wrap), and B
 *          diffuses into at least one neighbor that started empty.
 */
inline void test_gs_substep_signs_and_clamp() {
  std::vector<uint16_t> cA(GSWhiteBox::N, 65535), cB(GSWhiteBox::N, 0),
      nA(GSWhiteBox::N), nB(GSWhiteBox::N);
  const int seed = 4000; // an interior lattice node with a full neighbor ring
  cB[seed] = 65535;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::step(gs, cA.data(), cB.data(), nA.data(), nB.data());

  // a + (dA·0 - 1 + feed·0)·dt = 1 - 2.5 < 0 → clamps to 0 (not an unclamped
  // negative-float-to-uint16 wrap).
  HS_EXPECT_EQ((int)nA[seed], 0);
  // B diffuses outward: at least one initially-empty neighbor is now lit.
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && nB[nb] > 0)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0);
}

/**
 * @brief Verifies the explicit-Euler integrator does not diverge over many
 *        substeps at the worst-case stable slider setting.
 * @details At the Speed/dA/dB extremes the stability product dt·D·|λ|max =
 *          3·0.05·12 = 1.8 ≤ 2, so the scheme must stay bounded. A genuine
 *          instability (dt·D·|λ|max > 2) oscillates and the Q16 clamp pins the
 *          oscillating nodes to the 0/65535 rails — a saturate/oscillate blow-up
 *          that renders identically across runs and so slips past the
 *          determinism pass. After 256 steps from seeded nuclei almost no node
 *          should sit at the upper rail; assert that directly. (Whether B
 *          ultimately persists or decays is regime-dependent and not asserted —
 *          high diffusion legitimately dilutes the seeds toward the rest state.)
 */
inline void test_gs_evolution_stays_bounded() {
  std::vector<uint16_t> a(GSWhiteBox::N, 65535), b(GSWhiteBox::N, 0),
      sa(GSWhiteBox::N), sb(GSWhiteBox::N);
  for (int s : {500, 2500, 4500, 6500})
    b[s] = 65535;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.05f, 0.05f, 3.0f); // stable extreme
  uint16_t *cA = a.data(), *cB = b.data(), *nA = sa.data(), *nB = sb.data();
  for (int s = 0; s < 256; ++s) {
    GSWhiteBox::step(gs, cA, cB, nA, nB);
    std::swap(cA, nA);
    std::swap(cB, nB);
  }
  // No blow-up: a stable run leaves at most a handful of nodes at the upper
  // rail; an unstable oscillation would clamp a large fraction of the lattice
  // there. (< 5% is a generous tripwire — the stable run sits near zero.)
  int saturated = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (cB[i] == 65535)
      ++saturated;
  HS_EXPECT_LT(saturated, GSWhiteBox::N / 20);
}

// ---------------------------------------------------------------------------
// DreamBalls: preset-cycle / re-spawn white-box coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for DreamBalls' private preset-cycle bookkeeping
 *        (befriended in effects/DreamBalls.h).
 * @details spawn_sprite schedules its successor 288 frames out (a PeriodicTimer),
 *          but the smoke/determinism harness renders at most HS_SMOKE_FRAMES=120
 *          frames in CI — short of one period — so the re-spawn never fires under
 *          the generic passes: the preset advance, the active_bake_ ping-pong +
 *          rebake, and the reseed-on-change guard all stay dead. This seam drives
 *          spawn_sprite directly and reads the bake slot / preset index so those
 *          paths are pinned. <96,20> is used arbitrarily — the bookkeeping is
 *          resolution-independent.
 */
struct DreamBallsWhiteBox {
  using DB = DreamBalls<kDeviceW, kDeviceH>;
  static constexpr int kPresets = 4;

  static int active_bake(const DB &db) { return db.active_bake_; }
  static int last_preset_idx(const DB &db) { return db.last_preset_idx_; }
  static void spawn(DB &db, int idx) { db.spawn_sprite(idx); }
  // The literal solid-name pointer the given preset seeds into params on reseed.
  static const char *preset_name(const DB &db, int idx) {
    return db.preset_manager.entries[idx].params.solid_name;
  }
};

/**
 * @brief Drives spawn_sprite across a full preset cycle and asserts the bake-slot
 *        ping-pong, the modulo preset advance, and the reseed-on-change guard.
 * @details Calls spawn_sprite directly (no 288-frame wait) following the same idx
 *          progression the periodic callback uses: idx grows unbounded and the
 *          active preset is idx % 4. Each spawn must flip the bake slot (so a
 *          fading-out sprite keeps its own LUT) and, when the preset actually
 *          changes, reseed params to the new entry. A re-spawn of the SAME preset
 *          (the paused branch) must instead hold params so a live slider edit
 *          survives.
 */
inline void test_dreamballs_preset_cycle_bookkeeping() {
  using WB = DreamBallsWhiteBox;
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;

  WB::DB db;
  db.init(); // runs spawn_sprite(0)

  // init() spawned preset 0: it reseeded params (last_preset_idx_ -1 -> 0) and
  // flipped the bake slot once (0 -> 1).
  HS_EXPECT_EQ(WB::last_preset_idx(db), 0);
  HS_EXPECT_EQ(WB::active_bake(db), 1);
  HS_EXPECT_EQ(db.params.solid_name, WB::preset_name(db, 0));

  // Not-paused advance chain: the periodic callback re-invokes spawn_sprite(idx+1),
  // so idx grows unbounded and the preset is idx % 4. Drive a full cycle plus a
  // wrap; the bake slot must ping-pong every spawn and params must reseed to the
  // wrapped index each step.
  int expect_bake = WB::active_bake(db); // 1
  for (int idx = 1; idx <= 8; ++idx) {
    WB::spawn(db, idx);
    expect_bake ^= 1;
    const int safe = idx % WB::kPresets;
    HS_EXPECT_EQ(WB::active_bake(db), expect_bake);
    HS_EXPECT_EQ(WB::last_preset_idx(db), safe);
    HS_EXPECT_EQ(db.params.solid_name, WB::preset_name(db, safe));
  }

  // Paused-hold path: the callback re-spawns the SAME idx, so the reseed guard
  // (safe_idx == last_preset_idx_) holds and a live slider edit must survive the
  // re-spawn — while the bake slot still flips. last_preset_idx_ is now 0
  // (idx 8 % 4); re-spawn idx 8 again with a sentinel edit in place.
  const float sentinel = db.params.num_copies + 5.0f;
  db.params.num_copies = sentinel;
  const int held_idx = WB::last_preset_idx(db);
  expect_bake ^= 1;
  WB::spawn(db, 8); // 8 % 4 == held_idx: the paused re-spawn of the same preset
  HS_EXPECT_EQ(WB::last_preset_idx(db), held_idx); // preset unchanged
  HS_EXPECT_EQ(db.params.num_copies, sentinel);    // live edit preserved
  HS_EXPECT_EQ(WB::active_bake(db), expect_bake);  // bake slot still flipped
}

/**
 * @brief End-to-end check that the 288-frame re-spawn timer actually fires and
 *        honors the pause gate, by rendering past one period.
 * @details Covers the part the white-box driver cannot: the PeriodicTimer wiring
 *          and the animationsPaused()?idx:idx+1 hold-vs-advance decision in the
 *          scheduler lambda. Renders 300 frames (one period is 288; the next
 *          re-spawn is another 288 away, so exactly one fires). Unpaused, the
 *          preset advances 0 -> 1; paused, it re-spawns the same preset and holds.
 */
inline void test_dreamballs_respawn_fires_and_honors_pause() {
  using WB = DreamBallsWhiteBox;
  auto run = [](bool paused) {
    hs::random().seed(1337u);
    configure_arenas_default();
    Timeline().clear();
    global_timeline_t = 0;
    WB::DB db;
    db.init();
    db.setAnimationsPaused(paused);
    // The re-spawn timer is scheduled 288 frames out; render past it so exactly
    // one re-spawn fires (the following one is another 288 frames away).
    for (int f = 0; f < 300; ++f) {
      db.draw_frame();
      db.advance_display();
    }
    return WB::last_preset_idx(db);
  };
  HS_EXPECT_EQ(run(false), 1); // unpaused: re-spawn advanced preset 0 -> 1
  HS_EXPECT_EQ(run(true), 0);  // paused: re-spawn fired but held preset 0
}

// ---------------------------------------------------------------------------
// In-code-flagged numeric invariants with no oracle in the smoke harness
// ---------------------------------------------------------------------------
//
// The smoke/determinism passes only prove each effect renders something
// reproducibly; they cannot detect a wrong-but-still-non-black result. These
// white-box seams pin the fragile numeric invariants the effects flag in-code:
// a regression in any of them renders fine and slips past the generic harness.
//
// (MeshFeedback's SHAPE_FRAMES/PRESET_FRAMES coprimality — the fourth such
// invariant — is locked at compile time by a static_assert in the effect, so
// it needs no runtime case here: instantiating the effect in the smoke pass
// already checks it.)

/**
 * @brief White-box accessor for Comets' Lissajous-loop closing snap.
 * @details Befriended in effects/Comets.h. Reaches the private closing_domain()
 *          snap and the authored function table to verify every entry closes —
 *          path_fn(domain) == path_fn(0) — so the per-cycle drift reset never
 *          teleports the head. A wrong snap still renders a (discontinuous)
 *          curve, invisible to the smoke harness.
 */
struct CometsWhiteBox {
  /**
   * @brief Verifies every authored function table entry closes its loop.
   * @details For each entry the snapped endpoint must coincide with the t=0
   *          start (0,1,0), and the snap must stay positive (the floor-at-1
   *          guard keeps the head moving rather than freezing at path_fn(0)).
   */
  static void check_paths_close() {
    using C = Comets<kW, kH>;
    int idx = 0;
    for (const LissajousParams &cfg : C::functions) {
      const float cd = C::closing_domain(cfg);
      HS_EXPECT_GT(cd, 0.0f); // floor-at-1 keeps the head moving
      const Vector start = lissajous(cfg.m1, cfg.m2, cfg.a, 0.0f);
      const Vector end = lissajous(cfg.m1, cfg.m2, cfg.a, cd);
      const float gap = (end - start).magnitude();
      if (gap > 1e-3f)
        std::printf("  COMETS entry %d does not close: gap=%.6f\n", idx,
                    static_cast<double>(gap));
      HS_EXPECT_LT(gap, 1e-3f);
      ++idx;
    }
  }
};

/**
 * @brief White-box accessor for the Thrusters warp-decay endpoints.
 * @details Befriended in effects/Thrusters.h. Reaches the private warp_decay()
 *          curve to pin its shift-and-renormalized endpoints: a bare
 *          0.7*exp(-2t) would bottom out at ~0.095 and leave a residual wobble,
 *          which still renders and so passes the smoke harness.
 */
struct ThrustersWhiteBox {
  /**
   * @brief Verifies warp_decay peaks at 0.7 at t=0 and relaxes to exactly 0 at t=1.
   */
  static void check_warp_endpoints() {
    using T = Thrusters<kW, kH>;
    HS_EXPECT_NEAR(T::warp_decay(0.0f), 0.7f, 1e-6f); // peak at fire
    HS_EXPECT_NEAR(T::warp_decay(1.0f), 0.0f, 1e-6f); // full relaxation by end
  }
};

/**
 * @brief White-box accessor for RingShower's radius easing endpoints.
 * @details Befriended in effects/RingShower.h. Reaches the private Ring type to
 *          pin its age+1 convention: the ring must reach RADIUS_MAX on its final
 *          visible frame (age+1 == life) and render a non-zero first step rather
 *          than radius 0. An off-by-one in the convention still renders a ring.
 */
struct RingShowerWhiteBox {
  /**
   * @brief Verifies radius_at hits RADIUS_MAX on the final frame and is non-zero
   *        on the first.
   */
  static void check_radius_endpoints() {
    using RS = RingShower<kW, kH>;
    typename RS::Ring ring;
    ring.life = 50;
    ring.age = ring.life - 1; // final visible frame: age+1 == life -> t == 1
    HS_EXPECT_NEAR(ring.radius_at(), RS::Ring::RADIUS_MAX, 1e-5f);
    ring.age = 0; // first frame: one eased step in, not radius 0
    HS_EXPECT_GT(ring.radius_at(), 0.0f);
    HS_EXPECT_LT(ring.radius_at(), RS::Ring::RADIUS_MAX);
  }
};

/**
 * @brief White-box accessor for Dynamo's overlapping-wipe band ordering.
 * @details Befriended in effects/Dynamo.h. color() documents that a live Wipe-Dur
 *          change between two overlapping wipes can let the newer (front, index 0)
 *          boundary overtake the older one, transiently inverting band order — an
 *          acknowledged cosmetic gap the author asserts stays memory-safe and
 *          in-range. This pins that safety claim: it stages exactly that
 *          non-monotonic boundary order and sweeps the full angular span,
 *          asserting every color() call stays in bounds (every baked_palettes_
 *          access is HS_CHECK-guarded, so an OOB aborts here) and returns a finite
 *          alpha in [0, 1]. The smoke pass never reaches the inverted state.
 */
struct DynamoWhiteBox {
  /**
   * @brief Verifies color() is memory-safe and in-range under inverted bands.
   */
  static void check_overlapping_wipes_stay_in_range() {
    // Mirror smoke_one's global-state reset: Dynamo::init() bakes the LUT pool
    // from persistent_arena and schedules timers on the shared global timeline.
    hs::random().seed(1337u);
    configure_arenas_default();
    Timeline().clear();
    global_timeline_t = 0;

    Dynamo<kW, kH> effect;
    effect.init();

    // Two overlapping wipes -> two live boundaries (each pushed at the front at
    // angle 0; their Transitions are never stepped here, so they stay put until
    // we overwrite them below).
    effect.color_wipe();
    effect.color_wipe();
    HS_EXPECT_EQ(effect.palette_boundaries.size(), static_cast<size_t>(2));

    // Force the documented worst case: the newer band (index 0) has overtaken
    // the older (index 1) -> non-monotonic order. The chosen magnitudes also push
    // the scan past the first iteration into the second boundary and the
    // baked_palettes_[i+1] access for part of the sweep, exercising the bounds
    // path the safety claim rests on.
    effect.palette_boundaries[0] = 1.5f; // newer, overtaken ahead of the older
    effect.palette_boundaries[1] = 0.5f; // older, left behind
    HS_EXPECT_GT(effect.palette_boundaries[0], effect.palette_boundaries[1]);

    // palette_normal is Z_AXIS (set in the ctor, untouched since the timeline is
    // not stepped), so v = (sin theta, 0, cos theta) sweeps angle_between(v,
    // normal) across the full [0, PI] band span.
    constexpr int kSteps = 256;
    for (int i = 0; i <= kSteps; ++i) {
      float theta = PI_F * static_cast<float>(i) / static_cast<float>(kSteps);
      Vector v(std::sin(theta), 0.0f, std::cos(theta));
      Color4 c = effect.color(v, 0.5f);
      HS_EXPECT_TRUE(std::isfinite(c.alpha));
      HS_EXPECT_GE(c.alpha, 0.0f);
      HS_EXPECT_LE(c.alpha, 1.0f);
    }
  }
};

/**
 * @brief Module entry point for the effects test suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Runs the SH-decode check, then both smoke and determinism passes over
 * the full effect roster at the production and device resolutions.
 */
inline int run_effects_tests() {
  auto scope = hs_test::begin_module("effects");

  test_sh_decode_lm_valid_order();
  test_gs_q16_roundtrip();
  test_gs_rest_state_is_fixed_point();
  test_gs_substep_signs_and_clamp();
  test_gs_evolution_stays_bounded();
  test_dreamballs_preset_cycle_bookkeeping();
  test_dreamballs_respawn_fires_and_honors_pause();
  CometsWhiteBox::check_paths_close();
  ThrustersWhiteBox::check_warp_endpoints();
  RingShowerWhiteBox::check_radius_endpoints();
  DynamoWhiteBox::check_overlapping_wipes_stay_in_range();

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

