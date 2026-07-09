/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Effect smoke / robustness harness — exercises all registered effects.
 *
 * Two passes over the full effect roster:
 *   1. smoke_one  — construct/init/render/read-back, under native asserts.
 *   2. determinism_one — renders each effect twice under an injected, fixed
 *      per-frame clock (hs::set_mock_time) and asserts the two final frames are
 *      byte-identical.
 */
#pragma once

#include "core/engine/effects.h"
#include "core/render/canvas.h"
#include "core/engine/memory.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace hs_test {
namespace effects_tests {

/**
 * @brief Primary production render width in pixels.
 * @details The daydream simulator default and the device's full sphere; the
 * configuration effects are exercised at in deployment, so it is the
 * representative smoke target.
 */
constexpr int DEFAULT_W = 288;
/**
 * @brief Primary production render height in pixels.
 * @details Paired with DEFAULT_W for the full-sphere production resolution.
 */
constexpr int DEFAULT_H = 144;

/**
 * @brief Holosphere device render width in pixels.
 * @details The hardware target instantiates every effect at <96,20>. The native
 * suite is the only place that specialization runs under asserts (the device
 * forces NDEBUG; the CI WASM smoke runs assert-free), so a second roster pass at
 * this resolution exercises height-20-specific paths (PhiLUT<20> indexing,
 * small-aspect arena sizing, H_OFFSET interactions) that bypass every
 * assert-enabled layer otherwise. Both are in HS_WASM_RESOLUTIONS (wasm.cpp).
 */
constexpr int DEVICE_W = 96;
/**
 * @brief Holosphere device render height in pixels.
 * @details Paired with DEVICE_W for the <96,20> device specialization.
 */
constexpr int DEVICE_H = 20;

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
constexpr int DEFAULT_FRAMES = 8;

/**
 * @brief Resolves the per-effect frame count from the environment.
 * @return HS_SMOKE_FRAMES if set to a positive int, else DEFAULT_FRAMES.
 */
inline int smoke_frames() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *e = std::getenv("HS_SMOKE_FRAMES")) {
#pragma clang diagnostic pop
    int n = std::atoi(e);
    if (n > 0) return n;
  }
  return DEFAULT_FRAMES;
}

/**
 * @brief Forward declaration of the registered-but-unread param lint.
 * @param effect Effect instance whose editable params are probed.
 * @param name Effect name used in diagnostic output.
 */
inline void lint_dead_sliders(Effect &effect, const char *name);

// Sweep-wide "something lit up" counter: bumped by smoke_one() on a non-zero
// frame sum, asserted positive once per roster pass in run_effects_tests().
inline int g_nonblack_effects = 0;

// Per-effect non-black exemption. An effect whose ramp-up exceeds the frame
// window legitimately ends on an all-black frame; exempt it only below the frame
// count at which it first lights, so the assertion still fires at the longer
// windows (e.g. CI's HS_SMOKE_FRAMES=120) where output is mandatory.
//   RingShower: rings expand from zero radius; nothing is lit until ~frame 30,
//   so the default 8-frame local window is black by design (it lights by 120).
inline bool effect_may_be_dark(const char *name, int frames) {
  if (std::strcmp(name, "RingShower") == 0)
    return frames < 30;
  return false;
}

/**
 * @brief Drives one effect type through construct -> init -> render -> read-back.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param name Effect name used in the [ok] / diagnostic output.
 * @details Verifies the effect constructs, init's, renders smoke_frames()
 * frames, and reads back every pixel without tripping an assert/OOB/hang, and
 * that get_pixel is a stable pure accessor. Runs the dead-slider lint once on
 * the primary <DEFAULT_W,DEFAULT_H> pass.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void smoke_one(const char *name) {
  // Deterministic RNG so a given effect takes the same code path each run.
  hs::random().seed(1337u);
  // Fresh arena layout for every effect (each effect's init() may reconfigure).
  configure_arenas_default();
  // Timeline state is global/static and shared across effects: without this
  // reset, leftover events reference the previous (destroyed) effect instance.
  Timeline().clear();
  global_timeline_t = 0;

  E<W, H> effect;
  effect.init();

  HS_EXPECT_EQ(effect.width(), W);
  HS_EXPECT_EQ(effect.height(), H);

  const int frames = smoke_frames();
  for (int f = 0; f < frames; ++f) {
    effect.draw_frame();
    // Consume the queued frame, else the next Canvas ctor spin-waits forever.
    effect.advance_display();
  }

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

  // get_pixel is a pure accessor: re-reading the same frame yields the same sum.
  HS_EXPECT_EQ(acc, sum_buffer());

  if (acc > 0)
    ++g_nonblack_effects;
  std::printf("  [ok] %-20s rendered %d frames @ %dx%d (sum=%llu)\n", name,
              frames, W, H, static_cast<unsigned long long>(acc));

  // Per-effect "did it produce output": one effect regressing to all-black no
  // longer hides behind the roster-wide g_nonblack_effects aggregate.
  if (!effect_may_be_dark(name, frames)) {
    if (acc == 0)
      std::printf("  ALL-BLACK %-20s produced no lit pixel over %d frames "
                  "@ %dx%d\n",
                  name, frames, W, H);
    HS_EXPECT(acc > 0, "effect must produce non-black output");
  }

  // The dead-slider lint is resolution-independent; run it once on the primary pass.
  if constexpr (W == DEFAULT_W && H == DEFAULT_H)
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
    // Require the value near `target` AND strictly closer to it than to the
    // pre-write `cur`, catching a slow per-frame revert that 3 frames hide.
    const float eps = fmaxf(1e-3f, 1e-3f * range);
    const float now = def.get();
    const bool persisted =
        fabsf(now - target) <= eps && fabsf(now - target) < fabsf(now - cur);
    if (!persisted)
      std::printf("  DEAD SLIDER %s::%s — wrote %.4f, engine reverted to %.4f "
                  "(markAnimated / markReadonly / drive a private member)\n",
                  name, def.name, static_cast<double>(target),
                  static_cast<double>(now));
    HS_EXPECT(persisted, "editable param must persist across frames");
  }
}

/**
 * @brief Per-frame clock advance in milliseconds for the determinism pass (~30fps).
 * @details Identical across both runs, so any frame-to-frame animation driven by
 * hs::millis/micros or beatsin* is reproduced exactly rather than tracking the
 * wall clock.
 */
constexpr unsigned long FRAME_MS = 33;
/**
 * @brief Per-frame clock advance in microseconds, paired with FRAME_MS.
 */
constexpr unsigned long FRAME_US = 33000;

/**
 * @brief Renders one effect under the injected clock and copies the final buffer out.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param out Receives the final displayed frame, sized W*H pixels (row-major).
 * @param frames Number of frames to render before capture.
 * @details Resets every shared global the smoke path does (RNG seed, arenas,
 * Timeline) plus the generative-hue cursor and mock clock, so two calls start
 * from an identical state.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void render_capture(std::vector<Pixel> &out, int frames) {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;
  // Pin the global generative-hue cursor so both runs start identical (it
  // drifts across palette constructions by design; see GenerativePalette).
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
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
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param name Effect name used in the NONDETERMINISTIC diagnostic output.
 * @details The clock seam neutralizes wall-time, so a divergence here is real
 * nondeterminism (uninitialized read, stale global, address-dependent path) —
 * the defect class smoke coverage cannot see.
 */
/**
 * @brief Scrambles every output-affecting global that render_capture() resets.
 * @details Run between the two determinism captures so the second one must
 * RECOVER canonical output from a dirtied process state rather than merely
 * re-run from an identical pristine one. This makes each reset line in
 * render_capture load-bearing and self-testing: drop the RNG re-seed and an
 * effect that draws randomness diverges; drop the hue-seed reset and a
 * generative-palette effect diverges; drop the global_timeline_t reset and a
 * timeline-driven effect diverges. It also catches a global the effect READS
 * but never mutates that is missing from the reset list — the existing
 * pristine-twice check is blind to that class because both runs read the same
 * unchanged value. It does NOT catch a static seeded once and never reset:
 * in-process that value persists across both runs, so it stays out of reach (as
 * the determinism design notes above).
 */
inline void perturb_determinism_globals() {
  hs::random().seed(0xC0FFEEu);           // off the canonical seed(1337)
  global_timeline_t = 0x5EED;             // off zero
  GenerativePalette::reset_hue_seed(199); // off the zero hue cursor
}

template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void determinism_one(const char *name) {
  const int frames = smoke_frames();
  std::vector<Pixel> a, b;
  render_capture<E, W, H>(a, frames);
  perturb_determinism_globals();
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
  using GS = GSReactionDiffusion<DEVICE_W, DEVICE_H>;
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
    std::vector<float> fA(N), fB(N);
    gs.step_physics(cA, cB, nA, nB, fA.data(), fB.data());
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
  // rail; an unstable oscillation would clamp a large fraction there.
  int saturated = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (cB[i] == 65535)
      ++saturated;
  HS_EXPECT_LT(saturated, GSWhiteBox::N / 20);
}

// ---------------------------------------------------------------------------
// Belousov-Zhabotinsky reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for BZReactionDiffusion's private fixed-point and
 *        physics internals (befriended in effects/BZReactionDiffusion.h).
 * @details The matching seam to GSWhiteBox: the smoke/determinism harness cannot
 *          see a Q8 round-trip error, a sign/clamp slip in the Lotka-Volterra
 *          update, or a perturbation that wraps past the Q8 rail, so the
 *          conversions, one species step, the perturbation, and one fused physics
 *          substep are pinned directly. The lattice is the same fixed 7680-node
 *          graph as GS, independent of <W,H>, so the device resolution is used
 *          arbitrarily.
 */
struct BZWhiteBox {
  using BZ = BZReactionDiffusion<DEVICE_W, DEVICE_H>;
  static constexpr int N = BZ::RD_N;

  static uint8_t to_q8(float v) { return BZ::to_q8(v); }
  static float from_q8(uint8_t v) { return BZ::from_q8(v); }
  static void set_params(BZ &bz, float alpha, float D, float dt) {
    bz.params.alpha = alpha;
    bz.params.D = D;
    bz.params.dt = dt;
  }
  static uint8_t advance_species(const BZ &bz, float conc, float predator,
                                 float laplacian) {
    return bz.advance_species(conc, predator, laplacian);
  }
  static void perturb(uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    BZ::perturb_state(nA, nB, nC);
  }
  static int num_perturbations() { return BZ::NUM_PERTURBATIONS; }
  static void step(BZ &bz, const uint8_t *cA, const uint8_t *cB,
                   const uint8_t *cC, uint8_t *nA, uint8_t *nB, uint8_t *nC) {
    std::vector<float> fA(N), fB(N), fC(N);
    bz.step_physics(cA, cB, cC, nA, nB, nC, fA.data(), fB.data(), fC.data());
  }
  static void set_state(BZ &bz, uint8_t *A, uint8_t *B, uint8_t *C) {
    bz.state.A = A;
    bz.state.B = B;
    bz.state.C = C;
  }
  static void advance_substeps(BZ &bz, int steps, uint8_t *sA, uint8_t *sB,
                               uint8_t *sC) {
    std::vector<float> fA(N), fB(N), fC(N);
    bz.advance_substeps(
        steps,
        std::array<uint8_t *, 3>{bz.state.A, bz.state.B, bz.state.C},
        std::array<uint8_t *, 3>{sA, sB, sC}, [&](auto &cur, auto &nxt) {
          bz.step_physics(cur[0], cur[1], cur[2], nxt[0], nxt[1], nxt[2],
                          fA.data(), fB.data(), fC.data());
        });
  }
};

/**
 * @brief Verifies the Q8 fixed-point round-trip and the +0.5 rounding/clamp
 *        boundaries.
 * @details to_q8(from_q8(v)) must be the identity over every byte value, and
 *          to_q8 must clamp out-of-range floats and round to nearest (so 1.0 tops
 *          out at 255 with no overflow). A truncating or unclamped regression —
 *          which would bias the RD dynamics downward — fails here.
 */
inline void test_bz_q8_roundtrip() {
  HS_EXPECT_EQ(BZWhiteBox::to_q8(0.0f), (uint8_t)0);
  HS_EXPECT_EQ(BZWhiteBox::to_q8(1.0f), (uint8_t)255);
  HS_EXPECT_EQ(BZWhiteBox::to_q8(2.0f), (uint8_t)255);   // clamp high
  HS_EXPECT_EQ(BZWhiteBox::to_q8(-0.5f), (uint8_t)0);    // clamp low
  HS_EXPECT_NEAR(BZWhiteBox::from_q8(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(BZWhiteBox::from_q8(255), 1.0f, 1e-9f);
  int bad = 0;
  for (int v = 0; v <= 255; ++v)
    if (BZWhiteBox::to_q8(BZWhiteBox::from_q8((uint8_t)v)) != (uint8_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies advance_species has the right reaction/diffusion signs and
 *        that the Q8 clamp backstop holds even past the Euler stability bound.
 * @details advance_species is the single-species core of the BZ update:
 *          conc + (D·laplacian + conc·(1 − conc − α·predator))·dt, mapped through
 *          to_q8. Checks: an empty rest cell stays empty (no spurious growth);
 *          diffusion from higher neighbors grows an empty cell; logistic growth
 *          lifts a half-filled, predator-free cell; predation underflow clamps to
 *          0 (not a uint8 wrap to 255); and extreme over-/under-shoots — the
 *          documented to_q8 backstop the comment in the effect promises — clamp
 *          to the [0, 255] rails rather than wrapping.
 */
inline void test_bz_advance_species_signs_and_clamp() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 0.35f);

  // Empty rest cell: no diffusion, no reaction -> stays 0.
  HS_EXPECT_EQ((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, 0.0f), 0);

  // Diffusion from higher neighbors lifts an empty cell above 0.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, /*lap*/ 6.0f), 0);

  // Logistic growth: a half-filled, predator-free cell grows above its start.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.5f, 0.0f, 0.0f),
               (int)BZWhiteBox::to_q8(0.5f));

  // Predation drives a saturated cell negative; it must clamp to 0, not wrap.
  HS_EXPECT_EQ((int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1.0f, 0.0f),
               0);

  // Backstop past the Euler bound: a huge positive laplacian clamps to 255, a
  // huge predation clamps to 0 — to_q8 keeps every written state in range.
  HS_EXPECT_EQ((int)BZWhiteBox::advance_species(bz, 1.0f, 0.0f, /*lap*/ 1000.0f),
               255);
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1000.0f, 0.0f), 0);
}

/**
 * @brief Verifies perturb_state nudges nodes by a fixed Q8 amount and saturates
 *        at the 255 rail without wrapping.
 * @details Two RNG-agnostic invariants. (1) A fully-saturated field stays fully
 *          saturated: every nudge is a +PERTURB_AMOUNT saturating add, so a node
 *          already at 255 cannot wrap to a low value. (2) On a zero field, each
 *          touched entry is a small positive multiple of the nudge step and stays
 *          within [0, 255], and at least one entry is touched — whichever nodes
 *          the deterministic RNG happens to select.
 */
inline void test_bz_perturb_state_saturates_and_nudges() {
  // (1) Saturation / no-wrap: all rails stay at the rail.
  {
    std::vector<uint8_t> a(BZWhiteBox::N, 255), b(BZWhiteBox::N, 255),
        c(BZWhiteBox::N, 255);
    BZWhiteBox::perturb(a.data(), b.data(), c.data());
    int wrapped = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      if (a[i] != 255 || b[i] != 255 || c[i] != 255)
        ++wrapped;
    HS_EXPECT_EQ(wrapped, 0);
  }
  // (2) Zero field: touched entries are small positive multiples of the step.
  {
    std::vector<uint8_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    BZWhiteBox::perturb(a.data(), b.data(), c.data());
    int touched = 0, malformed = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      for (uint8_t v : {a[i], b[i], c[i]}) {
        if (v == 0)
          continue;
        ++touched;
        if (v % 3 != 0) // PERTURB_AMOUNT == 3; accumulations stay multiples of 3
          ++malformed;
      }
    HS_EXPECT_GT(touched, 0);
    HS_EXPECT_EQ(malformed, 0);
  }
}

/**
 * @brief Pins perturb_state's per-frame draw count on the shared RNG stream.
 * @details perturb_state advances hs::random() by exactly 2*NUM_PERTURBATIONS
 *          draws (idx + species per nudge). That count is part of the global
 *          determinism contract: every later effect's stream position depends on
 *          it, so a substep/perturbation retune that changes the draw count
 *          silently shifts all downstream effects. Reproduce the post-call stream
 *          position with a private generator stepped exactly that many times and
 *          require the global generator to be at the same position (next outputs
 *          equal), failing if the count drifts in either direction.
 */
inline void test_bz_perturb_state_draw_count_pinned() {
  const int expected_draws = 2 * BZWhiteBox::num_perturbations();

  constexpr uint64_t SEED = 1337u;
  hs::random().seed(SEED);
  std::vector<uint8_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
      c(BZWhiteBox::N, 0);
  BZWhiteBox::perturb(a.data(), b.data(), c.data());

  // A private generator from the same seed, advanced by the contracted count,
  // must now be at the same stream position as the global generator.
  hs::Pcg32 ref(SEED);
  for (int i = 0; i < expected_draws; ++i)
    (void)ref();
  HS_EXPECT_EQ(hs::random()(), ref());
  // Off-by-one in either direction would have landed at a different output.
  HS_EXPECT_EQ(hs::random()(), ref());
}

/**
 * @brief Verifies one fused physics substep diffuses a seeded species into its
 *        neighborhood with the right sign.
 * @details Seed species A fully at one interior node on an otherwise-empty field
 *          and run a single step. The seed's neighbors start empty but border a
 *          saturated node, so their graph-Laplacian is positive and A must
 *          diffuse into at least one of them; the seed node itself stays lit
 *          (it decays but does not vanish or wrap in one step). The stochastic
 *          perturbation runs too, but the seed-neighbor diffusion is independent
 *          of which nodes it nudges.
 */
inline void test_bz_substep_diffuses() {
  std::vector<uint8_t> cA(BZWhiteBox::N, 0), cB(BZWhiteBox::N, 0),
      cC(BZWhiteBox::N, 0);
  std::vector<uint8_t> nA(BZWhiteBox::N, 0), nB(BZWhiteBox::N, 0),
      nC(BZWhiteBox::N, 0);
  const int seed = 4000; // interior lattice node with a full neighbor ring
  cA[seed] = 255;

  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, 3.0f, 0.05f, 0.35f);
  BZWhiteBox::step(bz, cA.data(), cB.data(), cC.data(), nA.data(), nB.data(),
                   nC.data());

  HS_EXPECT_GT((int)nA[seed], 0); // the seed decays but does not vanish/wrap
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && nA[nb] > 0)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0); // A diffused into at least one empty neighbor
}

/**
 * @brief Verifies an odd substep count lands the final generation back in the
 *        persistent state buffers.
 * @details render() ping-pongs between the persistent state and scratch, then
 *          copies back only when the final generation ended up in scratch — the
 *          parity copy-back, which the even production STEPS_PER_FRAME never
 *          exercises. Drive a single (odd) substep on a seeded field: the
 *          diffused result must appear in the persistent A buffer (not just the
 *          scratch one), so a missing copy-back leaves the neighbors empty and
 *          fails here.
 */
inline void test_bz_odd_substep_lands_in_state() {
  std::vector<uint8_t> A(BZWhiteBox::N, 0), B(BZWhiteBox::N, 0),
      C(BZWhiteBox::N, 0);
  std::vector<uint8_t> sA(BZWhiteBox::N, 0), sB(BZWhiteBox::N, 0),
      sC(BZWhiteBox::N, 0);
  const int seed = 4000; // interior lattice node with a full neighbor ring
  A[seed] = 255;

  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, 3.0f, 0.05f, 0.35f);
  BZWhiteBox::set_state(bz, A.data(), B.data(), C.data());

  // One (odd) substep: the result must be copied back into A, not stranded in sA.
  BZWhiteBox::advance_substeps(bz, 1, sA.data(), sB.data(), sC.data());

  HS_EXPECT_GT((int)A[seed], 0); // seed decayed in place, landed back in A
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && A[nb] > 0)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0); // diffusion landed in the persistent buffer
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
  using DB = DreamBalls<DEVICE_W, DEVICE_H>;
  static constexpr int PRESETS = 4;

  static int active_bake(const DB &db) { return db.active_bake_; }
  static int last_preset_idx(const DB &db) { return db.last_preset_idx_; }
  static void spawn(DB &db, int idx) { db.spawn_sprite(idx); }
  // The literal solid-name pointer the given preset seeds into params on reseed.
  static const char *preset_name(const DB &db, int idx) {
    return db.preset_manager.get_entries()[idx].params.solid_name;
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
    const int safe = idx % WB::PRESETS;
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
    using C = Comets<DEFAULT_W, DEFAULT_H>;
    int idx = 0;
    for (const LissajousParams &cfg : C::functions) {
      const float cd = C::closing_domain(cfg);
      HS_EXPECT_GT(cd, 0.0f); // floor-at-1 keeps the head moving
      // Every authored entry must clear the floor (m2*domain >= PI rounds to >= 1
      // closing cycle) so the floor never silently rewrites an authored domain.
      HS_EXPECT_GE(cfg.m2 * cfg.domain, PI_F);
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
    using T = Thrusters<DEFAULT_W, DEFAULT_H>;
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
    using RS = RingShower<DEFAULT_W, DEFAULT_H>;
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
    // Dynamo::init() bakes from persistent_arena and schedules on the shared
    // global timeline, so reset the shared globals as smoke_one does.
    hs::random().seed(1337u);
    configure_arenas_default();
    Timeline().clear();
    global_timeline_t = 0;

    Dynamo<DEFAULT_W, DEFAULT_H> effect;
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
    constexpr int STEPS = 256;
    for (int i = 0; i <= STEPS; ++i) {
      float theta = PI_F * static_cast<float>(i) / static_cast<float>(STEPS);
      Vector v(std::sin(theta), 0.0f, std::cos(theta));
      Color4 c = effect.color(v, 0.5f);
      HS_EXPECT_TRUE(std::isfinite(c.alpha));
      HS_EXPECT_GE(c.alpha, 0.0f);
      HS_EXPECT_LE(c.alpha, 1.0f);
    }
  }
};

// ---------------------------------------------------------------------------
// Drift-prone effects: spawn-gap / emit-phase / pool-bound white-box pins.
//
// The smoke pass proves these render and the determinism pass proves they
// reproduce, but neither sees a per-frame accumulator that runs away, a phase
// that escapes its wrap interval, or a pool index that overruns its capacity —
// silent drift that still renders a plausible (wrong) frame. These seams reach
// the private state and assert the bound holds on every frame.
// ---------------------------------------------------------------------------

/**
 * @brief Shared per-effect setup mirroring smoke_one's global reset.
 */
inline void reset_effect_globals() {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;
}

/**
 * @brief White-box accessor for PetalFlow's spawn-gap accumulator and hue cursor
 *        (befriended in effects/PetalFlow.h).
 */
struct PetalFlowWhiteBox {
  using PF = PetalFlow<DEFAULT_W, DEFAULT_H>;
  static float gap(const PF &pf) { return pf.gap_accumulator; }
  static float spacing() { return PF::SPACING; }
  static float next_hue(const PF &pf) { return pf.next_hue; }
};

/**
 * @brief Verifies the spawn-gap accumulator drains every frame and the hue
 *        cursor stays wrapped.
 * @details check_spawn() integrates Speed*RHO_PER_SPEED into gap_accumulator and
 *          drains it by SPACING per spawn, so after any frame the residue must
 *          satisfy 0 <= gap < SPACING — a runaway (missing drain) or a negative
 *          residue both fail here. Run at the Speed slider top so the per-frame
 *          travel exceeds SPACING and the while-loop must emit several rings per
 *          frame. next_hue is advanced wrap(.,1) per spawn and must stay [0, 1).
 */
inline void test_petalflow_spawn_gap_bounded() {
  using WB = PetalFlowWhiteBox;
  reset_effect_globals();
  WB::PF pf;
  pf.init();
  pf.updateParameter("Speed", 20.0f); // slider top: per-frame travel >> SPACING

  const float spacing = WB::spacing();
  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    pf.draw_frame();
    pf.advance_display();
    const float gap = WB::gap(pf);
    HS_EXPECT_GE(gap, 0.0f);
    HS_EXPECT_LT(gap, spacing); // accumulator never runs away past one spacing
    const float hue = WB::next_hue(pf);
    HS_EXPECT_GE(hue, 0.0f);
    HS_EXPECT_LT(hue, 1.0f); // hue cursor stays wrapped
  }
}

/**
 * @brief White-box accessor for MindSplatter's per-emitter emit-phase and hue
 *        arrays (befriended in effects/MindSplatter.h).
 */
struct MindSplatterWhiteBox {
  using MS = MindSplatter<DEFAULT_W, DEFAULT_H>;
  static size_t num_emitters() { return MS::EmitSolid::NUM_VERTS; }
  static float emit_phase(const MS &ms, size_t i) { return ms.emit_phases[i]; }
  static float hue(const MS &ms, size_t i) { return ms.emitter_hues[i]; }
};

/**
 * @brief Verifies every per-emitter emission phase stays wrapped to [0, 2pi) and
 *        each hue stays [0, 1) across frames at the max angular rate.
 * @details Each emitter integrates Ang Spd into emit_phases[i] with fmodf(., 2pi)
 *          and advances emitter_hues[i] with fmodf(., 1); a dropped wrap lets the
 *          phase grow unbounded (fast_sinf range reduction then bands). Run at the
 *          Ang Spd slider top so the phase laps 2pi repeatedly.
 */
inline void test_mindsplatter_emit_phase_wrapped() {
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  WB::MS ms;
  ms.init();
  ms.updateParameter("Ang Spd", 1.0f); // slider top: phase laps fast

  const float two_pi = 2.0f * PI_F;
  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    ms.draw_frame();
    ms.advance_display();
    for (size_t i = 0; i < WB::num_emitters(); ++i) {
      const float ph = WB::emit_phase(ms, i);
      HS_EXPECT_GE(ph, 0.0f);
      HS_EXPECT_LT(ph, two_pi);
      const float h = WB::hue(ms, i);
      HS_EXPECT_GE(h, 0.0f);
      HS_EXPECT_LT(h, 1.0f);
    }
  }
}

/**
 * @brief White-box accessor for Flyby's noise-time and trig-phase accumulators
 *        (befriended in effects/Flyby.h).
 */
struct FlybyWhiteBox {
  using FB = Flyby<DEFAULT_W, DEFAULT_H>;
  static float noise_time(const FB &fb) { return fb.noise_time; }
  static float time_period() { return FB::TIME_PERIOD; }
  static float sin_phase(const FB &fb) { return fb.sin_phase; }
  static float drift_phase(const FB &fb) { return fb.drift_phase; }
};

/**
 * @brief Verifies Flyby's noise-time stays in [0, TIME_PERIOD) and both trig
 *        phases stay in [0, 2pi) across frames.
 * @details draw_frame wraps noise_time by fmodf(., TIME_PERIOD) and sin/drift
 *          phase by fmodf(., 2pi); a dropped wrap freezes the field (ULP) or
 *          bands fast_sinf range reduction. The preset Lerp drives Speed (and
 *          Drift scales the cos phase), so the accumulators must stay bounded.
 */
inline void test_flyby_phase_wrapped() {
  using WB = FlybyWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::FB fb;
  fb.init();

  const float period = WB::time_period();
  const float two_pi = 2.0f * PI_F;
  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fb.draw_frame();
    fb.advance_display();
    const float nt = WB::noise_time(fb);
    HS_EXPECT_GE(nt, 0.0f);
    HS_EXPECT_LT(nt, period);
    const float sp = WB::sin_phase(fb);
    HS_EXPECT_GE(sp, 0.0f);
    HS_EXPECT_LT(sp, two_pi);
    const float dp = WB::drift_phase(fb);
    HS_EXPECT_GE(dp, 0.0f);
    HS_EXPECT_LT(dp, two_pi);
  }
  hs::clear_mock_time();
}

/**
 * @brief White-box accessor for Liquid2D's noise-time and trig-phase
 *        accumulators (befriended in effects/Liquid2D.h).
 */
struct Liquid2DWhiteBox {
  using L2 = Liquid2D<DEFAULT_W, DEFAULT_H>;
  static float time_period() { return L2::TIME_PERIOD; }
  static float accumulated_time(const L2 &l2) { return l2.accumulated_time; }
  static float sin_phase(const L2 &l2) { return l2.sin_phase; }
  static float cos_phase(const L2 &l2) { return l2.cos_phase; }
  static float cycle_phase(const L2 &l2) { return l2.cycle_phase; }
  static void seed_accumulators(L2 &l2, float v) {
    l2.accumulated_time = v;
    l2.sin_phase = v;
    l2.cos_phase = v;
    l2.cycle_phase = v;
  }
};

/**
 * @brief Verifies Liquid2D's noise-time stays in [0, TIME_PERIOD) and all three
 *        trig phases stay in [0, 2pi) once wrapped.
 * @details draw_frame wraps accumulated_time by fmodf(., TIME_PERIOD) and the
 *          sin/cos/cycle phases by fmodf(., 2pi); a dropped wrap freezes the
 *          field (ULP) or bands fast_sinf range reduction. Seeding every
 *          accumulator past its ceiling makes a removed fmodf fail on frame 0.
 */
inline void test_liquid2d_phase_wrapped() {
  using WB = Liquid2DWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::L2 l2;
  l2.init();

  const float period = WB::time_period();
  const float two_pi = 2.0f * PI_F;
  WB::seed_accumulators(l2, period * 4.0f);

  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    l2.draw_frame();
    l2.advance_display();
    const float at = WB::accumulated_time(l2);
    HS_EXPECT_GE(at, 0.0f);
    HS_EXPECT_LT(at, period);
    const float sp = WB::sin_phase(l2);
    HS_EXPECT_GE(sp, 0.0f);
    HS_EXPECT_LT(sp, two_pi);
    const float cp = WB::cos_phase(l2);
    HS_EXPECT_GE(cp, 0.0f);
    HS_EXPECT_LT(cp, two_pi);
    const float cyp = WB::cycle_phase(l2);
    HS_EXPECT_GE(cyp, 0.0f);
    HS_EXPECT_LT(cyp, two_pi);
  }
  hs::clear_mock_time();
}

/**
 * @brief White-box accessor for FlowField's noise-time and particle pool
 *        (befriended in effects/FlowField.h).
 */
struct FlowFieldWhiteBox {
  using FF = FlowField<DEFAULT_W, DEFAULT_H>;
  static float noise_time(const FF &ff) { return ff.t; }
  static float time_period() { return FF::TIME_PERIOD; }
  static uint16_t active_count(const FF &ff) {
    return ff.particle_system.active_count;
  }
  static size_t pool_capacity(const FF &ff) {
    return ff.particle_system.pool.capacity();
  }
};

/**
 * @brief Verifies FlowField's noise-time stays in [0, TIME_PERIOD) and the
 *        particle pool never exceeds capacity.
 * @details draw_frame wraps t by fmodf(., TIME_PERIOD); the emitter back-fills
 *          the pool to capacity every frame, so active_count must hold at exactly
 *          capacity without overrun. Run at the Time Spd top so t advances fast.
 */
inline void test_flowfield_time_and_pool_bounded() {
  using WB = FlowFieldWhiteBox;
  reset_effect_globals();
  WB::FF ff;
  ff.init();
  ff.updateParameter("Time Spd", 0.05f); // slider top: t advances fast

  const float period = WB::time_period();
  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    ff.draw_frame();
    ff.advance_display();
    const float nt = WB::noise_time(ff);
    HS_EXPECT_GE(nt, 0.0f);
    HS_EXPECT_LT(nt, period);
    HS_EXPECT_LE(static_cast<size_t>(WB::active_count(ff)),
                 WB::pool_capacity(ff));
  }
}

/**
 * @brief White-box accessor for RingSpin's live ring count (befriended in
 *        effects/RingSpin.h).
 */
struct RingSpinWhiteBox {
  using RS = RingSpin<DEFAULT_W, DEFAULT_H>;
  static int num_rings(const RS &rs) { return rs.num_rings; }
  static int max_rings() { return RS::NUM_RINGS; }
  static void spawn(RS &rs) { rs.spawn_ring(&rs.baked_palettes[0]); }
};

/**
 * @brief Verifies spawn_ring clamps to the fixed pool and never overruns it.
 * @details init() fills the pool to NUM_RINGS; the guard (num_rings >= NUM_RINGS)
 *          must make every further spawn a no-op, holding num_rings pinned at the
 *          bound — an off-by-one would write past the rings[] allocation. Drive
 *          extra spawns directly past the bound to exercise the guard branch.
 */
inline void test_ringspin_pool_clamped() {
  using WB = RingSpinWhiteBox;
  reset_effect_globals();
  WB::RS rs;
  rs.init();
  HS_EXPECT_EQ(WB::num_rings(rs), WB::max_rings()); // init filled the pool
  for (int i = 0; i < WB::max_rings() + 4; ++i) {
    WB::spawn(rs); // each must be rejected by the >= NUM_RINGS guard
    HS_EXPECT_EQ(WB::num_rings(rs), WB::max_rings());
  }
}

/**
 * @brief Drives ShapeShifter through its 48-frame shape cut to cover the cycle
 *        wrap the default 8-frame smoke window never reaches.
 * @details draw_frame advances frame_count_ mod 48 and rotates current_shape on
 *          the wrap, so the cut path (and the post-cut shape's renderer) only
 *          executes past frame 48. Run two full periods under a fixed clock and
 *          require non-black output on every frame: the cut must neither blank
 *          the effect nor blow an assert on the new shape's first render. Uses
 *          only public draw_frame/get_pixel — no production seam.
 */
inline void test_shapeshifter_shape_cut_lifecycle() {
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  ShapeShifter<DEFAULT_W, DEFAULT_H> ss;
  ss.init();

  const int period = 48;
  const int frames = smoke_frames() < 2 * period ? 2 * period : smoke_frames();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    ss.draw_frame();
    ss.advance_display();
    uint64_t acc = 0;
    for (int y = 0; y < DEFAULT_H; ++y)
      for (int x = 0; x < DEFAULT_W; ++x) {
        const Pixel &p = ss.get_pixel(x, y);
        acc += static_cast<uint64_t>(p.r) + p.g + p.b;
      }
    HS_EXPECT_GT(acc, 0u); // each shape, including post-cut, must render
  }
  hs::clear_mock_time();
}

/**
 * @brief Verifies the segmented-mode full-frame gate on the real effect roster.
 * @details Effect::needs_full_frame() must report true for exactly the effects
 *          whose filter pipeline crosses segment bands (its trait-derived
 *          override surfaces the compile-time any_crosses_segments fold), and
 *          false otherwise so non-stateful effects keep the segment clipping
 *          win. Cross-segment today: MeshFeedback (Pixel::Feedback) and the two
 *          World::Trails effects (Dynamo, SplineFlow). This pins the gate end to
 *          end on constructed effects — the WASM driver reads exactly this query
 *          (targets/wasm/wasm.cpp setClip). See
 *          docs/segmented_stateful_effects_spec.md.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer (single-live guard) and
  // reconfigures the shared arenas/timeline in init(), so construct one at a
  // time with the same fresh setup smoke_one uses.
  auto reset = [] {
    hs::random().seed(1337u);
    configure_arenas_default();
    Timeline().clear();
    global_timeline_t = 0;
  };
  auto check = [&](Effect &fx, bool expected, const char *name) {
    fx.init();
    if (fx.needs_full_frame() != expected)
      std::printf("  [FAIL] %s needs_full_frame()=%d expected=%d\n", name,
                  fx.needs_full_frame(), expected);
    HS_EXPECT_EQ(fx.needs_full_frame(), expected);
  };

  // Cross-segment stateful effects -> full-frame.
  { reset(); MeshFeedback<DEFAULT_W, DEFAULT_H> fx; check(fx, true, "MeshFeedback"); }
  { reset(); Dynamo<DEFAULT_W, DEFAULT_H> fx;       check(fx, true, "Dynamo"); }
  { reset(); SplineFlow<DEFAULT_W, DEFAULT_H> fx;   check(fx, true, "SplineFlow"); }

  // Representative non-stateful effects -> keep band clipping (default false).
  { reset(); Voronoi<DEFAULT_W, DEFAULT_H> fx;  check(fx, false, "Voronoi"); }
  { reset(); RingSpin<DEFAULT_W, DEFAULT_H> fx; check(fx, false, "RingSpin"); }
}

/**
 * @brief Pins Voronoi's coarse-coherence equivalence: where a block's four
 *        corners share the canonical nearest-pair, every interior pixel resolves
 *        to that same pair under a full k=2 query.
 * @details Mirrors Voronoi::draw_frame's corner sampling (COHERENCE_BLOCK stride,
 *          clamped corners, order-independent {lo,hi} pair) over the public
 *          KDTree at the production resolution, using six maximally separated
 *          (octahedral) sites so no cell is smaller than a block — the regime in
 *          which the in-effect "bit-identical to the full query" fast path is
 *          exact. A sub-block cell dropped by all four corners (the documented
 *          high-density limitation) would break this invariant.
 */
inline void test_voronoi_coherence_equivalence() {
  constexpr int W = DEFAULT_W, H = DEFAULT_H;
  constexpr int B = 8; // Voronoi<W, H>::COHERENCE_BLOCK

  static uint8_t buf[64 * 1024];
  Arena arena(buf, sizeof(buf));

  const Vector sites[] = {
      Vector(1, 0, 0),  Vector(-1, 0, 0), Vector(0, 1, 0),
      Vector(0, -1, 0), Vector(0, 0, 1),  Vector(0, 0, -1),
  };
  constexpr int N = sizeof(sites) / sizeof(sites[0]);
  KDTree tree(arena, std::span<const Vector>(sites, N));

  struct Pair {
    uint16_t lo, hi;
    bool hasSecond;
  };
  auto classify = [&](const Vector &p) -> Pair {
    auto knn = tree.nearest(p, 2);
    uint16_t a = knn[0].original_index;
    bool hasSecond = knn.size() > 1;
    uint16_t b = hasSecond ? knn[1].original_index : a;
    return {std::min(a, b), std::max(a, b), hasSecond};
  };
  auto same = [](const Pair &x, const Pair &y) {
    return x.lo == y.lo && x.hi == y.hi && x.hasSecond == y.hasSecond;
  };

  constexpr int NBX = (W - 1) / B + 2;
  constexpr int NBY = (H - 1) / B + 2;
  static Pair corners[NBX * NBY];
  auto cx = [&](int j) { return std::min(j * B, W - 1); };
  auto cy = [&](int k) { return std::min(k * B, H - 1); };
  for (int k = 0; k < NBY; ++k)
    for (int j = 0; j < NBX; ++j)
      corners[k * NBX + j] = classify(pixel_to_vector<W, H>(cx(j), cy(k)));

  int fast_pixels = 0;
  for (int y = 0; y < H; ++y) {
    const int ky = y / B;
    for (int x = 0; x < W; ++x) {
      const int jx = x / B;
      const Pair &c00 = corners[ky * NBX + jx];
      const Pair &c10 = corners[ky * NBX + (jx + 1)];
      const Pair &c01 = corners[(ky + 1) * NBX + jx];
      const Pair &c11 = corners[(ky + 1) * NBX + (jx + 1)];
      // Seam/dense block -> the effect takes the full query; nothing to pin.
      if (!(same(c00, c10) && same(c00, c01) && same(c00, c11)))
        continue;
      HS_EXPECT_TRUE(same(classify(pixel_to_vector<W, H>(x, y)), c00));
      ++fast_pixels;
    }
  }
  // The fast path must actually be exercised, else the test pins nothing.
  HS_EXPECT_GT(fast_pixels, 0);
}

/**
 * @brief Minimal Effect owning a Canvas, so a mesh can be rasterized in a test.
 */
struct BudgetCanvasFx : public Effect {
  /**
   * @brief Constructs the fixture with a canvas of the given size.
   * @param w Canvas width in pixels.
   * @param h Canvas height in pixels.
   */
  BudgetCanvasFx(int w, int h) : Effect(w, h) {}
  /** @brief No-op per-frame hook; the mesh draw under test lights the canvas. */
  void draw_frame() override {}
};

/**
 * @brief Gates HankinSolids' hand-tuned scratch budgets against the real
 *        generate+classify+render+compaction peak of every simple solid at the
 *        device height (H=144).
 * @details HankinSolids::init sizes scratch_a=24 KB / scratch_b=32 KB for "the
 *          heaviest hankin mesh", but the effect's random morph never
 *          deterministically reaches every solid within a smoke window, so an
 *          under-sized budget would trap only on hardware. This reproduces the
 *          effect's per-solid arena discipline (load_shape → classify → draw_mesh
 *          → morph compaction) for all Platonic+Archimedean solids in
 *          headroomed arenas and asserts each scratch high-water fits its budget.
 */
inline void test_hankinsolids_arena_budget_covers_every_solid() {
  constexpr int W = 288, H = 144;
  constexpr size_t SCRATCH_A = 24 * 1024, SCRATCH_B = 32 * 1024;
  constexpr size_t MEASURE = 1024 * 1024; // headroom so a peak never traps here
  constexpr float ANGLE = PI_F / 4.0f;

  auto solids = Solids::Collections::get_simple_solids();
  for (size_t idx = 0; idx < solids.size(); ++idx) {
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * MEASURE, MEASURE, MEASURE);

    MeshPaletteBank palette_bank;
    palette_bank.bake_all(persistent_arena);

    MeshState mesh;
    CompiledHankin hankin;
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = Solids::finalize_solid(solids[idx].generate(a, b), a);
      hankin = CompiledHankin();
      MeshOps::compile_hankin(base, hankin, target, a);
      mesh.clear();
      MeshOps::update_hankin(hankin, mesh, target, ANGLE);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(mesh, scratch_arena_a, scratch_arena_b,
                                          persistent_arena);
    }

    // Render peak: transform into scratch_a, then Scan::Mesh::draw stacks a
    // FaceScratchBuffer on top (the scratch_a-binding path per init's comment).
    {
      ScratchScope a_guard(scratch_arena_a);
      Orientation<> orientation;
      OrientTransformer camera(orientation);
      MeshState rotated;
      MeshOps::transform(mesh, rotated, scratch_arena_a, camera);
      BudgetCanvasFx fx(W, H);
      Canvas canvas(fx);
      Pipeline<W, H> filters;
      auto frag = [](const Vector &, Fragment &f) {
        f.color = Color4(Pixel(1000, 1000, 1000), 1.0f);
      };
      Scan::Mesh::draw<W, H>(filters, canvas, rotated, frag, scratch_arena_a);
    }

    // Morph compaction peak: the CompiledHankin + palette bank survive into
    // scratch_b, the mesh into scratch_a, then persistent is reset — the same
    // Persist discipline start_morph_cycle uses to compact between solids.
    {
      Persist<CompiledHankin> ph(hankin, scratch_arena_b, persistent_arena);
      Persist<MeshState> pf(mesh, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank, scratch_arena_b, persistent_arena);
      persistent_arena.reset();
    }

    const size_t a_peak = scratch_arena_a.get_high_water_mark();
    const size_t b_peak = scratch_arena_b.get_high_water_mark();
    if (a_peak > SCRATCH_A || b_peak > SCRATCH_B)
      std::printf("  HankinSolids arena OVER BUDGET solid[%zu] '%s': "
                  "scratchA=%zu/%zu scratchB=%zu/%zu\n",
                  idx, solids[idx].name, a_peak, SCRATCH_A, b_peak, SCRATCH_B);
    HS_EXPECT_TRUE(a_peak <= SCRATCH_A);
    HS_EXPECT_TRUE(b_peak <= SCRATCH_B);
  }
}

/**
 * @brief Module entry point for the effects test suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Runs the SH-decode check, then both smoke and determinism passes over
 * the full effect roster at the production and device resolutions.
 */
inline int run_effects_tests() {
  hs_test::ModuleFixture fixture("effects");

  // SSOT anti-drift guard, mirroring the WASM startup check (targets/wasm/wasm.cpp):
  // the self-registering effect count (each header's REGISTER_EFFECT) must equal
  // the static HS_EFFECT_LIST roster, or an effect present in one and missing from
  // the other silently drops smoke coverage below. Active because the test build
  // defines HS_TEST_BUILD (see core/engine/effect_registry.h).
  HS_EXPECT_EQ(EffectRegistry::entries().size(),
               static_cast<size_t>(HS_EFFECT_COUNT));

  test_needs_full_frame_gate();
  test_voronoi_coherence_equivalence();
  test_sh_decode_lm_valid_order();
  test_gs_q16_roundtrip();
  test_gs_rest_state_is_fixed_point();
  test_gs_substep_signs_and_clamp();
  test_gs_evolution_stays_bounded();
  test_bz_q8_roundtrip();
  test_bz_advance_species_signs_and_clamp();
  test_bz_perturb_state_saturates_and_nudges();
  test_bz_perturb_state_draw_count_pinned();
  test_bz_substep_diffuses();
  test_bz_odd_substep_lands_in_state();
  test_dreamballs_preset_cycle_bookkeeping();
  test_dreamballs_respawn_fires_and_honors_pause();
  CometsWhiteBox::check_paths_close();
  ThrustersWhiteBox::check_warp_endpoints();
  RingShowerWhiteBox::check_radius_endpoints();
  DynamoWhiteBox::check_overlapping_wipes_stay_in_range();
  test_petalflow_spawn_gap_bounded();
  test_mindsplatter_emit_phase_wrapped();
  test_flyby_phase_wrapped();
  test_liquid2d_phase_wrapped();
  test_flowfield_time_and_pool_bounded();
  test_ringspin_pool_clamped();
  test_shapeshifter_shape_cut_lifecycle();
  test_hankinsolids_arena_budget_covers_every_solid();

  // Smoke every registered effect. The list is generated from the single-source
  // roster in core/engine/effects.h (HS_EFFECT_LIST), so it cannot drift from the
  // shipped set.
  g_nonblack_effects = 0;
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE
  // At least one effect must light up, catching a total regression-to-black.
  HS_EXPECT_GT(g_nonblack_effects, 0);

  // Second pass: cross-run determinism under the injected clock.
#define HS_DET_ONE(name) determinism_one<name>(#name);
  HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE

  // Repeat both passes at the Holosphere device resolution <96,20>, the only
  // place that specialization runs under native asserts (see DEVICE_W/DEVICE_H).
  std::printf("  -- device resolution %dx%d --\n", DEVICE_W, DEVICE_H);
  g_nonblack_effects = 0;
#define HS_SMOKE_ONE_DEV(name) smoke_one<name, DEVICE_W, DEVICE_H>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE_DEV)
#undef HS_SMOKE_ONE_DEV
  // The device <96,20> specialization is a distinct codepath; require it lights up too.
  HS_EXPECT_GT(g_nonblack_effects, 0);
#define HS_DET_ONE_DEV(name) determinism_one<name, DEVICE_W, DEVICE_H>(#name);
  HS_EFFECT_LIST(HS_DET_ONE_DEV)
#undef HS_DET_ONE_DEV

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test

