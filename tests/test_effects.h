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
#include "hardware/pov_segment_map.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <utility>
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
 * @brief Per-effect smoke frame count, resolved from HS_SMOKE_FRAMES.
 * @details Shared with the arena/stack budget gates (tests/test_fixture.h), so
 * they measure over the window this roster sweep renders. Frame count is a
 * minor cost lever here: the module is dominated by full-resolution software
 * raster (~71 ms/frame at 288x144 x 27 effects), which the QUICK tier skips
 * entirely (see effects_full_suite()).
 */
using hs_test::smoke_frames;

/**
 * @brief Selects the effects test depth tier from the environment.
 * @return true for the FULL suite (white-box correctness block + the
 * production-resolution 288x144 roster passes), false for the QUICK tier.
 * @details The full suite renders all 27 effects at the 288x144 production
 * resolution across a smoke pass, a determinism pass, and several full-frame
 * white-box tests — ~40 s, dominated by the ~71 ms/frame software raster of
 * 41,472-pixel frames (27 effects x N frames, back to back). The QUICK tier
 * (default) runs only the device-resolution <96,20> smoke + determinism passes,
 * ~1,920-pixel frames that cover every effect's construct/init/render/read-back
 * and cross-run determinism in ~2 s — the fast path for the local pre-commit
 * hook. CI opts into the full suite on every push/PR by setting HS_EFFECTS_FULL=1
 * (.github/workflows/ci.yml), so the full-resolution passes and the white-box
 * correctness block are the authoritative gate there, not locally. Set
 * HS_EFFECTS_FULL=1 to reproduce the CI depth in a local commit.
 */
inline bool effects_full_suite() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *e = std::getenv("HS_EFFECTS_FULL"))
#pragma clang diagnostic pop
    return std::atoi(e) != 0;
  return false;
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
// Compile-time streq, so the exemption name below is pinned to the roster.
constexpr bool roster_name_eq(const char *a, const char *b) {
  while (*a && *a == *b) {
    ++a;
    ++b;
  }
  return *a == *b;
}

// True when `name` is one of the HS_EFFECT_LIST roster class names.
constexpr bool is_roster_effect(const char *name) {
#define HS_ROSTER_NAME_MATCH(cls) || roster_name_eq(name, #cls)
  return false HS_EFFECT_LIST(HS_ROSTER_NAME_MATCH);
#undef HS_ROSTER_NAME_MATCH
}

inline bool effect_may_be_dark(const char *name, int frames) {
  // Renaming the effect class turns the exemption into a build error here
  // rather than a silently stale strcmp that drops the smoke assertion.
  static constexpr const char *kExemptDark = "RingShower";
  static_assert(is_roster_effect(kExemptDark),
                "all-black smoke exemption names a non-roster effect");
  if (std::strcmp(name, kExemptDark) == 0)
    return frames < 30;
  return false;
}

/**
 * @brief Resets the process-global effect state (RNG seed, arenas, timeline) to
 *        a clean per-effect baseline.
 * @details Every effect aliases the same static double buffer, global RNG,
 *          arenas, and timeline; without this reset leftover events reference
 *          the previous (destroyed) effect instance.
 */
inline void reset_effect_globals() {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
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
  reset_effect_globals();

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
 * mark_animated() and NOT flagged markReadonly() — must be genuinely
 * user-controllable, i.e. a value written through updateParameter() must persist
 * across frames. If the engine overwrites it every frame (a Mutation/Driver/Lerp
 * bound to the same member, or output-only telemetry), the slider is dead: the
 * author must drive a private member and mark_animated() it, or markReadonly()
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
                  "(mark_animated / markReadonly / drive a private member)\n",
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
 * @param frame_fold Optional out-param receiving an FNV-1a running checksum
 * folded over every displayed frame (not just the last). Lets a caller detect
 * mid-run nondeterminism that reconverges before the final frame, which the
 * final-buffer copy alone cannot see. Ignored when nullptr.
 * @details Resets every shared global the smoke path does (RNG seed, arenas,
 * Timeline) plus the generative-hue cursor and mock clock, so two calls start
 * from an identical state.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void render_capture(std::vector<Pixel> &out, int frames,
                           uint64_t *frame_fold = nullptr) {
  reset_effect_globals();
  // Pin the global generative-hue cursor so both runs start identical (it
  // drifts across palette constructions by design; see GenerativePalette).
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);

  uint64_t fold = 1469598103934665603ull; // FNV-1a offset basis
  const auto fold_byte = [&fold](uint8_t byte) {
    fold ^= byte;
    fold *= 1099511628211ull; // FNV-1a prime
  };

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    effect.draw_frame();
    effect.advance_display();
    if (frame_fold)
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel p = effect.get_pixel(x, y);
          fold_byte(p.r & 0xFF);
          fold_byte(p.r >> 8);
          fold_byte(p.g & 0xFF);
          fold_byte(p.g >> 8);
          fold_byte(p.b & 0xFF);
          fold_byte(p.b >> 8);
        }
  }
  if (frame_fold)
    *frame_fold = fold;

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
  uint64_t fold_a = 0, fold_b = 0;
  render_capture<E, W, H>(a, frames, &fold_a);
  perturb_determinism_globals();
  render_capture<E, W, H>(b, frames, &fold_b);
  hs::clear_mock_time();

  // Per-frame fold catches mid-run divergence that reconverges by the final
  // frame; the final-buffer comparison below alone would miss it.
  if (fold_a != fold_b)
    std::printf("  NONDETERMINISTIC %-20s per-frame checksum %llu != %llu over "
                "%d frames\n",
                name, static_cast<unsigned long long>(fold_a),
                static_cast<unsigned long long>(fold_b), frames);
  HS_EXPECT(fold_a == fold_b,
            "effect must render identically every frame across runs under a "
            "fixed clock");

  int first_diff = -1;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i].r != b[i].r || a[i].g != b[i].g || a[i].b != b[i].b) {
      first_diff = static_cast<int>(i);
      break;
    }

  if (first_diff >= 0) {
    const int x = first_diff % W, y = first_diff / W;
    std::printf(
        "  NONDETERMINISTIC %-20s pixel (%d,%d): runA (%d,%d,%d) != "
        "runB (%d,%d,%d) over %d frames\n",
        name, x, y, static_cast<int>(a[first_diff].r),
        static_cast<int>(a[first_diff].g), static_cast<int>(a[first_diff].b),
        static_cast<int>(b[first_diff].r), static_cast<int>(b[first_diff].g),
        static_cast<int>(b[first_diff].b), frames);
  }
  HS_EXPECT(first_diff < 0,
            "effect must render identically across runs under a fixed clock");
}

/** @brief Frames per segment in the clip-clear parity sweep. */
constexpr int PARITY_FRAMES = 16;
/** @brief Arm segments walked by the clip-clear parity sweep. */
constexpr int PARITY_SEGMENTS = 4;

/**
 * @brief Drives one effect under a moving segment clip with each clear scope and
 *        requires the displayed pixels to agree.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param name Effect name used in the diagnostic output.
 * @details The clip clear leaves the region outside the display band holding
 *          pixels from the frame that last wrote this buffer. Anything an effect
 *          reads back or carries forward from there shows up here as a mismatch
 *          against the whole-buffer clear.
 */
template <template <int, int> class E, int W = DEVICE_W, int H = DEVICE_H>
inline void clip_clear_parity_one(const char *name) {
  constexpr int S = H * 2;

  auto render = [&](int segment_id, bool full_clear) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    std::vector<Pixel> displayed;
    E<W, H> effect;
    effect.init();
    effect.force_full_buffer_clear = full_clear;
    const pov::SegmentMap map =
        pov::segment_map(segment_id, S, PARITY_SEGMENTS);
    for (int f = 0; f < PARITY_FRAMES; ++f) {
      const pov::SegmentClip clip =
          pov::segment_clip(map, (f & 1) == 0, S, PARITY_SEGMENTS, W);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      for (int y = clip.y0; y < clip.y1; ++y)
        for (int x = clip.x0; x < clip.x1; ++x)
          displayed.push_back(effect.get_pixel(x, y));
    }
    return displayed;
  };

  for (int segment_id = 0; segment_id < PARITY_SEGMENTS; ++segment_id) {
    const std::vector<Pixel> full = render(segment_id, true);
    const std::vector<Pixel> clipped = render(segment_id, false);
    HS_EXPECT_EQ(full.size(), clipped.size());

    size_t different = 0;
    for (size_t i = 0; i < full.size() && i < clipped.size(); ++i)
      if (full[i] != clipped[i])
        ++different;
    if (different)
      std::printf("  CLIP-CLEAR DRIFT %-20s segment %d: %zu of %zu displayed "
                  "pixels differ from the full-buffer clear\n",
                  name, segment_id, different, full.size());
    HS_EXPECT(different == 0,
              "clip clearing must not change any displayed pixel");
  }
  hs::clear_mock_time();
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
    HS_EXPECT_EQ(l * l + l + m, idx);  // round-trips the index
    HS_EXPECT_TRUE(m >= -l && m <= l); // order within the level's band
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

/** @brief Textbook associated Legendre P_l^m(x) in double, as test truth. */
inline double sh_reference_legendre(int l, int m, double x) {
  double pmm = 1.0;
  if (m > 0) {
    double somx2 = std::sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
    double fact = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;

  double pmmp1 = x * (2.0 * m + 1.0) * pmm;
  if (l == m + 1)
    return pmmp1;

  double pll = 0.0;
  for (int ll = m + 2; ll <= l; ll++) {
    pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

/**
 * @brief Verifies the Cartesian spherical harmonic matches the spherical form.
 * @details SHMath::spherical_harmonic drops the sin(phi)^|m| factor from the
 * Legendre term and recovers it from Re/Im((x + iz)^|m|). That identity holds
 * only for unit directions, so this sweeps a full (phi, theta) lattice against
 * a double-precision P_l^m(cos phi) * cos/sin(|m| theta), covering every (l, m)
 * the effect can reach.
 */
inline void test_sh_cartesian_matches_spherical() {
  double worst = 0.0;
  int worst_l = 0, worst_m = 0;

  for (int l = 0; l <= 4; ++l) {
    for (int m = -l; m <= l; ++m) {
      const float N = SHMath::normalization(l, m);
      const int abs_m = std::abs(m);
      for (int i = 0; i <= 40; ++i) {
        const double phi = PI_F * i / 40.0;
        const double sin_phi = std::sin(phi), cos_phi = std::cos(phi);
        for (int j = 0; j < 64; ++j) {
          const double theta = 2.0 * PI_F * j / 64.0;
          Vector p(static_cast<float>(sin_phi * std::cos(theta)),
                   static_cast<float>(cos_phi),
                   static_cast<float>(sin_phi * std::sin(theta)));

          const double azimuthal = (m > 0)   ? std::cos(abs_m * theta)
                                   : (m < 0) ? std::sin(abs_m * theta)
                                             : 1.0;
          const double want =
              N * sh_reference_legendre(l, abs_m, cos_phi) * azimuthal;
          const double err =
              std::fabs(SHMath::spherical_harmonic(l, m, p, N) - want);
          if (err > worst) {
            worst = err;
            worst_l = l;
            worst_m = m;
          }
        }
      }
    }
  }

  if (worst >= 1e-4)
    std::printf("  SH parity: worst error %g at l=%d m=%d\n", worst, worst_l,
                worst_m);
  HS_EXPECT(worst < 1e-4,
            "Cartesian harmonic must match the spherical form on unit vectors");
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
    gs.params.d_a = dA;
    gs.params.d_b = dB;
    gs.params.dt = dt;
  }
  static int dissolve_frame(const GS &gs) {
    return gs.transition.dissolve_frame;
  }
  static const uint16_t *b_field(const GS &gs) { return gs.state.B; }
  static constexpr int DISSOLVE_FRAMES = GS::DISSOLVE_FRAMES;
  static constexpr int MIN_GROW_FRAMES = GS::MIN_GROW_FRAMES;
  static constexpr int STABLE_HOLD_FRAMES = GS::STABLE_HOLD_FRAMES;

  static void step(GS &gs, const uint16_t *cA, const uint16_t *cB, uint16_t *nA,
                   uint16_t *nB) {
    // The production substep is float-resident; this seam keeps the tests'
    // Q16-in/Q16-out contract by converting at the edges, as render() does
    // once per frame.
    std::vector<float> fA(N), fB(N), gA(N), gB(N);
    for (int i = 0; i < N; ++i) {
      fA[i] = GS::from_q16(cA[i]);
      fB[i] = GS::from_q16(cB[i]);
    }
    gs.step_physics(fA.data(), fB.data(), gA.data(), gB.data());
    for (int i = 0; i < N; ++i) {
      nA[i] = GS::to_q16(gA[i]);
      nB[i] = GS::to_q16(gB[i]);
    }
  }

  // Float-domain substep: the Q16 edges of step() clamp on their own, so a test
  // that means to observe the substep's own bound has to read it undigested.
  static void step_float(GS &gs, const float *cA, const float *cB, float *nA,
                         float *nB) {
    gs.step_physics(cA, cB, nA, nB);
  }

  // Drives the shared substep driver over render()'s float generation buffers.
  static void advance_substeps(GS &gs, int steps, float *pA, float *pB,
                               float *sA, float *sB) {
    gs.advance_substeps(steps, std::array<float *, 2>{pA, pB},
                        std::array<float *, 2>{sA, sB},
                        [&](auto &cur, auto &nxt) {
                          gs.step_physics(cur[0], cur[1], nxt[0], nxt[1]);
                        });
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
  HS_EXPECT_EQ(GSWhiteBox::to_q16(2.0f), (uint16_t)65535); // clamp high
  HS_EXPECT_EQ(GSWhiteBox::to_q16(-0.5f), (uint16_t)0);    // clamp low
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(65535), 1.0f, 1e-9f);
  int bad = 0;
  for (int v = 0; v <= 65535; ++v)
    if (GSWhiteBox::to_q16(GSWhiteBox::from_q16((uint16_t)v)) != (uint16_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies the dissolve clears the sphere and reseeds a fresh reaction.
 * @details Drives a real effect to its stabilization transition, then pins the
 *          two things the dissolve must do: coverage falls monotonically enough
 *          to reach near-zero by the end of the window (a scattered rest
 *          node is refilled by its autocatalytic neighbours in one frame, so a
 *          dissolve that only converts the newly-crossed band never clears at
 *          all), and the field that comes back is the seed pattern, not the old
 *          one.
 */
inline void test_gs_dissolve_clears_and_reseeds() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  auto covered = [&]() {
    const uint16_t *b = GSWhiteBox::b_field(gs);
    int n = 0;
    for (int i = 0; i < GSWhiteBox::N; ++i)
      if (b[i] >= GSWhiteBox::to_q16(0.1f))
        ++n;
    return n;
  };

  // Run until the stalled field trips a dissolve; the default reaction settles
  // well inside this budget.
  const int budget =
      GSWhiteBox::MIN_GROW_FRAMES + GSWhiteBox::STABLE_HOLD_FRAMES + 400;
  int f = 0;
  for (; f < budget && GSWhiteBox::dissolve_frame(gs) < 0; ++f) {
    gs.draw_frame();
    gs.advance_display();
  }
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "stalled field did not start a dissolve");
  const int grown = covered();
  HS_EXPECT(grown > GSWhiteBox::N / 20,
            "default reaction grew no pattern to dissolve");

  // Step to the last frame before the window closes: the sphere must be all but
  // clear, and still mid-dissolve.
  for (int i = 0; i < GSWhiteBox::DISSOLVE_FRAMES - 1; ++i) {
    gs.draw_frame();
    gs.advance_display();
  }
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "dissolve ended before its window closed");
  HS_EXPECT(covered() < grown / 8,
            "dissolve left the sphere covered: cleared nodes healed");

  // The closing frame converts the remainder and reseeds.
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  HS_EXPECT(covered() > 0, "reseed left no nucleation sites");
}

/**
 * @brief Verifies editing the reaction constants starts a dissolve.
 * @details Feed/Kill/dA/dB define the reaction, so an edit clears and restarts
 *          it; Speed only sets the integration rate and must not.
 */
inline void test_gs_reaction_edit_starts_dissolve() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);

  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 1.0f); // Speed only
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);

  GSWhiteBox::set_params(gs, 0.03f, 0.06f, 0.02f, 0.01f, 1.0f); // Feed moved
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "a Feed edit did not start a dissolve");
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
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.05f, 0.05f,
                         3.0f); // stable extreme
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

/**
 * @brief Verifies the substep stays finite and inside [0, 1] at the joint
 *        feed/k corner of the slider box, past the Euler stability bound.
 * @details The dA/dB cap covers only the diffusion term; the reaction Jacobian's
 *          rows exceed 2 at feed = k = 0.1 with Speed at top, so the scheme is
 *          genuinely unstable there and the per-substep clamp is the only thing
 *          holding the field. Assert what that buys: after 256 substeps from
 *          seeded nuclei every node is still finite and on [0, 1], so the
 *          overshoot saturates rather than escaping into NaN/Inf and poisoning
 *          the palette lookup downstream. Read in float — step()'s Q16 edges
 *          would clamp the evidence away.
 */
inline void test_gs_reaction_corner_stays_bounded() {
  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f),
      sa(GSWhiteBox::N), sb(GSWhiteBox::N);
  for (int s : {500, 2500, 4500, 6500})
    b[s] = 1.0f;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.1f, 0.1f, 0.05f, 0.05f,
                         3.0f); // joint slider-box corner
  float *cA = a.data(), *cB = b.data(), *nA = sa.data(), *nB = sb.data();
  for (int s = 0; s < 256; ++s) {
    GSWhiteBox::step_float(gs, cA, cB, nA, nB);
    std::swap(cA, nA);
    std::swap(cB, nB);
  }
  int escaped = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (!std::isfinite(cA[i]) || !std::isfinite(cB[i]) || cA[i] < 0.0f ||
        cA[i] > 1.0f || cB[i] < 0.0f || cB[i] > 1.0f)
      ++escaped;
  HS_EXPECT_EQ(escaped, 0);
}

/**
 * @brief Verifies an odd substep count lands the final generation back in the
 *        driver's persistent buffers.
 * @details advance_substeps ping-pongs between the caller's persistent buffers
 *          and scratch, then copies back only when the final generation ended up
 *          in scratch — the parity copy-back, which the even production
 *          STEPS_PER_FRAME never exercises. Drive a single (odd) substep on a
 *          seeded field: the diffused result must appear in the persistent B
 *          buffer (not just the scratch one), so a missing copy-back leaves the
 *          neighbors at rest and fails here.
 */
inline void test_gs_odd_substep_lands_in_persistent() {
  std::vector<float> pA(GSWhiteBox::N, 1.0f), pB(GSWhiteBox::N, 0.0f),
      sA(GSWhiteBox::N, 0.0f), sB(GSWhiteBox::N, 0.0f);
  const int seed = 4000; // an interior lattice node with a full neighbor ring
  pB[seed] = 1.0f;

  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::advance_substeps(gs, 1, pA.data(), pB.data(), sA.data(),
                               sB.data());

  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && pB[nb] > 0.0f)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0); // diffusion landed in the persistent buffers
}

// ---------------------------------------------------------------------------
// Belousov-Zhabotinsky reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for BZReactionDiffusion's private fixed-point and
 *        physics internals (befriended in effects/BZReactionDiffusion.h).
 * @details The matching seam to GSWhiteBox: the smoke/determinism harness cannot
 *          see a Q16 round-trip error, a sign/clamp slip in the Lotka-Volterra
 *          update, or a perturbation that wraps past the Q16 rail, so the
 *          conversions, one species step, the perturbation, and one fused physics
 *          substep are pinned directly. The lattice is the same fixed 7680-node
 *          graph as GS, independent of <W,H>, so the device resolution is used
 *          arbitrarily.
 */
struct BZWhiteBox {
  using BZ = BZReactionDiffusion<DEVICE_W, DEVICE_H>;
  static constexpr int N = BZ::RD_N;

  static uint16_t to_q16(float v) { return BZ::to_q16(v); }
  static float from_q16(uint16_t v) { return BZ::from_q16(v); }
  static void set_params(BZ &bz, float alpha, float D, float dt) {
    bz.params.alpha = alpha;
    bz.params.D = D;
    bz.params.dt = dt;
  }
  static uint16_t advance_species(const BZ &bz, float conc, float predator,
                                  float laplacian) {
    return bz.advance_species(conc, predator, laplacian);
  }
  static void perturb(uint16_t *nA, uint16_t *nB, uint16_t *nC) {
    BZ::perturb_state(nA, nB, nC);
  }
  static int num_perturbations() { return BZ::NUM_PERTURBATIONS; }
  static int perturb_amount() { return BZ::PERTURB_AMOUNT; }
  static void step(BZ &bz, uint16_t *sA, uint16_t *sB, uint16_t *sC) {
    std::vector<float> fA(N), fB(N), fC(N);
    bz.step_physics(sA, sB, sC, fA.data(), fB.data(), fC.data());
  }
};

/**
 * @brief Verifies the Q16 fixed-point round-trip and the +0.5 rounding/clamp
 *        boundaries.
 * @details to_q16(from_q16(v)) must be the identity over every representable
 *          value, and to_q16 must clamp out-of-range floats and round to nearest
 *          (so 1.0 tops out at 65535 with no overflow). A truncating or unclamped
 *          regression — which would bias the RD dynamics downward — fails here.
 */
inline void test_bz_q16_roundtrip() {
  HS_EXPECT_EQ(BZWhiteBox::to_q16(0.0f), (uint16_t)0);
  HS_EXPECT_EQ(BZWhiteBox::to_q16(1.0f), (uint16_t)65535);
  HS_EXPECT_EQ(BZWhiteBox::to_q16(2.0f), (uint16_t)65535); // clamp high
  HS_EXPECT_EQ(BZWhiteBox::to_q16(-0.5f), (uint16_t)0);    // clamp low
  HS_EXPECT_NEAR(BZWhiteBox::from_q16(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(BZWhiteBox::from_q16(65535), 1.0f, 1e-9f);
  int bad = 0;
  for (int v = 0; v <= 65535; ++v)
    if (BZWhiteBox::to_q16(BZWhiteBox::from_q16((uint16_t)v)) != (uint16_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies the state resolution actually carries a slider-minimum
 *        diffusion step instead of rounding it away.
 * @details The finest per-substep move the sliders can ask for is a cold node
 *          beside one saturated neighbour at Diff's 0.001 floor: D·lap·dt =
 *          0.001·1·0.35 = 3.5e-4 of full scale. A store whose half-LSB exceeds
 *          that discards it, the node never changes, and the lattice degrades
 *          into uncoupled per-node ODEs while the reaction term keeps running.
 *          Drive advance_species at exactly that operating point and require the
 *          stored sample to move.
 */
inline void test_bz_min_diffusion_step_survives_quantization() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.001f, /*dt*/ 0.35f);
  // Cold node, one saturated neighbour out of RD_K: lap = 1, reaction term 0.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, /*lap*/ 1.0f),
               0);
  // The same gradient one LSB above the floor still resolves, so the response
  // is graded rather than a single lucky rounding boundary.
  const float lsb = BZWhiteBox::from_q16(1);
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, lsb, 0.0f, /*lap*/ 1.0f),
               (int)BZWhiteBox::to_q16(lsb));
}

/**
 * @brief Verifies advance_species has the right reaction/diffusion signs and
 *        that the Q16 clamp backstop holds even past the Euler stability bound.
 * @details advance_species is the single-species core of the BZ update:
 *          conc + (D·laplacian + conc·(1 − conc − α·predator))·dt, mapped through
 *          to_q16. Checks: an empty rest cell stays empty (no spurious growth);
 *          diffusion from higher neighbors grows an empty cell; logistic growth
 *          lifts a half-filled, predator-free cell; predation underflow clamps to
 *          0 (not a uint16 wrap to 65535); and extreme over-/under-shoots — the
 *          documented to_q16 backstop the comment in the effect promises — clamp
 *          to the [0, 65535] rails rather than wrapping.
 */
inline void test_bz_advance_species_signs_and_clamp() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 0.35f);

  // Empty rest cell: no diffusion, no reaction -> stays 0.
  HS_EXPECT_EQ((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, 0.0f), 0);

  // Diffusion from higher neighbors lifts an empty cell above 0.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.0f, 0.0f, /*lap*/ 6.0f),
               0);

  // Logistic growth: a half-filled, predator-free cell grows above its start.
  HS_EXPECT_GT((int)BZWhiteBox::advance_species(bz, 0.5f, 0.0f, 0.0f),
               (int)BZWhiteBox::to_q16(0.5f));

  // Predation drives a saturated cell negative; it must clamp to 0, not wrap.
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1.0f, 0.0f), 0);

  // Backstop past the Euler bound: a huge positive laplacian clamps to 65535, a
  // huge predation clamps to 0 — to_q16 keeps every written state in range.
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, 0.0f, /*lap*/ 1000.0f), 65535);
  HS_EXPECT_EQ(
      (int)BZWhiteBox::advance_species(bz, 1.0f, /*predator*/ 1000.0f, 0.0f),
      0);
}

/**
 * @brief Verifies perturb_state nudges nodes by a fixed Q16 amount and saturates
 *        at the 65535 rail without wrapping.
 * @details Two RNG-agnostic invariants. (1) A fully-saturated field stays fully
 *          saturated: every nudge is a +PERTURB_AMOUNT saturating add, so a node
 *          already at 65535 cannot wrap to a low value. (2) On a zero field, each
 *          touched entry is a small positive multiple of the nudge step and stays
 *          within [0, 65535], and at least one entry is touched — whichever nodes
 *          the deterministic RNG happens to select.
 */
inline void test_bz_perturb_state_saturates_and_nudges() {
  // (1) Saturation / no-wrap: all rails stay at the rail.
  {
    std::vector<uint16_t> a(BZWhiteBox::N, 65535), b(BZWhiteBox::N, 65535),
        c(BZWhiteBox::N, 65535);
    BZWhiteBox::perturb(a.data(), b.data(), c.data());
    int wrapped = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      if (a[i] != 65535 || b[i] != 65535 || c[i] != 65535)
        ++wrapped;
    HS_EXPECT_EQ(wrapped, 0);
  }
  // (2) Zero field: touched entries are small positive multiples of the step.
  {
    const int step = BZWhiteBox::perturb_amount();
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    BZWhiteBox::perturb(a.data(), b.data(), c.data());
    int touched = 0, malformed = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      for (uint16_t v : {a[i], b[i], c[i]}) {
        if (v == 0)
          continue;
        ++touched;
        if (v % step != 0) // accumulations stay multiples of the nudge step
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
  std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
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
 *        neighborhood with the right sign, in place.
 * @details Seed species A fully at one interior node on an otherwise-empty field
 *          and run a single step. The seed's neighbors start empty but border a
 *          saturated node, so their graph-Laplacian is positive and A must
 *          diffuse into at least one of them; the seed node itself stays lit
 *          (it decays but does not vanish or wrap in one step). The step writes
 *          the new generation over the state it read, so the same buffers carry
 *          the result — a Laplacian reading the half-written state instead of
 *          the float mirror would spread the seed asymmetrically along node
 *          order. The stochastic perturbation runs too, but the seed-neighbor
 *          diffusion is independent of which nodes it nudges.
 */
inline void test_bz_substep_diffuses() {
  std::vector<uint16_t> sA(BZWhiteBox::N, 0), sB(BZWhiteBox::N, 0),
      sC(BZWhiteBox::N, 0);
  const int seed = 4000; // interior lattice node with a full neighbor ring
  sA[seed] = 65535;

  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, 3.0f, 0.05f, 0.35f);
  BZWhiteBox::step(bz, sA.data(), sB.data(), sC.data());

  HS_EXPECT_GT((int)sA[seed], 0); // the seed decays but does not vanish/wrap
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && sA[nb] > 0)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0); // A diffused into at least one empty neighbor
}

// ---------------------------------------------------------------------------
// DreamBalls: preset-cycle / re-spawn white-box coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for DreamBalls' private preset-cycle bookkeeping
 *        (befriended in effects/DreamBalls.h).
 * @details spawn_sprite schedules its successor 320 frames out (a PeriodicTimer),
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
  // Not-paused step: advance the selector, then re-spawn (the scheduler's path).
  static void advance(DB &db) {
    db.preset_manager.next();
    db.spawn_sprite();
  }
  // Paused re-spawn of the current preset (no advance).
  static void respawn(DB &db) { db.spawn_sprite(); }
  // The literal solid-name pointer the given preset seeds into params on reseed.
  static const char *preset_name(const DB &db, int idx) {
    return db.preset_manager.get_entries()[idx].params.solid_name;
  }
  static const char *live_solid(const DB &db) { return db.params.solid_name; }
  static float &num_copies(DB &db) { return db.params.num_copies; }
};

/**
 * @brief Drives spawn_sprite across a full preset cycle and asserts the bake-slot
 *        ping-pong, the modulo preset advance, and the reseed-on-change guard.
 * @details Drives the advance without the 288-frame wait, following the same
 *          selector progression the periodic callback uses: each step calls
 *          Presets::next() then re-spawns, so the active preset walks step % 4.
 *          Each spawn must flip the bake slot (so a fading-out sprite keeps its
 *          own LUT) and, when the preset actually changes, reseed params to the
 *          new entry. A re-spawn of the SAME preset (the paused branch) must
 *          instead hold params so a live slider edit survives.
 */
inline void test_dreamballs_preset_cycle_bookkeeping() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init(); // runs spawn_sprite() at preset 0

  // init() spawned preset 0: it reseeded params (last_preset_idx_ -1 -> 0) and
  // flipped the bake slot once (0 -> 1).
  HS_EXPECT_EQ(WB::last_preset_idx(db), 0);
  HS_EXPECT_EQ(WB::active_bake(db), 1);
  HS_EXPECT_EQ(WB::live_solid(db), WB::preset_name(db, 0));

  // Not-paused advance chain: each step advances the selector then re-spawns, so
  // the preset is step % 4. Drive a full cycle plus a wrap; the bake slot must
  // ping-pong every spawn and params must reseed to the new index each step.
  int expect_bake = WB::active_bake(db); // 1
  for (int step = 1; step <= 8; ++step) {
    WB::advance(db);
    expect_bake ^= 1;
    const int safe = step % WB::PRESETS;
    HS_EXPECT_EQ(WB::active_bake(db), expect_bake);
    HS_EXPECT_EQ(WB::last_preset_idx(db), safe);
    HS_EXPECT_EQ(WB::live_solid(db), WB::preset_name(db, safe));
  }

  // Paused-hold path: re-spawn the SAME preset (no advance), so the reseed guard
  // (safe_idx == last_preset_idx_) holds and a live slider edit must survive the
  // re-spawn — while the bake slot still flips. last_preset_idx_ is now 0.
  const float sentinel = WB::num_copies(db) + 5.0f;
  WB::num_copies(db) = sentinel;
  const int held_idx = WB::last_preset_idx(db);
  expect_bake ^= 1;
  WB::respawn(db);
  HS_EXPECT_EQ(WB::last_preset_idx(db), held_idx); // preset unchanged
  HS_EXPECT_EQ(WB::num_copies(db), sentinel);      // live edit preserved
  HS_EXPECT_EQ(WB::active_bake(db), expect_bake);  // bake slot still flipped
}

/**
 * @brief End-to-end check that the 320-frame re-spawn timer actually fires and
 *        honors the pause gate, by rendering past one period.
 * @details Covers the part the white-box driver cannot: the PeriodicTimer wiring
 *          and the animations_paused()?idx:idx+1 hold-vs-advance decision in the
 *          scheduler lambda. Renders 340 frames (one period is 320; the next
 *          re-spawn is another 320 away, so exactly one fires). Unpaused, the
 *          preset advances 0 -> 1; paused, it re-spawns the same preset and holds.
 */
inline void test_dreamballs_respawn_fires_and_honors_pause() {
  using WB = DreamBallsWhiteBox;
  auto run = [](bool paused) {
    reset_effect_globals();
    WB::DB db;
    db.init();
    db.setAnimationsPaused(paused);
    // The re-spawn timer is scheduled 320 frames out; render past it so exactly
    // one re-spawn fires (the following one is another 320 frames away).
    for (int f = 0; f < 340; ++f) {
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
    for (const PresetEntry<LissajousParams> &entry : C::FUNCTIONS) {
      const LissajousParams &cfg = entry.params;
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
    reset_effect_globals();

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

/**
 * @brief White-box accessor for HopfFibration's S3-lift + stereographic
 *        projection (befriended in effects/HopfFibration.h).
 * @details The smoke/determinism harness only proves the effect renders and
 *          reproduces; it never pins hopf_project()'s numeric output. This seam
 *          sets the private per-frame cache (tumble sines/cosines, fold/flow/
 *          tumble-y phases) and a fiber's base coordinates, then calls the real
 *          projection so a test can compare it to the closed form.
 */
struct HopfWhiteBox {
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;
  static size_t fiber_count() { return HF::ACTUAL_FIBERS; }
  static Vector project(HF &fx, size_t i) { return fx.hopf_project(i); }
  static void set_fiber(HF &fx, size_t i, float azimuth, float polar) {
    fx.fibers[i] = Spherical(azimuth, polar);
  }
  static void set_cache(HF &fx, float cx, float sx, float cy, float sy,
                        float fold_base, float flow_rad, float ty_rad) {
    fx.cx = cx;
    fx.sx = sx;
    fx.cy = cy;
    fx.sy = sy;
    fx.fold_base = fold_base;
    fx.flow_rad = flow_rad;
    fx.ty_rad = ty_rad;
  }
  static void set_folding(HF &fx, float v) { fx.params.folding = v; }
  static void set_twist(HF &fx, float v) { fx.params.twist = v; }
};

/**
 * @brief Pins HopfFibration::hopf_project against a closed form and asserts its
 *        unit-direction/finite invariant under a nontrivial 4D tumble.
 * @details Identity tumble with zero folding and twist reduces fiber 0
 *          (beta == 0) to a plain S3 lift, whose stereographic image is the
 *          normalized (q0, q1, q2); the check uses the same fast trig the effect
 *          does so it pins the projection, not the trig approximation. The second
 *          pass exercises every fiber under active tumble/fold/twist and requires
 *          each result be finite, unit-length (normalized_or), and deterministic.
 */
inline void test_hopf_projection_math() {
  // init() bakes from persistent_arena and schedules on the shared timeline, so
  // reset the shared globals as smoke_one does.
  reset_effect_globals();
  using WB = HopfWhiteBox;

  HopfFibration<DEFAULT_W, DEFAULT_H> fx;
  fx.init();

  // Identity tumble + no folding/twist: fiber 0 has beta == 0, so q3 == 0 and the
  // projection collapses to the normalized (q0, q1, q2).
  WB::set_folding(fx, 0.0f);
  WB::set_twist(fx, 0.0f);
  WB::set_cache(fx, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);

  const float polar = 1.2f, azimuth = 0.7f;
  WB::set_fiber(fx, 0, azimuth, polar);
  const Vector v = WB::project(fx, 0);

  const float eta = polar * 0.5f;
  const float ce = fast_cosf(eta), se = fast_sinf(eta);
  const Vector expected =
      Vector(ce * fast_cosf(azimuth), ce * fast_sinf(azimuth), se).normalized();
  HS_EXPECT_NEAR(v.x, expected.x, 1e-4f);
  HS_EXPECT_NEAR(v.y, expected.y, 1e-4f);
  HS_EXPECT_NEAR(v.z, expected.z, 1e-4f);
  HS_EXPECT_NEAR(v.magnitude(), 1.0f, 1e-4f);

  // Every fiber projects to a finite unit direction under active tumble/fold/twist.
  const float ax = 0.9f, ay = 0.4f;
  WB::set_cache(fx, fast_cosf(ax), fast_sinf(ax), fast_cosf(ay), fast_sinf(ay),
                fast_sinf(ax * 0.5f) * 0.5f, 0.3f, ay);
  WB::set_folding(fx, 0.5f);
  WB::set_twist(fx, 2.0f);
  for (size_t i = 0; i < WB::fiber_count(); ++i) {
    const Vector p = WB::project(fx, i);
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-3f);
    const Vector again =
        WB::project(fx, i); // pure fn of the cache -> repeatable
    HS_EXPECT_NEAR(p.x, again.x, 1e-6f);
    HS_EXPECT_NEAR(p.y, again.y, 1e-6f);
    HS_EXPECT_NEAR(p.z, again.z, 1e-6f);
  }
}

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
 * @brief White-box accessor for DisplacementField's hue-table bake.
 */
struct DisplacementFieldWhiteBox {
  template <int W, int H>
  static void prepare_hue_table(DisplacementField<W, H> &effect,
                                const Color4 &color, float domain) {
    effect.prepare_hue_table(make_hue_rotate_base(color), domain);
  }

  template <int W, int H>
  static Pixel sample_hue_table(const DisplacementField<W, H> &effect,
                                float amount, float domain, bool cyclic) {
    return effect.sample_hue_table(amount, domain, cyclic);
  }

  template <int W, int H>
  static Pixel sample_hue_table_cached(DisplacementField<W, H> &effect,
                                       float amount, float domain, bool cyclic,
                                       const HueRotateBase &base,
                                       uint64_t *valid) {
    return effect.sample_hue_table_cached(amount, domain, cyclic, base, valid);
  }

  template <int W, int H>
  static int hue_table_size(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::HUE_TABLE_SIZE;
  }

  template <int W, int H>
  static Pixel hue_table_value(const DisplacementField<W, H> &effect,
                               int index) {
    return effect.hue_table[index];
  }

  template <int W, int H>
  static void clear_hue_table(DisplacementField<W, H> &effect) {
    for (int i = 0; i <= DisplacementField<W, H>::HUE_TABLE_SIZE; ++i)
      effect.hue_table[i] = Pixel(0, 0, 0);
  }

  template <int W, int H>
  static void configure_noise(DisplacementField<W, H> &effect, float rings,
                              float hue_scale) {
    effect.params.num_rings = rings;
    effect.params.hue_scale = hue_scale;
    effect.master_gain = 1.0f;
  }

  template <int W, int H>
  static void set_force_exact_hue(DisplacementField<W, H> &effect, bool exact) {
    effect.force_exact_hue = exact;
  }

  template <int W, int H>
  static int hue_table_uses(const DisplacementField<W, H> &effect) {
    return effect.hue_table_uses;
  }

  template <int W, int H>
  static Pixel hue_lut_value(const DisplacementField<W, H> &effect, int index) {
    // Slot 0 of the pooled bake: the first (only, in the one-ring tests)
    // drawn ring's hue LUT.
    return effect.hue_pool[index];
  }

  template <int W, int H>
  static Color4 current_ring_color(const DisplacementField<W, H> &effect,
                                   float ring_t) {
    return effect.palette.get(wrap_t(ring_t + effect.color_spin));
  }

  template <int W, int H>
  static int noise_lut_samples(const DisplacementField<W, H> &effect,
                               float sin_theta) {
    const float feature_scale = effect.params.scale1 + effect.params.scale2;
    return hs::clamp(
        static_cast<int>(ceilf(DisplacementField<W, H>::LUT_SAMPLES_PER_UNIT *
                               2.0f * PI_F * feature_scale * sin_theta)),
        DisplacementField<W, H>::LUT_MIN_SAMPLES, W);
  }
};

/**
 * @brief Verifies lazy table endpoints and outputs match eager construction.
 */
inline void test_displacement_field_lazy_hue_table_matches_eager() {
  reset_effect_globals();
  DisplacementField<DEVICE_W, DEVICE_H> effect;
  effect.init();
  GenerativePalette palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                            BrightnessProfile::FLAT);
  const Color4 color = palette.get(0.37f);
  const HueRotateBase base = make_hue_rotate_base(color);
  struct TableCase {
    float domain;
    float max_amount;
    bool cyclic;
  };
  const TableCase cases[] = {{0.4f, 0.4f, false},
                             {-0.4f, -0.4f, false},
                             {1.0f, 3.0f, true},
                             {-1.0f, -3.0f, true}};
  constexpr int SAMPLE_COUNT = 1024;
  const int table_size = DisplacementFieldWhiteBox::hue_table_size(effect);

  for (const TableCase &table_case : cases) {
    DisplacementFieldWhiteBox::prepare_hue_table(effect, color,
                                                 table_case.domain);
    std::vector<Pixel> endpoints(table_size + 1);
    for (int i = 0; i <= table_size; ++i)
      endpoints[i] = DisplacementFieldWhiteBox::hue_table_value(effect, i);
    std::vector<Pixel> expected(SAMPLE_COUNT + 1);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      expected[i] = DisplacementFieldWhiteBox::sample_hue_table(
          effect, amount, table_case.domain, table_case.cyclic);
    }

    DisplacementFieldWhiteBox::clear_hue_table(effect);
    uint64_t valid[2] = {0, 0};
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      Pixel actual = DisplacementFieldWhiteBox::sample_hue_table_cached(
          effect, amount, table_case.domain, table_case.cyclic, base, valid);
      HS_EXPECT_EQ(actual.r, expected[i].r);
      HS_EXPECT_EQ(actual.g, expected[i].g);
      HS_EXPECT_EQ(actual.b, expected[i].b);
    }
    HS_EXPECT_EQ(valid[0], ~uint64_t{0});
    HS_EXPECT_TRUE(valid[1] & uint64_t{1});
    for (int i = 0; i <= table_size; ++i) {
      Pixel actual = DisplacementFieldWhiteBox::hue_table_value(effect, i);
      HS_EXPECT_EQ(actual.r, endpoints[i].r);
      HS_EXPECT_EQ(actual.g, endpoints[i].g);
      HS_EXPECT_EQ(actual.b, endpoints[i].b);
    }
  }
}

/**
 * @brief Bounds dynamic and periodic hue tables over effect palette colors.
 */
inline void test_displacement_field_hue_table_fidelity() {
  reset_effect_globals();
  DisplacementField<DEVICE_W, DEVICE_H> effect;
  effect.init();

  constexpr float INV16 = 1.0f / 65535.0f;
  float default_delta_e = 0.0f;
  float cyclic_delta_e = 0.0f;
  int default_srgb_delta = 0;
  int cyclic_srgb_delta = 0;
  for (uint32_t seed = 0; seed < 12; ++seed) {
    GenerativePalette::reset_hue_seed(seed * 17u);
    GenerativePalette palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                              BrightnessProfile::FLAT);
    for (int color_index = 0; color_index < 48; ++color_index) {
      Color4 base = palette.get((color_index + 0.5f) / 48.0f);
      HueRotateBase exact_base = make_hue_rotate_base(base);
      for (int mode = 0; mode < 2; ++mode) {
        const bool cyclic = mode == 1;
        const float domain = cyclic ? 1.0f : 0.4f;
        const float max_amount = cyclic ? 3.0f : domain;
        DisplacementFieldWhiteBox::prepare_hue_table(effect, base, domain);
        for (int i = 0; i <= 1024; ++i) {
          const float amount = max_amount * i / 1024.0f;
          Pixel exact = hue_rotate(exact_base, amount).color;
          Pixel approx = DisplacementFieldWhiteBox::sample_hue_table(
              effect, amount, domain, cyclic);
          OKLab exact_lab = linear_rgb_to_oklab(
              exact.r * INV16, exact.g * INV16, exact.b * INV16);
          OKLab approx_lab = linear_rgb_to_oklab(
              approx.r * INV16, approx.g * INV16, approx.b * INV16);
          const float dl = exact_lab.L - approx_lab.L;
          const float da = exact_lab.a - approx_lab.a;
          const float db = exact_lab.b - approx_lab.b;
          const float delta_e = sqrtf(dl * dl + da * da + db * db);
          const int srgb_delta = std::max(
              std::abs(
                  static_cast<int>(linear_float_to_srgb8(exact.r * INV16)) -
                  linear_float_to_srgb8(approx.r * INV16)),
              std::max(std::abs(static_cast<int>(
                                    linear_float_to_srgb8(exact.g * INV16)) -
                                linear_float_to_srgb8(approx.g * INV16)),
                       std::abs(static_cast<int>(
                                    linear_float_to_srgb8(exact.b * INV16)) -
                                linear_float_to_srgb8(approx.b * INV16))));
          if (cyclic) {
            cyclic_delta_e = std::max(cyclic_delta_e, delta_e);
            cyclic_srgb_delta = std::max(cyclic_srgb_delta, srgb_delta);
          } else {
            default_delta_e = std::max(default_delta_e, delta_e);
            default_srgb_delta = std::max(default_srgb_delta, srgb_delta);
          }
        }
      }
    }
  }
  std::printf("  [hue table] default deltaE=%g sRGB8=%d; cyclic deltaE=%g "
              "sRGB8=%d\n",
              default_delta_e, default_srgb_delta, cyclic_delta_e,
              cyclic_srgb_delta);
  HS_EXPECT_LE(default_delta_e, 0.01f);
  HS_EXPECT_LE(default_srgb_delta, 10);
  HS_EXPECT_LE(cyclic_delta_e, 0.01f);
  HS_EXPECT_LE(cyclic_srgb_delta, 21);
  GenerativePalette::reset_hue_seed(0);
}

struct DisplacementHueFrame {
  std::vector<Pixel> pixels;
  int table_uses;
};

/**
 * @brief Renders a deterministic DisplacementField frame with table control.
 * @param exact Whether to bypass the hue table.
 * @return Final RGB16 framebuffer and table-use count.
 */
inline DisplacementHueFrame render_displacement_hue_frame(bool exact) {
  constexpr int W = 192;
  constexpr int H = 40;
  constexpr int FRAMES = 64;
  reset_effect_globals();
  global_timeline_t = 0;
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::set_force_exact_hue(effect, exact);
  for (int frame = 0; frame < FRAMES; ++frame) {
    hs::set_mock_time(frame * FRAME_MS, frame * FRAME_US);
    effect.draw_frame();
    effect.advance_display();
  }
  DisplacementHueFrame result{
      {}, DisplacementFieldWhiteBox::hue_table_uses(effect)};
  result.pixels.resize(W * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      result.pixels[y * W + x] = effect.get_pixel(x, y);
  hs::clear_mock_time();
  return result;
}

/**
 * @brief Bounds the table's final-frame delta against exact hue rotation.
 */
inline void test_displacement_field_hue_table_frame_fidelity() {
  DisplacementHueFrame exact = render_displacement_hue_frame(true);
  DisplacementHueFrame table = render_displacement_hue_frame(false);
  HS_EXPECT_EQ(exact.table_uses, 0);
  HS_EXPECT_GT(table.table_uses, 0);
  constexpr float INV16 = 1.0f / 65535.0f;
  float max_delta_e = 0.0f;
  int max_srgb_delta = 0;
  int changed = 0;
  for (size_t i = 0; i < exact.pixels.size(); ++i) {
    if (exact.pixels[i].r != table.pixels[i].r ||
        exact.pixels[i].g != table.pixels[i].g ||
        exact.pixels[i].b != table.pixels[i].b)
      ++changed;
    OKLab a = linear_rgb_to_oklab(exact.pixels[i].r * INV16,
                                  exact.pixels[i].g * INV16,
                                  exact.pixels[i].b * INV16);
    OKLab b = linear_rgb_to_oklab(table.pixels[i].r * INV16,
                                  table.pixels[i].g * INV16,
                                  table.pixels[i].b * INV16);
    const float dl = a.L - b.L;
    const float da = a.a - b.a;
    const float db = a.b - b.b;
    max_delta_e = std::max(max_delta_e, sqrtf(dl * dl + da * da + db * db));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].r * INV16)) -
            linear_float_to_srgb8(table.pixels[i].r * INV16)));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].g * INV16)) -
            linear_float_to_srgb8(table.pixels[i].g * INV16)));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].b * INV16)) -
            linear_float_to_srgb8(table.pixels[i].b * INV16)));
  }
  std::printf("  [hue frame] changed=%d/%zu uses=%d deltaE=%g sRGB8=%d\n",
              changed, exact.pixels.size(), table.table_uses, max_delta_e,
              max_srgb_delta);
  HS_EXPECT_LE(max_delta_e, 0.01f);
  HS_EXPECT_LE(max_srgb_delta, 3);
}

/**
 * @brief Verifies zero hue scale reproduces one exact zero rotation per ring.
 */
inline void test_displacement_field_zero_hue_scale_is_exact() {
  constexpr int W = 256;
  constexpr int H = 40;
  reset_effect_globals();
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::configure_noise(effect, 1.0f, 0.0f);
  effect.draw_frame();

  Color4 base = DisplacementFieldWhiteBox::current_ring_color(effect, 0.5f);
  Pixel expected = hue_rotate(make_hue_rotate_base(base), 0.0f).color;
  const int lut_n = DisplacementFieldWhiteBox::noise_lut_samples(effect, 1.0f);
  for (int i = 0; i <= lut_n; ++i) {
    Pixel actual = DisplacementFieldWhiteBox::hue_lut_value(effect, i);
    HS_EXPECT_EQ(actual.r, expected.r);
    HS_EXPECT_EQ(actual.g, expected.g);
    HS_EXPECT_EQ(actual.b, expected.b);
  }
  effect.advance_display();
  hs::clear_mock_time();
}

/**
 * @brief Verifies DisplacementField's clipped render tiles the full render:
 *        under a quadrant clip, every display-region pixel matches the
 *        full-canvas render bit-exactly, frame by frame.
 * @details Identical seeds and mock clock per run, so a divergence isolates
 *          the clip-only paths (the per-ring cap cull and the azimuth-chunk
 *          bake cull) dropping a reachable fragment or sampling a stale LUT
 *          entry. The frame window sits inside the ball phase, whose
 *          footprints drive both culls.
 */
inline void test_displacement_field_clip_tiles_full() {
  struct Quad {
    int x0, x1, y0, y1;
  };
  const Quad quads[] = {{0, DEFAULT_W / 2, 0, DEFAULT_H / 2},
                        {DEFAULT_W / 2, DEFAULT_W, DEFAULT_H / 2, DEFAULT_H}};
  const int frames = 60;

  auto fold_region = [&](bool clip, const Quad &q) -> uint64_t {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    DisplacementField<DEFAULT_W, DEFAULT_H> fx;
    fx.init();
    if (clip)
      fx.set_clip(q.y0, q.y1, q.x0, q.x1);
    uint64_t fold = 1469598103934665603ull; // FNV-1a offset basis
    for (int f = 0; f < frames; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      fx.draw_frame();
      fx.advance_display();
      for (int y = q.y0; y < q.y1; ++y)
        for (int x = q.x0; x < q.x1; ++x) {
          const Pixel p = fx.get_pixel(x, y);
          for (uint16_t c : {p.r, p.g, p.b}) {
            fold ^= c & 0xFF;
            fold *= 1099511628211ull; // FNV-1a prime
            fold ^= c >> 8;
            fold *= 1099511628211ull;
          }
        }
    }
    hs::clear_mock_time();
    return fold;
  };

  for (const Quad &q : quads)
    HS_EXPECT_EQ(fold_region(false, q), fold_region(true, q));
}

/**
 * @brief White-box accessor for MindSplatter's emitter and hole-kernel state.
 */
struct MindSplatterWhiteBox {
  using MS = MindSplatter<DEFAULT_W, DEFAULT_H>;

  template <typename T>
  using ObjectBytes = std::array<unsigned char, sizeof(T)>;

  template <int W, int H> struct ReplaySnapshot {
    using EffectType = MindSplatter<W, H>;
    using Particle = Animation::Particle<EffectType::TRAIL_LEN>;
    using Attractor = typename EffectType::ParticleSystem::Attractor;
    using Params = typename EffectType::Params;
    using PresetState = decltype(std::declval<EffectType>().presets);

    std::vector<ObjectBytes<Particle>> particles;
    std::vector<ObjectBytes<Attractor>> attractors;
    ObjectBytes<Params> params{};
    ObjectBytes<PresetState> presets{};
    ObjectBytes<Orientation<>> orientation{};
    ObjectBytes<MobiusParams> mobius{};
    std::array<Color4, BakedPalette::LUT_SIZE> palette{};
    std::array<float, EffectType::EmitSolid::NUM_VERTS> emitter_hues{};
    std::array<float, EffectType::EmitSolid::NUM_VERTS> emit_phases{};
    ClipRegion clip{};
    float friction = 0.0f;
    float gravity = 0.0f;
    uint16_t max_life = 0;
  };

  template <typename T> static ObjectBytes<T> object_bytes(const T &source) {
    static_assert(std::is_trivially_copyable_v<T>);
    ObjectBytes<T> bytes;
    std::memcpy(bytes.data(), &source, sizeof(T));
    return bytes;
  }

  template <typename T>
  static void restore_object(T &target, const ObjectBytes<T> &bytes) {
    static_assert(std::is_trivially_copyable_v<T>);
    std::memcpy(&target, bytes.data(), sizeof(T));
  }

  template <int W, int H>
  static ReplaySnapshot<W, H> capture(const MindSplatter<W, H> &ms) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    Snapshot snapshot;
    snapshot.particles.reserve(ms.particle_system.active());
    for (size_t i = 0; i < ms.particle_system.active(); ++i)
      snapshot.particles.push_back(object_bytes(ms.particle_system.pool[i]));
    snapshot.attractors.reserve(ms.particle_system.attractors.size());
    for (size_t i = 0; i < ms.particle_system.attractors.size(); ++i)
      snapshot.attractors.push_back(
          object_bytes(ms.particle_system.attractors[i]));
    snapshot.params = object_bytes(ms.params);
    snapshot.presets = object_bytes(ms.presets);
    snapshot.orientation = object_bytes(ms.orientation);
    snapshot.mobius = object_bytes(ms.mobius);
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
      const float t = static_cast<float>(i) /
                      static_cast<float>(BakedPalette::LUT_SIZE - 1);
      snapshot.palette[i] = ms.baked_palette_.get(t);
    }
    snapshot.emitter_hues = ms.emitter_hues;
    snapshot.emit_phases = ms.emit_phases;
    snapshot.clip = ms.clip();
    snapshot.friction = ms.particle_system.friction;
    snapshot.gravity = ms.particle_system.gravity;
    snapshot.max_life = ms.particle_system.max_life;
    static_assert(std::is_trivially_copyable_v<Particle>);
    return snapshot;
  }

  template <int W, int H>
  static void restore(MindSplatter<W, H> &ms,
                      const ReplaySnapshot<W, H> &snapshot) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    HS_CHECK(ms.particle_system.active() == 0);
    HS_CHECK(snapshot.particles.size() <= ms.particle_system.pool.capacity());
    HS_CHECK(snapshot.attractors.size() ==
             ms.particle_system.attractors.size());

    restore_object(ms.params, snapshot.params);
    restore_object(ms.presets, snapshot.presets);
    restore_object(ms.orientation, snapshot.orientation);
    restore_object(ms.mobius, snapshot.mobius);
    ms.emitter_hues = snapshot.emitter_hues;
    ms.emit_phases = snapshot.emit_phases;
    ms.set_clip(snapshot.clip.y_start, snapshot.clip.y_end,
                snapshot.clip.x_start, snapshot.clip.x_end);
    ms.set_margin(snapshot.clip.margin);
    ms.particle_system.friction = snapshot.friction;
    ms.particle_system.gravity = snapshot.gravity;
    ms.particle_system.max_life = snapshot.max_life;

    for (size_t i = 0; i < snapshot.attractors.size(); ++i)
      restore_object(ms.particle_system.attractors[i], snapshot.attractors[i]);
    for (const ObjectBytes<Particle> &bytes : snapshot.particles) {
      const size_t i = ms.particle_system.active();
      ms.particle_system.spawn(Vector(), Vector(), 0);
      restore_object(ms.particle_system.pool[i], bytes);
    }

    struct PaletteSource {
      const std::array<Color4, BakedPalette::LUT_SIZE> &entries;
      Color4 get(float t) const {
        const float scaled = t * static_cast<float>(BakedPalette::LUT_SIZE - 1);
        const int i = static_cast<int>(scaled + 0.5f);
        return entries[hs::clamp(i, 0, BakedPalette::LUT_SIZE - 1)];
      }
    };
    ms.baked_palette_.rebake(PaletteSource{snapshot.palette});
  }

  template <int W, int H> static void step_physics(MindSplatter<W, H> &ms) {
    Canvas canvas(ms);
    ms.particle_system.step(canvas);
    ms.params.active_count = static_cast<float>(ms.particle_system.active());
  }

  template <int W, int H>
  static bool same_snapshot(const ReplaySnapshot<W, H> &a,
                            const ReplaySnapshot<W, H> &b) {
    if (a.particles != b.particles || a.attractors != b.attractors ||
        a.params != b.params || a.presets != b.presets ||
        a.orientation != b.orientation || a.mobius != b.mobius)
      return false;
    if (std::memcmp(a.emitter_hues.data(), b.emitter_hues.data(),
                    sizeof(a.emitter_hues)) != 0 ||
        std::memcmp(a.emit_phases.data(), b.emit_phases.data(),
                    sizeof(a.emit_phases)) != 0 ||
        std::memcmp(&a.clip, &b.clip, sizeof(a.clip)) != 0 ||
        std::memcmp(&a.friction, &b.friction, sizeof(a.friction)) != 0 ||
        std::memcmp(&a.gravity, &b.gravity, sizeof(a.gravity)) != 0 ||
        a.max_life != b.max_life)
      return false;
    for (size_t i = 0; i < a.palette.size(); ++i) {
      if (a.palette[i].color != b.palette[i].color ||
          std::memcmp(&a.palette[i].alpha, &b.palette[i].alpha,
                      sizeof(a.palette[i].alpha)) != 0)
        return false;
    }
    return true;
  }

  static size_t num_emitters() { return MS::EmitSolid::NUM_VERTS; }
  static float emit_phase(const MS &ms, size_t i) { return ms.emit_phases[i]; }
  static float hue(const MS &ms, size_t i) { return ms.emitter_hues[i]; }
  static float event_horizon() { return MS::EVENT_HORIZON; }
  static float hole_alpha(const Vector &p) {
    return MS::octahedral_hole_alpha(p, fast_cosf(MS::EVENT_HORIZON));
  }
  static float reference_hole_alpha(const Vector &p) {
    const float cos_event_horizon = fast_cosf(MS::EVENT_HORIZON);
    float alpha = 1.0f;
    for (const Vector &axis : MS::AttractSolid::vertices) {
      const float cos_d = dot(p, axis);
      if (cos_d < cos_event_horizon)
        continue;
      const float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
      alpha *= quintic_kernel(d / MS::EVENT_HORIZON);
    }
    return alpha;
  }
  static Vector matrix_vertex(const Vector &v, const MobiusParams &mobius,
                              const Quaternion &orientation) {
    typename MS::RotationMatrix rotation(orientation);
    return rotation.apply(mobius_transform(v, mobius));
  }
  static Vector reference_vertex(const Vector &v, const MobiusParams &mobius,
                                 const Quaternion &orientation) {
    return rotate(mobius_transform(v, mobius), orientation);
  }
  static float normalized_color_seed(uint16_t seed) {
    return MS::normalize_color_seed(seed);
  }
  static float wrapped_color_t(float progress, float seed) {
    return MS::wrap_color_t(progress, seed);
  }
  static constexpr int trail_length() { return MS::TRAIL_LEN; }
  template <int W, int H>
  static void use_reference_orientation(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_orientation = enabled;
  }
  template <int W, int H>
  static void use_reference_color_seed_lookup(MindSplatter<W, H> &ms,
                                              bool enabled) {
    ms.reference_color_seed_lookup = enabled;
  }
  template <int W, int H>
  static void use_reference_vertex_pass(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_vertex_pass = enabled;
  }
  template <int W, int H>
  static void use_reference_hole_kernel(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_hole_kernel = enabled;
  }
  template <int W, int H>
  static void use_reference_palette_alpha(MindSplatter<W, H> &ms,
                                          bool enabled) {
    ms.reference_palette_alpha = enabled;
  }
  template <int W, int H>
  static bool palette_is_opaque(const MindSplatter<W, H> &ms) {
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
      const float t = static_cast<float>(i) /
                      static_cast<float>(BakedPalette::LUT_SIZE - 1);
      if (ms.baked_palette_.get(t).alpha != 1.0f)
        return false;
    }
    return true;
  }
  template <int W, int H>
  static std::array<Pixel, BakedPalette::LUT_SIZE>
  palette_colors(const MindSplatter<W, H> &ms) {
    std::array<Pixel, BakedPalette::LUT_SIZE> colors;
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
      const float t = static_cast<float>(i) /
                      static_cast<float>(BakedPalette::LUT_SIZE - 1);
      colors[i] = ms.baked_palette_.get(t).color;
    }
    return colors;
  }
  template <int W, int H> static void next_preset(MindSplatter<W, H> &ms) {
    ms.presets.next();
  }
  template <int W, int H>
  static size_t preset_index(const MindSplatter<W, H> &ms) {
    return ms.presets.current_index();
  }
  template <int W, int H>
  static void use_reference_signed_axis_physics(MindSplatter<W, H> &ms,
                                                bool enabled) {
    ms.particle_system.reference_signed_axis_physics = enabled;
  }
  template <int W, int H>
  static void use_clip_clear(MindSplatter<W, H> &ms, bool enabled) {
    ms.force_full_buffer_clear = !enabled;
  }
  template <int W, int H>
  static uint16_t active_particles(const MindSplatter<W, H> &ms) {
    return ms.particle_system.active();
  }
  template <int W, int H>
  static void draw_particles(MindSplatter<W, H> &ms, float opacity = 1.0f) {
    Canvas canvas(ms);
    ms.draw_particles(canvas, opacity);
  }
  template <int W, int H>
  static void draw_particles_reference(MindSplatter<W, H> &ms) {
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> sink;
    Canvas canvas(ms);
    ms.draw_particles_with(sink, canvas);
  }
  template <int W, int H>
  static size_t particle_capacity(const MindSplatter<W, H> &ms) {
    return ms.particle_system.pool.capacity();
  }
};

/** @brief Every MindSplatter palette rebuild stays exactly opaque. */
inline void test_mindsplatter_opaque_palette_invariant() {
  using MS = MindSplatter<DEVICE_W, DEVICE_H>;
  using WB = MindSplatterWhiteBox;
  std::array<Pixel, BakedPalette::LUT_SIZE> previous_colors{};
  bool saw_distinct_palette = false;
  bool have_previous = false;
  for (uint8_t seed : {0u, 1u, 63u, 127u, 191u, 255u}) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(seed);
    MS effect;
    effect.init();
    const auto rebuilt_colors = WB::palette_colors(effect);
    if (have_previous && rebuilt_colors != previous_colors)
      saw_distinct_palette = true;
    previous_colors = rebuilt_colors;
    have_previous = true;
    for (int preset = 0; preset < 4; ++preset) {
      HS_EXPECT_EQ(WB::preset_index(effect), static_cast<size_t>(preset));
      HS_EXPECT_TRUE(WB::palette_is_opaque(effect));
      HS_EXPECT_EQ(WB::palette_colors(effect), rebuilt_colors);
      WB::next_preset(effect);
    }
    HS_EXPECT_EQ(WB::preset_index(effect), static_cast<size_t>(0));
  }
  HS_EXPECT_TRUE(saw_distinct_palette);
}

/**
 * @brief Replays one frozen saturated state with exact particle and frame output.
 */
inline void test_mindsplatter_replay_snapshot_exact() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int WARMUP_FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  using Snapshot = WB::ReplaySnapshot<W, H>;

  Snapshot source;
  {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    MS effect;
    effect.init();
    for (int frame = 0; frame < WARMUP_FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
    }
    source = WB::capture(effect);
    HS_EXPECT_EQ(source.particles.size(), WB::particle_capacity(effect));
  }

  struct Replay {
    std::vector<Pixel> framebuffer;
    Snapshot post_step;
    bool restored_exactly = false;
  };
  auto replay = [&]() {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(WARMUP_FRAMES * FRAME_MS, WARMUP_FRAMES * FRAME_US);
    Replay result;
    MS effect;
    effect.init();
    WB::restore(effect, source);
    result.restored_exactly = WB::same_snapshot(source, WB::capture(effect));

    WB::step_physics(effect);
    effect.advance_display();
    result.post_step = WB::capture(effect);
    WB::draw_particles(effect);
    effect.advance_display();
    result.framebuffer.assign(effect.display_buffer(),
                              effect.display_buffer() + W * H);
    return result;
  };

  const Replay first = replay();
  const Replay second = replay();
  hs::clear_mock_time();
  HS_EXPECT_TRUE(first.restored_exactly);
  HS_EXPECT_TRUE(second.restored_exactly);
  HS_EXPECT_TRUE(WB::same_snapshot(first.post_step, second.post_step));
  HS_EXPECT_EQ(first.framebuffer.size(), second.framebuffer.size());
  size_t lit_pixels = 0;
  for (size_t i = 0; i < first.framebuffer.size(); ++i) {
    HS_EXPECT_EQ(first.framebuffer[i], second.framebuffer[i]);
    if (first.framebuffer[i].r | first.framebuffer[i].g |
        first.framebuffer[i].b)
      ++lit_pixels;
  }
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
}

/**
 * @brief The direct AA sink matches the generic pipeline for every quadrant of
 * a frozen saturated MindSplatter particle pool.
 */
inline void test_mindsplatter_saturated_quadrant_sink_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int WARMUP_FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  using Snapshot = WB::ReplaySnapshot<W, H>;

  Snapshot source;
  {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    MS source_effect;
    source_effect.init();
    for (int frame = 0; frame < WARMUP_FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      source_effect.draw_frame();
      source_effect.advance_display();
    }
    source = WB::capture(source_effect);
    HS_EXPECT_EQ(source.particles.size(), WB::particle_capacity(source_effect));
  }

  struct Quadrant {
    int x0, x1, y0, y1;
  };
  constexpr Quadrant quadrants[] = {
      {0, W / 2, 0, H / 2},
      {W / 2, W, 0, H / 2},
      {0, W / 2, H / 2, H},
      {W / 2, W, H / 2, H},
  };

  for (const Quadrant &quadrant : quadrants) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    MS effect;
    effect.init();
    WB::restore(effect, source);
    effect.set_clip(quadrant.y0, quadrant.y1, quadrant.x0, quadrant.x1);

    WB::draw_particles_reference(effect);
    effect.advance_display();
    const std::vector<Pixel> reference(effect.display_buffer(),
                                       effect.display_buffer() + W * H);

    WB::draw_particles(effect);
    effect.advance_display();
    const Pixel *const direct = effect.display_buffer();

    size_t lit_pixels = 0;
    for (int y = quadrant.y0; y < quadrant.y1; ++y) {
      for (int x = quadrant.x0; x < quadrant.x1; ++x) {
        const size_t i = static_cast<size_t>(y) * W + x;
        HS_EXPECT_EQ(direct[i], reference[i]);
        if (reference[i].r | reference[i].g | reference[i].b)
          ++lit_pixels;
      }
    }
    HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
  }
  hs::clear_mock_time();
}

/** @brief Opaque-LUT shading matches the generic path in a saturated quadrant. */
inline void test_mindsplatter_opaque_palette_framebuffer_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int WARMUP_FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  using Snapshot = WB::ReplaySnapshot<W, H>;

  Snapshot source;
  {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    MS effect;
    effect.init();
    for (int frame = 0; frame < WARMUP_FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
    }
    source = WB::capture(effect);
    HS_EXPECT_EQ(source.particles.size(), WB::particle_capacity(effect));
    HS_EXPECT_TRUE(WB::palette_is_opaque(effect));
  }

  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(WARMUP_FRAMES * FRAME_MS, WARMUP_FRAMES * FRAME_US);
    MS effect;
    effect.init();
    WB::restore(effect, source);
    effect.set_clip(0, H / 2, 0, W / 2);
    WB::use_reference_palette_alpha(effect, reference);
    WB::draw_particles(effect, 0.37f);
    effect.advance_display();
    std::vector<Pixel> quadrant;
    quadrant.reserve(static_cast<size_t>(W / 2) * (H / 2));
    for (int y = 0; y < H / 2; ++y)
      for (int x = 0; x < W / 2; ++x)
        quadrant.push_back(effect.get_pixel(x, y));
    return quadrant;
  };

  const std::vector<Pixel> reference = render(true);
  const std::vector<Pixel> opaque = render(false);
  hs::clear_mock_time();
  HS_EXPECT_EQ(reference.size(), opaque.size());
  size_t lit_pixels = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    HS_EXPECT_EQ(reference[i], opaque[i]);
    if (reference[i].r | reference[i].g | reference[i].b)
      ++lit_pixels;
  }
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
}

/** @brief Pins normalized color seeds and wrapped trail progress bit-exactly. */
inline void test_mindsplatter_normalized_color_seed_boundaries() {
  using WB = MindSplatterWhiteBox;
  const float progress_boundaries[] = {
      0.0f, std::nextafter(0.0f, 1.0f), 0.5f, std::nextafter(1.0f, 0.0f), 1.0f,
  };
  size_t samples = 0;
  for (uint32_t seed = 0; seed <= 65535; ++seed) {
    const float reference_seed = static_cast<float>(seed) / 65535.0f;
    const float normalized_seed =
        WB::normalized_color_seed(static_cast<uint16_t>(seed));
    HS_EXPECT_EQ(normalized_seed, reference_seed);
    for (float progress : progress_boundaries) {
      HS_EXPECT_EQ(WB::wrapped_color_t(progress, normalized_seed),
                   wrap_t(progress + reference_seed));
      ++samples;
    }
    const float seam = 1.0f - normalized_seed;
    for (float progress :
         {std::nextafter(seam, 0.0f), seam, std::nextafter(seam, 1.0f)}) {
      HS_EXPECT_EQ(WB::wrapped_color_t(progress, normalized_seed),
                   wrap_t(progress + reference_seed));
      ++samples;
    }
    for (int len = 2; len <= WB::trail_length(); ++len) {
      for (int i = 0; i < len; ++i) {
        const float progress =
            static_cast<float>(i) / static_cast<float>(len - 1);
        HS_EXPECT_EQ(WB::wrapped_color_t(progress, normalized_seed),
                     wrap_t(progress + reference_seed));
        ++samples;
      }
    }
  }
  HS_EXPECT_EQ(WB::normalized_color_seed(0), 0.0f);
  HS_EXPECT_EQ(WB::normalized_color_seed(65535), 1.0f);
  HS_EXPECT_EQ(WB::wrapped_color_t(1.0f, WB::normalized_color_seed(65535)),
               0.0f);
  HS_EXPECT_EQ(samples, static_cast<size_t>(18546688));
}

/**
 * @brief Bounds the precomputed orientation matrix against quaternion rotation.
 */
inline void test_mindsplatter_rotation_matrix_equivalence() {
  using WB = MindSplatterWhiteBox;
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;
  const MobiusParams transforms[] = {
      MobiusParams(), MobiusParams(1, 0, -1.2f, 0, 0, 0, 1, 0),
      MobiusParams(1, 0, -0.6f, 0.6f, 0, 0, 1, 0),
      MobiusParams(0.7f, 0.2f, -0.4f, 0.9f, 0.3f, -0.6f, 1.1f, 0.5f)};
  const Quaternion orientations[] = {
      Quaternion(),
      make_rotation(X_AXIS, PI_F * 0.5f),
      make_rotation(Y_AXIS, PI_F),
      make_rotation(Z_AXIS, -PI_F * 0.75f),
      Quaternion(0.3f, -0.4f, 0.5f, -0.7f).normalized(),
      -Quaternion(0.2f, 0.8f, -0.3f, 0.45f).normalized()};

  struct Tap {
    int x, y;
    uint16_t alpha;
  };
  auto taps = [](const PixelCoords &p) {
    std::array<Tap, 4> result{};
    size_t count = 0;
    Filter::Screen::AntiAlias<W, H> aa;
    aa.plot(p.x, p.y, Pixel(), 0.0f, 1.0f,
            [&](float x, float y, const Pixel &, float, float alpha) {
              result[count++] = {static_cast<int>(x), static_cast<int>(y),
                                 static_cast<uint16_t>(hs::clamp(
                                     alpha * 65535.0f + 0.5f, 0.0f, 65535.0f))};
            });
    return std::pair{result, count};
  };

  float max_component_error = 0.0f;
  float max_angular_error = 0.0f;
  float max_column_error = 0.0f;
  float max_row_error = 0.0f;
  int coverage_differences = 0;
  int q16_differences = 0;
  int max_q16_error = 0;
  size_t sample_count = 0;

  auto check = [&](const Vector &v, const MobiusParams &transform,
                   const Quaternion &orientation) {
    const Vector reference = WB::reference_vertex(v, transform, orientation);
    const Vector matrix = WB::matrix_vertex(v, transform, orientation);
    max_component_error =
        std::max(max_component_error,
                 std::max(std::abs(reference.x - matrix.x),
                          std::max(std::abs(reference.y - matrix.y),
                                   std::abs(reference.z - matrix.z))));
    const float chord = (reference - matrix).length();
    max_angular_error = std::max(
        max_angular_error, 2.0f * asinf(hs::clamp(chord * 0.5f, 0.0f, 1.0f)));

    const PixelCoords reference_pixel = vector_to_pixel<W, H>(reference);
    const PixelCoords matrix_pixel = vector_to_pixel<W, H>(matrix);
    float dx = std::abs(reference_pixel.x - matrix_pixel.x);
    dx = std::min(dx, static_cast<float>(W) - dx);
    max_column_error = std::max(max_column_error, dx);
    max_row_error =
        std::max(max_row_error, std::abs(reference_pixel.y - matrix_pixel.y));
    const auto reference_taps = taps(reference_pixel);
    const auto matrix_taps = taps(matrix_pixel);
    bool same_coverage = reference_taps.second == matrix_taps.second;
    if (same_coverage) {
      for (size_t j = 0; j < reference_taps.second; ++j)
        same_coverage &= reference_taps.first[j].x == matrix_taps.first[j].x &&
                         reference_taps.first[j].y == matrix_taps.first[j].y;
    }
    if (!same_coverage) {
      ++coverage_differences;
    } else {
      for (size_t j = 0; j < reference_taps.second; ++j) {
        const int delta =
            std::abs(static_cast<int>(reference_taps.first[j].alpha) -
                     static_cast<int>(matrix_taps.first[j].alpha));
        if (delta)
          ++q16_differences;
        max_q16_error = std::max(max_q16_error, delta);
      }
    }
    ++sample_count;
  };

  const Vector representative_vectors[] = {
      X_AXIS,
      Y_AXIS,
      Z_AXIS,
      -X_AXIS,
      -Y_AXIS,
      -Z_AXIS,
      Vector(1.0f, 1.0f, 1.0f).normalized(),
      Vector(-1.0f, 1.0f, -1.0f).normalized(),
  };
  hs::random().seed(0x6D617472);
  for (const MobiusParams &transform : transforms) {
    for (const Quaternion &orientation : orientations) {
      for (const Vector &v : representative_vectors)
        check(v, transform, orientation);
      for (int i = 0; i < 20000; ++i) {
        Vector v;
        do {
          const float vx = hs::rand_f(-1.0f, 1.0f);
          const float vy = hs::rand_f(-1.0f, 1.0f);
          const float vz = hs::rand_f(-1.0f, 1.0f);
          v = Vector(vx, vy, vz);
        } while (v.length() < 0.1f);
        v.normalize();
        check(v, transform, orientation);
      }
    }
  }

  std::printf("matrix samples=%zu component=%.9g angle=%.9g dx=%.9g dy=%.9g "
              "coverage=%d q16=%d max_q16=%d\n",
              sample_count, max_component_error, max_angular_error,
              max_column_error, max_row_error, coverage_differences,
              q16_differences, max_q16_error);
  HS_EXPECT_EQ(sample_count, static_cast<size_t>(480192));
  HS_EXPECT_LE(max_component_error, 5e-7f);
  HS_EXPECT_LE(max_angular_error, 5e-7f);
  HS_EXPECT_LE(coverage_differences, 256);
  // Subpixel coverage knife-edge: a ~2-ULP positional drift between the two
  // rotation formulas can flip an anti-alias weight, and the amount depends on
  // host libm rounding (Linux ~264 q16, Windows ~128). Geometric fidelity is
  // guarded by the component/angular bounds above; this only caps the residual.
  HS_EXPECT_LE(max_q16_error, 512);
}

/** @brief Bounds rendered output drift from the matrix orientation path. */
inline void test_mindsplatter_rotation_matrix_framebuffer_error() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    std::vector<Pixel> frames;
    frames.reserve(static_cast<size_t>(W) * H * FRAMES);
    MS effect;
    effect.init();
    WB::use_reference_orientation(effect, reference);
    for (int f = 0; f < FRAMES; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.push_back(effect.get_pixel(x, y));
    }
    return frames;
  };

  const std::vector<Pixel> reference = render(true);
  const std::vector<Pixel> matrix = render(false);
  hs::clear_mock_time();
  size_t different_pixels = 0;
  size_t coverage_differences = 0;
  int max_channel_error = 0;
  uint64_t total_channel_error = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    const Pixel a = reference[i];
    const Pixel b = matrix[i];
    if (a != b)
      ++different_pixels;
    const bool a_black = (a.r | a.g | a.b) == 0;
    const bool b_black = (b.r | b.g | b.b) == 0;
    if (a_black != b_black)
      ++coverage_differences;
    for (int delta :
         {std::abs(static_cast<int>(a.r) - static_cast<int>(b.r)),
          std::abs(static_cast<int>(a.g) - static_cast<int>(b.g)),
          std::abs(static_cast<int>(a.b) - static_cast<int>(b.b))}) {
      max_channel_error = std::max(max_channel_error, delta);
      total_channel_error += static_cast<uint64_t>(delta);
    }
  }
  std::printf("matrix framebuffer samples=%zu different=%zu coverage=%zu "
              "max_channel=%d total_channel=%llu\n",
              reference.size(), different_pixels, coverage_differences,
              max_channel_error,
              static_cast<unsigned long long>(total_channel_error));
  HS_EXPECT_EQ(coverage_differences, static_cast<size_t>(0));
  HS_EXPECT_LE(different_pixels, static_cast<size_t>(64));
  HS_EXPECT_LE(max_channel_error, 1);
  HS_EXPECT_LE(total_channel_error, static_cast<uint64_t>(64));
}

/** @brief The normalized-seed render matches the particle-pool lookup. */
inline void test_mindsplatter_color_seed_framebuffer_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    std::vector<Pixel> frames;
    frames.reserve(static_cast<size_t>(W) * H * FRAMES);
    MS effect;
    effect.init();
    WB::use_reference_color_seed_lookup(effect, reference);
    for (int f = 0; f < FRAMES; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.push_back(effect.get_pixel(x, y));
    }
    return frames;
  };

  const std::vector<Pixel> reference = render(true);
  const std::vector<Pixel> normalized = render(false);
  hs::clear_mock_time();
  size_t lit_pixels = 0;
  size_t different_pixels = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    if (reference[i].r | reference[i].g | reference[i].b)
      ++lit_pixels;
    if (reference[i] != normalized[i])
      ++different_pixels;
  }
  std::printf("color seed framebuffer samples=%zu lit=%zu different=%zu\n",
              reference.size(), lit_pixels, different_pixels);
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(different_pixels, static_cast<size_t>(0));
}

/** @brief Fusing the vertex pass into trail materialization is pixel exact. */
inline void test_mindsplatter_fused_vertex_framebuffer_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    std::vector<Pixel> frames;
    frames.reserve(static_cast<size_t>(W) * H * FRAMES);
    MS effect;
    effect.init();
    WB::use_reference_vertex_pass(effect, reference);
    for (int f = 0; f < FRAMES; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.push_back(effect.get_pixel(x, y));
    }
    return frames;
  };

  const std::vector<Pixel> reference = render(true);
  const std::vector<Pixel> fused = render(false);
  hs::clear_mock_time();
  size_t lit_pixels = 0;
  size_t different_pixels = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    if (reference[i].r | reference[i].g | reference[i].b)
      ++lit_pixels;
    if (reference[i] != fused[i])
      ++different_pixels;
  }
  std::printf("fused vertex framebuffer samples=%zu lit=%zu different=%zu\n",
              reference.size(), lit_pixels, different_pixels);
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(different_pixels, static_cast<size_t>(0));
}

/** @brief The multiply-only hole kernel matches the generic kernel exactly. */
inline void test_mindsplatter_hole_kernel_framebuffer_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    std::vector<Pixel> frames;
    frames.reserve(static_cast<size_t>(W) * H * FRAMES);
    MS effect;
    effect.init();
    WB::use_reference_hole_kernel(effect, reference);
    for (int f = 0; f < FRAMES; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.push_back(effect.get_pixel(x, y));
    }
    return frames;
  };

  const std::vector<Pixel> reference = render(true);
  const std::vector<Pixel> multiply_only = render(false);
  hs::clear_mock_time();
  HS_EXPECT_EQ(reference.size(), multiply_only.size());
  size_t lit_pixels = 0;
  size_t different_pixels = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    if (reference[i].r | reference[i].g | reference[i].b)
      ++lit_pixels;
    if (reference[i] != multiply_only[i])
      ++different_pixels;
  }
  std::printf("hole kernel framebuffer samples=%zu lit=%zu different=%zu\n",
              reference.size(), lit_pixels, different_pixels);
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(different_pixels, static_cast<size_t>(0));
}

/** @brief Clip clearing preserves every pixel displayed by the POV driver. */
inline void test_mindsplatter_clip_clear_display_parity() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int S = H * 2;
  constexpr int N = 4;
  constexpr int FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  struct Render {
    std::vector<Pixel> displayed;
    std::vector<uint16_t> active;
  };

  auto render = [&](int segment_id, bool clip_clear) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    Render result;
    result.displayed.reserve(static_cast<size_t>(W / 2) * (S / N) * FRAMES);
    result.active.reserve(FRAMES);
    MS effect;
    effect.init();
    WB::use_clip_clear(effect, clip_clear);
    const pov::SegmentMap map = pov::segment_map(segment_id, S, N);
    for (int f = 0; f < FRAMES; ++f) {
      const pov::SegmentClip clip =
          pov::segment_clip(map, (f & 1) == 0, S, N, W);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      result.active.push_back(WB::active_particles(effect));
      for (int y = clip.y0; y < clip.y1; ++y)
        for (int x = clip.x0; x < clip.x1; ++x)
          result.displayed.push_back(effect.get_pixel(x, y));
    }
    return result;
  };

  for (int segment_id = 0; segment_id < N; ++segment_id) {
    const Render full_clear = render(segment_id, false);
    const Render clip_clear = render(segment_id, true);
    HS_EXPECT_EQ(full_clear.active.size(), clip_clear.active.size());
    for (size_t i = 0; i < full_clear.active.size(); ++i)
      HS_EXPECT_EQ(full_clear.active[i], clip_clear.active[i]);
    HS_EXPECT_EQ(full_clear.displayed.size(), clip_clear.displayed.size());

    size_t lit_pixels = 0;
    size_t different_pixels = 0;
    size_t coverage_differences = 0;
    for (size_t i = 0; i < full_clear.displayed.size(); ++i) {
      const Pixel a = full_clear.displayed[i];
      const Pixel b = clip_clear.displayed[i];
      const bool a_black = (a.r | a.g | a.b) == 0;
      const bool b_black = (b.r | b.g | b.b) == 0;
      if (!a_black)
        ++lit_pixels;
      if (a != b)
        ++different_pixels;
      if (a_black != b_black)
        ++coverage_differences;
    }
    std::printf("clip clear segment=%d samples=%zu lit=%zu different=%zu "
                "coverage=%zu\n",
                segment_id, full_clear.displayed.size(), lit_pixels,
                different_pixels, coverage_differences);
    HS_EXPECT_EQ(full_clear.displayed.size(),
                 static_cast<size_t>(W / 2) * (S / N) * FRAMES);
    HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
    HS_EXPECT_EQ(different_pixels, static_cast<size_t>(0));
    HS_EXPECT_EQ(coverage_differences, static_cast<size_t>(0));
  }
  hs::clear_mock_time();
}

/** @brief Bounds full-lifetime render drift from signed-axis physics. */
inline void test_mindsplatter_signed_axis_framebuffer_error() {
  constexpr int W = DEVICE_W;
  constexpr int H = DEVICE_H;
  constexpr int FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  struct Render {
    std::vector<Pixel> frames;
    std::vector<uint16_t> active;
  };
  auto render = [&](bool reference) {
    reset_effect_globals();
    GenerativePalette::reset_hue_seed(0);
    hs::set_mock_time(0, 0);
    Render result;
    result.frames.reserve(static_cast<size_t>(W) * H * FRAMES);
    result.active.reserve(FRAMES);
    MS effect;
    effect.init();
    WB::use_reference_signed_axis_physics(effect, reference);
    for (int f = 0; f < FRAMES; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      result.active.push_back(WB::active_particles(effect));
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          result.frames.push_back(effect.get_pixel(x, y));
    }
    return result;
  };

  const Render reference = render(true);
  const Render specialized = render(false);
  hs::clear_mock_time();
  HS_EXPECT_EQ(reference.active.size(), specialized.active.size());
  for (size_t i = 0; i < reference.active.size(); ++i)
    HS_EXPECT_EQ(reference.active[i], specialized.active[i]);

  constexpr int CHECKPOINTS[] = {16, 80, 160};
  constexpr size_t MAX_DIFFERENT[] = {0, 64, 192};
  constexpr int MAX_CHANNEL[] = {0, 8, 512};
  constexpr uint64_t MAX_TOTAL[] = {0, 128, 2048};
  for (size_t checkpoint = 0; checkpoint < 3; ++checkpoint) {
    const int frame = CHECKPOINTS[checkpoint];
    const size_t offset = static_cast<size_t>(frame - 1) * W * H;
    size_t different_pixels = 0;
    size_t coverage_differences = 0;
    int max_channel_error = 0;
    uint64_t total_channel_error = 0;
    for (size_t i = 0; i < static_cast<size_t>(W) * H; ++i) {
      const Pixel a = reference.frames[offset + i];
      const Pixel b = specialized.frames[offset + i];
      if (a != b)
        ++different_pixels;
      const bool a_black = (a.r | a.g | a.b) == 0;
      const bool b_black = (b.r | b.g | b.b) == 0;
      if (a_black != b_black)
        ++coverage_differences;
      for (int delta : {
               std::abs(static_cast<int>(a.r) - static_cast<int>(b.r)),
               std::abs(static_cast<int>(a.g) - static_cast<int>(b.g)),
               std::abs(static_cast<int>(a.b) - static_cast<int>(b.b)),
           }) {
        max_channel_error = std::max(max_channel_error, delta);
        total_channel_error += static_cast<uint64_t>(delta);
      }
    }
    std::printf(
        "axis framebuffer frame=%d active=%u different=%zu coverage=%zu "
        "max_channel=%d total_channel=%llu\n",
        frame, reference.active[frame - 1], different_pixels,
        coverage_differences, max_channel_error,
        static_cast<unsigned long long>(total_channel_error));
    HS_EXPECT_EQ(coverage_differences, static_cast<size_t>(0));
    HS_EXPECT_LE(different_pixels, MAX_DIFFERENT[checkpoint]);
    HS_EXPECT_LE(max_channel_error, MAX_CHANNEL[checkpoint]);
    HS_EXPECT_LE(total_channel_error, MAX_TOTAL[checkpoint]);
  }
}

/**
 * @brief Pins the collapsed signed-axis hole kernel to the six-attractor loop.
 */
inline void test_mindsplatter_octahedral_hole_alpha_equivalence() {
  using WB = MindSplatterWhiteBox;
  constexpr float INV_SNORM16 = 1.0f / 32767.0f;
  auto check = [](const Vector &p) {
    HS_EXPECT_EQ(WB::hole_alpha(p), WB::reference_hole_alpha(p));
  };

  const Vector axes[] = {X_AXIS, -X_AXIS, Y_AXIS, -Y_AXIS, Z_AXIS, -Z_AXIS};
  for (const Vector &axis : axes) {
    Vector tangent = cross(axis, Y_AXIS);
    if (tangent.length() < 0.5f)
      tangent = cross(axis, X_AXIS);
    tangent = tangent.normalized();
    for (int i = 0; i <= 4096; ++i) {
      const float angle = (PI_F * 0.25f) * static_cast<float>(i) / 4096.0f;
      check(axis * cosf(angle) + tangent * sinf(angle));
    }
    const float horizon = WB::event_horizon();
    for (float angle :
         {std::nextafter(horizon, 0.0f), horizon,
          std::nextafter(horizon, std::numeric_limits<float>::infinity())})
      check(axis * cosf(angle) + tangent * sinf(angle));
  }

  hs::random().seed(0x0C7A);
  int outside_unit_sphere = 0;
  for (int i = 0; i < 100000; ++i) {
    Vector p;
    do {
      const float px = hs::rand_f(-1.0f, 1.0f);
      const float py = hs::rand_f(-1.0f, 1.0f);
      const float pz = hs::rand_f(-1.0f, 1.0f);
      p = Vector(px, py, pz);
    } while (p.length() < 0.1f);
    p = p.normalized();
    p = Vector(roundf(p.x * 32767.0f) * INV_SNORM16,
               roundf(p.y * 32767.0f) * INV_SNORM16,
               roundf(p.z * 32767.0f) * INV_SNORM16);
    if (dot(p, p) > 1.0f)
      ++outside_unit_sphere;
    check(p);
  }
  HS_EXPECT_GT(outside_unit_sphere, 1000);
}

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
  static Vector glitch_lens(const Vector &v) {
    return L2::apply_glitch_lens(v);
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
 * @brief Verifies Liquid2D::apply_glitch_lens maps unit directions to unit
 *        directions and returns the pole axis on its near-axis guard branch.
 * @details The lens is a hand-derived degree-3 rational sphere automorphism on
 *          the live per-pixel path; the positive-frame-sum smoke harness cannot
 *          catch a sign/coefficient slip, so pin |lens(v)| == 1 across a spread
 *          of directions plus the R^2 < 1e-6 pole return.
 */
inline void test_liquid2d_glitch_lens_unit_norm() {
  using WB = Liquid2DWhiteBox;
  const Vector dirs[] = {Vector(1, 0, 0),
                         Vector(0, 0, 1),
                         Vector(-1, 0, 0),
                         Vector(0, 0, -1),
                         Vector(1, 1, 1).normalized(),
                         Vector(-1, 2, -3).normalized(),
                         Vector(3, -1, 2).normalized(),
                         Vector(0.2f, -0.9f, 0.4f).normalized(),
                         Vector(-0.7f, 0.1f, 0.7f).normalized()};
  for (const Vector &v : dirs) {
    HS_EXPECT_NEAR(WB::glitch_lens(v).length(), 1.0f, 1e-3f);
  }

  // On-axis input (x = z = 0) trips the pole guard and returns the up vector.
  Vector pole = WB::glitch_lens(Vector(0, -1, 0));
  HS_EXPECT_NEAR(pole.x, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(pole.y, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(pole.z, 0.0f, 1e-6f);
}

/**
 * @brief White-box accessor for MobiusGrid's conformal-radius pole branch and
 *        counter-rotation singularity guard (befriended in effects/MobiusGrid.h).
 */
struct MobiusGridWhiteBox {
  using MG = MobiusGrid<DEFAULT_W, DEFAULT_H>;
  static float conformal_coord(float z, float phase) {
    return MG::conformal_coord(z, phase);
  }
  static Quaternion counter_rotation(const Vector &mid) {
    return MG::counter_rotation(mid);
  }
};

/**
 * @brief Pins MobiusGrid's two hand-derived numeric branches.
 * @details The longitude conformal-radius map R = sqrt((1+z)/(1-z)) is singular
 *          at the poles z = +/-1; the counter-rotation is undefined when the two
 *          Möbius-transformed poles cancel. A plausible frame renders even if
 *          either branch is wrong, so sweep z across the poles (coord must stay
 *          finite in [0,1], poles saturate to 1.0) and check the canceled
 *          midpoint leaves the rotation at identity while a well-defined one maps
 *          onto +Z.
 */
inline void test_mobiusgrid_conformal_and_counter_rotation() {
  using WB = MobiusGridWhiteBox;

  const float zs[] = {-1.0f, -0.999999f, -0.5f, 0.0f, 0.5f, 0.999999f, 1.0f};
  for (float z : zs) {
    for (float phase : {0.0f, 0.25f, 0.6f, 0.95f}) {
      float coord = WB::conformal_coord(z, phase);
      HS_EXPECT_TRUE(std::isfinite(coord));
      HS_EXPECT_GE(coord, 0.0f);
      HS_EXPECT_LE(coord, 1.0f);
    }
  }
  HS_EXPECT_NEAR(WB::conformal_coord(1.0f, 0.3f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::conformal_coord(-1.0f, 0.7f), 1.0f, 1e-6f);

  HS_EXPECT_TRUE(WB::counter_rotation(Vector(0.0f, 0.0f, 0.0f)) ==
                 Quaternion());

  Vector mid(0.3f, -0.7f, 0.4f);
  Vector r = rotate(mid.normalized(), WB::counter_rotation(mid));
  HS_EXPECT_NEAR(r.x, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.y, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.z, 1.0f, 1e-3f);
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
 * @brief Drives ShapeShifter at the Radius slider maximum through a full shape
 *        cycle, covering the antipode-folded (radius > 1) path of every shape.
 * @details The Radius range must stay within the angular-radius domain [0, 2]
 *          (past 2 the fold hands the SDFs a negative radius, which the Star
 *          constructor traps); pin the registered maximum and render four
 *          48-frame cut periods so every shape draws at the extreme.
 */
inline void test_shapeshifter_max_radius_survives_cycle() {
  reset_effect_globals();
  ShapeShifter<DEVICE_W, DEVICE_H> ss;
  ss.init();

  for (const auto &def : ss.getParameters())
    if (std::strcmp(def.name, "Radius") == 0)
      HS_EXPECT_LE(def.max, 2.0f);
  HS_EXPECT_TRUE(ss.updateParameter("Radius", 2.0f));

  for (int f = 0; f < 4 * 48 + 1; ++f) {
    ss.draw_frame();
    ss.advance_display();
  }
}

/**
 * @brief Verifies the segmented-mode full-frame gate on the real effect roster.
 * @details Effect::needs_full_frame() must report true for exactly the effects
 *          whose filter pipeline crosses segment bands (its trait-derived
 *          override surfaces the compile-time any_crosses_segments fold), and
 *          false otherwise so non-stateful effects keep the segment clipping
 *          win. Cross-segment today: MeshFeedback (Pixel::Feedback) and the
 *          World::Trails effect Dynamo. This pins the gate end to
 *          end on constructed effects — the WASM driver reads exactly this query
 *          (targets/wasm/wasm.cpp setClip). See
 *          docs/segmented_stateful_effects_spec.md.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer (single-live guard) and
  // reconfigures the shared arenas/timeline in init(), so construct one at a
  // time with the same fresh setup smoke_one uses.
  auto reset = [] { reset_effect_globals(); };
  auto check = [&](Effect &fx, bool expected, const char *name) {
    fx.init();
    if (fx.needs_full_frame() != expected)
      std::printf("  [FAIL] %s needs_full_frame()=%d expected=%d\n", name,
                  fx.needs_full_frame(), expected);
    HS_EXPECT_EQ(fx.needs_full_frame(), expected);
  };

  // Cross-segment stateful effects -> full-frame.
  {
    reset();
    MeshFeedback<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "MeshFeedback");
  }
  {
    reset();
    Dynamo<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "Dynamo");
  }

  // Representative non-stateful effects -> keep band clipping (default false).
  {
    reset();
    Voronoi<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "Voronoi");
  }
  {
    reset();
    RingSpin<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "RingSpin");
  }
}

/**
 * @brief White-box accessor for Voronoi's seeded sites and coarse-grid tuning
 *        constants. MAX_SITES and COHERENCE_BLOCK do not depend on W/H, so one
 *        instantiation supplies them; the block floor does, and is reached
 *        through coherence_block_min<W, H>().
 */
struct VoronoiWhiteBox {
  using VO = Voronoi<DEFAULT_W, DEFAULT_H>;
  static constexpr int MAX_SITES = VO::MAX_SITES;
  static constexpr int COHERENCE_BLOCK = VO::COHERENCE_BLOCK;

  /** @brief Adaptive block floor at render resolution W x H. */
  template <int W, int H> static constexpr int coherence_block_min() {
    return Voronoi<W, H>::COHERENCE_BLOCK_MIN;
  }

  /** @brief Number of currently seeded sites. */
  template <int W, int H> static size_t site_count(const Voronoi<W, H> &v) {
    return v.sites_buffer.size();
  }
  /** @brief Spin axis of seeded site @p i. */
  template <int W, int H>
  static Vector site_axis(const Voronoi<W, H> &v, size_t i) {
    return v.sites_buffer[i].axis;
  }
};

/**
 * @brief Verifies Voronoi seeds every spin axis through random_vector().
 */
inline void test_voronoi_axes_use_uniform_sampler() {
  using WB = VoronoiWhiteBox;
  reset_effect_globals();
  Voronoi<DEVICE_W, DEVICE_H> effect;
  effect.init();

  hs::random().seed(1337u);
  for (size_t i = 0; i < WB::site_count(effect); ++i)
    HS_EXPECT_VEC(WB::site_axis(effect, i), random_vector(), 0.0f);
}

/**
 * @brief Classifies the coarse-coherence corner grid at the effect's adaptive
 *        block size and measures block-candidate-union shading against a full
 *        k=1 query at every pixel.
 * @tparam W,H Render resolution.
 * @param sites Site positions on the unit sphere.
 * @param max_deficit Out: worst dot(p, true nearest) - dot(p, union nearest)
 *        over the mismatched pixels (0 when every pixel matches).
 * @return Fraction of pixels whose union-of-corner-pairs nearest matches the
 *         true nearest site.
 * @details Mirrors Voronoi::draw_frame's candidate path: the block edge B and
 *          the per-corner k=2 pair derivation match production, and each pixel
 *          takes an exact top-1 over the deduped union of its block's four
 *          corner pairs.
 */
template <int W, int H>
inline double voronoi_union_nearest_match(std::span<const Vector> sites,
                                          float &max_deficit) {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  KDTree tree(arena, sites);

  struct Pair {
    uint16_t lo, hi;
  };
  auto classify = [&](const Vector &p) -> Pair {
    auto knn = tree.nearest(p, 2);
    uint16_t a = knn[0].original_index;
    uint16_t b = knn.size() > 1 ? knn[1].original_index : a;
    return {std::min(a, b), std::max(a, b)};
  };

  // Mirror of Voronoi::draw_frame's block-edge derivation (effects/Voronoi.h).
  const float cell_px =
      (2.0f * H / PI_F) / sqrtf(static_cast<float>(sites.size()));
  const int B = std::clamp(static_cast<int>(cell_px),
                           VoronoiWhiteBox::coherence_block_min<W, H>(),
                           VoronoiWhiteBox::COHERENCE_BLOCK);

  const int nbx = (W - 1) / B + 2;
  const int nby = (H - 1) / B + 2;
  Pair *corners = static_cast<Pair *>(
      arena.allocate(size_t(nbx) * nby * sizeof(Pair), alignof(Pair)));
  auto cx = [&](int j) { return std::min(j * B, W - 1); };
  auto cy = [&](int k) { return std::min(k * B, H - 1); };
  for (int k = 0; k < nby; ++k)
    for (int j = 0; j < nbx; ++j)
      corners[k * nbx + j] = classify(pixel_to_vector<W, H>(cx(j), cy(k)));

  long matched = 0;
  max_deficit = 0.0f;
  for (int y = 0; y < H; ++y) {
    const int ky = y / B;
    for (int x = 0; x < W; ++x) {
      const int jx = x / B;
      uint16_t u[8];
      int un = 0;
      auto add = [&](uint16_t s) {
        for (int i = 0; i < un; ++i)
          if (u[i] == s)
            return;
        u[un++] = s;
      };
      for (const Pair *c :
           {&corners[ky * nbx + jx], &corners[ky * nbx + jx + 1],
            &corners[(ky + 1) * nbx + jx], &corners[(ky + 1) * nbx + jx + 1]}) {
        add(c->lo);
        add(c->hi);
      }
      const Vector p = pixel_to_vector<W, H>(x, y);
      float best = -2.0f;
      uint16_t bi = 0;
      for (int i = 0; i < un; ++i) {
        const float d = dot(p, sites[u[i]]);
        if (d > best) {
          best = d;
          bi = u[i];
        }
      }
      auto knn = tree.nearest(p, 1);
      // On an exact tie the union scan keeps the first candidate and the tree
      // keeps the other, so compare distance rather than site index.
      const float deficit = dot(p, knn[0].point) - best;
      if (knn[0].original_index == bi || deficit <= 0.0f)
        ++matched;
      else
        max_deficit = std::max(max_deficit, deficit);
    }
  }
  return static_cast<double>(matched) / (static_cast<double>(W) * H);
}

/**
 * @brief Pins Voronoi's block-candidate-union coverage across the adaptive
 *        block regime at both render resolutions: the union of a block's four
 *        corner pairs contains the true nearest site at every pixel in the
 *        low-density octahedral case, and at >= 99.9% of pixels (with a
 *        sub-visibility dot deficit on the rest) at the MAX_SITES Fibonacci
 *        spread.
 * @details The dense cases floor the adaptive block at COHERENCE_BLOCK_MIN. A
 *          block edge that outruns the cell pixel extent straddles whole cells
 *          and collapses the match fraction, so the short-canvas resolution —
 *          where MAX_SITES cells are sub-pixel and the floor drops to 1 — is
 *          checked for exact coverage.
 */
inline void test_voronoi_union_candidates_cover_nearest() {
  float deficit = 0.0f;

  const Vector octahedral[] = {
      Vector(1, 0, 0),  Vector(-1, 0, 0), Vector(0, 1, 0),
      Vector(0, -1, 0), Vector(0, 0, 1),  Vector(0, 0, -1),
  };
  const size_t octa_count = sizeof(octahedral) / sizeof(octahedral[0]);
  const std::span<const Vector> sparse(octahedral, octa_count);
  const double octa_match =
      voronoi_union_nearest_match<DEFAULT_W, DEFAULT_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match, 1.0);
  const double octa_match_dev =
      voronoi_union_nearest_match<DEVICE_W, DEVICE_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match_dev, 1.0);

  // Dense regime: seed MAX_SITES on a Fibonacci sphere exactly as
  // Voronoi::seed_sites places them, so the adaptive block floors at
  // COHERENCE_BLOCK_MIN.
  constexpr int N = VoronoiWhiteBox::MAX_SITES;
  static Vector fib[N];
  const float golden_angle = PI_F * (3.0f - sqrtf(5.0f));
  for (int i = 0; i < N; ++i) {
    const int span = N > 1 ? N - 1 : 1;
    const float y = 1.0f - (i / static_cast<float>(span)) * 2.0f;
    const float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
    const float theta = golden_angle * i;
    fib[i] = Vector(cosf(theta) * radius, y, sinf(theta) * radius);
  }
  const std::span<const Vector> dense(fib, N);
  const double fib_match =
      voronoi_union_nearest_match<DEFAULT_W, DEFAULT_H>(dense, deficit);
  HS_EXPECT_GE(fib_match, 0.999);
  HS_EXPECT_LE(deficit, 0.005f);

  const double fib_match_dev =
      voronoi_union_nearest_match<DEVICE_W, DEVICE_H>(dense, deficit);
  HS_EXPECT_EQ(fib_match_dev, 1.0);
  HS_EXPECT_EQ(deficit, 0.0f);
}

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

    // The effect's held graph-walk seed; dodecahedron is the largest Platonic.
    PolyMesh seed;
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed =
          Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), target);
    });

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
      MeshOps::classify_faces_by_topology(mesh, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Render peak: transform into scratch_a, then Scan::Mesh::draw stacks a
    // FaceScratchBuffer on top (the scratch_a-binding path per init's comment).
    {
      ScratchScope a_guard(scratch_arena_a);
      Orientation<> orientation;
      OrientTransformer camera(orientation);
      MeshState rotated;
      MeshOps::transform(mesh, rotated, scratch_arena_a, camera);
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      Pipeline<W, H> filters;
      auto frag = [](const Vector &, Fragment &f) {
        f.color = Color4(Pixel(1000, 1000, 1000), 1.0f);
      };
      Scan::Mesh::draw<W, H>(filters, canvas, rotated, frag, scratch_arena_a);
    }

    // Morph compaction peak: the CompiledHankin + palette bank survive into
    // scratch_b, the mesh + walk seed into scratch_a, then persistent is
    // reset — the same Persist discipline finish_morph_cycle uses to compact
    // between legs.
    {
      Persist<CompiledHankin> ph(hankin, scratch_arena_b, persistent_arena);
      Persist<MeshState> pf(mesh, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank, scratch_arena_b,
                                  persistent_arena);
      Persist<PolyMesh> ps(seed, scratch_arena_a, persistent_arena);
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
 * @brief White-box accessor for IslamicStars' private build-chain state
 *        (befriended in effects/IslamicStars.h).
 * @details The op-by-op recipe build only runs when the round-robin reaches a
 * non-null recipe entry (index 1), ~340 frames into a default-speed run — past
 * every generic smoke window. The probe pre-sets Trans Speed before init so
 * the whole crossing fits in ~100 frames, and reads build_active_/solid_idx
 * to pin that the build ran and completed. <96,20> keeps the raster cheap;
 * the build bookkeeping is resolution-independent.
 */
struct IslamicBuildProbe {
  using IS = IslamicStars<DEVICE_W, DEVICE_H>;
  template <int W, int H>
  static void set_trans_speed(IslamicStars<W, H> &e, float v) {
    e.params.trans_speed = v;
  }
  template <int W, int H>
  static bool build_active(const IslamicStars<W, H> &e) {
    return e.build_active_;
  }
  template <int W, int H> static int solid_idx(const IslamicStars<W, H> &e) {
    return e.solid_idx;
  }
  template <int W, int H> static int dual_bridges(const IslamicStars<W, H> &e) {
    return e.dual_bridges_built_;
  }
  template <int W, int H> static int gated_swaps(const IslamicStars<W, H> &e) {
    return e.gated_swaps_scheduled_;
  }
  template <int W, int H>
  static size_t persistent_budget(const IslamicStars<W, H> &e) {
    return e.device_persistent_budget_;
  }
  template <int W, int H> static int front_slot(IslamicStars<W, H> &e) {
    return e.carousel.front_index();
  }
  template <int W, int H>
  static const uint8_t *slot_palette(const IslamicStars<W, H> &e, int slot) {
    return e.slot_face_palette[slot];
  }
  template <int W, int H>
  static size_t slot_faces(IslamicStars<W, H> &e, int slot) {
    return e.carousel.slot(slot).topology.size();
  }
  static constexpr size_t bridge_scratch_a() {
    return IS::SPLIT_SCRATCH_A_BRIDGE;
  }
  static constexpr size_t default_scratch_a() {
    return IS::SPLIT_SCRATCH_A_DEFAULT;
  }
  static constexpr size_t bridge_scratch_b() {
    return IS::SPLIT_SCRATCH_B_BRIDGE;
  }
  static constexpr int sprite_fade_frames() { return IS::SPRITE_FADE_FRAMES; }
};

/**
 * @brief Verifies the first recipe seed uses the Sprite envelope's 16-frame
 *        fade-in and starts its build at the full-opacity boundary.
 */
inline void test_islamicstars_seed_sprite_fade_in() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();

  constexpr int FADE_FRAMES = IslamicBuildProbe::sprite_fade_frames();
  HS_EXPECT_EQ(FADE_FRAMES, 16);
  for (int frame = 1; frame < FADE_FRAMES; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
  }

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(IslamicBuildProbe::build_active(effect));
}

/**
 * @brief Drives IslamicStars across the second registry entry's complete
 *        op-by-op build at max trans speed: the build must activate and
 *        finish without a trap, the built shape's per-face colours must never
 *        change from finish_build through its still/ripple/fade display, and
 *        the shape after it must start cleanly and light pixels.
 */
inline void test_islamicstars_recipe_build_smoke() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  effect.init();

  // Entry 1 (dodecahedron_hk62_ambo_hk62) is the second spawn; the third
  // spawn requires its build to have completed.
  constexpr int MAX_FRAMES = 400;
  int frames = 0;
  int build_frames = 0;
  bool was_building = false;
  int built_slot = -1;
  int snap_solid = -1;
  int changed_after_build = 0;
  int constant_frames = 0;
  std::vector<uint8_t> built_pal;
  while (frames < MAX_FRAMES && IslamicBuildProbe::solid_idx(effect) < 2) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    const bool building = IslamicBuildProbe::build_active(effect);
    if (building)
      ++build_frames;
    // Snapshot the built shape's per-face colours the moment its build
    // completes; they must stay byte-identical through its still, ripple, and
    // fade phases (the next spawn retires the shape and may reuse its array).
    if (was_building && !building && built_slot < 0) {
      built_slot = IslamicBuildProbe::front_slot(effect);
      snap_solid = IslamicBuildProbe::solid_idx(effect);
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      built_pal.assign(pal,
                       pal + IslamicBuildProbe::slot_faces(effect, built_slot));
    } else if (built_slot >= 0 &&
               IslamicBuildProbe::solid_idx(effect) == snap_solid) {
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      for (size_t f = 0; f < built_pal.size(); ++f)
        if (pal[f] != built_pal[f])
          ++changed_after_build;
      ++constant_frames;
    }
    was_building = building;
  }
  HS_EXPECT_LT(frames, MAX_FRAMES);
  HS_EXPECT_GT(build_frames, 0);
  HS_EXPECT_TRUE(!IslamicBuildProbe::build_active(effect));
  HS_EXPECT_GE(built_slot, 0);
  HS_EXPECT_GT(built_pal.size(), (size_t)0);
  HS_EXPECT_GT(constant_frames, 0);
  HS_EXPECT_EQ(changed_after_build, 0);

  // The following shape renders lit frames.
  for (int f = 0; f < 12; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  uint64_t acc = 0;
  for (int y = 0; y < DEVICE_H; ++y)
    for (int x = 0; x < DEVICE_W; ++x) {
      const Pixel &p = effect.get_pixel(x, y);
      acc += static_cast<uint64_t>(p.r) + p.g + p.b;
    }
  HS_EXPECT_GT(acc, (uint64_t)0);
}

/**
 * @brief Drives IslamicStars through every registry entry, pinning the
 *        persistent arena against the effect's own budget.
 * @details The per-chain gate measures a build in isolation, which misses what
 *          the effect actually holds: a build follows whatever shape preceded
 *          it, and the roster's heaviest entry is 1082 faces. Cycling the whole
 *          roster is the only thing that exercises a build against that
 *          predecessor. Rendered small — the peak is mesh-driven, not
 *          canvas-driven — because the gate is about memory, not pixels. An
 *          Arena overrun traps, so an OOM fails this test by killing the run.
 */
inline void test_islamicstars_roster_cycle_fits_budget() {
  reset_effect_globals();
  // Production resolution, not device: the swept+compiled+draw scratch peak a
  // build leg reaches is mesh-driven and resolution-independent, but the effect
  // ships at 288x144, so the gate holds the real geometry. Trans Speed 8
  // compresses each stage so the whole roster cycles in ~1200 frames without
  // dropping any leg (never lower the resolution: that would shrink the peak).
  IslamicStars<288, 144> effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  effect.init();

  // Per-shape arena split (IslamicStars::spawn_shape): smooth kis/needle
  // bridge shapes use 129.5 KB / 74 KB scratch; ordinary recipe builds use
  // 116 KB / 72 KB; non-recipe generation retains 116 KB / 74 KB. On host the
  // scratch arenas are the exact device sizes and hard-capped, so an over-budget
  // leg traps -- completing the whole roster proves every shape's scratch fit
  // its own split. Persistent is host-inflated (soft), so it is checked per
  // frame against the effect's live per-shape device budget; host pointer
  // inflation makes the resident offset an upper bound on the device figure.
  auto solids = Solids::Collections::get_islamic_solids();
  const int entries = static_cast<int>(solids.size());
  int needle_idx = -1;
  for (int i = 0; i < entries; ++i)
    if (std::strstr(solids[i].name, "needle"))
      needle_idx = i;
  HS_EXPECT_GE(needle_idx, 0);

  constexpr int MAX_FRAMES = 20000;
  size_t a_peak = 0, b_peak = 0, persist_peak = 0;
  size_t na_peak = 0, nb_peak = 0, np_peak = 0; // needle-only peaks
  size_t worst_p = 0, worst_p_budget = 0;
  int worst_p_idx = -1; // shape at the worst persistent/budget ratio
  int frames = 0, shapes = 0, builds = 0;
  bool was_building = false;
  int last = IslamicBuildProbe::solid_idx(effect);
  while (frames < MAX_FRAMES && shapes <= entries) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    const int cur = IslamicBuildProbe::solid_idx(effect);
    // Palette variety of each finished build: distinct landed palettes on the
    // shape the real leg chain just landed.
    const bool building = IslamicBuildProbe::build_active(effect);
    if (was_building && !building) {
      const int front = IslamicBuildProbe::front_slot(effect);
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, front);
      const size_t nf = IslamicBuildProbe::slot_faces(effect, front);
      bool seen[MeshPaletteBank::N] = {};
      int distinct = 0;
      for (size_t f = 0; f < nf; ++f)
        if (pal[f] < MeshPaletteBank::N && !seen[pal[f]]) {
          seen[pal[f]] = true;
          ++distinct;
        }
      std::printf("  [built] %s: %d/%d palettes on %zu faces\n",
                  (cur >= 0 && cur < entries) ? solids[cur].name : "?",
                  distinct, MeshPaletteBank::N, nf);
    }
    was_building = building;
    const size_t p = persistent_arena.get_offset();
    const size_t a = scratch_arena_a.get_high_water_mark();
    const size_t b = scratch_arena_b.get_high_water_mark();
    const size_t p_budget = IslamicBuildProbe::persistent_budget(effect);
    a_peak = std::max(a_peak, a);
    b_peak = std::max(b_peak, b);
    persist_peak = std::max(persist_peak, p);
    if (p > worst_p) { // report the tightest persistent fit, not the raw max
      worst_p = p;
      worst_p_budget = p_budget;
      worst_p_idx = cur;
    }
    if (cur == needle_idx) {
      na_peak = std::max(na_peak, a);
      nb_peak = std::max(nb_peak, b);
      np_peak = std::max(np_peak, p);
    }
    // Per-shape persistent budget (device figure); scratch is trap-enforced.
    HS_EXPECT_LE(p, p_budget);
    if (building)
      ++builds;
    if (cur != last) {
      last = cur;
      ++shapes;
    }
  }

  const char *worst_name = (worst_p_idx >= 0 && worst_p_idx < entries)
                               ? solids[worst_p_idx].name
                               : "?";
  std::printf(
      "  [roster] %d shapes over %d frames, %d build frames: scratch_a "
      "peak=%zu B, scratch_b peak=%zu B, persistent peak=%zu B; tightest "
      "persistent %zu/%zu at %s\n",
      shapes, frames, builds, a_peak, b_peak, persist_peak, worst_p,
      worst_p_budget, worst_name);
  std::printf("  [roster] needle smooth-path peaks: scratch_a=%zu/%zu "
              "scratch_b=%zu/%zu persistent=%zu/%zu; gated_swaps=%d\n",
              na_peak, IslamicBuildProbe::bridge_scratch_a(), nb_peak,
              IslamicBuildProbe::bridge_scratch_b(), np_peak,
              DEVICE_GLOBAL_ARENA_SIZE - IslamicBuildProbe::bridge_scratch_a() -
                  IslamicBuildProbe::bridge_scratch_b(),
              IslamicBuildProbe::gated_swaps(effect));
  HS_EXPECT_GT(shapes, entries - 1);
  HS_EXPECT_GT(builds, 0);
  // Every dt/dtd chain takes the smooth path -- no gated (snapping) kis.
  HS_EXPECT_EQ(IslamicBuildProbe::gated_swaps(effect), 0);
  // needle actually reached its scratch_a-heavy split (proves the smooth path
  // ran, not a silently-dropped build).
  HS_EXPECT_GT(na_peak, 120u * 1024u);
}

/**
 * @brief Drives IslamicStars until every DUAL-bearing recipe has built its
 *        smooth three-leg dual bridge, pinning the scratch peaks against the
 *        effect's budget.
 * @details The roster gate at Trans Speed 8 compresses each build so far that a
 *          heavy shape's closing dual leg can be dropped before it runs; this
 *          drives at a modest speed so every bridge leg runs in full. Five
 *          registry recipes carry a DUAL (a needle, two gyros, a kis-gyro, and
 *          a truncate50d_ambo_dual). The bridge's leg 3 rebuilds the medial for
 *          its handoff centroids, whose scratch must not co-reside with the
 *          leg's own arrival mesh -- an over-budget leg traps in the host arena.
 */
inline void test_islamicstars_dual_bridge_fits_budget() {
  reset_effect_globals();
  IslamicStars<288, 144> effect;
  IslamicBuildProbe::set_trans_speed(effect, 2.0f);
  effect.init();

  constexpr int TARGET_BRIDGES = 5;
  constexpr int MAX_FRAMES = 40000;
  size_t a_peak = 0, b_peak = 0, persist_peak = 0;
  int frames = 0;
  // Scratch is hard-capped per-shape (spawn_shape's resplit), so a leg over its
  // split traps here; completing without a trap proves every full bridge fit.
  // Persistent is host-inflated, so it is checked per frame against the effect's
  // live per-shape device budget.
  while (frames < MAX_FRAMES &&
         IslamicBuildProbe::dual_bridges(effect) < TARGET_BRIDGES) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
    b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
    const size_t p = persistent_arena.get_offset();
    persist_peak = std::max(persist_peak, p);
    HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
  }
  std::printf(
      "  [dual-bridge] %d bridges over %d frames: scratch_a peak=%zu B, "
      "scratch_b peak=%zu B, persistent peak=%zu B; gated_swaps=%d\n",
      IslamicBuildProbe::dual_bridges(effect), frames, a_peak, b_peak,
      persist_peak, IslamicBuildProbe::gated_swaps(effect));
  HS_EXPECT_GE(IslamicBuildProbe::dual_bridges(effect), TARGET_BRIDGES);
  // Every full bridge runs the smooth path -- no gated (snapping) kis.
  HS_EXPECT_EQ(IslamicBuildProbe::gated_swaps(effect), 0);
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
  test_mindsplatter_opaque_palette_invariant();
  test_mindsplatter_replay_snapshot_exact();
  test_mindsplatter_saturated_quadrant_sink_parity();
  test_mindsplatter_opaque_palette_framebuffer_parity();
  test_mindsplatter_octahedral_hole_alpha_equivalence();
  test_mindsplatter_normalized_color_seed_boundaries();
  test_mindsplatter_rotation_matrix_equivalence();
  test_mindsplatter_rotation_matrix_framebuffer_error();
  test_mindsplatter_color_seed_framebuffer_parity();
  test_mindsplatter_fused_vertex_framebuffer_parity();
  test_mindsplatter_hole_kernel_framebuffer_parity();
  test_mindsplatter_clip_clear_display_parity();
  test_mindsplatter_signed_axis_framebuffer_error();

  // FULL tier only (HS_EFFECTS_FULL=1; CI on every push/PR). The white-box
  // correctness block and the 288x144 production-resolution roster passes below
  // are the ~40 s bulk of this module — dominated by full-frame software raster
  // (see effects_full_suite()). The QUICK tier (default, local pre-commit) skips
  // straight to the ~2 s device-resolution passes, which already cover every
  // effect's construct/init/render/read-back and cross-run determinism. So a
  // green local commit is NOT authoritative for the full-resolution paths or the
  // white-box invariants — CI is (same split as HS_SMOKE_FRAMES's cyclic-window
  // coverage; see the pre-commit hook header).
  if (effects_full_suite()) {
    test_needs_full_frame_gate();
    test_voronoi_axes_use_uniform_sampler();
    test_voronoi_union_candidates_cover_nearest();
    test_sh_decode_lm_valid_order();
    test_sh_cartesian_matches_spherical();
    test_gs_q16_roundtrip();
    test_gs_rest_state_is_fixed_point();
    test_gs_substep_signs_and_clamp();
    test_gs_evolution_stays_bounded();
    test_gs_reaction_corner_stays_bounded();
    test_gs_odd_substep_lands_in_persistent();
    test_gs_dissolve_clears_and_reseeds();
    test_gs_reaction_edit_starts_dissolve();
    test_bz_q16_roundtrip();
    test_bz_min_diffusion_step_survives_quantization();
    test_bz_advance_species_signs_and_clamp();
    test_bz_perturb_state_saturates_and_nudges();
    test_bz_perturb_state_draw_count_pinned();
    test_bz_substep_diffuses();
    test_dreamballs_preset_cycle_bookkeeping();
    test_dreamballs_respawn_fires_and_honors_pause();
    CometsWhiteBox::check_paths_close();
    ThrustersWhiteBox::check_warp_endpoints();
    RingShowerWhiteBox::check_radius_endpoints();
    DynamoWhiteBox::check_overlapping_wipes_stay_in_range();
    test_hopf_projection_math();
    test_petalflow_spawn_gap_bounded();
    test_displacement_field_lazy_hue_table_matches_eager();
    test_displacement_field_hue_table_fidelity();
    test_displacement_field_hue_table_frame_fidelity();
    test_displacement_field_zero_hue_scale_is_exact();
    test_displacement_field_clip_tiles_full();
    test_mindsplatter_emit_phase_wrapped();
    test_flyby_phase_wrapped();
    test_liquid2d_phase_wrapped();
    test_liquid2d_glitch_lens_unit_norm();
    test_mobiusgrid_conformal_and_counter_rotation();
    test_ringspin_pool_clamped();
    test_shapeshifter_shape_cut_lifecycle();
    test_shapeshifter_max_radius_survives_cycle();
    test_hankinsolids_arena_budget_covers_every_solid();
    test_islamicstars_seed_sprite_fade_in();
    test_islamicstars_recipe_build_smoke();
    test_islamicstars_roster_cycle_fits_budget();
    test_islamicstars_dual_bridge_fits_budget();

    // Full production-resolution roster passes (288x144): smoke, then cross-run
    // determinism under the injected clock.
    g_nonblack_effects = 0;
#define HS_SMOKE_ONE(name) smoke_one<name>(#name);
    HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE
    // At least one effect must light up, catching a total regression-to-black.
    HS_EXPECT_GT(g_nonblack_effects, 0);
#define HS_DET_ONE(name) determinism_one<name>(#name);
    HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE
  }

  // Device-resolution <96,20> roster passes — always run; this is the QUICK
  // tier's core and the only place that specialization runs under native asserts
  // (see DEVICE_W/DEVICE_H).
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

  // Every effect clip-clears unless it declares needs_full_frame, so the roster
  // is swept rather than the one effect the optimization started on.
  std::printf("  -- clip-clear display parity --\n");
#define HS_CLIP_PARITY_ONE(name)                                               \
  clip_clear_parity_one<name, DEVICE_W, DEVICE_H>(#name);
  HS_EFFECT_LIST(HS_CLIP_PARITY_ONE)
#undef HS_CLIP_PARITY_ONE

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test
