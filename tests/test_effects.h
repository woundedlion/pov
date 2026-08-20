/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Effect smoke / robustness harness — exercises all registered effects.
 *
 * Defines the per-effect sweep primitives:
 *   1. smoke_one  — construct/init/render/read-back, under native asserts.
 *   2. determinism_one — renders each effect twice under an injected, fixed
 *      per-frame clock (hs::set_mock_time) and asserts the two final frames are
 *      byte-identical.
 *   3. clip_clear_parity_one — renders under a moving segment clip with each
 *      clear scope and requires the displayed pixels to agree.
 * The roster sweeps that apply them are the effects_smoke module
 * (tests/test_effects_smoke.h); this module holds the white-box invariants.
 */
#pragma once

#include "core/engine/effects.h"
#include "core/render/canvas.h"
#include "core/render/sdf/volume.h"
#include "core/engine/memory.h"
#include "hardware/pov_segment_map.h"
#include "tests/pixel_test_util.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/mindsplatter_whitebox.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <vector>
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
 * @brief Small-aspect render width in pixels.
 * @details <96,20> is a simulator/test resolution, not a hardware target: every
 * PlatformIO env builds 288x144. The native suite is the only place the <96,20>
 * specialization runs under asserts (the device forces NDEBUG; the CI WASM
 * smoke runs assert-free), so a second roster pass at this resolution exercises
 * height-20-specific paths (PhiLUT<20> indexing, small-aspect arena sizing,
 * H_OFFSET interactions) that bypass every assert-enabled layer otherwise. Both
 * are in HS_RESOLUTIONS (core/control/registry.h).
 */
constexpr int SMALL_W = 96;
/**
 * @brief Small-aspect render height in pixels.
 * @details Paired with SMALL_W for the <96,20> specialization.
 */
constexpr int SMALL_H = 20;

/**
 * @brief Per-effect smoke frame count, resolved from HS_SMOKE_FRAMES.
 * @details Shared with the arena/stack budget gates (tests/test_fixture.h), so
 * they measure over the window this roster sweep renders. Frame count is a
 * minor cost lever here: the module is dominated by full-resolution software
 * raster (~71 ms/frame at 288x144, once per HS_EFFECT_LIST entry), which the
 * QUICK tier skips entirely (see effects_full_suite()).
 */
using hs_test::smoke_frames;

/**
 * @brief Selects the effects test depth tier from the environment.
 * @return true for the FULL suite (white-box correctness block + the
 * production-resolution 288x144 roster passes), false for the QUICK tier.
 * @details The full suite renders every roster effect at the 288x144 production
 * resolution across a smoke pass, a determinism pass, and several full-frame
 * white-box tests — dominated by the ~71 ms/frame software raster of
 * 41,472-pixel frames. The QUICK tier (default) runs only the small-aspect
 * <96,20> smoke + determinism passes, ~1,920-pixel frames that cover every
 * effect's construct/init/render/read-back and cross-run determinism in ~2 s —
 * the fast path for the local pre-commit hook. CI opts into the full suite on
 * every push/PR by setting HS_EFFECTS_FULL=1 (.github/workflows/ci.yml), so the
 * full-resolution passes and the white-box correctness block are the
 * authoritative gate there, not locally. Set HS_EFFECTS_FULL=1 to reproduce the
 * CI depth in a local commit. Read by both the effects and effects_smoke
 * modules, which split those two halves.
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

/**
 * @brief Forward declaration of the animated-param pause lint.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in diagnostic output.
 */
inline void lint_animated_pause(Effect &effect, const char *name);

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
 * @brief Resets the process-global effect state to a clean per-effect baseline.
 * @details Forwards to hs_test::reset_globals(). Every effect aliases the same
 *          static double buffer, global RNG, arenas, and timeline; without this
 *          reset leftover events reference the previous (destroyed) effect
 *          instance. The mock clock is released too, so a case that wants one
 *          pins it itself rather than inheriting whichever clock the previous
 *          case left behind.
 */
inline void reset_effect_globals() { hs_test::reset_globals(); }

/**
 * @brief Drives one effect type through construct -> init -> render -> read-back.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param name Effect name used in the [ok] / diagnostic output.
 * @details Verifies the effect constructs, init's, renders smoke_frames()
 * frames, and reads back every pixel without tripping an assert/OOB/hang, and
 * that get_pixel still aliases the displayed buffer once another frame has been
 * rendered and flipped in. Runs the dead-slider lint once on
 * the <SMALL_W,SMALL_H> pass, which both depth tiers execute. Finally requires
 * the effect to have overflowed no timeline event over its whole lifetime: a
 * drop is the one soft-degrade in the animation subsystem, and only this
 * per-effect delta attributes it.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void smoke_one(const char *name) {
  reset_effect_globals();
  const uint32_t dropped_before = Timeline::dropped_events();

  E<W, H> effect;
  effect.init();

  HS_EXPECT_EQ(effect.width(), W);
  HS_EXPECT_EQ(effect.height(), H);

  // A wider render margin costs rendered pixels, so it may only appear when a
  // filter's segment_margin asks for it; no roster pipeline asks past the
  // ClipRegion default today.
  HS_EXPECT_EQ(effect.clip().margin, ClipRegion{}.margin);

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

  std::printf("  [ok] %-20s rendered %d frames @ %dx%d (sum=%llu)\n", name,
              frames, W, H, static_cast<unsigned long long>(acc));

  if (!effect_may_be_dark(name, frames)) {
    if (acc == 0)
      std::printf("  ALL-BLACK %-20s produced no lit pixel over %d frames "
                  "@ %dx%d\n",
                  name, frames, W, H);
    HS_EXPECT(acc > 0, "effect must produce non-black output");
  }

  // The display read-back paths index display_buffer() directly (see
  // overrides_get_pixel), so get_pixel must still be its row-major view after
  // the flip advance_display performs.
  effect.draw_frame();
  effect.advance_display();
  const Pixel *displayed = effect.display_buffer();
  int unaliased = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (&effect.get_pixel(x, y) != displayed + y * W + x)
        ++unaliased;
  HS_EXPECT(unaliased == 0,
            "get_pixel must read the displayed buffer at the row-major offset");

  // The parameter lints are resolution-independent; run them once, on the pass
  // both depth tiers execute.
  if constexpr (W == SMALL_W && H == SMALL_H) {
    lint_dead_sliders(effect, name);
    lint_animated_pause(effect, name);
  }

  // A full timeline soft-drops: add() discards add_get()'s nullptr, so a chain
  // that re-arms itself from a .then() callback ends for good. The counter is
  // monotonic and process-wide, so this per-effect delta is the only thing that
  // attributes a drop to the effect that overflowed the table.
  const uint32_t dropped = Timeline::dropped_events() - dropped_before;
  if (dropped != 0)
    std::printf("  TIMELINE FULL %-20s dropped %u animation(s) over %d frames "
                "@ %dx%d (Timeline::MAX_EVENTS is %d)\n",
                name, static_cast<unsigned>(dropped), frames, W, H,
                Timeline::MAX_EVENTS);
  HS_EXPECT_EQ(dropped, 0u);
}

/**
 * @brief Build-time "registered-but-unread" lint for the live-art param system.
 * @param effect Effect instance whose editable params are probed.
 * @param name Effect name used in DEAD SLIDER diagnostic output.
 * @details Contract: a registered, editable param — one NOT flagged
 * animated and NOT flagged mark_readonly() — must be genuinely
 * user-controllable, i.e. a value written through updateParameter() must persist
 * across frames. If the engine overwrites it every frame (a Mutation/Driver/Lerp
 * bound to the same member, or output-only telemetry), the slider is dead: the
 * author must drive a private member and register_animated_param() it, or
 * mark_readonly()
 * pure telemetry. This is the build-time gate for the theme-4 dead-slider class
 * (it catches the per-frame overwrite mechanism behind MobiusGrid and the
 * preset-lerp group; a value that merely has no
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
    float target = (cur - def.min) > (def.max - cur) ? def.min + 0.25f * range
                                                     : def.min + 0.75f * range;
    // An integer target holds only whole numbers, so probe with a value it can
    // actually hold or every such param reads dead.
    if (def.is_integer()) {
      target = floorf(target);
      if (target == cur)
        continue;
    }
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
                  "(register_animated_param / mark_readonly / drive a private "
                  "member)\n",
                  name, def.name, static_cast<double>(target),
                  static_cast<double>(now));
    HS_EXPECT(persisted, "editable param must persist across frames");
  }
}

/**
 * @brief Verifies every advertised animated param remains fixed while paused.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in PAUSE LEAK diagnostic output.
 * @details Parameters are audited one at a time so independent valid controls
 *          are not combined into an invalid cross-field configuration. The
 *          aggregate audit spans at least 500 rendered frames.
 */
inline void lint_animated_pause(Effect &effect, const char *name) {
  std::vector<const char *> names;
  std::vector<float> original;
  std::vector<float> target;
  names.reserve(effect.getParameters().size());
  original.reserve(effect.getParameters().size());
  target.reserve(effect.getParameters().size());
  for (const auto &def : effect.getParameters()) {
    if (!def.animated)
      continue;
    // The write below lands on the requested value, so a schema-driven effect
    // whose rendered slot is canonicalized away from it only round-trips
    // through get_requested().
    const float current = def.get_requested();
    names.push_back(def.name);
    original.push_back(current);
    target.push_back(def.is_bool()
                         ? (current > 0.5f ? 0.0f : 1.0f)
                         : ((current - def.min) > (def.max - current)
                                ? def.min + 0.25f * (def.max - def.min)
                                : def.min + 0.75f * (def.max - def.min)));
  }
  const size_t count = names.size();
  if (count == 0)
    return;

  constexpr int PAUSE_AUDIT_FRAMES = 500;
  const int frames_per_param =
      std::max(4, static_cast<int>((PAUSE_AUDIT_FRAMES + count - 1) / count));
  bool leaked = false;
  for (size_t index = 0; index < count; ++index) {
    for (size_t restore = 0; restore < count; ++restore)
      effect.updateParameter(names[restore], original[restore]);
    HS_EXPECT_TRUE(effect.updateParameter(names[index], target[index]) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(effect.animations_paused());
    effect.draw_frame();
    effect.advance_display();
    const float held = effect.getParameters().find(names[index])->get();
    for (int frame = 0; frame < frames_per_param; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      const auto *def = effect.getParameters().find(names[index]);
      if (def->get() != held) {
        if (!leaked)
          std::printf("  PAUSE LEAK %s::%s moved %.4f -> %.4f at frame %d\n",
                      name, names[index], static_cast<double>(held),
                      static_cast<double>(def->get()), frame + 1);
        leaked = true;
      }
    }
  }
  HS_EXPECT(!leaked, "animated params must remain fixed on every paused frame");
}

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
 * Timeline, pole-LOD knob, scan counters) and pins the mock clock to the frame
 * cadence, so two calls start from an identical state.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void render_capture(std::vector<Pixel> &out, int frames,
                           uint64_t *frame_fold = nullptr) {
  reset_effect_globals();
  // Pre-init epoch, so construction sees the same clock on both runs.
  pin_frame_clock(0);

  uint64_t fold = hs_test::FNV1A64_BASIS;
  const auto fold_byte = [&fold](uint8_t byte) {
    fold = hs_test::fnv1a64_byte(fold, byte);
  };

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    pin_frame_clock(f);
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

  capture_frame<W, H>(effect, out);
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
  hs::random().seed(0xC0FFEEu); // off the canonical seed(1337)
  global_timeline_t = 0x5EED;   // off zero
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
/**
 * @brief Frames per segment for an effect that is still allowed to be dark at
 *        PARITY_FRAMES.
 * @details Comparing two all-black renders always agrees, so a slow-starting
 *          effect gets a sweep long enough to leave effect_may_be_dark()'s
 *          window instead of a pass that means nothing.
 */
constexpr int PARITY_FRAMES_SLOW = 64;
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
template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void clip_clear_parity_one(const char *name) {
  constexpr int S = H * 2;
  const int frames = effect_may_be_dark(name, PARITY_FRAMES)
                         ? PARITY_FRAMES_SLOW
                         : PARITY_FRAMES;

  auto render = [&](int segment_id, bool full_clear) {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    std::vector<Pixel> displayed;
    E<W, H> effect;
    effect.init();
    effect.force_full_buffer_clear = full_clear;
    const pov::SegmentMap map =
        pov::segment_map(segment_id, S, PARITY_SEGMENTS);
    for (int f = 0; f < frames; ++f) {
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

  size_t lit = 0, displayed_pixels = 0;
  for (int segment_id = 0; segment_id < PARITY_SEGMENTS; ++segment_id) {
    const std::vector<Pixel> full = render(segment_id, true);
    const std::vector<Pixel> clipped = render(segment_id, false);
    HS_EXPECT_EQ(full.size(), clipped.size());

    size_t different = 0;
    for (size_t i = 0; i < full.size() && i < clipped.size(); ++i) {
      if (full[i] != clipped[i])
        ++different;
      if (full[i].r | full[i].g | full[i].b)
        ++lit;
    }
    displayed_pixels += full.size();
    if (different)
      std::printf("  CLIP-CLEAR DRIFT %-20s segment %d: %zu of %zu displayed "
                  "pixels differ from the full-buffer clear\n",
                  name, segment_id, different, full.size());
    HS_EXPECT(different == 0,
              "clip clearing must not change any displayed pixel");
  }
  // Two all-black renders agree pixel for pixel, so the comparison above only
  // means something once the sweep has produced output. Sparse effects can miss
  // a single arm segment, so the requirement spans the whole sweep, and the
  // sweep length is chosen so no effect is still exempt at it — an exemption
  // here would make every assertion above vacuous.
  HS_EXPECT(!effect_may_be_dark(name, frames),
            "clip-clear parity must run past the all-black exemption window");
  if (lit == 0)
    std::printf("  CLIP-CLEAR DARK %-20s no lit pixel over %zu displayed "
                "pixels\n",
                name, displayed_pixels);
  HS_EXPECT(lit > 0, "clip-clear parity must compare a lit render");
  hs::clear_mock_time();
}

/** @brief Frames rendered per effect by the paused-render sweep. */
constexpr int PAUSED_FRAMES = 4;

/**
 * @brief Renders one effect with animations paused from before init() and
 *        requires the frame to light up.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param name Effect name used in the diagnostic output.
 * @details Mirrors the WASM load path (targets/wasm/engine_bindings.h), which
 * applies the retained pause before init() and re-engages it on every APPLIED
 * write to an animated param: every pausable event an effect schedules in init()
 * is then held at its first frame for as long as the user leaves the pause on.
 * A pause is not a reason to render nothing — ambient motion keeps running and
 * the sphere must stay lit.
 */
template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void paused_render_one(const char *name) {
  reset_effect_globals();

  E<W, H> effect;
  effect.setAnimationsPaused(true);
  effect.init();
  HS_EXPECT_TRUE(effect.animations_paused());
  for (int f = 0; f < PAUSED_FRAMES; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  uint64_t acc = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = effect.get_pixel(x, y);
      acc += static_cast<uint64_t>(p.r) + p.g + p.b;
    }

  if (!effect_may_be_dark(name, PAUSED_FRAMES)) {
    if (acc == 0)
      std::printf("  PAUSED-BLANK %-20s produced no lit pixel over %d paused "
                  "frames @ %dx%d\n",
                  name, PAUSED_FRAMES, W, H);
    HS_EXPECT(acc > 0, "effect must produce non-black output while paused");
  }
}

/**
 * @brief Roster sweep: every registered effect renders lit while paused.
 * @details The pause reaches an effect through whatever it schedules on the
 * timeline, so the gate belongs on the whole roster rather than on the one
 * effect a pause blanked.
 */
inline void test_every_effect_renders_while_paused() {
  std::printf("  -- paused render, %d frames --\n", PAUSED_FRAMES);
#define HS_PAUSED_ONE(name) paused_render_one<name>(#name);
  HS_EFFECT_LIST(HS_PAUSED_ONE)
#undef HS_PAUSED_ONE
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

// Highest degree the SH cases pin: one above the effect's MAX_DEGREE, so a
// widening lands on already-covered ground.
constexpr int SH_PIN_MAX_DEGREE = 5;

/** @brief One closed-form reduced Legendre polynomial. */
struct SHReducedLegendreRow {
  int l;
  int m;
  double c[SH_PIN_MAX_DEGREE + 1]; /**< Ascending-power coefficients in x. */
};

// P_l^m(x) / (1 - x²)^(m/2) as an explicit polynomial. Condon-Shortley
// P_l^m = (-1)^m (1 - x²)^(m/2) d^m/dx^m P_l, so each row is (-1)^m times the
// m-th derivative of the Legendre polynomial P_l — a closed form owing nothing
// to the recurrence under test. Covers every (l, m) with 0 <= m <= l <= 5.
inline constexpr SHReducedLegendreRow SH_REDUCED_LEGENDRE[] = {
    {0, 0, {1.0}},
    {1, 0, {0.0, 1.0}},
    {1, 1, {-1.0}},
    {2, 0, {-0.5, 0.0, 1.5}},
    {2, 1, {0.0, -3.0}},
    {2, 2, {3.0}},
    {3, 0, {0.0, -1.5, 0.0, 2.5}},
    {3, 1, {1.5, 0.0, -7.5}},
    {3, 2, {0.0, 15.0}},
    {3, 3, {-15.0}},
    {4, 0, {0.375, 0.0, -3.75, 0.0, 4.375}},
    {4, 1, {0.0, 7.5, 0.0, -17.5}},
    {4, 2, {-7.5, 0.0, 52.5}},
    {4, 3, {0.0, -105.0}},
    {4, 4, {105.0}},
    {5, 0, {0.0, 1.875, 0.0, -8.75, 0.0, 7.875}},
    {5, 1, {-1.875, 0.0, 26.25, 0.0, -39.375}},
    {5, 2, {0.0, -52.5, 0.0, 157.5}},
    {5, 3, {52.5, 0.0, -472.5}},
    {5, 4, {0.0, 945.0}},
    {5, 5, {-945.0}},
};
static_assert(sizeof(SH_REDUCED_LEGENDRE) / sizeof(SH_REDUCED_LEGENDRE[0]) ==
                  (SH_PIN_MAX_DEGREE + 1) * (SH_PIN_MAX_DEGREE + 2) / 2,
              "closed-form table must hold every (l, m) up to "
              "SH_PIN_MAX_DEGREE");

/** @brief Evaluates a closed-form table row at x, in double. */
inline double sh_closed_form_reduced(const SHReducedLegendreRow &row,
                                     double x) {
  double sum = 0.0;
  for (int k = SH_PIN_MAX_DEGREE; k >= 0; --k)
    sum = sum * x + row.c[k];
  return sum;
}

/** @brief Textbook associated Legendre P_l^m(x) in double, as test truth. */
inline double sh_reference_legendre(int l, int m, double x) {
  double somx2 = std::sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
  double value = 0.0;
  for (const SHReducedLegendreRow &row : SH_REDUCED_LEGENDRE)
    if (row.l == l && row.m == m)
      value = sh_closed_form_reduced(row, x);
  for (int i = 0; i < m; ++i)
    value *= somx2;
  return value;
}

/**
 * @brief Pins SHMath::reduced_legendre against the closed-form P_l^m table.
 * @details The hand-derived upward recurrence starts at a unit seed, steps to
 * P_{m+1}^m, then recurs in l, with the x-independent P_m^m factored out into
 * legendre_seed(); a wrong coefficient in any of the three shifts the result by
 * an O(1) relative amount. Sweeps every table row (so all three branches run for
 * degrees the effect can reach and one above) across x = cos(phi) in [-1, 1].
 */
inline void test_sh_reduced_legendre_matches_closed_form() {
  // Relative tolerance: the float recurrence tops out at 3.8e-6 measured under
  // -O2 -ffast-math, so 3e-5 leaves ~8x reassociation headroom while still
  // catching any coefficient error, which is O(1) relative.
  constexpr double REDUCED_LEGENDRE_REL_TOL = 3e-5;
  constexpr int SAMPLES = 200;

  double worst = 0.0;
  int worst_l = 0, worst_m = 0;
  for (const SHReducedLegendreRow &row : SH_REDUCED_LEGENDRE) {
    for (int i = 0; i <= SAMPLES; ++i) {
      const double x = -1.0 + 2.0 * i / SAMPLES;
      const double want = sh_closed_form_reduced(row, x);
      const double got =
          SHMath::legendre_seed(row.m) *
          SHMath::reduced_legendre(row.l, row.m, static_cast<float>(x));
      // Relative above unit magnitude; rows reach 945, so a plain absolute
      // bound would be either meaningless there or unreachable near zero.
      const double err = std::fabs(got - want) / std::max(1.0, std::fabs(want));
      if (err > worst) {
        worst = err;
        worst_l = row.l;
        worst_m = row.m;
      }
    }
  }

  if (worst >= REDUCED_LEGENDRE_REL_TOL)
    std::printf("  SH reduced Legendre: worst relative error %g at l=%d m=%d\n",
                worst, worst_l, worst_m);
  HS_EXPECT(worst < REDUCED_LEGENDRE_REL_TOL,
            "reduced Legendre recurrence must match the closed-form P_l^m");
}

/**
 * @brief Verifies the Cartesian spherical harmonic matches the spherical form.
 * @details SHMath::spherical_harmonic drops the sin(phi)^|m| factor from the
 * Legendre term and recovers it from Re/Im((x + iz)^|m|). That identity holds
 * only for unit directions, so this sweeps a full (phi, theta) lattice against
 * a double-precision P_l^m(cos phi) * cos/sin(|m| theta), covering every (l, m)
 * up to SH_PIN_MAX_DEGREE.
 */
inline void test_sh_cartesian_matches_spherical() {
  double worst = 0.0;
  int worst_l = 0, worst_m = 0;

  for (int l = 0; l <= SH_PIN_MAX_DEGREE; ++l) {
    for (int m = -l; m <= l; ++m) {
      const float N = SHMath::normalization(l, m);
      const float scale = SHMath::harmonic_scale(l, m);
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
              std::fabs(SHMath::spherical_harmonic(l, m, p, scale) - want);
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

/**
 * @brief White-box accessor for SphericalHarmonics' morph chain and shader inputs.
 * @details Befriended in effects/SphericalHarmonics.h. The morph legs are 64
 *          frames and the chain re-arms only from a Transition::then() callback,
 *          so a frozen chain renders a perfectly valid (static) sphere that the
 *          smoke and determinism sweeps both accept. This seam reads the mode
 *          bookkeeping and the palette/orientation the draw_frame shader consumed.
 */
struct SphericalHarmonicsWhiteBox {
  using SH = SphericalHarmonics<SMALL_W, SMALL_H>;
  using Field = SH::HarmonicField;

  static int current_idx(const SH &fx) { return fx.current_idx; }
  static int next_idx(const SH &fx) { return fx.next_idx; }
  static float morph_alpha(const SH &fx) { return fx.morph_alpha; }
  static Quaternion orientation(const SH &fx) { return fx.orientation.get(); }
  static const BakedPalette &palette(const SH &fx) { return fx.baked_palette; }
  static float amplitude(const SH &fx) { return fx.params.amplitude; }

  // Pinning both morph endpoints on one mode makes the blend an identity, so a
  // frame renders one pure harmonic whatever morph_alpha has reached. The
  // production roll excludes this state; only a test wants it.
  static void pin_mode(SH &fx, int idx) {
    fx.current_idx = idx;
    fx.next_idx = idx;
  }
};

/**
 * @brief Pins HarmonicField's write-through, blend endpoints, and local frame.
 * @details The field is the only channel between the harmonic math and the
 *          shader: distance() must report a constant "inside" so the whole
 *          sphere rasterizes, and the value must ride out through raw_dist.
 *          Nothing downstream can tell a frozen blend, a reversed blend, or a
 *          shape rotated the wrong way from a legitimate mode — every one of
 *          them still paints a plausible sphere.
 */
inline void test_sh_field_write_through_and_endpoints() {
  using Field = SphericalHarmonicsWhiteBox::Field;
  const Quaternion spin =
      make_rotation(Vector(0.3f, 0.8f, -0.5f).normalized(), 1.1f);
  const Quaternion identity;

  constexpr int LA = 1, MA = 0, LB = 3, MB = 2;
  Field mix_start(LA, MA, LB, MB, 0.0f, spin);
  Field mix_end(LA, MA, LB, MB, 1.0f, spin);
  Field pure_a(LA, MA, LA, MA, 0.0f, spin);
  Field pure_b(LB, MB, LB, MB, 0.0f, spin);
  Field pure_a_unrotated(LA, MA, LA, MA, 0.0f, identity);

  constexpr int PHI_STEPS = 24, THETA_STEPS = 32;
  int samples = 0, inside = 0, positives = 0, negatives = 0;
  double worst_start = 0.0, worst_end = 0.0, worst_frame = 0.0;
  for (int i = 0; i <= PHI_STEPS; ++i) {
    const float phi = PI_F * i / PHI_STEPS;
    for (int j = 0; j < THETA_STEPS; ++j) {
      const float theta = 2.0f * PI_F * j / THETA_STEPS;
      const Vector p(sinf(phi) * cosf(theta), cosf(phi),
                     sinf(phi) * sinf(theta));
      SDF::DistanceResult start, end, a, b, spun, unspun;
      mix_start.distance(p, start);
      mix_end.distance(p, end);
      pure_a.distance(p, a);
      pure_b.distance(p, b);
      pure_a_unrotated.distance(p, unspun);
      // Rotating the sample by the same quaternion the shape carries must land
      // back on the unrotated shape's value.
      pure_a.distance(rotate(p, spin), spun);

      ++samples;
      if (start.dist < 0.0f)
        ++inside;
      if (a.raw_dist > 0.02f)
        ++positives;
      if (a.raw_dist < -0.02f)
        ++negatives;
      worst_start =
          std::max<double>(worst_start, std::fabs(start.raw_dist - a.raw_dist));
      worst_end =
          std::max<double>(worst_end, std::fabs(end.raw_dist - b.raw_dist));
      worst_frame = std::max<double>(
          worst_frame, std::fabs(spun.raw_dist - unspun.raw_dist));
    }
  }

  // The whole sphere is covered: a positive distance anywhere would punch a hole.
  HS_EXPECT_EQ(inside, samples);
  // A dipole must actually change sign, or the polarity split below is vacuous.
  HS_EXPECT_GT(positives, 0);
  HS_EXPECT_GT(negatives, 0);

  if (worst_start > 0.0 || worst_end > 1e-5 || worst_frame > 1e-5)
    std::printf("  SH field: blend0 err=%g blend1 err=%g frame err=%g\n",
                worst_start, worst_end, worst_frame);
  // blend == 0 is the first mode and blend == 1 is the second, exactly.
  HS_EXPECT_EQ(worst_start, 0.0);
  HS_EXPECT_LT(worst_end, 1e-5);
  HS_EXPECT_LT(worst_frame, 1e-5);
}

/**
 * @brief Renders one frame of a mode-pinned SphericalHarmonics and inspects it.
 * @param idx Flat harmonic index pinned on both morph endpoints.
 * @param amplitude Value written to the Amplitude slider before the frame.
 * @param inspect Callback receiving (effect, field the shader saw, amplitude).
 */
template <typename FnT>
inline void sh_render_pinned_mode(int idx, float amplitude, FnT &&inspect) {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SH fx;
  fx.init();
  WB::pin_mode(fx, idx);
  HS_EXPECT_TRUE(fx.updateParameter("Amplitude", amplitude) ==
                 ParamSetResult::APPLIED);
  fx.draw_frame();
  fx.advance_display();

  auto [l, m] = SHMath::decode_lm(idx);
  WB::Field field(l, m, l, m, WB::morph_alpha(fx), WB::orientation(fx));
  inspect(fx, field, WB::amplitude(fx));
  hs::clear_mock_time();
}

/**
 * @brief Pins the diverging palette split and the ambient-occlusion shaping.
 * @details Reads the rendered frame back and classifies each pixel by the sign
 *          and magnitude the field carries there, so the assertions are on
 *          observed colors rather than on a second copy of the shader. Covers
 *          three contracts nothing else touches: a saturated positive lobe wears
 *          the palette's top color unmodified (which is also AO_AMBIENT +
 *          AO_RANGE == 1, i.e. no residual dimming at saturation), the negative
 *          lobe wears that same color with red and blue swapped and green
 *          dimmed, and below the AO falloff every pixel is scaled well under the
 *          palette entry its own magnitude selects.
 */
inline void test_sh_polarity_split_and_ao_shaping() {
  using WB = SphericalHarmonicsWhiteBox;
  // idx 2 is (l=1, m=0): one nodal circle, so both polarities cover a wide band.
  constexpr int DIPOLE_IDX = 2;

  // Amplitude 10 (the slider maximum) saturates the palette across most of both
  // lobes, leaving a wide margin above the |val| the classification thresholds
  // on even if the rasterizer's sample point differs from ours in the last bits.
  sh_render_pinned_mode(
      DIPOLE_IDX, 10.0f,
      [](const WB::SH &fx, const WB::Field &field, float amp) {
        constexpr float SATURATED =
            1.5f; // |val| * amp; palette saturates at 1.0
        Pixel pos(0, 0, 0), neg(0, 0, 0);
        int pos_n = 0, neg_n = 0, pos_split = 0, neg_split = 0;
        uint64_t pos_green = 0, neg_green = 0;
        for (int y = 0; y < SMALL_H; ++y)
          for (int x = 0; x < SMALL_W; ++x) {
            SDF::DistanceResult res;
            field.distance(pixel_to_vector<SMALL_W, SMALL_H>(x, y), res);
            const Pixel &px = fx.get_pixel(x, y);
            const bool positive = res.raw_dist > 0.0f;
            // A dipole's two lobes are antipodal mirrors, so their magnitude
            // distributions match and their green totals differ only by the
            // negative recolor's scale.
            (positive ? pos_green : neg_green) += px.g;

            const float mag = std::fabs(res.raw_dist) * amp;
            if (mag < SATURATED)
              continue;
            Pixel &slot = positive ? pos : neg;
            int &count = positive ? pos_n : neg_n;
            int &split = positive ? pos_split : neg_split;
            if (count++ == 0)
              slot = px;
            else if (px.r != slot.r || px.g != slot.g || px.b != slot.b)
              ++split;
          }

        std::printf("  SH polarity: pos=%d (%u,%u,%u) neg=%d (%u,%u,%u) "
                    "green %llu vs %llu\n",
                    pos_n, pos.r, pos.g, pos.b, neg_n, neg.r, neg.g, neg.b,
                    static_cast<unsigned long long>(pos_green),
                    static_cast<unsigned long long>(neg_green));
        HS_EXPECT_GT(pos_n, 0);
        HS_EXPECT_GT(neg_n, 0);
        // Past saturation the palette index, the seam weight and the AO factor are
        // all pinned, so every pixel of a lobe must carry one single color.
        HS_EXPECT_EQ(pos_split, 0);
        HS_EXPECT_EQ(neg_split, 0);

        // Positive lobe: the palette color, undimmed.
        const Pixel top = WB::palette(fx).get(1.0f).color;
        HS_EXPECT_EQ(pos.r, top.r);
        HS_EXPECT_EQ(pos.g, top.g);
        HS_EXPECT_EQ(pos.b, top.b);
        // Negative lobe: red and blue traded.
        HS_EXPECT_EQ(neg.r, pos.b);
        HS_EXPECT_EQ(neg.b, pos.r);
        // ...and the trade is visible, i.e. the two lobes are not the same color.
        HS_EXPECT_TRUE(pos.r != pos.b);
        // The green dimming is invisible at saturation (this palette's top color
        // has none), so pin it over the whole mirrored lobe instead.
        HS_EXPECT_GT(pos_green, 0u);
        HS_EXPECT_TRUE(neg_green * 10 < pos_green * 9);
      });

  // Amplitude 0.6 keeps every pixel below the AO falloff, so the whole frame
  // must read as a dimmed copy of the palette rather than the palette itself.
  sh_render_pinned_mode(
      DIPOLE_IDX, 0.6f,
      [](const WB::SH &fx, const WB::Field &field, float amp) {
        // Only sample where the palette entry is bright enough that a missing
        // occlusion factor would be unambiguous.
        constexpr uint32_t PALETTE_FLOOR = 4096;
        int probed = 0, undimmed = 0;
        for (int y = 0; y < SMALL_H; ++y)
          for (int x = 0; x < SMALL_W; ++x) {
            SDF::DistanceResult res;
            field.distance(pixel_to_vector<SMALL_W, SMALL_H>(x, y), res);
            const float mag = std::fabs(res.raw_dist) * amp;
            const Pixel want = WB::palette(fx).get(std::min(1.0f, mag)).color;
            const uint32_t want_energy =
                static_cast<uint32_t>(want.r) + want.g + want.b;
            if (want_energy < PALETTE_FLOOR)
              continue;
            const Pixel &px = fx.get_pixel(x, y);
            const uint32_t got_energy =
                static_cast<uint32_t>(px.r) + px.g + px.b;
            ++probed;
            // The occlusion factor cannot reach 0.8 at this amplitude.
            if (got_energy * 5 >= want_energy * 4)
              ++undimmed;
          }
        std::printf("  SH ambient occlusion: probed=%d undimmed=%d\n", probed,
                    undimmed);
        HS_EXPECT_GT(probed, 0);
        HS_EXPECT_EQ(undimmed, 0);
      });
}

/**
 * @brief Verifies the morph chain re-arms itself and keeps advancing modes.
 * @details start_morph() schedules a 64-frame Transition whose then() callback
 *          commits the target and calls start_morph() again; a chain that failed
 *          to re-arm would freeze the sphere on one mode while every smoke,
 *          determinism and math check still passed. Drives enough frames for at
 *          least three commits, so the callback must have re-armed twice.
 */
inline void test_sh_morph_chain_rearms() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SH fx;
  fx.init();

  const int seed = WB::current_idx(fx);
  HS_EXPECT_GT(seed, 0); // never the constant harmonic
  HS_EXPECT_TRUE(WB::next_idx(fx) != seed);

  constexpr int FRAMES = 260; // four 64-frame legs
  int commits = 0, held = seed, alpha_out_of_range = 0, self_blend = 0;
  int rearmed_at_zero = 0;
  float alpha_peak = 0.0f;
  std::vector<int> visited{seed};
  for (int f = 0; f < FRAMES; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();

    const float alpha = WB::morph_alpha(fx);
    alpha_peak = std::max(alpha_peak, alpha);
    if (!(alpha >= 0.0f && alpha <= 1.0f))
      ++alpha_out_of_range;
    if (WB::next_idx(fx) == WB::current_idx(fx))
      ++self_blend;

    const int now = WB::current_idx(fx);
    if (now != held) {
      ++commits;
      held = now;
      // A committed leg rewinds the blend and schedules the next one.
      if (alpha == 0.0f)
        ++rearmed_at_zero;
      if (std::find(visited.begin(), visited.end(), now) == visited.end())
        visited.push_back(now);
    }
  }
  hs::clear_mock_time();

  std::printf("  SH morph: %d commits over %d frames, %zu distinct modes, "
              "alpha peak %.3f\n",
              commits, FRAMES, visited.size(), static_cast<double>(alpha_peak));
  HS_EXPECT_GE(commits, 3); // the then() callback re-armed at least twice
  HS_EXPECT_EQ(rearmed_at_zero, commits);
  HS_EXPECT_GE(visited.size(), 3u);
  HS_EXPECT_EQ(alpha_out_of_range, 0);
  HS_EXPECT_EQ(self_blend, 0); // blending a mode into itself would freeze it
  HS_EXPECT_GT(alpha_peak, 0.9f);
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
 *          fixed 7680-node graph, independent of <W,H>, so the small-aspect
 *          resolution is used arbitrarily.
 */
struct GSWhiteBox {
  using GS = GSReactionDiffusion<SMALL_W, SMALL_H>;
  static constexpr int N = GS::RD_N;

  static uint16_t to_q16(float v) { return GS::to_q16(v); }
  static float from_q16(uint16_t v) { return GS::from_q16(v); }
  static Color4 palette_sample(const GS &gs, float t) {
    return gs.palette.get(t);
  }
  static void start_reaction(GS &gs) { gs.start_reaction(); }
  static void set_params(GS &gs, float feed, float k, float dA, float dB,
                         float dt) {
    gs.params.feed = feed;
    gs.params.k = k;
    gs.params.d_a = dA;
    gs.params.d_b = dB;
    gs.params.dt = dt;
  }
  static int dissolve_frame(const GS &gs) {
    return gs.transition.dissolve_frames;
  }
  static const uint16_t *b_field(const GS &gs) { return gs.state.B; }
  static const uint16_t *a_field(const GS &gs) { return gs.state.A; }
  static float dissolve_hash(int i, uint32_t seed) {
    return ::hash01(static_cast<uint32_t>(i), seed);
  }
  static void convert(GS &gs, float phase, uint32_t seed) {
    gs.transition.dissolve_seed = seed;
    gs.convert_below(phase);
  }
  static void set_node(GS &gs, int i, uint16_t a, uint16_t b) {
    gs.state.A[i] = a;
    gs.state.B[i] = b;
  }
  static constexpr float DISSOLVE_FADE_FRACTION = GS::DISSOLVE_FADE_FRACTION;
  static constexpr int DISSOLVE_FRAMES = GS::DISSOLVE_FRAMES;
  static constexpr int MIN_GROW_FRAMES = GS::MIN_GROW_FRAMES;
  static constexpr int STABLE_HOLD_FRAMES = GS::STABLE_HOLD_FRAMES;

  // refine_render_center's early-out certificates.
  static constexpr int POLE_BAND = GS::POLE_BAND;
  static constexpr float BULK_MIN_SPACING_FRAC = GS::BULK_MIN_SPACING_FRAC;
  static constexpr float BULK_CERTIFIED_D2 = GS::BULK_CERTIFIED_D2;
  static constexpr float POLE_CERTIFIED_D2 = GS::POLE_CERTIFIED_D2;

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

  static void step_float_reference(const GS &gs, const float *cA,
                                   const float *cB, float *nA, float *nB) {
    const float dt = gs.params.dt * GS::STEP_DT_SCALE;
    for (int i = 0; i < N; ++i) {
      float a = cA[i];
      float b = cB[i];
      float sumA = 0.0f, sumB = 0.0f;
      for (int k = 0; k < ReactionGraph::RD_K; ++k) {
        int ni = ReactionGraph::neighbors[i][k];
        sumA += cA[ni];
        sumB += cB[ni];
      }
      float lA = sumA - ReactionGraph::RD_K * a;
      float lB = sumB - ReactionGraph::RD_K * b;
      float abb = a * b * b;
      nA[i] = hs::clamp(
          a + (gs.params.d_a * lA - abb + gs.params.feed * (1.0f - a)) * dt,
          0.0f, 1.0f);
      nB[i] = hs::clamp(
          b + (gs.params.d_b * lB + abb - (gs.params.k + gs.params.feed) * b) *
                  dt,
          0.0f, 1.0f);
    }
  }

  struct ShaderError {
    int different = 0;
    int lit = 0;
    int coverage = 0;
    int hard = 0;
    int stencil_samples = 0;
    int stencil_changes = 0;
    int center_mismatches = 0;
    int max_channel = 0;
    uint64_t total_channel = 0;
  };

  template <int W, int H> static ShaderError shared_shader_error(GS &gs) {
    ScratchScope guard(scratch_arena_a);
    Vector *world_nodes = gs.orient_lattice();
    using Grid = Scan::Shader::SsaaGrid<W, H>;
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    Grid grid;
    ShaderError error;
    for (int y = 0; y < H; ++y) {
      grid.set_row(y);
      for (int x = 0; x < W; ++x) {
        Fragment frag;
        frag.pos = pixel_to_vector<W, H>(x, y);
        gs.seed_face_lut(frag);
        int seed = static_cast<int>(frag.v0);
        int shared_center =
            gs.refine_render_center(frag.pos, world_nodes, seed);
        if (shared_center != gs.refine_center(frag.pos, world_nodes, seed))
          ++error.center_mismatches;
        Pixel got = gs.shade_pixel(seed, frag.pos, world_nodes, grid, x);

        Pixel expected(0, 0, 0);
        for (int i = 0; i < 4; ++i) {
          Vector sample_v = grid.at(x, i);
          ++error.stencil_samples;
          if (gs.refine_center(sample_v, world_nodes, seed) != shared_center)
            ++error.stencil_changes;
          float b = gs.interpolate_b(sample_v, seed, world_nodes);
          if (b < GS::B_CULL_THRESHOLD)
            continue;
          float t = hs::clamp((b - GS::B_COLOR_FLOOR) * GS::B_COLOR_SCALE, 0.0f,
                              1.0f);
          Color4 sample = gs.palette.get(t);
          expected += sample.color * (sample.alpha * 0.25f);
        }
        int dr = std::abs(static_cast<int>(got.r) - expected.r);
        int dg = std::abs(static_cast<int>(got.g) - expected.g);
        int db = std::abs(static_cast<int>(got.b) - expected.b);
        int pixel_max = std::max(dr, std::max(dg, db));
        bool got_lit = got.r || got.g || got.b;
        bool expected_lit = expected.r || expected.g || expected.b;
        if (got != expected)
          ++error.different;
        if (got_lit || expected_lit)
          ++error.lit;
        if (got_lit != expected_lit)
          ++error.coverage;
        if (pixel_max > 4096)
          ++error.hard;
        error.max_channel = std::max(error.max_channel, pixel_max);
        error.total_channel += dr + dg + db;
      }
    }
    Pixel culled = gs.shade_pixel(-1, Vector(), world_nodes, grid, 0);
    if (culled != Pixel(0, 0, 0))
      ++error.different;
    return error;
  }
};

/** @brief Pins GS's opaque quarter-sample accumulation to Pixel arithmetic. */
inline void test_gs_opaque_quarter_accumulation_is_exact() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
    HS_EXPECT_EQ(GSWhiteBox::palette_sample(
                     gs, static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1))
                     .alpha,
                 1.0f);

  uint32_t sum = 0;
  Pixel accumulated(0, 0, 0);
  for (uint32_t channel = 0; channel <= 65535u; ++channel) {
    Pixel sample(static_cast<uint16_t>(channel), 0, 0);
    accumulated = Pixel(0, 0, 0);
    sum = 0;
    for (int i = 0; i < 4; ++i) {
      accumulated += sample * 0.25f;
      sum += (channel + 2u) >> 2;
    }
    HS_EXPECT_EQ(accumulated.r,
                 static_cast<uint16_t>(sum > 65535u ? 65535u : sum));
  }
}

/** @brief Verifies each reseeded reaction receives a freshly generated palette. */
inline void test_gs_reseed_generates_palette() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  std::array<Pixel, BakedPalette::LUT_SIZE> before;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
    before[i] = GSWhiteBox::palette_sample(gs, static_cast<float>(i) /
                                                   (BakedPalette::LUT_SIZE - 1))
                    .color;

  const size_t arena_offset = persistent_arena.get_offset();
  GSWhiteBox::start_reaction(gs);

  int changed = 0;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const Pixel after =
        GSWhiteBox::palette_sample(gs, static_cast<float>(i) /
                                           (BakedPalette::LUT_SIZE - 1))
            .color;
    changed += after != before[i];
  }
  HS_EXPECT_GT(changed, 0);
  HS_EXPECT_EQ(persistent_arena.get_offset(), arena_offset);
}

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
 * @brief Re-measures the lattice and requires refine_render_center's early-out
 *        certificates to bound the true minimum node spacing.
 * @details The early-out returns the seed unwalked whenever the query lies
 *          within sqrt(safe_d2) of it. That is the nearest node only while
 *          safe_d2 stays at or under a quarter of the seed's true
 *          nearest-neighbor distance squared, so BULK_CERTIFIED_D2 /
 *          POLE_CERTIFIED_D2 / BULK_MIN_SPACING_FRAC / POLE_BAND are all
 *          hand-measured properties of the generated lattice. The shipped
 *          static_asserts only cross-check them against each other: a change to
 *          RD_N, to the generator, or to the index-to-latitude mapping would
 *          leave every one of them silently wrong and the early-out returning a
 *          non-nearest node. This recomputes them from node() and fails if the
 *          shipped constants stop bounding the measurement.
 *
 *          The per-node minimum is exact, not sampled: node()'s y is strictly
 *          decreasing in the index and chord distance is at least |dy|, so
 *          scanning outward from each index and stopping once |dy| reaches the
 *          running best cannot skip a closer node.
 */
inline void test_gs_render_certificates_bound_lattice() {
  const int n = GSWhiteBox::N;
  std::vector<Vector> lattice(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    lattice[static_cast<size_t>(i)] = ReactionGraph::node(i);

  const float bulk_floor =
      GSWhiteBox::BULK_MIN_SPACING_FRAC * ReactionGraph::D_AVG;
  float bulk_min_d2 = 4.0f; // antipodal chord squared: the loosest possible
  float pole_min_d2 = 4.0f;
  int tight_band = 0; // band width the sub-bulk-spacing nodes actually need
  for (int i = 0; i < n; ++i) {
    const Vector &p = lattice[static_cast<size_t>(i)];
    float best = 4.0f;
    for (int step = -1; step <= 1; step += 2) {
      for (int j = i + step; j >= 0 && j < n; j += step) {
        const Vector &q = lattice[static_cast<size_t>(j)];
        float dy = p.y - q.y;
        if (dy * dy >= best)
          break;
        float dx = p.x - q.x, dz = p.z - q.z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best)
          best = d2;
      }
    }
    if (i < GSWhiteBox::POLE_BAND || i >= n - GSWhiteBox::POLE_BAND) {
      pole_min_d2 = std::min(pole_min_d2, best);
    } else {
      bulk_min_d2 = std::min(bulk_min_d2, best);
    }
    if (best < bulk_floor * bulk_floor)
      tight_band = std::max(tight_band, std::min(i, n - 1 - i) + 1);
  }

  const float bulk_frac = std::sqrt(bulk_min_d2) / ReactionGraph::D_AVG;
  std::printf("  [info] GS lattice: bulk min spacing %.6f*D_AVG, bulk d2/4 "
              "%.9f, pole d2/4 %.9f, tight band %d\n",
              bulk_frac, 0.25f * bulk_min_d2, 0.25f * pole_min_d2, tight_band);

  HS_EXPECT_LE(GSWhiteBox::BULK_MIN_SPACING_FRAC, bulk_frac);
  HS_EXPECT_LE(GSWhiteBox::BULK_CERTIFIED_D2, 0.25f * bulk_min_d2);
  HS_EXPECT_LE(GSWhiteBox::POLE_CERTIFIED_D2, 0.25f * pole_min_d2);
  // Every node packed tighter than the bulk minimum must fall inside the band
  // that gets the smaller certificate.
  HS_EXPECT_LE(tight_band, GSWhiteBox::POLE_BAND);
}

/**
 * @brief Bounds the shared stencil against independent SSAA refinement.
 * @details Compares production-resolution pixels after 64, 256, and 640
 * physics substeps. Near-tie Voronoi samples may change, but coverage and
 * high-amplitude errors remain confined to under 3% of lit pixels.
 */
inline void test_gs_shared_stencil_error_is_bounded() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr int probe_frames[] = {4, 16, 40};
  int next_probe = 0;
  for (int frame = 1; frame <= probe_frames[2]; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != probe_frames[next_probe])
      continue;
    auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(gs);
    std::printf("GS shared stencil frame=%d: different=%d lit=%d coverage=%d "
                "hard=%d center_changes=%d/%d center_mismatches=%d max=%d "
                "total=%llu\n",
                frame, error.different, error.lit, error.coverage, error.hard,
                error.stencil_changes, error.stencil_samples,
                error.center_mismatches, error.max_channel,
                static_cast<unsigned long long>(error.total_channel));
    HS_EXPECT_EQ(error.center_mismatches, 0);
    HS_EXPECT(error.different * 2 <= error.lit,
              "shared stencil changed over half of lit pixels");
    HS_EXPECT(error.coverage * 50 <= error.lit,
              "shared stencil moved over 2% of lit coverage");
    HS_EXPECT(error.hard * 100 <= error.lit * 3,
              "shared stencil put over 3% of lit pixels past 6.25% error");
    HS_EXPECT_LE(error.max_channel, 22000);
    HS_EXPECT(error.total_channel <=
                  static_cast<uint64_t>(error.lit) * 3u * 384u,
              "shared stencil mean channel error exceeds 3/512 full scale");
    if (next_probe < 2)
      ++next_probe;
  }
}

inline void test_gs_dissolve_frontier_fades_before_clear() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr float phase = 0.5f;
  constexpr uint32_t seed = 0x7b31d2a5u;
  int cleared = -1, fading = -1;
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    float h = GSWhiteBox::dissolve_hash(i, seed);
    if (cleared < 0 && h < phase)
      cleared = i;
    float fade = (h - phase) / GSWhiteBox::DISSOLVE_FADE_FRACTION;
    if (fading < 0 && fade > 0.25f && fade < 0.75f)
      fading = i;
  }
  HS_EXPECT(cleared >= 0 && fading >= 0,
            "dissolve probe did not find clear and fade nodes");
  GSWhiteBox::set_node(gs, cleared, 0, 65535);
  GSWhiteBox::set_node(gs, fading, 0, 65535);
  GSWhiteBox::convert(gs, phase, seed);
  HS_EXPECT_EQ(GSWhiteBox::a_field(gs)[cleared], (uint16_t)65535);
  HS_EXPECT_EQ(GSWhiteBox::b_field(gs)[cleared], (uint16_t)0);
  HS_EXPECT(GSWhiteBox::a_field(gs)[fading] > 0 &&
                GSWhiteBox::a_field(gs)[fading] < 65535,
            "dissolve fade did not ease A toward rest");
  HS_EXPECT(GSWhiteBox::b_field(gs)[fading] > 0 &&
                GSWhiteBox::b_field(gs)[fading] < 65535,
            "dissolve fade cleared B without an intermediate value");
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
 * @brief Pins the optimized float substep to the scalar equation bit-for-bit.
 */
inline void test_gs_substep_matches_scalar_reference() {
  std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
  std::vector<float> gotA(GSWhiteBox::N), gotB(GSWhiteBox::N);
  std::vector<float> refA(GSWhiteBox::N), refB(GSWhiteBox::N);
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    a[i] = static_cast<float>((i * 73) & 65535) * (1.0f / 65535.0f);
    b[i] = static_cast<float>((i * 151 + 17) & 65535) * (1.0f / 65535.0f);
  }
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.037f, 0.061f, 0.019f, 0.013f, 2.3f);
  GSWhiteBox::step_float(gs, a.data(), b.data(), gotA.data(), gotB.data());
  GSWhiteBox::step_float_reference(gs, a.data(), b.data(), refA.data(),
                                   refB.data());
  HS_EXPECT(
      std::memcmp(gotA.data(), refA.data(), gotA.size() * sizeof(float)) == 0,
      "optimized A substep differs from scalar reference");
  HS_EXPECT(
      std::memcmp(gotB.data(), refB.data(), gotB.size() * sizeof(float)) == 0,
      "optimized B substep differs from scalar reference");
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
 *        substeps at a high-diffusion stable setting.
 * @details At this high-diffusion setting the stability product dt·D·|λ|max =
 *          5·0.03·12 = 1.8 ≤ 2, so the scheme must stay bounded. A genuine
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
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.03f, 0.03f, 3.0f);
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
 *          graph as GS, independent of <W,H>, so the small-aspect resolution is
 *          used arbitrarily.
 */
struct BZWhiteBox {
  using BZ = BZReactionDiffusion<SMALL_W, SMALL_H>;
  using Grid = Scan::Shader::SsaaGrid<SMALL_W, SMALL_H>;
  static constexpr int N = BZ::RD_N;

  static Pixel palette_color(const BZ &bz, float t) {
    return bz.palette.get(t).color;
  }
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
  static void perturb(const BZ &bz, uint16_t *nA, uint16_t *nB, uint16_t *nC) {
    bz.perturb_state(nA, nB, nC);
  }
  static int num_perturbations() { return BZ::NUM_PERTURBATIONS; }
  static int perturb_amount() { return BZ::PERTURB_AMOUNT; }
  static void step(BZ &bz, uint16_t *sA, uint16_t *sB, uint16_t *sC) {
    std::vector<float> fA(N), fB(N), fC(N);
    bz.step_physics(sA, sB, sC, fA.data(), fB.data(), fC.data());
  }
  static void set_state(BZ &bz, uint16_t *a, uint16_t *b, uint16_t *c) {
    bz.state.A = a;
    bz.state.B = b;
    bz.state.C = c;
  }
  static int nearest_node(const Vector &v, const Vector *nodes) {
    int nearest = 0;
    float nearest_d2 = BZ::dist2(v, nodes[0]);
    for (int i = 1; i < N; ++i) {
      float d2 = BZ::dist2(v, nodes[i]);
      if (d2 < nearest_d2) {
        nearest = i;
        nearest_d2 = d2;
      }
    }
    return nearest;
  }
  struct CenterError {
    int pixels = 0;
    int refined_seeds = 0;
    int center_mismatches = 0;
  };

  /**
   * @brief Sweeps a <W,H> frame of cubemap-LUT seeds and compares the certified
   *        render-center early-out against the unconditional argmin walk.
   * @details Mirrors GSWhiteBox::shared_shader_error's center check. The seeds
   *          come from seed_face_lut, the way render() draws them, so a
   *          quantized seed that is not already the nearest node actually
   *          reaches the walk; refined_seeds counts those.
   */
  template <int W, int H> static CenterError render_center_error(BZ &bz) {
    ScratchScope guard(scratch_arena_a);
    Vector *world_nodes = bz.orient_lattice();
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    CenterError error;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        Fragment frag;
        frag.pos = pixel_to_vector<W, H>(x, y);
        bz.seed_face_lut(frag);
        int seed = static_cast<int>(frag.v0);
        int reference = BZ::refine_center(frag.pos, world_nodes, seed);
        if (reference != seed)
          ++error.refined_seeds;
        if (BZ::refine_render_center(frag.pos, world_nodes, seed) != reference)
          ++error.center_mismatches;
        ++error.pixels;
      }
    }
    return error;
  }
  static Pixel shade(const BZ &bz, int seed, const Vector &center,
                     const Vector *nodes, const Grid &grid, int x,
                     const Color4 &ca, const Color4 &cb, const Color4 &cc) {
    BZ::FloatRgb fa(ca.color), fb(cb.color), fc(cc.color);
    return bz.shade_pixel(seed, center, nodes, grid, x, fa, fb, fc);
  }
  static Pixel reference_shade(const BZ &bz, int seed, const Vector &center_rv,
                               const Vector *world_nodes, const Grid &grid,
                               int x, const Color4 &ca, const Color4 &cb,
                               const Color4 &cc) {
    int center = BZ::refine_center(center_rv, world_nodes, seed);
    Vector spos[BZ::RD_K + 1];
    uint16_t sa[BZ::RD_K + 1], sb[BZ::RD_K + 1], sc[BZ::RD_K + 1];
    spos[0] = world_nodes[center];
    sa[0] = bz.state.A[center];
    sb[0] = bz.state.B[center];
    sc[0] = bz.state.C[center];
    int k = 1;
    BZ::for_each_neighbor(center, [&](int ni) {
      spos[k] = world_nodes[ni];
      sa[k] = bz.state.A[ni];
      sb[k] = bz.state.B[ni];
      sc[k] = bz.state.C[ni];
      ++k;
    });

    constexpr float INV_SAMPLES = 1.0f / 4.0f;
    Pixel accum(0, 0, 0);
    for (int i = 0; i < 4; ++i) {
      Vector v = grid.at(x, i);
      float tw = 0, wa = 0, wb = 0, wc = 0;
      for (int j = 0; j < BZ::RD_K + 1; ++j)
        BZ::with_wendland_weight(BZ::dist2(v, spos[j]), [&](float w) {
          wa += sa[j] * w;
          wb += sb[j] * w;
          wc += sc[j] * w;
          tw += w;
        });
      if (tw <= BZ::KERNEL_MIN_TOTAL_WEIGHT)
        continue;
      float inv = BZ::Q16_INV / tw;
      float a = wa * inv, b = wb * inv, c = wc * inv;
      float total = a + b + c;
      if (total < BZ::SPECIES_EMPTY_EPS)
        continue;
      float color_inv = 1.0f / total;
      Pixel color(
          static_cast<uint16_t>(
              (ca.color.r * a + cb.color.r * b + cc.color.r * c) * color_inv +
              0.5f),
          static_cast<uint16_t>(
              (ca.color.g * a + cb.color.g * b + cc.color.g * c) * color_inv +
              0.5f),
          static_cast<uint16_t>(
              (ca.color.b * a + cb.color.b * b + cc.color.b * c) * color_inv +
              0.5f));
      accum += color * (hs::clamp(total, 0.0f, 1.0f) * INV_SAMPLES);
    }
    return accum;
  }
};

/** @brief Pins the three legacy BZ species colors. */
inline void test_bz_legacy_palette() {
  BZWhiteBox::BZ bz;
  HS_EXPECT_TRUE(BZWhiteBox::palette_color(bz, 0.0f) == Pixel(36844, 10770, 3));
  HS_EXPECT_TRUE(BZWhiteBox::palette_color(bz, 0.5f) == Pixel(0, 8112, 5753));
  HS_EXPECT_TRUE(BZWhiteBox::palette_color(bz, 1.0f) == Pixel(2059, 0, 9668));
}

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
 * @details Two RNG-agnostic invariants, at the dt = 1 top of the Speed slider
 *          where the nudge is the full PERTURB_AMOUNT. (1) A fully-saturated
 *          field stays fully saturated: every nudge is a saturating add, so a
 *          node already at 65535 cannot wrap to a low value. (2) On a zero
 *          field, each touched entry is a small positive multiple of the nudge
 *          step and stays within [0, 65535], and at least one entry is touched —
 *          whichever nodes the deterministic RNG happens to select.
 */
inline void test_bz_perturb_state_saturates_and_nudges() {
  BZWhiteBox::BZ bz;
  BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 1.0f);
  // (1) Saturation / no-wrap: all rails stay at the rail.
  {
    std::vector<uint16_t> a(BZWhiteBox::N, 65535), b(BZWhiteBox::N, 65535),
        c(BZWhiteBox::N, 65535);
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
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
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
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
 *          equal), failing if the count drifts in either direction. Checked at
 *          both ends of the Speed slider: the nudge magnitude scales with dt but
 *          the draw count must not, or a paused sphere would desynchronise every
 *          downstream effect.
 */
inline void test_bz_perturb_state_draw_count_pinned() {
  const int expected_draws = 2 * BZWhiteBox::num_perturbations();

  constexpr uint64_t SEED = 1337u;
  for (float dt : {1.0f, 0.0f}) {
    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, dt);
    hs::random().seed(SEED);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());

    // A private generator from the same seed, advanced by the contracted count,
    // must now be at the same stream position as the global generator.
    hs::Pcg32 ref(SEED);
    for (int i = 0; i < expected_draws; ++i)
      (void)ref();
    HS_EXPECT_EQ(hs::random()(), ref());
    // Off-by-one in either direction would have landed at a different output.
    HS_EXPECT_EQ(hs::random()(), ref());
  }
}

/**
 * @brief Verifies the stochastic nudge scales with the Speed slider and reaches
 *        zero where the integrator freezes.
 * @details The nudge is the only integrator term that was not multiplied by dt,
 *          so at dt = 0 — where advance_species round-trips its Q16 sample
 *          exactly — it was the one term still moving the state, walking every
 *          drawn node monotonically to the 65535 rail. Runs many perturbation
 *          passes over a zero field: at dt = 0 nothing may move, and at a
 *          mid-slider dt the touched entries must be multiples of the scaled
 *          step, strictly between zero and the full-rate step.
 */
inline void test_bz_perturb_scales_with_timestep() {
  constexpr int PASSES = 64;
  {
    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, /*dt*/ 0.0f);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    for (int p = 0; p < PASSES; ++p)
      BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int moved = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      if (a[i] != 0 || b[i] != 0 || c[i] != 0)
        ++moved;
    HS_EXPECT_EQ(moved, 0);
  }
  {
    constexpr float DT = 0.35f;
    const int full = BZWhiteBox::perturb_amount();
    const int step = static_cast<int>(full * DT);
    HS_EXPECT_GT(step, 0);
    HS_EXPECT_GT(full, step);

    BZWhiteBox::BZ bz;
    BZWhiteBox::set_params(bz, /*alpha*/ 3.0f, /*D*/ 0.05f, DT);
    std::vector<uint16_t> a(BZWhiteBox::N, 0), b(BZWhiteBox::N, 0),
        c(BZWhiteBox::N, 0);
    for (int p = 0; p < PASSES; ++p)
      BZWhiteBox::perturb(bz, a.data(), b.data(), c.data());
    int touched = 0, malformed = 0;
    for (int i = 0; i < BZWhiteBox::N; ++i)
      for (uint16_t v : {a[i], b[i], c[i]}) {
        if (v == 0)
          continue;
        ++touched;
        if (v % step != 0)
          ++malformed;
      }
    HS_EXPECT_GT(touched, 0);
    HS_EXPECT_EQ(malformed, 0);
  }
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

/**
 * @brief Pins the optimized BZ raster against its scalar sampling contract.
 * @details reference_shade composites the four SSAA samples the naive way: it
 *          normalizes each sample to concentrations, blends the palette into a
 *          uint16 Pixel, premultiplies by that sample's coverage, and adds the
 *          result into a uint16 accumulator. shade_pixel fuses the same algebra
 *          into float species coefficients and quantizes once at the end, so
 *          the two agree only up to the reference's intermediate rounding:
 *          four palette-blend roundings weighted by coverages that sum to at
 *          most 1 (<= 0.5 LSB), four premultiply roundings (<= 2.0 LSB), and
 *          the single final rounding both paths pay (<= 0.5 LSB). Any channel
 *          past that 3 LSB envelope is a formula difference, not round-off.
 */
inline void test_bz_raster_matches_reference() {
  using WhiteBox = BZWhiteBox;
  constexpr int MAX_ROUNDING_DRIFT = 3;
  auto drift = [](uint16_t p, uint16_t q) { return p > q ? p - q : q - p; };
  std::vector<Vector> nodes(WhiteBox::N);
  std::vector<uint16_t> a(WhiteBox::N), b(WhiteBox::N), c(WhiteBox::N);
  for (int i = 0; i < WhiteBox::N; ++i) {
    nodes[i] = ReactionGraph::node(i);
    a[i] = static_cast<uint16_t>((i * 4051u + 123u) & 0xffffu);
    b[i] = static_cast<uint16_t>((i * 7919u + 4567u) & 0xffffu);
    c[i] = static_cast<uint16_t>((i * 10429u + 8901u) & 0xffffu);
  }

  WhiteBox::BZ bz;
  WhiteBox::set_state(bz, a.data(), b.data(), c.data());
  if (!TrigLUT<SMALL_W, SMALL_H>::initialized)
    TrigLUT<SMALL_W, SMALL_H>::init();
  WhiteBox::Grid grid;
  const Color4 ca(Pixel(61123, 913, 17771), 1.0f);
  const Color4 cb(Pixel(2819, 59731, 1207), 1.0f);
  const Color4 cc(Pixel(4513, 7781, 62927), 1.0f);

  auto compare = [&] {
    int mismatches = 0;
    for (int y = 0; y < SMALL_H; y += 2) {
      grid.set_row(y);
      for (int x = 0; x < SMALL_W; x += 3) {
        Vector center = pixel_to_vector<SMALL_W, SMALL_H>(x, y);
        int seed = WhiteBox::nearest_node(center, nodes.data());
        Pixel expected = WhiteBox::reference_shade(
            bz, seed, center, nodes.data(), grid, x, ca, cb, cc);
        Pixel actual = WhiteBox::shade(bz, seed, center, nodes.data(), grid, x,
                                       ca, cb, cc);
        if (drift(actual.r, expected.r) > MAX_ROUNDING_DRIFT ||
            drift(actual.g, expected.g) > MAX_ROUNDING_DRIFT ||
            drift(actual.b, expected.b) > MAX_ROUNDING_DRIFT)
          ++mismatches;
      }
    }
    return mismatches;
  };

  int mismatches = compare();
  for (int i = 0; i < WhiteBox::N; ++i) {
    a[i] >>= 6;
    b[i] >>= 6;
    c[i] >>= 6;
  }
  mismatches += compare();
  HS_EXPECT_EQ(mismatches, 0);
}

/**
 * @brief Requires BZ's certified render-center early-out to agree with the
 *        unconditional argmin on production-resolution, LUT-seeded pixels.
 * @details shade_pixel centers its stencil with refine_render_center, which
 *          returns the seed unwalked inside BULK_CERTIFIED_D2 /
 *          POLE_CERTIFIED_D2. Only a cubemap-LUT seed — what render()'s vertex
 *          shader produces — can be a non-nearest node, so a sweep seeded from
 *          an exact argmin cannot tell the two refiners apart at all. A few
 *          frames of the orientation random walk move the lattice off its
 *          initial alignment first. refined_seeds must be nonzero or the
 *          comparison passed without ever exercising the walk.
 */
inline void test_bz_render_center_matches_reference() {
  hs_test::reset_globals();
  BZWhiteBox::BZ bz;
  bz.init();
  for (int frame = 0; frame < 4; ++frame) {
    bz.draw_frame();
    bz.advance_display();
  }
  auto error = BZWhiteBox::render_center_error<DEFAULT_W, DEFAULT_H>(bz);
  std::printf("BZ render centers: pixels=%d refined_seeds=%d mismatches=%d\n",
              error.pixels, error.refined_seeds, error.center_mismatches);
  HS_EXPECT_EQ(error.pixels, DEFAULT_W * DEFAULT_H);
  HS_EXPECT_GT(error.refined_seeds, 0);
  HS_EXPECT_EQ(error.center_mismatches, 0);
}

// ---------------------------------------------------------------------------
// DreamBalls: preset-cycle / re-spawn white-box coverage
// ---------------------------------------------------------------------------

/** @brief Test access to HankinSolids' morph-chain state. */
struct HankinPauseWhiteBox {
  using EffectT = HankinSolids<SMALL_W, SMALL_H>;

  static bool morph_active(const EffectT &effect) {
    return effect.pending_landing != nullptr;
  }
  static uint8_t node(const EffectT &effect) { return effect.node; }
};

/** @brief Verifies a manual Angle edit freezes an in-flight Conway morph. */
inline void test_hankinsolids_manual_pause_holds_morph() {
  reset_effect_globals();
  HankinPauseWhiteBox::EffectT effect;
  effect.init();

  int frames = 0;
  while (!HankinPauseWhiteBox::morph_active(effect) && frames < 100) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
  }
  HS_EXPECT_TRUE(HankinPauseWhiteBox::morph_active(effect));

  for (int i = 0; i < 5; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  const uint8_t held_node = HankinPauseWhiteBox::node(effect);
  HS_EXPECT_TRUE(effect.updateParameter("Angle", 0.7f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.animations_paused());

  for (int i = 0; i < 300; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(HankinPauseWhiteBox::morph_active(effect));
  HS_EXPECT_EQ(HankinPauseWhiteBox::node(effect), held_node);
  const auto *angle = effect.getParameters().find("Angle");
  HS_EXPECT_TRUE(angle != nullptr);
  HS_EXPECT_NEAR(angle->get(), 0.7f, 1e-6f);

  uint64_t energy = 0;
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief White-box accessor for DreamBalls' private preset-cycle bookkeeping
 *        (befriended in effects/DreamBalls.h).
 * @details spawn_sprite schedules its successor 320 frames out (a PeriodicTimer),
 *          but the smoke/determinism harness renders at most HS_SMOKE_FRAMES=120
 *          frames in CI — short of one period — so the re-spawn never fires under
 *          the generic passes: the preset advance, the active_bake ping-pong +
 *          rebake, and the reseed-on-change guard all stay dead. This seam drives
 *          spawn_sprite directly and reads the bake slot / preset index so those
 *          paths are pinned. <96,20> is used arbitrarily — the bookkeeping is
 *          resolution-independent.
 */
struct DreamBallsWhiteBox {
  using DB = DreamBalls<SMALL_W, SMALL_H>;
  static constexpr int PRESETS = 10;
  static constexpr size_t SOLID_COUNT = DB::SOLID_COUNT;
  static constexpr size_t MAX_SOLID_EDGES = Solids::MAX_SOLID_EDGES;
  static constexpr size_t SCRATCH_A_PEAK_BYTES = DB::SCRATCH_A_PEAK_BYTES;

  static int active_bake(const DB &db) { return db.active_bake; }
  static int last_preset_idx(const DB &db) { return db.last_preset_idx; }
  // Not-paused step: advance the selector, then re-spawn (the scheduler's path).
  static void advance(DB &db) {
    db.preset_manager.next();
    db.spawn_sprite();
  }
  // Re-spawn of the current preset (no advance).
  static void respawn(DB &db) { db.spawn_sprite(); }
  static DB::BaseMesh preset_mesh(const DB &db, int idx) {
    return db.preset_manager.get_entries()[idx].params.base_mesh;
  }
  static const DB::Params &preset_params(const DB &db, int idx) {
    return db.preset_manager.get_entries()[idx].params;
  }
  static const Palette *blood_stream_falloff(const DB &db) {
    return &db.blood_stream_falloff;
  }
  static bool solid_loaded(const DB &db, size_t idx) {
    return !db.loaded_solids[idx].mesh_state.vertices.is_empty();
  }
  static bool four_regular(const DB &db, size_t idx) {
    return db.loaded_solids[idx].four_regular;
  }
  static const auto &original_edges(const DB &db, size_t idx) {
    return db.loaded_solids[idx].original_edges;
  }
  static const auto &automatic_edges(const DB &db, size_t idx) {
    return db.loaded_solids[idx].automatic_edges;
  }
  static const auto &medial_edges(const DB &db, size_t idx) {
    const auto &solid = db.loaded_solids[idx];
    return solid.four_regular ? solid.medial_edges : solid.automatic_edges;
  }
  static size_t source_vertex_count(const DB &db, size_t idx) {
    return db.loaded_solids[idx].mesh_state.vertices.size();
  }
  static float under_gap_alpha(float edge_t, float gap) {
    return DB::under_gap_alpha(edge_t, gap);
  }
  static Vector woven_vertex(const DB &db, size_t solid, bool medial,
                             size_t vertex) {
    return DB::woven_vertex(db.loaded_solids[solid], medial, vertex);
  }
  static Vector tangent_axis(const Vector &normal) {
    return DB::tangent_axis(normal);
  }
  static Vector parallel_transport(const Vector &from, const Vector &to,
                                   const Vector &tangent) {
    return DB::parallel_transport(from, to, tangent);
  }
  static DB::BaseMesh live_mesh(const DB &db) { return db.params.base_mesh; }
  static DB::WeaveTopology live_weave_topology(const DB &db) {
    return db.params.weave_topology;
  }
  static float live_weave_gap(const DB &db) { return db.params.weave_gap; }
  static float &num_copies(DB &db) { return db.params.num_copies; }
};

/**
 * @brief Drives spawn_sprite across a full preset cycle and asserts the bake-slot
 *        ping-pong, the modulo preset advance, and the reseed-on-change guard.
 * @details Drives the advance without the 320-frame wait, following the same
 *          selector progression the periodic callback uses: each step calls
 *          Presets::next() then re-spawns, so the active preset walks modulo
 *          the preset count.
 *          Each spawn must flip the bake slot (so a fading-out sprite keeps its
 *          own LUT) and, when the preset actually changes, reseed params to the
 *          new entry. A re-spawn of the SAME preset must instead hold params
 *          so a live slider edit survives.
 */
inline void test_dreamballs_preset_cycle_bookkeeping() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init(); // runs spawn_sprite() at preset 0

  // init() spawned preset 0: it reseeded params (last_preset_idx -1 -> 0) and
  // flipped the bake slot once (0 -> 1).
  HS_EXPECT_EQ(WB::last_preset_idx(db), 0);
  HS_EXPECT_EQ(WB::active_bake(db), 1);
  HS_EXPECT_EQ(WB::live_mesh(db), WB::preset_mesh(db, 0));

  const auto &snub_cube = WB::preset_params(db, 4);
  HS_EXPECT_EQ(snub_cube.base_mesh, WB::DB::BaseMesh::SNUB_CUBE);
  HS_EXPECT_EQ(snub_cube.weave_topology, WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(snub_cube.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(snub_cube.num_copies, 4.534f, 1e-6f);
  HS_EXPECT_NEAR(snub_cube.offset_radius, 0.153f, 1e-6f);
  HS_EXPECT_NEAR(snub_cube.offset_speed, 2.025f, 1e-6f);
  HS_EXPECT_NEAR(snub_cube.warp_scale, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(snub_cube.alpha, 0.3f, 1e-6f);

  const auto &truncated_dodecahedron = WB::preset_params(db, 5);
  HS_EXPECT_EQ(truncated_dodecahedron.base_mesh,
               WB::DB::BaseMesh::TRUNCATED_DODECAHEDRON);
  HS_EXPECT_EQ(truncated_dodecahedron.weave_topology,
               WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(truncated_dodecahedron.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(truncated_dodecahedron.num_copies, 4.515f, 1e-6f);
  HS_EXPECT_NEAR(truncated_dodecahedron.offset_radius, 0.179f, 1e-6f);
  HS_EXPECT_NEAR(truncated_dodecahedron.offset_speed, 1.89f, 1e-6f);
  HS_EXPECT_NEAR(truncated_dodecahedron.warp_scale, 1.535f, 1e-6f);
  HS_EXPECT_NEAR(truncated_dodecahedron.alpha, 0.7f, 1e-6f);

  const auto &triakis_icosahedron = WB::preset_params(db, 6);
  HS_EXPECT_EQ(triakis_icosahedron.base_mesh,
               WB::DB::BaseMesh::TRIAKIS_ICOSAHEDRON);
  HS_EXPECT_EQ(triakis_icosahedron.weave_topology,
               WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(triakis_icosahedron.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron.num_copies, 4.515f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron.offset_radius, 0.131f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron.offset_speed, 1.89f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron.warp_scale, 1.535f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron.alpha, 0.7f, 1e-6f);

  const auto &triakis_icosahedron_six_copies = WB::preset_params(db, 7);
  HS_EXPECT_EQ(triakis_icosahedron_six_copies.base_mesh,
               WB::DB::BaseMesh::TRIAKIS_ICOSAHEDRON);
  HS_EXPECT_EQ(triakis_icosahedron_six_copies.weave_topology,
               WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.num_copies, 6.0f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.offset_radius, 0.078f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.offset_speed, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.warp_scale, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_six_copies.alpha, 0.3f, 1e-6f);

  const auto &disdyakis_triacontahedron = WB::preset_params(db, 8);
  HS_EXPECT_EQ(disdyakis_triacontahedron.base_mesh,
               WB::DB::BaseMesh::DISDYAKIS_TRIACONTAHEDRON);
  HS_EXPECT_EQ(disdyakis_triacontahedron.weave_topology,
               WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.num_copies, 6.0f, 1e-6f);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.offset_radius, 0.03f, 1e-6f);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.offset_speed, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.warp_scale, 1.795f, 1e-6f);
  HS_EXPECT_NEAR(disdyakis_triacontahedron.alpha, 0.3f, 1e-6f);

  const auto &triakis_icosahedron_compact = WB::preset_params(db, 9);
  HS_EXPECT_EQ(triakis_icosahedron_compact.base_mesh,
               WB::DB::BaseMesh::TRIAKIS_ICOSAHEDRON);
  HS_EXPECT_EQ(triakis_icosahedron_compact.weave_topology,
               WB::DB::WeaveTopology::AUTOMATIC);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.weave_gap, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.num_copies, 6.0f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.offset_radius, 0.03f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.offset_speed, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.warp_scale, 1.795f, 1e-6f);
  HS_EXPECT_NEAR(triakis_icosahedron_compact.alpha, 0.3f, 1e-6f);

  HS_EXPECT_TRUE(WB::preset_params(db, 0).palette ==
                 WB::blood_stream_falloff(db));
  HS_EXPECT_TRUE(WB::preset_params(db, 1).palette ==
                 WB::blood_stream_falloff(db));
  HS_EXPECT_TRUE(WB::preset_params(db, 2).palette == &Palettes::RICH_SUNSET);
  HS_EXPECT_TRUE(WB::preset_params(db, 3).palette == &Palettes::LAVENDER_LAKE);
  HS_EXPECT_TRUE(WB::preset_params(db, 4).palette == &Palettes::CORAL_BLUE);
  HS_EXPECT_TRUE(WB::preset_params(db, 5).palette == &Palettes::CORAL_BLUE);
  HS_EXPECT_TRUE(WB::preset_params(db, 6).palette == &Palettes::CORAL_BLUE);
  HS_EXPECT_TRUE(WB::preset_params(db, 7).palette == &Palettes::CORAL_BLUE);
  HS_EXPECT_TRUE(WB::preset_params(db, 8).palette == &Palettes::CORAL_BLUE);
  HS_EXPECT_TRUE(WB::preset_params(db, 9).palette == &Palettes::CORAL_BLUE);

  // Not-paused advance chain: each step advances the selector then re-spawns, so
  // the preset is step modulo the preset count. Drive two full cycles; the bake
  // slot must ping-pong and params must reseed to the new index each step.
  int expect_bake = WB::active_bake(db); // 1
  for (int step = 1; step <= 2 * WB::PRESETS; ++step) {
    WB::advance(db);
    expect_bake ^= 1;
    const int safe = step % WB::PRESETS;
    HS_EXPECT_EQ(WB::active_bake(db), expect_bake);
    HS_EXPECT_EQ(WB::last_preset_idx(db), safe);
    HS_EXPECT_EQ(WB::live_mesh(db), WB::preset_mesh(db, safe));
  }

  // Same-preset path: re-spawn without advancing, so the reseed guard
  // (safe_idx == last_preset_idx) holds and a live slider edit must survive the
  // re-spawn — while the bake slot still flips. last_preset_idx is now 0.
  const float sentinel = WB::num_copies(db) + 5.0f;
  WB::num_copies(db) = sentinel;
  const int held_idx = WB::last_preset_idx(db);
  expect_bake ^= 1;
  WB::respawn(db);
  HS_EXPECT_EQ(WB::last_preset_idx(db), held_idx); // preset unchanged
  HS_EXPECT_EQ(WB::num_copies(db), sentinel);      // live edit preserved
  HS_EXPECT_EQ(WB::active_bake(db), expect_bake);  // bake slot still flipped
}

/** @brief Verifies the Base Mesh dropdown covers every simple solid family. */
inline void test_dreamballs_base_mesh_selector() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  const auto *base_mesh = db.getParameters().find("Base Mesh");
  HS_EXPECT_TRUE(base_mesh != nullptr);
  if (!base_mesh)
    return;

  constexpr int EXPECTED_SOLIDS =
      static_cast<int>(Solids::PLATONIC_COUNT + Solids::ARCHIMEDEAN_COUNT +
                       Solids::CATALAN_COUNT);
  HS_EXPECT_EQ(base_mesh->option_count, EXPECTED_SOLIDS);
  HS_EXPECT_TRUE(base_mesh->animated);
  for (int i = 0; i < EXPECTED_SOLIDS; ++i)
    HS_EXPECT_TRUE(WB::solid_loaded(db, static_cast<size_t>(i)));
  HS_EXPECT_EQ(std::string_view(base_mesh->options[0]), "Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[17]), "Snub Dodecahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[18]), "Triakis Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[EXPECTED_SOLIDS - 1]),
               "Pentagonal Hexecontahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->export_options[EXPECTED_SOLIDS - 1]),
               "BaseMesh::PENTAGONAL_HEXECONTAHEDRON");

  HS_EXPECT_TRUE(db.updateParameter("Base Mesh",
                                    static_cast<float>(EXPECTED_SOLIDS - 1)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_mesh(db), WB::DB::BaseMesh::PENTAGONAL_HEXECONTAHEDRON);
  HS_EXPECT_TRUE(db.animations_paused());

  db.setAnimationsPaused(false);
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }
  uint64_t energy = 0;
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &pixel = db.get_pixel(x, y);
      energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief Renders the widest solid the dropdown reaches under medial topology and
 *        checks SCRATCH_A_PEAK_BYTES against the real frame peak.
 * @details No preset selects a MAX_SOLID_EDGES solid, and the smoke sweep stops
 *          well short of SPRITE_LIFE, so the widest woven staging the effect can
 *          bind is otherwise never drawn. Medial topology is forced so the
 *          staging takes one vertex per source edge and one framed edge per
 *          medial edge — the worst case the static_assert bounds.
 */
inline void test_dreamballs_max_edge_solid_render() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  size_t widest = 0;
  size_t widest_edges = 0;
  for (size_t i = 0; i < WB::SOLID_COUNT; ++i) {
    const size_t edges = WB::original_edges(db, i).size();
    if (edges > widest_edges) {
      widest_edges = edges;
      widest = i;
    }
  }
  HS_EXPECT_EQ(widest_edges, WB::MAX_SOLID_EDGES);

  HS_EXPECT_TRUE(db.updateParameter("Base Mesh", static_cast<float>(widest)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(db.updateParameter("Weave Topology", 2.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db), WB::DB::WeaveTopology::MEDIAL);
  HS_EXPECT_EQ(WB::medial_edges(db, widest).size(), 2 * widest_edges);

  db.setAnimationsPaused(false);
  scratch_arena_a.reset_high_water_mark();
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }

  // Three per-vertex buffers plus the framed vertex + edge-head mesh, at the
  // medial bound; a peak below this would mean the wide path never ran.
  const size_t staged_bytes = 6 * widest_edges * sizeof(Vector);
  HS_EXPECT_GE(scratch_arena_a.get_high_water_mark(), staged_bytes);
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), WB::SCRATCH_A_PEAK_BYTES);

  uint64_t energy = 0;
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &pixel = db.get_pixel(x, y);
      energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(energy, 0u);
}

/** @brief Verifies automatic, forced-medial, and defect weave topology. */
inline void test_dreamballs_weave_topology() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  const auto *topology = db.getParameters().find("Weave Topology");
  const auto *gap = db.getParameters().find("Weave Gap");
  HS_EXPECT_TRUE(topology != nullptr);
  HS_EXPECT_TRUE(gap != nullptr);
  if (!topology || !gap)
    return;

  HS_EXPECT_EQ(topology->option_count, 3);
  HS_EXPECT_EQ(std::string_view(topology->options[0]), "Automatic");
  HS_EXPECT_EQ(std::string_view(topology->options[1]), "Original with defects");
  HS_EXPECT_EQ(std::string_view(topology->options[2]), "Medial");
  HS_EXPECT_EQ(std::string_view(topology->export_options[1]),
               "WeaveTopology::ORIGINAL_WITH_DEFECTS");
  HS_EXPECT_EQ(std::string_view(topology->export_options[2]),
               "WeaveTopology::MEDIAL");
  HS_EXPECT_TRUE(topology->animated);
  HS_EXPECT_TRUE(gap->animated);

  constexpr int EXPECTED_SOLIDS =
      static_cast<int>(Solids::PLATONIC_COUNT + Solids::ARCHIMEDEAN_COUNT +
                       Solids::CATALAN_COUNT);
  int four_regular_count = 0;
  for (int i = 0; i < EXPECTED_SOLIDS; ++i) {
    const auto &original = WB::original_edges(db, i);
    const auto &automatic = WB::automatic_edges(db, i);
    const auto &medial = WB::medial_edges(db, i);
    const bool four_regular = WB::four_regular(db, i);
    const bool expected_four_regular =
        i == static_cast<int>(WB::DB::BaseMesh::OCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::CUBOCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::RHOMBICUBOCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::ICOSIDODECAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::RHOMBICOSIDODECAHEDRON);
    HS_EXPECT_EQ(four_regular, expected_four_regular);
    four_regular_count += four_regular ? 1 : 0;

    const size_t automatic_vertex_count =
        four_regular ? WB::source_vertex_count(db, i) : original.size();
    HS_EXPECT_EQ(automatic.size(),
                 four_regular ? original.size() : 2 * original.size());
    HS_EXPECT_EQ(medial.size(), 2 * original.size());

    std::vector<int> incoming(automatic_vertex_count, 0);
    std::vector<int> outgoing(automatic_vertex_count, 0);
    for (const auto &edge : automatic) {
      HS_EXPECT_LT(edge.u, automatic_vertex_count);
      HS_EXPECT_LT(edge.v, automatic_vertex_count);
      if (edge.u < automatic_vertex_count && edge.v < automatic_vertex_count) {
        outgoing[edge.u]++;
        incoming[edge.v]++;

        const Vector from = WB::woven_vertex(db, i, !four_regular, edge.u);
        const Vector to = WB::woven_vertex(db, i, !four_regular, edge.v);
        const Vector frame_u = WB::tangent_axis(from);
        const Vector offset = frame_u * 0.6f + cross(from, frame_u) * 0.8f;
        const Vector transported = WB::parallel_transport(from, to, offset);
        HS_EXPECT_NEAR(dot(transported, to), 0.0f, 2e-5f);
        HS_EXPECT_NEAR(dot(transported, transported), 1.0f, 2e-5f);
        HS_EXPECT_VEC(WB::parallel_transport(to, from, transported), offset,
                      2e-5f);
      }
    }
    for (size_t vertex = 0; vertex < automatic_vertex_count; ++vertex) {
      HS_EXPECT_EQ(incoming[vertex], 2);
      HS_EXPECT_EQ(outgoing[vertex], 2);
    }

    std::fill(incoming.begin(), incoming.end(), 0);
    std::fill(outgoing.begin(), outgoing.end(), 0);
    incoming.resize(original.size(), 0);
    outgoing.resize(original.size(), 0);
    for (const auto &edge : medial) {
      HS_EXPECT_LT(edge.u, original.size());
      HS_EXPECT_LT(edge.v, original.size());
      if (edge.u < original.size() && edge.v < original.size()) {
        outgoing[edge.u]++;
        incoming[edge.v]++;
      }
    }
    for (size_t vertex = 0; vertex < original.size(); ++vertex) {
      HS_EXPECT_EQ(incoming[vertex], 2);
      HS_EXPECT_EQ(outgoing[vertex], 2);
    }
  }
  HS_EXPECT_EQ(four_regular_count, 5);

  HS_EXPECT_NEAR(WB::under_gap_alpha(0.0f, 0.2f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(0.8f, 0.2f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(0.9f, 0.2f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(1.0f, 0.2f), 0.0f, 1e-6f);

  HS_EXPECT_TRUE(db.updateParameter("Base Mesh", 0.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(db.updateParameter("Weave Topology", 1.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(db.updateParameter("Weave Gap", 0.25f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db),
               WB::DB::WeaveTopology::ORIGINAL_WITH_DEFECTS);
  HS_EXPECT_NEAR(WB::live_weave_gap(db), 0.25f, 1e-6f);

  HS_EXPECT_TRUE(
      db.updateParameter("Base Mesh",
                         static_cast<float>(WB::DB::BaseMesh::OCTAHEDRON)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(db.updateParameter("Weave Topology", 2.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db), WB::DB::WeaveTopology::MEDIAL);

  db.setAnimationsPaused(false);
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }
  uint64_t energy = 0;
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &pixel = db.get_pixel(x, y);
      energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief Verifies pause freezes DreamBalls' sprite and next-preset clock.
 */
inline void test_dreamballs_respawn_fires_and_honors_pause() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();
  WB::DB db;
  db.init();
  const int held_bake = WB::active_bake(db);
  db.setAnimationsPaused(true);
  for (int f = 0; f < 340; ++f) {
    db.draw_frame();
    db.advance_display();
  }
  HS_EXPECT_EQ(WB::last_preset_idx(db), 0);
  HS_EXPECT_EQ(WB::active_bake(db), held_bake);

  db.setAnimationsPaused(false);
  for (int f = 0; f < 340; ++f) {
    db.draw_frame();
    db.advance_display();
  }
  HS_EXPECT_EQ(WB::last_preset_idx(db), 1);
  HS_EXPECT_EQ(WB::active_bake(db), held_bake ^ 1);
}

// ---------------------------------------------------------------------------
// MeshFeedback: draw_frame block-ordering coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for MeshFeedback's style / noise / preset state.
 * @details Befriended in effects/MeshFeedback.h. draw_frame runs four blocks
 *          whose ORDER is the contract — the preset switch leads apply_params()
 *          so the flush reads one preset's scalars, and the flush leads the mesh
 *          draw so the frame's own wireframe survives it. Reordering them is the
 *          cheapest possible refactor, changes the image, and is invisible to the
 *          roster sweeps: the PRESET_DWELL_FRAMES rotation never fires inside the
 *          smoke window at all.
 */
struct MeshFeedbackWhiteBox {
  using MF = MeshFeedback<SMALL_W, SMALL_H>;

  static const Feedback::Style &style(const MF &fx) { return fx.params.style; }
  static MF::BaseMesh base_mesh(const MF &fx) { return fx.params.base_mesh; }
  static MF::BaseMesh active_base_mesh(const MF &fx) {
    return fx.active_base_mesh;
  }
  static MF::BaseMesh preset_base_mesh(MF &fx, size_t index) {
    return fx.preset_params(index).base_mesh;
  }
  static size_t mesh_vertices(const MF &fx) { return fx.mesh.vertices.size(); }
  static size_t mesh_edges(const MF &fx) { return fx.edges.size(); }
  static const Animation::NoiseParams &noise(const MF &fx) {
    return fx.noise_params;
  }
  static size_t preset_index(const MF &fx) { return fx.getPresetIndex(); }
  static size_t mesh_storage_mark(const MF &fx) { return fx.mesh_storage_mark; }
  static void rebuild_mesh(MF &fx, MF::BaseMesh base_mesh) {
    fx.rebuild_mesh(base_mesh);
  }
};

/**
 * @brief Renders MeshFeedback for `frames` frames and copies out the last frame.
 * @param out Receives the displayed frame, SMALL_W * SMALL_H pixels row-major.
 * @param frames Number of frames to render.
 * @param feedback Value written to the Feedback toggle before the first frame.
 * @details Resets the same globals render_capture() does, so a feedback-on and a
 *          feedback-off run draw the identical wireframe on the identical frame.
 */
inline void meshfeedback_capture(std::vector<Pixel> &out, int frames,
                                 bool feedback) {
  reset_effect_globals();
  hs::set_mock_time(0, 0);

  MeshFeedbackWhiteBox::MF fx;
  fx.init();
  HS_EXPECT_TRUE(fx.updateParameter("Feedback", feedback ? 1.0f : 0.0f) ==
                 ParamSetResult::APPLIED);
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();
  }
  hs::clear_mock_time();

  capture_frame<SMALL_W, SMALL_H>(fx, out);
}

/**
 * @brief Verifies the feedback flush never decays the same frame's wireframe.
 * @details The flush composites the warped previous frame at alpha 1, i.e. it
 *          OVERWRITES the draw buffer; running it after the mesh draw would
 *          erase the frame's own wireframe and leave nothing but exponentially
 *          decaying history. Compares a feedback-on render against a
 *          feedback-off one of the same frame: with the flush leading, every
 *          channel of the on-render is at least the off-render's (the wireframe
 *          is laid over a brighter-or-equal background), and the trails make it
 *          strictly greater somewhere.
 */
inline void test_meshfeedback_flush_precedes_mesh_draw() {
  // Well past the empty-history early-out, short of the preset rotation.
  constexpr int FRAMES = 48;
  std::vector<Pixel> lit, bare;
  meshfeedback_capture(lit, FRAMES, true);
  meshfeedback_capture(bare, FRAMES, false);

  int decayed = 0, accumulated = 0;
  uint64_t bare_energy = 0;
  for (size_t i = 0; i < lit.size(); ++i) {
    const Pixel &a = lit[i], &b = bare[i];
    if (a.r < b.r || a.g < b.g || a.b < b.b)
      ++decayed;
    if (a.r > b.r || a.g > b.g || a.b > b.b)
      ++accumulated;
    bare_energy += static_cast<uint64_t>(b.r) + b.g + b.b;
  }

  std::printf("  MeshFeedback flush order: decayed=%d accumulated=%d "
              "wireframe energy=%llu\n",
              decayed, accumulated,
              static_cast<unsigned long long>(bare_energy));
  // A pair of black frames would agree trivially.
  HS_EXPECT_GT(bare_energy, 0u);
  HS_EXPECT_EQ(decayed, 0);
  HS_EXPECT_GT(accumulated, 0);
}

/**
 * @brief Drives the preset rotation and pins the switch-frame noise sync.
 * @details PRESET_DWELL_FRAMES is 241, so the rotation is out of reach of every
 *          roster sweep. Crosses two boundaries and requires the selector to
 *          step one entry each time — the second boundary only lands on frame
 *          482 if the first switch re-armed the dwell — and, the ordering
 *          contract, requires the bound NoiseParams to already carry the
 *          incoming preset's scalars when the switch frame ends: apply_params()
 *          runs after begin_automatic_transition(), so the flush that frame
 *          reads one preset's fade and noise, not two.
 */
inline void test_meshfeedback_preset_rotation_syncs_noise() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();
  hs::set_mock_time(0, 0);

  MF fx;
  fx.init();
  HS_EXPECT_EQ(WB::preset_index(fx), 0u);

  const auto in_sync = [&fx]() {
    const Feedback::Style &s = WB::style(fx);
    const Animation::NoiseParams &n = WB::noise(fx);
    return n.amplitude == s.amplitude && n.frequency == s.frequency &&
           n.speed == s.speed && n.scale == s.scale;
  };
  const auto matches_preset = [&fx](size_t idx) {
    const Feedback::Style &s = WB::style(fx);
    const Feedback::Style &p = MF::PRESETS[idx].params.style;
    return s.fade == p.fade && s.hue_shift == p.hue_shift &&
           s.scale == p.scale && s.amplitude == p.amplitude &&
           s.frequency == p.frequency && s.speed == p.speed &&
           s.space_fn == p.space_fn && s.color_fn == p.color_fn;
  };

  int switches = 0, desynced = 0, wrong_preset = 0;
  for (int f = 1; f <= 2 * MF::PRESET_DWELL_FRAMES + 1; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();

    const size_t expect_idx =
        static_cast<size_t>(f / MF::PRESET_DWELL_FRAMES) % MF::PRESETS.size();
    if (WB::preset_index(fx) != expect_idx)
      ++wrong_preset;
    if (f % MF::PRESET_DWELL_FRAMES == 0) {
      ++switches;
      // The live style is the incoming preset, and the noise the flush reads
      // already followed it within this same frame.
      HS_EXPECT_TRUE(matches_preset(expect_idx));
      HS_EXPECT_TRUE(in_sync());
    }
    if (!in_sync())
      ++desynced;
  }
  hs::clear_mock_time();

  std::printf("  MeshFeedback presets: %d switches, %d desynced frames, "
              "%d wrong index\n",
              switches, desynced, wrong_preset);
  HS_EXPECT_EQ(switches, 2);
  HS_EXPECT_EQ(wrong_preset, 0);
  HS_EXPECT_EQ(desynced, 0);
}

/** @brief Verifies MeshFeedback presets and GUI carry the selected solid. */
inline void test_meshfeedback_base_mesh_selector() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  const auto *base_mesh = effect.getParameters().find("Base Mesh");
  HS_EXPECT_TRUE(base_mesh != nullptr);
  if (!base_mesh)
    return;

  HS_EXPECT_EQ(base_mesh->option_count,
               static_cast<int>(Solids::BASE_MESH_COUNT));
  HS_EXPECT_TRUE(base_mesh->animated);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 0), MF::BaseMesh::ICOSAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 1), MF::BaseMesh::DODECAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 11),
               MF::BaseMesh::PENTAGONAL_HEXECONTAHEDRON);

  const auto selected = MF::BaseMesh::PENTAGONAL_HEXECONTAHEDRON;
  HS_EXPECT_TRUE(
      effect.updateParameter("Base Mesh", static_cast<float>(selected)) ==
      ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(WB::base_mesh(effect), selected);
  HS_EXPECT_EQ(WB::active_base_mesh(effect), selected);
  HS_EXPECT_GT(WB::mesh_vertices(effect), 0u);
  HS_EXPECT_GT(WB::mesh_edges(effect), 0u);
}

/**
 * @brief Verifies MeshFeedback's preset export covers exactly its Params.
 * @details A preset export emits one value per preset-flagged param into a
 *   PresetEntry<Params> aggregate, so a param backed by a member outside Params
 *   must be marked global or the pasted row will not compile.
 */
inline void test_meshfeedback_preset_export_arity() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  const auto *feedback = effect.getParameters().find("Feedback");
  HS_EXPECT_TRUE(feedback != nullptr);
  if (feedback)
    HS_EXPECT_FALSE(feedback->preset);

  // base_mesh + the six registered Style scalars.
  int preset_params = 0;
  for (const auto &def : effect.getParameters())
    if (def.preset)
      ++preset_params;
  HS_EXPECT_EQ(preset_params, 7);
}

/** @brief Verifies repeated MeshFeedback mesh changes reuse arena storage. */
inline void test_meshfeedback_mesh_rebuild_reuses_storage() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  const size_t storage_mark = WB::mesh_storage_mark(effect);
  std::array<size_t, Solids::BASE_MESH_COUNT> first_cycle_offsets{};
  for (size_t i = 0; i < Solids::BASE_MESH_COUNT; ++i) {
    WB::rebuild_mesh(effect, static_cast<MF::BaseMesh>(i));
    first_cycle_offsets[i] = persistent_arena.get_offset();
    HS_EXPECT_GT(first_cycle_offsets[i], storage_mark);
  }

  for (size_t i = 0; i < Solids::BASE_MESH_COUNT; ++i) {
    WB::rebuild_mesh(effect, static_cast<MF::BaseMesh>(i));
    HS_EXPECT_EQ(persistent_arena.get_offset(), first_cycle_offsets[i]);
  }
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

  using C = Comets<SMALL_W, SMALL_H>;
  static constexpr int wipe_frames() { return C::WIPE_FRAMES; }
  static void roll_palette_over(C &c) { c.update_palette(); }
  static int wipe_frames_remaining(const C &c) {
    return c.wipe.frames_remaining;
  }
  static const GenerativePalette::Snapshot &palette_start(const C &c) {
    return c.wipe.start;
  }
  static const GenerativePalette::Snapshot &palette_target(const C &c) {
    return c.wipe.target;
  }
  static const Quaternion &node_orientation(const C &c) {
    return c.node->orientation.get();
  }
  static size_t trail_length(const C &c) { return c.node->trail.length(); }
  static int cycle_period(const C &c) {
    return 2 * static_cast<int>(c.params.cycle_duration);
  }
};

/**
 * @brief Field-wise equality for two generative-palette snapshots.
 * @param a First snapshot.
 * @param b Second snapshot.
 * @return true when every authored key and axis range matches.
 * @details Compared field-wise rather than by memcmp: Snapshot has padding
 *          between its trailing bytes that carries no meaning.
 */
inline bool palette_snapshots_equal(const GenerativePalette::Snapshot &a,
                                    const GenerativePalette::Snapshot &b) {
  if (a.key_count != b.key_count || a.lightness_low != b.lightness_low ||
      a.lightness_high != b.lightness_high || a.chroma_low != b.chroma_low ||
      a.chroma_high != b.chroma_high ||
      a.lightness_curve != b.lightness_curve ||
      a.chroma_curve != b.chroma_curve)
    return false;
  for (size_t i = 0; i < a.keys.size(); ++i)
    if (a.keys[i].bytes != b.keys[i].bytes)
      return false;
  return true;
}

/**
 * @brief Pins Comets' mid-wipe rollover skip.
 * @details At the Cycle Dur floor the rollover timer's period is shorter than
 *          WIPE_FRAMES, so it can fire while a ColorWipe is still animating.
 *          The guard drops that rollover; without it the second wipe would
 *          overwrite palette_start/palette_target, which the live ColorWipe
 *          still holds references to. Both wipes still render, so the smoke
 *          pass cannot see the difference — assert the snapshots survive the
 *          mid-wipe call and that a rollover arms again once the wipe drains.
 */
inline void test_comets_rollover_skipped_mid_wipe() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  const GenerativePalette::Snapshot start = WB::palette_start(effect);
  const GenerativePalette::Snapshot target = WB::palette_target(effect);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  HS_EXPECT_TRUE(palette_snapshots_equal(start, WB::palette_start(effect)));
  HS_EXPECT_TRUE(palette_snapshots_equal(target, WB::palette_target(effect)));

  // step_wipe_rebake burns one armed frame plus WIPE_FRAMES stepped ones.
  for (int f = 0; f <= WB::wipe_frames(); ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), 0);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
}

/** @brief Manual Comets selection restarts the authored path deterministically. */
inline void test_comets_manual_preset_restarts_path() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();
  for (int f = 0; f < 7; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  constexpr size_t PRESET = 4;
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_TRUE(WB::node_orientation(effect) == Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  const Quaternion first_step = WB::node_orientation(effect);
  HS_EXPECT_EQ(WB::trail_length(effect), 1u);

  for (int f = 0; f < 7; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_TRUE(WB::node_orientation(effect) == Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(WB::node_orientation(effect) == first_step);

  // Manual selection pauses the automatic rollover timer as well as selecting
  // the requested path.
  for (int f = 0; f < WB::cycle_period(effect); ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), PRESET);
}

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
 *          asserting every color() call stays in bounds (every baked_palettes
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
    // baked_palettes[i+1] access for part of the sweep, exercising the bounds
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

  using D = Dynamo<DEFAULT_W, DEFAULT_H>;
  using Ring = Filter::World::Trails<D::TRAIL_CAPACITY>;
  static constexpr int trail_capacity() { return D::TRAIL_CAPACITY; }
  static void set_speed(D &d, float v) { d.params.speed = v; }
  static void set_trail_length(D &d, float v) { d.params.trail_length = v; }
  static float trail_ceiling(const D &d) { return d.params.trail_ceiling; }
  static uint32_t emitted_points(const D &d) { return d.emitted_points; }
  static uint32_t emission_points(const D &d) { return d.emission_points; }
  static size_t trail_points(D &d) { return d.filters.get<Ring>().size(); }
};

/**
 * @brief Pins Dynamo's trail-ring ceiling as the thing that keeps the ring from
 *        evicting.
 * @details trail_length_ceiling() caps the live trail at what the ring can hold
 *          for the current emission rate; past it the ring overruns and evicts
 *          points of arbitrary age (flush()'s compaction leaves it unordered),
 *          punching holes in the tail rather than shortening it — corruption
 *          that still renders, so the smoke pass never sees it. Drive the
 *          worst case both sliders allow and assert the requested trail really
 *          would have overrun while the ceiling keeps steady-state occupancy
 *          inside the ring.
 */
inline void test_dynamo_trail_ceiling_bounds_the_ring() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();

  constexpr float MAX_SPEED = 10.0f;  // "Speed" slider bound
  constexpr float MAX_TRAIL = 100.0f; // "Trail Len" slider bound
  WB::set_speed(effect, MAX_SPEED);
  WB::set_trail_length(effect, MAX_TRAIL);

  const auto capacity = static_cast<uint32_t>(WB::trail_capacity());
  for (int f = 0; f < 64; ++f) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_LT(WB::trail_points(effect), static_cast<size_t>(capacity));
  }

  // Points buffered per frame at this speed: emission_points per whole step,
  // and speed_accumulator carries at most one extra step into a frame.
  const uint32_t per_frame =
      WB::emission_points(effect) * (static_cast<uint32_t>(MAX_SPEED) + 1);
  HS_EXPECT_GT(per_frame, 0u);
  HS_EXPECT_GT(per_frame * static_cast<uint32_t>(MAX_TRAIL), capacity);
  HS_EXPECT_LE(per_frame * static_cast<uint32_t>(WB::trail_ceiling(effect)),
               capacity);
}

/**
 * @brief Pins what Dynamo's emitted_points counter actually measures.
 * @details draw_nodes() bumps it from inside the strand's fragment shader, so
 *          the trail-ring ceiling derived from it is only a valid bound while
 *          the rasterizer runs that shader exactly once per point it hands to
 *          the pipeline. Assert that directly: with a trail long enough that
 *          nothing ages out over the window, the ring's growth across a frame
 *          must equal the frame's emitted_points. A rasterizer that shaded a
 *          fragment it then dropped, or plotted one it never shaded, would
 *          silently rescale the bound.
 */
inline void test_dynamo_emitted_points_counts_ring_seeds() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();
  WB::set_trail_length(effect, 100.0f);

  for (int f = 0; f < 3; ++f) {
    const size_t before = WB::trail_points(effect);
    effect.draw_frame();
    effect.advance_display();
    const size_t seeded = WB::trail_points(effect) - before;
    HS_EXPECT_GT(WB::emitted_points(effect), 0u);
    HS_EXPECT_EQ(seeded, static_cast<size_t>(WB::emitted_points(effect)));
  }
}

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
  static size_t trim_start(size_t len, float alpha) {
    return HF::trail_trim_start(len, alpha);
  }
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

/**
 * @brief Pins HopfFibration's trail trim across every trail length and the
 *        whole visible alpha range.
 * @details render_trails() stages points [first, len) and rasterizes
 *          len - first fragments, so a trim that reached len would bind an
 *          empty polyline; the visibility gate ahead of it guarantees
 *          alpha >= MIN_VISIBLE_ALPHA, which clears the trim's own
 *          MIN_ENCODABLE_ALPHA floor and so bounds first at len - 2.
 *          The trim must also be exact — it drops a point only when that
 *          point's outgoing segment cannot paint — and monotone in alpha, since
 *          a brighter trail can never show less of its tail.
 */
inline void test_hopf_trail_trim_keeps_a_segment() {
  using WB = HopfWhiteBox;
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;
  // Ascending, and all at or above the gate render_trails() applies first.
  const float alphas[] = {MIN_VISIBLE_ALPHA, 0.005f, 0.01f, 0.05f, 0.4f, 1.0f};
  for (size_t len = 2; len <= static_cast<size_t>(HF::TRAIL_LEN); ++len) {
    size_t brighter = len;
    for (float alpha : alphas) {
      const size_t first = WB::trim_start(len, alpha);
      HS_EXPECT_LT(first + 1, len);
      HS_EXPECT_LE(first, brighter);
      // The kept segment paints; the one before it (if any) could not.
      HS_EXPECT_GE(static_cast<float>(first + 1) / (len - 1) * alpha,
                   MIN_ENCODABLE_ALPHA);
      if (first > 0)
        HS_EXPECT_LT(static_cast<float>(first) / (len - 1) * alpha,
                     MIN_ENCODABLE_ALPHA);
      brighter = first;
    }
  }
}

// ---------------------------------------------------------------------------
// Closed-form geometry pins: a ray-march cull sphere, a pixel-pitch star
// radius, and a trail-vertex scratch budget.
//
// Each is a constant the effect derives once and then trusts. The smoke pass
// renders all three and sees none of them: a cull sphere that understates its
// volume just silently clips surface, a star radius off by a factor renders a
// plausible frame at the wrong size, and a scratch estimate that under-counts
// only shows up as an arena overrun on some other resolution.
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for Raymarch's torus proportions (befriended in
 *        effects/Raymarch.h).
 */
struct RaymarchWhiteBox {
  using RM = Raymarch<DEFAULT_W, DEFAULT_H>;
  static constexpr float MAJOR = RM::MAJOR_K;
  static constexpr float MINOR = RM::MINOR_K;
  static constexpr float TWIST = RM::TWIST_K;
  static constexpr float VIS = RM::VIS_K;
  static constexpr float BOUNDS = RM::UNIT_BOUNDS;
};

/**
 * @brief Pins the constexpr Newton square root behind UNIT_BOUNDS against libm.
 * @details Eight fixed iterations from a unit seed. It runs only at compile
 *          time, so nothing else would notice it under-converging, and its one
 *          consumer is the cull-sphere radius — where a low answer culls real
 *          surface. A non-positive radicand short-circuits to 0.
 */
inline void test_raymarch_constexpr_sqrt_converges() {
  using WB = RaymarchWhiteBox;
  double worst = 0.0;
  for (int i = 1; i <= 80; ++i) {
    const float x = 0.05f * static_cast<float>(i); // 0.05 .. 4.0
    worst = std::max(worst, std::fabs(static_cast<double>(constexpr_sqrt(x)) -
                                      std::sqrt(static_cast<double>(x))));
  }
  HS_EXPECT_LT(worst, 1e-6);

  const float radicand = WB::MAJOR * WB::MAJOR + WB::TWIST * WB::TWIST;
  HS_EXPECT_NEAR(constexpr_sqrt(radicand), std::sqrt(radicand), 1e-7);

  static_assert(constexpr_sqrt(0.0f) == 0.0f);
  static_assert(constexpr_sqrt(-1.0f) == 0.0f);
  HS_EXPECT_EQ(constexpr_sqrt(0.0f), 0.0f);
}

/**
 * @brief Verifies UNIT_BOUNDS really bounds the twisted tube it culls against.
 * @details Raymarch hands scale*UNIT_BOUNDS (+ the AA pad) to Scan::Volume::draw
 *          as the ray-march cull sphere, so any surface outside it is dropped
 *          with no other tell. The twisted torus surface has a closed form —
 *          the tube circle of radius MINOR_K about the centerline
 *          (MAJOR_K, TWIST_K*sin(n*theta)) in the (xz-radius, y) half plane —
 *          so the whole surface is swept over the "Twist" slider's integer
 *          domain, checked against the production SDF's zero set, and required
 *          to sit inside the sphere. The sphere is also required to be tight: a
 *          slack radius is wasted ray steps at every vertex.
 */
inline void test_raymarch_unit_bounds_contains_twisted_tube() {
  using WB = RaymarchWhiteBox;
  constexpr double TWO_PI_DBL = 6.283185307179586;
  const double R = WB::MAJOR, r = WB::MINOR, A = WB::TWIST;
  const double bound = WB::BOUNDS;

  double max_off_surface = 0.0; // |SDF| at the sampled surface points
  double max_radius = 0.0;      // farthest surface point from the torus centre
  double min_on_shell = 1e30;   // least SDF value anywhere on the cull sphere

  // twist_n rounds the "Twist" slider, whose range is 0..8.
  for (int n = 0; n <= 8; ++n) {
    const SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{
        SDF::Torus{static_cast<float>(R), static_cast<float>(r)},
        SDF::Warp::Twist{n, static_cast<float>(A), static_cast<float>(R)}};

    for (int i = 0; i < 256; ++i) {
      const double th = TWO_PI_DBL * i / 256;
      const double ct = std::cos(th), st = std::sin(th);
      const double y_mid = A * std::sin(n * th);
      for (int j = 0; j < 256; ++j) {
        const double a = TWO_PI_DBL * j / 256;
        const double s = R + r * std::cos(a);
        const double y = y_mid + r * std::sin(a);
        const Vector p(static_cast<float>(s * ct), static_cast<float>(y),
                       static_cast<float>(s * st));
        max_off_surface =
            std::max(max_off_surface,
                     std::fabs(static_cast<double>(wv.raw_distance(p))));
        max_radius = std::max(max_radius, std::sqrt(s * s + y * y));
      }
    }

    // Sweep the cull sphere itself: no point of it may be inside the volume.
    for (int i = 0; i < 128; ++i) {
      const double th = TWO_PI_DBL * i / 128;
      for (int j = 0; j <= 64; ++j) {
        const double ph = 0.5 * TWO_PI_DBL * j / 64;
        const double sp = std::sin(ph);
        const Vector q(static_cast<float>(bound * sp * std::cos(th)),
                       static_cast<float>(bound * std::cos(ph)),
                       static_cast<float>(bound * sp * std::sin(th)));
        min_on_shell =
            std::min(min_on_shell, static_cast<double>(wv.raw_distance(q)));
      }
    }
  }

  HS_EXPECT_LT(max_off_surface, 2e-4); // the samples are the SDF's zero set
  HS_EXPECT_LE(max_radius, bound + 1e-6);
  HS_EXPECT_GT(min_on_shell, -1e-4);
  HS_EXPECT_NEAR(max_radius, bound, 1e-3);

  // The visible outer radius the per-vertex auto-size fits to the neighbour gap
  // is the ring rim, so it must sit inside the cull sphere.
  HS_EXPECT_LT(static_cast<double>(WB::VIS), bound);
}

/**
 * @brief White-box accessor for GnomonicStars' pixel-pitch star radius and its
 *        spiral cache (befriended in effects/GnomonicStars.h).
 */
struct GnomonicStarsWhiteBox {
  template <int W, int H> static constexpr float radius_px() {
    return GnomonicStars<W, H>::RADIUS_PX;
  }
  template <int W, int H> static int max_points(const GnomonicStars<W, H> &) {
    return GnomonicStars<W, H>::MAX_POINTS;
  }
  template <int W, int H>
  static void set_points(GnomonicStars<W, H> &fx, float n) {
    fx.params.points = n;
  }
  template <int W, int H>
  static int cached_points(const GnomonicStars<W, H> &fx) {
    return fx.cached_points;
  }
  template <int W, int H>
  static Vector cache_at(const GnomonicStars<W, H> &fx, int i) {
    return fx.spiral_cache[i];
  }
  template <int W, int H>
  static void poison_cache(GnomonicStars<W, H> &fx, int i, const Vector &v) {
    fx.spiral_cache[i] = v;
  }
};

/**
 * @brief Pins GnomonicStars' spiral-cache invalidation: the trig-heavy base
 *        lattice is rebuilt exactly when the clamped point count changes.
 * @details draw_frame() clamps "Points" into [1, MAX_POINTS] and rebuilds
 *          spiral_cache only on a change, so a cached_points left in step with
 *          the raw slider would draw the previous count's lattice, and one left
 *          behind would rebuild the trig every frame. The poisoned slot
 *          separates "not rebuilt" from "rebuilt to the same values", which
 *          re-reading the lattice alone cannot.
 */
inline void test_gnomonicstars_spiral_cache_invalidation() {
  using WB = GnomonicStarsWhiteBox;
  reset_effect_globals();
  GnomonicStars<SMALL_W, SMALL_H> fx;
  fx.init();

  auto step = [&](float points) {
    WB::set_points(fx, points);
    fx.draw_frame();
    fx.advance_display();
  };
  auto expect_lattice = [&](int n) {
    for (int i = 0; i < n; ++i) {
      const Vector want = fib_spiral(n, 0.5f, i);
      const Vector got = WB::cache_at(fx, i);
      HS_EXPECT_EQ(got.x, want.x);
      HS_EXPECT_EQ(got.y, want.y);
      HS_EXPECT_EQ(got.z, want.z);
    }
  };

  step(64.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 64);
  expect_lattice(64);

  // An unchanged count leaves the cache alone, so the poison survives.
  const Vector poison(0.0f, 0.0f, 1.0f);
  WB::poison_cache(fx, 7, poison);
  step(64.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 64);
  HS_EXPECT_EQ(WB::cache_at(fx, 7).z, poison.z);

  // A changed count rebuilds every slot, the poisoned one included.
  step(40.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 40);
  expect_lattice(40);

  // Sub-1 and over-capacity slider values clamp before they reach the cache.
  step(0.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 1);
  expect_lattice(1);
  step(static_cast<float>(WB::max_points(fx)) + 500.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), WB::max_points(fx));
}

/**
 * @brief Verifies RADIUS_PX is one pixel of azimuth at every build resolution.
 * @details Star sizes are authored as multiples of RADIUS_PX so they hold their
 *          pixel size across resolutions. Scan::Star turns its radius argument
 *          into SDF::Star's circumradius of radius*pi/2 radians, and one column
 *          spans 2*pi/W radians, so k*RADIUS_PX must come out as exactly k
 *          columns. The rendered pass then measures the lit set of a real star
 *          drawn on the equator, where a column is a full 2*pi/W of arc.
 */
inline void test_gnomonicstars_radius_px_spans_one_column() {
  using WB = GnomonicStarsWhiteBox;
  constexpr int W = DEFAULT_W, H = DEFAULT_H;
  const double column = 2.0 * static_cast<double>(PI_F) / W;
  const double small_column = 2.0 * static_cast<double>(PI_F) / SMALL_W;
  const Basis basis = make_basis(Quaternion(), X_AXIS);

  for (int k : {1, 2, 12}) {
    const SDF::Star shape(basis, k * WB::radius_px<W, H>(), 5, 0.0f);
    HS_EXPECT_NEAR(shape.circumradius, k * column, 1e-6);
    const SDF::Star small(basis, k * WB::radius_px<SMALL_W, SMALL_H>(), 5,
                          0.0f);
    HS_EXPECT_NEAR(small.circumradius, k * small_column, 1e-6);
  }

  constexpr int SPAN_PX = 12;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Scan::Star::draw<W, H, false>(
        pipe, c, basis, SPAN_PX * WB::radius_px<W, H>(), /*sides=*/5,
        [](const Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  size_t lit = 0;
  double max_arc = 0.0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      if (p.r == 0 && p.g == 0 && p.b == 0)
        continue;
      ++lit;
      const Vector v = pixel_to_vector<W, H>(x, y);
      max_arc = std::max(max_arc, static_cast<double>(std::acos(
                                      hs::clamp(dot(v, X_AXIS), -1.0f, 1.0f))));
    }

  HS_EXPECT_GT(lit, (size_t)0);
  // The tips reach SPAN_PX columns and nothing lands beyond them, allowing the
  // AA fringe plus a column of pixel quantization.
  HS_EXPECT_LE(max_arc, (SPAN_PX + 2) * column);
  HS_EXPECT_GE(max_arc, (SPAN_PX - 2) * column);
}

/**
 * @brief White-box accessor for Fishbowl's trail node (befriended in
 *        effects/Fishbowl.h).
 */
struct FishbowlWhiteBox {
  using EffectType = Fishbowl<SMALL_W, SMALL_H>;
  static size_t tween_vertices(const EffectType &fx) {
    size_t n = 0;
    deep_tween(fx.node->trail, [&](const Quaternion &, float) { ++n; });
    return n;
  }
  static Color4 sample_fire(const EffectType &fx, float t) {
    return fx.fire_palette.get(fx.duty_mod.modify(t));
  }
  static float palette_fill_scale(size_t trail_length) {
    return EffectType::palette_fill_scale(trail_length);
  }
  static float palette_phase_speed(float cycle_speed, float scale_factor,
                                   size_t trail_length) {
    return EffectType::palette_phase_speed(cycle_speed, scale_factor,
                                           trail_length);
  }
  static Color4 shade_trail(const EffectType &fx, float palette_t,
                            float age_t) {
    return fx.shade_trail(palette_t, age_t);
  }
  static bool needs_adaptive_midpoint(const Vector &a, const Vector &mid,
                                      const Vector &b) {
    return EffectType::needs_adaptive_midpoint(a, mid, b);
  }
};

/** @brief Verifies Fishbowl's sole preset and fire-palette duty window. */
inline void test_fishbowl_preset_and_fire_duty_cycle() {
  reset_effect_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  auto value = [&](const char *name) {
    const auto *def = fx.getParameters().find(name);
    HS_EXPECT(def != nullptr, "Fishbowl parameter is missing");
    return def ? def->get() : -1.0f;
  };

  HS_EXPECT_EQ(value("Alpha"), 1.0f);
  HS_EXPECT_EQ(value("Cycle Dur"), 80.0f);
  HS_EXPECT_EQ(value("Speed"), 0.235f);
  HS_EXPECT_EQ(value("Jitter Amp"), 2.88f);
  HS_EXPECT_EQ(value("Noise Scale"), 0.26974f);
  HS_EXPECT_EQ(value("Scale Factor"), 84.832001f);
  HS_EXPECT_EQ(value("Cycle Speed"), 0.672f);
  HS_EXPECT_EQ(value("Duty Cycle"), 0.5f);

  const Pixel black(0, 0, 0);
  const Color4 red = WB::sample_fire(fx, 0.20f);
  const Color4 yellow = WB::sample_fire(fx, 0.39f);
  HS_EXPECT_GT(red.color.r, red.color.g);
  HS_EXPECT_GT(red.color.r, red.color.b);
  HS_EXPECT_GT(yellow.color.g, red.color.g);
  HS_EXPECT_TRUE(WB::sample_fire(fx, 0.50f).color == black);
  HS_EXPECT_TRUE(WB::sample_fire(fx, 0.75f).color == black);
  const float dark_palette_t = 0.75f / value("Scale Factor");
  HS_EXPECT_EQ(WB::shade_trail(fx, dark_palette_t, 0.50f).alpha, 0.0f);

  HS_EXPECT_TRUE(fx.updateParameter("Duty Cycle", 0.25f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(WB::sample_fire(fx, 0.30f).color == black);

  HS_EXPECT_EQ(WB::palette_fill_scale(1), 1.0f / WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH / 2),
               static_cast<float>(WB::EffectType::TRAIL_LENGTH / 2) /
                   WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH), 1.0f);
  const float sample_index = 12.0f;
  const float short_coord = (sample_index / 32.0f) * WB::palette_fill_scale(32);
  const float long_coord = (sample_index / 96.0f) * WB::palette_fill_scale(96);
  HS_EXPECT_NEAR(short_coord, long_coord, 1e-7f);

  const float cycle_speed = value("Cycle Speed");
  const float scale_factor = value("Scale Factor");
  const float coordinate_shift = scale_factor / WB::EffectType::TRAIL_LENGTH;
  const float filling_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH - 1);
  const float saturated_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_NEAR(filling_speed, cycle_speed - coordinate_shift, 1e-7f);
  HS_EXPECT_NEAR(saturated_speed - coordinate_shift, filling_speed, 1e-7f);

  const Vector a = X_AXIS;
  const Vector b = Y_AXIS;
  HS_EXPECT_FALSE(WB::needs_adaptive_midpoint(a, slerp(a, b, 0.5f), b));
  HS_EXPECT_TRUE(WB::needs_adaptive_midpoint(a, Z_AXIS, b));
  HS_EXPECT_EQ(WB::EffectType::MAX_FRAGMENTS,
               2 * WB::EffectType::TRAIL_LENGTH *
                   WB::EffectType::ORIENTATION_SUBSTEPS);
}

/**
 * @brief Verifies the scratch-A split's predicted worst case really covers a
 *        saturated Fishbowl frame.
 * @details SCRATCH_A_BYTES is carved from the global arena against a closed-form
 *          worst case — the MAX_FRAGMENTS vertex buffer, the Multiline fragment
 *          buffer it binds, and rasterize's sub-step cache, all live at once.
 *          The static_assert only checks that estimate against the split, never
 *          against a real frame, so an estimate that under-counts would go
 *          unnoticed until the split shrank. Runs past TRAIL_LENGTH frames so
 *          the trail is full when the peak is read.
 */
inline void test_fishbowl_scratch_estimate_covers_peak() {
  reset_effect_globals();
  using WB = FishbowlWhiteBox;
  using EffectType = WB::EffectType;

  EffectType fx;
  fx.init();
  scratch_arena_a.reset_high_water_mark();

  size_t worst_vertices = 0;
  for (int i = 0; i < EffectType::TRAIL_LENGTH + 4; ++i) {
    fx.draw_frame();
    fx.advance_display();
    worst_vertices = std::max(worst_vertices, WB::tween_vertices(fx));
  }

  constexpr size_t PREDICTED =
      EffectType::MAX_FRAGMENTS * sizeof(typename EffectType::TrailVertex) +
      (EffectType::MAX_FRAGMENTS + 2) * sizeof(Fragment) +
      Plot::rasterize_scratch_a_bytes<SMALL_W>();
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), PREDICTED);
  // The trail fills, so the peak above is a saturated frame and not a warm-up.
  HS_EXPECT_GT(worst_vertices, (size_t)EffectType::TRAIL_LENGTH);
  HS_EXPECT_LE(worst_vertices, (size_t)EffectType::MAX_FRAGMENTS);
}

/**
 * @brief White-box accessor for RingSpin's ring pool (befriended in
 *        effects/RingSpin.h).
 */
struct RingSpinWhiteBox {
  using RS = RingSpin<DEFAULT_W, DEFAULT_H>;

  /** @brief Sub-frame count of ring @p i's current orientation step. */
  static int substeps(const RS &fx, int i) {
    return fx.rings[i].orientation.length();
  }
  /** @brief World-space axis of ring @p i's great circle at sub-frame @p s. */
  static Vector axis(const RS &fx, int i, int s) {
    return fx.rings[i].orientation.orient(fx.rings[i].normal, s).normalized();
  }
  /** @brief Half-width in radians of the trail-head stroke. */
  static float head_half_width(const RS &fx) {
    return 2.0f * (2.0f * PI_F / DEFAULT_W) * fx.params.thickness;
  }
};

/**
 * @brief Verifies every lit RingSpin pixel sits on one of its great circles.
 * @details RingSpin draws SDF::Ring at radius 1 (a great circle) about each
 *          ring's rotated plane normal, so on the first frame — where the trail
 *          holds a single orientation step — every lit pixel must lie within the
 *          stroke half-width of some ring's equator. This catches a basis,
 *          radius or trail-record regression that still renders a plausible
 *          frame. Alpha 0 must then blank the frame, pinning the trail-colour
 *          early-out.
 */
inline void test_ringspin_trail_hugs_its_great_circles() {
  reset_effect_globals();
  pin_frame_clock(0);
  using WB = RingSpinWhiteBox;
  WB::RS fx;
  fx.init();
  pin_frame_clock(1);
  fx.draw_frame();
  fx.advance_display();

  // Polar deviation from the equator that the stroke half-width admits. The
  // pipeline carries no Screen::AntiAlias, so there is no fringe to allow for;
  // the measured worst case is 0.011 against this 0.035 bound.
  const float band = std::sin(WB::head_half_width(fx));
  int lit = 0, off_circle = 0;
  float worst = 0.0f;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W; ++x) {
      const Pixel p = fx.get_pixel(x, y);
      if (p.r == 0 && p.g == 0 && p.b == 0)
        continue;
      ++lit;
      const Vector v = pixel_to_vector<DEFAULT_W, DEFAULT_H>(x, y);
      float nearest = 1.0f;
      for (int i = 0; i < WB::RS::NUM_RINGS; ++i)
        for (int s = 0; s < WB::substeps(fx, i); ++s)
          nearest = std::min(nearest, std::abs(dot(v, WB::axis(fx, i, s))));
      worst = std::max(worst, nearest);
      if (nearest > band)
        ++off_circle;
    }

  std::printf("  [info] RingSpin: lit=%d off-circle=%d worst |dot|=%.6f band "
              "%.6f\n",
              lit, off_circle, worst, band);
  HS_EXPECT_GT(lit, 0);
  HS_EXPECT_EQ(off_circle, 0);

  HS_EXPECT_TRUE(fx.updateParameter("Alpha", 0.0f) == ParamSetResult::APPLIED);
  fx.draw_frame();
  fx.advance_display();
  int still_lit = 0;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W; ++x) {
      const Pixel p = fx.get_pixel(x, y);
      if (p.r != 0 || p.g != 0 || p.b != 0)
        ++still_lit;
    }
  HS_EXPECT_EQ(still_lit, 0);
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
  static float live_spacing(const PF &pf) { return pf.spacing(); }
  static float move_dist(const PF &pf) { return pf.move_dist(); }
  static int spawns(const PF &pf) { return pf.spawns; }
  static int dropped_spawns(const PF &pf) { return pf.dropped_spawns; }
  static int max_rings() { return PF::MAX_RINGS; }
  static int rings_on_path() { return PF::RINGS_ON_PATH; }
  static int active_rings(const PF &pf) {
    int active = 0;
    for (int i = 0; i < PF::MAX_RINGS; ++i)
      if (pf.rings[i].active)
        ++active;
    return active;
  }
};

/**
 * @brief Verifies the spawn-gap accumulator drains every frame, the hue cursor
 *        stays wrapped, and the ring pool absorbs the worst slider corner.
 * @details check_spawn() integrates Speed*RHO_PER_SPEED into gap_accumulator and
 *          drains it by the live spacing per spawn, so after any frame the
 *          residue must satisfy 0 <= gap < spacing() — a runaway (missing drain)
 *          or a negative residue both fail here. Run at Speed_max x Density_max:
 *          the only corner where a frame's travel exceeds the live spacing, so
 *          the while-loop's multi-spawn branch runs, and the corner RINGS_ON_PATH
 *          is derived for, where the pool bound has no margin left — a dropped
 *          spawn there is a real defect, not a rounding allowance. next_hue is
 *          advanced wrap(.,1) per spawn and must stay [0, 1).
 */
inline void test_petalflow_spawn_gap_bounded() {
  using WB = PetalFlowWhiteBox;
  reset_effect_globals();
  WB::PF pf;
  pf.init();
  HS_EXPECT_TRUE(pf.updateParameter("Speed", 20.0f) == ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(pf.updateParameter("Density", 2.5f) ==
                 ParamSetResult::APPLIED);

  const float spacing = WB::live_spacing(pf);
  // Below this the while-loop emits at most one ring per frame, whatever the
  // frame count.
  HS_EXPECT_GT(WB::move_dist(pf), spacing);

  // Long enough for the path prefilled at the default density to be fully
  // replaced at the tightest one, so the pool bound is read at steady state
  // rather than mid-fill.
  const int frames = smoke_frames() < 128 ? 128 : smoke_frames();
  int worst_active = 0;
  int worst_burst = 0;
  int previous_spawns = WB::spawns(pf);
  for (int f = 0; f < frames; ++f) {
    pf.draw_frame();
    pf.advance_display();
    const float gap = WB::gap(pf);
    HS_EXPECT_GE(gap, 0.0f);
    HS_EXPECT_LT(gap, spacing); // accumulator never runs away past one spacing
    const float hue = WB::next_hue(pf);
    HS_EXPECT_GE(hue, 0.0f);
    HS_EXPECT_LT(hue, 1.0f); // hue cursor stays wrapped
    const int spawns = WB::spawns(pf);
    worst_burst = std::max(worst_burst, spawns - previous_spawns);
    previous_spawns = spawns;
    worst_active = std::max(worst_active, WB::active_rings(pf));
    HS_EXPECT_EQ(WB::dropped_spawns(pf), 0); // every spawn found a free slot
  }
  HS_EXPECT_GE(worst_burst, 2); // the multi-spawn branch actually ran
  HS_EXPECT_LE(worst_active, WB::rings_on_path());
  HS_EXPECT_LE(WB::rings_on_path(), WB::max_rings());
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
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  GenerativePalette palette(PaletteRecipes::profile(PaletteDomain::MIRROR,
                                                    PaletteHarmony::ANALOGOUS,
                                                    AxisCurve::CONSTANT, 0.0f));
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
 * @details Every bound below is the measured worst case over the sweep with
 * headroom: peak deltaE is 0.0015 default / 0.0053 cyclic against the 0.01
 * bound, and the paired peaks in encoded space are 8 / 19 sRGB8 codes against
 * 10 / 21. The sRGB8 pair is the looser gate because the encode is non-linear —
 * the same table-interpolation error spans more 8-bit codes where the transfer
 * curve is steep than deltaE weights it — so it is bounded rather than pinned to
 * the perceptual figure.
 */
inline void test_displacement_field_hue_table_fidelity() {
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();

  constexpr float INV16 = 1.0f / 65535.0f;
  float default_delta_e = 0.0f;
  float cyclic_delta_e = 0.0f;
  int default_srgb_delta = 0;
  int cyclic_srgb_delta = 0;
  for (uint32_t seed = 0; seed < 12; ++seed) {
    GenerativePalette palette(PaletteRecipes::profile(
        PaletteDomain::MIRROR, PaletteHarmony::ANALOGOUS, AxisCurve::CONSTANT,
        PaletteRecipes::hue_turns(seed * 17u)));
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

  size_t lit = 0, sampled_pixels = 0;
  auto fold_region = [&](bool clip, const Quad &q) -> uint64_t {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    DisplacementField<DEFAULT_W, DEFAULT_H> fx;
    fx.init();
    if (clip)
      fx.set_clip(q.y0, q.y1, q.x0, q.x1);
    uint64_t fold = hs_test::FNV1A64_BASIS;
    for (int f = 0; f < frames; ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      fx.draw_frame();
      fx.advance_display();
      for (int y = q.y0; y < q.y1; ++y)
        for (int x = q.x0; x < q.x1; ++x) {
          const Pixel p = fx.get_pixel(x, y);
          if (!clip) {
            if (p.r | p.g | p.b)
              ++lit;
            ++sampled_pixels;
          }
          for (uint16_t c : {p.r, p.g, p.b})
            fold = hs_test::fnv1a64_channel(fold, c);
        }
    }
    hs::clear_mock_time();
    return fold;
  };

  for (const Quad &q : quads)
    HS_EXPECT_EQ(fold_region(false, q), fold_region(true, q));

  // Two all-black folds agree, so the comparison above only means something
  // once the quadrants have produced output.
  if (lit == 0)
    std::printf("  CLIP-TILE DARK DisplacementField no lit pixel over %zu "
                "sampled pixels\n",
                sampled_pixels);
  HS_EXPECT(lit > 0, "clip tiling must compare a lit render");
}

/** @brief Verifies MindSplatter's Platonic emitter/dual-attractor selector. */
inline void test_mindsplatter_base_mesh_selector() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();

  MS effect;
  effect.init();

  const auto *base_mesh = effect.getParameters().find("Base Mesh");
  HS_EXPECT_TRUE(base_mesh != nullptr);
  if (!base_mesh)
    return;

  HS_EXPECT_EQ(base_mesh->option_count, 5);
  HS_EXPECT_TRUE(base_mesh->animated);
  HS_EXPECT_EQ(std::string_view(base_mesh->options[0]), "Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[4]), "Icosahedron");
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 0), MS::BaseMesh::CUBE);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 1), MS::BaseMesh::DODECAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 4), MS::BaseMesh::TETRAHEDRON);

  const auto select = [&](MS::BaseMesh mesh, size_t emitters,
                          size_t attractors) {
    HS_EXPECT_TRUE(
        effect.updateParameter("Base Mesh", static_cast<float>(mesh)) ==
        ParamSetResult::APPLIED);
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_EQ(WB::active_base_mesh(effect), mesh);
    HS_EXPECT_EQ(WB::active_emitters(effect), emitters);
    HS_EXPECT_EQ(WB::active_attractors(effect), attractors);
  };

  select(MS::BaseMesh::TETRAHEDRON, 4, 4);
  select(MS::BaseMesh::CUBE, 8, 6);
  select(MS::BaseMesh::OCTAHEDRON, 6, 8);
  select(MS::BaseMesh::DODECAHEDRON, 20, 12);
  select(MS::BaseMesh::ICOSAHEDRON, 12, 20);
}

/**
 * @brief Replays one frozen saturated state with exact particle and frame output.
 */
inline void test_mindsplatter_replay_snapshot_exact() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int WARMUP_FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  using Snapshot = WB::ReplaySnapshot<W, H>;
  static_assert(WB::trail_length() == 8);

  Snapshot source;
  {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    MS effect;
    effect.init();
    HS_EXPECT_EQ(WB::particle_capacity(effect), static_cast<size_t>(1672));
    for (int frame = 0; frame < WARMUP_FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
    }
    WB::fill_particle_capacity(effect);
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
 * @brief MindSplatter's single-pass direct-AA path preserves cached-path
 * coverage in every quadrant of a frozen saturated particle pool.
 * @details Single-pass stepping omits the cached endpoint normalization, so
 * interior sample phases, one fringe pixel, and accumulated channels can
 * differ. The production-resolution replay executable gates their error.
 */
inline void test_mindsplatter_saturated_quadrant_sink_parity() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int WARMUP_FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  using Snapshot = WB::ReplaySnapshot<W, H>;

  Snapshot source;
  {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    MS source_effect;
    source_effect.init();
    for (int frame = 0; frame < WARMUP_FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      source_effect.draw_frame();
      source_effect.advance_display();
    }
    WB::fill_particle_capacity(source_effect);
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
    size_t coverage_differences = 0;
    size_t changed_pixels = 0;
    int max_channel_error = 0;
    uint64_t total_channel_error = 0;
    for (int y = quadrant.y0; y < quadrant.y1; ++y) {
      for (int x = quadrant.x0; x < quadrant.x1; ++x) {
        const size_t i = static_cast<size_t>(y) * W + x;
        const bool lit =
            (reference[i].r | reference[i].g | reference[i].b) != 0;
        const bool direct_lit = (direct[i].r | direct[i].g | direct[i].b) != 0;
        coverage_differences += direct_lit != lit ? 1 : 0;
        bool changed = false;
        for (int delta : {std::abs(static_cast<int>(direct[i].r) -
                                   static_cast<int>(reference[i].r)),
                          std::abs(static_cast<int>(direct[i].g) -
                                   static_cast<int>(reference[i].g)),
                          std::abs(static_cast<int>(direct[i].b) -
                                   static_cast<int>(reference[i].b))}) {
          max_channel_error = std::max(max_channel_error, delta);
          total_channel_error += static_cast<uint64_t>(delta);
          changed = changed || delta != 0;
        }
        changed_pixels += changed ? 1 : 0;
        if (lit)
          ++lit_pixels;
      }
    }
    std::printf("sink parity quadrant x[%d,%d) y[%d,%d) lit=%zu changed=%zu "
                "coverage=%zu max_channel=%d total=%llu\n",
                quadrant.x0, quadrant.x1, quadrant.y0, quadrant.y1, lit_pixels,
                changed_pixels, coverage_differences, max_channel_error,
                static_cast<unsigned long long>(total_channel_error));
    HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
    HS_EXPECT_LE(coverage_differences, static_cast<size_t>(1));
  }
  hs::clear_mock_time();
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

/**
 * @brief Bounds rendered output drift from the matrix orientation path.
 * @details The matrix and the quaternion orient() spell the same rotation
 * differently and the shipping builds are -ffast-math, which reassociates each
 * spelling on its own, so the peak channel delta is bounded rather than pinned.
 * Coverage stays exact: both light the same pixels. The pixel and total-error
 * bounds stay tight, so a wrong matrix — which moves the whole render — is still
 * caught; only the peak carries the measured worst case (2 of 65535) with
 * headroom.
 */
inline void test_mindsplatter_rotation_matrix_framebuffer_error() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
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
  HS_EXPECT_LE(different_pixels, static_cast<size_t>(96));
  HS_EXPECT_LE(max_channel_error, 8);
  HS_EXPECT_LE(total_channel_error, static_cast<uint64_t>(128));
}

/** @brief Particle hue seeds advance in deterministic emission order. */
inline void test_mindsplatter_particle_gradients_follow_emission_order() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  MS effect;
  effect.init();
  const auto first = WB::trail_palette(effect, 0);
  const auto second = WB::trail_palette(effect, 32768);
  HS_EXPECT_TRUE(first != second);
  HS_EXPECT_NE(first.front(), first.back());
  HS_EXPECT_NE(second.front(), second.back());
  WB::step_physics(effect);
  HS_EXPECT_EQ(WB::active_particles(effect), WB::num_emitters(effect));
  for (size_t i = 0; i < WB::active_particles(effect); ++i)
    HS_EXPECT_EQ(WB::particle_color_seed(effect, i),
                 static_cast<uint16_t>(i << 8));
}

/**
 * @brief A manual preset outlasts the automatic blend it interrupts.
 * @details The automatic Lerp is frozen rather than finished while the manual
 *          selection holds animations paused, so it resumes on unpause.
 */
inline void test_mindsplatter_manual_preset_survives_unpause() {
  constexpr size_t MANUAL_PRESET = 5;
  constexpr int BLEND_FRAMES = 48;
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  MS effect;
  effect.init();

  WB::advance_preset(effect);
  for (int f = 0; f < 8; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  HS_EXPECT_TRUE(effect.selectPreset(MANUAL_PRESET));
  HS_EXPECT_TRUE(effect.animations_paused());
  effect.setAnimationsPaused(false);
  for (int f = 0; f < BLEND_FRAMES + 8; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  HS_EXPECT_EQ(effect.getPresetIndex(), MANUAL_PRESET);
  const auto &live = WB::live_params(effect);
  const auto &chosen = WB::preset_params(effect, MANUAL_PRESET);
  HS_EXPECT_TRUE(live.base_mesh == chosen.base_mesh);
  HS_EXPECT_NEAR(live.friction, chosen.friction, 1e-6f);
  HS_EXPECT_NEAR(live.well_strength, chosen.well_strength, 1e-6f);
  HS_EXPECT_NEAR(live.initial_speed, chosen.initial_speed, 1e-6f);
  HS_EXPECT_NEAR(live.angular_speed, chosen.angular_speed, 1e-6f);
  HS_EXPECT_NEAR(live.warp_scale, chosen.warp_scale, 1e-6f);
}

/** @brief Fusing the vertex pass into trail materialization is pixel exact. */
inline void test_mindsplatter_fused_vertex_framebuffer_parity() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
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

/**
 * @brief The multiply-only hole kernel matches the generic kernel.
 * @details Coverage is exact: both light the same pixels. Channel values carry a
 * tolerance because the shipping builds are -ffast-math, which reassociates the
 * branchless product and the max/acos spelling of the same falloff
 * independently, so the two agree to float rounding rather than bit-exactly.
 * The bounds are the measured worst case with headroom — 7 of 307200 samples
 * differ, by one 16-bit step. A wrong kernel changes the falloff shape across
 * the whole event horizon, which the differing-pixel count still catches.
 */
inline void test_mindsplatter_hole_kernel_framebuffer_parity() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  auto render = [&](bool reference) {
    reset_effect_globals();
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
  size_t coverage_differences = 0;
  int max_channel_error = 0;
  uint64_t total_channel_error = 0;
  for (size_t i = 0; i < reference.size(); ++i) {
    const Pixel a = reference[i];
    const Pixel b = multiply_only[i];
    if (a.r | a.g | a.b)
      ++lit_pixels;
    if (a != b)
      ++different_pixels;
    if (((a.r | a.g | a.b) == 0) != ((b.r | b.g | b.b) == 0))
      ++coverage_differences;
    for (int delta :
         {std::abs(static_cast<int>(a.r) - static_cast<int>(b.r)),
          std::abs(static_cast<int>(a.g) - static_cast<int>(b.g)),
          std::abs(static_cast<int>(a.b) - static_cast<int>(b.b))}) {
      max_channel_error = std::max(max_channel_error, delta);
      total_channel_error += static_cast<uint64_t>(delta);
    }
  }
  std::printf("hole kernel framebuffer samples=%zu lit=%zu different=%zu "
              "coverage=%zu max_channel=%d total_channel=%llu\n",
              reference.size(), lit_pixels, different_pixels,
              coverage_differences, max_channel_error,
              static_cast<unsigned long long>(total_channel_error));
  HS_EXPECT_GT(lit_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(coverage_differences, static_cast<size_t>(0));
  HS_EXPECT_LE(different_pixels, static_cast<size_t>(64));
  HS_EXPECT_LE(max_channel_error, 8);
  HS_EXPECT_LE(total_channel_error, static_cast<uint64_t>(64));
}

/** @brief Clip clearing preserves every pixel displayed by the POV driver. */
inline void test_mindsplatter_clip_clear_display_parity() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
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

/**
 * @brief Bounds full-lifetime render drift from signed-axis physics.
 * @details Particle count stays exact at every frame, and coverage stays exact
 * at every checkpoint: the two spellings agree on which particles live and which
 * pixels light. The per-checkpoint budgets grow with the frame index because the
 * drift compounds through the integrator. They also carry the -ffast-math the
 * shipping builds use, which reassociates the two spellings independently and
 * roughly doubles the mid-lifetime drift; the budgets are the measured worst
 * case across both math configurations with headroom.
 */
inline void test_mindsplatter_signed_axis_framebuffer_error() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int FRAMES = 160;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  struct Render {
    std::vector<Pixel> frames;
    std::vector<uint16_t> active;
  };
  auto render = [&](bool reference) {
    reset_effect_globals();
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
  constexpr size_t MAX_DIFFERENT[] = {0, 192, 192};
  constexpr int MAX_CHANNEL[] = {0, 256, 1152};
  constexpr uint64_t MAX_TOTAL[] = {0, 2048, 4608};
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
 * @brief Verifies every per-emitter emission phase stays wrapped to [0, 2pi)
 *        across frames at the max angular rate.
 * @details Each emitter integrates Ang Spd into emit_phases[i] with fmodf(., 2pi)
 *          so a dropped wrap lets the phase grow unbounded (fast_sinf range
 *          reduction then bands). The range assertions only bind once the phase
 *          has passed 2pi at least once, so the sweep runs at the Ang Spd slider
 *          top and counts the laps it observes.
 */
inline void test_mindsplatter_emit_phase_wrapped() {
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  WB::MS ms;
  ms.init();
  HS_EXPECT_TRUE(ms.updateParameter("Ang Spd", 1.0f) ==
                 ParamSetResult::APPLIED);

  const float two_pi = 2.0f * PI_F;
  const int frames = smoke_frames() < 64 ? 64 : smoke_frames();
  std::vector<float> previous(WB::num_emitters(ms), 0.0f);
  int laps = 0;
  for (int f = 0; f < frames; ++f) {
    ms.draw_frame();
    ms.advance_display();
    for (size_t i = 0; i < previous.size(); ++i) {
      const float ph = WB::emit_phase(ms, i);
      HS_EXPECT_GE(ph, 0.0f);
      HS_EXPECT_LT(ph, two_pi);
      if (ph < previous[i])
        ++laps;
      previous[i] = ph;
    }
  }
  HS_EXPECT_GT(laps, 0);
}

/**
 * @brief Verifies glitch_lens maps unit directions to unit directions and
 *        stays smooth through its equatorial fold.
 * @details Pins |lens(v)| == 1 across a spread of directions, its doubled-
 *          latitude topology, and matching one-sided derivatives at the fold.
 */
inline void test_shaderball_glitch_lens_unit_norm() {
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
    HS_EXPECT_NEAR(lenses::glitch_lens(v).length(), 1.0f, 1e-3f);
  }

  const Vector equator_x = lenses::glitch_lens(Vector(1, 0, 0));
  HS_EXPECT_NEAR(equator_x.x, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(equator_x.y, -1.0f, 1e-6f);
  HS_EXPECT_NEAR(equator_x.z, 0.0f, 1e-6f);

  const Vector north = lenses::glitch_lens(Vector(0, 1, 0));
  HS_EXPECT_NEAR(north.y, 1.0f, 1e-6f);
  const Vector south = lenses::glitch_lens(Vector(0, -1, 0));
  HS_EXPECT_NEAR(south.y, 1.0f, 1e-6f);

  constexpr float EPSILON = 1e-4f;
  const float radial = sqrtf(1.0f - EPSILON * EPSILON);
  const Vector center = lenses::glitch_lens(Vector(0.8f, 0.0f, 0.6f));
  const Vector above =
      lenses::glitch_lens(Vector(0.8f * radial, EPSILON, 0.6f * radial));
  const Vector below =
      lenses::glitch_lens(Vector(0.8f * radial, -EPSILON, 0.6f * radial));
  const Vector forward = above - center;
  const Vector backward = center - below;
  HS_EXPECT_NEAR(forward.x, backward.x, 1e-6f);
  HS_EXPECT_NEAR(forward.y, backward.y, 1e-6f);
  HS_EXPECT_NEAR(forward.z, backward.z, 1e-6f);
}

/** @brief White-box accessor for MobiusRings' singular numeric branches. */
struct MobiusRingsWhiteBox {
  using MR = MobiusRings<DEFAULT_W, DEFAULT_H>;
  static float conformal_coord(float z, float phase) {
    return MR::conformal_coord(z, phase);
  }
  static Quaternion counter_rotation(const Vector &mid) {
    return MR::counter_rotation(mid);
  }
};

inline void test_mobius_rings_conformal_and_counter_rotation() {
  using WB = MobiusRingsWhiteBox;

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
 * @brief Verifies ShapeShifter's slider contract and preset-row invariants.
 * @details A preset row is authored artistic data that is retuned freely, so its
 *          magnitudes (count, sides, amplitude, speed) are not pinned: a golden
 *          copy of them reds on every intentional retune and reports nothing. The
 *          invariants a retune must not break are pinned instead — every value
 *          inside the range register_param bound it to, the structural selections
 *          each row draws, the non-preset Alpha slider surviving a selection, the
 *          shape/falloff pairing, and every row landing on a distinct parameter
 *          vector.
 */
inline void test_shapeshifter_preset_defaults() {
  reset_effect_globals();
  ShapeShifter<DEFAULT_W, DEFAULT_H> ss;
  ss.init();

  auto value = [&](const char *name) {
    for (const auto &def : ss.getParameters())
      if (std::strcmp(def.name, name) == 0)
        return def.get();
    HS_EXPECT(false, "ShapeShifter parameter is missing");
    return -1.0f;
  };

  // register_param traps on a default outside its range, but a preset row is
  // assigned straight into params and never passes through it.
  auto expect_in_range = [&](const char *label) {
    HS_CONTEXT(label);
    for (const auto &def : ss.getParameters()) {
      HS_CONTEXT(def.name);
      const float v = def.get();
      HS_EXPECT_TRUE(std::isfinite(v));
      HS_EXPECT_GE(v, def.min);
      HS_EXPECT_LE(v, def.max);
      if (def.option_count > 0) {
        HS_EXPECT_EQ(v, std::floor(v));
        HS_EXPECT_LT(v, static_cast<float>(def.option_count));
      }
    }
  };

  expect_in_range("boot state");
  HS_EXPECT_EQ(value("Alpha"), 1.0f); // boots fully opaque
  HS_EXPECT_EQ(value("Shape"), 3.0f);
  HS_EXPECT_EQ(value("Spacing"), 1.0f);
  HS_EXPECT_EQ(value("Function"),
               static_cast<float>(
                   ShapeShifter<DEFAULT_W, DEFAULT_H>::PhaseFunction::SINE));
  HS_EXPECT_EQ(value("Opposite"), 0.0f);
  HS_EXPECT_EQ(value("Alpha Falloff"), 1.0f);

  const char *expected_export_order[] = {
      "Shape", "Count",    "Sides",         "Function", "Amplitude",
      "Speed", "Opposite", "Alpha Falloff", "Spacing"};
  size_t export_index = 0;
  for (const auto &def : ss.getParameters()) {
    if (!def.preset)
      continue;
    HS_EXPECT(export_index < std::size(expected_export_order),
              "ShapeShifter exports an unexpected parameter");
    if (export_index < std::size(expected_export_order))
      HS_EXPECT_EQ(std::string_view(def.name),
                   std::string_view(expected_export_order[export_index]));
    ++export_index;
  }
  HS_EXPECT_EQ(export_index, std::size(expected_export_order));

  const auto *alpha = ss.getParameters().find("Alpha");
  const auto *shape = ss.getParameters().find("Shape");
  const auto *falloff = ss.getParameters().find("Alpha Falloff");
  const auto *spacing = ss.getParameters().find("Spacing");
  const auto *count = ss.getParameters().find("Count");
  const auto *speed = ss.getParameters().find("Speed");
  HS_EXPECT(alpha != nullptr, "ShapeShifter Alpha parameter is missing");
  HS_EXPECT(shape != nullptr, "ShapeShifter Shape parameter is missing");
  HS_EXPECT(falloff != nullptr,
            "ShapeShifter Alpha Falloff parameter is missing");
  HS_EXPECT(spacing != nullptr, "ShapeShifter Spacing parameter is missing");
  HS_EXPECT(count != nullptr, "ShapeShifter Count parameter is missing");
  HS_EXPECT(speed != nullptr, "ShapeShifter Speed parameter is missing");
  if (count)
    HS_EXPECT_EQ(count->max, 288.0f);
  if (speed)
    HS_EXPECT_EQ(speed->max, 0.16f);
  if (alpha)
    HS_EXPECT_FALSE(alpha->preset);
  if (shape) {
    HS_EXPECT_TRUE(shape->is_enum());
    HS_EXPECT_EQ(std::string_view(shape->export_options[3]),
                 std::string_view("ShapeType::PLANAR_STAR"));
    HS_EXPECT_EQ(std::string_view(shape->export_options[4]),
                 std::string_view("ShapeType::SPHERICAL_STAR"));
  }
  if (falloff) {
    HS_EXPECT_TRUE(falloff->is_enum());
    HS_EXPECT_EQ(std::string_view(falloff->export_options[1]),
                 std::string_view("AlphaFalloff::TOWARD_EQUATOR"));
  }
  if (spacing) {
    HS_EXPECT_TRUE(spacing->is_enum());
    HS_EXPECT_EQ(std::string_view(spacing->export_options[1]),
                 std::string_view("RadiusSpacing::SCREEN_BALANCED"));
  }

  HS_EXPECT_TRUE(ss.updateParameter("Alpha", 0.37f) == ParamSetResult::APPLIED);

  // Structural selections: which primitive and falloff each row draws. These
  // pin the index -> row mapping the profile reports are keyed by; the
  // magnitudes each row sets are deliberately left free.
  const float expected_shapes[] = {3.0f, 1.0f, 3.0f, 2.0f, 3.0f,
                                   1.0f, 1.0f, 1.0f, 2.0f};
  const float expected_falloffs[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f,
                                     0.0f, 0.0f, 0.0f, 0.0f};
  std::vector<std::vector<float>> rows;
  for (size_t i = 0; i < std::size(expected_shapes); ++i) {
    HS_CONTEXT("preset", static_cast<int>(i));
    ss.profile_select_preset(i);
    expect_in_range("preset row");
    HS_EXPECT_EQ(value("Alpha"), 0.37f); // a preset never writes a non-preset
    HS_EXPECT_TRUE(ss.animations_paused());
    HS_EXPECT_EQ(value("Shape"), expected_shapes[i]);
    HS_EXPECT_EQ(value("Alpha Falloff"), expected_falloffs[i]);
    HS_EXPECT_EQ(value("Shape") == 3.0f, value("Alpha Falloff") == 1.0f);

    std::vector<float> row;
    for (const auto &def : ss.getParameters())
      if (def.preset)
        row.push_back(def.get());
    rows.push_back(row);
  }

  // Two rows that collapse onto the same parameter vector are one preset the
  // cycle visits twice, which no range or selection check above would show.
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = i + 1; j < rows.size(); ++j) {
      HS_CONTEXT("preset pair", static_cast<int>(i), static_cast<int>(j));
      HS_EXPECT(rows[i] != rows[j],
                "each ShapeShifter preset must be distinct");
    }
}

/**
 * @brief Renders every Shape and Function slider selection.
 * @details Each primitive is exercised at radii on both sides of the antipode
 * fold while the four phase functions advance through the same Plot pipeline.
 * Every selection renders one frame from an identical fresh state with the
 * preset timer paused, so the only input that moves between renders is the
 * slider: two selections folding to the same frame means the selection switch
 * did not dispatch on them.
 */
inline void test_shapeshifter_slider_selections_render() {
  using SS = ShapeShifter<SMALL_W, SMALL_H>;

  auto render = [](const char *slider, int selection) {
    reset_effect_globals();
    SS ss;
    ss.init();
    ss.setAnimationsPaused(true);
    HS_EXPECT_TRUE(ss.updateParameter(slider, static_cast<float>(selection)) ==
                   ParamSetResult::APPLIED);
    ss.draw_frame();
    ss.advance_display();

    uint64_t acc = 0;
    uint64_t fold = hs_test::FNV1A64_BASIS;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x) {
        const Pixel &pixel = ss.get_pixel(x, y);
        acc += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
        for (uint16_t channel : {pixel.r, pixel.g, pixel.b})
          fold = hs_test::fnv1a64_channel(fold, channel);
      }
    HS_EXPECT_GT(acc, 0u);
    return fold;
  };

  auto sweep = [&](const char *slider, int selections) {
    std::vector<uint64_t> folds;
    for (int selection = 0; selection < selections; ++selection)
      folds.push_back(render(slider, selection));
    for (int i = 0; i < selections; ++i)
      for (int j = i + 1; j < selections; ++j) {
        if (folds[i] == folds[j])
          std::printf("  SHAPESHIFTER %s selections %d and %d render the same "
                      "frame (fold %llu)\n",
                      slider, i, j, static_cast<unsigned long long>(folds[i]));
        HS_EXPECT(folds[i] != folds[j],
                  "each slider selection must render a distinct frame");
      }
  };

  sweep("Shape", SS::NUM_SHAPES);
  sweep("Function", SS::NUM_FUNCTIONS);
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
 *          (targets/wasm/engine_bindings.h setClip). See
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
  Voronoi<SMALL_W, SMALL_H> effect;
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
      voronoi_union_nearest_match<SMALL_W, SMALL_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match_dev, 1.0);

  // Dense regime: seed MAX_SITES on a Fibonacci sphere exactly as
  // Voronoi::seed_sites places them, so the adaptive block floors at
  // COHERENCE_BLOCK_MIN.
  constexpr int N = VoronoiWhiteBox::MAX_SITES;
  static Vector fib[N];
  for (int i = 0; i < N; ++i)
    fib[i] = fib_spiral(N, /*eps=*/0.5f, i);
  const std::span<const Vector> dense(fib, N);
  const double fib_match =
      voronoi_union_nearest_match<DEFAULT_W, DEFAULT_H>(dense, deficit);
  HS_EXPECT_GE(fib_match, 0.999);
  HS_EXPECT_LE(deficit, 0.005f);

  const double fib_match_dev =
      voronoi_union_nearest_match<SMALL_W, SMALL_H>(dense, deficit);
  HS_EXPECT_EQ(fib_match_dev, 1.0);
  HS_EXPECT_EQ(deficit, 0.0f);
}

/**
 * @brief Requires a Voronoi segment band to shade every pixel exactly as the
 *        full-canvas render does.
 * @details The coarse-coherence grid decides per block which sites reach a
 *          pixel's candidate union, so a grid whose phase follows the clip
 *          origin shades the same pixel differently depending on which band
 *          renders it — a discontinuity pinned to the segment seam, since
 *          Voronoi is neither full-frame nor persisting. The bands below start
 *          off a block boundary (the adaptive block is 6 px at the default site
 *          count), which an aligned split would hide.
 */
inline void test_voronoi_segment_render_matches_full_frame() {
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;

  auto render = [](int x0, int x1, int y0, int y1) {
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    effect.set_clip(y0, y1, x0, x1);
    effect.draw_frame();
    effect.advance_display();
    std::vector<Pixel> band;
    band.reserve(static_cast<size_t>(x1 - x0) * (y1 - y0));
    for (int y = y0; y < y1; ++y)
      for (int x = x0; x < x1; ++x)
        band.push_back(effect.get_pixel(x, y));
    return band;
  };

  const std::vector<Pixel> full = render(0, W, 0, H);

  struct Band {
    int x0, x1, y0, y1;
  };
  const Band bands[] = {
      {0, 100, 0, H}, {100, W, 0, H}, {0, W, 50, H}, {37, 205, 11, 93}};

  size_t lit = 0, compared_pixels = 0;
  for (const Band &b : bands) {
    const std::vector<Pixel> banded = render(b.x0, b.x1, b.y0, b.y1);
    HS_EXPECT_EQ(banded.size(),
                 static_cast<size_t>(b.x1 - b.x0) * (b.y1 - b.y0));
    size_t different = 0;
    size_t i = 0;
    for (int y = b.y0; y < b.y1; ++y)
      for (int x = b.x0; x < b.x1; ++x, ++i) {
        const Pixel &reference = full[static_cast<size_t>(y) * W + x];
        if (banded[i] != reference)
          ++different;
        if (reference.r | reference.g | reference.b)
          ++lit;
      }
    compared_pixels += banded.size();
    if (different)
      std::printf("  VORONOI SEAM band x[%d,%d) y[%d,%d): %zu of %zu pixels "
                  "differ from the full-canvas render\n",
                  b.x0, b.x1, b.y0, b.y1, different, banded.size());
    HS_EXPECT_EQ(different, static_cast<size_t>(0));
  }
  // Two all-black renders agree pixel for pixel, so the comparison above only
  // means something once the bands have produced output.
  if (lit == 0)
    std::printf("  VORONOI DARK no lit pixel over %zu compared pixels\n",
                compared_pixels);
  HS_EXPECT(lit > 0, "segment parity must compare a lit render");
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
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed =
          Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), target);
    });

    MeshState mesh;
    CompiledHankin hankin;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
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
  configure_arenas_default();
}

/**
 * @brief White-box accessor for IslamicStars' private build-chain state
 *        (befriended in effects/IslamicStars.h).
 * @details The op-by-op recipe build only runs when the round-robin reaches a
 * non-null recipe entry (index 1), ~340 frames into a default-speed run — past
 * every generic smoke window. The probe pre-sets Trans Speed before init so
 * the whole crossing fits in ~100 frames, and reads build_active/solid_idx
 * to pin that the build ran and completed. <96,20> keeps the raster cheap;
 * the build bookkeeping is resolution-independent.
 */
struct IslamicBuildProbe {
  using IS = IslamicStars<SMALL_W, SMALL_H>;
  template <int W, int H>
  static void set_trans_speed(IslamicStars<W, H> &e, float v) {
    e.params.trans_speed = v;
  }
  template <int W, int H>
  static bool build_active(const IslamicStars<W, H> &e) {
    return e.build_active;
  }
  template <int W, int H> static int solid_idx(const IslamicStars<W, H> &e) {
    return e.solid_idx;
  }
  template <int W, int H> static int dual_bridges(const IslamicStars<W, H> &e) {
    return e.dual_bridges_built;
  }
  template <int W, int H>
  static size_t persistent_budget(const IslamicStars<W, H> &e) {
    return e.device_persistent_budget;
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
  /**
   * @brief Spawns @p entry out of band, whatever the registry cycle is on.
   * @details Clears the timeline first: the pending scheduled spawn and the
   * outgoing shape's sprite would otherwise run against the injected shape.
   */
  template <int W, int H>
  static void spawn_entry(IslamicStars<W, H> &e, const Solids::Entry &entry) {
    e.timeline.clear();
    e.spawn_entry(entry);
  }
};

/** Whole-solid generator of the needle-ending recipe. */
inline PolyMesh generate_needle_recipe_solid(Arena &a, Arena &b) {
  return Solids::build_recipe(
      Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_RECIPE, a, b);
}

/** The needle-ending recipe as a spawnable entry. It is not in
 * islamic_registry, so the arena gate below builds it from the recipe constant
 * rather than finding it in the roster. */
inline constexpr Solids::Entry NEEDLE_ENTRY = {
    "truncatedIcosahedron_ambo_relax100_hk54_needle",
    generate_needle_recipe_solid, Solids::Category::Complex,
    &Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_RECIPE};

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
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &p = effect.get_pixel(x, y);
      acc += static_cast<uint64_t>(p.r) + p.g + p.b;
    }
  HS_EXPECT_GT(acc, (uint64_t)0);
}

/**
 * @brief Drives IslamicStars through every registry entry and then through the
 *        needle recipe, pinning the persistent arena against the effect's own
 *        budget.
 * @details The per-chain gate measures a build in isolation, which misses what
 *          the effect actually holds: a build follows whatever shape preceded
 *          it, and the roster's heaviest entry is 1082 faces. Cycling the whole
 *          roster is the only thing that exercises a build against that
 *          predecessor. Rendered small — the peak is mesh-driven, not
 *          canvas-driven — because the gate is about memory, not pixels. An
 *          Arena overrun traps, so an OOM fails this test by killing the run.
 *          The needle is the heaviest smooth-bridge shape and sets the
 *          scratch_a-heavy split, so it is measured separately.
 */
inline void test_islamicstars_roster_cycle_fits_budget() {
  reset_effect_globals();
  // Per-shape arena split (IslamicStars::spawn_entry): smooth kis/needle
  // bridge shapes use 129.5 KB / 74 KB scratch; ordinary recipe builds use
  // 116 KB / 72 KB; non-recipe generation retains 116 KB / 74 KB. On host the
  // scratch arenas are the exact device sizes and hard-capped, so an over-budget
  // leg traps -- completing the whole roster proves every shape's scratch fit
  // its own split. Persistent is host-inflated (soft), so it is checked per
  // frame against the effect's live per-shape device budget; host pointer
  // inflation makes the resident offset an upper bound on the device figure.
  // A resplit rebases each scratch high-water, so every peak below is the
  // measured shape's own.
  {
    // Production resolution, not device: the swept+compiled+draw scratch peak a
    // build leg reaches is mesh-driven and resolution-independent, but the
    // effect ships at 288x144, so the gate holds the real geometry. Trans Speed
    // 8 compresses each stage so the whole roster cycles in ~1200 frames
    // without dropping any leg (never lower the resolution: that would shrink
    // the peak).
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 8.0f);
    effect.init();

    auto solids = Solids::Collections::get_islamic_solids();
    const int entries = static_cast<int>(solids.size());
    constexpr int MAX_FRAMES = 20000;
    size_t a_peak = 0, b_peak = 0, persist_peak = 0;
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
      const size_t p_budget = IslamicBuildProbe::persistent_budget(effect);
      a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
      b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
      persist_peak = std::max(persist_peak, p);
      if (p > worst_p) { // report the tightest persistent fit, not the raw max
        worst_p = p;
        worst_p_budget = p_budget;
        worst_p_idx = cur;
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
    HS_EXPECT_GT(shapes, entries - 1);
    HS_EXPECT_GT(builds, 0);
  }

  // The needle build. Trans Speed 2, not the roster's 8: a compressed stage can
  // drop the closing bridge leg before it runs, which is the peak.
  reset_effect_globals();
  size_t na_peak = 0, nb_peak = 0, np_peak = 0;
  int needle_frames = 0;
  bool needle_built = false;
  {
    constexpr int NEEDLE_MAX_FRAMES = 4000;
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 2.0f);
    effect.init();
    IslamicBuildProbe::spawn_entry(effect, NEEDLE_ENTRY);
    bool was_building = false;
    while (needle_frames < NEEDLE_MAX_FRAMES) {
      effect.draw_frame();
      effect.advance_display();
      ++needle_frames;
      const bool building = IslamicBuildProbe::build_active(effect);
      const size_t p = persistent_arena.get_offset();
      na_peak = std::max(na_peak, scratch_arena_a.get_high_water_mark());
      nb_peak = std::max(nb_peak, scratch_arena_b.get_high_water_mark());
      np_peak = std::max(np_peak, p);
      HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
      if (building)
        was_building = true;
      else if (was_building) {
        needle_built = true;
        break;
      }
    }
  }
  std::printf("  [needle] smooth-path peaks over %d frames: scratch_a=%zu/%zu "
              "scratch_b=%zu/%zu persistent=%zu/%zu\n",
              needle_frames, na_peak, IslamicBuildProbe::bridge_scratch_a(),
              nb_peak, IslamicBuildProbe::bridge_scratch_b(), np_peak,
              DEVICE_GLOBAL_ARENA_SIZE - IslamicBuildProbe::bridge_scratch_a() -
                  IslamicBuildProbe::bridge_scratch_b());
  HS_EXPECT_TRUE(needle_built);
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
 *          bridges are reached in registry order: a gyro, a gyro-kis dt macro,
 *          then a kis-gyro whose dtd macro bridges twice before its own gyro
 *          bridges again. The bridge's leg 3 rebuilds the medial for
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
      "scratch_b peak=%zu B, persistent peak=%zu B\n",
      IslamicBuildProbe::dual_bridges(effect), frames, a_peak, b_peak,
      persist_peak);
  HS_EXPECT_GE(IslamicBuildProbe::dual_bridges(effect), TARGET_BRIDGES);
}

/**
 * @brief Module entry point for the effects white-box suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Per-effect invariants with an explicit oracle. The roster-wide smoke,
 * determinism and clip-clear parity sweeps live in the separate effects_smoke
 * module (tests/test_effects_smoke.h), which is IEEE-agnostic and so runs on the
 * fast-math axis this module is excluded from.
 */
template <typename EffectT>
inline void check_manual_preset_navigation(size_t expected_count) {
  reset_effect_globals();
  EffectT effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), expected_count);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_FALSE(effect.selectPreset(expected_count));
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_TRUE(effect.previousPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), expected_count - 1);
  HS_EXPECT_TRUE(effect.nextPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_TRUE(effect.animations_paused());
}

inline void test_manual_preset_navigation() {
  check_manual_preset_navigation<MindSplatter<SMALL_W, SMALL_H>>(8);
  check_manual_preset_navigation<DreamBalls<SMALL_W, SMALL_H>>(10);
  check_manual_preset_navigation<Comets<SMALL_W, SMALL_H>>(12);
  check_manual_preset_navigation<MeshFeedback<SMALL_W, SMALL_H>>(12);
  check_manual_preset_navigation<ShapeShifter<SMALL_W, SMALL_H>>(9);
}

inline int run_effects_tests() {
  hs_test::ModuleFixture fixture("effects");

  test_mindsplatter_base_mesh_selector();
  test_meshfeedback_base_mesh_selector();
  test_meshfeedback_preset_export_arity();
  test_meshfeedback_mesh_rebuild_reuses_storage();
  test_mindsplatter_replay_snapshot_exact();
  test_mindsplatter_saturated_quadrant_sink_parity();
  test_mindsplatter_octahedral_hole_alpha_equivalence();
  test_mindsplatter_rotation_matrix_equivalence();
  test_mindsplatter_rotation_matrix_framebuffer_error();
  test_mindsplatter_particle_gradients_follow_emission_order();
  test_mindsplatter_fused_vertex_framebuffer_parity();
  test_mindsplatter_hole_kernel_framebuffer_parity();
  test_mindsplatter_clip_clear_display_parity();
  test_mindsplatter_signed_axis_framebuffer_error();
  test_gs_opaque_quarter_accumulation_is_exact();
  test_gs_reseed_generates_palette();
  test_gs_render_certificates_bound_lattice();
  test_gs_shared_stencil_error_is_bounded();
  test_gs_dissolve_frontier_fades_before_clear();
  test_gs_substep_matches_scalar_reference();
  test_fishbowl_preset_and_fire_duty_cycle();
  test_shapeshifter_preset_defaults();
  test_shapeshifter_slider_selections_render();
  test_manual_preset_navigation();
  test_mindsplatter_manual_preset_survives_unpause();
  test_hankinsolids_manual_pause_holds_morph();
  test_every_effect_renders_while_paused();
  // Arena budgets both tiers run: a mesh or fragment-count change reds these,
  // and they cost under a second between them.
  test_fishbowl_scratch_estimate_covers_peak();
  test_hankinsolids_arena_budget_covers_every_solid();
  test_dreamballs_max_edge_solid_render();

  // FULL tier only (HS_EFFECTS_FULL=1; CI on every push/PR). The QUICK tier
  // (default, local pre-commit) skips the block below entirely, so a green local
  // commit is NOT authoritative for the full-resolution paths or the white-box
  // invariants — CI is (same split as HS_SMOKE_FRAMES's cyclic-window coverage;
  // see the pre-commit hook header).
  if (effects_full_suite()) {
    test_needs_full_frame_gate();
    test_voronoi_axes_use_uniform_sampler();
    test_voronoi_union_candidates_cover_nearest();
    test_voronoi_segment_render_matches_full_frame();
    test_sh_decode_lm_valid_order();
    test_sh_reduced_legendre_matches_closed_form();
    test_sh_cartesian_matches_spherical();
    test_sh_field_write_through_and_endpoints();
    test_sh_polarity_split_and_ao_shaping();
    test_sh_morph_chain_rearms();
    test_gs_q16_roundtrip();
    test_gs_rest_state_is_fixed_point();
    test_gs_substep_signs_and_clamp();
    test_gs_evolution_stays_bounded();
    test_gs_reaction_corner_stays_bounded();
    test_gs_dissolve_clears_and_reseeds();
    test_gs_reaction_edit_starts_dissolve();
    test_bz_legacy_palette();
    test_bz_q16_roundtrip();
    test_bz_min_diffusion_step_survives_quantization();
    test_bz_advance_species_signs_and_clamp();
    test_bz_perturb_state_saturates_and_nudges();
    test_bz_perturb_state_draw_count_pinned();
    test_bz_perturb_scales_with_timestep();
    test_bz_substep_diffuses();
    test_bz_raster_matches_reference();
    test_bz_render_center_matches_reference();
    test_dreamballs_preset_cycle_bookkeeping();
    test_dreamballs_base_mesh_selector();
    test_dreamballs_weave_topology();
    test_dreamballs_respawn_fires_and_honors_pause();
    test_meshfeedback_flush_precedes_mesh_draw();
    test_meshfeedback_preset_rotation_syncs_noise();
    CometsWhiteBox::check_paths_close();
    test_comets_rollover_skipped_mid_wipe();
    test_comets_manual_preset_restarts_path();
    ThrustersWhiteBox::check_warp_endpoints();
    RingShowerWhiteBox::check_radius_endpoints();
    DynamoWhiteBox::check_overlapping_wipes_stay_in_range();
    test_dynamo_trail_ceiling_bounds_the_ring();
    test_dynamo_emitted_points_counts_ring_seeds();
    test_ringspin_trail_hugs_its_great_circles();
    test_hopf_projection_math();
    test_hopf_trail_trim_keeps_a_segment();
    test_raymarch_constexpr_sqrt_converges();
    test_raymarch_unit_bounds_contains_twisted_tube();
    test_gnomonicstars_radius_px_spans_one_column();
    test_gnomonicstars_spiral_cache_invalidation();
    test_petalflow_spawn_gap_bounded();
    test_displacement_field_lazy_hue_table_matches_eager();
    test_displacement_field_hue_table_fidelity();
    test_displacement_field_hue_table_frame_fidelity();
    test_displacement_field_zero_hue_scale_is_exact();
    test_displacement_field_clip_tiles_full();
    test_mindsplatter_emit_phase_wrapped();
    test_shaderball_glitch_lens_unit_norm();
    test_mobius_rings_conformal_and_counter_rotation();
    test_islamicstars_seed_sprite_fade_in();
    test_islamicstars_recipe_build_smoke();
    test_islamicstars_roster_cycle_fits_budget();
    test_islamicstars_dual_bridge_fits_budget();
  }

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test
