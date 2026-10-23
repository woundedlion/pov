/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Per-effect sweep primitives (smoke_one, determinism_one,
 * clip_clear_parity_one) and the effects white-box suite.
 */
#pragma once

#include "core/animation/orientation.h"
#include "targets/effects.h"
#include "core/render/canvas.h"
#include "core/render/sdf/volume.h"
#include "core/memory.h"
#include "hardware/pov_segment_map.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/temporal_oracle.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace effects_tests {

/**
 * @brief Primary production render width in pixels.
 */
constexpr int DEFAULT_W = 288;
/**
 * @brief Primary production render height in pixels.
 * @details Paired with DEFAULT_W for the full-sphere production resolution.
 */
constexpr int DEFAULT_H = 144;

/**
 * @brief Per-effect smoke frame count, resolved from HS_SMOKE_FRAMES.
 */
using hs_test::smoke_frames;

/**
 * @brief Selects the effects test depth tier from the environment.
 * @return true for the FULL suite (the production-resolution roster passes
 * and the white-box cases in the FULL block), false for the QUICK tier.
 * @details HS_EFFECTS_FULL=1 selects the FULL tier; QUICK is the default.
 */
inline bool effects_full_suite() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *e = std::getenv("HS_EFFECTS_FULL"))
#pragma clang diagnostic pop
    return std::atoi(e) != 0;
  return false;
}

inline void lint_dead_sliders(Effect &effect, const char *name);

/**
 * @brief Forward declaration of the animated-param pause lint.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in diagnostic output.
 */
inline void lint_animated_pause(Effect &effect, const char *name);

/**
 * @brief Drives one effect type through construct -> init -> render -> read-back.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param name Effect name used in the [ok] / diagnostic output.
 * @details Renders smoke_frames() frames and reads back every pixel, checks
 * get_pixel aliases the displayed buffer, runs the parameter lints at
 * <SMALL_W,SMALL_H>, and requires no dropped timeline event. The pause lint
 * leaves the effect paused.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void smoke_one(const char *name) {
  HS_CONTEXT(name);
  reset_effect_globals();
  const uint32_t dropped_before = Timeline::dropped_events();
  pin_frame_clock(0);

  E<W, H> effect;
  effect.init();

  HS_EXPECT_EQ(effect.width(), W);
  HS_EXPECT_EQ(effect.height(), H);

  HS_EXPECT_EQ(effect.clip().margin, ClipRegion{}.margin);

  const int frames = smoke_frames();
  uint64_t previous_hash = 0;
  bool motion = false;
  for (int f = 0; f < frames; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    // Consume the queued frame before the next Canvas constructor watchdog expires.
    effect.advance_display();
    const uint64_t hash = frame_hash<W, H>(effect);
    if (f > frames / 2 && hash != previous_hash)
      motion = true;
    previous_hash = hash;
  }
  if (frames >= 4) {
    if (!motion)
      std::printf("  STATIC %-20s had no motion in the final %d frames\n", name,
                  frames - frames / 2 - 1);
    HS_EXPECT(motion, "effect must change output after warmup");
  }

  const uint64_t acc = frame_energy<W, H>(effect);

  std::printf("  [ok] %-20s rendered %d frames @ %dx%d (sum=%llu)\n", name,
              frames, W, H, static_cast<unsigned long long>(acc));

  {
    if (acc == 0)
      std::printf("  ALL-BLACK %-20s final frame (of %d) has no lit pixel "
                  "@ %dx%d\n",
                  name, frames, W, H);
    HS_EXPECT(acc > 0, "effect must produce non-black output");
  }

  // get_pixel must stay the row-major view of display_buffer() after a flip.
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

  if constexpr (W == SMALL_W && H == SMALL_H) {
    lint_dead_sliders(effect, name);
    lint_animated_pause(effect, name);
  }

  // The drop counter is process-wide and wraps at uint32_t; take the
  // per-effect delta.
  const uint32_t dropped = Timeline::dropped_events() - dropped_before;
  if (dropped != 0)
    std::printf("  TIMELINE FULL %-20s dropped %u animation(s) over %d frames "
                "@ %dx%d (Timeline::MAX_EVENTS is %d)\n",
                name, static_cast<unsigned>(dropped), frames, W, H,
                Timeline::MAX_EVENTS);
  HS_EXPECT_EQ(dropped, 0u);
}

#include "tests/effects/parameter_probe.h"

/**
 * @brief Runtime "registered-but-unread" lint for the live-art param system.
 * @param effect Effect instance whose editable params are probed.
 * @param name Effect name used in DEAD SLIDER diagnostic output.
 * @details A registered param that is neither animated nor mark_readonly()
 * must keep a value written through updateParameter() across frames. Detects
 * per-frame overwrites, not visual influence.
 */
inline void lint_dead_sliders(Effect &effect, const char *name) {
  for (const auto &def : effect.getParameters()) {
    HS_CONTEXT(def.name);
    if (def.animated || def.readonly)
      continue;
    const float range = def.max - def.min;
    if (range <= 0.0f)
      continue;
    const float cur = def.get();
    // An in-range target well clear of the current value, so a revert is
    // visible; a bool has only its flipped value.
    float target = parameter_probe_target(def, cur);
    // Dropdowns normalize to integral option indices.
    if (def.is_enum() || def.is_integer()) {
      target =
          cur < (def.min + def.max) * 0.5f ? ceilf(target) : floorf(target);
      if (target == cur)
        continue;
    }
    HS_EXPECT_EQ(effect.updateParameter(def.name, target),
                 ParamSetResult::APPLIED);
    for (int f = 0; f < 3; ++f) {
      effect.draw_frame();
      effect.advance_display();
    }
    // Require the value near `target` AND strictly closer to it than to the
    // pre-write `cur`, catching a slow per-frame revert. A bool must read back
    // exactly.
    const float eps = def.is_bool() ? 0.0f : fmaxf(1e-3f, 1e-3f * range);
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
    HS_EXPECT_EQ(effect.updateParameter(def.name, cur),
                 ParamSetResult::APPLIED);
  }
}

/**
 * @brief Checks paused rendering for advertised animated parameters.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in PAUSE LEAK diagnostic output.
 * @details Parameters are audited one at a time over PAUSE_AUDIT_FRAMES in
 *          total.
 */
inline void lint_animated_pause(Effect &effect, const char *name) {
  std::vector<const char *> names;
  std::vector<float> original;
  std::vector<float> target;
  names.reserve(effect.getParameters().size());
  original.reserve(effect.getParameters().size());
  target.reserve(effect.getParameters().size());
  for (const auto &def : effect.getParameters()) {
    if (!def.animated || def.readonly)
      continue;
    // The write below lands on the requested value, so a schema-driven effect
    // whose rendered slot is canonicalized away from it only round-trips
    // through get_requested().
    const float current = def.get_requested();
    names.push_back(def.name);
    original.push_back(current);
    target.push_back(parameter_probe_target(def, current));
  }
  const size_t count = names.size();
  if (count == 0)
    return;

  constexpr int PAUSE_AUDIT_FRAMES = 500;
  const int frames_per_param =
      std::max(4, static_cast<int>((PAUSE_AUDIT_FRAMES + count - 1) / count));
  bool leaked = false;
  for (size_t index = 0; index < count; ++index) {
    HS_CONTEXT(names[index]);
    for (size_t restore = 0; restore < count; ++restore)
      HS_EXPECT_EQ(effect.updateParameter(names[restore], original[restore]),
                   ParamSetResult::APPLIED);
    effect.setAnimationsPaused(false);
    HS_EXPECT_EQ(effect.updateParameter(names[index], target[index]),
                 ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(effect.animations_paused());
    const auto *written_def = effect.getParameters().find(names[index]);
    const float WRITTEN = written_def->get_requested();
    const float EPS =
        written_def->is_bool()
            ? 0.0f
            : fmaxf(1e-3f, 1e-3f * (written_def->max - written_def->min));
    effect.draw_frame();
    effect.advance_display();
    const auto *after_write = effect.getParameters().find(names[index]);
    const bool PERSISTED = fabsf(after_write->get_requested() - WRITTEN) <= EPS;
    if (!PERSISTED)
      std::printf(
          "  PAUSED WRITE LOST %s::%s wrote %.4f, engine reverted to %.4f\n",
          name, names[index], static_cast<double>(WRITTEN),
          static_cast<double>(after_write->get_requested()));
    HS_EXPECT(PERSISTED, "paused animated parameter writes must persist");
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
 * @param frame_fold Optional out-param receiving an FNV-1a checksum folded over
 * every displayed frame. Ignored when nullptr.
 * @param lit Optional output: whether any displayed frame contains a lit pixel.
 * @details Resets every shared global the smoke path does (RNG seed, arenas,
 * Timeline, pole-LOD knob, scan counters) and pins the mock clock to the frame
 * cadence, so two calls start from an identical state.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void render_capture(std::vector<Pixel> &out, int frames,
                           uint64_t *frame_fold = nullptr,
                           bool *lit = nullptr) {
  reset_effect_globals();
  // Pre-init epoch, so construction sees the same clock on both runs.
  pin_frame_clock(0);

  if (lit)
    *lit = false;
  uint64_t fold = hs_test::HASH_SEED64;

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    if (frame_fold || lit) {
      if (lit && count_lit_region<W, H>(effect) != 0)
        *lit = true;
      fold = frame_hash<W, H>(effect, fold);
    }
  }
  if (frame_fold)
    *frame_fold = fold;

  capture_frame<W, H>(effect, out);
}

/**
 * @brief Folds every displayed frame of an effect rendered with temporal
 *        sampling disabled.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param frames Frames to render.
 * @param clip Display band {y0, y1, x0, x1} set before the first frame, or
 *        null for the full canvas.
 * @return FNV-1a fold over every channel of every displayed frame.
 */
template <template <int, int> class E, int W, int H>
inline uint64_t untimed_frame_fold(int frames,
                                   const std::array<int, 4> *clip = nullptr) {
  reset_effect_globals();
  pin_frame_clock(0);
  Spacetime::set_flat_surface(true);
  uint64_t fold = hs_test::HASH_SEED64;
  {
    E<W, H> effect;
    effect.init();
    if (clip)
      effect.set_clip((*clip)[0], (*clip)[1], (*clip)[2], (*clip)[3]);
    for (int f = 0; f < frames; ++f) {
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      fold = frame_hash<W, H>(effect, fold);
    }
  }
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return fold;
}

/**
 * @brief Scrambles every output-affecting global that render_capture() resets.
 * @details Runs between captures to exercise recovery from dirty process state.
 */
inline void perturb_determinism_globals() {
  hs::random().seed(0xC0FFEEu);
  global_timeline_t = 0x5EED;
  Render::pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT + 0.5f;
  configure_arenas(DEFAULT_PERSISTENT_SIZE - 32, DEFAULT_SCRATCH_A_SIZE + 16,
                   DEFAULT_SCRATCH_B_SIZE + 16);
  for (Arena *arena : {&persistent_arena, &scratch_arena_a, &scratch_arena_b}) {
    const size_t bytes = arena->get_capacity();
    std::memset(arena->allocate(bytes), 0xA5, bytes);
  }
  HS_SCAN_METRIC(++hs::g_scan_metrics.pixels_tested);
}

/** @brief Frames per segment in the clip-clear parity sweep. */
constexpr int PARITY_FRAMES = 16;
/** @brief Arm segments walked by the clip-clear parity sweep. */
constexpr int PARITY_SEGMENTS = 4;

/**
 * @brief Cross-run determinism test: renders an effect twice and requires byte-identical frames.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels (defaults to DEFAULT_W).
 * @tparam H Render height in pixels (defaults to DEFAULT_H).
 * @param name Effect name used in the NONDETERMINISTIC diagnostic output.
 * @details The clock seam neutralizes wall-time, so a divergence is real
 * nondeterminism (uninitialized read or stale global).
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void determinism_one(const char *name) {
  HS_CONTEXT(name);
  const int frames = smoke_frames();
  std::vector<Pixel> a, b;
  uint64_t fold_a = 0, fold_b = 0;
  bool lit_a = false, lit_b = false;
  render_capture<E, W, H>(a, frames, &fold_a, &lit_a);
  perturb_determinism_globals();
  render_capture<E, W, H>(b, frames, &fold_b, &lit_b);
  hs::clear_mock_time();

  // Per-frame fold catches mid-run divergence that reconverges by the final
  // frame.
  if (fold_a != fold_b)
    std::printf("  NONDETERMINISTIC %-20s per-frame checksum %llu != %llu over "
                "%d frames\n",
                name, static_cast<unsigned long long>(fold_a),
                static_cast<unsigned long long>(fold_b), frames);
  HS_EXPECT(fold_a == fold_b,
            "effect must render identically every frame across runs under a "
            "fixed clock");

  HS_EXPECT_SIZE_OR_RETURN(b, a.size());

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
  HS_EXPECT_TRUE(lit_a);
  HS_EXPECT_TRUE(lit_b);
}

/**
 * @brief Drives one effect under a moving segment clip with each clear scope and
 *        requires the displayed pixels to agree.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param name Effect name used in the diagnostic output.
 * @param frames Frames rendered for each segment and clear scope.
 * @param probe Optional observer called after each displayed frame with the effect,
 *              segment, frame index, and full-clear flag.
 * @param require_lit_per_segment Require a nonblack reference in every segment.
 * @details The clip clear leaves the region outside the display band holding
 *          pixels from the frame that last wrote this buffer. Anything an effect
 *          reads back or carries forward from there shows up here as a mismatch
 *          against the whole-buffer clear.
 */
template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void clip_clear_parity_one(
    const char *name, int frames = PARITY_FRAMES,
    std::type_identity_t<std::function<void(const E<W, H> &, int, int, bool)>>
        probe = {},
    bool require_lit_per_segment = false) {
  HS_CONTEXT(name);
  constexpr int S = H * 2;

  auto render = [&](int segment_id, bool full_clear) {
    reset_effect_globals();
    pin_frame_clock(0);
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
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      if (probe)
        probe(effect, segment_id, f, full_clear);
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

    size_t different = 0, segment_lit = 0;
    for (size_t i = 0; i < full.size(); ++i) {
      if (full[i] != clipped[i])
        ++different;
      if (full[i].r | full[i].g | full[i].b)
        ++segment_lit;
    }
    lit += segment_lit;
    if (require_lit_per_segment)
      HS_EXPECT_GT(segment_lit, static_cast<size_t>(0));
    displayed_pixels += full.size();
    if (different)
      std::printf("  CLIP-CLEAR DRIFT %-20s segment %d: %zu of %zu displayed "
                  "pixels differ from the full-buffer clear\n",
                  name, segment_id, different, full.size());
    HS_EXPECT(different == 0,
              "clip clearing must not change any displayed pixel");
  }
  // Two all-black renders agree trivially; require output across the sweep.
  if (lit == 0)
    std::printf("  CLIP-CLEAR DARK %-20s no lit pixel over %zu displayed "
                "pixels\n",
                name, displayed_pixels);
  HS_EXPECT(lit > 0, "clip-clear parity must compare a lit render");
  hs::clear_mock_time();
}

/** @brief Adjacent-column colour steps at the arm boundaries and mid-arm. */
struct MeshSeamSteps {
  double mean_arm_ratios = 0.0; ///< Sum of mean arm-boundary ratios.
  double arm_ratios = 0.0; ///< Sum of the larger arm-boundary ratio per frame.
  double interior_ratios = 0.0; ///< Sum of the mid-arm ratios.
  int frames = 0;               ///< Consecutive frame pairs measured.

  /**
   * @brief Accumulates one pair of consecutive frames.
   * @param previous Previous displayed frame.
   * @param current Current displayed frame.
   * @param width Frame width.
   * @param height Frame height.
   */
  void add(const std::vector<Pixel> &previous,
           const std::vector<Pixel> &current, int width, int height) {
    const auto column_difference = [&](const auto &a, int x, const auto &b,
                                       int z) {
      double sum = 0.0;
      for (int y = 0; y < height; ++y) {
        const Pixel &p = a[size_t(y) * width + (x + width) % width];
        const Pixel &q = b[size_t(y) * width + (z + width) % width];
        sum += std::abs(int(p.r) - q.r) + std::abs(int(p.g) - q.g) +
               std::abs(int(p.b) - q.b);
      }
      return sum;
    };
    double arm = 0.0, interior = 0.0, interior_neighbors = 0.0;
    for (int boundary : {0, width / 2}) {
      double neighbors = 0.0;
      for (int j = 1; j <= 4; ++j) {
        neighbors += column_difference(previous, boundary - j - 1, previous,
                                       boundary - j);
        neighbors +=
            column_difference(current, boundary + j - 1, current, boundary + j);
      }
      const double seam =
          column_difference(previous, boundary - 1, current, boundary);
      if (neighbors > 0.0) {
        const double ratio = seam / (neighbors / 8.0);
        arm = std::max(arm, ratio);
        mean_arm_ratios += ratio * 0.5;
      }
      const int middle = boundary + width / 4;
      interior += column_difference(current, middle - 1, current, middle);
      for (int j : {-3, -2, 1, 2})
        interior_neighbors +=
            column_difference(current, middle - 1 + j, current, middle + j);
    }
    arm_ratios += arm;
    if (interior_neighbors > 0.0)
      interior_ratios += interior / (interior_neighbors / 4.0);
    ++frames;
  }
};

/**
 * @brief Checks boundary steps against their neighboring steps.
 * @param steps Accumulated frame differences.
 * @param label Effect name for the result.
 * @param mean_bound Bound on the mean arm-boundary ratio.
 * @param arm_bound Bound on the mean of each frame's larger arm ratio.
 * @param interior_bound Bound on the mean mid-arm ratio.
 */
inline void expect_mesh_seam_continuous(const MeshSeamSteps &steps,
                                        const char *label,
                                        double mean_bound = 1.06,
                                        double arm_bound = 1.15,
                                        double interior_bound = 1.06) {
  HS_CONTEXT(label);
  std::printf("  [%s] mean seam step %.3f x neighbors; mid-arm %.3f\n", label,
              steps.arm_ratios / steps.frames,
              steps.interior_ratios / steps.frames);
  std::printf("  mean across arms %.3f\n",
              steps.mean_arm_ratios / steps.frames);
  HS_EXPECT_GT(steps.frames, 0);
  HS_EXPECT_LE(steps.mean_arm_ratios / steps.frames, mean_bound);
  HS_EXPECT_LE(steps.arm_ratios / steps.frames, arm_bound);
  HS_EXPECT_LE(steps.interior_ratios / steps.frames, interior_bound);
}

/** @brief One arm side's copy of a face in a mesh deformation scan. */
struct FaceCopy {
  std::vector<uint16_t> indices;     ///< The face's vertex indices.
  std::vector<math::Vector> polygon; ///< Its scanned vertices, model space.
  float size = 0.0f;                 ///< Inradius the scan shaded it with.
};

/** @brief A deformation scan's face copies: first arm, then second, by face. */
using FaceCopies = std::array<std::map<size_t, FaceCopy>, 2>;

/// Largest extrapolation of a nonpolar face's contiguous capture branch.
inline constexpr float SEAM_EXTRAPOLATION = 1.0f;

/**
 * @brief Records the faces a deformation scan draws on one arm.
 * @tparam CAP Orientation history capacity.
 * @param canvas Canvas inside the scan.
 * @param mesh The scanned mesh, at the scan's captured vertices.
 * @param orientation Camera at the current endpoint.
 * @param copies Receives the arm's copies.
 * @return The arm, or -1 for a scan its render-x window does not hold to
 *         one arm.
 */
template <int CAP>
inline int record_face_copies(const Canvas &canvas, const MeshState &mesh,
                              const math::Orientation<CAP> &orientation,
                              FaceCopies &copies) {
  const Spacetime::TimeWindow &window = canvas.time_window();
  const Scan::ColumnBasis *time = canvas.column_basis();
  const int arm = window.arm(canvas.width());
  if (!window.active || arm < 0 || !time || !time->mesh.selected_faces)
    return -1;
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    if (time->mesh.selected_faces[f] != Scan::detail::FacePass::ARM_COPY)
      continue;
    FaceCopy &copy = copies[arm][f];
    copy = FaceCopy{};
    if (const Scan::FaceArmCopy *own = time->mesh.face_copy) {
      for (size_t k = 0; k < own->count; ++k) {
        copy.indices.push_back(own->sources[k]);
        copy.polygon.push_back(orientation.unorient(own->vertices[k]));
      }
      continue;
    }
    const uint16_t *face =
        mesh.get_faces_data() + mesh.get_face_offsets_data()[f];
    for (size_t k = 0; k < mesh.get_face_counts_data()[f]; ++k) {
      copy.indices.push_back(face[k]);
      copy.polygon.push_back(orientation.unorient(mesh.vertices[face[k]]));
    }
  }
  return arm;
}

/**
 * @brief Renders face copies as a still mesh through a camera.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam SetupFn Callable (size_t face, float size, FacePaletteShader &)
 *         shading a face as its effect does.
 * @param effect Effect whose display receives the render.
 * @param copies Face copies, by face.
 * @param faces The faces to draw.
 * @param camera Camera the copies are drawn through.
 * @param setup Shades each drawn face.
 * @param ids Whether to draw each face flat in its position in @p faces plus
 *        one, in the red channel, instead of shading it.
 * @return Displayed frame.
 */
template <int W, int H, typename EffectT, typename SetupFn>
inline std::vector<Pixel>
render_face_copies(EffectT &effect, const std::map<size_t, FaceCopy> &copies,
                   const std::vector<size_t> &faces,
                   const math::Quaternion &camera, SetupFn &&setup, bool ids) {
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(effect);
    ScratchScope guard(scratch_arena_a);
    size_t count = 0;
    for (size_t f : faces)
      count += copies.at(f).polygon.size();
    MeshState mesh;
    mesh.vertices.bind(scratch_arena_a, count);
    mesh.faces.bind(scratch_arena_a, count);
    mesh.face_counts.bind(scratch_arena_a, faces.size());
    mesh.face_offsets.bind(scratch_arena_a, faces.size());
    for (size_t f : faces) {
      const FaceCopy &copy = copies.at(f);
      mesh.face_offsets.push_back(static_cast<uint16_t>(mesh.faces.size()));
      mesh.face_counts.push_back(static_cast<uint8_t>(copy.polygon.size()));
      for (const math::Vector &v : copy.polygon) {
        mesh.faces.push_back(static_cast<uint16_t>(mesh.vertices.size()));
        mesh.vertices.push_back(math::rotate(v, camera));
      }
    }
    if (ids) {
      const auto id = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(static_cast<uint16_t>(f.v2 + 1.5f), 0, 0), 1);
      };
      Scan::Mesh::draw<W, H>(pipeline, canvas, mesh, id, scratch_arena_b);
    } else {
      FacePaletteShader shader;
      auto shade = [&](size_t j, float size) { setup(faces[j], size, shader); };
      Scan::Mesh::draw_specialized<W, H>(pipeline, canvas, mesh, shader,
                                         scratch_arena_b, nullptr, shade);
    }
  }
  effect.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(effect, out);
  return out;
}

/**
 * @brief Largest colour change a face's shading allows between two copies.
 * @tparam SetupFn As render_face_copies() takes it.
 * @param face Face index.
 * @param size The face's inradius in one copy.
 * @param shift Bound on how far the copies' inradius, and their edge depth
 *        at a pixel, differ.
 * @param setup Shades the face as its effect does.
 * @return 16-bit channel units: the shading at each edge depth against the
 *         shading within @p shift of that depth and of @p size.
 */
template <typename SetupFn>
inline int face_colour_bound(size_t face, float size, float shift,
                             SetupFn &&setup) {
  constexpr int DEPTHS = 64;
  constexpr float STEPS[] = {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f};
  FacePaletteShader one;
  setup(face, size, one);
  int bound = 0;
  for (float r : STEPS) {
    FacePaletteShader other;
    setup(face, fmaxf(size + r * shift, 0.0f), other);
    for (int k = 0; k <= DEPTHS; ++k) {
      Fragment a;
      a.v1 = -2.0f * size * k / DEPTHS;
      one(a.pos, a);
      for (float s : STEPS) {
        Fragment b;
        b.v1 = a.v1 + s * shift;
        other(b.pos, b);
        const Pixel &p = a.color.color, &q = b.color.color;
        bound = std::max({bound, std::abs(int(p.r) - q.r),
                          std::abs(int(p.g) - q.g), std::abs(int(p.b) - q.b)});
      }
    }
  }
  return bound;
}

/// Columns either side of each arm boundary whose colours a seam check reads.
inline constexpr int SEAM_CHECK_COLUMNS = 3;

/**
 * @brief Checks that one frame's copies before an arm boundary and the next
 *        frame's copies after it are one polygon there that shades alike.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam ToleranceFn Callable (uint16_t vertex) returning the bound on that
 *         vertex's gap between the copies (radians of chord).
 * @tparam SetupFn As render_face_copies() takes it.
 * @param effect Effect whose display renders the copies.
 * @param label Draw named in failure context.
 * @param seam The boundary's first column after it: 0 or W / 2.
 * @param late Frame N's copies on the arm the boundary ends.
 * @param early Frame N+1's copies on the arm the boundary starts.
 * @param camera Camera the copies' colours are compared through.
 * @param tolerance Per-vertex gap bound.
 * @param setup Shades each face as its effect does.
 * @return The number of faces compared.
 * @details Compares the faces the boundary's columns reach. A face's edge
 * depth and inradius move by at most twice its largest vertex gap, the
 * gnomonic stretch of a face within a quarter turn of its centre. Colours
 * are compared on the boundary columns' pixels a pixel clear of every edge.
 * Faces around a pole draw once and leave no copies.
 */
template <int W, int H, typename EffectT, typename ToleranceFn,
          typename SetupFn>
inline size_t
expect_seam_copies_match(EffectT &effect, const char *label, int seam,
                         const std::map<size_t, FaceCopy> &late,
                         const std::map<size_t, FaceCopy> &early,
                         const math::Quaternion &camera,
                         ToleranceFn &&tolerance, SetupFn &&setup) {
  HS_CONTEXT(label);
  HS_CONTEXT("seam column", seam);
  constexpr float PIXEL = math::TWO_PI_F / W;
  // Sign of x along the boundary's meridian.
  const float side = seam == 0 ? 1.0f : -1.0f;
  std::vector<size_t> faces;
  std::vector<int> colour_bounds;
  float worst_gap = 0.0f, worst_ratio = 0.0f, worst_size = 0.0f;
  for (const auto &[f, before] : late) {
    const auto after = early.find(f);
    if (after == early.end() || before.size <= 0.0f ||
        after->second.size <= 0.0f)
      continue;
    const auto seen = [&](const FaceCopy &copy) {
      std::vector<math::Vector> out;
      for (const auto &v : copy.polygon)
        out.push_back(math::rotate(v, camera));
      return out;
    };
    const std::vector<math::Vector> a = seen(before), b = seen(after->second);
    // A face meets the boundary when an edge crosses its meridian, or when
    // it lies within the checked columns on the boundary's side.
    const auto touches_seam = [&](const std::vector<math::Vector> &polygon) {
      float lo = 2.0f, hi = -2.0f, middle = 0.0f;
      for (size_t k = 0; k < polygon.size(); ++k) {
        const math::Vector &p = polygon[k];
        const math::Vector &q = polygon[(k + 1) % polygon.size()];
        lo = fminf(lo, p.z);
        hi = fmaxf(hi, p.z);
        middle += p.x;
        if (p.z * q.z <= 0.0f && p.z != q.z &&
            (p.x + (q.x - p.x) * p.z / (p.z - q.z)) * side > 0.0f)
          return true;
      }
      const float fringe = sinf(SEAM_CHECK_COLUMNS * PIXEL);
      return lo <= fringe && hi >= -fringe && middle * side > 0.0f;
    };
    if (!touches_seam(a) && !touches_seam(b))
      continue;
    // A face copy starts at its first vertex on its arm.
    const size_t n = before.indices.size();
    size_t offset = n;
    for (size_t o = 0; o < n && offset == n; ++o) {
      bool same = after->second.indices.size() == n;
      for (size_t k = 0; k < n && same; ++k)
        same = before.indices[k] == after->second.indices[(k + o) % n];
      if (same)
        offset = o;
    }
    HS_EXPECT_LT(offset, n);
    if (offset == n)
      continue;
    float face_tolerance = 0.0f;
    for (size_t k = 0; k < n; ++k) {
      const size_t j = (k + offset) % n;
      const float gap = sqrtf(
          math::distance_squared(before.polygon[k], after->second.polygon[j]));
      const float bound = tolerance(before.indices[k]);
      HS_EXPECT_LE(gap, bound);
      worst_gap = fmaxf(worst_gap, gap);
      worst_ratio = fmaxf(worst_ratio, gap / bound);
      face_tolerance = fmaxf(face_tolerance, bound);
    }
    const float size_gap = fabsf(before.size - after->second.size);
    HS_EXPECT_LE(size_gap, 2.0f * face_tolerance);
    worst_size = fmaxf(worst_size, size_gap / before.size);
    faces.push_back(f);
    colour_bounds.push_back(
        face_colour_bound(f, before.size, 2.0f * face_tolerance, setup));
  }
  const auto render = [&](const std::map<size_t, FaceCopy> &copies, bool ids) {
    return render_face_copies<W, H>(effect, copies, faces, camera, setup, ids);
  };
  const std::vector<Pixel> left = render(late, false);
  const std::vector<Pixel> right = render(early, false);
  const std::vector<Pixel> left_id = render(late, true);
  const std::vector<Pixel> right_id = render(early, true);
  const auto face_at = [&](int x, int y) {
    const int id = left_id[static_cast<size_t>(y) * W + x].r;
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx) {
        const int row = y + dy;
        if (row < 0 || row >= H)
          return 0;
        const size_t j = static_cast<size_t>(row) * W + (x + dx + W) % W;
        if (left_id[j].r != id || right_id[j].r != id)
          return 0;
      }
    return id;
  };
  int colour = 0, excess = 0, interior = 0;
  for (int y = 0; y < H; ++y)
    for (int c = -SEAM_CHECK_COLUMNS; c < SEAM_CHECK_COLUMNS; ++c) {
      const int x = (seam + c + W) % W;
      const int id = face_at(x, y);
      if (id == 0)
        continue;
      ++interior;
      const size_t i = static_cast<size_t>(y) * W + x;
      const int difference = std::max({std::abs(int(left[i].r) - right[i].r),
                                       std::abs(int(left[i].g) - right[i].g),
                                       std::abs(int(left[i].b) - right[i].b)});
      colour = std::max(colour, difference);
      excess = std::max(excess, difference - colour_bounds[id - 1]);
    }
  std::printf("  [%s, seam %d] %zu faces: vertex gap %.4f px "
              "(%.3f of bound), inradius %.5f, colour %.3f steps over %d "
              "interior pixels\n",
              label, seam, faces.size(), worst_gap / PIXEL, worst_ratio,
              worst_size, colour / 257.0, interior);
  HS_EXPECT_LE(excess, 0);
  return interior > 0 ? faces.size() : 0;
}

/**
 * @brief expect_seam_copies_match() at both arm boundaries.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam ToleranceFn As expect_seam_copies_match() takes it.
 * @tparam SetupFn As render_face_copies() takes it.
 * @param effect Effect whose display renders the copies.
 * @param label Draw named in failure context.
 * @param frames Frame N's and frame N+1's copies.
 * @param camera Camera at the frames' shared endpoint.
 * @param tolerance Per-vertex gap bound.
 * @param setup Shades each face as its effect does.
 * @return The number of faces compared, over both boundaries.
 */
template <int W, int H, typename EffectT, typename ToleranceFn,
          typename SetupFn>
inline size_t expect_seams_match(EffectT &effect, const char *label,
                                 const std::array<FaceCopies, 2> &frames,
                                 const math::Quaternion &camera,
                                 ToleranceFn &&tolerance, SetupFn &&setup) {
  return expect_seam_copies_match<W, H>(effect, label, W / 2, frames[0][0],
                                        frames[1][1], camera, tolerance,
                                        setup) +
         expect_seam_copies_match<W, H>(effect, label, 0, frames[0][1],
                                        frames[1][0], camera, tolerance, setup);
}

/**
 * @brief Rotation angle between two unit quaternions.
 * @param a First rotation.
 * @param b Second rotation.
 * @return Radians in [0, PI].
 */
inline float quaternion_angle(const math::Quaternion &a,
                              const math::Quaternion &b) {
  const float cosine = fabsf(math::dot(a, b));
  return 2.0f * acosf(fminf(1.0f, cosine));
}

/**
 * @brief Bound on how far two copies' capture times of one vertex differ.
 * @param motion The vertex's motion over either frame (radians).
 * @param turn Largest camera rotation over either frame (radians).
 * @param horizontal The vertex's distance from the polar axis.
 * @return Frame fraction: math::fast_atan2's peak error and one correction
 *         step's residual in each copy, and the column basis the
 *         extrapolating copy holds at its arm's end column.
 */
inline float seam_capture_slack(float motion, float turn, float horizontal) {
  constexpr float FAST_ATAN2_ERROR = 0.0038f;
  const float ratio = (motion + turn) / (math::PI_F * horizontal);
  return 2.0f * FAST_ATAN2_ERROR / math::PI_F +
         2.0f * (1.0f + SEAM_EXTRAPOLATION) * ratio * ratio +
         SEAM_EXTRAPOLATION * turn / (math::PI_F * horizontal);
}

/**
 * @brief Bound on a vertex's gap between frame N's late copy and frame N+1's
 *        early copy.
 * @tparam W Canvas width.
 * @tparam LateFn Callable (float u) giving the vertex on frame N's path.
 * @tparam EarlyFn Callable (float u) giving the vertex on frame N+1's path.
 * @param late Frame N's path, model space.
 * @param early Frame N+1's path, model space.
 * @param camera Camera at the frames' shared endpoint.
 * @param turn Largest camera rotation over either frame (radians).
 * @return The paths' largest gap within SEAM_EXTRAPOLATION of the shared
 *         endpoint, which is second order in the motion, plus the vertex's
 *         motion over seam_capture_slack(), one column interval, and any
 *         endpoint-held motion.
 */
template <int W, typename LateFn, typename EarlyFn>
inline float seam_vertex_tolerance(LateFn &&late, EarlyFn &&early,
                                   const math::Quaternion &camera, float turn) {
  constexpr int SAMPLES = 8;
  constexpr float FLOAT_ROUNDING = 1e-6f;
  const auto gap = [](const math::Vector &a, const math::Vector &b) {
    return sqrtf(math::distance_squared(a, b));
  };
  float paths = 0.0f;
  for (int s = -SAMPLES; s <= SAMPLES; ++s) {
    const float d = SEAM_EXTRAPOLATION * s / SAMPLES;
    paths = fmaxf(paths, gap(late(1.0f + d), early(d)));
  }
  const float before = gap(late(0.0f), late(1.0f));
  const float after = gap(early(0.0f), early(1.0f));
  const float motion = fmaxf(before, after);
  const bool held = fminf(before, after) <= Scan::HELD_VERTEX_ARC<W>;
  const float held_slack = held ? motion : 0.0f;
  const math::Vector seen = math::rotate(late(1.0f), camera);
  const float horizontal = sqrtf(seen.x * seen.x + seen.z * seen.z);
  return held_slack + paths +
         (seam_capture_slack(motion, turn, horizontal) + 2.0f / W) * motion +
         FLOAT_ROUNDING;
}

#include "tests/effects/spherical_harmonics.h"
#include "tests/effects/reaction_diffusion_gs.h"
#include "tests/effects/reaction_diffusion_bz.h"
#include "tests/effects/hankin.h"
#include "tests/effects/dreamballs_cases.h"
#include "tests/effects/mesh_feedback.h"
#include "tests/effects/numeric_invariants.h"
#include "tests/effects/effect_whitebox.h"
#include "tests/effects/galaxies_tests.h"
#include "tests/effects/displacement_field.h"
#include "tests/effects/shape_shifter.h"
#include "tests/effects/voronoi_cases.h"
#include "tests/effects/islamic_stars.h"
#include "tests/effects/column_time_parity.h"
/**
 * @brief Measures capture and scan work of the large scanned meshes under the
 *        device's segment clip: the top half of the rows and alternating
 *        half-width column runs.
 */
inline void measure_mesh_capture_work() {
  constexpr int W = 288, H = 144;
  const auto measure = [](auto &effect, const char *name, int frames,
                          auto draw) {
    Scan::MeshCaptureCounters total, peak;
    size_t pixels = 0, shaded = 0, pixel_peak = 0;
    int measured = 0;
    for (int frame = 0; frame < frames; ++frame) {
      const int first = frame % 2 ? W / 2 : 0;
      effect.set_clip(0, H / 2, first, first + W / 2);
      pin_frame_clock(frame);
      Scan::mesh_capture_counters = {};
      HS_SCAN_METRIC(hs::g_scan_metrics.reset());
      if (!draw(frame))
        break;
      ++measured;
      const auto &work = Scan::mesh_capture_counters;
      const auto add = [](size_t &sum, size_t &top, size_t value) {
        sum += value;
        top = std::max(top, value);
      };
      add(total.vertices, peak.vertices, work.vertices);
      add(total.evaluations, peak.evaluations, work.evaluations);
      add(total.captures, peak.captures, work.captures);
      add(total.shared_faces, peak.shared_faces, work.shared_faces);
      add(total.seam_faces, peak.seam_faces, work.seam_faces);
      add(total.cross_faces, peak.cross_faces, work.cross_faces);
      add(total.polar_faces, peak.polar_faces, work.polar_faces);
      add(total.swept_faces, peak.swept_faces, work.swept_faces);
      add(total.face_copies, peak.face_copies, work.face_copies);
      add(total.split_vertices, peak.split_vertices, work.split_vertices);
      add(total.scans, peak.scans, work.scans);
      add(total.face_setups, peak.face_setups, work.face_setups);
      HS_SCAN_METRIC(add(pixels, pixel_peak, hs::g_scan_metrics.pixels_tested));
      HS_SCAN_METRIC(shaded += hs::g_scan_metrics.shade_candidates);
    }
    std::printf(
        "  [capture work %s] frames=%d, total/peak per frame:\n"
        "    moving vertices %zu/%zu evaluations %zu/%zu captures %zu/%zu\n"
        "    shared faces %zu/%zu seam copies %zu/%zu cross %zu/%zu polar "
        "%zu/%zu (swept %zu/%zu) face copies %zu/%zu split vertices "
        "%zu/%zu\n"
        "    scans %zu/%zu face setups %zu/%zu pixels tested %zu/%zu "
        "shaded %zu\n",
        name, measured, total.vertices, peak.vertices, total.evaluations,
        peak.evaluations, total.captures, peak.captures, total.shared_faces,
        peak.shared_faces, total.seam_faces, peak.seam_faces, total.cross_faces,
        peak.cross_faces, total.polar_faces, peak.polar_faces,
        total.swept_faces, peak.swept_faces, total.face_copies,
        peak.face_copies, total.split_vertices, peak.split_vertices,
        total.scans, peak.scans, total.face_setups, peak.face_setups, pixels,
        pixel_peak, shaded);
  };
  for (const uint8_t node :
       {ConwayGraph::TRUNCATED_ICOSIDODECAHEDRON,
        ConwayGraph::TRUNCATED_ICOSAHEDRON, ConwayGraph::CUBOCTAHEDRON,
        ConwayGraph::TETRAHEDRON}) {
    reset_effect_globals();
    pin_frame_clock(0);
    HankinSolids<W, H> effect;
    effect.init();
    HankinSolidsWhiteBox::load_node(effect, node);
    HankinSolidsWhiteBox::rebuild(effect, 0.4f);
    const math::Vector axis = math::Vector(0.3f, 1.0f, 0.2f).normalized();
    char name[64];
    std::snprintf(name, sizeof(name), "Hankin %s",
                  Solids::simple_registry[node].name);
    measure(effect, name, 24, [&](int frame) {
      const auto &vertices = HankinSolidsWhiteBox::mesh_vertices(effect);
      std::vector<math::Vector> previous(vertices.begin(), vertices.end());
      HankinSolidsWhiteBox::rebuild(effect, 0.4f + 0.04f * (frame + 1));
      FaceCopies copies;
      HankinSolidsWhiteBox::draw_copies(
          effect, previous, math::make_rotation(axis, 0.02f * frame),
          math::make_rotation(axis, 0.02f * (frame + 1)), copies);
      return true;
    });
  }
  for (int shape : {6, 8}) {
    reset_effect_globals();
    pin_frame_clock(0);
    IslamicStars<W, H> effect;
    IslamicBuildProbe::set_trans_speed(effect, 4.0f);
    effect.init();
    IslamicBuildProbe::set_burst_size(effect, shape == 6 ? 8 : 1);
    IslamicBuildProbe::spawn_entry(effect,
                                   *IslamicBuildProbe::presets()[shape].params);
    measure(effect, shape == 6 ? "Islamic shape 6" : "Islamic shape 8", 120,
            [&](int) {
              effect.draw_frame();
              effect.advance_display();
              return IslamicBuildProbe::selected_entry(effect) ==
                     IslamicBuildProbe::presets()[shape].params;
            });
  }
  hs::clear_mock_time();
}

/// Frames of the device profile image's IslamicStars run.
inline constexpr int PROFILE_REPLAY_FRAMES = 280;

/**
 * @brief Replays the device profile image's IslamicStars run on the host and
 *        prints its deformation work per frame and per profile window.
 * @param preset Preset the image pins (HS_PROFILE_PRESET).
 * @param trans_speed The image's Trans Speed (HS_PROFILE_TRANS_SPEED).
 * @param frames Frames to draw, the first the full-width prime.
 * @param report Whether to print the work.
 * @details Mirrors the segmented driver on segment 0: the effect stream is
 * seeded from the effect's identity, the effect is built with its Trans
 * Speed, initialized, pinned and paused at @p preset, primed full width over
 * the segment's rows, then drawn one alternating half of them per frame.
 * Frame and window numbers match the image's `f` lines and windows. Run with
 * HS_MESH_CAPTURE_ONLY=3 (HS_REPLAY_PRESET, HS_REPLAY_TRANS_SPEED and
 * HS_REPLAY_FRAMES override 6, 4 and PROFILE_REPLAY_FRAMES).
 */
inline void replay_islamicstars_profile(size_t preset, float trans_speed,
                                        int frames, bool report = true) {
  constexpr int W = 288, H = 144, WINDOW = 16;
  reset_effect_globals();
  pin_frame_clock(0);
  hs::seed_effect_stream<IslamicStars<W, H>>();
  IslamicStars<W, H> effect;
  // The image compiles its Trans Speed into the constructor.
  IslamicBuildProbe::set_trans_speed(effect, trans_speed);
  effect.init();
  effect.profile_select_preset(preset);
  HS_EXPECT_TRUE(effect.updateParameter("Trans Speed", trans_speed) ==
                 ParamSetResult::APPLIED);
  Scan::MeshCaptureCounters window;
  [[maybe_unused]] size_t probes = 0;
  for (int frame = 1; frame <= frames; ++frame) {
    if (frame == 1)
      effect.set_clip(0, H / 2, 0, W);
    else
      effect.set_clip(0, H / 2, frame % 2 ? W / 2 : 0, frame % 2 ? W : W / 2);
    pin_frame_clock(frame);
    Scan::mesh_capture_counters = {};
    HS_SCAN_METRIC(hs::g_scan_metrics.reset());
    const size_t passes = IslamicBuildProbe::vertex_passes<W, H>();
    effect.draw_frame();
    effect.advance_display();
    const auto &w = Scan::mesh_capture_counters;
    HS_SCAN_METRIC(probes += hs::g_scan_metrics.pixels_tested);
    if (report && w.deformations)
      std::printf(
          "  f %d deform %zu moving %zu path evaluations %zu (scan %zu) "
          "captures %zu polar %zu swept %zu copies %zu\n",
          frame, w.deformations, w.vertices,
          IslamicBuildProbe::vertex_passes<W, H>() - passes, w.evaluations,
          w.captures, w.polar_faces, w.swept_faces, w.face_copies);
    window.deformations += w.deformations;
    window.vertices += w.vertices;
    window.evaluations += w.evaluations;
    window.polar_faces += w.polar_faces;
    window.swept_faces += w.swept_faces;
    window.face_copies += w.face_copies;
    if (report && frame % WINDOW == 0)
      std::printf("  window %d-%d: deform %zu polar %zu swept %zu copies %zu "
                  "moving %zu evaluations %zu probes %zu\n",
                  frame - WINDOW + 1, frame, window.deformations,
                  window.polar_faces, window.swept_faces, window.face_copies,
                  window.vertices, window.evaluations, probes);
    if (frame % WINDOW == 0) {
      window = {};
      probes = 0;
    }
  }
  if (report)
    std::printf("  scratch peak a %zu of %zu, b %zu of %zu\n",
                scratch_arena_a.get_lifetime_high_water_mark(),
                scratch_arena_a.get_capacity(),
                scratch_arena_b.get_lifetime_high_water_mark(),
                scratch_arena_b.get_capacity());
  hs::clear_mock_time();
}

/**
 * @brief Checks that one preset's profile run deforms within the effect's
 *        scratch arenas, which trap on overflow.
 * @param preset Preset the run pins.
 * @param frames Frames the run draws.
 */
inline void check_islamicstars_motion_fits_scratch(size_t preset, int frames) {
  HS_CONTEXT("preset", static_cast<int>(preset));
  scratch_arena_a.reset_high_water_mark();
  scratch_arena_b.reset_high_water_mark();
  replay_islamicstars_profile(preset, 4.0f, frames, false);
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(),
               scratch_arena_a.get_capacity());
  HS_EXPECT_LE(scratch_arena_b.get_high_water_mark(),
               scratch_arena_b.get_capacity());
  std::printf(
      "  [islamic scratch, preset %zu] a %zu of %zu, b %zu of %zu\n", preset,
      scratch_arena_a.get_high_water_mark(), scratch_arena_a.get_capacity(),
      scratch_arena_b.get_high_water_mark(), scratch_arena_b.get_capacity());
}

/**
 * @brief Preset 6, the scratch arenas' tightest, fits through its next
 *        shape's spawn at frame 151, where scratch A peaks.
 */
inline void test_islamicstars_peak_fits_scratch() {
  check_islamicstars_motion_fits_scratch(6, 160);
}

/**
 * @brief Every preset's whole profile run, its builds and its ripples,
 *        deforms within the effect's scratch arenas.
 */
inline void test_islamicstars_motion_fits_scratch() {
  for (size_t preset = 0; preset < std::size(IslamicBuildProbe::presets());
       ++preset)
    check_islamicstars_motion_fits_scratch(preset, PROFILE_REPLAY_FRAMES);
}

/** @brief Runs the focused IslamicStars white-box cases. */
inline void run_islamicstars_cases() {
  const auto run_case = [](auto test) {
    hs::random().seed(1337u);
    test();
  };
  run_case(test_islamicstars_seed_sprite_fade_in);
  run_case(test_islamicstars_opacity_shades_at_column_time);
  run_case(test_islamicstars_ripple_motion_is_seam_free);
  run_case(test_islamicstars_seam_copies_are_one_polygon);
  run_case(test_islamicstars_playback_motion_is_seam_free);
  run_case(test_islamicstars_build_legs_hand_their_path);
  run_case(test_islamicstars_opacity_history);
  run_case(test_islamicstars_burst_size_is_snapshotted_per_spawn);
  run_case(test_islamicstars_smooth_recipe_completion);
  run_case(test_islamicstars_default_cycle_frames);
  run_case(test_islamicstars_presets_cover_roster);
  run_case(test_islamicstars_manual_select_cuts_build);
  run_case(test_islamicstars_pause_holds_shape);

  run_case(test_islamicstars_prepared_build);
  run_case(test_islamicstars_prepared_ripples);
  run_case(test_islamicstars_endpoint_groups);
  run_case(test_islamicstars_peak_fits_scratch);
}

/**
 * @brief Module entry point for the effects white-box suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Per-effect invariants with an explicit oracle; excluded from the
 * fast-math axis.
 */
inline int run_effects_tests() {
  hs_test::ModuleFixture fixture("effects");
  const auto run_case = [](auto test) {
    hs::random().seed(1337u);
    test();
  };

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *islamic_only = std::getenv("HS_ISLAMICSTARS_ONLY");
  const char *mesh_only = std::getenv("HS_MESH_CAPTURE_ONLY");
#pragma clang diagnostic pop
  if (mesh_only && std::atoi(mesh_only) == 3) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    const char *preset = std::getenv("HS_REPLAY_PRESET");
    const char *speed = std::getenv("HS_REPLAY_TRANS_SPEED");
    const char *frames = std::getenv("HS_REPLAY_FRAMES");
#pragma clang diagnostic pop
    replay_islamicstars_profile(
        preset ? std::atoi(preset) : 6, speed ? std::atof(speed) : 4.0f,
        frames ? std::atoi(frames) : PROFILE_REPLAY_FRAMES);
    return fixture.result();
  }
  if (mesh_only && std::atoi(mesh_only) == 2) {
    run_case(test_hankinsolids_handoff_history);
    measure_mesh_capture_work();
    return fixture.result();
  }
  if (mesh_only && std::atoi(mesh_only) != 0) {
    run_case(test_hankinsolids_interlace_motion_is_seam_free);
    run_case(test_hankinsolids_seam_copies_are_one_polygon);
    run_case(test_hankinsolids_morph_motion_is_seam_free);
    run_case(test_islamicstars_ripple_motion_is_seam_free);
    run_case(test_islamicstars_seam_copies_are_one_polygon);
    run_case(test_islamicstars_playback_motion_is_seam_free);
    return fixture.result();
  }
  if (islamic_only && std::atoi(islamic_only) != 0) {
    run_islamicstars_cases();
    return fixture.result();
  }

  run_case(test_hopf_temporal_trail_oracle);
  run_case(test_dreamballs_temporal_source_oracle);
  run_case(test_dreamballs_column_work_is_bounded);
  run_case(DreamBallsWhiteBox::check_endpoint_error_bound);
  run_case(test_petalflow_temporal_source_oracle);
  run_case(test_meshfeedback_temporal_source_oracle);
  run_case(test_meshfeedback_stationary_timed_source);
  run_case(test_meshfeedback_trails_arm_cut_seam);
  run_case(test_parameter_probe_targets);
  run_case(test_float_dropdown_lint);
  run_case(test_ringspin_strobe_configuration_preserves_rendering);
  run_case(test_meshfeedback_base_mesh_selector);
  run_case(test_meshfeedback_preset_export_arity);
  run_case(test_meshfeedback_mesh_rebuild_reuses_storage);
  run_case(test_gs_palette_is_opaque);
  run_case(test_gs_seed_palettes_and_pigment_transport);
  run_case(test_gs_sparse_pigment_matches_dense);
  run_case(test_gs_direct_draw_matches_grid);
  run_case(test_gs_dot_kernel_matches_squared_distance);
  run_case(test_gs_shared_cube_projection);
  run_case(test_gs_shared_noise_palette_modifiers);
  run_case(test_gs_noise_bake_matches_reference);
  run_case(test_gs_reseed_generates_palette);
  run_case(test_gs_render_certificates_bound_lattice);
  run_case(test_gs_shared_stencil_error_is_bounded);
  run_case(test_gs_nearest_pigment_shader_matches_scalar_reference);
  run_case(test_gs_signed_coverage_and_concentration);
  run_case(test_gs_support_certificate_matches_display_geometry);
  run_case(test_gs_nearest_pigment_shader_fidelity);
  run_case(test_gs_partial_color_palette_rows);
  run_case(test_gs_reseed_invalidates_shimmer_cache);
  run_case(test_gs_dissolve_frontier_fades_before_clear);
  run_case(test_gs_substep_matches_scalar_reference);
  run_case(test_gs_inplace_frame_matches_jacobi);
  run_case(test_fishbowl_preset_and_fire_duty_cycle);
  run_case(test_fishbowl_temporal_trail_oracle);
  run_case(test_sh_preset_mode_mapping);
  run_case(test_sh_manual_preset_replaces_inflight_morph);
  run_case(test_shapeshifter_preset_defaults);
  run_case(test_shapeshifter_slider_selections_render);
  run_case(test_manual_preset_navigation);
  run_case(test_hankinsolids_manual_pause_holds_morph);
  run_case(test_hankinsolids_handoffs_shade_at_column_time);
  run_case(test_hankinsolids_handoff_history);
  run_case(test_hankinsolids_conway_crossfade_shades_at_column_time);
  run_case(test_hankinsolids_interlace_motion_is_seam_free);
  run_case(test_hankinsolids_seam_copies_are_one_polygon);
  run_case(test_hankinsolids_morph_draws_carry_previous_vertices);
  run_case(test_hankinsolids_morph_motion_is_seam_free);
  // Both tiers: arena budgets.
  run_case(test_fishbowl_scratch_estimate_covers_peak);
  run_case(test_hankinsolids_arena_budget_covers_every_solid);
  run_case(test_dreamballs_max_edge_solid_render);
  // Both tiers: Raymarch geometry and presets.
  run_case(test_raymarch_volume_random_walks_are_independent);
  run_case(test_raymarch_preset_and_placement_solids);
  run_case(test_raymarch_surface_frame_uv);
  run_case(test_raymarch_lighting_and_palette_at_column_time);
  run_case(test_raymarch_hue_noise_at_column_time);
  // Resolution-independent white-box math.
  run_case(test_gs_hot_flags_match_directed_graph);
  run_case(test_gs_pigment_codec_layout);
  run_case(test_gs_rest_state_is_fixed_point);
  run_case(test_gs_substep_signs_and_clamp);
  run_case(test_bz_q16_roundtrip);
  run_case(test_bz_advance_species_signs_and_clamp);
  run_case(test_bz_perturb_state_draw_count_pinned);
  run_case(test_hopf_projection_math);
  run_case(test_raymarch_constexpr_sqrt_converges);
  // Both tiers: additional white-box cases.
  run_case(test_needs_full_frame_gate);
  run_case(test_voronoi_axes_use_uniform_sampler);
  run_case(test_voronoi_segment_render_matches_full_frame);
  run_case(test_voronoi_instantaneous_render_pins);
  run_case(test_voronoi_endpoint_pairs_match_brute_force);
  run_case(test_voronoi_temporal_toggle);
  run_case(test_voronoi_reseed_resets_history);
  run_case(test_voronoi_retained_corners_match_search);
  run_case(test_sh_field_blend_endpoints_and_frame);
  run_case(test_sh_field_stays_inside_unit_range);
  run_case(test_sh_polarity_split_and_ao_shaping);
  run_case(test_sh_cycle_covers_every_mode);
  run_case(test_sh_instantaneous_pixel_golden);
  run_case(test_sh_columns_follow_sweep_time);
  run_case(test_sh_segmented_columns_match_single_engine);
  run_case(test_sh_manual_preset_resets_morph_history);
  run_case(test_sh_fidelity_toggle_keeps_history);
  run_case(test_sh_frame_sampler_rollover_and_sign);
  run_case(test_gs_reaction_edit_starts_dissolve);
  run_case(test_gs_stall_hold_starts_dissolve);
  run_case(test_gs_preset_pause_does_not_hold_lifecycle);
  run_case(test_gs_dissolve_absorbs_edits_and_stages_reseed);
  run_case(test_bz_authored_palette);
  run_case(test_bz_min_diffusion_step_survives_quantization);
  run_case(test_bz_perturb_state_saturates_and_nudges);
  run_case(test_bz_perturb_scales_with_timestep);
  run_case(test_bz_substep_diffuses);
  run_case(test_bz_substep_predator_wiring);
  run_case(test_bz_raster_matches_reference);
  run_case(test_bz_render_center_matches_reference);
  run_case(test_bz_untimed_frames_are_pinned);
  run_case(test_bz_slabs_shade_the_field_after_their_substep);
  run_case(test_bz_segmented_slabs_match_full_frame);
  run_case(test_bz_lattice_basis_matches_rotated_lattice);
  run_case(test_gs_untimed_frames_are_pinned);
  run_case(test_gs_lattice_basis_matches_rotated_lattice);
  run_case(test_gs_slabs_shade_the_field_after_their_substep);
  run_case(test_gs_slab_noise_keeps_the_previous_bake);
  run_case(test_gs_segmented_slabs_match_full_frame);
  run_case(test_dreamballs_preset_cycle_bookkeeping);
  CometsWhiteBox::check_paths_close();
  CometsWhiteBox::check_births_skip_earlier_columns();
  run_case(test_comets_rollover_skipped_mid_wipe);
  run_case(test_comets_wipe_shades_at_column_time);
  run_case(test_comets_wipe_history);
  ThrustersWhiteBox::check_warp_endpoints();
  ThrustersWhiteBox::check_fire_spawns_opposed_pair();
  ThrustersWhiteBox::check_collapsed_ring_falls_back_to_a_spin_axis();
  ThrustersWhiteBox::check_fifo_evicts_the_oldest_pair();
  ThrustersWhiteBox::check_expired_slots_retire_by_pair();
  RingShowerWhiteBox::check_radius_endpoints();
  DynamoWhiteBox::check_gap_reduction_uses_shortest_path();
  DynamoWhiteBox::check_overlapping_wipes_stay_in_range();
  DynamoWhiteBox::check_wipe_band_colors();
  run_case(test_dynamo_emitted_points_counts_ring_seeds);
  run_case(test_ringspin_trail_hugs_its_great_circles);
  run_case(test_hopf_trail_trim_keeps_a_segment);
  run_case(test_gnomonicstars_radius_px_covers_both_axes);
  run_case(test_gnomonicstars_spiral_cache_invalidation);
  run_case(test_gnomonicstars_stars_drift_from_last_draw);
  run_case(test_galaxies_spawn_along_spiral_with_orbital_velocity);
  run_case(test_galaxies_alternate_orbits_and_arm_spin);
  run_case(test_galaxies_arm_contrast_and_age_fade);
  run_case(test_galaxies_arm_sharpness);
  run_case(test_galaxies_star_color_mix);
  run_case(test_galaxies_core_fade);
  run_case(test_galaxies_white_star_stays_white);
  run_case(test_galaxies_emission_rate_and_capacity);
  run_case(test_galaxies_live_particle_count);
#if defined(_WIN32)
  run_case(test_column_time_scans_keep_fidelity_off_frames);
#endif
  run_case(test_column_time_scans_engage_with_fidelity_on);
  run_case(test_galaxies_black_hole_core_toggle);
  run_case(test_galaxies_packed_orbits_match_float);
  run_case(test_displacement_field_lazy_hue_table_matches_eager);
  run_case(test_displacement_field_hue_table_bind_clears_knots);
  run_case(test_displacement_field_ring_pool_abandoned_slot);
  run_case(test_displacement_field_ring_pool_capacity_release);
  run_case(test_displacement_field_footprint_covers_init);
  run_case(test_displacement_field_pause_holds_phase_only);
  run_case(test_displacement_field_default_cycle_frames);
  run_case(test_displacement_field_presets);
  run_case(test_displacement_field_edit_cancels_departure);
  run_case(test_displacement_field_startup_writes_preserve_entry);
  run_case(test_displacement_field_paused_balls_deform_render);
  run_case(test_displacement_field_paused_balls_keep_moving);
  run_case(test_displacement_field_knot_walk_unit_near_ball);
  run_case(test_displacement_field_zero_hue_scale_is_exact);
  run_case(test_displacement_field_octave_bake_tracks_noise);
  run_case(test_displacement_field_arm_boundaries_meet_endpoints);
  run_case(test_displacement_field_single_pass_reference);
  run_case(test_displacement_field_boundary_stencil);
  run_case(test_displacement_field_turn_follows_large_steps);
  run_case(test_displacement_field_alternate_groups_match_windows);
  run_case(test_glitch_lens_unit_norm);
  run_case(test_mobius_rings_conformal_and_counter_rotation);
  run_islamicstars_cases();

  // FULL tier only (HS_EFFECTS_FULL=1).
  if (effects_full_suite()) {
    run_case(test_voronoi_union_candidates_cover_nearest);
    run_case(test_voronoi_column_time_matches_brute_force);
    run_case(test_sh_pullback_matches_direct_colorize);
    run_case(test_sh_morph_chain_rearms);
    run_case(test_gs_evolution_stays_bounded);
    run_case(test_gs_reaction_corner_stays_bounded);
    run_case(test_gs_dissolve_clears_and_reseeds);
    run_case(test_gs_staged_reseed_matches_synchronous);
    run_case(test_dreamballs_base_mesh_selector);
    run_case(test_dreamballs_weave_topology);
    run_case(test_dreamballs_defect_weave_renders);
    run_case(test_dreamballs_respawn_fires_and_honors_pause);
    run_case(test_dreamballs_select_mid_fade_in_finishes_fade);
    run_case(test_dreamballs_select_mid_fade_out_respawns);
    run_case(test_meshfeedback_flush_precedes_mesh_draw);
    run_case(test_meshfeedback_preset_rotation_syncs_noise);
    run_case(test_comets_manual_preset_restarts_path);
    run_case(test_ash_cloud_value_cutout_gates_the_frame);
    run_case(test_dynamo_trail_ceiling_bounds_the_ring);
    run_case(test_raymarch_unit_bounds_contains_twisted_tube);
    run_case(test_petalflow_spawn_gap_bounded);
    run_case(test_petalflow_last_free_slot);
    run_case(test_displacement_field_hue_table_fidelity);
    run_case(test_displacement_field_hue_table_frame_fidelity);
    run_case(test_displacement_field_clip_tiles_full);
    run_case(test_displacement_field_ball_spans_and_lifecycle);
    run_case(test_islamicstars_recipe_build_smoke);
    run_case(test_islamicstars_roster_cycle_fits_budget);
    run_case(test_islamicstars_dual_bridge_fits_budget);
    run_case(test_islamicstars_motion_fits_scratch);
  } else {
    std::printf("  [TIER] full effects cases omitted; set HS_EFFECTS_FULL=1\n");
  }

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test
