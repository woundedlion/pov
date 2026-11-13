/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared effect rendering and oracle fixtures.
 */
#pragma once

#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "core/render/canvas.h"
#include "core/render/sdf/volume.h"
#include "core/memory.h"
#include "hardware/pov_segment_map.h"
#if HS_ENABLE_CHAIN_INTERPRETER
#include "workbench/shader/chain_host.h"
#endif
#include "tests/mesh_test_util.h"
#include "tests/fd_capture_util.h"
#include "tests/pixel_test_util.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/temporal_oracle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace effects_tests {

/** @brief Runs an effect case from the fixture random stream.
 * @param test Case entry point. */
inline void run_case(auto test) {
  reset_globals();
  test();
}

/** @brief Authored SphericalHarmonics preset mode indices. */
inline constexpr std::array<int, 24> SH_PRESET_MODES{
    6,  1,  2,  3,  4,  5,  7,  8,  9,  10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

/**
 * @brief Primary production render width in pixels.
 */
constexpr int DEFAULT_W = 288;
/**
 * @brief Primary production render height in pixels.
 * @details Paired with DEFAULT_W for the full-sphere production resolution.
 */
constexpr int DEFAULT_H = 144;

inline void lint_dead_sliders(Effect &effect, const char *name);

/** @brief Sets companion controls that make the probed parameter admissible. */
template <typename Effect>
inline void prepare_chain_parameter_probe(Effect &effect,
                                          std::string_view parameter) {
  const size_t separator = parameter.rfind('.');
  const std::string prefix(parameter.substr(0, separator + 1));
  const std::string_view field = parameter.substr(separator + 1);
  const auto set = [&](const char *companion, float value) {
    const std::string qualified = prefix + companion;
    HS_EXPECT_EQ(effect.updateParameter(qualified.c_str(), value),
                 ParamSetResult::APPLIED);
  };
  if (field == "hue-shift-mode")
    set("modulation-enabled", 0.0f);
  if (field == "pattern-scale" || field == "modulation-source" ||
      field == "modulation-evaluation") {
    set("modulation-mechanic", 0.0f);
    set("pattern-noise-basis", 0.0f);
    set("pattern-contrast", 0.0f);
    set("pattern-noise-variant", 0.0f);
    if (field != "modulation-evaluation")
      set("modulation-evaluation", 0.0f);
  }
  if (field == "modulation-mechanic")
    set("modulation-evaluation", 1.0f);
  if (field == "pattern-contrast" || field == "pattern-noise-basis" ||
      field == "pattern-noise-variant") {
    set("modulation-evaluation", 1.0f);
    if (field != "pattern-noise-variant")
      set("pattern-noise-variant", 1.0f);
  }
}

/**
 * @brief Checks that animated parameter edits pause automatic preset advancement.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in diagnostic output.
 */
template <typename FX>
inline void lint_animated_pause(FX &effect, const char *name);

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
  reset_globals();
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
 * @tparam FX Effect type, including its coordinated parameter API when present.
 * @param effect Effect instance whose automated params are probed.
 * @param name Effect name used in PAUSE LEAK diagnostic output.
 * @details Parameters are audited one at a time over PAUSE_AUDIT_FRAMES in
 *          total.
 */
template <typename FX>
inline void lint_animated_pause(FX &effect, const char *name) {
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

  // Modulation recipes restore whole, since their controls gate one another.
  std::vector<std::function<ParamSetResult()>> restore_recipes;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  if constexpr (requires {
                  effect.visit_document_modulations(
                      [](std::string_view, const auto &) {});
                }) {
    effect.visit_document_modulations([&](std::string_view resource,
                                          const auto &saved) {
      restore_recipes.push_back([&effect, resource = std::string(resource),
                                 saved] {
        return effect.edit_document_modulation(resource, [&](auto &candidate) {
          candidate = saved;
          return ParamSetResult::APPLIED;
        });
      });
    });
  }
#endif

  constexpr int PAUSE_AUDIT_FRAMES = 500;
  const int frames_per_param =
      std::max(4, static_cast<int>((PAUSE_AUDIT_FRAMES + count - 1) / count));
  bool leaked = false;
  for (size_t index = 0; index < count; ++index) {
    HS_CONTEXT(names[index]);
    bool batched = false;
#if HS_ENABLE_CHAIN_INTERPRETER
    if constexpr (requires {
                    effect.update_parameters(
                        std::span<const ShaderChainParameterWrite>{});
                  }) {
      std::vector<ShaderChainParameterWrite> restore;
      for (size_t item = 0; item < count; ++item)
        restore.push_back({names[item], original[item]});
      HS_EXPECT_EQ(effect.update_parameters(restore), ParamSetResult::APPLIED);
      batched = true;
    }
#endif
    if (!batched) {
      for (const auto &restore : restore_recipes)
        HS_EXPECT_EQ(restore(), ParamSetResult::APPLIED);
      for (size_t restore = 0; restore < count; ++restore)
        HS_EXPECT_EQ(effect.updateParameter(names[restore], original[restore]),
                     ParamSetResult::APPLIED);
    }
    if (batched) {
      prepare_chain_parameter_probe(effect, names[index]);
    }
    effect.setAnimationsPaused(false);
    // A recipe constraint can refuse the probe; halve it toward the original.
    float probe = target[index];
    ParamSetResult written = effect.updateParameter(names[index], probe);
    for (int halving = 0;
         !batched && written == ParamSetResult::INADMISSIBLE && halving < 8;
         ++halving) {
      probe = original[index] + 0.5f * (probe - original[index]);
      written = effect.updateParameter(names[index], probe);
    }
    HS_EXPECT_EQ(written, ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(effect.animations_paused());
    const auto *written_def = effect.getParameters().find(names[index]);
    const float written_value = written_def->get_requested();
    const float eps =
        written_def->is_bool()
            ? 0.0f
            : fmaxf(1e-3f, 1e-3f * (written_def->max - written_def->min));
    effect.draw_frame();
    effect.advance_display();
    const auto *after_write = effect.getParameters().find(names[index]);
    const bool persisted =
        fabsf(after_write->get_requested() - written_value) <= eps;
    if (!persisted)
      std::printf(
          "  PAUSED WRITE LOST %s::%s wrote %.4f, engine reverted to %.4f\n",
          name, names[index], static_cast<double>(written_value),
          static_cast<double>(after_write->get_requested()));
    HS_EXPECT(persisted, "paused animated parameter writes must persist");
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
 * @details Calls reset_globals() and pins the mock clock to the frame cadence,
 * so two calls start from an identical state.
 */
template <template <int, int> class E, int W = DEFAULT_W, int H = DEFAULT_H>
inline void render_capture(std::vector<Pixel> &out, int frames,
                           uint64_t *frame_fold = nullptr,
                           bool *lit = nullptr) {
  reset_globals();
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
 * @brief Records every displayed frame at the selected temporal sampling mode.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param frames Frames to render.
 * @param clip Display band set before the first frame, or null for the full
 *        canvas.
 * @param temporal Whether to sample columns at their emission times.
 * @return The recorded frames.
 */
template <template <int, int> class E, int W, int H>
inline RenderRun pinned_clock_run(int frames, const ClipWindow *clip = nullptr,
                                  bool temporal = false) {
  reset_globals();
  pin_frame_clock(0);
  Spacetime::set_flat_surface(!temporal);
  RenderRun run;
  {
    E<W, H> effect;
    effect.init();
    if (clip)
      effect.set_clip(clip->y0, clip->y1, clip->x0, clip->x1);
    for (int f = 0; f < frames; ++f) {
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      record_frame<W, H>(run, effect);
    }
  }
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return run;
}

/**
 * @brief Folds every displayed frame at the selected temporal sampling mode.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @tparam W Render width in pixels.
 * @tparam H Render height in pixels.
 * @param frames Frames to render.
 * @param clip Display band set before the first frame, or null for the full
 *        canvas.
 * @param temporal Whether to sample columns at their emission times.
 * @return Hash, RGB totals and brightness moments of the sequence.
 */
template <template <int, int> class E, int W, int H>
inline RenderFold pinned_clock_frame_fold(int frames,
                                          const ClipWindow *clip = nullptr,
                                          bool temporal = false) {
  return pinned_clock_run<E, W, H>(frames, clip, temporal).fold;
}

/**
 * @brief Checks a reaction-diffusion effect's renders with temporal sampling
 *        disabled: lit at the test and shipped sizes, and a banded clip
 *        equal to the full render inside the band.
 * @tparam E Effect class template, instantiated as E<W, H>.
 * @param name Effect label without spaces.
 */
template <template <int, int> class E>
inline void expect_untimed_rd_render_properties(const char *name) {
  constexpr ClipWindow BAND{3, 17, 39, 83};
  constexpr int FRAMES = 16, SHIPPED_FRAMES = 3;
  HS_CONTEXT(name);
  const RenderRun full = pinned_clock_run<E, SMALL_W, SMALL_H>(FRAMES);
  expect_lit("small", full);
  const RenderRun clipped =
      pinned_clock_run<E, SMALL_W, SMALL_H>(FRAMES, &BAND);
  expect_clip_matches_full<SMALL_W, SMALL_H>(name, full, clipped,
                                             [&](int) { return BAND; });
  const RenderRun shipped =
      pinned_clock_run<E, DEFAULT_W, DEFAULT_H>(SHIPPED_FRAMES);
  expect_lit("shipped", shipped);
  char label[32];
  std::snprintf(label, sizeof label, "%s-small", name);
  print_render_hash(label, "full", full.fold.hash);
  print_render_hash(label, "clipped", clipped.fold.hash);
  std::snprintf(label, sizeof label, "%s-shipped", name);
  print_render_hash(label, "full", shipped.fold.hash);
}

/**
 * @brief Scrambles every output-affecting global that render_capture() resets.
 * @details Runs between captures to exercise recovery from dirty process state.
 */
inline void perturb_determinism_globals() {
  Spacetime::set_flat_surface(true);
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
    reset_globals();
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

#include "tests/effects/mesh_seam_oracle.h"

/**
 * @brief Slab bands match the full frame after a simulation warm-up.
 * @tparam FX Reaction-diffusion effect type.
 * @param warmup Frames to advance before capturing comparisons.
 */
template <typename FX>
inline void expect_slab_bands_match_full_frame(int warmup) {
  constexpr int W = SMALL_W, H = SMALL_H, FRAMES = 3;
  auto render = [warmup](int x0, int x1) {
    reset_globals();
    pin_frame_clock(0);
    FX fx;
    fx.init();
    fx.set_clip(0, H, x0, x1);
    std::vector<Pixel> out;
    for (int f = 0; f < warmup + FRAMES; ++f) {
      pin_frame_clock(f);
      fx.draw_frame();
      fx.advance_display();
      if (f < warmup)
        continue;
      for (int y = 0; y < H; ++y)
        for (int x = x0; x < x1; ++x)
          out.push_back(fx.get_pixel(x, y));
    }
    hs::clear_mock_time();
    return out;
  };
  const std::vector<Pixel> full = render(0, W);
  size_t lit = 0;
  for (const Pixel &p : full)
    lit += !is_black(p);
  HS_EXPECT_GT(lit, static_cast<size_t>(0));
  const std::array<std::array<int, 2>, 4> BANDS{
      {{0, 30}, {20, 60}, {44, 52}, {70, W}}};
  for (const auto &[x0, x1] : BANDS) {
    HS_CONTEXT("band", x0, x1);
    const std::vector<Pixel> band = render(x0, x1);
    size_t next = 0;
    int mismatches = 0;
    for (int f = 0; f < FRAMES; ++f)
      for (int y = 0; y < H; ++y)
        for (int x = x0; x < x1; ++x)
          mismatches +=
              band[next++] != full[(static_cast<size_t>(f) * H + y) * W + x];
    HS_EXPECT_EQ(mismatches, 0);
  }
}

} // namespace effects_tests
} // namespace hs_test
