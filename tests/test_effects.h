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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string_view>
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
 * @brief Small-aspect render width in pixels.
 * @details The holosphere <96,20> resolution, run here under asserts.
 */
constexpr int SMALL_W = 96;
/**
 * @brief Small-aspect render height in pixels.
 * @details Paired with SMALL_W for the <96,20> specialization.
 */
constexpr int SMALL_H = 20;

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
 * @brief Resets the process-global effect state to a clean per-effect baseline.
 * @details Forwards to hs_test::reset_globals(), which also releases the mock
 *          clock.
 */
inline void reset_effect_globals() { hs_test::reset_globals(); }

/**
 * @brief Params of a choreographed effect's preset @p index.
 * @tparam E The effect type.
 * @param index Preset index.
 * @return The effect's own preset entry, or its initial params when it
 *         declares no preset(index).
 * @details ChoreographedEffect falls back to initial params for a
 * single-preset effect, and its own resolver is private, so a test reproduces
 * the fallback.
 */
template <typename E>
typename E::Params preset_params_or_initial(size_t index) {
  if constexpr (requires { E::preset(index).params; })
    return E::preset(index).params;
  else
    return E::initial_params();
}

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
    uint64_t hash = hs_test::FNV1A64_BASIS;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const Pixel &pixel = effect.get_pixel(x, y);
        for (uint16_t channel : {pixel.r, pixel.g, pixel.b})
          hash = hs_test::fnv1a64_channel(hash, channel);
      }
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
    if (def.animated || def.readonly)
      continue;
    const float range = def.max - def.min;
    if (range <= 0.0f)
      continue;
    const float cur = def.get();
    // An in-range target well clear of the current value, so a revert is
    // visible; a bool has only its flipped value.
    float target = parameter_probe_target(def, cur);
    // Round the probe to a value an integer target can hold.
    if (def.is_integer()) {
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
    // pre-write `cur`, catching a slow per-frame revert. A
    // bool must read back exactly.
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
    for (size_t restore = 0; restore < count; ++restore)
      effect.updateParameter(names[restore], original[restore]);
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
  uint64_t fold = hs_test::FNV1A64_BASIS;

  E<W, H> effect;
  effect.init();
  for (int f = 0; f < frames; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    if (frame_fold || lit)
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel p = effect.get_pixel(x, y);
          if (lit && (p.r || p.g || p.b))
            *lit = true;
          for (uint16_t channel : {p.r, p.g, p.b})
            fold = hs_test::fnv1a64_channel(fold, channel);
        }
  }
  if (frame_fold)
    *frame_fold = fold;

  capture_frame<W, H>(effect, out);
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
 * @details The clip clear leaves the region outside the display band holding
 *          pixels from the frame that last wrote this buffer. Anything an effect
 *          reads back or carries forward from there shows up here as a mismatch
 *          against the whole-buffer clear.
 */
template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void clip_clear_parity_one(const char *name) {
  constexpr int S = H * 2;
  const int frames = PARITY_FRAMES;

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

    size_t different = 0;
    for (size_t i = 0; i < full.size(); ++i) {
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
  // Two all-black renders agree trivially; require output across the sweep.
  if (lit == 0)
    std::printf("  CLIP-CLEAR DARK %-20s no lit pixel over %zu displayed "
                "pixels\n",
                name, displayed_pixels);
  HS_EXPECT(lit > 0, "clip-clear parity must compare a lit render");
  hs::clear_mock_time();
}

inline constexpr std::array<int, 24> SH_PRESET_MODES{
    6,  1,  2,  3,  4,  5,  7,  8,  9,  10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

/**
 * @brief White-box accessor for SphericalHarmonics' morph chain and shader inputs.
 * @details Befriended in effects/SphericalHarmonics.h.
 */
struct SphericalHarmonicsWhiteBox {
  using SH = SphericalHarmonics<SMALL_W, SMALL_H>;
  using Field = SH::HarmonicField;
  using Pipeline = SH::RenderPipeline;
  using PipelineFrame = Pipeline::Frame;

  static int current_idx(const SH &fx) { return fx.current_idx; }
  static int next_idx(const SH &fx) { return fx.next_idx; }
  static float morph_alpha(const SH &fx) { return fx.morph_alpha; }
  static math::Quaternion orientation(const SH &fx) {
    return fx.orientation.get();
  }
  static const BakedPalette &palette(const SH &fx) { return fx.baked_palette; }
  static float amplitude(const SH &fx) { return fx.params.amplitude; }
  static constexpr int max_mode_idx() { return SH::MAX_MODE_IDX; }

  static PipelineFrame prepare_pipeline(const SH &fx, int l1, int m1, int l2,
                                        int m2, float blend,
                                        const math::Quaternion &orientation,
                                        float amplitude) {
    return Pipeline::prepare({{l1, m1, l2, m2, blend, orientation},
                              {&fx.baked_palette.view(), amplitude}});
  }

  static Color4 shade_pipeline(const math::Vector &view,
                               const PipelineFrame &frame) {
    return Pipeline::evaluate(view, frame.ctx, frame.prepared);
  }

  static Color4 shade_legacy(const SH &fx, float value, float amplitude) {
    return SH::colorize_harmonic(value, {&fx.baked_palette.view(), amplitude});
  }

  // Pinning both morph endpoints on one mode makes the blend an identity, so a
  // frame renders one pure harmonic whatever morph_alpha has reached.
  static void pin_mode(SH &fx, int idx) {
    fx.current_idx = idx;
    fx.next_idx = idx;
  }
  static void set_next_idx(SH &fx, int idx) { fx.next_idx = idx; }
  static void step_timeline(SH &fx, Canvas &canvas) {
    fx.timeline.step(canvas);
  }
};

/**
 * @brief Pins HarmonicField's blend endpoints, polarity and local frame.
 */
inline void test_sh_field_write_through_and_endpoints() {
  using Field = SphericalHarmonicsWhiteBox::Field;
  const math::Quaternion spin =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 1.1f);
  const math::Quaternion identity;

  constexpr int LA = 1, MA = 0, LB = 3, MB = 2;
  Field mix_start(LA, MA, LB, MB, 0.0f, spin);
  Field mix_end(LA, MA, LB, MB, 1.0f, spin);
  Field pure_a(LA, MA, LA, MA, 0.0f, spin);
  Field pure_b(LB, MB, LB, MB, 0.0f, spin);
  Field pure_a_unrotated(LA, MA, LA, MA, 0.0f, identity);

  constexpr int PHI_STEPS = 24, THETA_STEPS = 32;
  int positives = 0, negatives = 0;
  double worst_start = 0.0, worst_end = 0.0, worst_frame = 0.0;
  for (int i = 0; i <= PHI_STEPS; ++i) {
    const float phi = math::PI_F * i / PHI_STEPS;
    for (int j = 0; j < THETA_STEPS; ++j) {
      const float theta = 2.0f * math::PI_F * j / THETA_STEPS;
      const math::Vector p(sinf(phi) * cosf(theta), cosf(phi),
                           sinf(phi) * sinf(theta));
      const float start = mix_start.sample(p);
      const float end = mix_end.sample(p);
      const float a = pure_a.sample(p);
      const float b = pure_b.sample(p);
      const float unspun = pure_a_unrotated.sample(p);
      // Rotating the sample by the same quaternion the shape carries must land
      // back on the unrotated shape's value.
      const float spun = pure_a.sample(math::rotate(p, spin));

      if (a > 0.02f)
        ++positives;
      if (a < -0.02f)
        ++negatives;
      worst_start =
          hs_test::fold_worst<double>(worst_start, std::fabs(start - a));
      worst_end = hs_test::fold_worst<double>(worst_end, std::fabs(end - b));
      worst_frame =
          hs_test::fold_worst<double>(worst_frame, std::fabs(spun - unspun));
    }
  }

  // A dipole must actually change sign, or the polarity split below is vacuous.
  HS_EXPECT_GT(positives, 0);
  HS_EXPECT_GT(negatives, 0);

  if (worst_start > 0.0 || worst_end > 1e-5 || worst_frame > 1e-5)
    std::printf("  SH field: blend0 err=%g blend1 err=%g frame err=%g\n",
                worst_start, worst_end, worst_frame);
  // Blend 0 is exact; blend 1 reaches the second mode within lerp rounding.
  HS_EXPECT_EQ(worst_start, 0.0);
  HS_EXPECT_LT(worst_end, 1e-5);
  HS_EXPECT_LT(worst_frame, 1e-5);
}

/** @brief Keeps every shipped harmonic morph inside SampleSphere's signed range. */
inline void test_sh_field_stays_inside_unit_range() {
  using WB = SphericalHarmonicsWhiteBox;
  const math::Quaternion orientation =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 1.1f);
  constexpr std::array<float, 5> BLENDS{0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  constexpr int PHI_STEPS = 48;
  constexpr int THETA_STEPS = 64;
  float peak = 0.0f;

  for (int first = 1; first <= WB::max_mode_idx(); ++first) {
    const int second = first % WB::max_mode_idx() + 1;
    auto [l1, m1] = SHMath::decode_lm(first);
    auto [l2, m2] = SHMath::decode_lm(second);
    for (float blend : BLENDS) {
      const WB::Field field(l1, m1, l2, m2, blend, orientation);
      for (int i = 0; i <= PHI_STEPS; ++i) {
        const float phi = math::PI_F * i / PHI_STEPS;
        for (int j = 0; j < THETA_STEPS; ++j) {
          const float theta = math::TWO_PI_F * j / THETA_STEPS;
          const math::Vector point(sinf(phi) * cosf(theta), cosf(phi),
                                   sinf(phi) * sinf(theta));
          peak = hs_test::fold_worst(peak, std::fabs(field.sample(point)));
        }
      }
    }
  }

  std::printf("  SH field bound: peak %.6f\n", static_cast<double>(peak));
  HS_EXPECT_LT(peak, 1.0f);
}

/** @brief Bounds color drift from the signed-to-unit carrier round trip. */
inline void test_sh_pullback_matches_legacy_shader() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  WB::SH fx;
  fx.init();

  const math::Quaternion orientation =
      math::make_rotation(math::Vector(-0.4f, 0.2f, 0.9f).normalized(), 0.83f);
  constexpr std::array<float, 3> BLENDS{0.0f, 0.5f, 1.0f};
  constexpr std::array<float, 3> AMPLITUDES{0.2f, 3.2f, 7.0f};
  constexpr int PHI_STEPS = 24;
  constexpr int THETA_STEPS = 32;
  int max_channel_delta = 0;
  int differing_pixels = 0;

  for (int first = 1; first <= WB::max_mode_idx(); ++first) {
    const int second = first % WB::max_mode_idx() + 1;
    auto [l1, m1] = SHMath::decode_lm(first);
    auto [l2, m2] = SHMath::decode_lm(second);
    for (float blend : BLENDS) {
      const WB::Field legacy_field(l1, m1, l2, m2, blend, orientation);
      for (float amplitude : AMPLITUDES) {
        const WB::PipelineFrame frame = WB::prepare_pipeline(
            fx, l1, m1, l2, m2, blend, orientation, amplitude);
        for (int i = 0; i <= PHI_STEPS; ++i) {
          const float phi = math::PI_F * i / PHI_STEPS;
          for (int j = 0; j < THETA_STEPS; ++j) {
            const float theta = math::TWO_PI_F * j / THETA_STEPS;
            const math::Vector point(sinf(phi) * cosf(theta), cosf(phi),
                                     sinf(phi) * sinf(theta));
            const Color4 legacy =
                WB::shade_legacy(fx, legacy_field.sample(point), amplitude);
            const Color4 pullback = WB::shade_pipeline(point, frame);
            const int red = std::abs(static_cast<int>(legacy.color.r) -
                                     static_cast<int>(pullback.color.r));
            const int green = std::abs(static_cast<int>(legacy.color.g) -
                                       static_cast<int>(pullback.color.g));
            const int blue = std::abs(static_cast<int>(legacy.color.b) -
                                      static_cast<int>(pullback.color.b));
            const int delta = std::max({red, green, blue});
            max_channel_delta = std::max(max_channel_delta, delta);
            differing_pixels += static_cast<int>(delta != 0);
            HS_EXPECT_EQ(legacy.alpha, pullback.alpha);
          }
        }
      }
    }
  }

  std::printf("  SH pullback parity: max channel delta %d, differing %d\n",
              max_channel_delta, differing_pixels);
  HS_EXPECT_LE(max_channel_delta, 1);
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
  HS_EXPECT_EQ(fx.updateParameter("Amplitude", amplitude),
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
 * @details Classifies each rendered pixel by the field's sign and magnitude: a
 *          saturated positive lobe wears the palette's top color unmodified, the
 *          negative lobe wears it with red and blue swapped and green dimmed,
 *          and below the AO falloff every pixel is dimmed under its palette
 *          entry.
 */
inline void test_sh_polarity_split_and_ao_shaping() {
  using WB = SphericalHarmonicsWhiteBox;
  // idx 2 is (l=1, m=0): one nodal circle, so both polarities cover a wide band.
  constexpr int DIPOLE_IDX = 2;

  // Amplitude 10 (the slider maximum) saturates the palette across most of both
  // lobes.
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
            const float val =
                field.sample(math::pixel_to_vector<SMALL_W, SMALL_H>(x, y));
            const Pixel &px = fx.get_pixel(x, y);
            const bool positive = val > 0.0f;
            // A dipole's two lobes are antipodal mirrors, so their magnitude
            // distributions match and their green totals differ only by the
            // negative recolor's scale.
            (positive ? pos_green : neg_green) += px.g;

            const float mag = std::fabs(val) * amp;
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
            const float val =
                field.sample(math::pixel_to_vector<SMALL_W, SMALL_H>(x, y));
            const float mag = std::fabs(val) * amp;
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
 *          commits the target and calls start_morph() again; three commits mean
 *          the callback re-armed twice.
 */
inline void test_sh_morph_chain_rearms() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SH fx;
  fx.init();

  const int seed = WB::current_idx(fx);
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
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
      const size_t expected_preset = static_cast<size_t>(
          std::find(SH_PRESET_MODES.begin(), SH_PRESET_MODES.end(), now) -
          SH_PRESET_MODES.begin());
      HS_EXPECT_EQ(fx.getPresetIndex(), expected_preset);
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

/** @brief Maps every runtime preset to its stable harmonic mode. */
inline void test_sh_preset_mode_mapping() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  WB::SH fx;
  fx.init();

  HS_EXPECT_EQ(fx.getParameters().size(), 1u);
  HS_EXPECT_TRUE(fx.getParameters().find("Amplitude") != nullptr);
  HS_EXPECT_TRUE(fx.getParameters().find("Debug BB") == nullptr);
  HS_EXPECT_EQ(fx.updateParameter("Debug BB", 1.0f),
               ParamSetResult::UNKNOWN_PARAM);

  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::current_idx(fx), SH_PRESET_MODES[0]);
  for (size_t preset = 0; preset < SH_PRESET_MODES.size(); ++preset) {
    HS_EXPECT_TRUE(fx.selectPreset(preset));
    HS_EXPECT_EQ(fx.getPresetIndex(), preset);
    HS_EXPECT_EQ(WB::current_idx(fx), SH_PRESET_MODES[preset]);
    HS_EXPECT_EQ(WB::next_idx(fx), SH_PRESET_MODES[preset]);
    HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);
  }
}

/** @brief Keeps a manual harmonic selection past the automatic leg it replaces. */
inline void test_sh_manual_preset_replaces_inflight_morph() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SH fx;
  fx.init();

  HS_EXPECT_EQ(fx.getPresetCount(), 24u);
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::current_idx(fx), 6);
  for (int frame = 0; frame < 8; ++frame) {
    hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                      static_cast<unsigned long>(frame) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();
  }
  HS_EXPECT_GT(WB::morph_alpha(fx), 0.0f);
  const int replaced_target = WB::next_idx(fx);

  const size_t selected_preset = replaced_target == 24 ? 22u : 23u;
  const int selected_mode = SH_PRESET_MODES[selected_preset];
  HS_EXPECT_TRUE(fx.selectPreset(selected_preset));
  HS_EXPECT_TRUE(fx.animations_paused());
  HS_EXPECT_EQ(fx.getPresetIndex(), selected_preset);
  HS_EXPECT_EQ(WB::current_idx(fx), selected_mode);
  HS_EXPECT_EQ(WB::next_idx(fx), selected_mode);
  HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);

  fx.draw_frame();
  fx.advance_display();
  HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);

  fx.setAnimationsPaused(false);
  for (int frame = 0; frame < 64; ++frame) {
    hs::set_mock_time(static_cast<unsigned long>(frame + 9) * FRAME_MS,
                      static_cast<unsigned long>(frame + 9) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();
  }
  hs::clear_mock_time();

  HS_EXPECT_EQ(WB::current_idx(fx), selected_mode);
  HS_EXPECT_TRUE(WB::current_idx(fx) != replaced_target);
  HS_EXPECT_EQ(fx.getPresetIndex(), selected_preset);
  HS_EXPECT_TRUE(WB::next_idx(fx) != selected_mode);
}

#include "tests/effects/reaction_diffusion_gs.h"
#include "tests/effects/reaction_diffusion_bz.h"
#include "tests/effects/hankin.h"
#include "tests/effects/mesh_feedback.h"
#include "tests/effects/numeric_invariants.h"
#include "tests/effects/geometry.h"
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

  run_case(test_parameter_probe_targets);
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
  run_case(test_gs_dissolve_frontier_fades_before_clear);
  run_case(test_gs_substep_matches_scalar_reference);
  run_case(test_gs_inplace_frame_matches_jacobi);
  run_case(test_fishbowl_preset_and_fire_duty_cycle);
  run_case(test_sh_preset_mode_mapping);
  run_case(test_sh_manual_preset_replaces_inflight_morph);
  run_case(test_shapeshifter_preset_defaults);
  run_case(test_shapeshifter_slider_selections_render);
  run_case(test_manual_preset_navigation);
  run_case(test_hankinsolids_manual_pause_holds_morph);
  // Both tiers: arena budgets.
  run_case(test_fishbowl_scratch_estimate_covers_peak);
  run_case(test_hankinsolids_arena_budget_covers_every_solid);
  run_case(test_dreamballs_max_edge_solid_render);
  // Both tiers: Raymarch geometry and presets.
  run_case(test_raymarch_volume_random_walks_are_independent);
  run_case(test_raymarch_preset_and_placement_solids);
  run_case(test_raymarch_surface_frame_uv);
  // Resolution-independent white-box math.
  run_case(test_gs_q16_roundtrip);
  run_case(test_gs_hot_flags_match_directed_graph);
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
  run_case(test_sh_field_write_through_and_endpoints);
  run_case(test_sh_field_stays_inside_unit_range);
  run_case(test_sh_polarity_split_and_ao_shaping);
  run_case(test_gs_reaction_edit_starts_dissolve);
  run_case(test_bz_legacy_palette);
  run_case(test_bz_min_diffusion_step_survives_quantization);
  run_case(test_bz_perturb_state_saturates_and_nudges);
  run_case(test_bz_perturb_scales_with_timestep);
  run_case(test_bz_substep_diffuses);
  run_case(test_bz_raster_matches_reference);
  run_case(test_bz_render_center_matches_reference);
  run_case(test_dreamballs_preset_cycle_bookkeeping);
  CometsWhiteBox::check_paths_close();
  run_case(test_comets_rollover_skipped_mid_wipe);
  ThrustersWhiteBox::check_warp_endpoints();
  ThrustersWhiteBox::check_fire_spawns_opposed_pair();
  ThrustersWhiteBox::check_collapsed_ring_falls_back_to_a_spin_axis();
  ThrustersWhiteBox::check_fifo_evicts_the_oldest_pair();
  ThrustersWhiteBox::check_expired_slots_retire_by_pair();
  RingShowerWhiteBox::check_radius_endpoints();
  DynamoWhiteBox::check_overlapping_wipes_stay_in_range();
  run_case(test_dynamo_emitted_points_counts_ring_seeds);
  run_case(test_ringspin_trail_hugs_its_great_circles);
  run_case(test_hopf_trail_trim_keeps_a_segment);
  run_case(test_gnomonicstars_radius_px_covers_both_axes);
  run_case(test_gnomonicstars_spiral_cache_invalidation);
  run_case(test_displacement_field_lazy_hue_table_matches_eager);
  run_case(test_displacement_field_zero_hue_scale_is_exact);
  run_case(test_displacement_field_octave_bake_tracks_noise);
  run_case(test_glitch_lens_unit_norm);
  run_case(test_mobius_rings_conformal_and_counter_rotation);
  run_case(test_islamicstars_seed_sprite_fade_in);
  run_case(test_islamicstars_burst_size_is_snapshotted_per_spawn);
  run_case(test_islamicstars_smooth_recipe_completion);

  // FULL tier only (HS_EFFECTS_FULL=1).
  if (effects_full_suite()) {
    run_case(test_voronoi_union_candidates_cover_nearest);
    run_case(test_sh_pullback_matches_legacy_shader);
    run_case(test_sh_morph_chain_rearms);
    run_case(test_gs_evolution_stays_bounded);
    run_case(test_gs_reaction_corner_stays_bounded);
    run_case(test_gs_dissolve_clears_and_reseeds);
    run_case(test_gs_staged_reseed_matches_synchronous);
    run_case(test_dreamballs_base_mesh_selector);
    run_case(test_dreamballs_weave_topology);
    run_case(test_dreamballs_respawn_fires_and_honors_pause);
    run_case(test_meshfeedback_flush_precedes_mesh_draw);
    run_case(test_meshfeedback_preset_rotation_syncs_noise);
    run_case(test_comets_manual_preset_restarts_path);
    run_case(test_ash_cloud_value_cutout_gates_the_frame);
    run_case(test_dynamo_trail_ceiling_bounds_the_ring);
    run_case(test_raymarch_unit_bounds_contains_twisted_tube);
    run_case(test_petalflow_spawn_gap_bounded);
    run_case(test_displacement_field_hue_table_fidelity);
    run_case(test_displacement_field_hue_table_frame_fidelity);
    run_case(test_displacement_field_clip_tiles_full);
    run_case(test_displacement_field_ball_spans_and_lifecycle);
    run_case(test_islamicstars_recipe_build_smoke);
    run_case(test_islamicstars_roster_cycle_fits_budget);
    run_case(test_islamicstars_dual_bridge_fits_budget);
  } else {
    std::printf("  [TIER] full effects cases omitted; set HS_EFFECTS_FULL=1\n");
  }

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test
