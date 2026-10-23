/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Pixel parity of every filter-pipeline effect against recorded frame hashes.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include "core/render/spacetime.h"
#include "hardware/pov_segment_map.h"
#include "targets/effects.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace plot_parity_tests {

/** @brief Render width of every parity scenario. */
constexpr int PARITY_W = 288;
/** @brief Render height of every parity scenario. */
constexpr int PARITY_H = 144;
/** @brief Frames folded into each scenario's hash. */
constexpr int PARITY_FRAMES = 48;
/** @brief Segments the clipped scenario's display is divided into. */
constexpr int PARITY_SEGMENTS = 4;
/** @brief Segment the clipped scenario renders. */
constexpr int PARITY_SEGMENT_ID = 1;

/** @brief Frame-hash folds of one effect across the parity scenarios. */
struct ParityHashes {
  const char *effect;       /**< Effect name. */
  uint64_t temporal;        /**< Full frames, temporal sampling enabled. */
  uint64_t frame_end;       /**< Full frames, temporal sampling disabled. */
  uint64_t segment_clipped; /**< One display segment's clip per frame. */
};

/**
 * @brief Fingerprint of the floating-point environment the hashes were
 *        recorded under.
 * @details Folds libm transcendental results and the fast-math build flag, so
 * a platform whose libm rounds differently compares no hashes.
 * @return FNV-1a fold of the probe results.
 */
inline uint64_t float_environment_fingerprint() {
  uint64_t hash = HASH_SEED64;
  auto fold = [&](float value) {
    uint32_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    hash = fnv1a64_channel(hash, static_cast<uint16_t>(bits));
    hash = fnv1a64_channel(hash, static_cast<uint16_t>(bits >> 16));
  };
  for (int i = 0; i < 4096; ++i) {
    const float t = static_cast<float>(i) / 4096.0f;
    const float angle = (t * 2.0f - 1.0f) * 12.0f;
    fold(sinf(angle));
    fold(cosf(angle));
    fold(tanf(angle * 0.1f));
    fold(acosf(t * 2.0f - 1.0f));
    fold(asinf(t * 2.0f - 1.0f));
    fold(atan2f(sinf(angle), cosf(angle) + 0.5f));
    fold(expf(angle * 0.5f));
    fold(logf(t + 0.001f));
    fold(powf(t + 0.001f, 2.2f));
    fold(static_cast<float>(std::sin(static_cast<double>(angle))));
    fold(static_cast<float>(std::cos(static_cast<double>(angle))));
    fold(static_cast<float>(std::log(static_cast<double>(t) + 0.001)));
  }
#if HS_TEST_FAST_MATH
  hash = fnv1a64_channel(hash, 1);
#endif
  return hash;
}

/** @brief Fingerprint of the environment that recorded PARITY_HASHES. */
inline constexpr uint64_t RECORDED_FINGERPRINT = 0xa4066a8d264bc108ull;

/** @brief Recorded temporal, frame-end and segment-clipped frame hashes. */
inline constexpr ParityHashes PARITY_HASHES[] = {
    {"Comets", 0x8b8536ade45aba11ull, 0x7d05d3cd68d1cbb9ull,
     0x7c59ab4460ae24e3ull},
    {"DisplacementField", 0xc4bc19266261803dull, 0x47fa576d9a94eb76ull,
     0xc0cc816e5cbabb92ull},
    {"DreamBalls", 0x4949b1bb2f407d5aull, 0x448c45bc7689b9e1ull,
     0xb1ee4826d631b513ull},
    {"Dynamo", 0x3abd50e13f110816ull, 0x3abd50e13f110816ull,
     0x3abd50e13f110816ull},
    {"Fishbowl", 0xe7715a41d1d9b309ull, 0xe470fb4fe3b1eed2ull,
     0x4ea8cf3a5cefad6full},
    {"Galaxies", 0xe5a998917a1e8e25ull, 0x4637663771062f8full,
     0xfeb6f3baf57701b6ull},
    {"GnomonicStars", 0x81bb46212d3052a6ull, 0x16a3021ea493af52ull,
     0x19e77fedacd5f287ull},
    {"HankinSolids", 0x9cb8fc1dd116511dull, 0x530957c44028e83bull,
     0x0cdcd78b958faf19ull},
    {"HopfFibration", 0x6e1f434fa1b2d41bull, 0xde39c32a03287cc3ull,
     0x9382a2af48e7c07cull},
    {"IslamicStars", 0x96b8c46b2f41db54ull, 0x66b52924497c1972ull,
     0x8ef51aa75425b232ull},
    {"MeshFeedback", 0x6275b4b7567b384aull, 0x20674c39f14b693cull,
     0x6275b4b7567b384aull},
    {"MindSplatter", 0x0e9ed49566641abdull, 0x2d2ed377b25b0f51ull,
     0x6271120babdf94f3ull},
    {"MobiusRings", 0x67d2b1e0df5fd6feull, 0x00f0e15d48db5d71ull,
     0x00b3fed29e618219ull},
    {"PetalFlow", 0x04889c42ab9c5e65ull, 0x1312633d0088b40full,
     0xa6d67e4fad3537daull},
    {"RingShower", 0x7624bf028d7c940aull, 0x7624bf028d7c940aull,
     0x685b21e41b9424d0ull},
    {"RingSpin", 0x6c0a0a7a0d9b5f37ull, 0x67e9653ae68abed2ull,
     0xbaff0e1259f62c26ull},
    {"ShapeShifter", 0x776048d2aefa8d00ull, 0x47d2347a023fe8dfull,
     0xdb96e15259c022f0ull},
    {"Thrusters", 0x04886484f4367d8full, 0x04886484f4367d8full,
     0xae12ff3819824aafull},
};

/**
 * @brief Renders PARITY_FRAMES frames of one effect and folds every displayed
 *        frame into one hash.
 * @tparam E Effect class template.
 * @param temporal Temporal sampling selection for the whole run.
 * @param segment_clipped Render one display segment's clip each frame, where
 *        the segmented driver clips the effect.
 * @param lit Cleared when no displayed frame has a lit pixel.
 * @return The folded hash.
 */
template <template <int, int> class E>
inline uint64_t render_fold(bool temporal, bool segment_clipped, bool &lit) {
  reset_globals();
  Spacetime::set_flat_surface(!temporal);
  pin_frame_clock(0);
  E<PARITY_W, PARITY_H> effect;
  effect.init();
  constexpr int S = PARITY_H * 2;
  const pov::SegmentMap map =
      pov::segment_map(PARITY_SEGMENT_ID, S, PARITY_SEGMENTS);
  const bool clip_applies = pov::segment_clip_applies(effect.needs_full_frame(),
                                                      effect.persists_pixels());
  uint64_t fold = HASH_SEED64;
  bool any_lit = false;
  for (int f = 0; f < PARITY_FRAMES; ++f) {
    if (segment_clipped && clip_applies) {
      const pov::SegmentClip clip =
          pov::segment_clip(map, (f & 1) == 0, S, PARITY_SEGMENTS, PARITY_W);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
    }
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    fold = frame_hash<PARITY_W, PARITY_H>(effect, fold);
    any_lit = any_lit || count_lit_region<PARITY_W, PARITY_H>(effect) != 0;
  }
  lit = lit && any_lit;
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return fold;
}

/**
 * @brief Renders every parity scenario of one effect.
 * @tparam E Effect class template.
 * @param name Effect name.
 * @return The effect's hashes.
 */
template <template <int, int> class E>
inline ParityHashes render_hashes(const char *name) {
  bool lit = true;
  const ParityHashes hashes{name, render_fold<E>(true, false, lit),
                            render_fold<E>(false, false, lit),
                            render_fold<E>(true, true, lit)};
  HS_CONTEXT(name);
  HS_EXPECT(lit, "every parity scenario must display a lit pixel");
  return hashes;
}

/**
 * @brief Every filter-pipeline effect matches its recorded frame hashes.
 * @details Hashes compare only under the recorded floating-point environment.
 * Elsewhere the scenarios still render, the computed table is printed and the
 * comparison is a skipped case, except under glibc's libm, which renders
 * without comparing.
 */
inline void test_filter_pipeline_effects_match_recorded_frames() {
  const ParityHashes computed[] = {
      render_hashes<Comets>("Comets"),
      render_hashes<DisplacementField>("DisplacementField"),
      render_hashes<DreamBalls>("DreamBalls"),
      render_hashes<Dynamo>("Dynamo"),
      render_hashes<Fishbowl>("Fishbowl"),
      render_hashes<Galaxies>("Galaxies"),
      render_hashes<GnomonicStars>("GnomonicStars"),
      render_hashes<HankinSolids>("HankinSolids"),
      render_hashes<HopfFibration>("HopfFibration"),
      render_hashes<IslamicStars>("IslamicStars"),
      render_hashes<MeshFeedback>("MeshFeedback"),
      render_hashes<MindSplatter>("MindSplatter"),
      render_hashes<MobiusRings>("MobiusRings"),
      render_hashes<PetalFlow>("PetalFlow"),
      render_hashes<RingShower>("RingShower"),
      render_hashes<RingSpin>("RingSpin"),
      render_hashes<ShapeShifter>("ShapeShifter"),
      render_hashes<Thrusters>("Thrusters"),
  };
  const uint64_t fingerprint = float_environment_fingerprint();
  const bool recorded_environment = fingerprint == RECORDED_FINGERPRINT;
  bool mismatch = !recorded_environment;
  if (recorded_environment) {
    HS_EXPECT_EQ(std::size(computed), std::size(PARITY_HASHES));
    for (size_t i = 0; i < std::size(computed) && i < std::size(PARITY_HASHES);
         ++i) {
      const ParityHashes &want = PARITY_HASHES[i];
      const ParityHashes &got = computed[i];
      HS_CONTEXT(got.effect);
      HS_EXPECT(std::strcmp(want.effect, got.effect) == 0,
                "recorded effect order");
      HS_EXPECT_EQ(got.temporal, want.temporal);
      HS_EXPECT_EQ(got.frame_end, want.frame_end);
      HS_EXPECT_EQ(got.segment_clipped, want.segment_clipped);
      mismatch = mismatch || got.temporal != want.temporal ||
                 got.frame_end != want.frame_end ||
                 got.segment_clipped != want.segment_clipped;
    }
  } else {
    std::printf("  plot_parity: float environment 0x%016llxull differs from "
                "the recorded 0x%016llxull; hashes not compared\n",
                static_cast<unsigned long long>(fingerprint),
                static_cast<unsigned long long>(RECORDED_FINGERPRINT));
#if defined(__GLIBC__) && !defined(HS_TEST_FAST_MATH)
    std::printf("  plot_parity: hashes are recorded against another libm; "
                "glibc renders without comparing\n");
#else
    skip_case(__func__, "recorded frame hashes need the recorded float "
                        "environment");
#endif
  }
  if (mismatch) {
    std::printf("  RECORDED_FINGERPRINT = 0x%016llxull;\n",
                static_cast<unsigned long long>(fingerprint));
    for (const ParityHashes &h : computed)
      std::printf("    {\"%s\", 0x%016llxull, 0x%016llxull, 0x%016llxull},\n",
                  h.effect, static_cast<unsigned long long>(h.temporal),
                  static_cast<unsigned long long>(h.frame_end),
                  static_cast<unsigned long long>(h.segment_clipped));
  }
}

/**
 * @brief Runs every plot parity test case.
 * @return The module's failure count.
 */
inline int run_plot_parity_tests() {
  hs_test::ModuleFixture fixture("plot_parity");
  test_filter_pipeline_effects_match_recorded_frames();
  return fixture.result();
}

} // namespace plot_parity_tests
} // namespace hs_test
