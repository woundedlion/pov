/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Render properties of the plot and scan filter-pipeline effects.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

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

/** @brief Largest relative energy gap between temporal and frame-end runs. */
constexpr double TEMPORAL_ENERGY_BAND = 0.3;
/** @brief Largest brightness-centroid gap between temporal and frame-end
 *  runs, in pixels. */
constexpr double TEMPORAL_CENTROID_BAND = 4.0;
/** @brief Largest share of displayed pixels a chunk-batched effect's clipped
 *  render may differ on: overlapping trails may stack in another order. */
constexpr double BATCHED_REORDER_SHARE = 0.001;

/** @brief Every displayed frame of one parity run. */
struct ParityRun {
  std::vector<Pixel> frames; /**< PARITY_FRAMES whole frames, row-major. */
  RenderFold fold;           /**< Hash and fingerprint of the whole run. */
  bool lit = false;     /**< Whether any displayed frame has a lit pixel. */
  bool clipped = false; /**< Whether the segment clip applied. */
};

/**
 * @brief Display-segment clip the clipped run applies on frame @p f.
 * @param f Frame index.
 * @return The clip window.
 */
inline pov::SegmentClip parity_clip(int f) {
  constexpr int S = PARITY_H * 2;
  const pov::SegmentMap map =
      pov::segment_map(PARITY_SEGMENT_ID, S, PARITY_SEGMENTS);
  return pov::segment_clip(map, (f & 1) == 0, S, PARITY_SEGMENTS, PARITY_W);
}

/**
 * @brief Renders PARITY_FRAMES frames of one effect and keeps every displayed
 *        frame.
 * @tparam E Effect class template.
 * @param temporal Temporal sampling selection for the whole run.
 * @param segment_clipped Render one display segment's clip each frame, where
 *        the segmented driver clips the effect.
 * @return The run's frames and fold.
 */
template <template <int, int> class E>
inline ParityRun render_run(bool temporal, bool segment_clipped) {
  reset_globals();
  Spacetime::set_flat_surface(!temporal);
  pin_frame_clock(0);
  E<PARITY_W, PARITY_H> effect;
  effect.init();
  ParityRun run;
  run.clipped =
      segment_clipped && pov::segment_clip_applies(effect.needs_full_frame(),
                                                   effect.persists_pixels());
  run.frames.reserve(static_cast<size_t>(PARITY_FRAMES) * PARITY_W * PARITY_H);
  for (int f = 0; f < PARITY_FRAMES; ++f) {
    if (run.clipped) {
      const pov::SegmentClip clip = parity_clip(f);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
    }
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    fold_displayed_frame<PARITY_W, PARITY_H>(run.fold, effect);
    run.lit = run.lit || count_lit_region<PARITY_W, PARITY_H>(effect) != 0;
    for (int y = 0; y < PARITY_H; ++y)
      for (int x = 0; x < PARITY_W; ++x)
        run.frames.push_back(effect.get_pixel(x, y));
  }
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return run;
}

/**
 * @brief Counts pixels inside each frame's displayed window where two runs
 *        differ.
 * @param full Run whose frames are whole.
 * @param clipped Run under test.
 * @param lit Set to the number of lit pixels compared.
 * @param displayed Set to the number of pixels compared.
 * @return Differing pixels.
 * @details Without a clip the whole frame is the window.
 */
inline size_t displayed_differences(const ParityRun &full,
                                    const ParityRun &clipped, size_t &lit,
                                    size_t &displayed) {
  size_t different = 0;
  lit = 0;
  displayed = 0;
  for (int f = 0; f < PARITY_FRAMES; ++f) {
    pov::SegmentClip window{0, PARITY_W, 0, PARITY_H};
    if (clipped.clipped)
      window = parity_clip(f);
    const size_t base = static_cast<size_t>(f) * PARITY_W * PARITY_H;
    for (int y = window.y0; y < window.y1; ++y)
      for (int x = window.x0; x < window.x1; ++x) {
        const size_t i = base + static_cast<size_t>(y) * PARITY_W + x;
        different += full.frames[i] != clipped.frames[i];
        lit += (full.frames[i].r | full.frames[i].g | full.frames[i].b) != 0;
        ++displayed;
      }
  }
  return different;
}

/**
 * @brief Checks one effect's render properties across the parity scenarios.
 * @tparam E Effect class template.
 * @param name Effect name.
 * @param batched Whether the effect batches trail edges in a ColumnChunk, whose
 *        flush points, and so the stacking of overlapping trails, follow the
 *        clip.
 * @details Every scenario displays a lit pixel. The segment-clipped run equals
 *          the temporal run inside each frame's displayed window, which is a
 *          determinism check where no clip applies; a batched effect may
 *          differ on at most BATCHED_REORDER_SHARE of the displayed pixels. The temporal and frame-end
 *          runs agree in energy and brightness centroid within coarse bands.
 *          With HS_PRINT_RENDER_HASHES set, each run's hash is printed for
 *          `tools/render_diff.sh`.
 */
template <template <int, int> class E>
inline void check_parity_properties(const char *name, bool batched = false) {
  HS_CONTEXT(name);
  const ParityRun temporal = render_run<E>(true, false);
  const ParityRun frame_end = render_run<E>(false, false);
  const ParityRun clipped = render_run<E>(true, true);
  HS_EXPECT(temporal.lit && frame_end.lit && clipped.lit,
            "every parity scenario must display a lit pixel");

  size_t lit = 0;
  size_t displayed = 0;
  const size_t different =
      displayed_differences(temporal, clipped, lit, displayed);
  if (different)
    std::printf("  CLIP DRIFT %-20s %zu displayed pixels differ from the "
                "unclipped render\n",
                name, different);
  const size_t allowed =
      batched ? static_cast<size_t>(BATCHED_REORDER_SHARE *
                                    static_cast<double>(displayed))
              : 0;
  HS_EXPECT_LE(different, allowed);
  HS_EXPECT_GT(lit, size_t{0});

  const RenderFingerprint &t = temporal.fold.fingerprint;
  const RenderFingerprint &e = frame_end.fold.fingerprint;
  if (e.energy() != 0 && t.energy() != 0) {
    const double energy_gap = static_cast<double>(t.energy()) / e.energy() - 1;
    const double row_gap = static_cast<double>(t.row) / t.energy() -
                           static_cast<double>(e.row) / e.energy();
    const double column_gap = static_cast<double>(t.column) / t.energy() -
                              static_cast<double>(e.column) / e.energy();
    HS_EXPECT_LE(std::abs(energy_gap), TEMPORAL_ENERGY_BAND);
    HS_EXPECT_LE(std::abs(row_gap), TEMPORAL_CENTROID_BAND);
    HS_EXPECT_LE(std::abs(column_gap), TEMPORAL_CENTROID_BAND);
  }

  if (env_var("HS_PRINT_RENDER_HASHES"))
    std::printf("render-hash %s temporal=%016llx frame_end=%016llx "
                "clipped=%016llx\n",
                name, static_cast<unsigned long long>(temporal.fold.hash),
                static_cast<unsigned long long>(frame_end.fold.hash),
                static_cast<unsigned long long>(clipped.fold.hash));
}

/**
 * @brief Each filter-pipeline effect is lit, clip-consistent and agrees
 *        coarsely between temporal and frame-end sampling.
 */
inline void test_filter_pipeline_effect_properties() {
  check_parity_properties<Comets>("Comets");
  check_parity_properties<DisplacementField>("DisplacementField");
  check_parity_properties<DreamBalls>("DreamBalls");
  check_parity_properties<Dynamo>("Dynamo");
  check_parity_properties<Fishbowl>("Fishbowl", true);
  check_parity_properties<Galaxies>("Galaxies");
  check_parity_properties<GnomonicStars>("GnomonicStars");
  check_parity_properties<HankinSolids>("HankinSolids");
  check_parity_properties<HopfFibration>("HopfFibration", true);
  check_parity_properties<IslamicStars>("IslamicStars");
  check_parity_properties<MeshFeedback>("MeshFeedback");
  check_parity_properties<MindSplatter>("MindSplatter", true);
  check_parity_properties<MobiusRings>("MobiusRings");
  check_parity_properties<PetalFlow>("PetalFlow");
  check_parity_properties<Raymarch>("Raymarch");
  check_parity_properties<RingShower>("RingShower");
  check_parity_properties<RingSpin>("RingSpin");
  check_parity_properties<ShapeShifter>("ShapeShifter");
  check_parity_properties<Thrusters>("Thrusters");
}

/**
 * @brief Runs every plot parity test case.
 * @return The module's failure count.
 */
inline int run_plot_parity_tests() {
  hs_test::ModuleFixture fixture("plot_parity");
  HS_RUN(test_filter_pipeline_effect_properties);
  return fixture.result();
}

} // namespace plot_parity_tests
} // namespace hs_test
