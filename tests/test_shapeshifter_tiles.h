/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Clipped-tile parity for the ShapeShifter oracle: a mosaic of segment renders
 * must reproduce the unclipped frame pixel for pixel. The candidate Flower's
 * band split cuts at clip-dependent points, so its tiles are exempt. Exact
 * parity runs with Plot::PlanarChords' pole-run split off.
 *
 * Parity holds only under IEEE: under -ffast-math the clipped and unclipped
 * planar samplers reassociate differently and can move a whole splat.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/test_shapeshifter_oracle.h"

namespace hs_test {
namespace shapeshifter_tiles_tests {

using namespace hs_test::shapeshifter_oracle_tests;

inline void copy_clip(OracleFrame &destination, const OracleFrame &source,
                      const OracleClip &clip) {
  for (int y = clip.y0; y < clip.y1; ++y)
    for (int x = clip.x0; x < clip.x1; ++x)
      destination.pixels[static_cast<size_t>(y) * ORACLE_W + x] =
          source.at(x, y);
}

inline constexpr OracleClip QUADRANTS[] = {
    {0, ORACLE_H / 2, 0, ORACLE_W / 2},
    {0, ORACLE_H / 2, ORACLE_W / 2, ORACLE_W},
    {ORACLE_H / 2, ORACLE_H, 0, ORACLE_W / 2},
    {ORACLE_H / 2, ORACLE_H, ORACLE_W / 2, ORACLE_W}};

template <typename Render>
inline void expect_mosaic_matches(OracleState state, Render render,
                                  std::span<const OracleClip> clips,
                                  bool approximate = false) {
  const OracleFrame full = capture_frame(state, render);
  OracleFrame tiled;
  tiled.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
  for (const OracleClip &clip : clips) {
    state.clip = clip;
    copy_clip(tiled, capture_frame(state, render), clip);
  }
  if (approximate) {
    const uint64_t ENERGY = frame_energy(full);
    HS_EXPECT_GT(ENERGY, uint64_t{0});
    HS_EXPECT_LT(std::fabs(static_cast<double>(frame_energy(tiled)) - ENERGY) /
                     ENERGY,
                 0.035);
    size_t uncovered = 0;
    for (int y = 3; y < ORACLE_H - 3; ++y)
      for (int x = 0; x < ORACLE_W; ++x) {
        if (pixel_is_bright(full.at(x, y)))
          uncovered += !candidate_covers_neighborhood(tiled, x, y);
      }
    HS_EXPECT_EQ(uncovered, size_t{0});
    return;
  }
  const FrameErrorStats error = compare_buffers(full, tiled);
  HS_EXPECT_GT(frame_energy(full), uint64_t{0});
  HS_EXPECT_TRUE(error.exact());
  HS_EXPECT_EQ(error.different_pixels, size_t{0});
  HS_EXPECT_EQ(error.total_absolute_error, uint64_t{0});
}

template <typename Render>
inline void
expect_segment_tiles_reconstruct_full_frame(Render render,
                                            bool band_split_flower) {
  const auto matrix = shape_function_matrix();
  for (int shape = 0; shape < 5; ++shape) {
    OracleState state = matrix[shape * 4 + shape % 4];
    state.orientation = math::Quaternion();
    expect_mosaic_matches(state, render, QUADRANTS,
                          band_split_flower &&
                              state.shape == OracleEffect::ShapeType::FLOWER);
  }
}

inline void test_segment_tiles_reconstruct_full_frame() {
  expect_segment_tiles_reconstruct_full_frame(reference_renderer(), false);
  expect_segment_tiles_reconstruct_full_frame(candidate_renderer(), true);

  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.count = 144;
  state.sides = 7;
  state.phase = 0.37f;
  state.alpha = 0.274f;
  state.orientation =
      math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
  expect_mosaic_matches(state, candidate_renderer(), QUADRANTS);
  state.phase = 0.125f;
  state.orientation = math::make_rotation(math::X_AXIS, math::Y_AXIS);
  expect_mosaic_matches(state, candidate_renderer(), QUADRANTS);

  // Dense screen-balanced contours at both poles, with pole-run splitting off.
  state.count = 288;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  state.phase = 0.249f;
  expect_mosaic_matches(state, candidate_renderer(), QUADRANTS);
  state.orientation =
      math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized();
  expect_mosaic_matches(state, candidate_renderer(), QUADRANTS);
}

/**
 * @brief Pins the star cap's azimuthal cull against narrow column clips.
 * @details Full-height W/4 columns leave only the azimuthal bound deciding
 * visibility.
 */
inline void test_star_azimuthal_cull_spans_narrow_columns() {
  const OracleClip columns[] = {{0, ORACLE_H, 0, ORACLE_W / 4},
                                {0, ORACLE_H, ORACLE_W / 4, ORACLE_W / 2},
                                {0, ORACLE_H, ORACLE_W / 2, 3 * ORACLE_W / 4},
                                {0, ORACLE_H, 3 * ORACLE_W / 4, ORACLE_W}};
  const std::array<math::Quaternion, 4> orientations = {{
      math::Quaternion(),
      math::make_rotation(math::X_AXIS, math::Z_AXIS),
      math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized(),
      math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
  }};
  const std::array<float, 4> phases = {{0.0f, 0.37f, 0.5f, 0.83f}};

  for (size_t i = 0; i < orientations.size(); ++i) {
    OracleState state;
    state.shape = OracleEffect::ShapeType::PLANAR_STAR;
    state.function = OracleEffect::PhaseFunction::SINE;
    state.count = 144;
    state.sides = 7;
    state.phase = phases[i];
    state.alpha = 0.274f;
    state.orientation = orientations[i];

    expect_mosaic_matches(state, candidate_renderer(), columns);
  }
}

inline void test_split_pole_runs_tiles_within_energy_budget() {
  Plot::g_planar_chords_split_pole_runs = true;
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.count = 288;
  state.sides = 7;
  state.alpha = 0.274f;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  state.phase = 0.249f;
  state.orientation = math::make_rotation(math::X_AXIS, math::Y_AXIS);
  expect_mosaic_matches(state, candidate_renderer(), QUADRANTS, true);
}

/**
 * @brief Module entry point for the clipped-tile parity sweeps.
 * @return Module result code from hs_test::end_module (0 on success).
 */
inline int run_shapeshifter_tiles_tests() {
  ModuleFixture fixture("shapeshifter_tiles");
  struct PoleRunScope {
    bool saved = Plot::g_planar_chords_split_pole_runs;
    ~PoleRunScope() { Plot::g_planar_chords_split_pole_runs = saved; }
  } scope;
  Plot::g_planar_chords_split_pole_runs = false;
  test_segment_tiles_reconstruct_full_frame();
  test_star_azimuthal_cull_spans_narrow_columns();
  test_split_pole_runs_tiles_within_energy_budget();
  return fixture.result();
}

} // namespace shapeshifter_tiles_tests
} // namespace hs_test
