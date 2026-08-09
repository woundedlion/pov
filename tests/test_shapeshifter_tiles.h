/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Clipped-tile parity for the ShapeShifter oracle: a mosaic of segment renders
 * must reproduce the unclipped frame pixel for pixel.
 *
 * Separate from the shapeshifter_oracle module because this is the one property
 * there that holds only under IEEE. An active clip switches Plot::rasterize to
 * the planar sampler rebuilt from the cull span; that construction is
 * bit-identical to the unclipped one under IEEE, but the two are inlined into
 * different contexts and reassociate separately under -ffast-math
 * -fno-finite-math-only, which moves a sample across a pixel boundary and with
 * it a whole splat. The divergence is bounded (6.9e-5 of frame energy) but not
 * by anything a channel tolerance can separate from a real geometry defect, so
 * this module is excluded from the fast-math axis and the rest of the oracle —
 * the reference-versus-candidate visual budgets, which compare two renders from
 * the same binary — keeps running there.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

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

template <typename Render>
inline void expect_segment_tiles_reconstruct_full_frame(Render render) {
  const OracleClip clips[] = {{0, ORACLE_H / 2, 0, ORACLE_W / 2},
                              {0, ORACLE_H / 2, ORACLE_W / 2, ORACLE_W},
                              {ORACLE_H / 2, ORACLE_H, 0, ORACLE_W / 2},
                              {ORACLE_H / 2, ORACLE_H, ORACLE_W / 2, ORACLE_W}};
  const auto matrix = exhaustive_matrix();

  for (int shape = 0; shape < 5; ++shape) {
    OracleState full_state = matrix[shape * 4 + shape % 4];
    full_state.orientation = Quaternion();
    OracleFrame full = capture_frame(full_state, render);
    OracleFrame tiled;
    tiled.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);

    for (const OracleClip &clip : clips) {
      OracleState segment_state = full_state;
      segment_state.clip = clip;
      OracleFrame segment = capture_frame(segment_state, render);
      copy_clip(tiled, segment, clip);
    }

    FrameErrorStats error = compare_buffers(full, tiled);
    HS_EXPECT_GT(frame_energy(full), static_cast<uint64_t>(0));
    HS_EXPECT_TRUE(error.exact());
    HS_EXPECT_EQ(error.different_pixels, static_cast<size_t>(0));
    HS_EXPECT_EQ(error.total_absolute_error, static_cast<uint64_t>(0));
  }
}

inline void test_segment_tiles_reconstruct_full_frame() {
  expect_segment_tiles_reconstruct_full_frame(reference_renderer());
  expect_segment_tiles_reconstruct_full_frame(candidate_renderer());

  const OracleClip clips[] = {{0, ORACLE_H / 2, 0, ORACLE_W / 2},
                              {0, ORACLE_H / 2, ORACLE_W / 2, ORACLE_W},
                              {ORACLE_H / 2, ORACLE_H, 0, ORACLE_W / 2},
                              {ORACLE_H / 2, ORACLE_H, ORACLE_W / 2, ORACLE_W}};
  auto expect_tiled = [&](OracleState state) {
    const OracleFrame full = capture_frame(state, candidate_renderer());
    OracleFrame tiled;
    tiled.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
    for (const OracleClip &clip : clips) {
      state.clip = clip;
      copy_clip(tiled, capture_frame(state, candidate_renderer()), clip);
    }
    const FrameErrorStats error = compare_buffers(full, tiled);
    HS_EXPECT_GT(frame_energy(full), uint64_t{0});
    HS_EXPECT_TRUE(error.exact());
    HS_EXPECT_EQ(error.total_absolute_error, uint64_t{0});
  };
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.count = 144;
  state.sides = 7;
  state.phase = 0.37f;
  state.alpha = 0.274f;
  state.orientation = Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
  expect_tiled(state);
  state.phase = 0.125f;
  state.orientation = make_rotation(X_AXIS, Y_AXIS);
  expect_tiled(state);
}

/**
 * @brief Pins the star cap's azimuthal cull against narrow column clips.
 * @details Four full-height W/4 columns tile the canvas, so the y-band half of
 * the cull passes everything and only the azimuthal bound decides visibility.
 * A quarter-width column shrinks the cap half-width to a quarter turn, where the
 * bound rejects most of a dense star stack; any shape it drops that reaches the
 * column shows up as a mismatch against the unclipped frame.
 */
inline void test_star_azimuthal_cull_spans_narrow_columns() {
  const OracleClip columns[] = {{0, ORACLE_H, 0, ORACLE_W / 4},
                                {0, ORACLE_H, ORACLE_W / 4, ORACLE_W / 2},
                                {0, ORACLE_H, ORACLE_W / 2, 3 * ORACLE_W / 4},
                                {0, ORACLE_H, 3 * ORACLE_W / 4, ORACLE_W}};
  const std::array<Quaternion, 4> orientations = {{
      Quaternion(),
      make_rotation(X_AXIS, Z_AXIS),
      Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized(),
      Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
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

    const OracleFrame full = capture_frame(state, candidate_renderer());
    OracleFrame tiled;
    tiled.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
    for (const OracleClip &column : columns) {
      state.clip = column;
      copy_clip(tiled, capture_frame(state, candidate_renderer()), column);
    }

    const FrameErrorStats error = compare_buffers(full, tiled);
    HS_EXPECT_GT(frame_energy(full), uint64_t{0});
    HS_EXPECT_TRUE(error.exact());
    HS_EXPECT_EQ(error.different_pixels, size_t{0});
    HS_EXPECT_EQ(error.total_absolute_error, uint64_t{0});
  }
}

/**
 * @brief Module entry point for the clipped-tile parity sweeps.
 * @return Module result code from hs_test::end_module (0 on success).
 */
inline int run_shapeshifter_tiles_tests() {
  ModuleFixture fixture("shapeshifter_tiles");
  test_segment_tiles_reconstruct_full_frame();
  test_star_azimuthal_cull_spans_narrow_columns();
  return fixture.result();
}

} // namespace shapeshifter_tiles_tests
} // namespace hs_test
