/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * MindSplatter white-box invariants: emitter/attractor mesh selection, the
 * golden geometry replay, the fused-vertex, hole-kernel, signed-axis and
 * rotation-matrix framebuffer parities, clip/clear display parity, preset
 * timeline bookkeeping and emission-phase wrapping.
 *
 * Split out of the effects module so the roster's heaviest single-effect block
 * is its own CTest and can be sharded independently. The cases share the
 * effects module's render fixtures and resolution constants, hence the
 * test_effects.h include and the using-directive below.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "tests/mindsplatter_whitebox.h"
#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mindsplatter_tests {

using namespace hs_test::effects_tests;

struct FrameRun {
  std::vector<Pixel> pixels;
  std::vector<uint16_t> active;
};

struct FrameDiff {
  size_t lit = 0;
  size_t different = 0;
  size_t coverage = 0;
  int max_channel = 0;
  uint64_t total_channel = 0;
};

template <int W, int H, int FRAMES, typename Toggle>
std::array<FrameRun, 2> render_toggled(Toggle toggle) {
  std::array<FrameRun, 2> runs;
  for (int reference = 0; reference < 2; ++reference) {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    auto &run = runs[reference];
    run.pixels.reserve(static_cast<size_t>(W) * H * FRAMES);
    run.active.reserve(FRAMES);
    MindSplatter<W, H> effect;
    effect.init();
    toggle(effect, reference != 0);
    for (int frame = 0; frame < FRAMES; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame) * FRAME_MS,
                        static_cast<unsigned long>(frame) * FRAME_US);
      effect.draw_frame();
      effect.advance_display();
      run.active.push_back(MindSplatterWhiteBox::active_particles(effect));
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          run.pixels.push_back(effect.get_pixel(x, y));
    }
  }
  hs::clear_mock_time();
  return runs;
}

inline FrameDiff diff_frames(std::span<const Pixel> reference,
                             std::span<const Pixel> candidate) {
  HS_EXPECT_EQ(reference.size(), candidate.size());
  FrameDiff diff;
  for (size_t i = 0; i < reference.size(); ++i) {
    const Pixel a = reference[i];
    const Pixel b = candidate[i];
    const bool a_black = (a.r | a.g | a.b) == 0;
    const bool b_black = (b.r | b.g | b.b) == 0;
    diff.lit += !a_black;
    diff.different += a != b;
    diff.coverage += a_black != b_black;
    for (int delta :
         {std::abs(static_cast<int>(a.r) - static_cast<int>(b.r)),
          std::abs(static_cast<int>(a.g) - static_cast<int>(b.g)),
          std::abs(static_cast<int>(a.b) - static_cast<int>(b.b))}) {
      diff.max_channel = std::max(diff.max_channel, delta);
      diff.total_channel += static_cast<uint64_t>(delta);
    }
  }
  return diff;
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
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 1), MS::BaseMesh::CUBE);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 2), MS::BaseMesh::CUBE);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 4), MS::BaseMesh::TETRAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 5), MS::BaseMesh::DODECAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 6), MS::BaseMesh::CUBE);

  const auto select = [&](MS::BaseMesh mesh, size_t emitters,
                          size_t attractors) {
    HS_EXPECT_EQ(effect.updateParameter("Base Mesh", static_cast<float>(mesh)),
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

/** @brief The whitebox search and replay retain the selected particle geometry. */
inline void test_mindsplatter_whitebox_geometry_replay() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  std::vector<unsigned char> state;

  {
    reset_effect_globals();
    MS effect;
    effect.init();
    WB::select_preset(effect, 5);
    effect.setAnimationsPaused(true);
    WB::step_state_without_render(effect);
    HS_EXPECT_EQ(WB::active_base_mesh(effect), MS::BaseMesh::DODECAHEDRON);
    HS_EXPECT_EQ(WB::active_emitters(effect), static_cast<size_t>(20));
    HS_EXPECT_EQ(WB::active_attractors(effect), static_cast<size_t>(12));
    state = WB::serialize_render(WB::capture(effect));
  }

  reset_effect_globals();
  MS replay;
  replay.init();
  WB::restore_render(replay, state);
  HS_EXPECT_EQ(WB::active_base_mesh(replay), MS::BaseMesh::DODECAHEDRON);
  HS_EXPECT_EQ(WB::active_emitters(replay), static_cast<size_t>(20));
  HS_EXPECT_EQ(WB::active_attractors(replay), static_cast<size_t>(12));
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
 * @brief MindSplatter's single-pass direct-AA path matches the cached
 * AntiAlias-sink reference in every quadrant of a frozen saturated particle
 * pool.
 * @details Single-pass stepping omits the cached endpoint normalization, so
 * interior sample phases, one fringe pixel, and accumulated channels can
 * differ; the bounds are the measured worst case (exact at 96x20 under IEEE)
 * with headroom. The production-resolution replay executable compares the
 * same candidate against draw_particles_replay_reference, the shipping filter
 * stack under reference screen stepping, not against this sink.
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
    HS_EXPECT_LE(changed_pixels, static_cast<size_t>(4));
    HS_EXPECT_LE(max_channel_error, 64);
    HS_EXPECT_LE(total_channel_error, static_cast<uint64_t>(128));
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
  const auto runs = render_toggled<W, H, FRAMES>([](MS &effect, bool value) {
    WB::use_reference_orientation(effect, value);
  });
  const FrameRun &matrix = runs[0];
  const FrameRun &reference = runs[1];
  const FrameDiff diff = diff_frames(reference.pixels, matrix.pixels);
  std::printf("matrix framebuffer samples=%zu different=%zu coverage=%zu "
              "max_channel=%d total_channel=%llu\n",
              reference.pixels.size(), diff.different, diff.coverage,
              diff.max_channel,
              static_cast<unsigned long long>(diff.total_channel));
  HS_EXPECT_EQ(diff.coverage, static_cast<size_t>(0));
  HS_EXPECT_LE(diff.different, static_cast<size_t>(96));
  HS_EXPECT_LE(diff.max_channel, 8);
  HS_EXPECT_LE(diff.total_channel, static_cast<uint64_t>(128));
}

/** @brief Particle hue seeds advance in deterministic emission order. */
inline void test_mindsplatter_particle_gradients_follow_emission_order() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  MS effect;
  effect.init();
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
  HS_EXPECT_EQ(live.base_mesh, chosen.base_mesh);
  HS_EXPECT_NEAR(live.friction, chosen.friction, 1e-6f);
  HS_EXPECT_NEAR(live.well_strength, chosen.well_strength, 1e-6f);
  HS_EXPECT_NEAR(live.initial_speed, chosen.initial_speed, 1e-6f);
  HS_EXPECT_NEAR(live.angular_speed, chosen.angular_speed, 1e-6f);
  HS_EXPECT_NEAR(live.warp_scale, chosen.warp_scale, 1e-6f);
}

inline void test_mindsplatter_full_timeline_retries_transition() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();
  MS effect;
  effect.init();
  float sink = 0.0f;
  WB::saturate_timeline(effect, sink);
  const uint32_t dropped_before = Timeline::dropped_events();

  constexpr int INITIAL_DWELL =
      MS::PRESET_DWELL_FRAMES + MS::PRESET_SEGUE.frames;
  for (int f = 0; f < INITIAL_DWELL; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);
  HS_EXPECT_FALSE(WB::transition_active(effect));
  HS_EXPECT_EQ(WB::preset_index(effect), size_t{0});

  // The rejection restarts the dwell, so the retry costs one drop per dwell
  // rather than one per frame.
  for (int f = 1; f < MS::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);

  WB::clear_timeline(effect);
  WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(WB::preset_index(effect), size_t{1});
}

/** @brief Fusing the vertex pass into trail materialization is pixel exact. */
inline void test_mindsplatter_fused_vertex_framebuffer_parity() {
  constexpr int W = SMALL_W;
  constexpr int H = SMALL_H;
  constexpr int FRAMES = 16;
  using MS = MindSplatter<W, H>;
  using WB = MindSplatterWhiteBox;
  const auto runs = render_toggled<W, H, FRAMES>([](MS &effect, bool value) {
    WB::use_reference_vertex_pass(effect, value);
  });
  const FrameDiff diff = diff_frames(runs[1].pixels, runs[0].pixels);
  std::printf("fused vertex framebuffer samples=%zu lit=%zu different=%zu\n",
              runs[1].pixels.size(), diff.lit, diff.different);
  HS_EXPECT_GT(diff.lit, static_cast<size_t>(0));
  HS_EXPECT_EQ(diff.different, static_cast<size_t>(0));
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
  const auto runs = render_toggled<W, H, FRAMES>([](MS &effect, bool value) {
    WB::use_reference_hole_kernel(effect, value);
  });
  const FrameDiff diff = diff_frames(runs[1].pixels, runs[0].pixels);
  std::printf("hole kernel framebuffer samples=%zu lit=%zu different=%zu "
              "coverage=%zu max_channel=%d total_channel=%llu\n",
              runs[1].pixels.size(), diff.lit, diff.different, diff.coverage,
              diff.max_channel,
              static_cast<unsigned long long>(diff.total_channel));
  HS_EXPECT_GT(diff.lit, static_cast<size_t>(0));
  HS_EXPECT_EQ(diff.coverage, static_cast<size_t>(0));
  HS_EXPECT_LE(diff.different, static_cast<size_t>(64));
  HS_EXPECT_LE(diff.max_channel, 8);
  HS_EXPECT_LE(diff.total_channel, static_cast<uint64_t>(64));
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
  const auto runs = render_toggled<W, H, FRAMES>([](MS &effect, bool value) {
    WB::use_reference_signed_axis_physics(effect, value);
  });
  const FrameRun &specialized = runs[0];
  const FrameRun &reference = runs[1];
  HS_EXPECT_EQ(reference.active.size(), specialized.active.size());
  for (size_t i = 0; i < reference.active.size(); ++i)
    HS_EXPECT_EQ(reference.active[i], specialized.active[i]);

  constexpr int CHECKPOINTS[] = {16, 80, 160};
  constexpr size_t MAX_DIFFERENT[] = {1, 192, 192};
  constexpr int MAX_CHANNEL[] = {1, 256, 1152};
  constexpr uint64_t MAX_TOTAL[] = {1, 2048, 4608};
  for (size_t checkpoint = 0; checkpoint < 3; ++checkpoint) {
    const int frame = CHECKPOINTS[checkpoint];
    const size_t offset = static_cast<size_t>(frame - 1) * W * H;
    const FrameDiff diff =
        diff_frames(std::span(reference.pixels).subspan(offset, W * H),
                    std::span(specialized.pixels).subspan(offset, W * H));
    std::printf(
        "axis framebuffer frame=%d active=%u different=%zu coverage=%zu "
        "max_channel=%d total_channel=%llu\n",
        frame, reference.active[frame - 1], diff.different, diff.coverage,
        diff.max_channel, static_cast<unsigned long long>(diff.total_channel));
    HS_EXPECT_EQ(diff.coverage, static_cast<size_t>(0));
    HS_EXPECT_LE(diff.different, MAX_DIFFERENT[checkpoint]);
    HS_EXPECT_LE(diff.max_channel, MAX_CHANNEL[checkpoint]);
    HS_EXPECT_LE(diff.total_channel, MAX_TOTAL[checkpoint]);
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

/** @brief Attractor caps match independent quintic radial samples. */
inline void test_mindsplatter_attractor_hole_alpha_equivalence() {
  using MS = MindSplatter<SMALL_W, SMALL_H>;
  using WB = MindSplatterWhiteBox;
  reset_effect_globals();

  MS effect;
  effect.init();
  constexpr MS::BaseMesh MESHES[] = {
      MS::BaseMesh::TETRAHEDRON, MS::BaseMesh::OCTAHEDRON,
      MS::BaseMesh::DODECAHEDRON, MS::BaseMesh::ICOSAHEDRON};
  for (MS::BaseMesh mesh : MESHES) {
    HS_EXPECT_EQ(effect.updateParameter("Base Mesh", static_cast<float>(mesh)),
                 ParamSetResult::APPLIED);
    effect.draw_frame();
    effect.advance_display();

    for (size_t i = 0; i < WB::active_attractors(effect); ++i) {
      const Vector attractor = WB::attractor_position(effect, i);
      Vector tangent = cross(attractor, Y_AXIS);
      if (tangent.length() < 0.5f)
        tangent = cross(attractor, X_AXIS);
      tangent = tangent.normalized();
      constexpr float FRACTIONS[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.25f};
      constexpr float EXPECTED[] = {0.0f, 0.103515625f, 0.5f, 0.896484375f,
                                    1.0f};
      for (size_t sample = 0; sample < std::size(FRACTIONS); ++sample) {
        const float angle = WB::event_horizon() * FRACTIONS[sample];
        const Vector p = attractor * cosf(angle) + tangent * sinf(angle);
        HS_EXPECT_NEAR(WB::attractor_hole_alpha(effect, p), EXPECTED[sample],
                       1e-3f);
      }
    }
  }
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
  HS_EXPECT_EQ(ms.updateParameter("Ang Spd", 1.0f), ParamSetResult::APPLIED);

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
 * @brief Module entry point for the MindSplatter white-box block.
 * @return Module result code from hs_test::end_module (0 on success).
 */
inline int run_mindsplatter_tests() {
  ModuleFixture fixture("mindsplatter");

  test_mindsplatter_base_mesh_selector();
  test_mindsplatter_whitebox_geometry_replay();
  test_mindsplatter_replay_snapshot_exact();
  test_mindsplatter_saturated_quadrant_sink_parity();
  test_mindsplatter_octahedral_hole_alpha_equivalence();
  test_mindsplatter_attractor_hole_alpha_equivalence();
  test_mindsplatter_rotation_matrix_equivalence();
  test_mindsplatter_rotation_matrix_framebuffer_error();
  test_mindsplatter_particle_gradients_follow_emission_order();
  test_mindsplatter_fused_vertex_framebuffer_parity();
  test_mindsplatter_hole_kernel_framebuffer_parity();
  test_mindsplatter_clip_clear_display_parity();
  test_mindsplatter_signed_axis_framebuffer_error();
  test_mindsplatter_manual_preset_survives_unpause();
  test_mindsplatter_full_timeline_retries_transition();

  // FULL tier only (HS_EFFECTS_FULL=1), matching the effects module's split.
  if (effects_full_suite())
    test_mindsplatter_emit_phase_wrapped();

  return fixture.result();
}

} // namespace mindsplatter_tests
} // namespace hs_test
