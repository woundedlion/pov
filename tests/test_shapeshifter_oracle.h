/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "effects/ShapeShifter.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shapeshifter_oracle_tests {

constexpr int ORACLE_W = 288;
constexpr int ORACLE_H = 144;
constexpr uint32_t BRIGHT_ENERGY = 12288;
constexpr uint32_t COVERAGE_ENERGY = 512;
constexpr uint32_t HIGH_CHANNEL_ERROR = 4096;
constexpr double MAX_MEAN_ABSOLUTE_ERROR = 350.0;
constexpr double MAX_ROOT_MEAN_SQUARED_ERROR = 1250.0;
constexpr uint32_t MAX_CHANNEL_ERROR = 32768;
constexpr double MAX_ENERGY_DRIFT = 0.006;
constexpr size_t MAX_HIGH_ERROR_PIXELS = 1600;
using OracleEffect = ShapeShifter<ORACLE_W, ORACLE_H>;

/** @brief Display bounds used by one oracle render. */
struct OracleClip {
  int y0 = 0;
  int y1 = ORACLE_H;
  int x0 = 0;
  int x1 = ORACLE_W;
};

/** @brief Fully deterministic input state for a production-resolution render. */
struct OracleState {
  OracleEffect::ShapeType shape = OracleEffect::ShapeType::PLANAR_POLYGON;
  OracleEffect::PhaseFunction function = OracleEffect::PhaseFunction::SINE;
  int count = 7;
  int sides = 5;
  float phase = 0.0f;
  float alpha = 0.625f;
  bool opposite = false;
  Quaternion orientation;
  OracleClip clip;
};

/** @brief Row-major production-resolution framebuffer capture. */
struct OracleFrame {
  std::vector<Pixel> pixels;
  uint32_t backstop_hits =
      0; /**< plot() steps_cache backstop trips during this capture. */

  const Pixel &at(int x, int y) const {
    return pixels[static_cast<size_t>(y) * ORACLE_W + x];
  }
};

/** @brief Exact and aggregate channel-error statistics for two captures. */
struct FrameErrorStats {
  size_t compared_pixels = 0;
  size_t different_pixels = 0;
  size_t different_channels = 0;
  uint32_t max_absolute_error = 0;
  uint64_t total_absolute_error = 0;
  uint64_t total_squared_error = 0;
  int first_different_pixel = -1;

  bool exact() const { return different_pixels == 0; }

  double mean_absolute_error() const {
    const size_t channels = compared_pixels * 3;
    return channels == 0 ? 0.0
                         : static_cast<double>(total_absolute_error) / channels;
  }

  double root_mean_squared_error() const {
    const size_t channels = compared_pixels * 3;
    return channels == 0
               ? 0.0
               : std::sqrt(static_cast<double>(total_squared_error) / channels);
  }
};

/** @brief Compares two same-sized captures channel by channel. */
inline FrameErrorStats compare_buffers(const OracleFrame &reference,
                                       const OracleFrame &candidate) {
  HS_CHECK(reference.pixels.size() == candidate.pixels.size());
  FrameErrorStats stats;
  stats.compared_pixels = reference.pixels.size();

  for (size_t i = 0; i < reference.pixels.size(); ++i) {
    const Pixel &a = reference.pixels[i];
    const Pixel &b = candidate.pixels[i];
    const uint16_t channels_a[] = {a.r, a.g, a.b};
    const uint16_t channels_b[] = {b.r, b.g, b.b};
    bool pixel_differs = false;
    for (int channel = 0; channel < 3; ++channel) {
      const uint32_t error = channels_a[channel] > channels_b[channel]
                                 ? channels_a[channel] - channels_b[channel]
                                 : channels_b[channel] - channels_a[channel];
      if (error == 0)
        continue;
      pixel_differs = true;
      ++stats.different_channels;
      stats.max_absolute_error = std::max(stats.max_absolute_error, error);
      stats.total_absolute_error += error;
      stats.total_squared_error += static_cast<uint64_t>(error) * error;
    }
    if (!pixel_differs)
      continue;
    if (stats.first_different_pixel < 0)
      stats.first_different_pixel = static_cast<int>(i);
    ++stats.different_pixels;
  }
  return stats;
}

/** @brief Test-only access to ShapeShifter's pinned state and shipping renderer. */
struct ShapeShifterWhiteBox {
  static void configure(OracleEffect &effect, const OracleState &state) {
    effect.params.alpha = state.alpha;
    effect.params.shape = static_cast<float>(state.shape);
    effect.params.count = static_cast<float>(state.count);
    effect.params.sides = static_cast<float>(state.sides);
    effect.params.function = static_cast<float>(state.function);
    effect.params.speed = 0.0f;
    effect.params.opposite = state.opposite;
    effect.phase = state.phase;
    effect.orientation.set(state.orientation);
    effect.set_clip(state.clip.y0, state.clip.y1, state.clip.x0, state.clip.x1);
  }

  static void render_reference(OracleEffect &effect, Canvas &canvas) {
    effect.plot_filters.prepare(canvas);
    effect.draw_all_reference(canvas);
  }

  static void render_candidate(OracleEffect &effect, Canvas &canvas) {
    effect.plot_filters.prepare(canvas);
    effect.draw_all(canvas);
  }

  static float advance_phase(OracleEffect &effect, float speed,
                             float amplitude) {
    effect.phase = 0.0f;
    effect.params.speed = speed;
    effect.params.amplitude = amplitude;
    effect.advance_phase();
    return effect.phase;
  }

  static float phase_direction(OracleEffect &effect, bool opposite,
                               float radius) {
    effect.params.opposite = opposite;
    return effect.phase_direction(radius);
  }

  static void next_preset(OracleEffect &effect) { effect.next_preset(); }
};

/** @brief Captures one renderer callable from a fresh deterministic effect. */
template <typename Render>
inline OracleFrame capture_frame(const OracleState &state, Render &&render) {
  reset_globals();
  GenerativePalette::reset_hue_seed(0);
  hs::set_mock_time(0, 0);
  hs::g_scan_metrics.reset();

  OracleFrame frame;
  {
    OracleEffect effect;
    effect.init();
    ShapeShifterWhiteBox::configure(effect, state);
    {
      Canvas canvas(effect);
      std::forward<Render>(render)(effect, canvas);
    }
    effect.advance_display();
    frame.backstop_hits = hs::g_scan_metrics.plot_backstop_hits;
    frame.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
    for (int y = 0; y < ORACLE_H; ++y)
      for (int x = 0; x < ORACLE_W; ++x)
        frame.pixels[static_cast<size_t>(y) * ORACLE_W + x] =
            effect.get_pixel(x, y);
  }
  Timeline().clear();
  hs::clear_mock_time();
  return frame;
}

struct RenderComparison {
  OracleFrame reference;
  OracleFrame candidate;
  FrameErrorStats error;
};

/** @brief Captures and compares two renderer callables from identical state. */
template <typename ReferenceRender, typename CandidateRender>
inline RenderComparison compare_renders(const OracleState &state,
                                        ReferenceRender &&reference_render,
                                        CandidateRender &&candidate_render) {
  RenderComparison comparison;
  comparison.reference =
      capture_frame(state, std::forward<ReferenceRender>(reference_render));
  comparison.candidate =
      capture_frame(state, std::forward<CandidateRender>(candidate_render));
  comparison.error =
      compare_buffers(comparison.reference, comparison.candidate);
  return comparison;
}

inline auto reference_renderer() {
  return [](OracleEffect &effect, Canvas &canvas) {
    ShapeShifterWhiteBox::render_reference(effect, canvas);
  };
}

inline auto candidate_renderer() {
  return [](OracleEffect &effect, Canvas &canvas) {
    ShapeShifterWhiteBox::render_candidate(effect, canvas);
  };
}

inline uint64_t frame_energy(const OracleFrame &frame) {
  uint64_t energy = 0;
  for (const Pixel &pixel : frame.pixels)
    energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
  return energy;
}

inline std::array<OracleState, 16> exhaustive_matrix() {
  using Function = OracleEffect::PhaseFunction;
  using Shape = OracleEffect::ShapeType;
  const Shape shapes[] = {Shape::PLANAR_POLYGON, Shape::SPHERICAL_POLYGON,
                          Shape::FLOWER, Shape::STAR};
  const Function functions[] = {Function::SINE, Function::TRIANGLE,
                                Function::SAWTOOTH, Function::SQUARE};
  const int counts[] = {1, 2, 7, 75};
  const int sides[] = {3, 5, 9, 16};
  const float phases[] = {0.0f, 0.249f, 0.5f, 0.999f};
  const Quaternion orientations[] = {
      Quaternion(), make_rotation(X_AXIS, Y_AXIS),
      make_rotation(X_AXIS, -Y_AXIS), make_rotation(X_AXIS, Z_AXIS)};

  std::array<OracleState, 16> matrix;
  size_t index = 0;
  for (int shape = 0; shape < 4; ++shape)
    for (int function = 0; function < 4; ++function) {
      OracleState state;
      state.shape = shapes[shape];
      state.function = functions[function];
      state.count = counts[function];
      state.sides = sides[shape];
      state.phase = phases[function];
      state.orientation = orientations[(shape + function) % 4];
      matrix[index++] = state;
    }
  return matrix;
}

inline void test_buffer_comparator_statistics() {
  OracleFrame reference;
  reference.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
  OracleFrame candidate = reference;

  FrameErrorStats exact = compare_buffers(reference, candidate);
  HS_EXPECT_TRUE(exact.exact());
  HS_EXPECT_EQ(exact.compared_pixels, static_cast<size_t>(ORACLE_W) * ORACLE_H);
  HS_EXPECT_EQ(exact.different_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(exact.different_channels, static_cast<size_t>(0));
  HS_EXPECT_EQ(exact.max_absolute_error, static_cast<uint32_t>(0));
  HS_EXPECT_EQ(exact.total_absolute_error, static_cast<uint64_t>(0));
  HS_EXPECT_EQ(exact.total_squared_error, static_cast<uint64_t>(0));
  HS_EXPECT_EQ(exact.first_different_pixel, -1);
  HS_EXPECT_EQ(exact.mean_absolute_error(), 0.0);
  HS_EXPECT_EQ(exact.root_mean_squared_error(), 0.0);

  constexpr size_t PIXEL_INDEX = 17;
  candidate.pixels[PIXEL_INDEX] = Pixel(7, 11, 0);
  FrameErrorStats changed = compare_buffers(reference, candidate);
  HS_EXPECT_FALSE(changed.exact());
  HS_EXPECT_EQ(changed.different_pixels, static_cast<size_t>(1));
  HS_EXPECT_EQ(changed.different_channels, static_cast<size_t>(2));
  HS_EXPECT_EQ(changed.max_absolute_error, static_cast<uint32_t>(11));
  HS_EXPECT_EQ(changed.total_absolute_error, static_cast<uint64_t>(18));
  HS_EXPECT_EQ(changed.total_squared_error, static_cast<uint64_t>(170));
  HS_EXPECT_EQ(changed.first_different_pixel, static_cast<int>(PIXEL_INDEX));
}

inline void test_reference_matrix_is_exact_and_nonblack() {
  for (const OracleState &state : exhaustive_matrix()) {
    RenderComparison comparison =
        compare_renders(state, reference_renderer(), reference_renderer());
    HS_EXPECT_TRUE(comparison.error.exact());
    HS_EXPECT_EQ(comparison.error.different_pixels, static_cast<size_t>(0));
    HS_EXPECT_EQ(comparison.error.total_absolute_error,
                 static_cast<uint64_t>(0));
    HS_EXPECT_EQ(comparison.reference.pixels.size(),
                 static_cast<size_t>(ORACLE_W) * ORACLE_H);
    HS_EXPECT_GT(frame_energy(comparison.reference), static_cast<uint64_t>(0));
  }
}

inline bool pixel_is_bright(const Pixel &pixel) {
  return static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > BRIGHT_ENERGY;
}

inline bool pixel_has_coverage(const Pixel &pixel) {
  return static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > COVERAGE_ENERGY;
}

inline bool candidate_covers_neighborhood(const OracleFrame &candidate, int x,
                                          int y) {
  for (int dy = -1; dy <= 1; ++dy) {
    const int sample_y = y + dy;
    if (sample_y < 0 || sample_y >= ORACLE_H)
      continue;
    for (int dx = -1; dx <= 1; ++dx) {
      const int sample_x = (x + dx + ORACLE_W) % ORACLE_W;
      if (pixel_has_coverage(candidate.at(sample_x, sample_y)))
        return true;
    }
  }
  return false;
}

inline void test_candidate_matrix_stays_within_visual_budget() {
  for (const OracleState &state : exhaustive_matrix()) {
    RenderComparison comparison =
        compare_renders(state, reference_renderer(), candidate_renderer());
    const uint64_t reference_energy = frame_energy(comparison.reference);
    const int64_t energy_delta =
        static_cast<int64_t>(frame_energy(comparison.candidate)) -
        static_cast<int64_t>(reference_energy);
    const double energy_ratio =
        reference_energy == 0
            ? 0.0
            : std::fabs(static_cast<double>(energy_delta)) / reference_energy;
    size_t uncovered_bright_pixels = 0;
    size_t high_error_pixels = 0;
    for (size_t pixel = 0; pixel < comparison.reference.pixels.size();
         ++pixel) {
      const Pixel &a = comparison.reference.pixels[pixel];
      const Pixel &b = comparison.candidate.pixels[pixel];
      if (pixel_is_bright(a)) {
        const int x = static_cast<int>(pixel % ORACLE_W);
        const int y = static_cast<int>(pixel / ORACLE_W);
        if (!candidate_covers_neighborhood(comparison.candidate, x, y))
          ++uncovered_bright_pixels;
      }
      const uint32_t max_error = std::max({a.r > b.r ? a.r - b.r : b.r - a.r,
                                           a.g > b.g ? a.g - b.g : b.g - a.g,
                                           a.b > b.b ? a.b - b.b : b.b - a.b});
      if (max_error > HIGH_CHANNEL_ERROR)
        ++high_error_pixels;
    }
    HS_EXPECT_LT(comparison.error.mean_absolute_error(),
                 MAX_MEAN_ABSOLUTE_ERROR);
    HS_EXPECT_LT(comparison.error.root_mean_squared_error(),
                 MAX_ROOT_MEAN_SQUARED_ERROR);
    HS_EXPECT_LT(comparison.error.max_absolute_error, MAX_CHANNEL_ERROR);
    HS_EXPECT_LT(energy_ratio, MAX_ENERGY_DRIFT);
    HS_EXPECT_LT(high_error_pixels, MAX_HIGH_ERROR_PIXELS);
    HS_EXPECT_EQ(uncovered_bright_pixels, size_t{0});
    HS_EXPECT_EQ(comparison.reference.backstop_hits, uint32_t{0});
    HS_EXPECT_EQ(comparison.candidate.backstop_hits, uint32_t{0});
  }
}

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

  for (int shape = 0; shape < 4; ++shape) {
    OracleState full_state = matrix[shape * 4 + shape];
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
    HS_EXPECT_TRUE(error.exact());
    HS_EXPECT_EQ(error.different_pixels, static_cast<size_t>(0));
    HS_EXPECT_EQ(error.total_absolute_error, static_cast<uint64_t>(0));
  }
}

inline void test_segment_tiles_reconstruct_full_frame() {
  expect_segment_tiles_reconstruct_full_frame(reference_renderer());
  expect_segment_tiles_reconstruct_full_frame(candidate_renderer());
}

inline void test_amplitude_preserves_sweep_velocity() {
  OracleEffect effect;
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 1.0f),
                 0.04f, 1e-6f);
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 4.0f),
                 0.01f, 1e-6f);
}

inline void test_opposite_halves_direction() {
  {
    OracleEffect effect;
    effect.init();
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, false, 0.5f),
                 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, false, 1.5f),
                 -1.0f);
    HS_EXPECT_TRUE(effect.updateParameter("Opposite", 1.0f) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, true, 0.5f),
                 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, true, 1.5f),
                 1.0f);
  }
  Timeline().clear();

  OracleState state;
  state.shape = OracleEffect::ShapeType::STAR;
  state.function = OracleEffect::PhaseFunction::SAWTOOTH;
  state.count = 6;
  state.sides = 7;
  state.phase = 0.23f;
  const OracleFrame unchecked = capture_frame(state, candidate_renderer());
  state.opposite = true;
  const OracleFrame checked = capture_frame(state, candidate_renderer());
  HS_EXPECT_GT(compare_buffers(unchecked, checked).different_pixels,
               static_cast<size_t>(0));
}

inline void test_preset_transition_snaps() {
  {
    OracleEffect effect;
    effect.init();
    ShapeShifterWhiteBox::next_preset(effect);

    auto value = [&](const char *name) {
      for (const auto &def : effect.getParameters())
        if (std::strcmp(def.name, name) == 0)
          return def.get();
      HS_EXPECT(false, "ShapeShifter parameter is missing");
      return -1.0f;
    };

    HS_EXPECT_EQ(value("Alpha"), 0.5f);
    HS_EXPECT_EQ(value("Shape"), 2.793f);
    HS_EXPECT_EQ(value("Count"), 43.327999f);
    HS_EXPECT_EQ(value("Sides"), 6.562f);
    HS_EXPECT_EQ(value("Function"), 0.0f);
    HS_EXPECT_EQ(value("Amplitude"), 1.0f);
    HS_EXPECT_EQ(value("Speed"), 0.0142f);
    HS_EXPECT_EQ(value("Opposite"), 0.0f);
  }
  Timeline().clear();
}

inline int run_shapeshifter_oracle_tests() {
  ModuleFixture fixture("shapeshifter_oracle");
  test_buffer_comparator_statistics();
  test_reference_matrix_is_exact_and_nonblack();
  test_candidate_matrix_stays_within_visual_budget();
  test_segment_tiles_reconstruct_full_frame();
  test_amplitude_preserves_sweep_velocity();
  test_opposite_halves_direction();
  test_preset_transition_snaps();
  return fixture.result();
}

} // namespace shapeshifter_oracle_tests
} // namespace hs_test
