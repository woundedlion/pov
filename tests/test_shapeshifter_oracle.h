/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <type_traits>
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

/** Candidate-vs-reference visual budget. Every bound is the maximum measured
 * over the oracle matrix and the high-count star presets times a headroom
 * factor, and every case prints achieved-vs-budget. Whole-frame statistics
 * average over the entire canvas and take the tight factor; the worst single
 * pixel is a discrete extremum that hops to a neighbouring edge pixel when the
 * rounding shifts, so it takes the loose one. Mean, RMS and the star high-error
 * count peak on the square-wave planar star at count=75; the worst pixel and
 * the non-star high-error count on the square-wave flower at count=75. */
constexpr double WHOLE_FRAME_HEADROOM = 1.2;
constexpr double WORST_PIXEL_HEADROOM = 1.5;
constexpr double MEASURED_MEAN_ABSOLUTE_ERROR = 279.1;
constexpr double MEASURED_ROOT_MEAN_SQUARED_ERROR = 947.7;
constexpr uint32_t MEASURED_CHANNEL_ERROR = 15768;
constexpr double MEASURED_ENERGY_DRIFT = 0.00397;
constexpr double MEASURED_HIGH_COUNT_STAR_ENERGY_DRIFT = 0.0120;
constexpr size_t MEASURED_HIGH_ERROR_PIXELS = 452;
constexpr size_t MEASURED_STAR_HIGH_ERROR_PIXELS = 1456;

constexpr double MAX_MEAN_ABSOLUTE_ERROR =
    MEASURED_MEAN_ABSOLUTE_ERROR * WHOLE_FRAME_HEADROOM;
constexpr double MAX_ROOT_MEAN_SQUARED_ERROR =
    MEASURED_ROOT_MEAN_SQUARED_ERROR * WHOLE_FRAME_HEADROOM;
constexpr uint32_t MAX_CHANNEL_ERROR =
    static_cast<uint32_t>(MEASURED_CHANNEL_ERROR * WORST_PIXEL_HEADROOM);
constexpr double MAX_ENERGY_DRIFT =
    MEASURED_ENERGY_DRIFT * WHOLE_FRAME_HEADROOM;
constexpr double MAX_HIGH_COUNT_STAR_ENERGY_DRIFT =
    MEASURED_HIGH_COUNT_STAR_ENERGY_DRIFT * WHOLE_FRAME_HEADROOM;
constexpr size_t MAX_HIGH_ERROR_PIXELS =
    static_cast<size_t>(MEASURED_HIGH_ERROR_PIXELS * WHOLE_FRAME_HEADROOM);
constexpr size_t MAX_STAR_HIGH_ERROR_PIXELS =
    static_cast<size_t>(MEASURED_STAR_HIGH_ERROR_PIXELS * WHOLE_FRAME_HEADROOM);
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
  OracleEffect::AlphaFalloff alpha_falloff =
      OracleEffect::AlphaFalloff::TOWARD_EQUATOR;
  OracleEffect::RadiusSpacing spacing = OracleEffect::RadiusSpacing::UNIFORM;
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
  static constexpr size_t PRESET_COUNT = OracleEffect::PRESETS.size();
  static_assert(PRESET_COUNT == 9,
                "shipping preset roster changed; update the oracle's "
                "per-preset expectations");

  static void configure(OracleEffect &effect, const OracleState &state) {
    effect.alpha = state.alpha;
    effect.params.shape = state.shape;
    effect.params.count = static_cast<float>(state.count);
    effect.params.sides = static_cast<float>(state.sides);
    effect.params.function = state.function;
    effect.params.speed = 0.0f;
    effect.params.opposite = state.opposite;
    effect.params.alpha_falloff = state.alpha_falloff;
    effect.params.spacing = state.spacing;
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

  static float alpha_falloff_at(OracleEffect::AlphaFalloff falloff,
                                float radius_t, int count) {
    return OracleEffect::alpha_falloff_at(falloff, radius_t, count);
  }

  static float screen_balanced_radius_t(float radius_t) {
    return OracleEffect::screen_balanced_radius_t(radius_t);
  }

  static std::vector<int> folded_draw_order(int count) {
    std::vector<int> order;
    for (int ordinal = 0; ordinal < count; ++ordinal)
      order.push_back(OracleEffect::folded_draw_index(count, ordinal));
    return order;
  }

  static int shape_count() { return OracleEffect::NUM_SHAPES; }

  static const char *shape_option(int index) {
    return OracleEffect::SHAPE_OPTIONS[index];
  }

  static const char *function_option(int index) {
    return OracleEffect::FUNCTION_OPTIONS[index];
  }

  static const char *shape_export_option(int index) {
    return OracleEffect::SHAPE_EXPORT_OPTIONS[index];
  }

  static const char *spacing_option(int index) {
    return OracleEffect::SPACING_OPTIONS[index];
  }

  static const char *spacing_export_option(int index) {
    return OracleEffect::SPACING_EXPORT_OPTIONS[index];
  }

  static const OracleEffect::Params &preset_params(size_t index) {
    return OracleEffect::PRESETS[index].params;
  }

  static OracleEffect::ShapeType preset_shape(size_t index) {
    return preset_params(index).shape;
  }

  static OracleEffect::AlphaFalloff preset_alpha_falloff(size_t index) {
    return preset_params(index).alpha_falloff;
  }

  static OracleEffect::RadiusSpacing preset_spacing(size_t index) {
    return preset_params(index).spacing;
  }

  static void next_preset(OracleEffect &effect) { effect.next_preset(); }

  static void step_timeline(OracleEffect &effect, Canvas &canvas) {
    effect.timeline.step(canvas);
  }

  static float preset_opacity(const OracleEffect &effect) {
    return effect.preset_opacity;
  }

  static size_t preset_index(const OracleEffect &effect) {
    return effect.presets.current_index();
  }

  static int preset_frames() { return OracleEffect::PRESET_FRAMES; }
  static int preset_segue_frames() { return OracleEffect::PRESET_SEGUE_FRAMES; }
};

/** @brief Captures one renderer callable from a fresh deterministic effect. */
template <typename Render>
inline OracleFrame capture_frame(const OracleState &state, Render &&render) {
  reset_globals();
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

inline uint64_t row_energy(const OracleFrame &frame, int y) {
  uint64_t energy = 0;
  for (int x = 0; x < ORACLE_W; ++x) {
    const Pixel &pixel = frame.at(x, y);
    energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
  }
  return energy;
}

inline std::array<OracleState, 20> exhaustive_matrix() {
  using Function = OracleEffect::PhaseFunction;
  using Shape = OracleEffect::ShapeType;
  const Shape shapes[] = {Shape::PLANAR_POLYGON, Shape::SPHERICAL_POLYGON,
                          Shape::FLOWER, Shape::PLANAR_STAR,
                          Shape::SPHERICAL_STAR};
  const Function functions[] = {Function::SINE, Function::TRIANGLE,
                                Function::SAWTOOTH, Function::SQUARE};
  const int counts[] = {1, 2, 7, 75};
  const int sides[] = {3, 5, 9, 16, 7};
  const float phases[] = {0.0f, 0.249f, 0.5f, 0.999f};
  const Quaternion orientations[] = {
      Quaternion(), make_rotation(X_AXIS, Y_AXIS),
      make_rotation(X_AXIS, -Y_AXIS), make_rotation(X_AXIS, Z_AXIS)};

  std::array<OracleState, 20> matrix;
  size_t index = 0;
  for (int shape = 0; shape < 5; ++shape)
    for (int function = 0; function < 4; ++function) {
      OracleState state;
      state.shape = shapes[shape];
      const bool star = state.shape == Shape::PLANAR_STAR ||
                        state.shape == Shape::SPHERICAL_STAR;
      state.alpha_falloff = star ? OracleEffect::AlphaFalloff::TOWARD_EQUATOR
                                 : OracleEffect::AlphaFalloff::CONSTANT_HALF;
      state.function = functions[function];
      state.count = state.shape == Shape::SPHERICAL_STAR
                        ? std::min(counts[function], 7)
                        : counts[function];
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

/**
 * @brief Verifies the reference render is reproducible and non-black over the
 *        whole matrix.
 * @details Both sides of the comparison are reference_renderer(), so this
 *          establishes determinism, not correctness: the capture harness
 *          reproduces a state bit-identically, and the reference frame carries
 *          energy. Those are the preconditions the candidate-vs-reference
 *          budget cases rest on — a nondeterministic harness or an all-black
 *          reference would make every error bound below vacuous.
 */
inline void test_reference_matrix_is_deterministic_and_nonblack() {
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

/** @brief Prints one case's achieved statistics against its budget. */
inline void report_visual_budget(const OracleState &state,
                                 const FrameErrorStats &error,
                                 double energy_ratio, size_t high_error_pixels,
                                 double max_energy_drift,
                                 size_t max_high_error_pixels) {
  std::printf(
      "  [%s %s count=%d sides=%d] MAE %.1f/%.1f; RMSE %.1f/%.1f; max channel "
      "%u/%u; energy drift %.3f%%/%.3f%%; high-error px %zu/%zu\n",
      ShapeShifterWhiteBox::shape_option(static_cast<int>(state.shape)),
      ShapeShifterWhiteBox::function_option(static_cast<int>(state.function)),
      state.count, state.sides, error.mean_absolute_error(),
      MAX_MEAN_ABSOLUTE_ERROR, error.root_mean_squared_error(),
      MAX_ROOT_MEAN_SQUARED_ERROR, error.max_absolute_error, MAX_CHANNEL_ERROR,
      100.0 * energy_ratio, 100.0 * max_energy_drift, high_error_pixels,
      max_high_error_pixels);
}

inline void expect_candidate_within_visual_budget(
    const OracleState &state, double max_energy_drift = MAX_ENERGY_DRIFT,
    size_t max_high_error_pixels = MAX_HIGH_ERROR_PIXELS) {
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
  for (size_t pixel = 0; pixel < comparison.reference.pixels.size(); ++pixel) {
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
  report_visual_budget(state, comparison.error, energy_ratio, high_error_pixels,
                       max_energy_drift, max_high_error_pixels);
  HS_EXPECT_LT(comparison.error.mean_absolute_error(), MAX_MEAN_ABSOLUTE_ERROR);
  HS_EXPECT_LT(comparison.error.root_mean_squared_error(),
               MAX_ROOT_MEAN_SQUARED_ERROR);
  HS_EXPECT_LT(comparison.error.max_absolute_error, MAX_CHANNEL_ERROR);
  HS_EXPECT_LT(energy_ratio, max_energy_drift);
  HS_EXPECT_LT(high_error_pixels, max_high_error_pixels);
  HS_EXPECT_EQ(uncovered_bright_pixels, size_t{0});
  HS_EXPECT_EQ(comparison.reference.backstop_hits, uint32_t{0});
  HS_EXPECT_EQ(comparison.candidate.backstop_hits, uint32_t{0});
}

inline void test_candidate_matrix_stays_within_visual_budget() {
  for (const OracleState &state : exhaustive_matrix()) {
    const size_t max_high_error_pixels =
        state.shape == OracleEffect::ShapeType::PLANAR_STAR ||
                state.shape == OracleEffect::ShapeType::SPHERICAL_STAR
            ? MAX_STAR_HIGH_ERROR_PIXELS
            : MAX_HIGH_ERROR_PIXELS;
    expect_candidate_within_visual_budget(state, MAX_ENERGY_DRIFT,
                                          max_high_error_pixels);
  }
}

inline void test_star_projection_policies_render_different_edges() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SQUARE;
  state.count = 3;
  state.sides = 5;
  state.phase = 0.17f;
  state.orientation = Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
  const OracleFrame planar = capture_frame(state, candidate_renderer());
  state.shape = OracleEffect::ShapeType::SPHERICAL_STAR;
  const OracleFrame spherical = capture_frame(state, candidate_renderer());
  HS_EXPECT_GT(compare_buffers(planar, spherical).different_pixels,
               size_t{100});
}

inline void test_screen_balanced_spacing_is_opt_in_for_every_shape() {
  using Shape = OracleEffect::ShapeType;
  constexpr std::array<Shape, 5> SHAPES = {
      Shape::PLANAR_POLYGON, Shape::SPHERICAL_POLYGON, Shape::FLOWER,
      Shape::PLANAR_STAR, Shape::SPHERICAL_STAR};
  for (Shape shape : SHAPES) {
    OracleState state;
    state.shape = shape;
    state.function = OracleEffect::PhaseFunction::TRIANGLE;
    state.count = 17;
    state.sides = 7;
    state.phase = 0.23f;
    state.orientation = Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
    const OracleFrame uniform = capture_frame(state, candidate_renderer());
    state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
    const OracleFrame balanced = capture_frame(state, candidate_renderer());
    HS_EXPECT_GT(compare_buffers(uniform, balanced).different_pixels,
                 size_t{100});
  }
}

inline void test_star_options_and_shipping_presets_are_planar() {
  using Falloff = OracleEffect::AlphaFalloff;
  using Shape = OracleEffect::ShapeType;
  using Spacing = OracleEffect::RadiusSpacing;
  HS_EXPECT_EQ(ShapeShifterWhiteBox::shape_count(), 5);
  HS_EXPECT_TRUE(
      std::strcmp(ShapeShifterWhiteBox::shape_option(3), "Planar Star") == 0);
  HS_EXPECT_TRUE(std::strcmp(ShapeShifterWhiteBox::shape_option(4),
                             "Spherical Star") == 0);
  HS_EXPECT_TRUE(std::strcmp(ShapeShifterWhiteBox::shape_export_option(3),
                             "ShapeType::PLANAR_STAR") == 0);
  HS_EXPECT_TRUE(std::strcmp(ShapeShifterWhiteBox::shape_export_option(4),
                             "ShapeType::SPHERICAL_STAR") == 0);
  HS_EXPECT_TRUE(std::strcmp(ShapeShifterWhiteBox::spacing_option(1),
                             "Screen Balanced") == 0);
  HS_EXPECT_TRUE(std::strcmp(ShapeShifterWhiteBox::spacing_export_option(1),
                             "RadiusSpacing::SCREEN_BALANCED") == 0);
  for (size_t index : {size_t{0}, size_t{2}, size_t{4}}) {
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_shape(index), Shape::PLANAR_STAR);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_alpha_falloff(index),
                 Falloff::TOWARD_EQUATOR);
  }
  HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_spacing(0),
               Spacing::SCREEN_BALANCED);
  for (size_t index = 1; index < ShapeShifterWhiteBox::PRESET_COUNT; ++index)
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_spacing(index), Spacing::UNIFORM);
  for (size_t index :
       {size_t{1}, size_t{3}, size_t{5}, size_t{6}, size_t{7}, size_t{8}})
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_alpha_falloff(index),
                 Falloff::CONSTANT_HALF);
}

inline void test_high_count_star_preset_stays_within_visual_budget() {
  using Function = OracleEffect::PhaseFunction;
  using Shape = OracleEffect::ShapeType;
  const std::array<Quaternion, 6> orientations = {{
      Quaternion(),
      make_rotation(X_AXIS, Y_AXIS),
      make_rotation(X_AXIS, -Y_AXIS),
      make_rotation(X_AXIS, Z_AXIS),
      Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(),
      Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
  }};
  const std::array<float, 6> phases = {
      {0.0f, 0.125f, 0.249f, 0.5f, 0.75f, 0.999f}};
  for (size_t i = 0; i < orientations.size(); ++i) {
    OracleState state;
    state.shape = Shape::PLANAR_STAR;
    state.function = Function::SINE;
    state.count = 144;
    state.sides = 7;
    state.phase = phases[i];
    state.alpha = 0.274f;
    state.orientation = orientations[i];
    expect_candidate_within_visual_budget(
        state, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT, MAX_STAR_HIGH_ERROR_PIXELS);
  }
}

inline void test_screen_balanced_star_preset_stays_within_visual_budget() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.alpha_falloff = OracleEffect::AlphaFalloff::TOWARD_EQUATOR;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  state.count = 208;
  state.sides = 7;
  state.phase = 0.249f;
  state.alpha = 0.274f;
  state.orientation = Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized();
  expect_candidate_within_visual_budget(state, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT,
                                        MAX_STAR_HIGH_ERROR_PIXELS);
}

inline void test_high_count_star_preset_covers_north_pole() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.count = 144;
  state.sides = 7;
  state.phase = 0.125f;
  state.alpha = 0.274f;
  state.orientation = make_rotation(X_AXIS, Y_AXIS);
  const RenderComparison comparison =
      compare_renders(state, reference_renderer(), candidate_renderer());
  const uint64_t reference_pole_energy = row_energy(comparison.reference, 0);
  const uint64_t candidate_pole_energy = row_energy(comparison.candidate, 0);
  size_t uncovered_pole_pixels = 0;
  for (int x = 0; x < ORACLE_W; ++x) {
    if (!pixel_is_bright(comparison.reference.at(x, 0)))
      continue;
    bool covered = false;
    for (int dx = -1; dx <= 1; ++dx) {
      const int candidate_x = (x + dx + ORACLE_W) % ORACLE_W;
      covered |= pixel_has_coverage(comparison.candidate.at(candidate_x, 0));
    }
    if (!covered)
      ++uncovered_pole_pixels;
  }
  HS_EXPECT_GT(reference_pole_energy, uint64_t{0});
  HS_EXPECT_GE(candidate_pole_energy * 10, reference_pole_energy * 9);
  HS_EXPECT_EQ(uncovered_pole_pixels, size_t{0});
}

inline uint32_t minimum_energy_around(const OracleFrame &frame,
                                      const Vector &position, float radius_px) {
  const PixelCoords center = vector_to_pixel<ORACLE_W, ORACLE_H>(position);
  uint32_t minimum = UINT32_MAX;
  for (int y = 0; y < ORACLE_H; ++y) {
    for (int x = 0; x < ORACLE_W; ++x) {
      float dx = std::fabs(static_cast<float>(x) - center.x);
      dx = std::min(dx, static_cast<float>(ORACLE_W) - dx);
      const float dy = static_cast<float>(y) - center.y;
      if (dx * dx + dy * dy > radius_px * radius_px)
        continue;
      const Pixel &pixel = frame.at(x, y);
      minimum =
          std::min(minimum, static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b);
    }
  }
  return minimum;
}

inline void test_high_count_planar_star_caps_cover_chart_centers() {
  constexpr uint32_t MIN_CAP_ENERGY = 60000;
  const std::array<Quaternion, 4> orientations = {{
      Quaternion(),
      Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(),
      Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
      Quaternion(0.51f, 0.63f, -0.28f, 0.51f).normalized(),
  }};
  for (const auto [count, spacing] :
       {std::pair{144, OracleEffect::RadiusSpacing::UNIFORM},
        std::pair{208, OracleEffect::RadiusSpacing::SCREEN_BALANCED}}) {
    for (const Quaternion &orientation : orientations) {
      OracleState state;
      state.shape = OracleEffect::ShapeType::PLANAR_STAR;
      state.function = OracleEffect::PhaseFunction::SINE;
      state.spacing = spacing;
      state.count = count;
      state.sides = 7;
      state.phase = 0.125f;
      state.alpha = 0.274f;
      state.orientation = orientation;
      const OracleFrame frame = capture_frame(state, candidate_renderer());
      const Vector near = rotate(X_AXIS, orientation);
      const Vector far = -near;
      HS_EXPECT_GE(minimum_energy_around(frame, near, 1.5f), MIN_CAP_ENERGY);
      HS_EXPECT_GE(minimum_energy_around(frame, far, 1.5f), MIN_CAP_ENERGY);
    }
  }
}

inline int first_covered_north_row(const OracleFrame &frame) {
  for (int y = 0; y < 32; ++y)
    for (int x = 0; x < ORACLE_W; ++x) {
      const Pixel &pixel = frame.at(x, y);
      if (static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > 0)
        return y;
    }
  return 32;
}

inline void test_high_count_spherical_star_contours_reach_display_north() {
  for (float phase : {0.0f, 0.125f, 0.249f, 0.5f, 0.75f, 0.999f}) {
    OracleState state;
    state.shape = OracleEffect::ShapeType::SPHERICAL_STAR;
    state.function = OracleEffect::PhaseFunction::SINE;
    state.count = 144;
    state.sides = 7;
    state.phase = phase;
    state.alpha = 0.274f;
    state.orientation = Quaternion();
    const OracleFrame frame = capture_frame(state, candidate_renderer());
    HS_EXPECT_LE(first_covered_north_row(frame), 1);
  }
}

inline void test_amplitude_preserves_sweep_velocity() {
  OracleEffect effect;
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 1.0f),
                 0.04f, 1e-6f);
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 4.0f),
                 0.01f, 1e-6f);
}

/** @brief Contour midpoint the render path samples the baked palette at. */
inline float midpoint_radius_t(int index, int count) {
  return (static_cast<float>(index) + 0.5f) / static_cast<float>(count);
}

/** @brief Discrete per-contour alpha ramp, in index space. */
inline float midpoint_alpha_reference(int index, int count) {
  const int steps_to_equator = (count - 1) / 2;
  if (steps_to_equator == 0)
    return 1.0f;
  const int opposite_index = count - index - 1;
  const int distance_from_pole =
      index < opposite_index ? index : opposite_index;
  const float equator_alpha = 2.0f / static_cast<float>(count);
  return 1.0f - (1.0f - equator_alpha) *
                    static_cast<float>(distance_from_pole) /
                    static_cast<float>(steps_to_equator);
}

inline void test_alpha_falloff_fades_to_equator() {
  using Falloff = OracleEffect::AlphaFalloff;
  auto toward_equator = [](int index, int count) {
    return ShapeShifterWhiteBox::alpha_falloff_at(
        Falloff::TOWARD_EQUATOR, midpoint_radius_t(index, count), count);
  };
  HS_EXPECT_EQ(toward_equator(0, 1), 1.0f);
  HS_EXPECT_NEAR(toward_equator(0, 6), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(1, 6), 2.0f / 3.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(2, 6), 1.0f / 3.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(3, 6), 1.0f / 3.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(4, 6), 2.0f / 3.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(5, 6), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(toward_equator(2, 5), 2.0f / 5.0f, 1e-6f);
  for (int index = 0; index < 6; ++index)
    HS_EXPECT_EQ(ShapeShifterWhiteBox::alpha_falloff_at(
                     Falloff::CONSTANT_HALF, midpoint_radius_t(index, 6), 6),
                 0.5f);

  for (int count = 1; count <= OracleEffect::MAX_SHAPES; ++count) {
    for (int index = 0; index < count; ++index) {
      const float radius_t = midpoint_radius_t(index, count);
      HS_EXPECT_NEAR(ShapeShifterWhiteBox::alpha_falloff_at(
                         Falloff::TOWARD_EQUATOR, radius_t, count),
                     midpoint_alpha_reference(index, count), 1e-6f);
      HS_EXPECT_EQ(ShapeShifterWhiteBox::alpha_falloff_at(
                       Falloff::CONSTANT_HALF, radius_t, count),
                   0.5f);
    }
  }
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
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
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

inline void test_folded_shapes_draw_from_equator_to_both_poles() {
  HS_EXPECT_EQ(ShapeShifterWhiteBox::folded_draw_order(6),
               (std::vector<int>{3, 2, 4, 1, 5, 0}));
  HS_EXPECT_EQ(ShapeShifterWhiteBox::folded_draw_order(5),
               (std::vector<int>{2, 3, 1, 4, 0}));
}

inline void test_screen_balanced_spacing_follows_sampling_envelope() {
  constexpr int COUNT = 208;
  constexpr float DENSITY_FLOOR = 0.5f;
  constexpr float DENSITY_INTEGRAL = 1.1278247916f;
  auto mapped = [](int index) {
    const float radius_t =
        (static_cast<float>(index) + 0.5f) / static_cast<float>(COUNT);
    return ShapeShifterWhiteBox::screen_balanced_radius_t(radius_t);
  };

  const float equator_equivalent_count =
      1.0f / (mapped(COUNT / 2 - 1) - mapped(COUNT / 2 - 2));
  const float pole_equivalent_count = 1.0f / (mapped(1) - mapped(0));
  const float density_scale = COUNT * PI_F / (2.0f * DENSITY_INTEGRAL);
  HS_EXPECT_NEAR(equator_equivalent_count, density_scale, 0.2f);
  HS_EXPECT_NEAR(equator_equivalent_count, OracleEffect::MAX_SHAPES, 1.8f);
  HS_EXPECT_NEAR(pole_equivalent_count, density_scale * DENSITY_FLOOR, 0.2f);
  HS_EXPECT_NEAR(pole_equivalent_count, ORACLE_H, 1.0f);
  HS_EXPECT_NEAR(mapped(0) * PI_F, 2.0f * DENSITY_INTEGRAL / COUNT, 1e-6f);

  for (int index = 0; index < COUNT; ++index) {
    const float radius_t =
        (static_cast<float>(index) + 0.5f) / static_cast<float>(COUNT);
    const float opposite_radius_t = 1.0f - radius_t;
    HS_EXPECT_NEAR(
        ShapeShifterWhiteBox::screen_balanced_radius_t(radius_t) +
            ShapeShifterWhiteBox::screen_balanced_radius_t(opposite_radius_t),
        1.0f, 1e-6f);
  }
}

inline void test_preset_transition_snaps() {
  {
    OracleEffect effect;
    effect.init();
    HS_EXPECT_TRUE(effect.updateParameter("Alpha", 0.42f) ==
                   ParamSetResult::APPLIED);
    ShapeShifterWhiteBox::next_preset(effect);

    auto value = [&](const char *name) {
      for (const auto &def : effect.getParameters())
        if (std::strcmp(def.name, name) == 0)
          return def.get();
      HS_EXPECT(false, "ShapeShifter parameter is missing");
      return -1.0f;
    };

    auto option_index = [](auto field) {
      return static_cast<float>(
          static_cast<std::underlying_type_t<decltype(field)>>(field));
    };
    const auto &preset = ShapeShifterWhiteBox::preset_params(1);

    HS_EXPECT_EQ(value("Alpha"), 0.42f);
    HS_EXPECT_EQ(value("Shape"), option_index(preset.shape));
    HS_EXPECT_EQ(value("Count"), preset.count);
    HS_EXPECT_EQ(value("Spacing"), option_index(preset.spacing));
    HS_EXPECT_EQ(value("Sides"), preset.sides);
    HS_EXPECT_EQ(value("Function"), option_index(preset.function));
    HS_EXPECT_EQ(value("Amplitude"), preset.amplitude);
    HS_EXPECT_EQ(value("Speed"), preset.speed);
    HS_EXPECT_EQ(value("Opposite"), preset.opposite ? 1.0f : 0.0f);
    HS_EXPECT_EQ(value("Alpha Falloff"), option_index(preset.alpha_falloff));
  }
  Timeline().clear();
}

/** @brief Pins ShapeShifter's eight-frame fade-out and fade-in preset seam. */
inline void test_preset_transition_fades_through_black_in_16_frames() {
  {
    OracleEffect effect;
    effect.init();
    Canvas canvas(effect);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_segue_frames(), 16);

    for (int frame = 0; frame < ShapeShifterWhiteBox::preset_frames() - 8;
         ++frame)
      ShapeShifterWhiteBox::step_timeline(effect, canvas);
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 1.0f, 1e-6f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{0});

    for (int frame = 0; frame < 7; ++frame)
      ShapeShifterWhiteBox::step_timeline(effect, canvas);
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 0.125f, 1e-6f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{0});

    ShapeShifterWhiteBox::step_timeline(effect, canvas);
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 0.0f, 1e-6f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});

    ShapeShifterWhiteBox::step_timeline(effect, canvas);
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 0.125f, 1e-6f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});

    for (int frame = 0; frame < 7; ++frame)
      ShapeShifterWhiteBox::step_timeline(effect, canvas);
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 1.0f, 1e-6f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});
  }
  Timeline().clear();
}

inline int run_shapeshifter_oracle_tests() {
  ModuleFixture fixture("shapeshifter_oracle");
  test_buffer_comparator_statistics();
  test_reference_matrix_is_deterministic_and_nonblack();
  test_candidate_matrix_stays_within_visual_budget();
  test_star_projection_policies_render_different_edges();
  test_screen_balanced_spacing_is_opt_in_for_every_shape();
  test_star_options_and_shipping_presets_are_planar();
  test_high_count_star_preset_stays_within_visual_budget();
  test_screen_balanced_star_preset_stays_within_visual_budget();
  test_high_count_star_preset_covers_north_pole();
  test_high_count_planar_star_caps_cover_chart_centers();
  test_high_count_spherical_star_contours_reach_display_north();
  test_amplitude_preserves_sweep_velocity();
  test_alpha_falloff_fades_to_equator();
  test_opposite_halves_direction();
  test_folded_shapes_draw_from_equator_to_both_poles();
  test_screen_balanced_spacing_follows_sampling_envelope();
  test_preset_transition_snaps();
  test_preset_transition_fades_through_black_in_16_frames();
  return fixture.result();
}

} // namespace shapeshifter_oracle_tests
} // namespace hs_test
