/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Shape rendering comparisons and visual budgets.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "effects/ShapeShifter.h"
#include "tests/pixel_test_util.h"
#include "tests/test_radial_surface.h"
#include "tests/test_radial_surface_mask.h"
#include "tests/temporal_oracle.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::shapeshifter_oracle_tests {

constexpr int ORACLE_W = 288;
constexpr int ORACLE_H = 144;
constexpr uint32_t BRIGHT_ENERGY = 12288;
constexpr uint32_t COVERAGE_ENERGY = 512;
constexpr uint32_t HIGH_CHANNEL_ERROR = 4096;

/** Measured oracle-matrix energy drift and high-count star error baselines,
 * with headroom multipliers. */
constexpr double WHOLE_FRAME_HEADROOM = 1.2;
constexpr double WORST_PIXEL_HEADROOM = 1.5;
constexpr double MEASURED_MEAN_ABSOLUTE_ERROR = 187.5;
constexpr double MEASURED_ROOT_MEAN_SQUARED_ERROR = 529.5;
constexpr uint32_t MEASURED_CHANNEL_ERROR = 10482;
constexpr double MEASURED_ENERGY_DRIFT = 0.00415;
constexpr double MEASURED_PRESET_ENERGY_DRIFT = 0.03337;
constexpr double MEASURED_HIGH_COUNT_STAR_ENERGY_DRIFT = 0.01179;
constexpr size_t MEASURED_STAR_HIGH_ERROR_PIXELS = 236;

constexpr double MAX_MEAN_ABSOLUTE_ERROR =
    MEASURED_MEAN_ABSOLUTE_ERROR * WHOLE_FRAME_HEADROOM;
constexpr double MAX_ROOT_MEAN_SQUARED_ERROR =
    MEASURED_ROOT_MEAN_SQUARED_ERROR * WHOLE_FRAME_HEADROOM;
constexpr uint32_t MAX_CHANNEL_ERROR =
    static_cast<uint32_t>(MEASURED_CHANNEL_ERROR * WORST_PIXEL_HEADROOM);
constexpr double MAX_ENERGY_DRIFT =
    MEASURED_ENERGY_DRIFT * WHOLE_FRAME_HEADROOM;
constexpr double MAX_PRESET_ENERGY_DRIFT =
    MEASURED_PRESET_ENERGY_DRIFT * WHOLE_FRAME_HEADROOM;
constexpr double MAX_HIGH_COUNT_STAR_ENERGY_DRIFT =
    MEASURED_HIGH_COUNT_STAR_ENERGY_DRIFT * WHOLE_FRAME_HEADROOM;
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
  size_t preset_index = 0;
  bool force_contours = false;
  int count = 7;
  int sides = 5;
  float phase = 0.0f;
  float amplitude = 1.0f;
  float alpha = 0.625f;
  bool opposite = false;
  math::Quaternion orientation;
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

/** @brief Test-only access to ShapeShifter's pinned state and shipping renderer. */
struct ShapeShifterWhiteBox {
  /** @brief Kernel family used by footprint oracles. */
  using Surface = OracleEffect::Surface::Family;
  /** @brief Contour tables of the oracle effect. */
  using Layout = Plot::ContourLayout<OracleEffect::MAX_SHAPES>;
  static constexpr size_t PRESET_COUNT = OracleEffect::PRESETS.size();
  static constexpr float NYQUIST_PHASE_STEP = OracleEffect::NYQUIST_PHASE_STEP;
  static constexpr float AMPLITUDE_MAX = OracleEffect::AMPLITUDE_MAX;
  static constexpr size_t SURFACE_SCRATCH_A_BYTES =
      OracleEffect::SURFACE_SCRATCH_A_BYTES;
  static constexpr size_t SURFACE_SCRATCH_B_BYTES =
      OracleEffect::SURFACE_SCRATCH_B_BYTES;

  static_assert(PRESET_COUNT >= 2,
                "the preset transition test steps from preset 0 to 1");

  static void prepare_count(OracleEffect &effect, int count) {
    effect.prepare_count(count);
  }

  static void force_contours(OracleEffect &effect, bool enabled) {
    effect.test_force_contours = enabled;
  }

  static void configure(OracleEffect &effect, const OracleState &state) {
    force_contours(effect, state.force_contours);
    if (effect.getPresetIndex() != state.preset_index)
      effect.profile_select_preset(state.preset_index);
    effect.alpha = state.alpha;
    effect.params.shape = state.shape;
    effect.params.count = static_cast<float>(state.count);
    effect.params.sides = static_cast<float>(state.sides);
    effect.params.function = state.function;
    effect.params.amplitude = state.amplitude;
    effect.params.speed = 0.0f;
    effect.params.opposite = state.opposite;
    effect.params.alpha_falloff = state.alpha_falloff;
    effect.params.spacing = state.spacing;
    effect.phase = state.phase;
    effect.orientation.set(state.orientation);
    effect.set_clip(state.clip.y0, state.clip.y1, state.clip.x0, state.clip.x1);
  }

  template <typename F>
  static void dispatch_plot_reference(OracleEffect &effect, Canvas &canvas,
                                      const math::Basis &basis,
                                      OracleEffect::ShapeType shape,
                                      float radius, int sides,
                                      const F &fragment_shader,
                                      float shape_phase) {
    using ShapeType = OracleEffect::ShapeType;
    PipelineRef pipeline(effect.plot_filters, canvas);
    switch (shape) {
    case ShapeType::PLANAR_POLYGON:
      Plot::Polygon<Plot::PlanarProjection>::draw<ORACLE_W, ORACLE_H>(
          pipeline, canvas, basis, radius, sides, fragment_shader, shape_phase);
      break;
    case ShapeType::SPHERICAL_POLYGON:
      Plot::Polygon<Plot::GeodesicProjection>::draw<ORACLE_W, ORACLE_H>(
          pipeline, canvas, basis, radius, sides, fragment_shader, shape_phase);
      break;
    case ShapeType::FLOWER:
      Plot::Flower::draw<ORACLE_W, ORACLE_H>(pipeline, canvas, basis, radius,
                                             sides, fragment_shader, {},
                                             shape_phase);
      break;
    case ShapeType::PLANAR_STAR:
      Plot::Star<Plot::PlanarProjection>::draw<ORACLE_W, ORACLE_H>(
          pipeline, canvas, basis, radius, sides, fragment_shader, shape_phase);
      break;
    case ShapeType::SPHERICAL_STAR:
      Plot::Star<Plot::GeodesicProjection>::draw_continuous<ORACLE_W, ORACLE_H>(
          pipeline, canvas, basis, radius, sides, fragment_shader, shape_phase);
      break;
    }
  }

  static void draw_all_reference(OracleEffect &effect, Canvas &canvas) {
    using ShapeType = OracleEffect::ShapeType;
    const int count = hs::clamp(static_cast<int>(effect.params.count), 1,
                                OracleEffect::DRAW_LIMIT);
    if (count != effect.prepared_count ||
        effect.params.spacing != effect.prepared_spacing)
      effect.prepare_count(count);
    const auto &palette = effect.selected_palette();
    const int sides = hs::clamp(static_cast<int>(effect.params.sides),
                                static_cast<int>(OracleEffect::SIDES_MIN),
                                static_cast<int>(OracleEffect::SIDES_MAX));
    const ShapeType shape = effect.selected_shape();
    const auto function = effect.selected_function();
    const math::Basis basis =
        math::make_basis(effect.orientation.get(), math::X_AXIS);
    const float global_alpha = effect.alpha * effect.preset_opacity;

    const bool continuous_star = shape == ShapeType::SPHERICAL_STAR;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const int i = continuous_star ? count - ordinal - 1
                                    : effect.layout.draw_index(ordinal);
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      const float radius = 2.0f * effect.layout.spaced_radius_t(i);
      const float direction =
          continuous_star
              ? Layout::phase_sign_star(effect.params.opposite, radius)
              : Layout::phase_sign_contour(effect.params.opposite, radius);
      const float contour_phase =
          direction * effect.params.amplitude *
          effect.evaluate(function, radius_t + effect.phase);
      const Color4 color = palette.get(radius_t);
      auto shader = [&](const math::Vector &, Fragment &fragment) {
        fragment.color = color;
        fragment.color.alpha *= global_alpha;
      };
      dispatch_plot_reference(effect, canvas, basis, shape, radius, sides,
                              shader, contour_phase);
    }
    if (shape == ShapeType::PLANAR_STAR)
      effect.draw_planar_star_pole_caps(canvas, basis, count, sides, palette);
  }

  static void render_reference(OracleEffect &effect, Canvas &canvas) {
    effect.plot_filters.prepare(canvas);
    draw_all_reference(effect, canvas);
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

  static float phase_sign_contour(const OracleEffect &effect, float radius) {
    return Layout::phase_sign_contour(effect.params.opposite, radius);
  }

  static float phase_sign_star(const OracleEffect &effect, float radius) {
    return Layout::phase_sign_star(effect.params.opposite, radius);
  }

  static float alpha_falloff_at(OracleEffect::AlphaFalloff falloff,
                                float radius_t, int count) {
    return OracleEffect::alpha_falloff_at(falloff, radius_t, count);
  }

  static float screen_balanced_radius_t(float radius_t) {
    return Plot::screen_balanced_radius_t(radius_t);
  }

  static std::vector<int> folded_draw_order(int count) {
    std::vector<int> order;
    for (int ordinal = 0; ordinal < count; ++ordinal)
      order.push_back(Layout::folded_draw_index(count, ordinal));
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

  /** @brief Steps the timeline and the preset choreography one frame. */
  static void step_choreography(OracleEffect &effect, Canvas &canvas) {
    effect.timeline.step(canvas);
    effect.step_choreography();
  }

  static OracleEffect::ShapeType shape(const OracleEffect &effect) {
    return effect.params.shape;
  }

  static float preset_opacity(const OracleEffect &effect) {
    return effect.preset_opacity;
  }

  static size_t preset_index(const OracleEffect &effect) {
    return effect.getPresetIndex();
  }

  /** @brief The effect's orientation.
   * @param effect Effect.
   * @return Its orientation. */
  static math::Orientation<> &orientation(OracleEffect &effect) {
    return effect.orientation;
  }

  /** @brief Sets the waveform phase and its advance over the frame.
   * @param effect Effect.
   * @param phase Phase at the frame's end, turns.
   * @param step Advance over the frame, turns. */
  static void set_twist(OracleEffect &effect, float phase, float step) {
    effect.phase = phase;
    effect.twist_step = step;
  }

  static void replay_step(OracleEffect &e, Canvas &canvas) {
    e.timeline.step(canvas);
    e.step_choreography();
    e.advance_phase();
  }
  static float phase(const OracleEffect &e) { return e.phase; }
  static float twist_step(const OracleEffect &e) { return e.twist_step; }
  static void surface_kernel(OracleEffect &e, Canvas &canvas, bool timed,
                             Scan::SurfaceAreaStats *stats = nullptr,
                             float displacement = 0, bool mask = true) {
    e.prepare_count(288);
    const Scan::SurfaceFrame frame = e.surface_frame(e.selected_palette());
    if (mask)
      OracleEffect::Surface::draw<true, true>(canvas, frame, timed, stats,
                                              displacement);
    else
      OracleEffect::Surface::draw<true, false>(canvas, frame, timed, stats,
                                               displacement);
  }
  static int preset_frames() { return OracleEffect::PRESET_FRAMES; }
  static constexpr int preset_dwell_frames() {
    return OracleEffect::PRESET_DWELL_FRAMES;
  }
  static constexpr int preset_segue_frames() {
    return OracleEffect::DEPARTURE.frames;
  }
};

inline OracleState preset_state(size_t i) {
  const auto &preset = ShapeShifterWhiteBox::preset_params(i);
  OracleState state;
  state.shape = preset.shape;
  state.function = preset.function;
  state.alpha_falloff = preset.alpha_falloff;
  state.spacing = preset.spacing;
  state.count = static_cast<int>(preset.count);
  state.sides = static_cast<int>(preset.sides);
  state.amplitude = preset.amplitude;
  state.opposite = preset.opposite;
  state.phase = 0.37f;
  state.alpha = 0.274f;
  return state;
}

inline float preset_phase_step(size_t i) {
  const auto &preset = ShapeShifterWhiteBox::preset_params(i);
  return fminf(preset.speed / preset.amplitude,
               ShapeShifterWhiteBox::NYQUIST_PHASE_STEP);
}

inline OracleState dense_star_state() {
  OracleState state;
  state.preset_index = 1;
  state.force_contours = true;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.alpha_falloff = OracleEffect::AlphaFalloff::TOWARD_EQUATOR;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  state.count = 288;
  state.sides = 7;
  state.phase = 0.249f;
  state.alpha = 0.274f;
  state.orientation =
      math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized();
  return state;
}

/** @brief Captures one renderer callable from a fresh deterministic effect.
 * @details Temporal sampling is off: these oracles judge the shapes' still
 * geometry, which the timed route draws at each column's time. */
template <typename Render>
inline OracleFrame render_oracle_frame(const OracleState &state,
                                       Render &&render) {
  reset_globals();
  const ScopedTemporalSampling sampling(false);
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
    hs_test::capture_frame<ORACLE_W, ORACLE_H>(effect, frame.pixels);
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
  comparison.reference = render_oracle_frame(
      state, std::forward<ReferenceRender>(reference_render));
  comparison.candidate = render_oracle_frame(
      state, std::forward<CandidateRender>(candidate_render));
  comparison.error =
      compare_buffers(comparison.reference.pixels, comparison.candidate.pixels);
  return comparison;
}

/** @brief Renders the reference shape path. */
inline auto reference_renderer() {
  return [](OracleEffect &effect, Canvas &canvas) {
    ShapeShifterWhiteBox::render_reference(effect, canvas);
  };
}

/** @brief Renders the candidate shape path. */
inline auto candidate_renderer() {
  return [](OracleEffect &effect, Canvas &canvas) {
    ShapeShifterWhiteBox::render_candidate(effect, canvas);
  };
}

/** @brief Sums linear RGB channels across the frame.
 * @param frame Rendered oracle frame.
 * @return Total linear RGB energy. */
inline uint64_t frame_energy(const OracleFrame &frame) {
  return hs_test::frame_energy(frame.pixels);
}

/** @brief Sums linear RGB channels along one frame row.
 * @param frame Rendered oracle frame.
 * @param y Row index.
 * @return Total linear RGB energy on the row. */
inline uint64_t row_energy(const OracleFrame &frame, int y) {
  uint64_t energy = 0;
  for (int x = 0; x < ORACLE_W; ++x) {
    const Pixel &pixel = frame.at(x, y);
    energy += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
  }
  return energy;
}

/**
 * @brief One state per (shape, phase function) pair, shape-major: entry
 *        shape * 4 + function.
 */
inline std::array<OracleState, 20> shape_function_matrix() {
  using Function = OracleEffect::PhaseFunction;
  using Shape = OracleEffect::ShapeType;
  const Shape SHAPES[] = {Shape::PLANAR_POLYGON, Shape::SPHERICAL_POLYGON,
                          Shape::FLOWER, Shape::PLANAR_STAR,
                          Shape::SPHERICAL_STAR};
  const Function FUNCTIONS[] = {Function::SINE, Function::TRIANGLE,
                                Function::SAWTOOTH, Function::SQUARE};
  const int COUNTS[] = {1, 2, 7, 75};
  const int SIDES[] = {3, 5, 9, 16, 7};
  const float PHASES[] = {0.0f, 0.249f, 0.5f, 0.999f};
  const math::Quaternion orientations[] = {
      math::Quaternion(), math::make_rotation(math::X_AXIS, math::Y_AXIS),
      math::make_rotation(math::X_AXIS, -math::Y_AXIS),
      math::make_rotation(math::X_AXIS, math::Z_AXIS)};

  std::array<OracleState, 20> matrix;
  size_t index = 0;
  for (int shape = 0; shape < 5; ++shape)
    for (int function = 0; function < 4; ++function) {
      OracleState state;
      state.shape = SHAPES[shape];
      const bool star = state.shape == Shape::PLANAR_STAR ||
                        state.shape == Shape::SPHERICAL_STAR;
      state.alpha_falloff = star ? OracleEffect::AlphaFalloff::TOWARD_EQUATOR
                                 : OracleEffect::AlphaFalloff::CONSTANT_HALF;
      state.function = FUNCTIONS[function];
      state.count = state.shape == Shape::SPHERICAL_STAR
                        ? std::min(COUNTS[function], 7)
                        : COUNTS[function];
      state.sides = SIDES[shape];
      state.phase = PHASES[function];
      state.orientation = orientations[(shape + function) % 4];
      matrix[index++] = state;
    }
  return matrix;
}

/** @brief Tests whether a pixel exceeds the brightness threshold. */
inline bool pixel_is_bright(const Pixel &pixel) {
  return static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > BRIGHT_ENERGY;
}

/** @brief Tests whether a pixel has visible coverage. */
inline bool pixel_has_coverage(const Pixel &pixel) {
  return static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > COVERAGE_ENERGY;
}

/** @brief Checks candidate coverage near a reference pixel. */
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
                                 size_t max_high_error_pixels, double max_mae,
                                 double max_rmse, uint32_t max_channel) {
  std::printf(
      "  [%s %s count=%d sides=%d] MAE %.1f/%.1f; RMSE %.1f/%.1f; max channel "
      "%u/%u; energy drift %.3f%%/%.3f%%; high-error px %zu/%zu\n",
      ShapeShifterWhiteBox::shape_option(static_cast<int>(state.shape)),
      ShapeShifterWhiteBox::function_option(static_cast<int>(state.function)),
      state.count, state.sides, error.mean_absolute_error(), max_mae,
      error.root_mean_squared_error(), max_rmse, error.max_absolute_error,
      max_channel, 100.0 * energy_ratio, 100.0 * max_energy_drift,
      high_error_pixels, max_high_error_pixels);
}

/** @brief Asserts the candidate image stays within the visual budget. */
inline void expect_candidate_within_visual_budget(const OracleState &state,
                                                  double max_energy_drift,
                                                  size_t max_high_error_pixels,
                                                  double max_mae,
                                                  double max_rmse,
                                                  uint32_t max_channel) {
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
    if (max_channel_gap(a, b) > HIGH_CHANNEL_ERROR)
      ++high_error_pixels;
  }
  report_visual_budget(state, comparison.error, energy_ratio, high_error_pixels,
                       max_energy_drift, max_high_error_pixels, max_mae,
                       max_rmse, max_channel);
  HS_EXPECT_GT(reference_energy, uint64_t{0});
  HS_EXPECT_LT(comparison.error.mean_absolute_error(), max_mae);
  HS_EXPECT_LT(comparison.error.root_mean_squared_error(), max_rmse);
  HS_EXPECT_LT(comparison.error.max_absolute_error, max_channel);
  HS_EXPECT_LT(energy_ratio, max_energy_drift);
  HS_EXPECT_LT(high_error_pixels, max_high_error_pixels);
  HS_EXPECT_EQ(uncovered_bright_pixels, size_t{0});
  HS_EXPECT_EQ(comparison.reference.backstop_hits, uint32_t{0});
  HS_EXPECT_EQ(comparison.candidate.backstop_hits, uint32_t{0});
}

/** @brief Lowest RGB sum within a pixel radius of a direction.
 * @param frame Rendered oracle frame.
 * @param position Unit direction at the search center.
 * @param radius_px Search radius in pixels, wrapping horizontally.
 * @return Lowest pixel RGB sum, or UINT32_MAX when the disk contains no pixel. */
inline uint32_t minimum_energy_around(const OracleFrame &frame,
                                      const math::Vector &position,
                                      float radius_px) {
  const math::PixelCoords center =
      math::vector_to_pixel<ORACLE_W, ORACLE_H>(position);
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

/** @brief Finds the northernmost row with any nonzero energy. */
inline int first_covered_north_row(const OracleFrame &frame) {
  for (int y = 0; y < 32; ++y)
    for (int x = 0; x < ORACLE_W; ++x) {
      const Pixel &pixel = frame.at(x, y);
      if (static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > 0)
        return y;
    }
  return 32;
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

/** @brief Controls a timed or frozen contour render. */
struct TimedContourRenderOptions {
  float step = 0.0f;          ///< Waveform phase advance over the frame.
  bool timed = true;          ///< Uses the production temporal draw.
  float fraction = 1.0f;      ///< Frozen frame fraction.
  float displacement = 0.0f;  ///< Frozen-reference shift in rows and columns.
  uint32_t *walked = nullptr; ///< Receives column-walk setup count.
};

/** @brief Renders production contours or their frozen geometry.
 * @param state Source geometry and rendering state.
 * @param options Phase motion and frozen-reference controls.
 * @return Row-major normalized mean RGB channels. */
inline std::vector<float>
render_timed_contours(const OracleState &state,
                      const TimedContourRenderOptions &options) {
  const auto &[step, timed, fraction, displacement, walked] = options;
  reset_globals();
  const ScopedTemporalSampling sampling(true);
  hs::set_mock_time(0, 0);
  std::vector<float> luminance;
  {
    OracleEffect effect;
    effect.init();
    ShapeShifterWhiteBox::configure(effect, state);
    const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
    const math::Quaternion start_pose =
        math::make_rotation(math::X_AXIS, 0.27f);
    const math::Quaternion end_pose =
        math::make_rotation(axis, 0.03f) * start_pose;
    ShapeShifterWhiteBox::set_twist(effect, state.phase, step);
    if (timed) {
      ShapeShifterWhiteBox::orientation(effect).set(start_pose);
      ShapeShifterWhiteBox::orientation(effect).push(end_pose);
    } else {
      ShapeShifterWhiteBox::set_twist(
          effect, state.phase - step * (1.0f - fraction), 0.0f);
      ShapeShifterWhiteBox::orientation(effect).set(start_pose);
      ShapeShifterWhiteBox::orientation(effect).push(end_pose);
      const math::Quaternion pose =
          ShapeShifterWhiteBox::orientation(effect).sample(fraction);
      const int column = Spacetime::surface_column(fraction, ORACLE_W);
      const math::Vector meridian = math::pixel_to_vector<ORACLE_W, ORACLE_H>(
          static_cast<float>(column), ORACLE_H / 2.0f);
      const math::Vector across =
          math::cross(math::Y_AXIS, meridian).normalized();
      ShapeShifterWhiteBox::orientation(effect).set(
          math::make_rotation(across,
                              displacement * math::RADIANS_PER_ROW<ORACLE_H>) *
          math::make_rotation(
              math::Y_AXIS, displacement * math::RADIANS_PER_COLUMN<ORACLE_W>) *
          pose);
    }
    {
      Canvas canvas(effect);
      Spacetime::geometry_counters = {};
      ShapeShifterWhiteBox::render_candidate(effect, canvas);
      if (walked)
        *walked = Spacetime::geometry_counters.column_draw_setups;
    }
    effect.advance_display();
    luminance = capture_mean_channel<ORACLE_W, ORACLE_H>(effect);
  }
  Timeline().clear();
  hs::clear_mock_time();
  return luminance;
}
} // namespace hs_test::shapeshifter_oracle_tests
