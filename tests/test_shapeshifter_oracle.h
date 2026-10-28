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
#include <fstream>
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

namespace hs_test {
namespace shapeshifter_oracle_tests {

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

/** @brief Image-dump path prefix, or null when disabled.
 * @return Value of HS_SHAPESHIFTER_IMAGES. */
static inline const char *shapeshifter_image_prefix() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  return std::getenv("HS_SHAPESHIFTER_IMAGES");
#pragma clang diagnostic pop
}

/** @brief Writes an oracle frame as a gamma-encoded P6 image.
 * @param frame Frame to write.
 * @param path Output image path. */
static inline void write_oracle_ppm(const OracleFrame &frame,
                                    const std::string &path) {
  std::ofstream file(path, std::ios::binary);
  HS_EXPECT_TRUE(file.good());
  if (!file)
    return;
  file << "P6\n" << ORACLE_W << " " << ORACLE_H << "\n255\n";
  for (const Pixel &pixel : frame.pixels)
    for (auto value : {pixel.r, pixel.g, pixel.b})
      file.put(char(int(255.0f * sqrtf(value / 65535.0f) + 0.5f)));
}

/** @brief Test-only access to ShapeShifter's pinned state and shipping renderer. */
struct ShapeShifterWhiteBox {
  /** @brief Kernel family used by footprint oracles. */
  using Surface = OracleEffect::Surface;
  static constexpr size_t PRESET_COUNT = OracleEffect::PRESETS.size();

  static_assert(PRESET_COUNT >= 2,
                "the preset transition test steps from preset 0 to 1");

  static void prepare_count(OracleEffect &effect, int count) {
    effect.prepare_count(count);
  }

  static void configure(OracleEffect &effect, const OracleState &state) {
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
                                    : effect.folded_draw_indices[ordinal];
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      const float radius = 2.0f * effect.spaced_radius_t[i];
      const float direction = continuous_star
                                  ? effect.star_phase_direction(radius)
                                  : effect.phase_direction(radius);
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

  static float phase_direction(const OracleEffect &effect, float radius) {
    return effect.phase_direction(radius);
  }

  static float star_phase_direction(const OracleEffect &effect, float radius) {
    return effect.star_phase_direction(radius);
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
    if (mask)
      e.template draw_surface_all<true, true>(canvas, e.selected_palette(),
                                              timed, stats, displacement);
    else
      e.template draw_surface_all<true, false>(canvas, e.selected_palette(),
                                               timed, stats, displacement);
  }
  /** @brief Renders the matched N24 family with a 32 by 32 area quadrature.
   * @param e Prepared effect pose and palette.
   * @param canvas Destination.
   * @param timed Uses each column's authored time. */
  static void surface_full_area(OracleEffect &e, Canvas &canvas, bool timed) {
    using Surface = OracleEffect::Surface;
    constexpr int N = 24;
    constexpr int SAMPLES = 32;
    e.prepare_count(288);
    const auto &palette = e.selected_palette();
    ScratchScope scope_a(scratch_arena_a), scope_b(scratch_arena_b);
    const auto columns = timed ? Scan::prepare_column_basis<288, 144>(
                                     scratch_arena_a, e.orientation, true)
                               : Scan::ColumnBasis{};
    const auto basis = math::make_basis(e.orientation.get(), math::X_AXIS);
    Surface near, far;
    near.bind(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, N);
    far.bind(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, N);
    auto *knots = scratch_arena_b.allocate_n<Surface::Knot>(N + 1);
    auto *colors = scratch_arena_b.allocate_n<Scan::SurfaceColor>(N + 1);
    for (int i = 0; i <= N; ++i) {
      const auto COLOR = palette.get(1 - .5f * i / N);
      const float OPACITY = COLOR.alpha * e.alpha * e.preset_opacity;
      colors[i] = {COLOR.color.r * OPACITY, COLOR.color.g * OPACITY,
                   COLOR.color.b * OPACITY};
    }
    float sine = 0, cosine = 1;
    auto prepare = [&](int column) {
      const float TIME =
          timed ? Spacetime::SamplingSurface<288>::time(column) : 1;
      const float ANGLE =
          math::TWO_PI_F * (e.phase - e.twist_step * (1 - TIME));
      sine = sinf(ANGLE);
      cosine = cosf(ANGLE);
    };
    auto ensure = [&](int side) {
      for (int i = 0; i <= N; ++i) {
        const float Q = .5f * i / N;
        const float RADIUS = math::PI_F * e.screen_balanced_radius_t(Q);
        const float S = sinf(math::TWO_PI_F * Q), C = cosf(math::TWO_PI_F * Q);
        const float WAVE = side ? sine * C - cosine * S : sine * C + cosine * S;
        const float TURN =
            (side && !e.params.opposite ? -1 : 1) * e.params.amplitude * WAVE;
        float angle = TURN - floorf(TURN / math::TWO_PI_F) * math::TWO_PI_F;
        const float SIGN = angle > math::PI_F ? -1.0f : 1.0f;
        if (angle > math::PI_F)
          angle -= math::PI_F;
        float turn_sine, turn_cosine;
        math::poly_sincosf_0_pi(angle, turn_sine, turn_cosine);
        knots[i] = {RADIUS * turn_cosine * SIGN, RADIUS * turn_sine * SIGN,
                    TURN};
      }
      (side ? far : near).set_knots({knots, N + 1});
    };
    auto shade = [&](const Surface::Hit &hit) {
      const auto A = colors[hit.interval], B = colors[hit.interval + 1];
      return Scan::SurfaceColor{A.r + (B.r - A.r) * hit.fraction,
                                A.g + (B.g - A.g) * hit.fraction,
                                A.b + (B.b - A.b) * hit.fraction};
    };
    auto area = [&](const Surface &surface, auto p, float, float,
                    auto &footprint, auto &, auto &, auto *) {
      const auto query = footprint();
      const auto dx = query.dx, dy = query.dy;
      Scan::SurfaceColor color;
      for (int y = 0; y < SAMPLES; ++y)
        for (int x = 0; x < SAMPLES; ++x) {
          const float U = (x + .5f) / SAMPLES - .5f,
                      V = (y + .5f) / SAMPLES - .5f;
          const Surface::Point Q{p.x + U * dx.x + V * dy.x,
                                 p.y + U * dx.y + V * dy.y};
          const auto HIT = surface.sample<false, true, false>(
              Q, Scan::radial_atan2(Q.y, Q.x));
          if (HIT.found()) {
            const auto C = shade(HIT);
            color.r += C.r / (SAMPLES * SAMPLES);
            color.g += C.g / (SAMPLES * SAMPLES);
            color.b += C.b / (SAMPLES * SAMPLES);
          }
        }
      return color;
    };
    auto lower_bound = [](float) { return 0; };
    Scan::draw_radial_surface<288, 144>(canvas, near, far, area, basis, columns,
                                        prepare, ensure, lower_bound, shade);
  }
  static int preset_frames() { return OracleEffect::PRESET_FRAMES; }
  static int preset_dwell_frames() { return OracleEffect::PRESET_DWELL_FRAMES; }
  static int preset_segue_frames() { return OracleEffect::DEPARTURE.frames; }
};

/** @brief Captures one renderer callable from a fresh deterministic effect.
 * @details Temporal sampling is off: these oracles judge the shapes' still
 * geometry, which the timed route draws at each column's time. */
template <typename Render>
inline OracleFrame capture_frame(const OracleState &state, Render &&render) {
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

/** @brief Frame energy. */
inline uint64_t frame_energy(const OracleFrame &frame) {
  return hs_test::frame_energy(frame.pixels);
}

/** @brief Row energy. */
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
  const Shape shapes[] = {Shape::PLANAR_POLYGON, Shape::SPHERICAL_POLYGON,
                          Shape::FLOWER, Shape::PLANAR_STAR,
                          Shape::SPHERICAL_STAR};
  const Function functions[] = {Function::SINE, Function::TRIANGLE,
                                Function::SAWTOOTH, Function::SQUARE};
  const int counts[] = {1, 2, 7, 75};
  const int sides[] = {3, 5, 9, 16, 7};
  const float phases[] = {0.0f, 0.249f, 0.5f, 0.999f};
  const math::Quaternion orientations[] = {
      math::Quaternion(), math::make_rotation(math::X_AXIS, math::Y_AXIS),
      math::make_rotation(math::X_AXIS, -math::Y_AXIS),
      math::make_rotation(math::X_AXIS, math::Z_AXIS)};

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

/** @brief Image comparison reports the expected error statistics. */
inline void test_buffer_comparator_statistics() {
  OracleFrame reference;
  reference.pixels.resize(static_cast<size_t>(ORACLE_W) * ORACLE_H);
  OracleFrame candidate = reference;

  FrameErrorStats exact = compare_buffers(reference.pixels, candidate.pixels);
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
  FrameErrorStats changed = compare_buffers(reference.pixels, candidate.pixels);
  HS_EXPECT_FALSE(changed.exact());
  HS_EXPECT_EQ(changed.different_pixels, static_cast<size_t>(1));
  HS_EXPECT_EQ(changed.different_channels, static_cast<size_t>(2));
  HS_EXPECT_EQ(changed.max_absolute_error, static_cast<uint32_t>(11));
  HS_EXPECT_EQ(changed.total_absolute_error, static_cast<uint64_t>(18));
  HS_EXPECT_EQ(changed.total_squared_error, static_cast<uint64_t>(170));
  HS_EXPECT_EQ(changed.first_different_pixel, static_cast<int>(PIXEL_INDEX));
  HS_EXPECT_EQ(changed.lit_reference_pixels, static_cast<size_t>(0));
  HS_EXPECT_EQ(changed.coverage_mismatches, static_cast<size_t>(1));

  // Accumulating over split spans indexes across the whole comparison.
  const std::span<const Pixel> reference_span(reference.pixels);
  const std::span<const Pixel> candidate_span(candidate.pixels);
  const FrameErrorStats head =
      compare_buffers(reference_span.first(PIXEL_INDEX - 1),
                      candidate_span.first(PIXEL_INDEX - 1));
  const FrameErrorStats split =
      compare_buffers(reference_span.subspan(PIXEL_INDEX - 1),
                      candidate_span.subspan(PIXEL_INDEX - 1), head);
  HS_EXPECT_EQ(split.compared_pixels, changed.compared_pixels);
  HS_EXPECT_EQ(split.total_squared_error, changed.total_squared_error);
  HS_EXPECT_EQ(split.first_different_pixel, static_cast<int>(PIXEL_INDEX));
}

/**
 * @brief Verifies the reference render is reproducible and non-black over the
 *        whole matrix.
 * @details Both sides are reference_renderer(), so this pins determinism and
 *          nonzero energy, not correctness.
 */
inline void test_reference_matrix_is_deterministic_and_nonblack() {
  for (const OracleState &state : shape_function_matrix()) {
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

/** @brief Candidate matrix stays within visual budget. */
inline void test_candidate_matrix_stays_within_visual_budget() {
  struct Budget {
    double mae;
    double rmse;
    uint32_t channel;
    size_t pixels;
  };
  // Per-case native baselines with aggregate and peak headroom.
  constexpr Budget budgets[] = {
      {0.4, 5, 480, 1},         {2.8, 36, 3800, 1},   {9, 87, 9700, 9},
      {30, 172, 10910, 17},     {0.3, 1.5, 90, 1},    {0.4, 3.5, 210, 1},
      {6.2, 82, 7720, 7},       {54, 232, 11340, 44}, {0.4, 6, 650, 1},
      {2.9, 52, 5010, 1},       {9.7, 102, 6190, 7},  {188, 806, 23780, 541},
      {1.7, 15, 960, 1},        {12, 106, 3810, 1},   {39, 294, 10860, 111},
      {351, 1192, 20450, 1845}, {1.2, 14, 1250, 1},   {9.5, 101, 4550, 1},
      {23, 185, 11160, 25},     {15, 135, 8530, 9}};
  const auto matrix = shape_function_matrix();
  static_assert(std::size(budgets) == std::tuple_size_v<decltype(matrix)>);
  for (size_t index = 0; index < matrix.size(); ++index) {
    const auto &budget = budgets[index];
    expect_candidate_within_visual_budget(matrix[index], MAX_ENERGY_DRIFT,
                                          budget.pixels, budget.mae,
                                          budget.rmse, budget.channel);
  }
}

/** @brief Star projection policies render different edges. */
inline void test_star_projection_policies_render_different_edges() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SQUARE;
  state.count = 3;
  state.sides = 5;
  state.phase = 0.17f;
  state.orientation =
      math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
  const OracleFrame planar = capture_frame(state, candidate_renderer());
  state.shape = OracleEffect::ShapeType::SPHERICAL_STAR;
  const OracleFrame spherical = capture_frame(state, candidate_renderer());
  HS_EXPECT_GT(
      compare_buffers(planar.pixels, spherical.pixels).different_pixels,
      size_t{100});
}

/** @brief Screen balanced spacing is opt in for every shape. */
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
    state.orientation =
        math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized();
    const OracleFrame uniform = capture_frame(state, candidate_renderer());
    state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
    const OracleFrame balanced = capture_frame(state, candidate_renderer());
    HS_EXPECT_GT(
        compare_buffers(uniform.pixels, balanced.pixels).different_pixels,
        size_t{100});
  }
}

/** @brief Star options and shipping presets are planar. */
inline void test_star_options_and_shipping_presets_are_planar() {
  using Falloff = OracleEffect::AlphaFalloff;
  using Shape = OracleEffect::ShapeType;
  HS_EXPECT_EQ(ShapeShifterWhiteBox::shape_count(), 5);
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::shape_option(3), "Planar Star");
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::shape_option(4), "Spherical Star");
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::shape_export_option(3),
                  "ShapeType::PLANAR_STAR");
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::shape_export_option(4),
                  "ShapeType::SPHERICAL_STAR");
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::spacing_option(1), "Screen Balanced");
  HS_EXPECT_STREQ(ShapeShifterWhiteBox::spacing_export_option(1),
                  "RadiusSpacing::SCREEN_BALANCED");
  size_t planar_stars = 0;
  for (size_t index = 0; index < ShapeShifterWhiteBox::PRESET_COUNT; ++index) {
    const Shape shape = ShapeShifterWhiteBox::preset_shape(index);
    const bool planar_star = shape == Shape::PLANAR_STAR;
    HS_EXPECT_NE(shape, Shape::SPHERICAL_STAR);
    const Falloff falloff = ShapeShifterWhiteBox::preset_alpha_falloff(index);
    HS_EXPECT_EQ(planar_star, falloff == Falloff::TOWARD_EQUATOR);
    if (planar_star)
      ++planar_stars;
  }
  HS_EXPECT_GT(planar_stars, size_t{0});
}

/** @brief High count uniform star stays within visual budget. */
inline void test_high_count_uniform_star_stays_within_visual_budget() {
  using Function = OracleEffect::PhaseFunction;
  using Shape = OracleEffect::ShapeType;
  const std::array<math::Quaternion, 6> orientations = {{
      math::Quaternion(),
      math::make_rotation(math::X_AXIS, math::Y_AXIS),
      math::make_rotation(math::X_AXIS, -math::Y_AXIS),
      math::make_rotation(math::X_AXIS, math::Z_AXIS),
      math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(),
      math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
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
        state, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT, MAX_STAR_HIGH_ERROR_PIXELS,
        MAX_MEAN_ABSOLUTE_ERROR, MAX_ROOT_MEAN_SQUARED_ERROR,
        MAX_CHANNEL_ERROR);
  }
}

/** @brief Screen-balanced star strokes stay within their visual budget. */
inline void test_screen_balanced_star_stroke_stays_within_visual_budget() {
  OracleState state;
  state.preset_index = 1;
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
  expect_candidate_within_visual_budget(
      state, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT, MAX_STAR_HIGH_ERROR_PIXELS,
      MAX_MEAN_ABSOLUTE_ERROR, MAX_ROOT_MEAN_SQUARED_ERROR, MAX_CHANNEL_ERROR);
}

/** @brief Frozen temporal star strokes preserve the still stroke look. */
inline void test_timed_dense_star_stroke_stays_within_visual_budget() {
  OracleState state;
  state.preset_index = 1;
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
  const auto comparison = compare_renders(
      state, candidate_renderer(), [](OracleEffect &effect, Canvas &canvas) {
        const ScopedTemporalSampling sampling(true);
        ShapeShifterWhiteBox::render_candidate(effect, canvas);
      });
  const auto energy = frame_energy(comparison.reference);
  const double drift =
      std::fabs(double(frame_energy(comparison.candidate)) - double(energy)) /
      double(energy);
  std::printf(
      "  frozen temporal star stroke look: MAE=%.3f RMSE=%.3f energy_drift=%.6f reference=%f candidate=%f\n",
      comparison.error.mean_absolute_error(),
      comparison.error.root_mean_squared_error(), drift, double(energy),
      double(frame_energy(comparison.candidate)));
  if (const char *prefix = shapeshifter_image_prefix()) {
    write_oracle_ppm(comparison.reference,
                     std::string(prefix) + "-reference.ppm");
    write_oracle_ppm(comparison.candidate,
                     std::string(prefix) + "-temporal-frozen.ppm");
  }
  HS_EXPECT_LT(comparison.error.mean_absolute_error(), MAX_MEAN_ABSOLUTE_ERROR);
  HS_EXPECT_LT(comparison.error.root_mean_squared_error(),
               MAX_ROOT_MEAN_SQUARED_ERROR);
  HS_EXPECT_LT(drift, 0.014);
  const auto pole_reference = row_energy(comparison.reference, 0);
  const auto pole_candidate = row_energy(comparison.candidate, 0);
  HS_EXPECT_GE(pole_candidate * 10, pole_reference * 9);
  std::printf(
      "  frozen star stroke row-zero energy: reference=%llu candidate=%llu\n",
      static_cast<unsigned long long>(pole_reference),
      static_cast<unsigned long long>(pole_candidate));
}

/** @brief High count uniform star covers north pole. */
inline void test_high_count_uniform_star_covers_north_pole() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SINE;
  state.count = 144;
  state.sides = 7;
  state.phase = 0.125f;
  state.alpha = 0.274f;
  state.orientation = math::make_rotation(math::X_AXIS, math::Y_AXIS);
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

/** @brief Minimum energy around. */
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

/** @brief High count planar star caps cover chart centers. */
inline void test_high_count_planar_star_caps_cover_chart_centers() {
  constexpr uint32_t MIN_CAP_ENERGY = 60000;
  const std::array<math::Quaternion, 4> orientations = {{
      math::Quaternion(),
      math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(),
      math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(),
      math::Quaternion(0.51f, 0.63f, -0.28f, 0.51f).normalized(),
  }};
  for (const auto [count, spacing] :
       {std::pair{144, OracleEffect::RadiusSpacing::UNIFORM},
        std::pair{288, OracleEffect::RadiusSpacing::SCREEN_BALANCED}}) {
    for (const math::Quaternion &orientation : orientations) {
      OracleState state;
      state.preset_index = 1;
      state.shape = OracleEffect::ShapeType::PLANAR_STAR;
      state.function = OracleEffect::PhaseFunction::SINE;
      state.spacing = spacing;
      state.count = count;
      state.sides = 7;
      state.phase = 0.125f;
      state.alpha = 0.274f;
      state.orientation = orientation;
      const OracleFrame frame = capture_frame(state, candidate_renderer());
      const math::Vector near = math::rotate(math::X_AXIS, orientation);
      const math::Vector far = -near;
      HS_EXPECT_GE(minimum_energy_around(frame, near, 1.5f), MIN_CAP_ENERGY);
      HS_EXPECT_GE(minimum_energy_around(frame, far, 1.5f), MIN_CAP_ENERGY);
    }
  }
}

/** @brief Finds the northernmost row with visible coverage. */
inline int first_covered_north_row(const OracleFrame &frame) {
  for (int y = 0; y < 32; ++y)
    for (int x = 0; x < ORACLE_W; ++x) {
      const Pixel &pixel = frame.at(x, y);
      if (static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > 0)
        return y;
    }
  return 32;
}

/** @brief High count spherical star contours reach display north. */
inline void test_high_count_spherical_star_contours_reach_display_north() {
  for (float phase : {0.0f, 0.125f, 0.249f, 0.5f, 0.75f, 0.999f}) {
    OracleState state;
    state.shape = OracleEffect::ShapeType::SPHERICAL_STAR;
    state.function = OracleEffect::PhaseFunction::SINE;
    state.count = 144;
    state.sides = 7;
    state.phase = phase;
    state.alpha = 0.274f;
    state.orientation = math::Quaternion();
    const OracleFrame frame = capture_frame(state, candidate_renderer());
    HS_EXPECT_LE(first_covered_north_row(frame), 1);
  }
}

/** @brief Amplitude preserves sweep velocity. */
inline void test_amplitude_preserves_sweep_velocity() {
  OracleEffect effect;
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 1.0f),
                 0.04f, 1e-6f);
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.04f, 4.0f),
                 0.01f, 1e-6f);
}

/** @brief Phase step stops at the nyquist limit. */
inline void test_phase_step_stops_at_the_nyquist_limit() {
  OracleEffect effect;
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.049f, 0.1f),
                 0.49f, 1e-6f);
  HS_EXPECT_NEAR(ShapeShifterWhiteBox::advance_phase(effect, 0.16f, 0.1f), 0.5f,
                 1e-6f);
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

/** @brief Alpha falloff fades to equator. */
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

/** @brief Opposite halves direction. */
inline void test_opposite_halves_direction() {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    HS_EXPECT_EQ(effect.updateParameter("Opposite", 0.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, 1.5f), -1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::star_phase_direction(effect, 0.5f),
                 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::star_phase_direction(effect, 1.5f),
                 1.0f);
    HS_EXPECT_EQ(effect.updateParameter("Opposite", 1.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_direction(effect, 1.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::star_phase_direction(effect, 0.5f),
                 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::star_phase_direction(effect, 1.5f),
                 -1.0f);
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
  HS_EXPECT_GT(
      compare_buffers(unchecked.pixels, checked.pixels).different_pixels,
      static_cast<size_t>(0));
}

/** @brief Folded shapes draw from equator to both poles. */
inline void test_folded_shapes_draw_from_equator_to_both_poles() {
  HS_EXPECT_EQ(ShapeShifterWhiteBox::folded_draw_order(6),
               (std::vector<int>{3, 2, 4, 1, 5, 0}));
  HS_EXPECT_EQ(ShapeShifterWhiteBox::folded_draw_order(5),
               (std::vector<int>{2, 3, 1, 4, 0}));
}

/** @brief Screen balanced spacing follows sampling envelope. */
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
  const float density_scale = COUNT * math::PI_F / (2.0f * DENSITY_INTEGRAL);
  HS_EXPECT_NEAR(equator_equivalent_count, density_scale, 0.2f);
  HS_EXPECT_NEAR(equator_equivalent_count, OracleEffect::MAX_SHAPES, 1.8f);
  HS_EXPECT_NEAR(pole_equivalent_count, density_scale * DENSITY_FLOOR, 0.2f);
  HS_EXPECT_NEAR(pole_equivalent_count, ORACLE_H, 1.0f);
  HS_EXPECT_NEAR(mapped(0) * math::PI_F, 2.0f * DENSITY_INTEGRAL / COUNT,
                 1e-6f);

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

/** @brief Preset transition adopts next params. */
inline void test_preset_transition_adopts_next_params() {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    Canvas canvas(effect);
    HS_EXPECT_EQ(effect.updateParameter("Alpha", 0.42f),
                 ParamSetResult::APPLIED);
    for (int frame = 0; frame < ShapeShifterWhiteBox::preset_frames(); ++frame)
      ShapeShifterWhiteBox::step_choreography(effect, canvas);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});

    auto value = [&](const char *name) {
      const auto *def = hs_test::required_param(effect, name);
      return def ? def->get() : -1.0f;
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

/**
 * @brief Pins ShapeShifter's preset cadence with an eight-frame fade out,
 *        adoption in the dark, and an eight-frame fade back in.
 */
inline void test_preset_transition_fades_through_black_in_16_frames() {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    Canvas canvas(effect);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_segue_frames(), 16);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_dwell_frames() +
                     ShapeShifterWhiteBox::preset_segue_frames(),
                 ShapeShifterWhiteBox::preset_frames());
    const auto FIRST = ShapeShifterWhiteBox::preset_shape(0);
    const auto SECOND = ShapeShifterWhiteBox::preset_shape(1);
    HS_EXPECT_NE(static_cast<int>(FIRST), static_cast<int>(SECOND));
    auto step = [&] {
      ShapeShifterWhiteBox::step_choreography(effect, canvas);
    };

    for (int frame = 0; frame < ShapeShifterWhiteBox::preset_dwell_frames();
         ++frame)
      step();
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 1.0f, 1e-6f);
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(FIRST));

    for (int frame = 1; frame < 8; ++frame) {
      step();
      HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect),
                     1.0f - static_cast<float>(frame) / 8.0f, 1e-6f);
      HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                   static_cast<int>(FIRST));
    }

    step();
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 0.0f, 1e-6f);
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(SECOND));
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});

    for (int frame = 1; frame <= 8; ++frame) {
      step();
      HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect),
                     static_cast<float>(frame) / 8.0f, 1e-6f);
    }
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(SECOND));
  }
  Timeline().clear();
}

/** @brief Count slider matches draw limit. */
inline void test_count_slider_matches_draw_limit() {
  reset_globals();
  {
    ShapeShifter<96, 20> effect;
    effect.init();
    const auto *count = effect.getParameters().find("Count");
    HS_EXPECT_TRUE(count != nullptr);
    if (!count)
      return;
    HS_EXPECT_EQ(count->min, 1.0f);
    HS_EXPECT_EQ(count->max, 40.0f);
    HS_EXPECT_EQ(effect.updateParameter("Count", 288.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(count->get(), 40.0f);
  }
  Timeline().clear();
}

/** @brief Renders a timed contour stack or its frozen frame.
 * @param state Shape, count, sides and amplitude.
 * @param step Waveform phase advance over the frame, turns.
 * @param timed Uses the production temporal draw.
 * @param fraction Frozen frame fraction.
 * @param displacement Frozen-reference shift in rows and columns.
 * @param walked Receives the production draw's column-walk setups, when
 * given.
 * @return Row-major linear luminance.
 * @details The orientation turns about a tilted axis over the frame and the
 * waveform phase advances by @p step. The frozen frame holds the phase and
 * the orientation of @p fraction, turned @p displacement columns about the
 * vertical and @p displacement rows along the meridians of @p fraction's
 * columns, and draws through the same timed route with no motion. Chord and
 * point samples are single samples at the point route's tolerance, held to a
 * reference displaced along both. */
inline std::vector<float> render_timed_contours(const OracleState &state,
                                                float step, bool timed,
                                                float fraction,
                                                float displacement,
                                                uint32_t *walked = nullptr) {
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
          ShapeShifterWhiteBox::orientation(effect).get(fraction);
      int column = 0;
      for (int c = 1; c < ORACLE_W / 2; ++c)
        if (fabsf(Spacetime::SamplingSurface<ORACLE_W>::time(c) - fraction) <
            fabsf(Spacetime::SamplingSurface<ORACLE_W>::time(column) -
                  fraction))
          column = c;
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
    luminance = capture_luminance<ORACLE_W, ORACLE_H>(effect);
  }
  Timeline().clear();
  hs::clear_mock_time();
  return luminance;
}

/** @brief Timed contours match frozen columns for every shape: geodesic
 *         contours through the column walk, chart-straight ones through
 *         moving chords. */
inline void test_timed_contours_match_frozen_columns() {
  using ShapeType = OracleEffect::ShapeType;
  struct Case {
    ShapeType shape;
    int count;
    int sides;
    float amplitude;
    float step;
    const char *name;
  };
  const Case cases[] = {
      {ShapeType::SPHERICAL_POLYGON, 24, 4, 2.0f, 0.01f, "spherical polygon"},
      {ShapeType::SPHERICAL_STAR, 12, 5, 1.0f, 0.016f, "spherical star"},
      {ShapeType::PLANAR_STAR, 48, 7, 1.0f, 0.016f, "planar star"},
      {ShapeType::PLANAR_STAR, 12, 5, 1.0f, 0.016f, "sparse planar star"},
      {ShapeType::PLANAR_POLYGON, 16, 5, 1.0f, 0.016f, "planar polygon"},
      {ShapeType::FLOWER, 24, 3, 1.0f, 0.0186f, "flower"}};
  for (const Case &test : cases) {
    OracleState state;
    state.shape = test.shape;
    state.count = test.count;
    state.sides = test.sides;
    state.amplitude = test.amplitude;
    state.phase = 0.37f;
    uint32_t walked = 0;
    render_timed_contours(state, test.step, true, 1.0f, 0.0f, &walked);
    const bool geodesic = test.shape == ShapeType::SPHERICAL_POLYGON ||
                          test.shape == ShapeType::SPHERICAL_STAR;
    HS_CONTEXT("shape", static_cast<int>(test.shape));
    if (geodesic)
      HS_EXPECT_GT(walked, 0u);
    hs_test::seam_compare_render<ORACLE_W, ORACLE_H>(
        [&](bool timed, float fraction, float displacement) {
          return render_timed_contours(state, test.step, timed, fraction,
                                       displacement);
        },
        test.name);
  }
}

/** @brief Shipping presets match their frozen columns at the authored speed. */
inline void test_timed_presets_match_frozen_columns() {
  for (size_t i = 0; i < ShapeShifterWhiteBox::PRESET_COUNT; ++i) {
    if (i == 0)
      continue;

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
    const float step = fminf(preset.speed / preset.amplitude, 0.5f);
    const std::string name = "preset " + std::to_string(i + 1);
    HS_CONTEXT("preset", static_cast<int>(i + 1));
    hs_test::seam_compare_render<ORACLE_W, ORACLE_H>(
        [&](bool timed, float fraction, float displacement) {
          return render_timed_contours(state, step, timed, fraction,
                                       displacement);
        },
        name.c_str());
    state.orientation =
        math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized();
    const auto look = compare_renders(
        state, candidate_renderer(), [](OracleEffect &effect, Canvas &canvas) {
          const ScopedTemporalSampling sampling(true);
          ShapeShifterWhiteBox::render_candidate(effect, canvas);
        });
    const double energy = double(frame_energy(look.reference));
    const double drift =
        (double(frame_energy(look.candidate)) - energy) / energy;
    std::printf(
        "  frozen style preset %zu: MAE=%.3f RMSE=%.3f max=%u energy_drift=%.6f\n",
        i + 1, look.error.mean_absolute_error(),
        look.error.root_mean_squared_error(), look.error.max_absolute_error,
        drift);
    HS_EXPECT_LT(std::fabs(drift), MAX_PRESET_ENERGY_DRIFT);
    HS_EXPECT_LT(look.error.mean_absolute_error(), MAX_MEAN_ABSOLUTE_ERROR);
    HS_EXPECT_LT(look.error.root_mean_squared_error(),
                 MAX_ROOT_MEAN_SQUARED_ERROR);
    HS_EXPECT_LT(look.error.max_absolute_error, MAX_CHANNEL_ERROR);
    if (const char *prefix = shapeshifter_image_prefix()) {
      const std::string path =
          std::string(prefix) + "-p" + std::to_string(i + 1);
      write_oracle_ppm(look.reference, path + "-still.ppm");
      write_oracle_ppm(look.candidate, path + "-timed.ppm");
    }
  }
}

/** @brief Clipping moving flowers preserves visible pixels and skips chart samples. */
inline void test_timed_flower_clip_culls_anchor_intervals() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::FLOWER;
  state.count = 70;
  state.sides = 3;
  state.phase = 0.37f;
  const auto full = render_timed_contours(state, 0.0186f, true, 1.0f, 0.0f);
  const auto full_counts = Spacetime::geometry_counters;
  state.clip = {0, 72, 0, 144};
  const auto clipped = render_timed_contours(state, 0.0186f, true, 1.0f, 0.0f);
  const auto clipped_counts = Spacetime::geometry_counters;
  HS_EXPECT_GT(clipped_counts.chord_fast_samples, 0u);
  HS_EXPECT_GT(clipped_counts.chord_intervals_culled, 0u);
  HS_EXPECT_GT(clipped_counts.polar_column_samples, 0u);
  HS_EXPECT_LT(clipped_counts.polar_column_samples,
               full_counts.polar_column_samples);
  for (int y = 0; y < 72; ++y)
    for (int x = 0; x < 144; ++x)
      HS_EXPECT_NEAR(clipped[y * ORACLE_W + x], full[y * ORACLE_W + x], 1e-5f);
  std::printf(
      "  flower clip: intervals=%u culled=%u coherent=%u screen=%u projections=%u column_samples=%u; full projections=%u\n",
      clipped_counts.chord_intervals, clipped_counts.chord_intervals_culled,
      clipped_counts.chord_fast_samples, clipped_counts.chord_screen_samples,
      clipped_counts.point_piece_evaluations, clipped_counts.column_samples,
      full_counts.point_piece_evaluations);
  std::printf(
      "  flower solves: anchors=%u point_pieces=%u candidates=%u\n",
      clipped_counts.chord_capture_solves, clipped_counts.point_capture_solves,
      clipped_counts.chord_fast_samples + clipped_counts.chord_screen_samples);
  state.preset_index = 1;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.count = 288;
  state.sides = 7;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  render_timed_contours(state, 0.016f, true, 1.0f, 0.0f);
  const auto star_counts = Spacetime::geometry_counters;
  HS_EXPECT_GT(star_counts.polar_column_samples, 0u);
  std::printf(
      "  polar visits/captures: flower=%u/%u full=%u/%u star=%u/%u\n",
      clipped_counts.polar_column_visits, clipped_counts.polar_geometry_samples,
      full_counts.polar_column_visits, full_counts.polar_geometry_samples,
      star_counts.polar_column_visits, star_counts.polar_geometry_samples);
  std::printf("  polar geometry captures: flower=%u full=%u star=%u\n",
              clipped_counts.polar_geometry_samples,
              full_counts.polar_geometry_samples,
              star_counts.polar_geometry_samples);
  std::printf("  polar clip samples: flower=%u full=%u star=%u\n",
              clipped_counts.polar_column_samples,
              full_counts.polar_column_samples,
              star_counts.polar_column_samples);
  std::printf("  dense star clip: intervals=%u culled=%u coherent=%u screen=%u "
              "projections=%u column_samples=%u\n",
              star_counts.chord_intervals, star_counts.chord_intervals_culled,
              star_counts.chord_fast_samples, star_counts.chord_screen_samples,
              star_counts.point_piece_evaluations, star_counts.column_samples);
}

/** @brief Chart-edge rotation bounds include bowed interiors and the equator. */
inline void test_chart_edge_rotation_radius() {
  using Chords = Plot::PlanarChords<ORACLE_W, ORACLE_H>;
  /** @brief One straight source-chart segment. */
  struct Edge {
    math::PixelCoords a, b; ///< Chart endpoints, radians.
  };
  const Edge edges[] = {{{0.0f, 0.0f}, {0.0f, 0.0f}},
                        {{0.2f, 0.0f}, {0.8f, 0.0f}},
                        {{3.0f * cosf(math::PI_F / 6.0f), -1.5f},
                         {3.0f * cosf(math::PI_F / 6.0f), 1.5f}},
                        {{1.1f, -1.9052559f}, {1.1f, 1.9052559f}},
                        {{-2.5f, 0.1f}, {2.5f, 0.1f}}};
  for (const auto &edge : edges) {
    const float bound = Chords::rotation_radius(edge.a, edge.b);
    double maximum = 0.0;
    for (int sample = 0; sample <= 200; ++sample) {
      const double t = sample / 200.0;
      const double x = edge.a.x + (edge.b.x - edge.a.x) * t;
      const double z = edge.a.y + (edge.b.y - edge.a.y) * t;
      const double radius = std::hypot(x, z);
      const double physical = std::abs(std::sin(radius));
      maximum = std::max(maximum, physical);
      HS_EXPECT_LE(physical, double(bound) + 2e-6);
    }
    HS_EXPECT_GE(bound, 0.0f);
    HS_EXPECT_LE(bound, 1.0f);
    HS_EXPECT_LT(double(bound) - maximum, 0.0001);
  }
  HS_EXPECT_GT(Chords::rotation_radius(edges[2].a, edges[2].b), 0.5f);
}

/** @brief Polar captures match exact rotated chart segments at arm cuts. */
inline void test_polar_walk_matches_rotated_chart() {
  using Walk = Plot::PolarWalk<ORACLE_W, ORACLE_H>;
  const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
  float worst = 0.0f;
  for (bool north : {true, false}) {
    for (int column : {0, 1, 71, 143, 144, 145, 215, 287}) {
      const float lambda = column * math::RADIANS_PER_COLUMN<ORACLE_W>;
      auto source = [&](float u) {
        const float x = 8.0f * cosf(lambda - 0.15f) * (1.0f - u) +
                        9.0f * cosf(lambda + 0.15f) * u;
        const float z = 8.0f * sinf(lambda - 0.15f) * (1.0f - u) +
                        9.0f * sinf(lambda + 0.15f) * u;
        const float rho = hypotf(x, z);
        const float phi = rho * math::RADIANS_PER_ROW<ORACLE_H>;
        return math::Vector(x * sinf(phi) / rho,
                            (north ? 1.0f : -1.0f) * cosf(phi),
                            z * sinf(phi) / rho);
      };
      auto exact = [&](float u, float t) {
        return math::rotate(source(u),
                            math::make_rotation(axis, 0.04f * (t - 1.0f)));
      };
      auto capture = [&](float u) {
        return Plot::PolarAnchor::capture(Walk::project(exact(u, 0.0f), north),
                                          Walk::project(exact(u, 0.5f), north),
                                          Walk::project(exact(u, 1.0f), north));
      };
      const auto a = capture(0.0f), b = capture(1.0f);
      const float t = Spacetime::SamplingSurface<ORACLE_W>::time(column);
      const float dx = cosf(lambda), dz = sinf(lambda);
      float u = 0.0f;
      HS_EXPECT_TRUE(Walk::crossing(a.at(t), b.at(t), dx, dz, u));
      float lo = 0.0f, hi = 1.0f;
      const auto first = exact(0.0f, t);
      const float start_side = first.x * dz - first.z * dx;
      for (int iteration = 0; iteration < 30; ++iteration) {
        const float mid = (lo + hi) * 0.5f;
        const auto q = exact(mid, t);
        if ((q.x * dz - q.z * dx) * start_side > 0.0f)
          lo = mid;
        else
          hi = mid;
      }
      const auto q = exact((lo + hi) * 0.5f, t);
      const float expected = atan2f(hypotf(q.x, q.z), north ? q.y : -q.y) *
                             math::ROWS_PER_RADIAN<ORACLE_H>;
      const auto from = a.at(t), to = b.at(t);
      const float actual =
          hypotf(from.x + u * (to.x - from.x), from.y + u * (to.y - from.y));
      worst = fmaxf(worst, fabsf(actual - expected));
      HS_EXPECT_LT(fabsf(actual - expected), 0.05f);
    }
  }
  std::printf("  polar exact intersection max error: %.6f px\n", worst);
}

inline void write_shapeshifter_preset_hashes() {
  for (bool temporal : {false, true})
    for (size_t preset = 0; preset < ShapeShifterWhiteBox::PRESET_COUNT;
         ++preset) {
      reset_globals();
      hs::seed_effect_stream(hs::effect_seed_key<OracleEffect>("ShapeShifter"));
      hs::set_mock_time(0, 0);
      uint64_t hash = 14695981039346656037ull;
      {
        OracleEffect effect;
        effect.init();
        effect.profile_select_preset(preset);
        effect.set_clip(0, 72, 0, 144);
        const ScopedTemporalSampling sampling(temporal);
        for (int frame = 0; frame < 8; ++frame) {
          hs::set_mock_time(frame * 62500ul / 1000, frame * 62500ul);
          effect.draw_frame();
          effect.advance_display();
          for (int y = 0; y < 72; ++y)
            for (int x = 0; x < 144; ++x) {
              const auto p = effect.get_pixel(x, y);
              for (const uint16_t channel : {p.r, p.g, p.b})
                hash = hs_test::fnv1a64_channel(hash, channel);
            }
        }
      }
      Timeline().clear();
      HS_EXPECT_NE(hash, uint64_t{0});
      std::printf(
          "SHAPESHIFTER_PRESET_HASH preset=%zu temporal=%d frames=8 hash=0x%016llx\n",
          preset + 1, int(temporal), static_cast<unsigned long long>(hash));
    }
  hs::clear_mock_time();
}

#include "tests/test_shapeshifter_surface.h"

/** @brief Runs shapeshifter oracle tests. */
inline int run_shapeshifter_oracle_tests() {
  ModuleFixture fixture("shapeshifter_oracle");
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (std::getenv("HS_SHAPESHIFTER_SURFACE_CENSUS")) {
    test_surface_natural_census();
    return fixture.result();
  }
  if (std::getenv("HS_SURFACE_DARK_DIAGNOSTICS")) {
    HS_RUN(test_surface_footprint_darkening);
    return fixture.result();
  }
  if (std::getenv("HS_SURFACE_LOBE_DIAGNOSTICS")) {
    test_surface_lobe_outliers();
    return fixture.result();
  }
  if (std::getenv("HS_SHAPESHIFTER_PRESET_HASHES")) {
    write_shapeshifter_preset_hashes();
    return fixture.result();
  }
#pragma clang diagnostic pop
  if (shapeshifter_image_prefix()) {
    test_surface_images();
    return fixture.result();
  }
  HS_RUN(radial_surface_tests::run_radial_surface_checks);
  HS_RUN(radial_surface_mask_tests::test_radial_surface_mask);
  HS_RUN(radial_surface_mask_tests::test_radial_surface_area);
  HS_RUN(test_surface_seams);
  HS_RUN(test_surface_clips);
  HS_RUN(test_surface_mask_parity);
  HS_RUN(test_surface_preset_dispatch);
  HS_RUN(test_surface_max_amplitude);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (std::getenv("HS_SHAPESHIFTER_SURFACE_ONLY"))
    return fixture.result();
#pragma clang diagnostic pop
  HS_RUN(test_chart_edge_rotation_radius);
  HS_RUN(test_polar_walk_matches_rotated_chart);
  HS_RUN(test_timed_contours_match_frozen_columns);
  HS_RUN(test_timed_presets_match_frozen_columns);
  HS_RUN(test_timed_flower_clip_culls_anchor_intervals);
  HS_RUN(test_count_slider_matches_draw_limit);
  HS_RUN(test_buffer_comparator_statistics);
  HS_RUN(test_reference_matrix_is_deterministic_and_nonblack);
  HS_RUN(test_candidate_matrix_stays_within_visual_budget);
  HS_RUN(test_star_projection_policies_render_different_edges);
  HS_RUN(test_screen_balanced_spacing_is_opt_in_for_every_shape);
  HS_RUN(test_star_options_and_shipping_presets_are_planar);
  HS_RUN(test_high_count_uniform_star_stays_within_visual_budget);
  HS_RUN(test_screen_balanced_star_stroke_stays_within_visual_budget);
  HS_RUN(test_timed_dense_star_stroke_stays_within_visual_budget);
  HS_RUN(test_high_count_uniform_star_covers_north_pole);
  HS_RUN(test_high_count_planar_star_caps_cover_chart_centers);
  HS_RUN(test_high_count_spherical_star_contours_reach_display_north);
  HS_RUN(test_amplitude_preserves_sweep_velocity);
  HS_RUN(test_phase_step_stops_at_the_nyquist_limit);
  HS_RUN(test_alpha_falloff_fades_to_equator);
  HS_RUN(test_opposite_halves_direction);
  HS_RUN(test_folded_shapes_draw_from_equator_to_both_poles);
  HS_RUN(test_screen_balanced_spacing_follows_sampling_envelope);
  HS_RUN(test_preset_transition_adopts_next_params);
  HS_RUN(test_preset_transition_fades_through_black_in_16_frames);
  return fixture.result();
}

} // namespace shapeshifter_oracle_tests
} // namespace hs_test
