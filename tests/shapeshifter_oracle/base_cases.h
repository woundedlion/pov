/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Shapeshifter oracle module cases. */

#include "tests/shapeshifter_oracle/context.h"

namespace hs_test::shapeshifter_oracle_tests {

/** @brief Odd contour tables fit their budget from an unaligned arena offset. */
inline void test_contour_layout_storage_alignment() {
  using Layout = Plot::ContourLayout<3>;
  alignas(Layout::RadiusTrig) std::array<uint8_t, Layout::storage_bytes()>
      storage{};
  Arena arena(storage.data(), storage.size());
  arena.allocate(1, alignof(uint8_t));
  Layout layout;
  layout.init_storage(arena);
  HS_EXPECT_LE(arena.get_offset(), Layout::storage_bytes());
  layout.prepare(3, false);
  HS_EXPECT_EQ(Layout::folded_draw_index(3, 0), 1);
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

/** @brief Candidate matrix stays within visual budget. */
inline void test_candidate_matrix_stays_within_visual_budget() {
  struct Budget {
    double mae;
    double rmse;
    uint32_t channel;
    size_t pixels;
  };
  // Per-case native baselines with aggregate and peak headroom.
  constexpr Budget BUDGETS[] = {
      {0.4, 5, 480, 1},         {2.8, 36, 3800, 1},   {9, 87, 9700, 9},
      {30, 172, 10910, 17},     {0.3, 1.5, 90, 1},    {0.4, 3.5, 210, 1},
      {6.2, 82, 7720, 7},       {54, 232, 11340, 44}, {0.4, 6, 650, 1},
      {2.9, 52, 5010, 1},       {9.7, 102, 6190, 7},  {188, 806, 23780, 541},
      {1.7, 15, 960, 1},        {12, 106, 3810, 1},   {39, 294, 10860, 111},
      {351, 1192, 20450, 1845}, {1.2, 14, 1250, 1},   {9.5, 101, 4550, 1},
      {23, 185, 11160, 25},     {15, 135, 8530, 9}};
  const auto matrix = shape_function_matrix();
  static_assert(std::size(BUDGETS) == std::tuple_size_v<decltype(matrix)>);
  for (size_t index = 0; index < matrix.size(); ++index) {
    const auto &budget = BUDGETS[index];
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
  const OracleFrame planar = render_oracle_frame(state, candidate_renderer());
  state.shape = OracleEffect::ShapeType::SPHERICAL_STAR;
  const OracleFrame spherical =
      render_oracle_frame(state, candidate_renderer());
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
    const OracleFrame uniform =
        render_oracle_frame(state, candidate_renderer());
    state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
    const OracleFrame balanced =
        render_oracle_frame(state, candidate_renderer());
    HS_EXPECT_GT(
        compare_buffers(uniform.pixels, balanced.pixels).different_pixels,
        size_t{100});
  }
}

/** @brief The roster has five shapes and no shipping preset selects the
 * spherical star. */
inline void test_shipping_presets_exclude_spherical_star() {
  HS_EXPECT_EQ(ShapeShifterWhiteBox::shape_count(), 5);
  for (size_t index = 0; index < ShapeShifterWhiteBox::PRESET_COUNT; ++index)
    HS_EXPECT_NE(ShapeShifterWhiteBox::preset_shape(index),
                 OracleEffect::ShapeType::SPHERICAL_STAR);
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
  const std::array<float, 6> PHASES = {
      {0.0f, 0.125f, 0.249f, 0.5f, 0.75f, 0.999f}};
  for (size_t i = 0; i < orientations.size(); ++i) {
    OracleState state;
    state.shape = Shape::PLANAR_STAR;
    state.function = Function::SINE;
    state.count = 144;
    state.sides = 7;
    state.phase = PHASES[i];
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
  const OracleState state = dense_star_state();
  expect_candidate_within_visual_budget(
      state, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT, MAX_STAR_HIGH_ERROR_PIXELS,
      MAX_MEAN_ABSOLUTE_ERROR, MAX_ROOT_MEAN_SQUARED_ERROR, MAX_CHANNEL_ERROR);
}

/** @brief Frozen temporal star strokes preserve the still stroke. */
inline void test_timed_dense_star_stroke_stays_within_visual_budget() {
  const OracleState state = dense_star_state();
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
      "  frozen temporal star stroke: MAE=%.3f RMSE=%.3f energy_drift=%.6f reference=%f candidate=%f\n",
      comparison.error.mean_absolute_error(),
      comparison.error.root_mean_squared_error(), drift, double(energy),
      double(frame_energy(comparison.candidate)));
  HS_EXPECT_LT(comparison.error.mean_absolute_error(), MAX_MEAN_ABSOLUTE_ERROR);
  HS_EXPECT_LT(comparison.error.root_mean_squared_error(),
               MAX_ROOT_MEAN_SQUARED_ERROR);
  HS_EXPECT_LT(drift, MAX_HIGH_COUNT_STAR_ENERGY_DRIFT);
  const auto pole_reference = row_energy(comparison.reference, 0);
  const auto pole_candidate = row_energy(comparison.candidate, 0);
  HS_EXPECT_GT(pole_reference, uint64_t{0});
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
      state.force_contours = true;
      state.shape = OracleEffect::ShapeType::PLANAR_STAR;
      state.function = OracleEffect::PhaseFunction::SINE;
      state.spacing = spacing;
      state.count = count;
      state.sides = 7;
      state.phase = 0.125f;
      state.alpha = 0.274f;
      state.orientation = orientation;
      const OracleFrame frame =
          render_oracle_frame(state, candidate_renderer());
      const math::Vector near = math::rotate(math::X_AXIS, orientation);
      const math::Vector far = -near;
      HS_EXPECT_GE(minimum_energy_around(frame, near, 1.5f), MIN_CAP_ENERGY);
      HS_EXPECT_GE(minimum_energy_around(frame, far, 1.5f), MIN_CAP_ENERGY);
    }
  }
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
    const OracleFrame frame = render_oracle_frame(state, candidate_renderer());
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

/** @brief Opposite flips contour phase across the equator and star phase across hemispheres. */
inline void test_opposite_halves_direction() {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    HS_EXPECT_EQ(effect.updateParameter("Opposite", 0.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_contour(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_contour(effect, 1.5f), -1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_star(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_star(effect, 1.5f), 1.0f);
    HS_EXPECT_EQ(effect.updateParameter("Opposite", 1.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_contour(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_contour(effect, 1.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_star(effect, 0.5f), 1.0f);
    HS_EXPECT_EQ(ShapeShifterWhiteBox::phase_sign_star(effect, 1.5f), -1.0f);
  }
  Timeline().clear();

  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.function = OracleEffect::PhaseFunction::SAWTOOTH;
  state.count = 6;
  state.sides = 7;
  state.phase = 0.23f;
  const OracleFrame unchecked =
      render_oracle_frame(state, candidate_renderer());
  state.opposite = true;
  const OracleFrame checked = render_oracle_frame(state, candidate_renderer());
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
  const float DENSITY_SCALE = COUNT * math::PI_F / (2.0f * DENSITY_INTEGRAL);
  HS_EXPECT_NEAR(equator_equivalent_count, DENSITY_SCALE, 0.2f);
  HS_EXPECT_NEAR(equator_equivalent_count, OracleEffect::MAX_SHAPES, 1.8f);
  HS_EXPECT_NEAR(pole_equivalent_count, DENSITY_SCALE * DENSITY_FLOOR, 0.2f);
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
      return required_param_value(effect, name);
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
 * @brief A preset change fades out over half the segue, adopts the next
 *        preset in the dark, and fades back in over the other half.
 */
inline void test_preset_transition_fades_through_black() {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    Canvas canvas(effect);
    constexpr int HALF = ShapeShifterWhiteBox::preset_segue_frames() / 2;
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_dwell_frames() +
                     ShapeShifterWhiteBox::preset_segue_frames(),
                 ShapeShifterWhiteBox::preset_frames());
    const auto first = ShapeShifterWhiteBox::preset_shape(0);
    const auto second = ShapeShifterWhiteBox::preset_shape(1);
    HS_EXPECT_NE(static_cast<int>(first), static_cast<int>(second));
    auto step = [&] {
      ShapeShifterWhiteBox::step_choreography(effect, canvas);
    };

    for (int frame = 0; frame < ShapeShifterWhiteBox::preset_dwell_frames();
         ++frame)
      step();
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 1.0f, 1e-6f);
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(first));

    for (int frame = 1; frame < HALF; ++frame) {
      step();
      HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect),
                     1.0f - static_cast<float>(frame) / HALF, 1e-6f);
      HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                   static_cast<int>(first));
    }

    step();
    HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect), 0.0f, 1e-6f);
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(second));
    HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{1});

    for (int frame = 1; frame <= HALF; ++frame) {
      step();
      HS_EXPECT_NEAR(ShapeShifterWhiteBox::preset_opacity(effect),
                     static_cast<float>(frame) / HALF, 1e-6f);
    }
    HS_EXPECT_EQ(static_cast<int>(ShapeShifterWhiteBox::shape(effect)),
                 static_cast<int>(second));
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
    HS_EXPECT_EQ(count->max,
                 static_cast<float>(ShapeShifter<96, 20>::DRAW_LIMIT));
    HS_EXPECT_EQ(effect.updateParameter("Count", 288.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(count->get(),
                 static_cast<float>(ShapeShifter<96, 20>::DRAW_LIMIT));
  }
  Timeline().clear();
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
  const Case CASES[] = {
      {ShapeType::SPHERICAL_POLYGON, 24, 4, 2.0f, 0.01f, "spherical polygon"},
      {ShapeType::SPHERICAL_STAR, 12, 5, 1.0f, 0.016f, "spherical star"},
      {ShapeType::PLANAR_STAR, 48, 7, 1.0f, 0.016f, "planar star"},
      {ShapeType::PLANAR_STAR, 12, 5, 1.0f, 0.016f, "sparse planar star"},
      {ShapeType::PLANAR_POLYGON, 16, 5, 1.0f, 0.016f, "planar polygon"},
      {ShapeType::FLOWER, 24, 3, 1.0f, 0.0186f, "flower"}};
  for (const Case &test : CASES) {
    OracleState state;
    state.shape = test.shape;
    state.count = test.count;
    state.sides = test.sides;
    state.amplitude = test.amplitude;
    state.phase = 0.37f;
    uint32_t walked = 0;
    render_timed_contours(state, {.step = test.step,
                                  .timed = true,
                                  .fraction = 1.0f,
                                  .displacement = 0.0f,
                                  .walked = &walked});
    const bool geodesic = test.shape == ShapeType::SPHERICAL_POLYGON ||
                          test.shape == ShapeType::SPHERICAL_STAR;
    HS_CONTEXT("shape", static_cast<int>(test.shape));
    if (geodesic) {
      HS_EXPECT_GT(walked, 0u);
    } else {
      HS_EXPECT_GT(Spacetime::geometry_counters.chord_intervals, 0u);
      HS_EXPECT_EQ(walked, 0u);
    }
    hs_test::seam_compare_render<ORACLE_W, ORACLE_H>(
        [&](bool timed, float fraction, float displacement) {
          return render_timed_contours(state, {.step = test.step,
                                               .timed = timed,
                                               .fraction = fraction,
                                               .displacement = displacement});
        },
        test.name);
  }
}

/** @brief Preset 5's polar-walked strokes keep every lit pixel continuous
 *         across a small phase shift. */
inline void test_timed_preset_five_preserves_polar_anchor() {
  constexpr size_t PRESET = 4;
  constexpr float LIT = 0.025f;
  OracleState state = preset_state(PRESET);
  state.alpha = 1.0f;
  state.phase = 0.002375f;
  const TimedContourRenderOptions options{.step = preset_phase_step(PRESET)};
  const auto before = render_timed_contours(state, options);
  HS_EXPECT_GT(Spacetime::geometry_counters.polar_column_visits, 0u);
  state.phase += 0.00002f;
  const auto after = render_timed_contours(state, options);
  size_t lit = 0;
  for (size_t pixel = 0; pixel < before.size(); ++pixel) {
    if (before[pixel] <= LIT)
      continue;
    ++lit;
    HS_CONTEXT("pixel", pixel);
    HS_EXPECT_NEAR(after[pixel], before[pixel], 0.002f);
  }
  HS_EXPECT_GT(lit, 0u);
}

/** @brief Contour-route presets match their frozen columns at the authored speed.
 * @details Flower contours use L1 error relative to mean frame energy. */
inline void test_timed_presets_match_frozen_columns() {
  constexpr double MAX_FLOWER_RELATIVE_ERROR = 0.27;
  size_t contour_presets = 0;
  for (size_t i = 0; i < ShapeShifterWhiteBox::PRESET_COUNT; ++i) {
    if (ShapeShifterWhiteBox::preset_is_surface(i))
      continue;
    ++contour_presets;
    OracleState state = preset_state(i);
    const float step = preset_phase_step(i);
    const std::string name = "preset " + std::to_string(i + 1);
    HS_CONTEXT("preset", static_cast<int>(i + 1));
    hs_test::seam_compare_render<ORACLE_W, ORACLE_H>(
        [&](bool timed, float fraction, float displacement) {
          return render_timed_contours(state, {.step = step,
                                               .timed = timed,
                                               .fraction = fraction,
                                               .displacement = displacement});
        },
        name.c_str());
    state.orientation =
        math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized();
    const auto frozen_style = compare_renders(
        state, candidate_renderer(), [](OracleEffect &effect, Canvas &canvas) {
          const ScopedTemporalSampling sampling(true);
          ShapeShifterWhiteBox::render_candidate(effect, canvas);
        });
    const double reference_energy =
        double(frame_energy(frozen_style.reference));
    const double candidate_energy =
        double(frame_energy(frozen_style.candidate));
    const double drift =
        (candidate_energy - reference_energy) / reference_energy;
    const double mean_energy = (reference_energy + candidate_energy) * 0.5;
    const double relative_error =
        static_cast<double>(frozen_style.error.total_absolute_error) /
        mean_energy;
    std::printf(
        "  frozen style preset %zu: MAE=%.3f RMSE=%.3f max=%u energy_drift=%.6f relative_error=%.6f\n",
        i + 1, frozen_style.error.mean_absolute_error(),
        frozen_style.error.root_mean_squared_error(),
        frozen_style.error.max_absolute_error, drift, relative_error);
    HS_EXPECT_LT(std::fabs(drift), MAX_PRESET_ENERGY_DRIFT);
    if (state.shape == OracleEffect::ShapeType::FLOWER) {
      HS_EXPECT_LT(relative_error, MAX_FLOWER_RELATIVE_ERROR);
    } else {
      HS_EXPECT_LT(frozen_style.error.mean_absolute_error(),
                   MAX_MEAN_ABSOLUTE_ERROR);
    }
    HS_EXPECT_LT(frozen_style.error.root_mean_squared_error(),
                 MAX_ROOT_MEAN_SQUARED_ERROR);
    HS_EXPECT_LT(frozen_style.error.max_absolute_error, MAX_CHANNEL_ERROR);
  }
  HS_EXPECT_GT(contour_presets, 0u);
}

/** @brief Clipping moving flowers preserves visible pixels and skips chart samples. */
inline void test_timed_flower_clip_culls_anchor_intervals() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::FLOWER;
  state.count = 70;
  state.sides = 3;
  state.phase = 0.37f;
  const auto full = render_timed_contours(
      state,
      {.step = 0.0186f, .timed = true, .fraction = 1.0f, .displacement = 0.0f});
  const auto full_counts = Spacetime::geometry_counters;
  state.clip = {0, ORACLE_H / 2, 0, ORACLE_W / 2};
  const auto clipped = render_timed_contours(
      state,
      {.step = 0.0186f, .timed = true, .fraction = 1.0f, .displacement = 0.0f});
  const auto clipped_counts = Spacetime::geometry_counters;
  HS_EXPECT_GT(clipped_counts.chord_fast_samples, 0u);
  HS_EXPECT_GT(clipped_counts.chord_intervals_culled, 0u);
  HS_EXPECT_GT(clipped_counts.polar_column_samples, 0u);
  HS_EXPECT_LT(clipped_counts.polar_column_samples,
               full_counts.polar_column_samples);
  for (int y = 0; y < ORACLE_H / 2; ++y)
    for (int x = 0; x < ORACLE_W / 2; ++x)
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
}

/** @brief Mid-arm clips preserve fast flower motion through clip corners. */
inline void test_timed_flower_mid_arm_clip_preserves_motion() {
  using Chords = Plot::PlanarChords<ORACLE_W, ORACLE_H>;
  const ClipRegion clip{ORACLE_H / 4, 3 * ORACLE_H / 4, ORACLE_W / 4,
                        3 * ORACLE_W / 4};
  const auto band = Plot::ClipBand<ORACLE_W, ORACLE_H>::of(clip);
  const float left = band.x_center - band.x_half;
  const math::PixelCoords from0{left + 3.0f, band.row_lo - 1.1f};
  const math::PixelCoords from1{left - 2.1f, band.row_lo + 3.0f};
  const math::PixelCoords to0{from0.x + 0.01f, from0.y};
  const math::PixelCoords to1{from1.x + 0.01f, from1.y};
  const auto [first, last] =
      Chords::sample_swept_band_range(band, 1.0f, from0, to0, from1, to1, {});
  int visible = 0;
  for (int frame = 0; frame <= 100; ++frame) {
    const float t = frame / 100.0f;
    const float x = from0.x + (from1.x - from0.x) * t;
    const float y = from0.y + (from1.y - from0.y) * t;
    if (x >= left && x <= band.x_center + band.x_half && y >= band.row_lo &&
        y <= band.row_hi) {
      ++visible;
      HS_EXPECT_LE(first, 0.0f);
      HS_EXPECT_GE(last, 0.0f);
    }
  }
  HS_EXPECT_GT(visible, 0);
  const float right = band.x_center + band.x_half;
  for (const auto ends :
       {std::pair{math::PixelCoords{left + 3.0f, band.row_lo - 3.0f},
                  math::PixelCoords{left + 3.0f, band.row_hi + 3.0f}},
        std::pair{math::PixelCoords{left - 3.0f, band.row_lo + 3.0f},
                  math::PixelCoords{right + 3.0f, band.row_lo + 3.0f}}}) {
    const auto wide = Chords::sample_swept_band_range(
        band, 1.0f, ends.first, ends.first, ends.second, ends.second, {});
    HS_EXPECT_LE(wide.first, 0.0f);
    HS_EXPECT_GE(wide.second, 0.0f);
  }
  OracleState state;
  state.shape = OracleEffect::ShapeType::FLOWER;
  state.count = 70;
  state.sides = 3;
  for (float step : {0.0186f, 0.045f}) {
    for (float phase : {0.0f, 0.37f, 0.8f, 1.4f}) {
      state.phase = phase;
      state.clip = {0, ORACLE_H, 0, ORACLE_W};
      const auto full = render_timed_contours(state, {.step = step});
      state.clip = {ORACLE_H / 4, 3 * ORACLE_H / 4, ORACLE_W / 4,
                    3 * ORACLE_W / 4};
      const auto clipped = render_timed_contours(state, {.step = step});
      double energy = 0.0;
      double error = 0.0;
      for (int y = state.clip.y0; y < state.clip.y1; ++y)
        for (int x = state.clip.x0; x < state.clip.x1; ++x) {
          const size_t i = static_cast<size_t>(y) * ORACLE_W + x;
          energy += full[i];
          error += std::abs(double(clipped[i]) - full[i]);
        }
      HS_EXPECT_GT(energy, 0.0);
      HS_EXPECT_NEAR(error, 0.0, 1e-5);
    }
  }
}

/** @brief A clipped timed dense star samples its polar columns. */
inline void test_timed_dense_star_clip_samples_polar_columns() {
  OracleState state;
  state.clip = {0, ORACLE_H / 2, 0, ORACLE_W / 2};
  state.phase = .37f;
  state.preset_index = 1;
  state.force_contours = true;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.count = 288;
  state.sides = 7;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  render_timed_contours(
      state,
      {.step = 0.016f, .timed = true, .fraction = 1.0f, .displacement = 0.0f});
  const auto star_counts = Spacetime::geometry_counters;
  HS_EXPECT_GT(star_counts.polar_column_samples, 0u);
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
    Plot::ChartPoint a, b; ///< Chart endpoints, radians.
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
} // namespace hs_test::shapeshifter_oracle_tests
