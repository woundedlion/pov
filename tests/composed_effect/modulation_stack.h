/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Ordered modulation-stack kernel equivalence tests.
 */

#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

/** @brief A modulation layer with editable sources observes the affine distance. */
inline void test_modulation_stack_observes_warp_distance() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  using Product = FX::RenderPipeline;
  using EditableLayer = Stage::Modulate<Color::Modulate<
      ModulationProvider<FX::Binding, "modulate2", Pattern::SourceKind::NOISE,
                         Pattern::Evaluation::DIRECT, true>>>;
  using Observer = decltype([]<size_t... I>(std::index_sequence<I...>) {
    return Pullback::Pipeline<FX::Binding, Product::StageAt<I>...,
                              EditableLayer>{};
  }(std::make_index_sequence<Product::STAGE_COUNT - 1>{}));
  static_assert(Observer::MODULATION_STACK &&
                !Observer::ELIDES_FINAL_WARP_PATH);
  reset_globals();
  FX effect;
  effect.init();
  auto state = effect.frame_for_test();
  state.resources.template get<"outer_warp">().phase = 0.125f;
  state.params.template get<"modulate">().modulation.enabled = false;
  state.params.template get<"modulate1">().modulation.enabled = false;
  auto &mask = state.params.template get<"modulate2">().modulation;
  mask.source = Pattern::SourceKind::DISPLACEMENT;
  mask.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  mask.enabled = true;
  const auto frame = Observer::prepare(state);
  const auto column = Observer::bind_modulation_stack(frame);
  const auto scan =
      Observer::scan_evaluator([&]() -> const auto & { return column; });
  mask.enabled = false;
  const auto unmasked = Observer::prepare(state);
  bool observed_distance = false;
  for (int x = 0; x < SMALL_W; x += 3) {
    const auto view = math::pixel_to_vector<SMALL_W, SMALL_H>(x, SMALL_H / 3);
    const auto expected = Observer::evaluate(view, frame.ctx, frame.prepared);
    const auto actual = Observer::shade(view, frame, column);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_EQ(actual.alpha, expected.alpha);
    for (int repeat = 0; repeat < 2; ++repeat) {
      const auto cached = Observer::shade(view, frame, scan());
      HS_EXPECT_EQ(cached.color, expected.color);
      HS_EXPECT_EQ(cached.alpha, expected.alpha);
    }
    const auto base = Observer::evaluate(view, unmasked.ctx, unmasked.prepared);
    observed_distance |= actual.alpha > 0.5f * base.alpha;
  }
  HS_EXPECT_TRUE(observed_distance);
}

/** @brief Ordered mechanics share direct samples only for equivalent scalar fields. */
inline void test_modulation_stack_direct_sample_counts() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = FX::RenderPipeline::prepare(effect.frame_for_test());
  constexpr size_t PALETTE = FX::RenderPipeline::NODE_COUNT - 4;
  auto palette = std::get<PALETTE>(frame.prepared);
  palette.previous = palette.current;
  auto source = *std::get<PALETTE + 1>(frame.prepared).current;
  static int samples;
  source.evaluate = [](const Modulation::Prepared &field,
                       const Pattern::SphereDirection &) {
    ++samples;
    return field.config.seed % 2 == 0 ? -0.8f : 0.35f;
  };
  const auto equivalent = source;
  auto distinct = source;
  ++distinct.config.seed;
  const FieldSample field{0.47f, 0.61f, math::Z_AXIS, 0.0f};
  Color::ModulationState hue{&source,
                             &equivalent,
                             Pattern::SourceKind::NOISE,
                             {ColorModulation::Mechanic::HUE_ROTATION, 0.21f},
                             true};
  auto lift = hue;
  lift.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.45f};
  auto mask = hue;
  mask.mechanic = {ColorModulation::Mechanic::OPACITY, 0.65f};
  const std::array<const Color::ModulationState *, 3> layers{&hue, &lift,
                                                             &mask};
  for (int change = 0; change < 3; ++change) {
    if (change == 1)
      mask.current = &distinct;
    if (change == 2)
      mask.enabled = false;
    for (float fraction : {0.0f, 0.37f, 1.0f}) {
      palette.fraction = fraction;
      auto expected = Color::apply_generated_palette(field, palette.current);
      auto previous = Color::apply_generated_palette(field, palette.previous);
      for (const auto *layer : layers) {
        Color::apply_modulation_in_place(previous, *layer, true);
        Color::apply_modulation_in_place(expected, *layer, false);
      }
      if (fraction <= 0.0f)
        expected = previous;
      else if (fraction < 1.0f)
        blend_color_samples_in_place(previous, expected, fraction);
      const Color4 reference = finish_color_sample(expected);
      samples = 0;
      const Color::ModulationStack<3> stack(palette, true, layers);
      const auto actual = stack.shade(field);
      HS_EXPECT_EQ(samples, change == 1 && fraction > 0.0f ? 2 : 1);
      HS_EXPECT_EQ(actual.color, reference.color);
      HS_EXPECT_EQ(actual.alpha, reference.alpha);
    }
  }
}

/** @brief Only endpoints read by an active spatial mechanic require a source. */
inline void test_modulation_stack_required_endpoints() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = FX::RenderPipeline::prepare(effect.frame_for_test());
  constexpr size_t PALETTE = FX::RenderPipeline::NODE_COUNT - 4;
  const auto base = std::get<PALETTE>(frame.prepared).current;

  std::array<float, Pattern::SphereLutView::SIZE> linear;
  linear.fill(32.0f);
  Modulation::Prepared cube{}, direct{};
  cube.cube = {nullptr, true, linear.data()};
  cube.evaluate = Modulation::sample_cube;
  direct.evaluate = [](const Modulation::Prepared &,
                       const Pattern::SphereDirection &) { return -0.5f; };
  const Color::ModulationState disabled{
      nullptr,
      nullptr,
      Pattern::SourceKind::NOISE,
      {ColorModulation::Mechanic::HUE_ROTATION, 0.5f},
      false};
  const Color::ModulationState path{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.2f},
      true};
  const FieldSample field{0.47f, 0.61f, math::Z_AXIS, 0.3f};

  for (float fraction : {0.0f, 0.37f, 1.0f}) {
    Color::TemporalGeneratedPaletteState palette{base, base, fraction};
    Color::ModulationState layer{&direct,
                                 &cube,
                                 Pattern::SourceKind::NOISE,
                                 {ColorModulation::Mechanic::OPACITY, 0.65f},
                                 true};
    if (fraction == 0.0f) {
      palette.current.palette = nullptr;
      layer.current = nullptr;
    }
    if (fraction == 1.0f) {
      palette.previous.palette = nullptr;
      layer.previous = nullptr;
    }
    const std::array<const Color::ModulationState *, 3> layers{&disabled, &path,
                                                               &layer};
    const auto endpoint = [&](bool previous) {
      auto sample = Color::apply_generated_palette(
          field, previous ? palette.previous : palette.current);
      for (const auto *op : layers)
        Color::apply_modulation_in_place(sample, *op, previous);
      return sample;
    };
    auto expected = endpoint(fraction == 0.0f);
    if (fraction > 0.0f && fraction < 1.0f) {
      auto previous = endpoint(true);
      blend_color_samples_in_place(previous, expected, fraction);
    }
    const auto reference = finish_color_sample(expected);
    const auto actual =
        Color::ModulationStack<3>(palette, true, layers).shade(field);
    HS_EXPECT_EQ(actual.color, reference.color);
    HS_EXPECT_EQ(actual.alpha, reference.alpha);
  }
}

/** @brief Independent byte cubes can differ at any one corner of a shared tap. */
inline void test_modulation_stack_cube_corners() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = FX::RenderPipeline::prepare(effect.frame_for_test());
  constexpr size_t PALETTE = FX::RenderPipeline::NODE_COUNT - 4;
  auto palette = std::get<PALETTE>(frame.prepared);
  palette.previous = palette.current;
  palette.fraction = 0.37f;
  auto lift = std::get<PALETTE + 2>(frame.prepared);
  auto mask = std::get<PALETTE + 3>(frame.prepared);
  auto retained = *mask.current;
  auto current = retained;
  std::array<int8_t, Pattern::SphereLutView::SIZE> first{}, second{};
  retained.cube = {first.data(), true};
  current.cube = {second.data(), true};
  lift.previous = lift.current = mask.previous = &retained;
  mask.current = &current;
  lift.paired_cube_ready = mask.paired_cube_ready = true;
  const std::array<const Color::ModulationState *, 2> layers{&lift, &mask};
  const FieldSample field{0.47f, 0.61f, math::Z_AXIS, 0.0f};
  const auto tap = Pattern::prepare_sphere_lut_sample(
      Pattern::SphereDirection::from_unit(field.sphere));
  constexpr int STRIDE = Pattern::SphereLutView::FACE_STEPS;
  for (int corner : {0, 1, STRIDE, STRIDE + 1}) {
    second.fill(0);
    second[tap.offset + corner] = 100;
    const auto initial = Color::apply_generated_palette(field, palette.current);
    auto a = initial;
    auto b = initial;
    for (const auto *layer : layers) {
      Color::apply_modulation_in_place(a, *layer, true);
      Color::apply_modulation_in_place(b, *layer, false);
    }
    HS_EXPECT_TRUE(a.color.alpha != b.color.alpha);
    blend_color_samples_in_place(a, b, palette.fraction);
    const auto expected = finish_color_sample(b);
    const Color::ModulationStack<2> stack(palette, true, layers);
    const auto actual = stack.shade(field);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_EQ(actual.alpha, expected.alpha);
  }

  // A shared cube can feed a signed mechanic even when its peer only feeds lift.
  first.fill(-110);
  second.fill(-110);
  lift.current = &retained;
  lift.previous = &current;
  auto hue = mask;
  hue.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.41f};
  hue.current = hue.previous = &current;
  const std::array<const Color::ModulationState *, 3> inputs{&lift, &mask,
                                                             &hue};
  std::array<unsigned, 3> order{0, 1, 2};
  do {
    const std::array<const Color::ModulationState *, 3> ordered{
        inputs[order[0]], inputs[order[1]], inputs[order[2]]};
    auto previous = Color::apply_generated_palette(field, palette.previous);
    auto next = Color::apply_generated_palette(field, palette.current);
    for (const auto *layer : ordered) {
      Color::apply_modulation_in_place(previous, *layer, true);
      Color::apply_modulation_in_place(next, *layer, false);
    }
    blend_color_samples_in_place(previous, next, palette.fraction);
    const auto expected = finish_color_sample(next);
    const auto actual =
        Color::ModulationStack<3>(palette, true, ordered).shade(field);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_EQ(actual.alpha, expected.alpha);
  } while (std::next_permutation(order.begin(), order.end()));

  // Repeated mechanics must accumulate each endpoint before temporal blending.
  std::array<int8_t, Pattern::SphereLutView::SIZE> third{}, fourth{};
  third.fill(85);
  fourth.fill(-45);
  first.fill(105);
  second.fill(-80);
  auto other_retained = retained, other_current = current;
  other_retained.cube = {third.data(), true};
  other_current.cube = {fourth.data(), true};
  lift.previous = &retained;
  lift.current = &current;
  mask.previous = &other_retained;
  mask.current = &other_current;
  const auto compare_repeated = [&]<ColorModulation::Mechanic M>() {
    constexpr Color::LayerLayout kind{Color::LayerSampling::BYTE_CUBE, M};
    lift.mechanic = {M, 0.6f};
    mask.mechanic = {M, 0.4f};
    const auto stack =
        Color::ModulationStack<2, kind, kind>(palette, true, layers);
    auto previous = Color::apply_generated_palette(field, palette.previous);
    auto next = Color::apply_generated_palette(field, palette.current);
    for (const auto *layer : layers) {
      Color::apply_modulation_in_place(previous, *layer, true);
      Color::apply_modulation_in_place(next, *layer, false);
    }
    blend_color_samples_in_place(previous, next, palette.fraction);
    const auto expected = finish_color_sample(next);
    const auto actual = stack.shade(field);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_EQ(actual.alpha, expected.alpha);
  };
  compare_repeated
      .template operator()<ColorModulation::Mechanic::HUE_ROTATION>();
  compare_repeated.template
  operator()<ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT>();
  compare_repeated.template operator()<ColorModulation::Mechanic::OPACITY>();
}

/** @brief Direct and displacement layouts match the ordered kernel at both endpoints. */
inline void test_modulation_stack_direct_or_path() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  using M = ColorModulation::Mechanic;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = FX::RenderPipeline::prepare(effect.frame_for_test());
  constexpr size_t PALETTE = FX::RenderPipeline::NODE_COUNT - 4;
  auto palette = std::get<PALETTE>(frame.prepared);
  palette.previous = palette.current;
  auto direct = std::get<PALETTE + 1>(frame.prepared);
  auto cube = std::get<PALETTE + 2>(frame.prepared);
  auto previous = *cube.current;
  auto current = previous;
  std::array<int8_t, Pattern::SphereLutView::SIZE> a{}, b{};
  a.fill(-95);
  b.fill(61);
  previous.cube = {a.data(), true};
  current.cube = {b.data(), true};
  cube.previous = &previous;
  cube.current = &current;
  auto compare = [&]<M First, M Second, bool PathFirst>() {
    constexpr Color::LayerLayout mixed{Color::LayerSampling::DIRECT_OR_PATH,
                                       First};
    constexpr Color::LayerLayout baked{Color::LayerSampling::BYTE_CUBE, Second};
    const std::array<const Color::ModulationState *, 2> layers =
        PathFirst
            ? std::array<const Color::ModulationState *, 2>{&direct, &cube}
            : std::array<const Color::ModulationState *, 2>{&cube, &direct};
    direct.mechanic = {First, First == M::HUE_ROTATION ? 2.0f : 0.75f};
    cube.mechanic = {Second, 0.63f};
    for (const auto source :
         {Pattern::SourceKind::NOISE, Pattern::SourceKind::DISPLACEMENT}) {
      direct.source = source;
      for (float fraction : {0.0f, 0.43f, 1.0f}) {
        palette.fraction = fraction;
        const Color::ModulationStack<2, PathFirst ? mixed : baked,
                                     PathFirst ? baked : mixed>
            typed(palette, true, layers);
        const Color::ModulationStack<2> ordered(palette, true, layers);
        for (float path : {0.0f, 0.25f, 1.0f, 2.8f}) {
          const FieldSample field{0.47f, 0.61f, math::Z_AXIS, path};
          const auto expected = ordered.shade(field);
          const auto actual = typed.shade(field);
          HS_EXPECT_EQ(actual.color, expected.color);
          HS_EXPECT_EQ(actual.alpha, expected.alpha);
        }
      }
    }
  };
  compare
      .template operator()<M::HUE_ROTATION, M::POSITIVE_LIGHTNESS_LIFT, true>();
  compare.template
  operator()<M::POSITIVE_LIGHTNESS_LIFT, M::POSITIVE_LIGHTNESS_LIFT, false>();
  compare.template operator()<M::OPACITY, M::OPACITY, true>();
}

/** @brief Typed spatial layouts preserve the opacity blend with displacement. */
inline void test_modulation_stack_spatial_path_mix() {
  using namespace Pullback;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  using M = ColorModulation::Mechanic;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = FX::RenderPipeline::prepare(effect.frame_for_test());
  constexpr size_t PALETTE = FX::RenderPipeline::NODE_COUNT - 4;
  auto palette = std::get<PALETTE>(frame.prepared);
  palette.previous = palette.current;
  palette.fraction = 1.0f;
  auto direct = std::get<PALETTE + 1>(frame.prepared);
  auto cube = std::get<PALETTE + 2>(frame.prepared);
  auto prepared = *cube.current;
  std::array<int8_t, Pattern::SphereLutView::SIZE> bytes{};
  bytes.fill(61);
  prepared.cube = {bytes.data(), true};
  cube.previous = cube.current = &prepared;

  const auto compare =
      [&]<Color::LayerSampling Sampling>(Color::ModulationState layer) {
        constexpr Color::LayerLayout layout{Sampling, M::OPACITY};
        layer.enabled = true;
        layer.mechanic = {M::OPACITY, 1.0f};
        layer.path_mix = {0.6f, 0.8f, 0.7f, false};
        const Color::ModulationStack<1, layout> typed(palette, true, {&layer});
        const Color::ModulationStack<1> ordered(palette, true, {&layer});
        const FieldSample start{0.47f, 0.61f, math::Z_AXIS, 0.0f};
        const FieldSample end{0.47f, 0.61f, math::Z_AXIS, 1.0f};
        const auto expected_start = ordered.shade(start);
        const auto expected_end = ordered.shade(end);
        HS_EXPECT_GT(expected_end.alpha, expected_start.alpha);
        const auto actual_start = typed.shade(start);
        const auto actual_end = typed.shade(end);
        HS_EXPECT_EQ(actual_start.color, expected_start.color);
        HS_EXPECT_EQ(actual_start.alpha, expected_start.alpha);
        HS_EXPECT_EQ(actual_end.color, expected_end.color);
        HS_EXPECT_EQ(actual_end.alpha, expected_end.alpha);
      };
  compare.template operator()<Color::LayerSampling::DIRECT>(direct);
  compare.template operator()<Color::LayerSampling::BYTE_CUBE>(cube);
}

} // namespace hs_test::composed_effect_tests
