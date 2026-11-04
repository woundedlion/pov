/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file catalog_policies.h
 * @brief Authored fixed-policy parity with catalog operator models. */

#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/operators/table.h"
#include "tests/test_harness.h"

namespace hs_test::catalog_policy_tests {

namespace Op = Pullback::Interp::Op;

/** @brief Copies matching numeric controls across catalog parameter layouts.
 * @tparam To Destination family.
 * @tparam From Source family.
 * @param to Destination controls.
 * @param from Source controls. */
template <typename To, typename From>
void copy_numeric(To &to, const From &from) {
  for (const auto &target : To::FIELDS)
    for (const auto &source : From::FIELDS)
      if (std::string_view(target.id) == source.id)
        to.*target.member = from.*source.member;
}

/** @brief Frame binding for one independently keyed authored stage.
 * @tparam Family Parameter family.
 * @tparam Kind Resource kind. */
template <typename Family, Pullback::ResourceKind Kind,
          Pullback::ResourceKey Key = "probe">
struct StageFrame {
  using Resource = Pullback::ParameterResource<Key, Family, Kind>;
  using Params = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<Resource>>;
  using Frame = Pullback::FrameIntervalState<Params>;
  using Binding = Pullback::Binding<Frame>;
};

/** @brief Compares every authored warp selector against its catalog model.
 * @tparam Family Authored family.
 * @tparam Model Catalog operator model. */
template <typename Family, typename Model> void check_warp() {
  using Fixture = StageFrame<Family, Pullback::ResourceKind::WARP>;
  using Policy =
      typename Pullback::WarpPolicyFor<Family, typename Fixture::Binding,
                                       "probe", true>::Type;
  using Resources =
      typename Pullback::ComposedDetail::PolicyResources<Policy>::Type;
  static_assert(
      std::is_same_v<Resources, Pullback::ComposedDetail::ResourceList<
                                    typename Fixture::Resource>>);
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"probe">();
  auto &clock = frame.resources.template get<"probe">();
  clock.phase = 0.23f;
  if constexpr (requires { params.strength; })
    params.strength = 0.02f;
  if constexpr (requires { params.scale; })
    params.scale = 0.7f;
  if constexpr (requires { params.frequency; })
    params.frequency = 1.7f;
  if constexpr (requires { params.radial_scale; })
    params.radial_scale = 1.3f;
  if constexpr (requires { params.radial_phase; })
    params.radial_phase = 0.2f;
  if constexpr (requires { params.angular_phase; })
    params.angular_phase = -0.1f;
  if constexpr (requires { params.translation_x; }) {
    params.translation_x = 2.0f;
    params.translation_y = -1.0f;
    params.lattice_period = 3.0f;
    params.rotation_rate = 1.2f;
    clock.rotation = 0.4f;
  }
  if constexpr (requires { params.repeat_cells; })
    params.repeat_cells = 3.0f;
  constexpr int BASIS_COUNT = requires { params.basis; } ? 3 : 1;
  constexpr int ENVELOPE_COUNT = requires { params.envelope; } ? 3 : 1;
  constexpr int MODE_COUNT = requires { params.mode; } ? 2 : 1;
  constexpr int HARMONIC_COUNT = requires { params.harmonic; } ? 16 : 1;
  constexpr int INTEGRATOR_COUNT = requires { params.integrator; } ? 3 : 1;
  for (int basis = 0; basis < BASIS_COUNT; ++basis)
    for (int envelope = 0; envelope < ENVELOPE_COUNT; ++envelope)
      for (int mode = 0; mode < MODE_COUNT; ++mode)
        for (int harmonic = 0; harmonic < HARMONIC_COUNT; ++harmonic)
          for (int integrator = 0; integrator < INTEGRATOR_COUNT;
               ++integrator) {
            typename Model::Params expected_params{};
            copy_numeric(expected_params, params);
            if constexpr (requires { params.basis; })
              params.basis = expected_params.basis =
                  static_cast<uint8_t>(basis);
            if constexpr (requires { params.envelope; })
              params.envelope = expected_params.envelope =
                  static_cast<uint8_t>(envelope);
            if constexpr (requires { params.mode; })
              params.mode = expected_params.mode = static_cast<uint8_t>(mode);
            if constexpr (requires { params.harmonic; })
              params.harmonic = expected_params.harmonic =
                  static_cast<uint8_t>(harmonic);
            if constexpr (requires { params.integrator; })
              params.integrator = expected_params.integrator =
                  static_cast<uint8_t>(integrator);
            typename Model::State state{};
            state.phase = clock.phase;
            if constexpr (requires { state.rotation; })
              state.rotation = clock.rotation;
            if constexpr (requires { state.noise; }) {
              Pullback::init_effect_noise(state.noise, 31);
              clock.noise = &state.noise;
            }
            const auto actual_prepared = Policy::prepare(frame);
            const auto expected_prepared = [&] {
              if constexpr (requires {
                              Model::prepare({}, expected_params, state,
                                             Pullback::TranslationLattice{});
                            })
                return Model::prepare({}, expected_params, state,
                                      Pullback::TranslationLattice{});
              else
                return Model::prepare({}, expected_params, state);
            }();
            for (const math::Complex point :
                 {math::Complex{0.17f, 0.23f}, math::Complex{-0.8f, 0.4f},
                  math::Complex{0x1p21f, 0.2f}}) {
              Pullback::PlaneSample input{point, {}, math::Y_AXIS, 0.3f};
              input.provenance.value_weight = 0.4f;
              input.provenance.fade_edge_distance = 0.03f;
              const auto actual = Policy::apply(input.coords, input.provenance,
                                                frame, actual_prepared);
              const auto expected =
                  Model::run(input, {}, expected_params, expected_prepared);
              HS_EXPECT_NEAR(actual.coords.re, expected.coords.re, 1e-6f);
              HS_EXPECT_NEAR(actual.coords.im, expected.coords.im, 1e-6f);
              HS_EXPECT_NEAR(actual.path_length + input.path_length,
                             expected.path_length, 1e-6f);
            }
          }
}

/** @brief Frame binding of an affine warp keyed "probe" before a lattice
    source keyed "source". */
struct LatticeAffineFrame {
  using Warp = Pullback::ParameterResource<"probe", Pullback::AffineParams,
                                           Pullback::ResourceKind::WARP>;
  using Source =
      Pullback::ParameterResource<"source", Pullback::LatticeSourceParams,
                                  Pullback::ResourceKind::SOURCE>;
  using Params = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<Warp, Source>>;
  using Frame = Pullback::FrameIntervalState<Params>;
  using Binding = Pullback::Binding<Frame>;
};

/** @brief Whole windings of the source lattice preserve repeat-relative
    displacement across phase wraps. */
inline void check_affine_lattice_windings() {
  using namespace Pullback;
  using Fixture = LatticeAffineFrame;
  using Tracked =
      typename WarpPolicyFor<AffineParams, typename Fixture::Binding, "probe",
                             true>::Type;
  using Untracked =
      typename WarpPolicyFor<AffineParams, typename Fixture::Binding, "probe",
                             false>::Type;
  static_assert(Tracked::template PROVIDER_VALID<typename Fixture::Binding>);
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"probe">();
  auto &clock = frame.resources.template get<"probe">();
  frame.params.template get<"source">().lattice_cell_scale = 2.5f;
  params.lattice_period = 7.0f;
  for (const float winding : {-2.51f, -1.5f, -1.49f, 1.49f, 1.5f, 2.51f}) {
    params.translation_x = winding;
    params.translation_y = -winding;
    for (const float phase : {0.0f, 0.125f, 0.5f, 0.875f}) {
      clock.phase = phase;
      const auto prepared = Tracked::prepare(frame);
      const auto actual = Tracked::apply({}, {}, frame, prepared);
      const float expected_x = -phase * std::round(winding) * 0.4f;
      HS_EXPECT_NEAR(actual.coords.re, expected_x, 3e-7f);
      HS_EXPECT_NEAR(actual.coords.im, -expected_x, 3e-7f);
      const float remainder = std::remainder(expected_x, 0.4f);
      HS_EXPECT_NEAR(actual.path_length,
                     std::sqrt(2.0f * remainder * remainder), 3e-7f);
      const auto skipped =
          Untracked::apply({}, {}, frame, Untracked::prepare(frame));
      HS_EXPECT_NEAR(skipped.coords.re, actual.coords.re, 0.0f);
      HS_EXPECT_NEAR(skipped.path_length, 0.0f, 0.0f);
    }
  }
  params.translation_x = 1.7f;
  params.translation_y = -2.2f;
  params.scale_x = 1.3f;
  params.scale_y = 0.7f;
  params.shear = 0.3f;
  params.repeat_cells = 2.0f;
  clock.rotation = 0.7f;
  const math::Complex point{0.13f, -0.21f};
  clock.phase = 0.999999f;
  const auto before = Tracked::apply(point, {}, frame, Tracked::prepare(frame));
  clock.phase = 0.000001f;
  const auto after = Tracked::apply(point, {}, frame, Tracked::prepare(frame));
  HS_EXPECT_NEAR(after.coords.re - before.coords.re, 1.6f, 6e-6f);
  HS_EXPECT_NEAR(after.coords.im - before.coords.im, -1.6f, 6e-6f);
  HS_EXPECT_NEAR(before.path_length, after.path_length, 3e-6f);
  HS_EXPECT_TRUE(valid(params));
  params.lattice_period = 0.0f;
  HS_EXPECT_FALSE(valid(params));
  for (const float direction : {-1.0f, 1.0f}) {
    typename Fixture::Frame from{}, to{};
    auto &a = from.params.template get<"probe">();
    auto &b = to.params.template get<"probe">();
    a.translation_x = 1.0f;
    b.translation_x = -2.0f;
    auto &start = from.resources.template get<"probe">();
    auto &end = to.resources.template get<"probe">();
    start.phase = direction > 0.0f ? 0.999f : 0.001f;
    end.phase = direction > 0.0f ? 0.001f : 0.999f;
    start.rotation = math::TWO_PI_F - 0.001f;
    end.rotation = 0.001f;
    FrameInterval interval(from, to);
    for (const float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
      const auto sampled = interval.at(fraction);
      const auto &sampled_clock = sampled.resources.template get<"probe">();
      HS_EXPECT_NEAR(math::wrap_t(sampled_clock.phase),
                     math::wrap_t(start.phase + direction * 0.002f * fraction),
                     2e-7f);
      const auto expected =
          math::wrap_t((start.rotation + 0.002f * fraction) / math::TWO_PI_F);
      HS_EXPECT_NEAR(math::wrap_t(sampled_clock.rotation / math::TWO_PI_F),
                     expected, 2e-7f);
      HS_EXPECT_NEAR(sampled.params.template get<"probe">().translation_x,
                     fraction < 1.0f ? 1.0f : -2.0f, 0.0f);
    }
  }
}

/** @brief The composed affine stage before a lattice source matches the
    catalog model given that source's lattice. */
inline void check_affine_lattice_parity() {
  using namespace Pullback;
  using Fixture = LatticeAffineFrame;
  using Policy = typename WarpPolicyFor<AffineParams, typename Fixture::Binding,
                                        "probe", true>::Type;
  static_assert(WarpProvider<typename Fixture::Binding, "probe", AffineParams,
                             true>::repeat_reduction() ==
                Warp::RepeatReduction::RECTANGULAR);
  static_assert(
      WarpProvider<
          typename StageFrame<AffineParams, ResourceKind::WARP>::Binding,
          "probe", AffineParams, true>::repeat_reduction() ==
      Warp::RepeatReduction::UNWRAPPED);
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"probe">();
  auto &source = frame.params.template get<"source">();
  auto &clock = frame.resources.template get<"probe">();
  source.lattice_cell_scale = 3.0f;
  params.translation_x = 2.0f;
  params.translation_y = -1.0f;
  params.repeat_cells = 3.0f;
  params.lattice_period = 5.0f;
  params.rotation_rate = 1.2f;
  params.scale_x = 1.4f;
  params.shear = 0.2f;
  clock.phase = 0.23f;
  clock.rotation = 0.4f;
  const Op::WarpAffineV4::State state{clock.phase, clock.rotation};
  const auto expected_prepared = Op::WarpAffineV4::prepare(
      {}, params, state, Source::translation_lattice(source));
  const auto actual_prepared = Policy::prepare(frame);
  for (const math::Complex point :
       {math::Complex{0.17f, 0.23f}, math::Complex{-0.8f, 0.4f},
        math::Complex{3.1f, -2.2f}}) {
    PlaneSample input{point, {}, math::Y_AXIS, 0.0f};
    const auto actual =
        Policy::apply(input.coords, input.provenance, frame, actual_prepared);
    const auto expected =
        Op::WarpAffineV4::run(input, {}, params, expected_prepared);
    HS_EXPECT_NEAR(actual.coords.re, expected.coords.re, 0.0f);
    HS_EXPECT_NEAR(actual.coords.im, expected.coords.im, 0.0f);
    HS_EXPECT_NEAR(actual.path_length, expected.path_length, 0.0f);
    HS_EXPECT_LE(actual.path_length, std::sqrt(2.0f) * 0.5f + 1e-6f);
  }
}

/** @brief Every planar authored selector matches the catalog's observable step. */
inline void test_authored_warp_catalog_parity() {
  check_warp<Pullback::Warp::AuthoredWaveShearParams, Op::WarpWaveShear>();
  check_warp<Pullback::Warp::AuthoredVectorNoiseParams, Op::WarpVectorNoise>();
  check_warp<Pullback::Warp::AuthoredPolarParams, Op::WarpPolarChart>();
  check_warp<Pullback::Warp::AuthoredCurlFlowParams, Op::WarpCurlFlow>();
  check_warp<Pullback::AffineParams, Op::WarpAffineV4>();
  check_affine_lattice_windings();
  check_affine_lattice_parity();
}

/** @brief Compares every authored source selector against its catalog model.
 * @tparam Family Source family.
 * @tparam Model Catalog operator model. */
template <typename Family, typename Model> void check_source() {
  using Fixture = StageFrame<Family, Pullback::ResourceKind::SOURCE, "source">;
  using Policy =
      typename Pullback::SourcePolicyFor<Family,
                                         typename Fixture::Binding>::Type;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  auto &clock = frame.resources.template get<"source">();
  clock.primary = 0.4f;
  clock.secondary = 0.7f;
  clock.angle = 0.2f;
  clock.noise_time = 0.23f;
  constexpr int BASIS_COUNT = requires { params.basis; } ? 3 : 1;
  constexpr int KIND_COUNT = requires { params.kind; } ? 3 : 1;
  for (int basis = 0; basis < BASIS_COUNT; ++basis)
    for (int kind = 0; kind < KIND_COUNT; ++kind) {
      typename Model::Params expected_params{};
      copy_numeric(expected_params, params);
      if constexpr (requires { expected_params.weight_mode; }) {
        expected_params.weight_mode = 0;
        expected_params.coverage_mode = 0;
      }
      if constexpr (requires { params.basis; })
        params.basis = expected_params.basis = static_cast<uint8_t>(basis);
      if constexpr (requires { params.kind; })
        params.kind = expected_params.kind = static_cast<uint8_t>(kind);
      typename Model::State state{};
      if constexpr (requires { state.primary; }) {
        state.primary = clock.primary;
        state.secondary = clock.secondary;
        state.angle = clock.angle;
      }
      if constexpr (requires { state.noise; }) {
        Pullback::init_effect_noise(state.noise, 31);
        clock.noise = &state.noise;
        state.phase = clock.noise_time;
      }
      const auto actual_prepared = Policy::prepare(frame);
      const auto expected_prepared = Model::prepare({}, expected_params, state);
      for (const math::Complex point :
           {math::Complex{0.17f, 0.23f}, math::Complex{-0.8f, 0.4f},
            math::Complex{0x1p21f, 0.2f}}) {
        typename Model::Input input{};
        if constexpr (requires { input.coords; }) {
          input.coords = point;
          input.sphere = math::Vector(point.re, 0.6f, point.im).normalized();
        } else
          input.dir = math::Vector(point.re, 0.6f, point.im).normalized();
        const float actual = Policy::sample(input, frame, actual_prepared);
        const auto expected =
            Model::run(input, {}, expected_params, expected_prepared);
        HS_EXPECT_NEAR(0.5f + 0.5f * actual, expected.value, 1e-6f);
      }
    }
}

/** @brief Every authored source selector matches the catalog's sampled value. */
inline void test_authored_source_catalog_parity() {
  check_source<Pullback::RingsSourceParams, Op::SampleRings>();
  check_source<Pullback::Source::AuthoredTessellationSourceParams,
               Op::SampleTessellation>();
  check_source<Pullback::Source::AuthoredProjectedNoiseSourceParams,
               Op::SampleProjectedNoise>();
}

/** @brief Runtime signal and coverage selectors preserve the crossing contract. */
inline void test_authored_sample_crossing() {
  using Family = Pullback::Source::AuthoredRingsSampleParams;
  using Fixture = StageFrame<Family, Pullback::ResourceKind::SOURCE, "source">;
  using Binding = typename Fixture::Binding;
  using Provider = Pullback::SourceProvider<Binding, Family>;
  using Source = typename Pullback::SourcePolicyFor<Family, Binding>::Type;
  using Weight = Pullback::Weight::Authored<Provider>;
  using Coverage = Pullback::ProjectionCoverage::Authored<Provider>;
  using Crossing = Pullback::Stage::Sample<Source, Weight, Coverage>;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  frame.resources.template get<"source">().primary = 0.6f;
  params.edge_width = 0.1f;
  Pullback::PlaneSample input{{0.17f, 0.23f}, {}, math::Y_AXIS, 0.3f};
  input.provenance.value_weight = 0.4f;
  input.provenance.domain_coverage = 0.7f;
  constexpr float DISTANCES[] = {-0.1f,  0.0f, 0.025f, 0.05f,
                                 0.075f, 0.1f, 0.2f};
  constexpr float FADES[] = {0.0f, 0.0f, 0.15625f, 0.5f, 0.84375f, 1.0f, 1.0f};
  for (uint8_t weight = 0; weight < 2; ++weight)
    for (uint8_t coverage = 0; coverage < 5; ++coverage)
      for (size_t point = 0; point < std::size(DISTANCES); ++point) {
        params.weight_mode = weight;
        params.coverage_mode = coverage;
        input.provenance.fade_edge_distance = DISTANCES[point];
        const auto prepared = Crossing::template prepare<Binding>(frame);
        const auto actual =
            Crossing::template run<Binding>(input, frame, prepared);
        const float raw = Pattern::rings(
            Pattern::PlanePosition{projections::stereo_pattern_args(
                input.coords, params.pattern_freq)},
            Pattern::prepare(0.6f, 0.0f, 0.0f));
        const float expected_value = 0.5f + 0.5f * raw * (weight ? 0.4f : 1.0f);
        const float expected_coverage[] = {1.0f, 0.4f, 0.16f, FADES[point],
                                           FADES[point] * FADES[point]};
        HS_EXPECT_NEAR(actual.value, expected_value, 1e-6f);
        HS_EXPECT_NEAR(actual.coverage, 0.7f * expected_coverage[coverage],
                       1e-6f);
        HS_EXPECT_NEAR(actual.path_length, 0.3f, 0.0f);
      }
  params.weight_mode = 2;
  HS_EXPECT_FALSE(Pullback::valid(params));
  params.weight_mode = 1;
  params.coverage_mode = 5;
  HS_EXPECT_FALSE(Pullback::valid(params));
  params.coverage_mode = 4;
  HS_EXPECT_TRUE(Pullback::valid(params));
  Pullback::Source::AuthoredGridSampleParams grid{};
  grid.secondary_rate = 1.8f;
  HS_EXPECT_TRUE(Pullback::valid(grid));
}

/** @brief Static squared edge fade squares coverage without changing source weights. */
inline void test_edge_fade_squared_curve() {
  struct Frame {
    float width;
  };
  struct Provider {
    static float edge_width(const Frame &frame) { return frame.width; }
  };
  using Policy = Pullback::ProjectionCoverage::EdgeFadeSquared<Provider>;
  Pullback::ProjectionProvenance provenance{};
  provenance.value_weight = 0.25f;
  for (const float width : {0.0f, 0.1f})
    for (const float distance : {-0.1f, 0.0f, 0.05f, 0.1f, 0.2f}) {
      provenance.fade_edge_distance = distance;
      const float expected = width == 0.0f ? (distance > 0.0f ? 1.0f : 0.0f)
                             : distance <= 0.0f ? 0.0f
                             : distance >= 0.1f ? 1.0f
                                                : 0.25f;
      HS_EXPECT_NEAR(Policy::apply(provenance, Frame{width}), expected, 1e-6f);
    }
}

/** @brief Ridge and cosine bands preserve their authored transfer controls. */
inline void test_authored_transfer_catalog_parity() {
  using Fixture =
      StageFrame<Pullback::SmoothBandsParams, Pullback::ResourceKind::VALUE>;
  using Provider =
      Pullback::ValueProvider<typename Fixture::Binding,
                              Pullback::SmoothBandsParams, "probe">;
  using Policy = Pullback::Transfer::SmoothBands<Provider>;
  typename Fixture::Frame frame{};
  for (const float count : {1.0f, 7.0f, 32.0f})
    for (const float phase : {0.0f, 0.3f, 5.9f}) {
      auto &params = frame.params.template get<"probe">();
      params.band_count = count;
      params.band_phase = phase;
      for (const float value : {0.0f, 0.19f, 0.5f, 0.91f, 1.0f}) {
        const Pullback::FieldSample input{value, 0.7f, math::Y_AXIS, 0.0f};
        const auto expected =
            Op::TransferSmoothBands::run(input, {}, params, {});
        HS_EXPECT_NEAR(Policy::apply(value, frame), expected.value, 1e-6f);
        HS_EXPECT_NEAR(Pullback::Transfer::Ridge::apply(value, frame),
                       Op::TransferRidge::run(input, {}, {}, {}).value, 1e-6f);
      }
    }
}

} // namespace hs_test::catalog_policy_tests
