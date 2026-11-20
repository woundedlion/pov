/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_sources.h
 * @brief Spherical source carriers, composed policies and instance clocks. */

#include "core/render/pullback/sphere_tiling_edges_baked.h"
#include "tests/shader_chain_fixture.h"

namespace hs_test::shader_chain_tests {

/** @brief Composed frame for one spherical source resource.
 * @tparam Family Source parameter family. */
template <typename Family> struct SphericalSourceFrame {
  using Resource =
      PB::ParameterResource<"source", Family, PB::ResourceKind::SOURCE>;
  using Params = PB::ComposedDetail::ParameterSet<
      PB::ComposedDetail::ResourceList<Resource>>;
  using Frame = PB::FrameIntervalState<Params>;
  using Binding = PB::Binding<Frame>;
  using Policy = typename PB::SourcePolicyFor<Family, Binding>::Type;
  using Stage = PB::Stage::SampleSphere<Policy>;
  static_assert(
      std::is_same_v<typename PB::ComposedDetail::PolicyResources<Policy>::Type,
                     PB::ComposedDetail::ResourceList<Resource>>);
};

/** @brief Spherical noise adapters preserve the canonical source and carrier. */
inline void test_shader_chain_spherical_noise_policy() {
  using Model = In::Op::SampleSphericalNoise;
  using Fixture = SphericalSourceFrame<Model::Params>;
  using Binding = Fixture::Binding;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  params.noise_scale = 3.5f;
  params.noise_contrast = 1.7f;
  params.noise_time_rate = -1.0f / 128.0f;
  Model::State state;
  Model::init(state,
              {"source", Model::ID, In::instance_hash("source", Model::ID)});
  state.phase = 0.37f;
  auto &clock = frame.resources.template get<"source">();
  clock.noise = &state.noise;
  clock.noise_time = state.phase;
  const auto oracle_noise = make_oracle_noise("source", Model::ID, 1.0f);
  for (uint8_t basis = 0;
       basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3); ++basis) {
    params.basis = basis;
    const auto prepared = Model::prepare({}, params, state);
    const auto typed = Fixture::Stage::template prepare<Binding>(frame);
    const auto oracle = Pattern::SphereNoise::prepare(
        {3.5f, -1.0f / 128.0f, static_cast<math::NoiseBasis>(basis), 1.7f,
         Pattern::NoiseVariant::OCTAVE_CONTOUR},
        {0.37f}, oracle_noise);
    for (const auto &view : sweep_views()) {
      const float unit_raw = Pattern::SphereNoise::sample(
          Pattern::SphereDirection::from(view), oracle);
      for (float scale : {1.0f, std::sqrt(1.003f)}) {
        const PB::SphereSample input{scale * view, 0.375f};
        const auto actual = Model::run(input, {}, params, prepared);
        const auto composed =
            Fixture::Stage::template run<Binding>(input, frame, typed);
        const float raw = Pattern::SphereNoise::sample(
            Pattern::SphereDirection::from(input.dir), oracle);
        HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * raw, 1.0e-7f);
        HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * unit_raw, 2.0e-5f);
        HS_EXPECT_TRUE(field_identical(actual, composed));
        HS_EXPECT_EQ(actual.coverage, 1.0f);
        HS_EXPECT_EQ(actual.path_length, 0.375f);
        HS_EXPECT_EQ(actual.sphere.x, input.dir.x);
        HS_EXPECT_EQ(actual.sphere.y, input.dir.y);
        HS_EXPECT_EQ(actual.sphere.z, input.dir.z);
      }
    }
  }
}

/** @brief Spherical cell selectors agree between typed and runtime crossings. */
inline void test_shader_chain_spherical_cells_policy() {
  using Model = In::Op::SampleSphericalCells;
  using Fixture = SphericalSourceFrame<Model::Params>;
  using Binding = Fixture::Binding;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  params.count = 17.0f;
  params.seed = 391.0f;
  for (uint8_t layout = 0;
       layout <= static_cast<uint8_t>(Pattern::SphereCellLayout::SCATTERED);
       ++layout) {
    params.layout = layout;
    for (uint8_t output = 0;
         output <= static_cast<uint8_t>(Pattern::SphereCellOutput::CELL_VALUE);
         ++output) {
      params.output = output;
      const auto prepared = Model::prepare({}, params, {});
      const auto typed = Fixture::Stage::template prepare<Binding>(frame);
      for (const auto &view : sweep_views()) {
        const float unit_raw = Pattern::SphereCells::sample(
            Pattern::SphereDirection::from(view), prepared);
        for (float scale : {1.0f, std::sqrt(1.003f)}) {
          const PB::SphereSample input{scale * view, 1.25f};
          const auto actual = Model::run(input, {}, params, prepared);
          const auto composed =
              Fixture::Stage::template run<Binding>(input, frame, typed);
          HS_EXPECT_TRUE(field_identical(actual, composed));
          const float raw = Pattern::SphereCells::sample(
              Pattern::SphereDirection::from(input.dir), prepared);
          HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * raw, 1.0e-7f);
          if (output !=
              static_cast<uint8_t>(Pattern::SphereCellOutput::CELL_VALUE))
            HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * unit_raw, 3.0e-6f);
          HS_EXPECT_TRUE(actual.value >= 0.0f && actual.value <= 1.0f);
          HS_EXPECT_EQ(actual.coverage, 1.0f);
          HS_EXPECT_EQ(actual.path_length, 1.25f);
          HS_EXPECT_EQ(actual.sphere.x, input.dir.x);
          HS_EXPECT_EQ(actual.sphere.y, input.dir.y);
          HS_EXPECT_EQ(actual.sphere.z, input.dir.z);
        }
      }
    }
  }
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  const PB::SphereSample INPUT{math::X_AXIS, 2.0f};
  const auto prepared = Model::prepare({}, params, {});
  HS_EXPECT_NEAR(Model::run(INPUT, {}, params, prepared).value, 1.0f, 1.0e-6f);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::BORDERS);
  params.line_thickness = 0.1f;
  params.line_softness = 0.2f;
  const float angle =
      std::numbers::pi_v<float> / 4.0f -
      0.2f / std::sqrt(6.0f / (4.0f * std::numbers::pi_v<float>));
  const math::Vector border{std::cos(angle), std::sin(angle), 0.0f};
  const auto border_prepared = Model::prepare({}, params, {});
  const auto border_typed = Fixture::Stage::template prepare<Binding>(frame);
  for (float scale : {1.0f, std::sqrt(1.003f)}) {
    const PB::SphereSample drifted{scale * border, 2.0f};
    const auto actual = Model::run(drifted, {}, params, border_prepared);
    const auto composed =
        Fixture::Stage::template run<Binding>(drifted, frame, border_typed);
    HS_EXPECT_NEAR(actual.value, 0.5f, 2.0e-6f);
    HS_EXPECT_NEAR(composed.value, 0.5f, 2.0e-6f);
    HS_EXPECT_EQ(actual.sphere.x, drifted.dir.x);
    HS_EXPECT_EQ(composed.sphere.x, drifted.dir.x);
    HS_EXPECT_EQ(actual.path_length, 2.0f);
    HS_EXPECT_EQ(composed.path_length, 2.0f);
  }
}

/** @brief Baked edge samples match freshly built registry geometry. */
inline void test_shader_chain_baked_tiling_edges() {
  using namespace PB::Source;
  const SphereTilingEdgeBand bands[] = {
      {1.0f, 0.5f}, {0.0f, 1.398054f}, {10.0f, 10.0f}};
  constexpr const char *BAKED_NAMES[] = {
      "truncatedIcosahedron", "truncatedIcosahedron_truncate50d_ambo_dual",
      "dodecahedron_hk54_ambo_hk72"};
  unsigned solid = 0;
  for (const auto registry : Solids::all_registries()) {
    for (const auto &entry : registry) {
      const bool expected = std::ranges::any_of(BAKED_NAMES, [&](auto name) {
        return std::string_view(entry.name) == name;
      });
      const auto baked = baked_sphere_tiling_edges(solid);
      HS_EXPECT_EQ(supports_baked_sphere_tiling_edges(solid), expected);
      HS_EXPECT_EQ(baked.edges.empty(), !expected);
      if (expected) {
        const SphereTilingGeometry geometry{solid};
        const auto runtime = geometry.edge_view();
        HS_EXPECT_EQ(baked.vertices.size(), runtime.vertices.size());
        HS_EXPECT_EQ(baked.edges.size(), runtime.edges.size());
        for (int i = 0; i < 4096; ++i) {
          const auto direction = math::fib_spiral(4096, 0.47f, i);
          for (const auto &band : bands) {
            const auto a = baked.sample(direction, band);
            const auto b = runtime.sample(direction, band);
            HS_EXPECT_NEAR(a.coverage, b.coverage, 1.0e-4f);
            if (b.coverage >= 1.0e-4f)
              HS_EXPECT_NEAR(a.value, b.value, 1.0e-4f);
          }
        }
      }
      ++solid;
    }
  }
}

/** @brief Tiling maps palette position and coverage along cube edges and faces. */
inline void test_shader_chain_spherical_tiling() {
  using Model = In::Op::SampleSphericalTiling;
  using Geometry = PB::Source::SphereTilingGeometry;
  using Output = PB::Source::TilingOutput;
  constexpr float DEG = math::PI_F / 180.0f;
  In::FrameContext context;
  context.angular_sample_radius = 0.125f * DEG;
  const auto face_point = [](float depth) {
    const float angle = (45.0f - depth) * DEG;
    return PB::SphereSample{{std::cos(angle), std::sin(angle), 0.0f}, 0.25f};
  };
  Model::Params params;
  params.solid = static_cast<uint8_t>(Solids::BaseMesh::CUBE);
  params.line_thickness = 2.0f;
  params.line_softness = 2.0f;
  Model::State state;
  const auto prepared = Model::prepare({}, params, state);
  const auto seam = face_point(0.0f);
  const auto center = face_point(45.0f);
  const auto edge_field = Model::run(seam, {}, params, prepared);
  HS_EXPECT_NEAR(edge_field.value, 0.5f, 2.0e-6f);
  HS_EXPECT_EQ(edge_field.coverage, 1.0f);
  HS_EXPECT_EQ(edge_field.path_length, seam.path_length);
  const float half_arc = std::atan(std::sqrt(0.5f));
  for (float position : {0.25f, 0.75f}) {
    const float arc = (2.0f * position - 1.0f) * half_arc;
    const auto quarter =
        std::cos(arc) * seam.dir + std::sin(arc) * math::Z_AXIS;
    const auto inside = Model::run({quarter, 0.0f}, {}, params, prepared);
    HS_EXPECT_NEAR(inside.value, position, 2.0e-6f);
    for (float side : {-1.0f, 1.0f}) {
      const auto across = math::Vector(side, -side, 0.0f).normalized();
      const auto fringe_point =
          std::cos(2.0f * DEG) * quarter + std::sin(2.0f * DEG) * across;
      const auto fringe =
          Model::run({fringe_point, 0.0f}, {}, params, prepared);
      HS_EXPECT_NEAR(fringe.value, position, 2.0e-6f);
      HS_EXPECT_NEAR(fringe.coverage, 0.5f, 2.0e-6f);
    }
  }
  HS_EXPECT_EQ(Model::run(face_point(4.0f), {}, params, prepared).coverage,
               0.0f);

  params.output = static_cast<uint8_t>(Output::FACES);
  params.line_softness = 20.0f / 44.0f;
  const auto faces = Model::prepare({}, params, state);
  HS_EXPECT_TRUE(faces.geometry == prepared.geometry);
  const auto cutout = Model::run(seam, {}, params, faces);
  HS_EXPECT_EQ(cutout.value, 0.0f);
  HS_EXPECT_NEAR(cutout.coverage, 0.0f, 2.0e-6f);
  const auto fringe = Model::run(face_point(2.0f), {}, params, faces);
  HS_EXPECT_NEAR(fringe.value, 1.0f / 44.0f, 2.0e-6f);
  HS_EXPECT_NEAR(fringe.coverage, 0.5f, 2.0e-6f);
  HS_EXPECT_NEAR(Model::run(face_point(1.5f), {}, params, faces).coverage,
                 0.25f, 5.0e-6f);
  HS_EXPECT_NEAR(Model::run(face_point(2.5f), {}, params, faces).coverage,
                 0.75f, 5.0e-6f);
  const auto middle = Model::run(face_point(23.0f), {}, params, faces);
  HS_EXPECT_NEAR(middle.value, 0.5f, 2.0e-6f);
  HS_EXPECT_EQ(middle.coverage, 1.0f);
  const auto face_center = Model::run(center, {}, params, faces);
  HS_EXPECT_NEAR(face_center.value, 1.0f, 2.0e-6f);
  HS_EXPECT_EQ(face_center.coverage, 1.0f);

  params.line_thickness = 0.0f;
  params.line_softness = 5.0f / 45.0f;
  const auto zero_width = Model::prepare({}, params, state);
  const auto seam_field = Model::run(seam, context, params, zero_width);
  HS_EXPECT_NEAR(seam_field.value, 0.0f, 2.0e-6f);
  HS_EXPECT_NEAR(seam_field.coverage, 0.125f, 2.0e-6f);
  const auto zero_width_middle =
      Model::run(face_point(22.5f), context, params, zero_width);
  HS_EXPECT_NEAR(zero_width_middle.value, 0.5f, 2.0e-6f);
  HS_EXPECT_EQ(zero_width_middle.coverage, 1.0f);

  const auto zero_width_fringe =
      Model::run(face_point(0.25f), context, params, zero_width);
  HS_EXPECT_NEAR(zero_width_fringe.coverage, 0.5f, 1.0e-5f);
  params.line_softness = 10.0f / 45.0f;
  const auto softer_zero = Model::prepare({}, params, state);
  HS_EXPECT_NEAR(Model::run(seam, context, params, softer_zero).coverage,
                 0.0625f, 2.0e-6f);

  params.line_thickness = 0.002f;
  params.line_softness = 5.0f / 44.999f;
  const auto near_zero = Model::prepare({}, params, state);
  HS_EXPECT_NEAR(Model::run(seam, context, params, near_zero).coverage,
                 0.123008f, 2.0e-6f);

  params.line_thickness = 0.5f;
  params.line_softness = 20.0f / 44.75f;
  const auto thin = Model::prepare({}, params, state);
  const auto inside_edge = face_point(0.0625f);
  const auto within_fringe = face_point(1.25f);
  HS_EXPECT_EQ(Model::run(inside_edge, context, params, thin).coverage, 0.0f);
  const auto thin_fringe = Model::run(within_fringe, context, params, thin);
  HS_EXPECT_NEAR(thin_fringe.coverage, 0.5f, 5.0e-6f);
  const auto farther_in = face_point(3.25f);
  HS_EXPECT_EQ(Model::run(farther_in, context, params, thin).coverage, 1.0f);
  params.line_softness = 40.0f / 44.75f;
  const auto softer = Model::prepare({}, params, state);
  HS_EXPECT_EQ(Model::run(inside_edge, context, params, softer).coverage, 0.0f);
  const auto softer_fringe = Model::run(within_fringe, context, params, softer);
  HS_EXPECT_NEAR(softer_fringe.coverage, 0.25f, 5.0e-6f);
  HS_EXPECT_EQ(softer_fringe.value, thin_fringe.value);
  HS_EXPECT_NEAR(Model::run(farther_in, context, params, softer).coverage,
                 0.75f, 5.0e-6f);
  HS_EXPECT_EQ(Model::run(face_point(5.25f), context, params, softer).coverage,
               1.0f);
  params.line_thickness = 0.0f;

  const auto &geometry = *prepared.geometry;
  const uint16_t front = geometry.face_at(math::X_AXIS);
  HS_EXPECT_TRUE(front != Geometry::NO_FACE);
  HS_EXPECT_TRUE(geometry.face_at(math::Y_AXIS) != front);

  for (uint8_t solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    params.solid = solid;
    const auto tiled = Model::prepare({}, params, state);
    for (int sample = 0; sample < 192; ++sample) {
      const auto direction = math::fib_spiral(192, 0.5f, sample);
      const uint16_t owner = tiled.geometry->face_at(direction);
      HS_EXPECT_TRUE(owner != Geometry::NO_FACE);
      const auto field = Model::run({direction, 0.0f}, context, params, tiled);
      HS_EXPECT_TRUE(field.value >= 0.0f && field.value <= 1.0f);
      HS_EXPECT_TRUE(field.coverage > 0.0f && field.coverage <= 1.0f);
    }
  }

  params.solid = static_cast<uint8_t>(Solids::BaseMesh::CUBOCTAHEDRON);
  const auto square_center = math::X_AXIS;
  const auto triangle_center = math::Vector{1.0f, 1.0f, 1.0f}.normalized();
  const auto shared_seam = math::Vector{2.0f, 1.0f, 1.0f}.normalized();
  const float square_depth = std::atan(std::sqrt(0.5f));
  const float triangle_depth = std::acos(4.0f / std::sqrt(18.0f));
  const auto from_seam = [&](const math::Vector &center, float depth) {
    const auto tangent =
        (center - math::dot(center, shared_seam) * shared_seam).normalized();
    return PB::SphereSample{
        std::cos(depth) * shared_seam + std::sin(depth) * tangent, 0.0f};
  };
  for (float thickness : {0.0f, 2.0f}) {
    params.line_thickness = thickness;
    params.line_softness = 2.0f;
    const auto mixed = Model::prepare({}, params, state);
    const float half_width = 0.5f * thickness * DEG;
    const auto square = from_seam(
        square_center, half_width + 0.1f * (square_depth - half_width));
    const auto triangle = from_seam(
        triangle_center, half_width + 0.1f * (triangle_depth - half_width));
    HS_EXPECT_NEAR(Model::run(square, {}, params, mixed).coverage, 0.5f,
                   2.0e-5f);
    HS_EXPECT_NEAR(Model::run(triangle, {}, params, mixed).coverage, 0.5f,
                   2.0e-5f);
    HS_EXPECT_EQ(
        Model::run({triangle_center, 0.0f}, context, params, mixed).coverage,
        1.0f);

    params.line_softness = 10.0f;
    const auto full = Model::prepare({}, params, state);
    for (float fraction : {0.25f, 0.75f, 1.0f}) {
      const auto square_sample = from_seam(
          square_center, half_width + fraction * (square_depth - half_width));
      const auto triangle_sample =
          from_seam(triangle_center,
                    half_width + fraction * (triangle_depth - half_width));
      HS_EXPECT_NEAR(Model::run(square_sample, {}, params, full).coverage,
                     fraction, 2.0e-5f);
      HS_EXPECT_NEAR(Model::run(triangle_sample, {}, params, full).coverage,
                     fraction, 2.0e-5f);
    }
  }

  params.line_thickness = 0.0f;
  params.line_softness = 2.0f;
  context.angular_sample_radius = DEG;
  const auto mixed = Model::prepare({}, params, state);
  const auto square = from_seam(square_center, 0.0005f * DEG);
  const auto triangle = from_seam(triangle_center, 0.0005f * DEG);
  HS_EXPECT_NE(mixed.geometry->face_at(square.dir),
               mixed.geometry->face_at(triangle.dir));
  const float square_coverage =
      Model::run(square, context, params, mixed).coverage;
  const float triangle_coverage =
      Model::run(triangle, context, params, mixed).coverage;
  HS_EXPECT_TRUE(square_coverage > 0.0f && square_coverage < 1.0f);
  HS_EXPECT_NEAR(square_coverage, triangle_coverage, 1.0e-4f);
}

/** @brief A spherical tiling preview blends coverage and the palette position. */
inline void test_shader_chain_spherical_tiling_morph() {
  using WB = ShaderChainWhiteBox;
  using Model = In::Op::SampleSphericalTiling;
  using Seed = Solids::BaseMesh;
  reset_globals();
  auto effect = std::make_unique<WB::FX>();
  effect->init();
  const In::ChainEntryRequest chain[] = {
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID},
  };
  if (!expect_compiled(effect->set_chain(chain)))
    return;
  auto &program = WB::program(*effect);
  Model::Params from, to;
  from.solid = static_cast<uint8_t>(Seed::CUBE);
  to.solid = static_cast<uint8_t>(Seed::OCTAHEDRON);
  from.line_softness = 6.0f;
  to.line_softness = 8.0f;
  const ShaderChainParameterWrite from_writes[] = {
      {"source.solid", static_cast<float>(from.solid)},
      {"source.line-softness", from.line_softness},
  };
  const ShaderChainParameterWrite to_writes[] = {
      {"source.solid", static_cast<float>(to.solid)},
      {"source.line-softness", to.line_softness},
  };
  const std::span<const ShaderChainParameterWrite> a = from_writes;
  const std::span<const ShaderChainParameterWrite> b = to_writes;
  Model::State state;
  const auto ctx = WB::frame_context(*effect);
  const auto prepared_a = Model::prepare(ctx, from, state);
  const auto prepared_b = Model::prepare(ctx, to, state);
  HS_EXPECT_NE(prepared_a.geometry, prepared_b.geometry);
  HS_EXPECT_EQ(Model::prepare(ctx, from, state).geometry, prepared_a.geometry);
  const math::Vector views[] = {
      math::Vector(1.0f, 1.0f, .19f).normalized(),
      math::Vector(1.0f, 1.0f, .5f).normalized(),
      math::Vector(1.0f, .5f, 0.0f).normalized(),
      math::Vector(1.0f, .5f, .5f).normalized(),
  };
  for (bool reverse : {false, true}) {
    bool overlapping = false;
    bool exclusive = false;
    for (uint16_t frame : {0, 1, 4}) {
      const ShaderChainSourceMorphWrite morph{
          "source", reverse ? b : a, reverse ? a : b, {frame, 4, true, true}};
      const auto result = effect->update_parameters({}, true, &morph);
      HS_EXPECT_EQ(result, ParamSetResult::APPLIED);
      if (result != ParamSetResult::APPLIED)
        return;
      program.prepare(ctx);
      const float mix = math::ease_in_out_sin(static_cast<float>(frame) / 4.0f);
      for (const auto &view : views) {
        const PB::SphereSample point{view, 0.0f};
        const auto first = Model::run(point, ctx, from, prepared_a);
        const auto second = Model::run(point, ctx, to, prepared_b);
        const auto outgoing = reverse ? second : first;
        const auto incoming = reverse ? first : second;
        auto field = frame == 4 ? incoming : outgoing;
        if (frame == 1) {
          overlapping |= first.coverage > 0.0f && second.coverage > 0.0f &&
                         std::abs(first.value - second.value) > 0.1f;
          exclusive |= (first.coverage > 0.0f) != (second.coverage > 0.0f);
          const float old_alpha = (1.0f - mix) * outgoing.coverage;
          const float new_alpha = mix * incoming.coverage;
          field.coverage = old_alpha + new_alpha;
          field.value =
              field.coverage > 0.0f
                  ? (old_alpha * outgoing.value + new_alpha * incoming.value) /
                        field.coverage
                  : 0.0f;
        }
        const auto expected = run_entry<Color4>(program, ctx, 1, field);
        const auto actual = program.evaluate(view, ctx);
        HS_EXPECT_EQ(actual.color, expected.color);
        HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
      }
    }
    HS_EXPECT_TRUE(overlapping);
    HS_EXPECT_TRUE(exclusive);
  }
}

/** @brief A spherical source reaches cutout and palette without a projection. */
inline void test_shader_chain_spherical_source_pipeline() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest CHAIN[] = {
      {"source", In::Op::SampleSphericalCells::ID},
      {"cutout", In::Op::CoverageValueCutout::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  const auto result = program.compile(CHAIN);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  auto &params = param_as<In::Op::SampleSphericalCells::Params>(program, 0);
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  auto &cutout = param_as<In::Op::ValueCutoutChainParams>(program, 1);
  cutout.cutout_threshold = 0.75f;
  cutout.cutout_softness = 0.0f;
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  const PB::SphereSample CENTER{math::X_AXIS, 0.25f};
  const PB::SphereSample edge{math::Vector(1, 1, 0).normalized(), 0.25f};
  const auto center_field = run_entry<PB::FieldSample>(program, ctx, 0, CENTER);
  const auto edge_field = run_entry<PB::FieldSample>(program, ctx, 0, edge);
  HS_EXPECT_NEAR(center_field.value, 1.0f, 1.0e-6f);
  HS_EXPECT_TRUE(edge_field.value < 0.75f);
  const auto visible =
      run_entry<PB::FieldSample>(program, ctx, 1, center_field);
  const auto hidden = run_entry<PB::FieldSample>(program, ctx, 1, edge_field);
  HS_EXPECT_EQ(visible.coverage, 1.0f);
  HS_EXPECT_EQ(hidden.coverage, 0.0f);
  HS_EXPECT_EQ(hidden.path_length, 0.25f);
  HS_EXPECT_EQ(hidden.sphere.x, edge.dir.x);
  HS_EXPECT_EQ(hidden.sphere.y, edge.dir.y);
  const auto color = program.evaluate(CENTER.dir, ctx);
  const auto cut = program.evaluate(edge.dir, ctx);
  HS_EXPECT_EQ(color.alpha, 1.0f);
  HS_EXPECT_EQ(cut.alpha, 0.0f);
  HS_EXPECT_TRUE(carrier_live(color));
}

/** @brief Spherical noise keeps its instance clock through edits and restore. */
inline void test_shader_chain_spherical_noise_lifecycle() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  using Model = In::Op::SampleSphericalNoise;
  const In::ChainEntryRequest CHAIN[] = {
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  auto result = program.compile(CHAIN);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  param_as<Model::Params>(program, 0).noise_time_rate = 1.0f / 128.0f;
  for (int frame = 0; frame < 3; ++frame)
    program.advance();
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 3.0f / 128.0f);
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  SweepColors before;
  snapshot_render(program, ctx, before);
  const auto &runtime = program.ops()[0].op->runtime;
  const auto saved = runtime.capture_state(program.state_block(0));
  program.advance();
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 4.0f / 128.0f);
  HS_EXPECT_TRUE(runtime.restore_state(program.state_block(0), saved));
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 3.0f / 128.0f);
  program.prepare(ctx);
  SweepColors restored;
  snapshot_render(program, ctx, restored);
  expect_sweep_identical(restored, before);
  const auto seed = state_as<Model::State>(program, 0).noise_seed;
  const In::ChainEntryRequest EDITED[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(EDITED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).phase, 3.0f / 128.0f);
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).noise_seed, seed);
  const In::ChainEntryRequest RENAMED[] = {
      {"other", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(RENAMED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 0.0f);
  HS_EXPECT_TRUE(state_as<Model::State>(program, 0).noise_seed != seed);
}

} // namespace hs_test::shader_chain_tests
