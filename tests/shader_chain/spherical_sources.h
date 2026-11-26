/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_sources.h
 * @brief Spherical source carriers, composed policies and instance clocks. */

#include "core/render/pullback/composed_tiling.h"
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

/** @brief Build arenas large enough for every registry solid's geometry. */
struct TilingBuildArenas {
  static constexpr size_t BLOCK_BYTES = 2 * 1024 * 1024; ///< Geometry block.
  static constexpr size_t SCRATCH_BYTES = 512 * 1024;    ///< Each scratch.
  std::vector<uint8_t> block = std::vector<uint8_t>(BLOCK_BYTES); ///< Block.
  std::vector<uint8_t> storage =
      std::vector<uint8_t>(2 * SCRATCH_BYTES);    ///< Scratch storage.
  Arena scratch_a{storage.data(), SCRATCH_BYTES}; ///< Scratch A.
  Arena scratch_b{storage.data() + SCRATCH_BYTES,
                  SCRATCH_BYTES}; ///< Scratch B.

  /** @brief Builds a solid's edge-only geometry into the block.
   * @param solid Registry solid ordinal.
   * @param indices Edge indices to build.
   * @return Geometry viewing the block. */
  PB::Source::SphereTilingGeometry
  build(unsigned solid, PB::Source::SphereTilingGeometry::Indices indices =
                            PB::Source::SphereTilingGeometry::Indices::BOTH) {
    return PB::Source::SphereTilingGeometry::build(solid, block, scratch_a,
                                                   scratch_b, indices);
  }
};

/** @brief Thin and wide show bands, with parameters drawing each. */
struct TilingShowBand {
  float thickness; ///< Full edge width in degrees.
  float softness;  ///< Fringe width in degrees.
};

/// Filigree's thin band and a band wider than the narrow index serves.
inline constexpr TilingShowBand TILING_SHOW_BANDS[] = {{0.0f, 1.2350903f},
                                                       {6.0f, 4.0f}};

/** @brief The edge store admits exactly the solids whose larger build fits,
 *         and each admitted solid's thin build has its narrow index. */
inline void test_shader_chain_tiling_edge_store_admission() {
  using Store = PB::Source::SphereTilingEdgeStore;
  using Indices = PB::Source::SphereTilingGeometry::Indices;
  auto arenas = std::make_unique<TilingBuildArenas>();
  for (unsigned solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    HS_CONTEXT("solid", static_cast<int>(solid));
    const auto thin = arenas->build(solid, Indices::NARROW);
    const bool narrow = !thin.edge_view().narrow.cells.empty();
    const size_t thin_bytes = thin.bytes();
    const size_t bytes =
        std::max(thin_bytes, arenas->build(solid, Indices::COARSE).bytes());
    HS_EXPECT_EQ(bytes, size_t{PB::Source::SPHERE_TILING_EDGE_BYTES[solid]});
    HS_EXPECT_EQ(Store::admits(solid),
                 bytes + Store::ROUNDING_SLACK <= Store::BYTES);
    if (Store::admits(solid))
      HS_EXPECT_EQ(narrow, PB::Source::SphereTilingGeometry::narrow_indexed(
                               thin.edge_view().edges.size()));
    PB::Source::SphericalTilingSourceParams params;
    params.solid = static_cast<uint8_t>(solid);
    HS_EXPECT_EQ(PB::Source::valid(params), Store::admits(solid));
    params.output = static_cast<uint8_t>(PB::Source::TilingOutput::FACES);
    HS_EXPECT_FALSE(PB::Source::valid(params));
  }
  HS_EXPECT_FALSE(Store::admits(Solids::NUM_ENTRIES));
}

/** @brief Every admitted solid builds in the edge store with engine scratch,
 *         holding only the index its band reads, and draws edges at thin and
 *         wide show bands. */
inline void test_shader_chain_tiling_edge_store_builds_admitted() {
  using namespace PB::Source;
  reset_globals();
  auto store = std::make_unique<SphereTilingEdgeStore>();
  for (unsigned solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    if (!SphereTilingEdgeStore::admits(solid))
      continue;
    HS_CONTEXT("solid", static_cast<int>(solid));
    for (const auto &show : TILING_SHOW_BANDS) {
      HS_CONTEXT("thickness", static_cast<int>(show.thickness));
      SphericalTilingSourceParams params;
      params.solid = static_cast<uint8_t>(solid);
      params.line_thickness = show.thickness;
      params.line_softness = show.softness;
      store->refresh(params);
      store->refresh(params);
      HS_EXPECT_FALSE(store->refresh(params));
      HS_EXPECT_LE(store->bytes(), size_t{SPHERE_TILING_EDGE_BYTES[solid]});
      HS_EXPECT_EQ(scratch_arena_a.get_offset(), size_t{0});
      HS_EXPECT_EQ(scratch_arena_b.get_offset(), size_t{0});
      HS_EXPECT_EQ(store->built_solid(), params.solid);
      const auto view = store->edges();
      const SphereTilingEdgeBand band(show.thickness, show.softness);
      const bool narrow =
          band.outer <= SphereTilingNarrowIndex::OUTER &&
          SphereTilingGeometry::narrow_indexed(view.edges.size());
      HS_EXPECT_EQ(view.narrow.cells.empty(), !narrow);
      HS_EXPECT_EQ(view.offsets.empty(), narrow);
      int covered = 0;
      for (int i = 0; i < 4096; ++i)
        covered +=
            view.sample(math::fib_spiral(4096, 0.47f, i), band).coverage > 0.0f;
      HS_EXPECT_GT(covered, 0);
    }
  }
}

/** @brief Pruned coarse and narrow indices find the edge coverage an exhaustive
 *         scan of every edge finds. */
inline void test_shader_chain_tiling_indices_match_exhaustive_scan() {
  using namespace PB::Source;
  auto arenas = std::make_unique<TilingBuildArenas>();
  const SphereTilingEdgeBand thin{0.0f, 1.398054f};
  const SphereTilingEdgeBand wide{6.0f, 4.0f};
  for (unsigned solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    if (!SphereTilingEdgeStore::admits(solid))
      continue;
    HS_CONTEXT("solid", static_cast<int>(solid));
    const auto view = arenas->build(solid).edge_view();
    std::vector<uint16_t> all(view.edges.size());
    for (size_t id = 0; id < all.size(); ++id)
      all[id] = static_cast<uint16_t>(id);
    const std::vector<uint32_t> offsets = {
        0,
        static_cast<uint32_t>(all.size()),
        static_cast<uint32_t>(2 * all.size()),
        static_cast<uint32_t>(3 * all.size()),
        static_cast<uint32_t>(4 * all.size()),
        static_cast<uint32_t>(5 * all.size()),
        static_cast<uint32_t>(6 * all.size())};
    std::vector<uint16_t> ids;
    for (int face = 0; face < 6; ++face)
      ids.insert(ids.end(), all.begin(), all.end());
    const SphereTilingEdgeView exhaustive{view.vertices, view.edges, offsets,
                                          ids, 1};
    auto coarse = view;
    coarse.narrow = {};
    for (int i = 0; i < 4096; ++i) {
      const auto direction = math::fib_spiral(4096, 0.31f, i);
      HS_EXPECT_NEAR(coarse.sample(direction, wide).coverage,
                     exhaustive.sample(direction, wide).coverage, 1.0e-6f);
      HS_EXPECT_NEAR(view.sample(direction, thin).coverage,
                     exhaustive.sample(direction, thin).coverage, 1.0e-4f);
    }
  }
}

/** @brief Coarse and narrow tables of a view, rebuilt by scanning every edge
 *         for every cell. */
struct ExhaustiveTilingTables {
  std::vector<uint32_t> offsets{0}; ///< Coarse cell offsets.
  std::vector<uint16_t> ids;        ///< Coarse candidates.
  std::vector<uint16_t> cells;      ///< Narrow cells.
  std::vector<uint16_t> candidates; ///< Narrow multi-edge candidates.

  /** @brief Scans a view's edges.
   * @param view Built edge view supplying vertices, edges and resolutions. */
  explicit ExhaustiveTilingTables(
      const PB::Source::SphereTilingEdgeView &view) {
    using namespace PB::Source;
    using namespace PB::Source::tiling_edges_detail;
    const auto direction = Pattern::sphere_face_direction;
    const int steps = view.steps;
    const float step = 2.0f / steps;
    for (int cell = 0; cell < 6 * steps * steps; ++cell) {
      const int side = cell / (steps * steps);
      const float u = -1.0f + step * (cell % steps);
      const float v = -1.0f + step * ((cell / steps) % steps);
      const auto center = direction(side, u + 0.5f * step, v + 0.5f * step);
      float largest = 0.0f;
      for (float x : {u, u + step})
        for (float y : {v, v + step})
          largest = std::max(
              largest, math::distance_squared(center, direction(side, x, y)));
      const float radius = chord_angle(largest) + 1.0e-4f;
      float nearest = 4.0f;
      for (const auto &edge : view.edges)
        nearest =
            std::min(nearest, edge_chord_squared(edge, center, view.vertices));
      const float limit = chord_squared(
          std::min(math::PI_F, chord_angle(nearest) + 2.0f * radius));
      for (size_t id = 0; id < view.edges.size(); ++id)
        if (edge_chord_squared(view.edges[id], center, view.vertices) <= limit)
          ids.push_back(static_cast<uint16_t>(id));
      offsets.push_back(static_cast<uint32_t>(ids.size()));
    }
    if (view.narrow.cells.empty())
      return;
    constexpr int STEPS = SphereTilingNarrowIndex::STEPS;
    constexpr float STEP = 2.0f / STEPS;
    constexpr float PAD = 1.0e-4f;
    const float limit_sine = std::sin(SphereTilingNarrowIndex::OUTER + PAD);
    for (int face = 0; face < 6; ++face)
      for (int y = 0; y < STEPS; ++y)
        for (int x = 0; x < STEPS; ++x) {
          const float u = -1.0f + x * STEP;
          const float v = -1.0f + y * STEP;
          const auto center = direction(face, u + 0.5f * STEP, v + 0.5f * STEP);
          const math::Vector corners[] = {direction(face, u, v),
                                          direction(face, u + STEP, v),
                                          direction(face, u, v + STEP),
                                          direction(face, u + STEP, v + STEP)};
          float radius = 0.0f;
          for (const auto &corner : corners)
            radius = std::max(radius, math::distance_squared(center, corner));
          const float limit = chord_squared(
              chord_angle(radius) + SphereTilingNarrowIndex::OUTER + PAD);
          const size_t begin = candidates.size();
          for (size_t id = 0; id < view.edges.size(); ++id) {
            const auto &edge = view.edges[id];
            if (edge_chord_squared(edge, center, view.vertices) > limit)
              continue;
            bool front = true, back = true, interior = edge.inv_length > 0.0f;
            for (const auto &corner : corners) {
              const float sine = math::dot(edge.normal, corner);
              front &= sine > limit_sine;
              back &= sine < -limit_sine;
              interior &= math::dot(edge.from_a, corner) > PAD &&
                          math::dot(edge.from_b, corner) > PAD;
            }
            if (!front && !back)
              candidates.push_back(static_cast<uint16_t>(
                  id | (interior ? SphereTilingNarrowIndex::INTERIOR : 0)));
          }
          const size_t count = candidates.size() - begin;
          if (count == 0) {
            cells.push_back(0);
          } else if (count == 1) {
            cells.push_back(candidates.back() | SphereTilingNarrowIndex::LAST);
            candidates.pop_back();
          } else {
            cells.push_back(static_cast<uint16_t>(begin + 1));
            candidates.back() |= SphereTilingNarrowIndex::LAST;
          }
        }
  }
};

/** @brief Every admitted solid's edge-driven tables equal the exhaustive
 *         scan's. */
inline void test_shader_chain_tiling_tables_match_exhaustive_scan() {
  using namespace PB::Source;
  auto arenas = std::make_unique<TilingBuildArenas>();
  for (unsigned solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    if (!SphereTilingEdgeStore::admits(solid))
      continue;
    HS_CONTEXT("solid", static_cast<int>(solid));
    const auto view = arenas->build(solid).edge_view();
    const ExhaustiveTilingTables expected(view);
    HS_EXPECT_TRUE(std::ranges::equal(view.offsets, expected.offsets));
    HS_EXPECT_TRUE(std::ranges::equal(view.candidate_ids, expected.ids));
    HS_EXPECT_TRUE(std::ranges::equal(view.narrow.cells, expected.cells));
    HS_EXPECT_TRUE(
        std::ranges::equal(view.narrow.candidates, expected.candidates));
    using Indices = SphereTilingGeometry::Indices;
    const auto wide = arenas->build(solid, Indices::COARSE).edge_view();
    HS_EXPECT_TRUE(std::ranges::equal(wide.offsets, expected.offsets));
    HS_EXPECT_TRUE(std::ranges::equal(wide.candidate_ids, expected.ids));
    HS_EXPECT_TRUE(wide.narrow.cells.empty());
    if (!expected.cells.empty()) {
      const auto thin = arenas->build(solid, Indices::NARROW).edge_view();
      HS_EXPECT_TRUE(thin.offsets.empty());
      HS_EXPECT_TRUE(std::ranges::equal(thin.narrow.cells, expected.cells));
      HS_EXPECT_TRUE(
          std::ranges::equal(thin.narrow.candidates, expected.candidates));
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

/** @brief Restoring the retired tiling output option 2 selects Faces, and no
 *         other out-of-range output value is admitted. */
inline void test_shader_chain_spherical_tiling_legacy_output_restore() {
  using WB = ShaderChainWhiteBox;
  using Model = In::Op::SampleSphericalTiling;
  using Output = PB::Source::TilingOutput;
  constexpr float RETIRED_OUTPUT = 2.0f;
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
  const auto output = [&] {
    return param_as<Model::Params>(program, 0).output;
  };
  HS_EXPECT_EQ(output(), static_cast<uint8_t>(Output::EDGES));
  const ChainSnapshot saved = effect->snapshot();
  const auto with_output = [&](float value) {
    ChainSnapshot snapshot = saved;
    bool found = false;
    for (auto &parameter : snapshot.parameters)
      if (parameter.name == "source.output") {
        parameter.value = value;
        found = true;
      }
    HS_EXPECT_TRUE(found);
    return snapshot;
  };
  HS_EXPECT_EQ(effect->restore_snapshot(with_output(RETIRED_OUTPUT)),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(output(), static_cast<uint8_t>(Output::FACES));
  HS_EXPECT_EQ(effect->restore_snapshot(with_output(RETIRED_OUTPUT + 1.0f)),
               ChainSnapshotRestoreResult::INVALID_VALUE);
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
