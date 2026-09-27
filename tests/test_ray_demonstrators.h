/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/ray/events.h"
#include "core/render/ray/shade.h"
#include "core/render/ray/query.h"
#include "core/render/pullback/ray.h"
#include "core/render/sdf/framework.h"
#include "core/render/sdf/lattice_field.h"
#include "core/render/sdf/volume.h"
#include "tests/test_harness.h"

namespace hs_test::ray_demonstrator_tests {

inline void test_framework_generic_event_rendering() {
  SDF::TriangularFramework geometry;
  Raycast::TraceLimits limits;
  limits.max_candidates = 64;
  limits.max_layers = 32;
  const Raycast::Ray RAY{{0.0f, 0.0f, 0.3f}, {0.0f, 0.0f, -1.0f}, {0.0f, 2.3f}};
  SDF::FrameworkEvents events(geometry, RAY);
  std::array<Raycast::Contribution, 8> output;
  size_t count = 0;
  const auto RESULT = Raycast::trace_events(
      events, RAY.interval, limits, [&](const Raycast::Contribution &hit) {
        output[count++] = hit;
        return true;
      });
  HS_EXPECT_EQ(RESULT.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(count, size_t{3});
  for (size_t i = 0; i < count; ++i) {
    HS_EXPECT_NEAR(output[i].t, 0.3f + static_cast<float>(i), 1e-6f);
    HS_EXPECT_EQ(output[i].coverage, 1.0f);
    HS_EXPECT_FALSE(output[i].verified);
  }
  const math::Vector DIRECTION = math::Vector{0.7f, 0.4f, 0.3f}.normalized();
  const Raycast::Ray OBLIQUE{{0.0f, 0.0f, 0.0f}, DIRECTION, {0.0f, 5.0f}};
  SDF::FrameworkEvents slanted(geometry, OBLIQUE, {0.05f, 0.2f});
  for (size_t i = 0; i < SDF::FrameworkEvents::STREAM_COUNT; ++i)
    HS_EXPECT_TRUE(slanted.active(i));
  float previous = -1.0f;
  count = 0;
  const auto SLANTED = Raycast::trace_events(
      slanted, OBLIQUE.interval, limits, [&](const Raycast::Contribution &hit) {
        HS_EXPECT_GE(hit.t, previous);
        previous = hit.t;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(SLANTED.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_GT(count, size_t{1});
}

struct CloseStreams {
  static constexpr size_t STREAM_COUNT = 3;
  std::array<std::array<float, 2>, STREAM_COUNT> times = {
      {{{0.2f, 0.20003f}}, {{0.20001f, 0.9f}}, {{0.20002f, 1.0f}}}};
  std::array<size_t, STREAM_COUNT> indices{};
  bool active(size_t i) const { return indices[i] < times[i].size(); }
  float distance(size_t i) const { return times[i][indices[i]]; }
  Raycast::Contribution candidate(size_t i) const {
    Raycast::Contribution result;
    result.t = distance(i);
    result.merge_identity = i == 1 ? 0 : i;
    result.material = i == 2 ? 1 : 0;
    result.coverage = i == 1 ? 0.8f : 0.3f;
    return result;
  }
  void advance(size_t i) { ++indices[i]; }
};

inline void test_repeated_stream_grouping_preserves_order_and_endpoints() {
  CloseStreams events;
  Raycast::TraceLimits limits;
  std::array<Raycast::Contribution, 8> output;
  size_t count = 0;
  const auto RESULT = Raycast::trace_events(
      events, {0.0f, 1.0f}, limits, [&](const Raycast::Contribution &hit) {
        output[count++] = hit;
        return true;
      });
  HS_EXPECT_EQ(RESULT.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(RESULT.counters.candidates, 6);
  HS_EXPECT_EQ(count, size_t{4});
  HS_EXPECT_NEAR(output[0].t, 0.2f, 1e-7f);
  HS_EXPECT_NEAR(output[0].coverage, 0.8f, 1e-7f);
  HS_EXPECT_EQ(output[1].material, uint32_t{1});
  HS_EXPECT_EQ(output[count - 1].t, 1.0f);
  for (size_t i = 1; i < count; ++i)
    HS_EXPECT_GE(output[i].t, output[i - 1].t);
}

inline void test_lattice_volume_camera_demonstrators() {
  Raycast::PreparedCamera camera;
  camera.center = {{0.3f, 0.2f, 0.0f, 0.0f}};
  camera.interval = {0.0f, 1.0f};
  SDF::WireLattice<3> cubic;
  Raycast::DomainQuery3 query3{cubic, camera};
  Raycast::TraceLimits limits;
  limits.max_queries = 256;
  limits.max_steps = 256;
  limits.max_refinements = 256;
  limits.position_tolerance = 1e-5f;
  const math::Vector DIRECTION{0.0f, -1.0f, 0.0f};
  const auto CUBIC =
      Raycast::surface_search(query3, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_TRUE(CUBIC.has_surface);
  HS_EXPECT_TRUE(CUBIC.contribution.verified);
  HS_EXPECT_NEAR(CUBIC.contribution.t, 0.15f, 1e-5f);

  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.center[3] = 0.03f;
  SDF::WireLattice<4> hypercubic;
  Raycast::DomainQuery4 query4{hypercubic, camera};
  const auto HYPERCUBIC =
      Raycast::surface_search(query4, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_TRUE(HYPERCUBIC.has_surface);
  HS_EXPECT_NEAR(HYPERCUBIC.contribution.t, 0.16f, 1e-5f);
  camera.center[1] = 0.0f;
  const auto EXIT = Raycast::surface_search(
      query4, camera.ray(math::Vector{0.0f, 1.0f, 0.0f}), {}, limits);
  HS_EXPECT_TRUE(EXIT.has_surface);
  HS_EXPECT_NEAR(EXIT.contribution.t, 0.04f, 1e-5f);
  camera.center[3] = 0.1f;
  const auto EMPTY_SLICE =
      Raycast::surface_search(query4, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_FALSE(EMPTY_SLICE.has_surface);
  HS_EXPECT_EQ(EMPTY_SLICE.status, Raycast::TraceStatus::RANGE_COMPLETE);
}

struct WhitePalette {
  Color4 get(float) const { return {{50000, 30000, 10000}, 1.0f}; }
};

template <typename Shape> struct VolumeFrame {
  const Shape *shape;
  Raycast::PreparedCamera camera;
  math::Vector object_center;
  math::Quaternion inverse_rotation;
  Raycast::Appearance appearance;
  mutable int preparations = 0;
};

template <typename Shape> struct VolumeRenderer {
  struct Prepared {
    Raycast::PreparedCamera camera;
    math::Vector center;
    math::Quaternion inverse_rotation;
  };
  static Prepared prepare(const VolumeFrame<Shape> &frame) {
    ++frame.preparations;
    return {frame.camera, frame.object_center, frame.inverse_rotation};
  }
  static Raycast::ShadedTrace trace(const math::Vector &direction,
                                    const VolumeFrame<Shape> &frame,
                                    const Prepared &prepared) {
    Raycast::VolumeQuery query{
        *frame.shape, Raycast::QueryCapabilities{true, false, true, 0}};
    Raycast::PlacedQuery placed{query, prepared.center,
                                prepared.inverse_rotation, 1.0f};
    Raycast::DomainQuery3 domain{placed, prepared.camera};
    Raycast::TraceLimits limits;
    limits.max_queries = 256;
    limits.max_steps = 256;
    limits.max_refinements = 256;
    return Raycast::shade_surface(domain, prepared.camera.ray(direction), {},
                                  limits, frame.appearance);
  }
  static Color4 shade(const math::Vector &direction,
                      const VolumeFrame<Shape> &frame,
                      const Prepared &prepared) {
    return trace(direction, frame, prepared).color;
  }
};

template <typename Shape> struct VolumeBinding {
  using FrameState = VolumeFrame<Shape>;
  using Instrumentation = Pullback::NoInstrumentation;
};

template <typename Shape>
void check_placed_volume_stage(const Shape &shape,
                               const BakedPalette &palette) {
  using Renderer = VolumeRenderer<Shape>;
  using Pipeline =
      Pullback::Pipeline<VolumeBinding<Shape>, Pullback::RayStage<Renderer>>;
  VolumeFrame<Shape> frame{};
  frame.shape = &shape;
  frame.object_center = {100.0f, -50.0f, 30.0f};
  frame.camera.center = {{100.0f, -50.0f, 27.0f, 0.0f}};
  frame.camera.interval = {0.0f, 7.0f};
  frame.appearance.inv_far = 1.0f / 7.0f;
  frame.appearance.near_start = -1.0f;
  frame.appearance.depth_palette = &palette;
  frame.appearance.feature_palette = &palette;
  std::array<uint64_t, 2> hashes{};
  for (int pose = 0; pose < 2; ++pose) {
    frame.inverse_rotation =
        math::make_rotation(math::Vector{1.0f, 0.0f, 0.0f}, -pose * 0.8f);
    const auto PREPARED = Pipeline::prepare_stages(frame);
    HS_EXPECT_EQ(frame.preparations, pose + 1);
    uint64_t hash = hs_test::FNV1A64_BASIS;
    int hits = 0;
    int misses = 0;
    for (int y = -6; y <= 6; ++y) {
      for (int x = -6; x <= 6; ++x) {
        const math::Vector DIRECTION =
            math::Vector{x * 0.08f, y * 0.08f, 1.0f}.normalized();
        const Color4 COLOR = Pipeline::evaluate(DIRECTION, frame, PREPARED);
        const auto TRACE =
            Renderer::trace(DIRECTION, frame, std::get<0>(PREPARED));
        HS_EXPECT_EQ(COLOR.color.r, TRACE.color.color.r);
        HS_EXPECT_EQ(COLOR.alpha, TRACE.color.alpha);
        HS_EXPECT_TRUE(Raycast::finite(COLOR.alpha));
        HS_EXPECT_GE(COLOR.alpha, 0.0f);
        HS_EXPECT_LE(COLOR.alpha, 1.0f);
        if (TRACE.trace.has_surface) {
          ++hits;
          HS_EXPECT_TRUE(TRACE.trace.contribution.verified);
          HS_EXPECT_GT(COLOR.alpha, 0.0f);
        } else {
          ++misses;
          HS_EXPECT_EQ(COLOR.alpha, 0.0f);
        }
        hash = hs_test::fnv1a64_channel(hash, COLOR.color.r);
      }
    }
    HS_EXPECT_GT(hits, 10);
    HS_EXPECT_GT(misses, 10);
    HS_EXPECT_EQ(frame.camera.center[0], 100.0f);
    HS_EXPECT_EQ(frame.camera.center[2], 27.0f);
    HS_EXPECT_EQ(frame.preparations, pose + 1);
    hashes[pose] = hash;
  }
  HS_EXPECT_NE(hashes[0], hashes[1]);
}

inline void test_torus_and_warped_volume_spherical_stage() {
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, WhitePalette{});
  const SDF::Torus TORUS{1.0f, 0.25f};
  check_placed_volume_stage(TORUS, palette.view());
  const SDF::WarpedVolume WARPED{TORUS, SDF::Warp::Twist{3, 0.15f, 1.0f},
                                 0.001f};
  check_placed_volume_stage(WARPED, palette.view());
}

struct OffSliceBall {
  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }
  Raycast::QuerySample sample(const math::Vec4 &p) const {
    const float W = p[3] - 1.001f;
    const float VALUE =
        sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + W * W) - 1.0f;
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, 0};
  }
};

inline void test_verified_filter_off_slice_geometry_and_projected_normal() {
  Raycast::PreparedCamera camera;
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.interval = {0.0f, 0.1f};
  OffSliceBall ball;
  Raycast::DomainQuery4 query{ball, camera};
  constexpr std::array<math::Vector, 4> DIRECTIONS = {
      math::Vector{1.0f, 0.0f, 0.0f}, math::Vector{-1.0f, 0.0f, 0.0f},
      math::Vector{0.0f, 1.0f, 0.0f}, math::Vector{0.0f, -1.0f, 0.0f}};
  const auto FILTERED =
      Raycast::verified_filter(DIRECTIONS, [&](const auto &dir) {
        Raycast::ShadedTrace sample;
        sample.trace = Raycast::surface_search(query, camera.ray(dir), {}, {});
        sample.color = {{65535, 65535, 65535}, 1.0f};
        return sample;
      });
  HS_EXPECT_FALSE(FILTERED.trace.has_surface);
  HS_EXPECT_EQ(FILTERED.color.alpha, 0.0f);
  math::Vector projected{1.0f, 1.0f, 1.0f};
  HS_EXPECT_FALSE(camera.project_normal({{0.0f, 0.0f, 0.0f, 1.0f}}, projected));
  HS_EXPECT_EQ(projected.magnitude(), 0.0f);
  const auto PARTIAL =
      Raycast::verified_filter(DIRECTIONS, [](const auto &dir) {
        Raycast::ShadedTrace sample;
        if (dir.x > 0.0f) {
          sample.trace.status = Raycast::TraceStatus::SURFACE;
          sample.trace.has_surface = true;
          sample.trace.contribution.verified = true;
          sample.color = {{40000, 20000, 10000}, 0.8f};
        } else if (dir.y < 0.0f) {
          sample.trace.status = Raycast::TraceStatus::UNRESOLVED;
        }
        return sample;
      });
  HS_EXPECT_EQ(PARTIAL.trace.status, Raycast::TraceStatus::UNRESOLVED);
  HS_EXPECT_TRUE(PARTIAL.trace.has_surface);
  HS_EXPECT_TRUE(PARTIAL.trace.contribution.verified);
  HS_EXPECT_NEAR(PARTIAL.trace.contribution.coverage, 0.25f, 1e-7f);
  HS_EXPECT_NEAR(PARTIAL.color.alpha, 0.2f, 1e-7f);
  HS_EXPECT_EQ(PARTIAL.color.color.r, uint16_t{40000});
}

inline int run_ray_demonstrator_tests() {
  const auto MODULE = hs_test::begin_module("ray_demonstrators");
  test_framework_generic_event_rendering();
  test_repeated_stream_grouping_preserves_order_and_endpoints();
  test_lattice_volume_camera_demonstrators();
  test_torus_and_warped_volume_spherical_stage();
  test_verified_filter_off_slice_geometry_and_projected_normal();
  return hs_test::end_module(MODULE);
}

} // namespace hs_test::ray_demonstrator_tests
