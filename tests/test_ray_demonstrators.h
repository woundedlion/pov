/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Analytic and periodic ray demonstrator rendering contracts.
 */
#pragma once

#include "core/render/ray/events.h"
#include "core/render/ray/shade.h"
#include "core/render/ray/query.h"
#include "core/render/pullback/ray.h"
#include "core/render/sdf/framework.h"
#include "core/render/sdf/lattice_field.h"
#include "core/render/sdf/volume.h"
#include "core/render/sdf/affine_lattice.h"
#include "core/render/sdf/periodic_shells.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"

namespace hs_test::ray_demonstrator_tests {

/** @brief Pins octet crossing coverage against ray line distance. */
inline void test_octet_crossing_coverage_against_ray_line_distance() {
  int compared = 0;
  int hits = 0;
  for (int sample = 0; sample < 600; ++sample) {
    SDF::OctetFramework geometry;
    geometry.cell_size = sample % 3 == 0   ? .25f
                         : sample % 3 == 1 ? 1.5f
                                           : 10.0f;
    geometry.wire_radius = geometry.cell_size * .055f;
    geometry.origin = {.3f, -.2f, .7f};
    const Raycast::Ray RAY{{sinf(sample * .43f) * 5, cosf(sample * .71f) * 5,
                            sinf(sample * .19f) * 5},
                           math::Vector{cosf(sample * .31f),
                                        sinf(sample * .57f),
                                        cosf(sample * .23f)}
                               .normalized(),
                           {.1f, 16.0f}};
    const Raycast::Footprint FOOTPRINT{sample % 5 == 0 ? 0.0f : .02f, .3f};
    const auto PLANES = geometry.plane_families();
    SDF::OctetEvents events(geometry, RAY, FOOTPRINT);
    SDF::FrameworkPlaneEvents reference(RAY, FOOTPRINT);
    reference.initialize(PLANES, geometry.origin, geometry.valid());
    std::array<bool, SDF::OctetEvents::STREAM_COUNT> owns_pair{};
    for (size_t i = 0; i < PLANES.size(); ++i) {
      for (size_t j = i + 1; j < PLANES.size(); ++j) {
        const float A = math::dot(RAY.direction, PLANES[i].normal);
        const float B = math::dot(RAY.direction, PLANES[j].normal);
        const size_t OWNER = fabsf(A) >= fabsf(B) ? i : j;
        owns_pair[OWNER] = owns_pair[OWNER] || (OWNER == i ? A : B) != 0;
      }
    }
    for (size_t stream = 0; stream < PLANES.size(); ++stream)
      reference.cursors[stream].active &= owns_pair[stream];
    for (size_t stream = 0; stream < PLANES.size(); ++stream) {
      HS_EXPECT_EQ(events.active(stream), reference.active(stream));
      for (int crossing = 0; crossing < 8 && events.active(stream);
           ++crossing) {
        ++compared;
        const auto HIT = events.candidate(stream);
        HS_EXPECT_NEAR(HIT.t, reference.distance(stream),
                       2e-5f * std::max(1.0f, HIT.t));
        const auto POINT = RAY.at(HIT.t) - geometry.origin;
        float best = INFINITY;
        uint32_t feature = 0;
        uint32_t pair = 0;
        for (size_t i = 0; i < PLANES.size(); ++i)
          for (size_t j = i + 1; j < PLANES.size(); ++j, ++pair) {
            const float A = math::dot(RAY.direction, PLANES[i].normal);
            const float B = math::dot(RAY.direction, PLANES[j].normal);
            const size_t OWNER = fabsf(A) >= fabsf(B) ? i : j;
            if (OWNER != stream)
              continue;
            const size_t OTHER = OWNER == i ? j : i;
            const float U = math::dot(POINT, PLANES[OTHER].normal);
            const float RESIDUAL =
                U - PLANES[OTHER].spacing * roundf(U / PLANES[OTHER].spacing);
            const auto CROSS = math::cross(
                RAY.direction, math::cross(PLANES[i].normal, PLANES[j].normal));
            const float DISTANCE =
                fabsf(RESIDUAL * (OWNER == i ? A : B)) / CROSS.magnitude();
            if (DISTANCE < best) {
              best = DISTANCE;
              feature = pair;
            }
          }
        const float WIDTH = FOOTPRINT.at(HIT.t);
        const float COVERAGE =
            WIDTH > 0 ? std::clamp(.5f - (best - geometry.wire_radius) / WIDTH,
                                   0.0f, 1.0f)
            : best <= geometry.wire_radius ? 1.0f
                                           : 0.0f;
        HS_EXPECT_NEAR(HIT.coverage, COVERAGE, 3e-4f);
        if (COVERAGE > 0) {
          ++hits;
          HS_EXPECT_EQ(HIT.feature, feature);
        }
        events.advance(stream);
        reference.advance(stream);
        HS_EXPECT_EQ(events.active(stream), reference.active(stream));
      }
    }
  }
  HS_EXPECT_TRUE(compared > 1000);
  HS_EXPECT_TRUE(hits > compared / 20);
}

/** @brief Pins framework generic event rendering. */
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
        if (count < output.size())
          output[count] = hit;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(RESULT.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(count, size_t{3});
  for (size_t i = 0; i < count && i < output.size(); ++i) {
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

/** @brief Pins repeated stream grouping preserves order and endpoints. */
inline void test_repeated_stream_grouping_preserves_order_and_endpoints() {
  CloseStreams events;
  Raycast::TraceLimits limits;
  std::array<Raycast::Contribution, 8> output;
  size_t count = 0;
  const auto RESULT = Raycast::trace_events(
      events, {0.0f, 1.0f}, limits, [&](const Raycast::Contribution &hit) {
        if (count < output.size())
          output[count] = hit;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(RESULT.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(RESULT.counters.candidates, 6);
  HS_EXPECT_EQ(count, size_t{4});
  if (count == 0 || count > output.size())
    return;
  HS_EXPECT_NEAR(output[0].t, 0.2f, 1e-7f);
  HS_EXPECT_NEAR(output[0].coverage, 0.8f, 1e-7f);
  HS_EXPECT_EQ(output[1].material, uint32_t{1});
  HS_EXPECT_EQ(output[count - 1].t, 1.0f);
  for (size_t i = 1; i < count; ++i)
    HS_EXPECT_GE(output[i].t, output[i - 1].t);
}

/** @brief Pins lattice volume camera demonstrators. */
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

struct ConstantPalette {
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
  frame.appearance.palette = &palette;
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

/** @brief Pins torus and warped volume spherical stage. */
inline void test_torus_and_warped_volume_spherical_stage() {
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, ConstantPalette{});
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

/** @brief Pins verified filter off slice geometry and projected normal. */
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

/** @brief Pins affine lattice world metric and periods. */
inline void test_affine_lattice_world_metric_and_periods() {
  const SDF::AffineLattice GEOMETRY{2, .55f, 1.4f};
  const math::Vec4 POINT{{.13f, -.24f, .42f, .37f}};
  const auto WORLD = GEOMETRY.point(POINT);
  const auto RECOVERED = GEOMETRY.inverse(WORLD);
  for (int k = 0; k < 4; ++k) {
    HS_EXPECT_NEAR(RECOVERED[k], POINT[k], 1e-6f);
    auto translated = POINT;
    translated[k] += 3;
    const auto WRAPPED = GEOMETRY.wrap(GEOMETRY.point(translated));
    const auto EXPECTED = GEOMETRY.wrap(WORLD);
    for (int j = 0; j < 4; ++j)
      HS_EXPECT_NEAR(WRAPPED[j], EXPECTED[j], 1e-6f);
  }
  Raycast::PreparedCamera camera;
  camera.center = {{.1f, .2f, .08f, 0}};
  camera.interval = {0, 8};
  SDF::AffineLatticeEvents events(camera, {1, 0, 0}, GEOMETRY, .05f, {.1f, 0});
  const auto HIT = events.candidate(0);
  const float EXPECTED =
      std::clamp(.5f - (.08f - .05f) / (.1f * HIT.t), 0.0f, 1.0f);
  HS_EXPECT_NEAR(HIT.coverage, EXPECTED, 1e-5f);
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.center[3] = .06f;
  SDF::AffineLatticeEvents slice(camera, {1, 0, 0}, GEOMETRY, .05f, {.1f, 0});
  const auto SLICE_HIT = slice.candidate(0);
  HS_EXPECT_NEAR(
      SLICE_HIT.coverage,
      std::clamp(.5f - (.1f - .05f) / (.1f * SLICE_HIT.t), 0.0f, 1.0f), 1e-5f);
}

/** @brief Pins periodic shell roots and slices. */
inline void test_periodic_shell_roots_and_slices() {
  SDF::PeriodicShells geometry{2, .3f};
  const auto SPHERE = geometry.intersect({{-2, 0, 0, 0}}, {{1, 0, 0, 0}}, 3);
  HS_EXPECT_TRUE(SPHERE.hit);
  HS_EXPECT_NEAR(SPHERE.near, 1.4f, 1e-6f);
  HS_EXPECT_NEAR(SPHERE.far, 2.6f, 1e-6f);
  const auto SLICE = geometry.intersect({{-2, 0, 0, .36f}}, {{1, 0, 0, 0}}, 4);
  HS_EXPECT_NEAR(SLICE.near, 1.52f, 1e-6f);
  HS_EXPECT_NEAR(SLICE.far, 2.48f, 1e-6f);
  HS_EXPECT_FALSE(geometry.intersect({{-2, 0, 0, .7f}}, {{1, 0, 0, 0}}, 4).hit);
  const auto INSIDE = geometry.intersect({{0, 0, 0, 0}}, {{1, 0, 0, 0}}, 4);
  HS_EXPECT_NEAR(INSIDE.near, -.6f, 1e-6f);
  HS_EXPECT_NEAR(INSIDE.far, .6f, 1e-6f);
}

/** @brief Pins affine cached metric against ray line oracle. */
inline void test_affine_cached_metric_against_ray_line_oracle() {
  int compared = 0;
  int hits = 0;
  for (int sample = 0; sample < 240; ++sample) {
    Raycast::PreparedCamera camera;
    camera.domain = sample % 2 ? Raycast::SamplingDomain::SLICE_4D
                               : Raycast::SamplingDomain::SPATIAL_3D;
    camera.center = {{sinf(sample * .13f), cosf(sample * .31f),
                      sinf(sample * .17f), sample % 2 ? .27f : 0}};
    math::rotate_plane(camera.embedding, 0, 1, sample * .19f);
    if (sample % 2)
      math::rotate_plane(camera.embedding, 1, 3, sample * .21f);
    const SDF::AffineLattice GEOMETRY{sample % 3 == 0 ? .25f : 2.0f,
                                      sinf(sample * .37f),
                                      .5f + (sample % 7) * .25f};
    const math::Vector VIEW =
        math::Vector{cosf(sample * .51f), sinf(sample * .41f), .43f}
            .normalized();
    const Raycast::Footprint FOOTPRINT{.12f, .1f};
    SDF::AffineLatticeEvents events(camera, VIEW, GEOMETRY,
                                    .055f * GEOMETRY.cell_size, FOOTPRINT);
    const auto DIRECTION =
        GEOMETRY.inverse(camera.embedding.apply({{VIEW.x, VIEW.y, VIEW.z, 0}}));
    for (int plane = 0; plane < events.dimensions; ++plane) {
      bool expected_active = false;
      for (int free = 0; free < events.dimensions; ++free) {
        int owner = -1;
        float speed = 0;
        for (int other = 0; other < events.dimensions; ++other) {
          if (other != free && fabsf(DIRECTION[other]) > speed) {
            speed = fabsf(DIRECTION[other]);
            owner = other;
          }
        }
        math::Vec4 axis{};
        axis[free] = 1;
        const auto U = GEOMETRY.point(axis);
        double uu = 0, ud = 0;
        for (int k = 0; k < events.dimensions; ++k) {
          uu += static_cast<double>(U[k]) * U[k];
          ud += static_cast<double>(U[k]) * events.ambient_direction[k];
        }
        expected_active |= owner == plane && uu - ud * ud > 1e-12;
      }
      expected_active &= DIRECTION[plane] != 0;
      HS_EXPECT_EQ(events.active(plane), expected_active);
      for (int crossing = 0; crossing < 4 && events.active(plane); ++crossing) {
        ++compared;
        const auto HIT = events.candidate(plane);
        double best = INFINITY;
        for (int free = 0; free < events.dimensions; ++free) {
          if (events.owner[free] != plane)
            continue;
          math::Vec4 axis{};
          axis[free] = 1;
          const auto U = GEOMETRY.point(axis);
          double uu = 0, ud = 0, dd = 0;
          for (int k = 0; k < events.dimensions; ++k) {
            uu += static_cast<double>(U[k]) * U[k];
            ud += static_cast<double>(U[k]) * events.ambient_direction[k];
            dd += static_cast<double>(events.ambient_direction[k]) *
                  events.ambient_direction[k];
          }
          if (uu - ud * ud <= 1e-12)
            continue;
          const int COUNT = events.dimensions == 4 ? 9 : 3;
          for (int neighbor = 0; neighbor < COUNT; ++neighbor) {
            int digits = neighbor;
            math::Vec4 residual{};
            for (int k = 0; k < events.dimensions; ++k) {
              if (k == free || k == plane)
                continue;
              const float P = events.origin[k] + HIT.t * events.direction[k];
              residual[k] = P - roundf(P) + static_cast<float>(digits % 3 - 1);
              digits /= 3;
            }
            const auto R = GEOMETRY.point(residual);
            double rr = 0, ru = 0, rd = 0;
            for (int k = 0; k < events.dimensions; ++k) {
              rr += static_cast<double>(R[k]) * R[k];
              ru += static_cast<double>(R[k]) * U[k];
              rd += static_cast<double>(R[k]) * events.ambient_direction[k];
            }
            best = std::min(
                best, std::max(0.0, rr - (dd * ru * ru - 2 * ud * ru * rd +
                                          uu * rd * rd) /
                                             (uu * dd - ud * ud)));
          }
        }
        const float EXPECTED =
            std::clamp(.5f - (static_cast<float>(sqrt(best)) - events.radius) /
                                 FOOTPRINT.at(HIT.t),
                       0.0f, 1.0f);
        HS_EXPECT_NEAR(HIT.coverage, EXPECTED, 3e-5f);
        hits += EXPECTED > 0;
        events.advance(plane);
        HS_EXPECT_EQ(events.active(plane), expected_active);
      }
    }
  }
  HS_EXPECT_TRUE(compared > 1000);
  HS_EXPECT_TRUE(hits > 10);
}

/** @brief Pins periodic shell traversal budgets. */
inline void test_periodic_shell_traversal_budgets() {
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, ConstantPalette{});
  const Raycast::Appearance APPEARANCE{.1f, 0, 1, &palette.view()};
  Raycast::PreparedCamera camera;
  camera.center = {{-.5f, 0, 0, 0}};
  camera.interval = {0, 1};
  Raycast::TraceLimits limits;
  const auto COMPLETE = SDF::shade_periodic_shells(camera, {1, 0, 0}, 1, .3f,
                                                   {}, limits, APPEARANCE);
  HS_EXPECT_EQ(COMPLETE.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(COMPLETE.trace.counters.layers, 2);
  HS_EXPECT_FALSE(COMPLETE.trace.has_surface);
  HS_EXPECT_NEAR(COMPLETE.trace.contribution.t, .8f, 1e-6f);
  limits.max_layers = 1;
  const auto PARTIAL = SDF::shade_periodic_shells(camera, {1, 0, 0}, 1, .3f, {},
                                                  limits, APPEARANCE);
  HS_EXPECT_EQ(PARTIAL.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(PARTIAL.trace.counters.layers, 1);
  HS_EXPECT_TRUE(PARTIAL.color.alpha > 0);
  HS_EXPECT_NEAR(PARTIAL.trace.contribution.t, .2f, 1e-6f);
  limits.max_steps = 0;
  const auto EMPTY = SDF::shade_periodic_shells(camera, {1, 0, 0}, 1, .3f, {},
                                                limits, APPEARANCE);
  HS_EXPECT_EQ(EMPTY.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(EMPTY.trace.counters.layers, 0);
  limits = {};
  camera.interval = {.5f, .9f};
  const auto CLIPPED = SDF::shade_periodic_shells(camera, {1, 0, 0}, 1, .3f, {},
                                                  limits, APPEARANCE);
  HS_EXPECT_EQ(CLIPPED.trace.counters.layers, 1);
  HS_EXPECT_NEAR(CLIPPED.trace.contribution.t, .8f, 1e-6f);
  camera.center = {{.5f, 0, 0, 0}};
  camera.interval = {0, 1};
  const auto REVERSE = SDF::shade_periodic_shells(camera, {-1, 0, 0}, 1, .3f,
                                                  {}, limits, APPEARANCE);
  HS_EXPECT_EQ(REVERSE.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(REVERSE.trace.counters.layers, 2);
  camera.center = {{-.5f, .305f, 0, 0}};
  const auto FILTERED = SDF::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {.04f, 0}, limits, APPEARANCE);
  HS_EXPECT_EQ(FILTERED.trace.counters.layers, 1);
  HS_EXPECT_FALSE(FILTERED.trace.contribution.verified);
  HS_EXPECT_NEAR(FILTERED.trace.contribution.coverage, .25f, 1e-5f);
  const auto UNFILTERED = SDF::shade_periodic_shells(camera, {1, 0, 0}, 1, .3f,
                                                     {}, limits, APPEARANCE);
  HS_EXPECT_EQ(UNFILTERED.trace.counters.layers, 0);

  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.center = {{-.5f, .245f, 0, .18f}};
  const Raycast::Footprint FOOTPRINT{.08f, 0};
  const auto SLICE_FILTERED = SDF::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, APPEARANCE);
  const float AMBIENT_DISTANCE = sqrtf(.245f * .245f + .18f * .18f);
  const float PROJECTED_DISTANCE =
      (AMBIENT_DISTANCE - .3f) * AMBIENT_DISTANCE / .245f;
  HS_EXPECT_EQ(SLICE_FILTERED.trace.counters.layers, 1);
  HS_EXPECT_FALSE(SLICE_FILTERED.trace.contribution.verified);
  HS_EXPECT_NEAR(SLICE_FILTERED.trace.contribution.coverage,
                 .5f - PROJECTED_DISTANCE / FOOTPRINT.at(.5f), 1e-5f);
  HS_EXPECT_NEAR(SLICE_FILTERED.trace.contribution.normal.y, 1, 1e-6f);

  camera.center = {{-.5f, .235f, 0, .18f}};
  const auto SLICE_HIT = SDF::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, APPEARANCE);
  const float HALF_CHORD = sqrtf(.3f * .3f - .18f * .18f - .235f * .235f);
  HS_EXPECT_EQ(SLICE_HIT.trace.counters.layers, 2);
  HS_EXPECT_TRUE(SLICE_HIT.trace.contribution.verified);
  HS_EXPECT_NEAR(SLICE_HIT.trace.contribution.t, .5f + HALF_CHORD, 1e-6f);
  HS_EXPECT_NEAR(SLICE_HIT.trace.contribution.coverage,
                 .5f + .5f * HALF_CHORD * HALF_CHORD /
                           (.24f * FOOTPRINT.at(.5f + HALF_CHORD)),
                 1e-5f);

  camera.center = {{-.5f, 0, 0, .305f}};
  const auto OFF_SLICE = SDF::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, APPEARANCE);
  HS_EXPECT_EQ(OFF_SLICE.trace.counters.layers, 0);
}

/** @brief Pins prepared shells match sphere roots. */
inline void test_prepared_shells_match_sphere_roots() {
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, ConstantPalette{});
  const Raycast::Appearance APPEARANCE{.1f, 0, 1, &palette.view()};
  static_assert(sizeof(SDF::PreparedPeriodicShells) <= 28);
  for (int sample = 0; sample < 120; ++sample) {
    Raycast::PreparedCamera camera;
    camera.domain = sample % 2 ? Raycast::SamplingDomain::SLICE_4D
                               : Raycast::SamplingDomain::SPATIAL_3D;
    camera.center = {{.04f, -.03f, .02f, sample % 2 ? .06f : 0}};
    camera.interval = {0, .49f};
    math::rotate_plane(camera.embedding, 0, 2, sample * .37f);
    if (sample % 2)
      math::rotate_plane(camera.embedding, 1, 3, sample * .21f);
    const SDF::PeriodicShells GEOMETRY{1, .12f + (sample % 6) * .05f};
    const auto PREPARED = SDF::prepare_periodic_shells(
        camera, GEOMETRY.cell_size, GEOMETRY.shell_radius, {});
    HS_EXPECT_TRUE(PREPARED.valid);
    const math::Vector VIEW =
        math::Vector{cosf(sample * .51f), sinf(sample * .41f), .43f}
            .normalized();
    const auto AMBIENT = camera.embedding.apply({{VIEW.x, VIEW.y, VIEW.z, 0}});
    const auto ROOTS =
        GEOMETRY.intersect(camera.center, AMBIENT, sample % 2 ? 4 : 3);
    const auto HIT =
        SDF::shade_periodic_shells(PREPARED, camera, VIEW, {}, APPEARANCE);
    HS_EXPECT_TRUE(ROOTS.hit);
    HS_EXPECT_EQ(HIT.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
    HS_EXPECT_EQ(HIT.trace.counters.layers, 1);
    HS_EXPECT_NEAR(HIT.trace.contribution.t, ROOTS.far, 1e-6f);
    HS_EXPECT_TRUE(HIT.trace.contribution.verified);
    HS_EXPECT_TRUE(HIT.trace.contribution.has_normal);
    HS_EXPECT_NEAR(HIT.trace.contribution.normal.magnitude(), 1.0f, 2e-6f);
  }
  Raycast::PreparedCamera invalid;
  invalid.embedding.m[0][0] = 2;
  const auto PREPARED = SDF::prepare_periodic_shells(invalid, 1, .3f, {});
  HS_EXPECT_FALSE(PREPARED.valid);
  HS_EXPECT_EQ(
      SDF::shade_periodic_shells(PREPARED, invalid, {1, 0, 0}, {}, APPEARANCE)
          .trace.status,
      Raycast::TraceStatus::INVALID_QUERY);
}

struct GradientPalette {
  Color4 get(float t) const {
    return {{static_cast<uint16_t>(9000 + 50000 * t),
             static_cast<uint16_t>(60000 - 45000 * t), 21000},
            1.0f};
  }
};

struct ShellMarchCase {
  float cell, radius, far, near_fade, angular, radial_start;
};

struct ShellRandom {
  uint32_t state;
  float operator()(float low, float high) {
    state = state * 1664525u + 1013904223u;
    return low + (high - low) * static_cast<float>(state >> 8) * 0x1p-24f;
  }
};

inline Raycast::PreparedCamera shell_camera(const ShellMarchCase &settings,
                                            ShellRandom &uniform,
                                            Raycast::SamplingDomain domain) {
  Raycast::PreparedCamera camera;
  camera.domain = domain;
  camera.center = {{uniform(-2, 2), uniform(-2, 2), uniform(-2, 2), 0}};
  const int dimensions = domain == Raycast::SamplingDomain::SLICE_4D ? 4 : 3;
  if (dimensions == 4)
    camera.center[3] = uniform(-2, 2);
  camera.radial_start = settings.radial_start;
  camera.interval = {0, settings.far};
  for (int a = 0; a < dimensions; ++a)
    for (int b = a + 1; b < dimensions; ++b)
      math::rotate_plane(camera.embedding, a, b, uniform(0, 6.3f));
  return camera;
}

inline SDF::ShellSample trace_single_owner_shell(
    const SDF::PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance,
    SDF::ShellLayerStorage &) {
  return SDF::trace_periodic_shells_3d(prepared, camera, direction, limits,
                                       appearance);
}

/** @brief Compares shell march composites against per-cell traversal over seeded rays. */
template <size_t Count, typename CameraSetup, typename CheckPrepared,
          typename Trace>
inline void compare_shell_march(const ShellMarchCase (&cases)[Count],
                                uint32_t seed, int frames,
                                CameraSetup camera_setup,
                                CheckPrepared check_prepared, Trace trace) {
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, GradientPalette{});
  ShellRandom uniform{seed};
  int compared = 0, hits = 0, mismatched = 0;
  for (const auto &settings : cases) {
    for (int frame = 0; frame < frames; ++frame) {
      const auto camera = camera_setup(settings, uniform);
      const auto prepared = SDF::prepare_periodic_shells(
          camera, settings.cell, settings.radius,
          {settings.angular, settings.radial_start});
      check_prepared(prepared);
      const Raycast::Appearance appearance{
          1 / settings.far, 0, 1 / settings.near_fade, &palette.view()};
      const Raycast::TraceLimits limits;
      SDF::ShellLayerStorage layers;
      for (int ray = 0; ray < 200; ++ray) {
        const math::Vector direction =
            math::Vector{uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)}
                .normalized();
        const auto cells = SDF::shade_periodic_shells(
            prepared, camera, direction, limits, appearance);
        const auto march =
            trace(prepared, camera, direction, limits, appearance, layers);
        const Pixel expected = cells.color.color * cells.color.alpha;
        HS_EXPECT_EQ(march.status, cells.trace.status);
        // Tangent-root rounding can choose either side of the silhouette jump.
        mismatched += !(abs(march.color.r - expected.r) <= 24 &&
                        abs(march.color.g - expected.g) <= 24 &&
                        abs(march.color.b - expected.b) <= 24);
        ++compared;
        hits += cells.trace.counters.layers > 0;
      }
    }
  }
  HS_EXPECT_TRUE(hits > compared / 20);
  HS_EXPECT_TRUE(mismatched <= compared / 2000);
}

/**
 * @brief The 3D layer march preserves traversal status and bounds composite error.
 * @details Channels differ by at most 24 codes except for at most 1/2000 rays
 *          at tangent silhouettes. Neighbor-reaching frames are declined.
 */
inline void test_shell_layer_march_matches_cell_traversal() {
  const ShellMarchCase CASES[] = {{.78625f, .1f, 10.736f, 2, .0218f, 0},
                                  {.4645f, .15f, 5.836f, .6f, .0218f, 0},
                                  {1, .28f, 4, .5f, 0, 0},
                                  {1, .25f, 3, .3f, .01f, .7f},
                                  {.6f, .2f, 6, 1, .01f, 0}};
  compare_shell_march(
      CASES, 0x5e11u, 24,
      [](const ShellMarchCase &settings, ShellRandom &random) {
        return shell_camera(settings, random,
                            Raycast::SamplingDomain::SPATIAL_3D);
      },
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_TRUE(prepared.single_owner);
      },
      trace_single_owner_shell);

  Raycast::PreparedCamera camera;
  camera.interval = {0, 16};
  HS_EXPECT_FALSE(SDF::prepare_periodic_shells(camera, .25f, .32f, {.0218f, 0})
                      .single_owner);
  HS_EXPECT_FALSE(
      SDF::prepare_periodic_shells(camera, 1, .3f, {}).single_owner);
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  HS_EXPECT_FALSE(
      SDF::prepare_periodic_shells(camera, 1, .1f, {}).single_owner);
}

/**
 * @brief The 3D neighbor march preserves traversal status and bounds composite error.
 * @details Channels differ by at most 24 codes except for at most 1/2000 rays
 *          at tangent silhouettes, including neighbor-reaching spheres.
 */
inline void test_shell_neighbor_march_matches_cell_traversal() {
  const ShellMarchCase CASES[] = {{.4645f, .15f, 7.5f, 1, .0218f, 0},
                                  {.7287f, .15f, 10.85f, 1, .0218f, 0},
                                  {1, .4f, 6, 1, .01f, 0}};
  compare_shell_march(
      CASES, 0x3d7e11u, 16,
      [](const ShellMarchCase &settings, ShellRandom &random) {
        return shell_camera(settings, random,
                            Raycast::SamplingDomain::SPATIAL_3D);
      },
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_FALSE(prepared.single_owner);
        HS_EXPECT_TRUE(prepared.march);
      },
      SDF::trace_periodic_shells_march<3>);
}

/**
 * @brief Shell marches take radial_start from the prepared footprint.
 */
inline void test_shell_march_uses_footprint_radial_start() {
  const ShellMarchCase OWNER[] = {{1, .1f, 4, .5f, .03f, 7}};
  const ShellMarchCase NEIGHBORS[] = {{1, .1f, 5, .5f, .04f, 12}};
  const auto camera_setup = [](const ShellMarchCase &settings,
                               ShellRandom &random) {
    auto camera =
        shell_camera(settings, random, Raycast::SamplingDomain::SPATIAL_3D);
    camera.radial_start = 0;
    return camera;
  };
  compare_shell_march(
      OWNER, 0xf007u, 16, camera_setup,
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_TRUE(prepared.single_owner);
      },
      trace_single_owner_shell);
  compare_shell_march(
      NEIGHBORS, 0xf008u, 16, camera_setup,
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_FALSE(prepared.single_owner);
        HS_EXPECT_TRUE(prepared.march);
      },
      SDF::trace_periodic_shells_march<3>);
  compare_shell_march(
      NEIGHBORS, 0xf009u, 16,
      [](const ShellMarchCase &settings, ShellRandom &random) {
        auto camera =
            shell_camera(settings, random, Raycast::SamplingDomain::SLICE_4D);
        camera.radial_start = 0;
        return camera;
      },
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_TRUE(prepared.march);
      },
      SDF::trace_periodic_shells_march<4>);
}

/**
 * @brief The 4D-slice march preserves traversal status and bounds composite error.
 * @details Channels differ by at most 24 codes except for at most 1/2000 rays
 *          at tangent silhouettes.
 */
inline void test_shell_slice_march_matches_cell_traversal() {
  const ShellMarchCase CASES[] = {{1, .15f, 16, 2, .0218f, 0},
                                  {1, .3f, 6, .5f, .0218f, 0},
                                  {.6f, .2f, 8, 1, .01f, .4f},
                                  {1, .45f, 4, .5f, 0, 0}};
  compare_shell_march(
      CASES, 0x4d511u, 24,
      [](const ShellMarchCase &settings, ShellRandom &random) {
        return shell_camera(settings, random,
                            Raycast::SamplingDomain::SLICE_4D);
      },
      [](const SDF::PreparedPeriodicShells &prepared) {
        HS_EXPECT_TRUE(prepared.march);
        HS_EXPECT_FALSE(prepared.single_owner);
      },
      SDF::trace_periodic_shells_march<4>);

  Raycast::PreparedCamera camera;
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.interval = {0, 16};
  HS_EXPECT_FALSE(
      SDF::prepare_periodic_shells(camera, .5f, .3f, {.0218f, 0}).march);
}

inline int run_ray_demonstrator_tests() {
  hs_test::ModuleFixture fixture("ray_demonstrators");
  test_octet_crossing_coverage_against_ray_line_distance();
  test_affine_lattice_world_metric_and_periods();
  test_affine_cached_metric_against_ray_line_oracle();
  test_periodic_shell_roots_and_slices();
  test_periodic_shell_traversal_budgets();
  test_prepared_shells_match_sphere_roots();
  test_shell_layer_march_matches_cell_traversal();
  test_shell_slice_march_matches_cell_traversal();
  test_shell_neighbor_march_matches_cell_traversal();
  test_shell_march_uses_footprint_radial_start();
  test_framework_generic_event_rendering();
  test_repeated_stream_grouping_preserves_order_and_endpoints();
  test_lattice_volume_camera_demonstrators();
  test_torus_and_warped_volume_spherical_stage();
  test_verified_filter_off_slice_geometry_and_projected_normal();
  return fixture.result();
}

} // namespace hs_test::ray_demonstrator_tests
