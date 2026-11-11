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
#include "core/render/sdf/lattice_trace.h"
#include "core/render/sdf/periodic_shells.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"
#include "tests/pixel_test_util.h"

#include <array>
#include <utility>
#include <vector>

namespace hs_test::ray_demonstrator_tests {

template <typename Source> struct PaletteFixture {
  alignas(
      Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()> buffer;
  Arena arena{buffer.data(), buffer.size()};
  BakedPaletteStorage palette;

  PaletteFixture() { palette.bake(arena, Source{}); }
};

/** @brief Palette sampling requires the appearance to cover the ray interval. */
inline void test_appearance_palette_interval() {
  Raycast::Appearance appearance;
  appearance.inv_far = 0.1f;
  HS_EXPECT_TRUE(appearance.valid_for({0.0f, 10.0f}));
  HS_EXPECT_FALSE(appearance.valid_for({0.0f, 11.0f}));
  appearance.inv_far = 0.0f;
  HS_EXPECT_TRUE(appearance.valid_for({0.0f, 11.0f}));
  appearance.inv_far = -0.1f;
  HS_EXPECT_FALSE(appearance.valid_for({0.0f, 10.0f}));
}

/** @brief Plane crossings evaluated in world-distance units. */
struct PlaneCrossingReference {
  struct Cursor {
    float next = 0.0f;
    float step = 0.0f;
    bool active = false;
  };
  std::array<Cursor, 4> cursors{};

  PlaneCrossingReference(const Raycast::Ray &ray,
                         const std::array<SDF::FrameworkPlane, 4> &families,
                         const math::Vector &origin) {
    const math::Vector start = ray.at(ray.interval.near) - origin;
    for (size_t i = 0; i < cursors.size(); ++i) {
      const float speed = math::dot(ray.direction, families[i].normal);
      if (speed == 0.0f)
        continue;
      const float position = math::dot(start, families[i].normal);
      const float spacing = families[i].spacing;
      const float cell = position / spacing;
      const float plane = speed > 0.0f ? ceilf(cell) : floorf(cell);
      auto &cursor = cursors[i];
      cursor.next = ray.interval.near + (plane * spacing - position) / speed;
      cursor.step = spacing / fabsf(speed);
      cursor.active = Raycast::finite(cursor.next) &&
                      Raycast::finite(cursor.step) && cursor.step > 0.0f;
    }
  }

  bool active(size_t index) const { return cursors[index].active; }
  float distance(size_t index) const { return cursors[index].next; }
  void advance(size_t index) {
    auto &cursor = cursors[index];
    const float next = cursor.next + cursor.step;
    cursor.active = Raycast::finite(next) && next > cursor.next;
    cursor.next = next;
  }
};

/** @brief Octet crossing coverage matches an independent ray-to-line distance reference. */
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
    const Raycast::Ray ray{{sinf(sample * .43f) * 5, cosf(sample * .71f) * 5,
                            sinf(sample * .19f) * 5},
                           math::Vector{cosf(sample * .31f),
                                        sinf(sample * .57f),
                                        cosf(sample * .23f)}
                               .normalized(),
                           {.1f, 16.0f}};
    const Raycast::Footprint footprint{sample % 5 == 0 ? 0.0f : .02f, .3f};
    const auto planes = geometry.plane_families();
    SDF::OctetEvents events(geometry, ray, footprint);
    PlaneCrossingReference reference(ray, planes, geometry.origin);
    HS_EXPECT_TRUE(geometry.valid());
    std::array<bool, SDF::OctetEvents::STREAM_COUNT> owns_pair{};
    for (size_t i = 0; i < planes.size(); ++i) {
      for (size_t j = i + 1; j < planes.size(); ++j) {
        const float a = math::dot(ray.direction, planes[i].normal);
        const float b = math::dot(ray.direction, planes[j].normal);
        const size_t owner = fabsf(a) >= fabsf(b) ? i : j;
        owns_pair[owner] = owns_pair[owner] || (owner == i ? a : b) != 0;
      }
    }
    for (size_t stream = 0; stream < planes.size(); ++stream)
      reference.cursors[stream].active &= owns_pair[stream];
    for (size_t stream = 0; stream < planes.size(); ++stream) {
      HS_EXPECT_EQ(events.active(stream), reference.active(stream));
      for (int crossing = 0; crossing < 8 && events.active(stream);
           ++crossing) {
        ++compared;
        const auto hit = events.candidate(stream);
        HS_EXPECT_NEAR(hit.t, reference.distance(stream),
                       2e-5f * std::max(1.0f, hit.t));
        const auto point = ray.at(hit.t) - geometry.origin;
        float best = INFINITY;
        uint32_t feature = 0;
        uint32_t pair = 0;
        for (size_t i = 0; i < planes.size(); ++i)
          for (size_t j = i + 1; j < planes.size(); ++j, ++pair) {
            const float a = math::dot(ray.direction, planes[i].normal);
            const float b = math::dot(ray.direction, planes[j].normal);
            const size_t owner = fabsf(a) >= fabsf(b) ? i : j;
            if (owner != stream)
              continue;
            const size_t other = owner == i ? j : i;
            const float u = math::dot(point, planes[other].normal);
            const float residual =
                u - planes[other].spacing * roundf(u / planes[other].spacing);
            const auto cross = math::cross(
                ray.direction, math::cross(planes[i].normal, planes[j].normal));
            const float distance =
                fabsf(residual * (owner == i ? a : b)) / cross.magnitude();
            if (distance < best) {
              best = distance;
              feature = pair;
            }
          }
        const float width = footprint.at(hit.t);
        const float coverage =
            width > 0 ? std::clamp(.5f - (best - geometry.wire_radius) / width,
                                   0.0f, 1.0f)
            : best <= geometry.wire_radius ? 1.0f
                                           : 0.0f;
        HS_EXPECT_NEAR(hit.coverage, coverage, 3e-4f);
        if (coverage > 0) {
          ++hits;
          HS_EXPECT_EQ(hit.feature, feature);
        }
        events.advance(stream);
        reference.advance(stream);
        HS_EXPECT_EQ(events.active(stream), reference.active(stream));
      }
    }
  }
  HS_EXPECT_GT(compared, 1000);
  HS_EXPECT_GT(hits, compared / 20);
}

/** @brief Framework events emit ordered crossings with expected distance, coverage and verification state. */
inline void test_framework_generic_event_rendering() {
  SDF::TriangularFramework geometry;
  Raycast::TraceLimits limits;
  limits.max_candidates = 64;
  limits.max_layers = 32;
  const Raycast::Ray RAY{{0.0f, 0.0f, 0.3f}, {0.0f, 0.0f, -1.0f}, {0.0f, 2.3f}};
  SDF::FrameworkEvents events(geometry, RAY);
  std::array<Raycast::Contribution, 8> output;
  size_t count = 0;
  const auto result = Raycast::trace_events(
      events, RAY.interval, limits, [&](const Raycast::Contribution &hit) {
        if (count < output.size())
          output[count] = hit;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(count, size_t{3});
  for (size_t i = 0; i < count && i < output.size(); ++i) {
    HS_EXPECT_NEAR(output[i].t, 0.3f + static_cast<float>(i), 1e-6f);
    HS_EXPECT_EQ(output[i].coverage, 1.0f);
    HS_EXPECT_FALSE(output[i].verified);
  }
  const math::Vector direction = math::Vector{0.7f, 0.4f, 0.3f}.normalized();
  const Raycast::Ray oblique{{0.0f, 0.0f, 0.0f}, direction, {0.0f, 5.0f}};
  SDF::FrameworkEvents slanted(geometry, oblique, {0.05f, 0.2f});
  for (size_t i = 0; i < SDF::FrameworkEvents::STREAM_COUNT; ++i)
    HS_EXPECT_TRUE(slanted.active(i));
  float previous = -1.0f;
  count = 0;
  const auto slanted_value = Raycast::trace_events(
      slanted, oblique.interval, limits, [&](const Raycast::Contribution &hit) {
        HS_EXPECT_GE(hit.t, previous);
        previous = hit.t;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(slanted_value.status, Raycast::TraceStatus::RANGE_COMPLETE);
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
    result.merge_identity = i == 1 || i == 2 ? 0 : i;
    result.material = i == 2 ? 1 : 0;
    result.coverage = i == 1 ? 0.8f : 0.3f;
    return result;
  }
  void advance(size_t i) { ++indices[i]; }
};

/** @brief Repeated event groups preserve sorted output and terminal endpoints. */
inline void test_repeated_stream_grouping_preserves_order_and_endpoints() {
  CloseStreams events;
  Raycast::TraceLimits limits;
  std::array<Raycast::Contribution, 8> output;
  size_t count = 0;
  const auto result = Raycast::trace_events(
      events, {0.0f, 1.0f}, limits, [&](const Raycast::Contribution &hit) {
        if (count < output.size())
          output[count] = hit;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(result.counters.candidates, 6);
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

/** @brief Lattice cameras resolve analytic entry and exit distances and reject empty 4D slices. */
inline void test_lattice_volume_camera_demonstrators() {
  Raycast::PreparedCamera camera;
  camera.center = {{0.3f, 0.2f, 0.0f, 0.0f}};
  camera.interval = {0.0f, 1.0f};
  SDF::WireLattice<3> cubic;
  Raycast::DomainQuery3<decltype(cubic)> query3{cubic, camera};
  Raycast::TraceLimits limits;
  limits.max_queries = 256;
  limits.max_steps = 256;
  limits.max_refinements = 256;
  limits.position_tolerance = 1e-5f;
  const math::Vector DIRECTION{0.0f, -1.0f, 0.0f};
  const auto cubic_value =
      Raycast::surface_search(query3, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_TRUE(cubic_value.has_surface);
  HS_EXPECT_TRUE(cubic_value.contribution.verified);
  HS_EXPECT_NEAR(cubic_value.contribution.t, 0.15f, 1e-5f);

  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.center[3] = 0.03f;
  SDF::WireLattice<4> hypercubic;
  Raycast::DomainQuery4<decltype(hypercubic)> query4{hypercubic, camera};
  const auto hypercubic_value =
      Raycast::surface_search(query4, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_TRUE(hypercubic_value.has_surface);
  HS_EXPECT_NEAR(hypercubic_value.contribution.t, 0.16f, 1e-5f);
  camera.center[1] = 0.0f;
  const auto exit = Raycast::surface_search(
      query4, camera.ray(math::Vector{0.0f, 1.0f, 0.0f}), {}, limits);
  HS_EXPECT_TRUE(exit.has_surface);
  HS_EXPECT_NEAR(exit.contribution.t, 0.04f, 1e-5f);
  camera.center[3] = 0.1f;
  const auto empty_slice =
      Raycast::surface_search(query4, camera.ray(DIRECTION), {}, limits);
  HS_EXPECT_FALSE(empty_slice.has_surface);
  HS_EXPECT_EQ(empty_slice.status, Raycast::TraceStatus::RANGE_COMPLETE);
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
    Raycast::DomainQuery3<decltype(placed)> domain{placed, prepared.camera};
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
    const auto prepared = Pipeline::prepare_stages(frame);
    HS_EXPECT_EQ(frame.preparations, pose + 1);
    uint64_t hash = hs_test::HASH_SEED64;
    int hits = 0;
    int misses = 0;
    for (int y = -6; y <= 6; ++y) {
      for (int x = -6; x <= 6; ++x) {
        const math::Vector direction =
            math::Vector{x * 0.08f, y * 0.08f, 1.0f}.normalized();
        const Color4 color = Pipeline::evaluate(direction, frame, prepared);
        const auto trace =
            Renderer::trace(direction, frame, std::get<0>(prepared));
        HS_EXPECT_EQ(color.color.r, trace.color.color.r);
        HS_EXPECT_EQ(color.alpha, trace.color.alpha);
        HS_EXPECT_TRUE(Raycast::finite(color.alpha));
        HS_EXPECT_GE(color.alpha, 0.0f);
        HS_EXPECT_LE(color.alpha, 1.0f);
        if (trace.trace.has_surface) {
          ++hits;
          HS_EXPECT_TRUE(trace.trace.contribution.verified);
          HS_EXPECT_GT(color.alpha, 0.0f);
        } else {
          ++misses;
          HS_EXPECT_EQ(color.alpha, 0.0f);
        }
        hash = hs_test::fnv1a64_channel(hash, color.color.r);
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

/** @brief Placed torus stages preserve direct trace shading and respond to pose changes. */
inline void test_torus_and_warped_volume_spherical_stage() {
  PaletteFixture<ConstantPalette> fixture;
  auto &palette = fixture.palette;
  const SDF::Torus TORUS{1.0f, 0.25f};
  check_placed_volume_stage(TORUS, palette.view());
  const SDF::WarpedVolume warped{TORUS, SDF::Warp::Twist{3, 0.15f, 1.0f},
                                 0.001f};
  check_placed_volume_stage(warped, palette.view());
}

struct OffSliceBall {
  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }
  Raycast::QuerySample sample(const math::Vec4 &p) const {
    const float w = p[3] - 1.001f;
    const float value =
        sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + w * w) - 1.0f;
    return {value, fabsf(value), value == 0.0f, 0, 0};
  }
};

/** @brief Verified filtering suppresses off-slice geometry and preserves partial-hit coverage and failure state. */
inline void test_verified_filter_off_slice_geometry_and_projected_normal() {
  Raycast::PreparedCamera camera;
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.interval = {0.0f, 0.1f};
  OffSliceBall ball;
  Raycast::DomainQuery4<decltype(ball)> query{ball, camera};
  constexpr std::array<math::Vector, 4> DIRECTIONS = {
      math::Vector{1.0f, 0.0f, 0.0f}, math::Vector{-1.0f, 0.0f, 0.0f},
      math::Vector{0.0f, 1.0f, 0.0f}, math::Vector{0.0f, -1.0f, 0.0f}};
  const auto filtered =
      Raycast::verified_filter(DIRECTIONS, [&](const auto &dir) {
        Raycast::ShadedTrace sample;
        sample.trace = Raycast::surface_search(query, camera.ray(dir), {}, {});
        sample.color = {{65535, 65535, 65535}, 1.0f};
        return sample;
      });
  HS_EXPECT_FALSE(filtered.trace.has_surface);
  HS_EXPECT_EQ(filtered.color.alpha, 0.0f);
  math::Vector projected{1.0f, 1.0f, 1.0f};
  HS_EXPECT_FALSE(camera.project_normal({{0.0f, 0.0f, 0.0f, 1.0f}}, projected));
  HS_EXPECT_EQ(projected.magnitude(), 0.0f);
  const auto partial =
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
  HS_EXPECT_EQ(partial.trace.status, Raycast::TraceStatus::UNRESOLVED);
  HS_EXPECT_TRUE(partial.trace.has_surface);
  HS_EXPECT_TRUE(partial.trace.contribution.verified);
  HS_EXPECT_NEAR(partial.trace.contribution.coverage, 0.25f, 1e-7f);
  HS_EXPECT_NEAR(partial.color.alpha, 0.2f, 1e-7f);
  HS_EXPECT_EQ(partial.color.color.r, uint16_t{40000});
}

/** @brief Periodic-shell roots match analytic sphere chords and reject disjoint slices. */
inline void test_periodic_shell_roots_and_slices() {
  SDF::PeriodicShells geometry{2, .3f};
  const auto sphere = geometry.intersect({{-2, 0, 0, 0}}, {{1, 0, 0, 0}}, 3);
  HS_EXPECT_TRUE(sphere.hit);
  HS_EXPECT_NEAR(sphere.near, 1.4f, 1e-6f);
  HS_EXPECT_NEAR(sphere.far, 2.6f, 1e-6f);
  const auto slice = geometry.intersect({{-2, 0, 0, .36f}}, {{1, 0, 0, 0}}, 4);
  HS_EXPECT_NEAR(slice.near, 1.52f, 1e-6f);
  HS_EXPECT_NEAR(slice.far, 2.48f, 1e-6f);
  HS_EXPECT_FALSE(geometry.intersect({{-2, 0, 0, .7f}}, {{1, 0, 0, 0}}, 4).hit);
  const auto inside = geometry.intersect({{0, 0, 0, 0}}, {{1, 0, 0, 0}}, 4);
  HS_EXPECT_NEAR(inside.near, -.6f, 1e-6f);
  HS_EXPECT_NEAR(inside.far, .6f, 1e-6f);
}

/** @brief Shell traversal clips intervals, respects budgets and projects slice antialias coverage. */
inline void test_periodic_shell_traversal_budgets() {
  PaletteFixture<ConstantPalette> fixture;
  auto &palette = fixture.palette;
  const Raycast::Appearance appearance{.1f, 0, 1, &palette.view()};
  Raycast::PreparedCamera camera;
  camera.center = {{-.5f, 0, 0, 0}};
  camera.interval = {0, 1};
  Raycast::TraceLimits limits;
  const auto complete = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(complete.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(complete.trace.counters.layers, 2);
  HS_EXPECT_FALSE(complete.trace.has_surface);
  HS_EXPECT_NEAR(complete.trace.contribution.t, .8f, 1e-6f);
  limits.max_layers = 1;
  const auto partial = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(partial.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(partial.trace.counters.layers, 1);
  HS_EXPECT_GT(partial.color.alpha, 0);
  HS_EXPECT_NEAR(partial.trace.contribution.t, .2f, 1e-6f);
  limits.max_steps = 0;
  const auto empty = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(empty.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(empty.trace.counters.layers, 0);
  limits = {};
  camera.interval = {.5f, .9f};
  const auto clipped = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(clipped.trace.counters.layers, 1);
  HS_EXPECT_NEAR(clipped.trace.contribution.t, .8f, 1e-6f);
  camera.center = {{.5f, 0, 0, 0}};
  camera.interval = {0, 1};
  const auto reverse = SDF::detail::reference::shade_periodic_shells(
      camera, {-1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(reverse.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(reverse.trace.counters.layers, 2);
  camera.center = {{-.5f, .305f, 0, 0}};
  const auto filtered = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {.04f, 0}, limits, appearance);
  HS_EXPECT_EQ(filtered.trace.counters.layers, 1);
  HS_EXPECT_FALSE(filtered.trace.contribution.verified);
  HS_EXPECT_NEAR(filtered.trace.contribution.coverage, .25f, 1e-5f);
  const auto unfiltered = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, {}, limits, appearance);
  HS_EXPECT_EQ(unfiltered.trace.counters.layers, 0);

  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  camera.center = {{-.5f, .245f, 0, .18f}};
  const Raycast::Footprint FOOTPRINT{.08f, 0};
  const auto slice_filtered = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, appearance);
  const float ambient_distance = sqrtf(.245f * .245f + .18f * .18f);
  const float projected_distance =
      (ambient_distance - .3f) * ambient_distance / .245f;
  HS_EXPECT_EQ(slice_filtered.trace.counters.layers, 1);
  HS_EXPECT_FALSE(slice_filtered.trace.contribution.verified);
  HS_EXPECT_NEAR(slice_filtered.trace.contribution.coverage,
                 .5f - projected_distance / FOOTPRINT.at(.5f), 1e-5f);
  HS_EXPECT_NEAR(slice_filtered.trace.contribution.normal.y, 1, 1e-6f);

  camera.center = {{-.5f, .235f, 0, .18f}};
  const auto slice_hit = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, appearance);
  const float half_chord = sqrtf(.3f * .3f - .18f * .18f - .235f * .235f);
  HS_EXPECT_EQ(slice_hit.trace.counters.layers, 2);
  HS_EXPECT_TRUE(slice_hit.trace.contribution.verified);
  HS_EXPECT_NEAR(slice_hit.trace.contribution.t, .5f + half_chord, 1e-6f);
  HS_EXPECT_NEAR(slice_hit.trace.contribution.coverage,
                 .5f + .5f * half_chord * half_chord /
                           (.24f * FOOTPRINT.at(.5f + half_chord)),
                 1e-5f);

  camera.center = {{-.5f, 0, 0, .305f}};
  const auto off_slice = SDF::detail::reference::shade_periodic_shells(
      camera, {1, 0, 0}, 1, .3f, FOOTPRINT, limits, appearance);
  HS_EXPECT_EQ(off_slice.trace.counters.layers, 0);
}

/** @brief Prepared shell hits match analytic roots with unit normals and reject invalid embeddings. */
inline void test_prepared_shells_match_sphere_roots() {
  PaletteFixture<ConstantPalette> fixture;
  auto &palette = fixture.palette;
  const Raycast::Appearance appearance{.1f, 0, 1, &palette.view()};
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
    const SDF::PeriodicShells geometry{1, .12f + (sample % 6) * .05f};
    const auto prepared = SDF::prepare_periodic_shells(
        camera, geometry.cell_size, geometry.shell_radius, {});
    HS_EXPECT_TRUE(prepared.valid);
    const math::Vector view =
        math::Vector{cosf(sample * .51f), sinf(sample * .41f), .43f}
            .normalized();
    const auto ambient = camera.embedding.apply({{view.x, view.y, view.z, 0}});
    const auto roots =
        geometry.intersect(camera.center, ambient, sample % 2 ? 4 : 3);
    const auto hit = SDF::detail::reference::shade_periodic_shells(
        prepared, camera, view, {}, appearance);
    HS_EXPECT_TRUE(roots.hit);
    HS_EXPECT_EQ(hit.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
    HS_EXPECT_EQ(hit.trace.counters.layers, 1);
    HS_EXPECT_NEAR(hit.trace.contribution.t, roots.far, 1e-6f);
    HS_EXPECT_TRUE(hit.trace.contribution.verified);
    HS_EXPECT_TRUE(hit.trace.contribution.has_normal);
    HS_EXPECT_NEAR(hit.trace.contribution.normal.magnitude(), 1.0f, 2e-6f);
  }
  Raycast::PreparedCamera invalid;
  invalid.embedding.m[0][0] = 2;
  const auto prepared = SDF::prepare_periodic_shells(invalid, 1, .3f, {});
  HS_EXPECT_FALSE(prepared.valid);
  HS_EXPECT_EQ(SDF::detail::reference::shade_periodic_shells(
                   prepared, invalid, {1, 0, 0}, {}, appearance)
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

inline Raycast::PreparedCamera shell_camera(const ShellMarchCase &settings,
                                            hs::Pcg32 &random,
                                            Raycast::SamplingDomain domain) {
  Raycast::PreparedCamera camera;
  camera.domain = domain;
  const float x = rand_uniform(random, -2.0f, 2.0f);
  const float y = rand_uniform(random, -2.0f, 2.0f);
  const float z = rand_uniform(random, -2.0f, 2.0f);
  camera.center = {{x, y, z, 0}};
  const int dimensions = domain == Raycast::SamplingDomain::SLICE_4D ? 4 : 3;
  if (dimensions == 4)
    camera.center[3] = rand_uniform(random, -2.0f, 2.0f);
  camera.radial_start = settings.radial_start;
  camera.interval = {0, settings.far};
  for (int a = 0; a < dimensions; ++a)
    for (int b = a + 1; b < dimensions; ++b) {
      const float angle = rand_uniform(random, 0.0f, 6.3f);
      math::rotate_plane(camera.embedding, a, b, angle);
    }
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
  PaletteFixture<GradientPalette> fixture;
  auto &palette = fixture.palette;
  hs::Pcg32 random{seed};
  int compared = 0, hits = 0, mismatched = 0;
  for (const auto &settings : cases) {
    for (int frame = 0; frame < frames; ++frame) {
      const auto camera = camera_setup(settings, random);
      const auto prepared = SDF::prepare_periodic_shells(
          camera, settings.cell, settings.radius,
          {settings.angular, settings.radial_start});
      check_prepared(prepared);
      const Raycast::Appearance appearance{
          1 / settings.far, 0, 1 / settings.near_fade, &palette.view()};
      const Raycast::TraceLimits limits;
      SDF::ShellLayerStorage layers;
      for (int ray = 0; ray < 200; ++ray) {
        const float x = rand_uniform(random, -1.0f, 1.0f);
        const float y = rand_uniform(random, -1.0f, 1.0f);
        const float z = rand_uniform(random, -1.0f, 1.0f);
        const math::Vector direction = math::Vector{x, y, z}.normalized();
        const auto cells = SDF::detail::reference::shade_periodic_shells(
            prepared, camera, direction, limits, appearance);
        const auto march =
            trace(prepared, camera, direction, limits, appearance, layers);
        const Pixel expected = cells.color.color * cells.color.alpha;
        HS_EXPECT_EQ(march.status, cells.trace.status);
        // Tangent-root rounding can choose either side of the silhouette jump.
        mismatched += max_channel_gap(march.color, expected) > 24;
        ++compared;
        hits += cells.trace.counters.layers > 0;
      }
    }
  }
  HS_EXPECT_GT(hits, compared / 20);
  HS_EXPECT_LE(mismatched, compared / 2000);
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
      [](const ShellMarchCase &settings, hs::Pcg32 &random) {
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
      [](const ShellMarchCase &settings, hs::Pcg32 &random) {
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
                               hs::Pcg32 &random) {
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
      [](const ShellMarchCase &settings, hs::Pcg32 &random) {
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
      [](const ShellMarchCase &settings, hs::Pcg32 &random) {
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

/** @brief One covered sphere contribution along a ray. */
struct ShellContribution {
  float t, coverage;
};

/**
 * @brief Every shell contribution to a ray, from all lattice centers in a box
 *        around the ray segment, sorted by distance.
 */
inline std::vector<ShellContribution>
shell_neighborhood(const Raycast::PreparedCamera &camera,
                   const math::Vector &direction, float cell, float radius,
                   Raycast::Footprint footprint) {
  const int dimensions =
      camera.domain == Raycast::SamplingDomain::SLICE_4D ? 4 : 3;
  const auto &e = camera.embedding.m;
  const math::Vec4 d =
      camera.embedding.apply({{direction.x, direction.y, direction.z, 0}});
  const math::Vec4 o = camera.point4(direction * camera.radial_start);
  const float near = camera.interval.near, far = camera.interval.far;
  const float r = radius * cell;
  const float margin = r + footprint.at(far) + cell;
  int low[4] = {}, high[4] = {};
  for (int k = 0; k < dimensions; ++k) {
    const float a = o[k] + near * d[k], b = o[k] + far * d[k];
    low[k] = static_cast<int>(floorf((fminf(a, b) - margin) / cell));
    high[k] = static_cast<int>(ceilf((fmaxf(a, b) + margin) / cell));
  }
  const auto projected = [&](const float (&g)[4]) {
    float length2 = 0;
    for (int j = 0; j < 3; ++j) {
      float component = 0;
      for (int k = 0; k < dimensions; ++k)
        component += e[k][j] * g[k];
      length2 += component * component;
    }
    return sqrtf(length2);
  };
  std::vector<ShellContribution> found;
  int index[4] = {low[0], low[1], low[2], low[3]};
  for (;;) {
    float g0[4] = {}, a = 0, b = 0, c = -r * r;
    for (int k = 0; k < dimensions; ++k) {
      g0[k] = o[k] - index[k] * cell;
      a += d[k] * d[k];
      b += g0[k] * d[k];
      c += g0[k] * g0[k];
    }
    const float discriminant = b * b - a * c;
    if (discriminant >= 0) {
      const float half_chord = sqrtf(discriminant) / a;
      const float center = -b / a;
      for (const float T : {center - half_chord, center + half_chord}) {
        if (T < near || T > far ||
            (T == center + half_chord && half_chord == 0))
          continue;
        float g[4] = {};
        for (int k = 0; k < dimensions; ++k)
          g[k] = g0[k] + T * d[k];
        const float width = footprint.at(T);
        const float length = projected(g);
        const float depth =
            length > 0 ? half_chord * half_chord / (2 * length) : 0;
        found.push_back(
            {T, width > 0 ? hs::clamp(.5f + depth / width, 0.0f, 1.0f) : 1.0f});
      }
    } else if (footprint.angular_radius > 0) {
      const float t = -b / a;
      float g[4] = {}, length2 = 0;
      for (int k = 0; k < dimensions; ++k) {
        g[k] = g0[k] + t * d[k];
        length2 += g[k] * g[k];
      }
      const float length = sqrtf(length2);
      const float projected_value = projected(g);
      const float width = footprint.at(t);
      if (t >= near && t <= far && width > 0 && projected_value > 0) {
        const float distance = (length - r) * length / projected_value;
        const float coverage = hs::clamp(.5f - distance / width, 0.0f, 1.0f);
        if (coverage > 0)
          found.push_back({t, coverage});
      }
    }
    int k = 0;
    for (; k < dimensions && ++index[k] > high[k]; ++k)
      index[k] = low[k];
    if (k == dimensions)
      break;
  }
  std::sort(found.begin(), found.end(),
            [](const ShellContribution &x, const ShellContribution &y) {
              return x.t < y.t;
            });
  return found;
}

/**
 * @brief Front-to-back composite of the first @p count contributions.
 * @return Premultiplied color, and whether the layer budget or saturation
 *         stopped it before @p count.
 */
inline std::pair<Pixel, bool>
composite_shells(const std::vector<ShellContribution> &found, size_t count,
                 const Raycast::Appearance &appearance,
                 const Raycast::TraceLimits &limits) {
  LayerComposite composite;
  for (size_t i = 0; i < count; ++i) {
    if (static_cast<int>(i) >= limits.max_layers ||
        static_cast<int>(i) >= limits.max_candidates)
      return {composite.premultiplied(), true};
    appearance.composite(composite, found[i].t, found[i].coverage);
    if (composite.saturated())
      return {composite.premultiplied(), true};
  }
  return {composite.premultiplied(), false};
}

/**
 * @brief A filtered ray between two sphere rows gathers both rows' silhouettes
 *        through the lattice trace.
 */
inline void test_wide_shell_footprint_gathers_neighbor_rows() {
  PaletteFixture<GradientPalette> fixture;
  auto &palette = fixture.palette;
  SDF::ShellLayerStorage layers;
  SDF::LatticeTrace::Settings settings;
  settings.geometry = SDF::LatticeTrace::Geometry::SHELLS;
  settings.cell_size = 1;
  settings.shell_radius = .32f;
  settings.far_distance = 8;
  settings.near_fade = 1;
  settings.pixel_half_angle = .1f;
  settings.center = {{-2, .5f, 0, 0}};
  settings.palette = &palette.view();
  settings.shell_layers = &layers;
  const auto prepared = SDF::LatticeTrace::prepare(settings);
  HS_EXPECT_TRUE(prepared.valid);
  HS_EXPECT_FALSE(prepared.periodic_shells.march);
  const math::Vector VIEW{1, 0, 0};
  const auto found =
      shell_neighborhood(prepared.camera, VIEW, 1, .32f, prepared.footprint);
  HS_EXPECT_EQ(found.size(), size_t{10});
  LayerComposite expected;
  for (const auto &hit : found)
    expected.add(prepared.appearance.color(hit.t) * .1f,
                 hit.coverage * prepared.appearance.opacity(hit.t));
  const auto traced = SDF::LatticeTrace::shade<false>(VIEW, prepared);
  HS_EXPECT_EQ(traced.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_LE(max_channel_gap(traced.color, expected.premultiplied()), 2);
  HS_EXPECT_EQ(SDF::detail::reference::shade_periodic_shells(
                   prepared.periodic_shells, prepared.camera, VIEW,
                   prepared.limits, prepared.appearance)
                   .trace.status,
               Raycast::TraceStatus::INVALID_QUERY);
}

/**
 * @brief The wide shell march gathers a silhouette whose center lies one
 *        layer behind the interval start.
 */
inline void test_wide_shell_march_gathers_layers_behind_start() {
  PaletteFixture<GradientPalette> fixture;
  auto &palette = fixture.palette;
  const math::Vector view = math::Vector{1, 1, 1}.normalized();
  const math::Vector START{-.49f, -.3415f, -.3415f};
  const Raycast::Footprint FOOTPRINT{.1f, 8};
  Raycast::PreparedCamera camera;
  camera.radial_start = FOOTPRINT.radial_start;
  camera.interval = {0, 8};
  const math::Vector center = START - view * camera.radial_start;
  camera.center = {{center.x, center.y, center.z, 0}};
  const auto prepared =
      SDF::prepare_periodic_shells(camera, 1, .32f, FOOTPRINT);
  HS_EXPECT_FALSE(prepared.march);
  const Raycast::Appearance appearance{1 / 8.0f, 0, 50, &palette.view()};
  Raycast::TraceLimits limits;
  limits.max_steps = 256;
  const auto found = shell_neighborhood(camera, view, 1, .32f, FOOTPRINT);
  HS_EXPECT_FALSE(found.empty());
  if (found.empty())
    return;
  HS_EXPECT_NEAR(found.front().t, .1f, .01f);
  SDF::ShellLayerStorage layers;
  const auto traced = SDF::trace_periodic_shells_march<3, true>(
      prepared, camera, view, limits, appearance, layers);
  const auto expected =
      composite_shells(found, found.size(), appearance, limits);
  HS_EXPECT_LE(max_channel_gap(traced.color, expected.first), 2);
}

/**
 * @brief The wide shell march matches a brute-force neighborhood composite in
 *        3D and in a 4D slice.
 * @details No ray may differ by more than 24 codes; budget-cut rays match
 *          a distance-ordered prefix of the neighborhood.
 */
inline void test_wide_shell_march_matches_neighborhood() {
  PaletteFixture<GradientPalette> fixture;
  auto &palette = fixture.palette;
  const ShellMarchCase CASES[] = {{1, .32f, 8, 1, .1f, 0},
                                  {.25f, .32f, 6, 1, .0218f, .5f},
                                  {.5f, .2f, 10, 2, .06f, 0},
                                  {1, .15f, 12, 1, .2f, 1},
                                  {1, .2f, 4, .02f, .1f, 12}};
  hs::Pcg32 random{0x91de5u};
  for (const auto domain : {Raycast::SamplingDomain::SPATIAL_3D,
                            Raycast::SamplingDomain::SLICE_4D}) {
    int compared = 0, hits = 0, mismatched = 0, truncated = 0;
    for (const auto &settings : CASES) {
      for (int frame = 0; frame < 6; ++frame) {
        const auto camera = shell_camera(settings, random, domain);
        const Raycast::Footprint footprint{settings.angular,
                                           settings.radial_start};
        const auto prepared = SDF::prepare_periodic_shells(
            camera, settings.cell, settings.radius, footprint);
        HS_EXPECT_FALSE(prepared.march);
        const Raycast::Appearance appearance{
            1 / settings.far, 0, 1 / settings.near_fade, &palette.view()};
        Raycast::TraceLimits limits;
        limits.max_steps = 256;
        SDF::ShellLayerStorage layers;
        for (int ray = 0; ray < 40; ++ray) {
          const float x = rand_uniform(random, -1.0f, 1.0f);
          const float y = rand_uniform(random, -1.0f, 1.0f);
          const float z = rand_uniform(random, -1.0f, 1.0f);
          const math::Vector direction = math::Vector{x, y, z}.normalized();
          const auto found = shell_neighborhood(
              camera, direction, settings.cell, settings.radius, footprint);
          const auto traced =
              domain == Raycast::SamplingDomain::SLICE_4D
                  ? SDF::trace_periodic_shells_march<4, true>(
                        prepared, camera, direction, limits, appearance, layers)
                  : SDF::trace_periodic_shells_march<3, true>(
                        prepared, camera, direction, limits, appearance,
                        layers);
          ++compared;
          hits += !found.empty();
          const auto expected =
              composite_shells(found, found.size(), appearance, limits);
          if (traced.status == Raycast::TraceStatus::BUDGET_EXHAUSTED) {
            truncated += !expected.second;
            bool prefix = false;
            for (size_t count = 0; count <= found.size() && !prefix; ++count)
              prefix = max_channel_gap(
                           traced.color,
                           composite_shells(found, count, appearance, limits)
                               .first) <= 2;
            mismatched += !prefix;
            continue;
          }
          HS_EXPECT_EQ(traced.status,
                       expected.second ? Raycast::TraceStatus::SATURATED
                                       : Raycast::TraceStatus::RANGE_COMPLETE);
          mismatched += max_channel_gap(traced.color, expected.first) > 24;
        }
      }
    }
    HS_EXPECT_GT(hits, compared / 2);
    HS_EXPECT_LE(truncated, compared / 20);
    HS_EXPECT_LE(mismatched, compared / 2000);
  }
}

/** @brief Sphere shading darkens color without weakening occlusion. */
inline void test_shell_shading_preserves_opacity() {
  PaletteFixture<ConstantPalette> fixture;
  auto &palette = fixture.palette;
  const Raycast::Appearance appearance{0, 0, 100, &palette.view()};
  for (const auto domain : {Raycast::SamplingDomain::SPATIAL_3D,
                            Raycast::SamplingDomain::SLICE_4D}) {
    for (const float OFFSET : {0.0f, .09f, .12f, .144f}) {
      Raycast::PreparedCamera camera;
      camera.domain = domain;
      camera.center = domain == Raycast::SamplingDomain::SPATIAL_3D
                          ? math::Vec4{{-.5f, OFFSET, 0, 0}}
                          : math::Vec4{{-.5f, 0, 0, OFFSET}};
      camera.interval = {0, 1};
      const float half_chord = sqrtf(.15f * .15f - OFFSET * OFFSET);
      const float BRIGHTNESS = 1 - .9f * OFFSET * OFFSET / (.15f * .15f);
      const Pixel expected = appearance.color(.5f - half_chord) * BRIGHTNESS;
      const auto prepared = SDF::prepare_periodic_shells(camera, 1, .15f, {});
      SDF::ShellLayerStorage storage;
      const auto narrow =
          domain == Raycast::SamplingDomain::SPATIAL_3D
              ? SDF::trace_periodic_shells_march<3, false, true>(
                    prepared, camera, math::X_AXIS, {}, appearance, storage)
              : SDF::trace_periodic_shells_march<4, false, true>(
                    prepared, camera, math::X_AXIS, {}, appearance, storage);
      const auto wide =
          domain == Raycast::SamplingDomain::SPATIAL_3D
              ? SDF::trace_periodic_shells_march<3, true, true>(
                    prepared, camera, math::X_AXIS, {}, appearance, storage)
              : SDF::trace_periodic_shells_march<4, true, true>(
                    prepared, camera, math::X_AXIS, {}, appearance, storage);
      HS_EXPECT_EQ(narrow.status, Raycast::TraceStatus::SATURATED);
      HS_EXPECT_EQ(wide.status, Raycast::TraceStatus::SATURATED);
      HS_EXPECT_LE(max_channel_gap(narrow.color, expected), 2);
      HS_EXPECT_LE(max_channel_gap(wide.color, expected), 2);
      if (domain == Raycast::SamplingDomain::SPATIAL_3D) {
        const auto owner = SDF::trace_periodic_shells_3d<true>(
            prepared, camera, math::X_AXIS, {}, appearance);
        HS_EXPECT_EQ(owner.status, Raycast::TraceStatus::SATURATED);
        HS_EXPECT_LE(max_channel_gap(owner.color, expected), 2);
      }
    }
  }
}

inline int run_ray_demonstrator_tests() {
  hs_test::ModuleFixture fixture("ray_demonstrators");
  HS_RUN(test_appearance_palette_interval);
  HS_RUN(test_octet_crossing_coverage_against_ray_line_distance);
  HS_RUN(test_periodic_shell_roots_and_slices);
  HS_RUN(test_shell_shading_preserves_opacity);
  HS_RUN(test_periodic_shell_traversal_budgets);
  HS_RUN(test_prepared_shells_match_sphere_roots);
  HS_RUN(test_shell_layer_march_matches_cell_traversal);
  HS_RUN(test_shell_slice_march_matches_cell_traversal);
  HS_RUN(test_shell_neighbor_march_matches_cell_traversal);
  HS_RUN(test_shell_march_uses_footprint_radial_start);
  HS_RUN(test_wide_shell_footprint_gathers_neighbor_rows);
  HS_RUN(test_wide_shell_march_gathers_layers_behind_start);
  HS_RUN(test_wide_shell_march_matches_neighborhood);
  HS_RUN(test_framework_generic_event_rendering);
  HS_RUN(test_repeated_stream_grouping_preserves_order_and_endpoints);
  HS_RUN(test_lattice_volume_camera_demonstrators);
  HS_RUN(test_torus_and_warped_volume_spherical_stage);
  HS_RUN(test_verified_filter_off_slice_geometry_and_projected_normal);
  return fixture.result();
}

} // namespace hs_test::ray_demonstrator_tests
