/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/sdf/framework.h"
#include "core/render/sdf/lattice_field.h"
#include "core/render/sdf/periodic_surface.h"
#include "core/render/ray/events.h"
#include "tests/test_harness.h"

namespace hs_test::sdf_pattern_tests {

inline void test_lattice_world_metric_and_fourth_axis() {
  SDF::WireLattice<3> cubic;
  HS_EXPECT_TRUE(cubic.valid());
  HS_EXPECT_NEAR(cubic.distance({0.25f, 0.1f, 0.0f}), 0.05f, 1e-6f);
  HS_EXPECT_NEAR(cubic.distance({0.25f, 1e-5f, 0.0f}),
                 1e-5f - cubic.wire_radius, 1e-7f);
  HS_EXPECT_NEAR(cubic.distance({0.0f, 0.0f, 0.0f}), -cubic.wire_radius, 1e-7f);
  cubic.cell_size = 2.0f;
  cubic.wire_radius = 0.1f;
  cubic.origin = {{3.0f, -2.0f, 1.0f, 0.0f}};
  HS_EXPECT_NEAR(cubic.distance({3.5f, -1.8f, 1.0f}), 0.1f, 1e-6f);
  math::rotate_plane(cubic.rotation, 0, 1, 0.7f);
  HS_EXPECT_TRUE(cubic.valid());
  const math::Vec4 ROTATED = cubic.rotation.apply({{0.5f, 0.2f, 0.0f, 0.0f}});
  const math::Vector P{ROTATED[0] + 3.0f, ROTATED[1] - 2.0f, ROTATED[2] + 1.0f};
  HS_EXPECT_NEAR(cubic.distance(P), 0.1f, 1e-6f);
  HS_EXPECT_NEAR(cubic.normal(P).magnitude(), 1.0f, 1e-6f);
  cubic.rotation.m[0][0] *= 2.0f;
  HS_EXPECT_FALSE(cubic.valid());

  SDF::WireLattice<4> hypercubic;
  HS_EXPECT_NEAR(hypercubic.distance({{0.4f, 0.0f, 0.0f, 0.0f}}), -0.05f,
                 1e-7f);
  HS_EXPECT_NEAR(hypercubic.distance({{0.4f, 0.0f, 0.0f, 0.3f}}), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(hypercubic.distance({{0.1f, 0.2f, 0.3f, 0.4f}}),
                 sqrtf(0.14f) - 0.05f, 1e-6f);
  HS_EXPECT_TRUE(hypercubic.capabilities().interior_clearance);
}

inline void test_lattice_reference_and_lipschitz() {
  SDF::WireLattice<4> lattice;
  for (int sample = 0; sample < 180; ++sample) {
    math::Vec4 p;
    for (int d = 0; d < 4; ++d)
      p[d] = sinf(static_cast<float>(sample * (d + 1)) * 0.73f) * 2.0f;
    float best = 100.0f;
    for (int axis = 0; axis < 4; ++axis) {
      float squared = 0.0f;
      for (int d = 0; d < 4; ++d) {
        if (d == axis)
          continue;
        float nearest = 100.0f;
        for (int plane = -3; plane <= 3; ++plane)
          nearest = std::min(nearest, fabsf(p[d] - plane));
        squared += nearest * nearest;
      }
      best = std::min(best, sqrtf(squared));
    }
    HS_EXPECT_NEAR(lattice.distance(p), best - lattice.wire_radius, 1e-6f);
    math::Vec4 q = p;
    q[3] += 0.03f;
    HS_EXPECT_LE(fabsf(lattice.distance(p) - lattice.distance(q)), 0.030001f);
  }
}

inline void test_framework_geometry_and_plane_streams() {
  SDF::TriangularFramework framework;
  HS_EXPECT_TRUE(framework.valid());
  const float H = SDF::TriangularFramework::TRIANGLE_HEIGHT;
  HS_EXPECT_NEAR(framework.distance({0.5f, H, 0.37f}), -0.05f, 1e-6f);
  HS_EXPECT_NEAR(framework.distance({0.25f, 0.5f * H, 0.0f}), -0.05f, 1e-6f);
  HS_EXPECT_GT(framework.distance({0.5f, H / 3.0f, 0.5f}), 0.2f);
  const math::Vector P{0.23f, 0.19f, 0.17f};
  HS_EXPECT_NEAR(framework.distance(P),
                 framework.distance(P + math::Vector{0.5f, H, 1.0f}), 1e-6f);
  const auto FAMILIES = framework.plane_families();
  HS_EXPECT_NEAR(math::dot(FAMILIES[0].normal, FAMILIES[1].normal), 0.5f,
                 1e-6f);
  Raycast::Ray ray{{0.2f, 0.1f, 0.3f}, {0.0f, 0.0f, -1.0f}, {0.0f, 4.0f}};
  SDF::FrameworkEvents events(framework, ray);
  HS_EXPECT_FALSE(events.active(0));
  HS_EXPECT_FALSE(events.active(1));
  HS_EXPECT_FALSE(events.active(2));
  HS_EXPECT_TRUE(events.active(3));
  HS_EXPECT_NEAR(events.distance(3), 0.3f, 1e-6f);
  HS_EXPECT_FALSE(events.candidate(3).verified);
  events.advance(3);
  HS_EXPECT_NEAR(events.distance(3), 1.3f, 1e-6f);
  framework.origin = {2.0f, -1.0f, 0.7f};
  HS_EXPECT_NEAR(framework.distance(P + framework.origin),
                 SDF::TriangularFramework{}.distance(P), 1e-6f);
  framework.cell_size = 0.0f;
  HS_EXPECT_FALSE(framework.valid());
}

inline float octet_segment_reference(const math::Vector &p) {
  constexpr float HALF_CUBE = 0.7071067811865475f;
  const std::array<math::Vector, 6> DIRECTIONS{
      {{1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1}}};
  float best = INFINITY;
  for (int x = -2; x <= 2; ++x)
    for (int y = -2; y <= 2; ++y)
      for (int z = -2; z <= 2; ++z) {
        if ((x + y + z) % 2 != 0)
          continue;
        const math::Vector VERTEX =
            math::Vector{static_cast<float>(x), static_cast<float>(y),
                         static_cast<float>(z)} *
            HALF_CUBE;
        for (const auto &direction : DIRECTIONS) {
          const math::Vector EDGE = direction * HALF_CUBE;
          const float ALONG =
              std::clamp(math::dot(p - VERTEX, EDGE), 0.0f, 1.0f);
          best = std::min(best, (p - VERTEX - EDGE * ALONG).magnitude());
        }
      }
  return best;
}

inline void test_octet_fcc_geometry_and_symmetry() {
  SDF::OctetFramework octet;
  HS_EXPECT_TRUE(octet.valid());
  const auto PLANES = octet.plane_families();
  for (size_t i = 0; i < PLANES.size(); ++i) {
    HS_EXPECT_NEAR(PLANES[i].normal.magnitude(), 1.0f, 1e-6f);
    HS_EXPECT_NEAR(PLANES[i].spacing, sqrtf(2.0f / 3.0f), 1e-6f);
    for (size_t j = i + 1; j < PLANES.size(); ++j)
      HS_EXPECT_NEAR(math::dot(PLANES[i].normal, PLANES[j].normal),
                     -1.0f / 3.0f, 1e-6f);
  }
  constexpr float H = 0.7071067811865475f;
  for (int zero = 0; zero < 3; ++zero)
    for (int a : {-1, 1})
      for (int b : {-1, 1}) {
        const math::Vector END = zero == 0   ? math::Vector{0, a * H, b * H}
                                 : zero == 1 ? math::Vector{a * H, 0, b * H}
                                             : math::Vector{a * H, b * H, 0};
        HS_EXPECT_NEAR(END.magnitude(), octet.cell_size, 1e-6f);
        for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
          HS_EXPECT_NEAR(octet.distance(END * fraction), -octet.wire_radius,
                         1e-6f);
      }
  HS_EXPECT_NEAR(octet.distance({H, 0, 0}), 0.5f - octet.wire_radius, 1e-6f);
  HS_EXPECT_NEAR(octet.distance({H / 2, H / 2, H / 2}),
                 H / 2 - octet.wire_radius, 1e-6f);
  for (int i = 0; i < 100; ++i) {
    const math::Vector P{0.7f * sinf(i * 0.73f), 0.7f * cosf(i * 0.47f),
                         0.7f * sinf(i * 0.31f)};
    const float DISTANCE = octet.distance(P);
    HS_EXPECT_NEAR(DISTANCE, octet_segment_reference(P) - octet.wire_radius,
                   2e-6f);
    HS_EXPECT_NEAR(DISTANCE, octet.distance({P.z, P.x, P.y}), 1e-6f);
    HS_EXPECT_NEAR(DISTANCE, octet.distance({-P.x, P.y, P.z}), 1e-6f);
    HS_EXPECT_NEAR(DISTANCE, octet.distance(P + math::Vector{H, H, 0}), 1e-6f);
    HS_EXPECT_NEAR(octet.normal(P).magnitude(), 1.0f, 1e-6f);
    const math::Vector Q = P + math::Vector{0.01f, -0.03f, 0.02f};
    HS_EXPECT_LE(fabsf(DISTANCE - octet.distance(Q)),
                 (P - Q).magnitude() + 1e-6f);
    SDF::OctetFramework scaled;
    scaled.cell_size *= 3;
    scaled.wire_radius *= 3;
    scaled.origin = {2, -3, 4};
    HS_EXPECT_NEAR(scaled.distance(P * 3 + scaled.origin), DISTANCE * 3, 3e-6f);
    for (size_t plane = 0; plane < PLANES.size(); ++plane) {
      const auto N = PLANES[plane].normal;
      const auto ON_PLANE = P - N * math::dot(P, N);
      const auto SAMPLE = octet.plane_sample(plane, ON_PLANE);
      HS_EXPECT_GE(SAMPLE.field + 2e-6f, octet.distance(ON_PLANE));
      HS_EXPECT_FALSE(SAMPLE.boundary);
    }
  }
  HS_EXPECT_FALSE(octet.sample({H, 0, 0}).boundary);
  octet.cell_size = 0;
  HS_EXPECT_FALSE(octet.valid());
  octet.cell_size = 1;
  octet.wire_radius = NAN;
  HS_EXPECT_FALSE(octet.valid());
}

inline void test_octet_events_ties_limits_and_invalid_inputs() {
  const SDF::OctetFramework OCTET;
  const Raycast::Ray RAY{{0, 0, 0}, {1, 0, 0}, {0, 1.5f}};
  SDF::OctetEvents events(OCTET, RAY);
  Raycast::TraceLimits limits;
  size_t count = 0;
  const auto RESULT = Raycast::trace_events(
      events, RAY.interval, limits, [&](const Raycast::Contribution &hit) {
        HS_EXPECT_NEAR(hit.t, static_cast<float>(count) * sqrtf(2.0f), 1e-6f);
        HS_EXPECT_EQ(hit.coverage, 1.0f);
        HS_EXPECT_FALSE(hit.verified);
        ++count;
        return true;
      });
  HS_EXPECT_EQ(RESULT.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_FALSE(RESULT.has_surface);
  HS_EXPECT_EQ(RESULT.counters.candidates, 8);
  HS_EXPECT_EQ(count, size_t{2});
  limits.max_candidates = 3;
  events = SDF::OctetEvents(OCTET, RAY);
  const auto LIMITED = Raycast::trace_events(events, RAY.interval, limits,
                                             [](const auto &) { return true; });
  HS_EXPECT_EQ(LIMITED.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(LIMITED.counters.candidates, 3);
  HS_EXPECT_EQ(LIMITED.counters.layers, 1);
  limits.max_candidates = 64;
  limits.max_layers = 1;
  events = SDF::OctetEvents(OCTET, RAY);
  const auto LAYERS = Raycast::trace_events(events, RAY.interval, limits,
                                            [](const auto &) { return true; });
  HS_EXPECT_EQ(LAYERS.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(LAYERS.counters.layers, 1);
  const Raycast::Ray EDGE_RAY{
      {0, 0, 0}, math::Vector{1, 1, 0}.normalized(), {0, 3}};
  SDF::OctetEvents edge_events(OCTET, EDGE_RAY);
  HS_EXPECT_FALSE(edge_events.active(1));
  HS_EXPECT_FALSE(edge_events.active(2));
  edge_events.advance(0);
  HS_EXPECT_NEAR(edge_events.distance(0), 1.0f, 1e-6f);
  auto invalid = RAY;
  invalid.direction.x = 0;
  SDF::OctetEvents bad_ray(OCTET, invalid);
  auto bad_geometry = OCTET;
  bad_geometry.cell_size = -1;
  SDF::OctetEvents bad_shape(bad_geometry, RAY);
  for (size_t i = 0; i < SDF::OctetEvents::STREAM_COUNT; ++i) {
    HS_EXPECT_FALSE(bad_ray.active(i));
    HS_EXPECT_FALSE(bad_shape.active(i));
  }
  const Raycast::Ray OBLIQUE{
      {0.2f, -0.3f, 0.4f}, math::Vector{-1, 2, -3}.normalized(), {0.3f, 4}};
  SDF::OctetEvents slanted(OCTET, OBLIQUE, {0.05f, 0.1f});
  float previous = OBLIQUE.interval.near;
  limits.max_layers = 32;
  const auto SLANTED = Raycast::trace_events(slanted, OBLIQUE.interval, limits,
                                             [&](const auto &hit) {
                                               HS_EXPECT_GE(hit.t, previous);
                                               HS_EXPECT_GE(hit.coverage, 0.0f);
                                               HS_EXPECT_LE(hit.coverage, 1.0f);
                                               HS_EXPECT_FALSE(hit.verified);
                                               previous = hit.t;
                                               return true;
                                             });
  HS_EXPECT_EQ(SLANTED.status, Raycast::TraceStatus::RANGE_COMPLETE);
}

template <typename Surface> void check_periodic_surface() {
  Surface surface;
  surface.period = 2.3f;
  surface.iso = 0.27f;
  HS_EXPECT_TRUE(surface.valid());
  HS_EXPECT_TRUE(surface.capabilities().interior_clearance);
  constexpr float DELTA = 0.0003f;
  for (int i = 0; i < 90; ++i) {
    const math::Vector P{0.073f * i, sinf(i * 0.7f), cosf(i * 0.4f)};
    const math::Vector G = surface.gradient(P);
    HS_EXPECT_LE(G.magnitude(), surface.lipschitz() + 1e-5f);
    HS_EXPECT_NEAR(G.x,
                   (surface.field(P + math::Vector{DELTA, 0.0f, 0.0f}) -
                    surface.field(P - math::Vector{DELTA, 0.0f, 0.0f})) /
                       (2.0f * DELTA),
                   0.012f);
    HS_EXPECT_NEAR(surface.field(P),
                   surface.field(P + math::Vector{surface.period, 0.0f, 0.0f}),
                   1e-5f);
    const math::Vector Q = P + math::Vector{0.11f, -0.03f, 0.07f};
    HS_EXPECT_LE(fabsf(surface.distance(P) - surface.distance(Q)),
                 (P - Q).magnitude() + 1e-6f);
    const auto SAMPLE = surface.sample(P);
    HS_EXPECT_NEAR(SAMPLE.clearance, fabsf(surface.distance(P)), 1e-7f);
    Surface scaled = surface;
    scaled.period *= 3.0f;
    HS_EXPECT_NEAR(scaled.distance(P * 3.0f), surface.distance(P) * 3.0f,
                   3e-6f);
    Surface placed = surface;
    placed.origin = {3.0f, -5.0f, 7.0f};
    HS_EXPECT_NEAR(placed.distance(P + placed.origin), surface.distance(P),
                   2e-6f);
  }
  surface.period = 0.0f;
  HS_EXPECT_FALSE(surface.valid());
}

inline void test_periodic_surface_bounds_and_gradients() {
  check_periodic_surface<SDF::CosineSurface>();
  check_periodic_surface<SDF::GyroidSurface>();
  SDF::CosineSurface cosine;
  HS_EXPECT_NEAR(cosine.field({0.0f, 0.0f, 0.0f}), 3.0f, 1e-6f);
  HS_EXPECT_NEAR(cosine.field({0.5f, 0.5f, 0.5f}), -3.0f, 1e-6f);
  HS_EXPECT_EQ(cosine.normal({0.0f, 0.0f, 0.0f}).magnitude(), 0.0f);
  SDF::GyroidSurface gyroid;
  HS_EXPECT_EQ(gyroid.field({0.0f, 0.0f, 0.0f}), 0.0f);
  HS_EXPECT_TRUE(gyroid.sample({0.0f, 0.0f, 0.0f}).boundary);
}

inline int run_sdf_pattern_tests() {
  const auto MODULE = hs_test::begin_module("sdf_patterns");
  test_lattice_world_metric_and_fourth_axis();
  test_lattice_reference_and_lipschitz();
  test_framework_geometry_and_plane_streams();
  test_octet_fcc_geometry_and_symmetry();
  test_octet_events_ties_limits_and_invalid_inputs();
  test_periodic_surface_bounds_and_gradients();
  return hs_test::end_module(MODULE);
}

} // namespace hs_test::sdf_pattern_tests
