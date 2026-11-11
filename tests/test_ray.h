/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Ray tracing convergence, verification, query bounds, and status contracts.
 */
#pragma once

#include "core/render/ray/query.h"
#include "core/render/ray/march.h"
#include "core/render/sdf/volume.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::ray_tests {

struct Sphere {
  float distance(const math::Vector &p) const {
    return sqrtf(math::dot(p, p)) - 1.0f;
  }
};

struct SliceBall {
  float offset = 1.001f;
  Raycast::QueryCapabilities capabilities() const { return {true, true, true}; }
  Raycast::QuerySample sample(const math::Vec4 &p) const {
    const float w = p[3] - offset;
    const float field =
        sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + w * w) - 1.0f;
    return {field, fabsf(field), field == 0.0f};
  }
};

struct ConstantQuery {
  float field = 1.0f;
  float clearance = 0.0f;
  Raycast::QueryCapabilities capabilities() const { return {true, true, true}; }
  Raycast::QuerySample sample(const math::Vector &) const {
    return {field, clearance, false};
  }
};

struct ThinNeighbors {
  float distance(const math::Vector &p) const {
    return std::min(fabsf(p.x - 1.0f) - 0.0002f, fabsf(p.x - 1.001f) - 0.0002f);
  }
};

/** @brief Camera transforms and normal projection obey the embedding contract and reject invalid state. */
inline void test_camera() {
  Raycast::PreparedCamera camera;
  camera.center = {{4.0f, 5.0f, 6.0f, 0.0f}};
  camera.radial_start = 2.0f;
  camera.interval = {1.0f, 8.0f};
  HS_EXPECT_TRUE(camera.valid());
  auto invalid_spatial = camera;
  invalid_spatial.center[3] = 1.0f;
  HS_EXPECT_FALSE(invalid_spatial.valid());
  invalid_spatial = camera;
  math::rotate_plane(invalid_spatial.embedding, 0, 3, 0.6f);
  HS_EXPECT_FALSE(invalid_spatial.valid());
  const auto ray = camera.ray(math::Vector(1, 0, 0));
  const auto point = camera.point3(ray.at(3.0f));
  HS_EXPECT_NEAR(point.x, 9.0f, 1e-6f);
  HS_EXPECT_NEAR(point.y, 5.0f, 1e-6f);
  HS_EXPECT_NEAR((Raycast::Footprint{0.1f, 2.0f}.at(3.0f)), 0.5f, 1e-6f);
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  math::rotate_plane(camera.embedding, 0, 3, 0.6f);
  HS_EXPECT_TRUE(camera.valid());
  const auto direction = camera.point4(math::Vector(1, 0, 0));
  HS_EXPECT_NEAR(direction[3], sinf(0.6f), 1e-6f);
  math::Vector normal;
  HS_EXPECT_TRUE(camera.project_normal({{1, 0, 0, 0}}, normal));
  HS_EXPECT_NEAR(normal.x, 1.0f, 1e-6f);
  camera.embedding = math::Mat4::identity();
  float inverse_length = 9.0f;
  HS_EXPECT_TRUE(camera.project_normal_dimension<3>({{3, 4, 0, 0}}, normal,
                                                    inverse_length));
  HS_EXPECT_NEAR(inverse_length, 0.2f, 1e-6f);
  HS_EXPECT_NEAR(normal.x, 0.6f, 1e-6f);
  HS_EXPECT_NEAR(normal.y, 0.8f, 1e-6f);
  HS_EXPECT_FALSE(camera.project_normal_dimension<4>({{0, 0, 0, 1}}, normal,
                                                     inverse_length));
  HS_EXPECT_EQ(inverse_length, 0.0f);
  HS_EXPECT_EQ(normal.magnitude(), 0.0f);
  HS_EXPECT_TRUE(!camera.project_normal({{0, 0, 0, 1}}, normal));
  camera.embedding.m[0][0] = 2.0f;
  HS_EXPECT_TRUE(!camera.valid());
  camera.embedding = math::Mat4::identity();
  camera.radial_start = -1.0f;
  HS_EXPECT_TRUE(!camera.valid());
  camera.radial_start = 0.0f;
  camera.interval = {1, 1};
  HS_EXPECT_TRUE(!camera.valid());
}

/** @brief Sphere searches resolve entry, exit and on-surface starts while rejecting unsupported interior starts. */
inline void test_surface_boundaries() {
  const Sphere sphere;
  const Raycast::VolumeQuery query{
      sphere, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::Ray ENTRY{
      math::Vector(-3, 0, 0), math::Vector(1, 0, 0), {0, 6}};
  auto result = Raycast::surface_search(query, ENTRY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 2.0f, 1e-4f);
  HS_EXPECT_TRUE(result.contribution.verified);
  HS_EXPECT_EQ(result.counters.layers, 1);
  auto exit = ENTRY;
  exit.origin = math::Vector(0, 0, 0);
  result = Raycast::surface_search(query, exit, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.0f, 1e-4f);
  const Raycast::VolumeQuery exterior_only{sphere,
                                           Raycast::QueryCapabilities{}};
  result = Raycast::surface_search(exterior_only, exit, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::UNSUPPORTED_START);
  exit.origin = math::Vector(1, 0, 0);
  result = Raycast::surface_search(query, exit, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_EQ(result.contribution.t, 0.0f);
  HS_EXPECT_EQ(result.counters.queries, 1);
  auto miss = ENTRY;
  miss.origin = math::Vector(-3, 2, 0);
  result = Raycast::surface_search(query, miss, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_TRUE(!result.has_surface);
}

/** @brief Search failures preserve query and step budgets and reject invalid samples. */
inline void test_bounded_failures() {
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 6}};
  auto result = Raycast::surface_search(ConstantQuery{}, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::UNRESOLVED);
  HS_EXPECT_TRUE(!result.has_surface);
  HS_EXPECT_EQ(result.counters.queries,
               Raycast::TraceLimits{}.max_refinements + 1);
  HS_EXPECT_EQ(result.counters.steps,
               Raycast::TraceLimits{}.max_refinements + 1);
  Raycast::TraceLimits limits;
  limits.max_queries = 1;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.queries, 1);
  limits.max_queries = 96;
  limits.max_steps = 1;
  result = Raycast::surface_search(ConstantQuery{1, 0.1f}, RAY, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.steps, 1);
  result = Raycast::surface_search(ConstantQuery{NAN, 1}, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  result = Raycast::surface_search(ConstantQuery{1, -1}, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  limits.position_tolerance = 0.0f;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
}

/** @brief An off-slice ball produces no surface while an intersecting slice matches its analytic root. */
inline void test_slice_no_phantom() {
  Raycast::PreparedCamera camera;
  const SliceBall ball;
  const Raycast::DomainQuery4<SliceBall> invalid_query{ball, camera};
  const Raycast::Ray SPATIAL_RAY{
      math::Vector(-2, 0, 0), math::Vector(1, 0, 0), {0, 4}};
  HS_EXPECT_EQ(Raycast::surface_search(invalid_query, SPATIAL_RAY, {}).status,
               Raycast::TraceStatus::INVALID_QUERY);
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  const Raycast::DomainQuery4<SliceBall> query{ball, camera};
  const Raycast::Ray RAY{math::Vector(-2, 0, 0), math::Vector(1, 0, 0), {0, 4}};
  Raycast::TraceLimits limits;
  limits.max_steps = 512;
  limits.max_queries = 1024;
  const auto result = Raycast::surface_search(query, RAY, {0.1f, 0}, limits);
  HS_EXPECT_TRUE(!result.has_surface);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  const SliceBall CUT_BALL{0.5f};
  const Raycast::DomainQuery4<SliceBall> cut_query{CUT_BALL, camera};
  const auto cut = Raycast::surface_search(cut_query, RAY, {}, limits);
  HS_EXPECT_EQ(cut.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(cut.contribution.t, 2.0f - sqrtf(0.75f), 1e-4f);
}

/** @brief Placed and warped shapes resolve their first surface and reject invalid transforms. */
inline void test_placement_and_shapes() {
  const Sphere sphere;
  const Raycast::VolumeQuery query{
      sphere, Raycast::QueryCapabilities{true, true, true}};
  Raycast::PlacedQuery placed{query, math::Vector(5, 0, 0), math::Quaternion(),
                              2.0f};
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 10}};
  auto result = Raycast::surface_search(placed, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 3.0f, 1e-4f);
  HS_EXPECT_NEAR(placed.sample(math::Vector()).field, 3.0f, 1e-6f);
  HS_EXPECT_NEAR(placed.sample(math::Vector(5, 0, 0)).field, -2.0f, 1e-6f);
  placed.scale = 0.0f;
  HS_EXPECT_TRUE(!placed.valid());
  result = Raycast::surface_search(placed, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  placed.scale = 1.0f;
  placed.inverse_rotation = math::Quaternion(0, 0, 0, 0);
  HS_EXPECT_TRUE(!placed.valid());
  const SDF::Torus TORUS{2.0f, 0.25f};
  const Raycast::VolumeQuery TORUS_QUERY{
      TORUS, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::PreparedCamera camera;
  const Raycast::DomainQuery3<decltype(TORUS_QUERY)> domain_query{TORUS_QUERY,
                                                                  camera};
  result = Raycast::surface_search(domain_query, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.75f, 1e-4f);
  const SDF::WarpedVolume warped{TORUS, SDF::Warp::Twist{3, 0.1f, 2.0f}, 0.01f};
  const Raycast::VolumeQuery warp_query{warped, Raycast::QueryCapabilities{}};
  const Raycast::DomainQuery3<decltype(warp_query)> warp_domain{warp_query,
                                                                camera};
  result = Raycast::surface_search(warp_domain, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.75f, 1e-4f);
  const Raycast::Ray displaced_ray{
      math::Vector(), math::Vector(sqrtf(3.0f) * 0.5f, 0.0f, 0.5f), {0, 10}};
  result = Raycast::surface_search(warp_domain, displaced_ray, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  const float hit = result.contribution.t;
  HS_EXPECT_TRUE(fabsf(hit - 1.75f) > 1e-3f);
  HS_EXPECT_TRUE(fabsf(warped.raw_distance(displaced_ray.at(hit))) <= 1e-3f);
  for (int k = 0; k < 64; ++k)
    HS_EXPECT_TRUE(warped.distance(displaced_ray.at(hit * k / 64.0f)) > 0.0f);
}

/** @brief Search resolves the first in-range thin boundary within the configured position tolerance. */
inline void test_first_boundary_and_tolerances() {
  const ThinNeighbors neighbors;
  const Raycast::VolumeQuery query{
      neighbors, Raycast::QueryCapabilities{true, true, true}};
  Raycast::Ray ray{math::Vector(), math::Vector(1, 0, 0), {0, 2}};
  Raycast::TraceLimits limits;
  limits.position_tolerance = 1e-5f;
  auto result = Raycast::surface_search(query, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 0.9998f, 1e-5f);
  ray.interval.far = 0.9f;
  result = Raycast::surface_search(query, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  ray.interval = {1.0005f, 2.0f};
  result = Raycast::surface_search(query, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.0008f, 1e-5f);
  const Sphere sphere;
  Raycast::VolumeQuery sphere_query{
      sphere, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::PlacedQuery tiny{sphere_query, math::Vector(0.0005f, 0, 0),
                                  math::Quaternion(), 0.0001f};
  ray.interval = {0, 0.001f};
  limits.position_tolerance = 1e-8f;
  result = Raycast::surface_search(tiny, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 0.0004f, 1e-8f);
  sphere_query.guarantees.error = 0.001f;
  result = Raycast::surface_search(sphere_query, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  sphere_query.guarantees.error = 0;
  ray = {math::Vector(-2, 1, 0), math::Vector(1, 0, 0), {0, 4}};
  limits = {};
  result = Raycast::surface_search(sphere_query, ray, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
}

/** @brief Zero layer and refinement budgets stop work and invalid footprints are rejected. */
inline void test_limits_and_nonfinite() {
  HS_EXPECT_TRUE(!Raycast::finite(NAN));
  HS_EXPECT_TRUE(!Raycast::finite(INFINITY));
  HS_EXPECT_TRUE(Raycast::finite(0.0f));
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 6}};
  Raycast::TraceLimits limits;
  limits.max_layers = 0;
  auto result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.counters.queries, 0);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  limits.max_layers = 1;
  limits.max_refinements = 0;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.counters.refinements, 0);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::UNRESOLVED);
  result = Raycast::surface_search(ConstantQuery{}, RAY, {-1, 0});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
}

/** @brief Nested query wrappers validate the underlying query once across repeated samples. */
inline void test_nested_query_validation_once() {
  struct CountedQuery : ConstantQuery {
    mutable int validations = 0;
    bool valid() const {
      ++validations;
      return true;
    }
  };
  const CountedQuery QUERY{{1.0f, 0.1f}};
  const Raycast::PreparedCamera camera;
  const Raycast::DomainQuery3<CountedQuery> domain{QUERY, camera};
  const Raycast::PlacedQuery placed{domain, math::Vector(), math::Quaternion(),
                                    1.0f};
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 4}};
  const auto result = Raycast::surface_search(placed, RAY, {});
  HS_EXPECT_GT(result.counters.queries, 1);
  HS_EXPECT_EQ(QUERY.validations, 1);
}

/**
 * @brief Clearance within query error reaches the interval end through
 * bounded refinement.
 */
inline void test_uncertain_clearance_progress() {
  struct UncertainQuery : ConstantQuery {
    Raycast::QueryCapabilities capabilities() const {
      return {true, true, true, 0.0001f};
    }
  };
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 0.001f}};
  Raycast::TraceLimits limits;
  limits.position_tolerance = 0.0002f;
  const auto result = Raycast::surface_search(UncertainQuery{{1.0f, 0.00005f}},
                                              RAY, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_TRUE(!result.has_surface);
  HS_EXPECT_GT(result.counters.refinements, 0);
  HS_EXPECT_LT(result.counters.refinements, limits.max_refinements);
}

/** @brief Invalid clearance and unresolved zero-width signs never produce a surface. */
inline void test_unverified_boundaries_remain_unresolved() {
  struct OverstatedClearance {
    Raycast::QueryCapabilities capabilities() const {
      return {true, true, true};
    }
    Raycast::QuerySample sample(const math::Vector &p) const {
      return {p.x < 0.5f ? 1.0f : -1.0f, 1.0f, false};
    }
  };
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 6}};
  auto result = Raycast::surface_search(OverstatedClearance{}, RAY, {});
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  HS_EXPECT_FALSE(result.has_surface);
  Raycast::TraceLimits limits;
  limits.position_tolerance = 0.0002f;
  const Raycast::Ray SHORT_RAY{
      math::Vector(), math::Vector(1, 0, 0), {0, 0.001f}};
  result = Raycast::surface_search(ConstantQuery{0, 0}, SHORT_RAY, {}, limits);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::UNRESOLVED);
  HS_EXPECT_FALSE(result.has_surface);
}

inline int run_ray_tests() {
  hs_test::ModuleFixture fixture("ray");
  HS_RUN(test_camera);
  HS_RUN(test_unverified_boundaries_remain_unresolved);
  HS_RUN(test_surface_boundaries);
  HS_RUN(test_bounded_failures);
  HS_RUN(test_slice_no_phantom);
  HS_RUN(test_placement_and_shapes);
  HS_RUN(test_nested_query_validation_once);
  HS_RUN(test_limits_and_nonfinite);
  HS_RUN(test_first_boundary_and_tolerances);
  HS_RUN(test_uncertain_clearance_progress);
  return fixture.result();
}

} // namespace hs_test::ray_tests
