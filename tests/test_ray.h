/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
    const float W = p[3] - offset;
    const float FIELD =
        sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + W * W) - 1.0f;
    return {FIELD, fabsf(FIELD), FIELD == 0.0f};
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

inline void test_camera() {
  Raycast::PreparedCamera camera;
  camera.center = {{4.0f, 5.0f, 6.0f, 0.0f}};
  camera.radial_start = 2.0f;
  camera.interval = {1.0f, 8.0f};
  HS_EXPECT_TRUE(camera.valid());
  const auto RAY = camera.ray(math::Vector(1, 0, 0));
  const auto POINT = camera.point3(RAY.at(3.0f));
  HS_EXPECT_NEAR(POINT.x, 9.0f, 1e-6f);
  HS_EXPECT_NEAR(POINT.y, 5.0f, 1e-6f);
  HS_EXPECT_NEAR((Raycast::Footprint{0.1f, 2.0f}.at(3.0f)), 0.5f, 1e-6f);
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  math::rotate_plane(camera.embedding, 0, 3, 0.6f);
  HS_EXPECT_TRUE(camera.valid());
  const auto DIRECTION = camera.point4(math::Vector(1, 0, 0));
  HS_EXPECT_NEAR(DIRECTION[3], sinf(0.6f), 1e-6f);
  math::Vector normal;
  HS_EXPECT_TRUE(camera.project_normal({{1, 0, 0, 0}}, normal));
  HS_EXPECT_NEAR(normal.x, 1.0f, 1e-6f);
  camera.embedding = math::Mat4::identity();
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

inline void test_surface_boundaries() {
  const Sphere SPHERE;
  const Raycast::VolumeQuery QUERY{
      SPHERE, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::Ray ENTRY{
      math::Vector(-3, 0, 0), math::Vector(1, 0, 0), {0, 6}};
  auto result = Raycast::surface_search(QUERY, ENTRY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 2.0f, 1e-4f);
  HS_EXPECT_TRUE(result.contribution.verified);
  HS_EXPECT_EQ(result.counters.layers, 1);
  auto exit = ENTRY;
  exit.origin = math::Vector(0, 0, 0);
  result = Raycast::surface_search(QUERY, exit, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.0f, 1e-4f);
  const Raycast::VolumeQuery EXTERIOR_ONLY{SPHERE,
                                           Raycast::QueryCapabilities{}};
  result = Raycast::surface_search(EXTERIOR_ONLY, exit, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::UNSUPPORTED_START);
  exit.origin = math::Vector(1, 0, 0);
  result = Raycast::surface_search(QUERY, exit, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_EQ(result.contribution.t, 0.0f);
  HS_EXPECT_EQ(result.counters.queries, 1);
  auto miss = ENTRY;
  miss.origin = math::Vector(-3, 2, 0);
  result = Raycast::surface_search(QUERY, miss, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_TRUE(!result.has_surface);
}

inline void test_bounded_failures() {
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 6}};
  auto result = Raycast::surface_search(ConstantQuery{}, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::UNRESOLVED);
  HS_EXPECT_TRUE(!result.has_surface);
  HS_EXPECT_EQ(result.counters.queries, 2);
  Raycast::TraceLimits limits;
  limits.max_queries = 1;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.queries, 1);
  limits.max_queries = 96;
  limits.max_steps = 1;
  result = Raycast::surface_search(ConstantQuery{1, 0.1f}, RAY, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.steps, 1);
  result = Raycast::surface_search(ConstantQuery{NAN, 1}, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
  result = Raycast::surface_search(ConstantQuery{1, -1}, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
  limits.position_tolerance = 0.0f;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
}

inline void test_slice_no_phantom() {
  Raycast::PreparedCamera camera;
  camera.domain = Raycast::SamplingDomain::SLICE_4D;
  const SliceBall BALL;
  const Raycast::DomainQuery4<SliceBall> QUERY{BALL, camera};
  const Raycast::Ray RAY{math::Vector(-2, 0, 0), math::Vector(1, 0, 0), {0, 4}};
  Raycast::TraceLimits limits;
  limits.max_steps = 512;
  limits.max_queries = 1024;
  const auto RESULT = Raycast::surface_search(QUERY, RAY, {0.1f, 0}, limits);
  HS_EXPECT_TRUE(!RESULT.has_surface);
  HS_EXPECT_TRUE(RESULT.status == Raycast::TraceStatus::RANGE_COMPLETE);
  const SliceBall CUT_BALL{0.5f};
  const Raycast::DomainQuery4<SliceBall> CUT_QUERY{CUT_BALL, camera};
  const auto CUT = Raycast::surface_search(CUT_QUERY, RAY, {}, limits);
  HS_EXPECT_TRUE(CUT.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(CUT.contribution.t, 2.0f - sqrtf(0.75f), 1e-4f);
}

inline void test_placement_and_shapes() {
  const Sphere SPHERE;
  const Raycast::VolumeQuery QUERY{
      SPHERE, Raycast::QueryCapabilities{true, true, true}};
  Raycast::PlacedQuery PLACED{QUERY, math::Vector(5, 0, 0), math::Quaternion(),
                              2.0f};
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 10}};
  auto result = Raycast::surface_search(PLACED, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 3.0f, 1e-4f);
  HS_EXPECT_NEAR(PLACED.sample(math::Vector()).field, 3.0f, 1e-6f);
  HS_EXPECT_NEAR(PLACED.sample(math::Vector(5, 0, 0)).field, -2.0f, 1e-6f);
  PLACED.scale = 0.0f;
  HS_EXPECT_TRUE(!PLACED.valid());
  result = Raycast::surface_search(PLACED, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
  PLACED.scale = 1.0f;
  PLACED.inverse_rotation = math::Quaternion(0, 0, 0, 0);
  HS_EXPECT_TRUE(!PLACED.valid());
  const SDF::Torus TORUS{2.0f, 0.25f};
  const Raycast::VolumeQuery TORUS_QUERY{
      TORUS, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::PreparedCamera CAMERA;
  const Raycast::DomainQuery3<decltype(TORUS_QUERY)> DOMAIN_QUERY{TORUS_QUERY,
                                                                  CAMERA};
  result = Raycast::surface_search(DOMAIN_QUERY, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.75f, 1e-4f);
  const SDF::WarpedVolume WARPED{TORUS, SDF::Warp::Twist{3, 0.1f, 2.0f}, 0.01f};
  const Raycast::VolumeQuery WARP_QUERY{WARPED, Raycast::QueryCapabilities{}};
  const Raycast::DomainQuery3<decltype(WARP_QUERY)> WARP_DOMAIN{WARP_QUERY,
                                                                CAMERA};
  result = Raycast::surface_search(WARP_DOMAIN, RAY, {});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.75f, 1e-4f);
}

inline void test_first_boundary_and_tolerances() {
  const ThinNeighbors NEIGHBORS;
  const Raycast::VolumeQuery QUERY{
      NEIGHBORS, Raycast::QueryCapabilities{true, true, true}};
  Raycast::Ray ray{math::Vector(), math::Vector(1, 0, 0), {0, 2}};
  Raycast::TraceLimits limits;
  limits.position_tolerance = 1e-5f;
  auto result = Raycast::surface_search(QUERY, ray, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 0.9998f, 1e-5f);
  ray.interval.far = 0.9f;
  result = Raycast::surface_search(QUERY, ray, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::RANGE_COMPLETE);
  ray.interval = {1.0005f, 2.0f};
  result = Raycast::surface_search(QUERY, ray, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 1.0008f, 1e-5f);
  const Sphere SPHERE;
  Raycast::VolumeQuery sphere_query{
      SPHERE, Raycast::QueryCapabilities{true, true, true}};
  const Raycast::PlacedQuery TINY{sphere_query, math::Vector(0.0005f, 0, 0),
                                  math::Quaternion(), 0.0001f};
  ray.interval = {0, 0.001f};
  limits.position_tolerance = 1e-8f;
  result = Raycast::surface_search(TINY, ray, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::SURFACE);
  HS_EXPECT_NEAR(result.contribution.t, 0.0004f, 1e-8f);
  sphere_query.guarantees.error = 0.001f;
  result = Raycast::surface_search(sphere_query, ray, {}, limits);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
  sphere_query.guarantees.error = 0;
  ray = {math::Vector(-2, 1, 0), math::Vector(1, 0, 0), {0, 4}};
  limits = {};
  result = Raycast::surface_search(sphere_query, ray, {}, limits);
  HS_EXPECT_TRUE(result.status != Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
}

inline void test_limits_and_nonfinite() {
  HS_EXPECT_TRUE(!Raycast::finite(NAN));
  HS_EXPECT_TRUE(!Raycast::finite(INFINITY));
  HS_EXPECT_TRUE(Raycast::finite(0.0f));
  const Raycast::Ray RAY{math::Vector(), math::Vector(1, 0, 0), {0, 6}};
  Raycast::TraceLimits limits;
  limits.max_layers = 0;
  auto result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.counters.queries, 0);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED);
  limits.max_layers = 1;
  limits.max_refinements = 0;
  result = Raycast::surface_search(ConstantQuery{}, RAY, {}, limits);
  HS_EXPECT_EQ(result.counters.refinements, 0);
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::UNRESOLVED);
  result = Raycast::surface_search(ConstantQuery{}, RAY, {-1, 0});
  HS_EXPECT_TRUE(result.status == Raycast::TraceStatus::INVALID_QUERY);
}

inline void test_nested_query_validation_once() {
  struct CountedQuery : ConstantQuery {
    mutable int validations = 0;
    bool valid() const {
      ++validations;
      return true;
    }
  };
  const CountedQuery query{{1.0f, 0.1f}};
  const Raycast::PreparedCamera camera;
  const Raycast::DomainQuery3<CountedQuery> domain{query, camera};
  const Raycast::PlacedQuery placed{domain, math::Vector(), math::Quaternion(),
                                    1.0f};
  const Raycast::Ray ray{math::Vector(), math::Vector(1, 0, 0), {0, 4}};
  const auto result = Raycast::surface_search(placed, ray, {});
  HS_EXPECT_GT(result.counters.queries, 1);
  HS_EXPECT_EQ(query.validations, 1);
}

inline int run_ray_tests() {
  hs_test::ModuleFixture fixture("ray");
  test_camera();
  test_surface_boundaries();
  test_bounded_failures();
  test_slice_no_phantom();
  test_placement_and_shapes();
  test_nested_query_validation_once();
  test_limits_and_nonfinite();
  test_first_boundary_and_tolerances();
  return fixture.result();
}

} // namespace hs_test::ray_tests
