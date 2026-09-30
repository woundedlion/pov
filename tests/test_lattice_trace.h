/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/sdf/lattice.h"
#include "render/sdf/lattice_trace.h"
#include "tests/test_harness.h"

namespace hs_test::lattice_trace_tests {

struct DepthPalette {
  Color4 get(float t) const {
    return {Pixel(uint16_t(1000 + 50000 * t), uint16_t(60000 - 45000 * t),
                  uint16_t(12000 + 30000 * t)),
            1.0f};
  }
};

inline math::Vector direction(int sample) {
  if (sample == 0)
    return math::X_AXIS;
  if (sample == 1)
    return math::Y_AXIS;
  if (sample == 2)
    return math::Z_AXIS;
  const float z = (sample - 2.5f) / 24.0f - 1.0f;
  const float radius = sqrtf(1.0f - z * z);
  const float phase = sample * 2.39996323f;
  return {radius * cosf(phase), radius * sinf(phase), z};
}

template <bool SLICE, uint8_t SHELLS>
void check_cubic(const SDF::Lattice::PreparedShading &prepared,
                 const math::Vector &view, int &lit) {
  const auto composite =
      SDF::Lattice::composite_crossings<SLICE, SHELLS>(view, prepared);
  const auto actual = composite.finish();
  SDF::Lattice::Events<SLICE, SHELLS> events(view, prepared.lattice);
  const auto expected = Raycast::shade_events(
      events, {0, prepared.lattice.far_distance}, {}, prepared.appearance);
  HS_EXPECT_NEAR(actual.alpha, expected.color.alpha, 3e-4f);
  const Pixel color = composite.premultiplied();
  const Pixel expected_color = expected.color.color * expected.color.alpha;
  HS_EXPECT_NEAR(color.r, expected_color.r, 2);
  HS_EXPECT_NEAR(color.g, expected_color.g, 2);
  HS_EXPECT_NEAR(color.b, expected_color.b, 2);
  lit += actual.alpha > 0.0f;
}

inline void test_cubic_compositor_matches_event_shading() {
  alignas(std::max_align_t) uint8_t storage[4096];
  Arena arena(storage, sizeof(storage));
  BakedPaletteStorage palette;
  palette.bake(arena, DepthPalette{});
  SDF::Lattice::CrossingList crossings;
  int lit = 0;
  for (const auto domain :
       {SDF::Lattice::Domain::THREE_D, SDF::Lattice::Domain::FOUR_D_SLICE}) {
    for (const float cell : {.4f, 1.7f}) {
      SDF::Lattice::Settings settings;
      settings.mode = domain;
      settings.cell_size = cell;
      settings.wire_radius = .12f;
      auto embedding = math::Mat4::identity();
      math::rotate_plane(embedding, 0, 1, .31f);
      if (domain == SDF::Lattice::Domain::FOUR_D_SLICE)
        math::rotate_plane(embedding, 0, 3, .73f);
      for (const auto shells :
           {SDF::Lattice::ShellCount::ONE, SDF::Lattice::ShellCount::TWO,
            SDF::Lattice::ShellCount::THREE}) {
        settings.shells = shells;
        const SDF::Lattice::PreparedShading prepared{
            SDF::Lattice::prepare(settings, {{.231f, -.437f, .719f, .383f}},
                                  embedding, 7.0f, .018f),
            {1.0f / 7.0f, 0.0f, 2.0f, &palette.view(), .8f},
            &crossings};
        for (int i = 0; i < 51; ++i) {
          const auto view = direction(i);
          check_cubic<false, 0>(prepared, view, lit);
          check_cubic<true, 1>(prepared, view, lit);
          check_cubic<true, 2>(prepared, view, lit);
          check_cubic<true, 3>(prepared, view, lit);
        }
      }
    }
  }
  HS_EXPECT_GT(lit, 100);
}

inline void test_octet_preparation_matches_world_events() {
  namespace Trace = SDF::LatticeTrace;
  alignas(std::max_align_t) uint8_t storage[4096];
  Arena arena(storage, sizeof(storage));
  BakedPaletteStorage palette;
  palette.bake(arena, DepthPalette{});
  Trace::CrossingStorage crossings;
  Trace::Settings settings;
  settings.palette = &palette.view();
  settings.crossings = &crossings;
  settings.pixel_half_angle = .018f;
  int lit = 0;
  for (const auto domain : {Raycast::SamplingDomain::SPATIAL_3D,
                            Raycast::SamplingDomain::SLICE_4D}) {
    settings.domain = domain;
    const bool SLICE = domain == Raycast::SamplingDomain::SLICE_4D;
    settings.center = {{.231f, -.437f, .719f, SLICE ? .383f : 0.0f}};
    settings.embedding = math::Mat4::identity();
    math::rotate_plane(settings.embedding, 0, 1, .31f);
    if (SLICE)
      math::rotate_plane(settings.embedding, 0, 3, .73f);
    for (float cell : {.4f, 1.7f}) {
      settings.cell_size = cell;
      settings.wire_radius = .055f * cell;
      for (float radial : {0.0f, .7f}) {
        settings.radial_start = radial;
        const auto prepared = Trace::prepare(settings);
        HS_EXPECT_TRUE(prepared.valid);
        for (int i = 0; i < 51; ++i) {
          const auto view = direction(i);
          const auto ray = prepared.camera.ray(view);
          const auto ambient =
              prepared.camera.embedding.apply({{view.x, view.y, view.z, 0.0f}});
          Raycast::ShadedTrace expected;
          if (SLICE) {
            SDF::OctetEvents4 events(prepared.octet4,
                                     prepared.camera.point4(ray.origin),
                                     ambient, ray.interval, prepared.footprint);
            expected = Raycast::shade_events(
                events, ray.interval, prepared.limits, prepared.appearance);
          } else {
            const Raycast::Ray world{prepared.camera.point3(ray.origin),
                                     {ambient[0], ambient[1], ambient[2]},
                                     ray.interval};
            SDF::OctetEvents events(prepared.octet, world, prepared.footprint);
            expected = Raycast::shade_events(
                events, ray.interval, prepared.limits, prepared.appearance);
          }
          const auto actual = SLICE ? Trace::shade<true>(view, prepared)
                                    : Trace::shade<false>(view, prepared);
          const Pixel premultiplied =
              expected.color.color * expected.color.alpha;
          HS_EXPECT_EQ(actual.status, expected.trace.status);
          HS_EXPECT_NEAR(actual.color.r, premultiplied.r, 2);
          HS_EXPECT_NEAR(actual.color.g, premultiplied.g, 2);
          HS_EXPECT_NEAR(actual.color.b, premultiplied.b, 2);
          lit += actual.color != Pixel{};
        }
      }
    }
  }
  HS_EXPECT_GT(lit, 10);
}

inline void run_lattice_trace_cases() {
  test_cubic_compositor_matches_event_shading();
  test_octet_preparation_matches_world_events();
}

} // namespace hs_test::lattice_trace_tests
