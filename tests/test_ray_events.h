/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Ray events test-module entry. */

#include "tests/ray_events/lattice_trace.h"
#include "tests/ray_events/base_cases.h"

namespace hs_test::ray_event_tests {

inline int run_ray_event_tests() {
  hs_test::ModuleFixture fixture("ray_events");
  HS_RUN(test_octet_prepared_projection);
  HS_RUN(test_octet_4d_canonical_trace);
  HS_RUN(test_octet_4d_canonical_boundary_directions);
  HS_RUN(test_octet_trace_sorted_matches_events);
  HS_RUN(test_lattice_view_directions_are_unit_length);
  HS_RUN(test_lattice_prepare_validates_wire_radius_per_geometry);
  HS_RUN(test_cubic_compositor_matches_event_shading);
  HS_RUN(test_octet_preparation_matches_world_events);
  HS_RUN(test_octet_distant_center_periodicity);
  HS_RUN(test_single_group_capacity);
  HS_RUN(test_failure_status_survives_flush);
  HS_RUN(test_event_stream_contracts);
  HS_RUN(test_verified_filter_contracts);
  return fixture.result();
}
} // namespace hs_test::ray_event_tests
