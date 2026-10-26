/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_ray.h"
#include "tests/test_ray_demonstrators.h"
#include "tests/test_sdf_patterns.h"
#include "tests/test_ray_events.h"
#include "tests/test_hyper_lattice.h"
#include "tests/test_pullback.h"
#include "tests/test_shader_chain.h"
#include "tests/test_filter.h"
#include "tests/test_plot_parity.h"
#include "tests/test_plot_scan.h"
#include "tests/test_canvas.h"
#include "tests/test_scan.h"
#include "tests/test_mesh_raster.h"

namespace hs_test::module_runners {

int run_ray_tests() { return hs_test::ray_tests::run_ray_tests(); }
int run_ray_demonstrator_tests() {
  return hs_test::ray_demonstrator_tests::run_ray_demonstrator_tests();
}
int run_sdf_pattern_tests() {
  return hs_test::sdf_pattern_tests::run_sdf_pattern_tests();
}
int run_ray_event_tests() {
  return hs_test::ray_event_tests::run_ray_event_tests();
}
int run_hyper_lattice_tests() {
  return hs_test::hyper_lattice_tests::run_hyper_lattice_tests();
}
int run_pullback_tests() {
  return hs_test::pullback_tests::run_pullback_tests();
}
int run_shader_chain_tests() {
  return hs_test::shader_chain_tests::run_shader_chain_tests();
}
int run_filter_tests() { return hs_test::filter_tests::run_filter_tests(); }
int run_plot_parity_tests() {
  return hs_test::plot_parity_tests::run_plot_parity_tests();
}
int run_plot_scan_tests() {
  return hs_test::plot_scan_tests::run_plot_scan_tests();
}
int run_canvas_tests() { return hs_test::canvas_tests::run_canvas_tests(); }
int run_scan_tests() { return hs_test::scan_tests::run_scan_tests(); }
int run_mesh_raster_tests() {
  return hs_test::mesh_raster_tests::run_mesh_raster_tests();
}

} // namespace hs_test::module_runners
