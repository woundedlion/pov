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
#include "tests/test_pullback.h"
#include "tests/test_filter.h"
#include "tests/test_plot_scan.h"
#include "tests/test_canvas.h"
#include "tests/test_scan.h"
#include "tests/test_mesh_raster.h"

namespace hs_test::module_runners {

HS_FORWARD_MODULE(ray_tests, run_ray_tests)
HS_FORWARD_MODULE(ray_demonstrator_tests, run_ray_demonstrator_tests)
HS_FORWARD_MODULE(sdf_pattern_tests, run_sdf_pattern_tests)
HS_FORWARD_MODULE(ray_event_tests, run_ray_event_tests)
HS_FORWARD_MODULE(pullback_tests, run_pullback_tests)
HS_FORWARD_MODULE(filter_tests, run_filter_tests)
HS_FORWARD_MODULE(plot_scan_tests, run_plot_scan_tests)
HS_FORWARD_MODULE(canvas_tests, run_canvas_tests)
HS_FORWARD_MODULE(scan_tests, run_scan_tests)
HS_FORWARD_MODULE(mesh_raster_tests, run_mesh_raster_tests)

} // namespace hs_test::module_runners
