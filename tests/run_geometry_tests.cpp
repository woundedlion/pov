/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_sdf.h"
#include "tests/test_conway.h"
#include "tests/test_partition_seam.h"
#include "tests/test_hankin.h"
#include "tests/test_geometry.h"
#include "tests/test_spherical_field.h"
#include "tests/test_spherical_harmonics.h"
#include "tests/test_mesh.h"
#include "tests/test_solids.h"
#include "tests/test_reaction_graph.h"

namespace hs_test::module_runners {

HS_FORWARD_MODULE(sdf_tests, run_sdf_tests)
HS_FORWARD_MODULE(conway_tests, run_conway_tests)
HS_FORWARD_MODULE(partition_seam_tests, run_partition_seam_tests)
HS_FORWARD_MODULE(hankin_tests, run_hankin_tests)
HS_FORWARD_MODULE(geometry_tests, run_geometry_tests)
HS_FORWARD_MODULE(spherical_field_tests, run_spherical_field_tests)
HS_FORWARD_MODULE(spherical_harmonics_tests, run_spherical_harmonics_tests)
HS_FORWARD_MODULE(mesh_tests, run_mesh_tests)
HS_FORWARD_MODULE(solids_tests, run_solids_tests)
HS_FORWARD_MODULE(reaction_graph_tests, run_reaction_graph_tests)

} // namespace hs_test::module_runners
