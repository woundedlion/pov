/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

#include "tests/test_sdf.h"
#include "tests/test_conway.h"
#include "tests/test_conway_morph.h"
#include "tests/test_conway_continuity.h"
#include "tests/test_partition_seam.h"
#include "tests/test_conway_soak.h"
#include "tests/test_opchain_probe.h"
#include "tests/test_opchain_arena_survey.h"
#include "tests/test_hankin.h"
#include "tests/test_geometry.h"
#include "tests/test_spherical_field.h"
#include "tests/test_spherical_harmonics.h"
#include "tests/test_mesh.h"
#include "tests/test_solids.h"
#include "tests/test_reaction_graph.h"

namespace hs_test::module_runners {

int run_sdf_tests() { return hs_test::sdf_tests::run_sdf_tests(); }
int run_conway_tests() { return hs_test::conway_tests::run_conway_tests(); }
int run_conway_morph_tests() {
  return hs_test::conway_morph_tests::run_conway_morph_tests();
}
int run_conway_continuity_tests() {
  return hs_test::conway_continuity_tests::run_conway_continuity_tests();
}
int run_partition_seam_tests() {
  return hs_test::partition_seam_tests::run_partition_seam_tests();
}
int run_conway_soak_tests() {
  return hs_test::conway_soak_tests::run_conway_soak_tests();
}
int run_opchain_probe_tests() {
  return hs_test::opchain_probe_tests::run_opchain_probe_tests();
}
int run_opchain_arena_survey_tests() {
  return hs_test::opchain_arena_survey_tests::run_opchain_arena_survey_tests();
}
int run_hankin_tests() { return hs_test::hankin_tests::run_hankin_tests(); }
int run_geometry_tests() {
  return hs_test::geometry_tests::run_geometry_tests();
}
int run_spherical_field_tests() {
  return hs_test::spherical_field_tests::run_spherical_field_tests();
}
int run_spherical_harmonics_tests() {
  return hs_test::spherical_harmonics_tests::run_spherical_harmonics_tests();
}
int run_mesh_tests() { return hs_test::mesh_tests::run_mesh_tests(); }
int run_solids_tests() { return hs_test::solids_tests::run_solids_tests(); }
int run_reaction_graph_tests() {
  return hs_test::reaction_graph_tests::run_reaction_graph_tests();
}

} // namespace hs_test::module_runners
