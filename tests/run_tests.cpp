#include <cstdio>
#include <cstdio>
#include "tests/test_3dmath.h"
#include "tests/test_memory.h"
#include "tests/test_spatial.h"
#include "tests/test_static_circular_buffer.h"
#include "tests/test_sdf.h"
#include "tests/test_conway.h"
#include "tests/test_hankin.h"
#include "tests/test_geometry.h"
#include "tests/test_mesh.h"

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  int failures = 0;
  failures += hs_test::math3d::run_3dmath_tests();
  failures += hs_test::mem::run_memory_tests();
  failures += hs_test::spatial::run_spatial_tests();
  failures += hs_test::scb::run_static_circular_buffer_tests();
  failures += hs_test::sdf::run_sdf_tests();
  failures += hs_test::conway_tests::run_conway_tests();
  failures += hs_test::hankin_tests::run_hankin_tests();
  failures += hs_test::geometry::run_geometry_tests();
  failures += hs_test::mesh_tests::run_mesh_tests();
  return failures;
}
