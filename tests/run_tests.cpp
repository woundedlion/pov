#include "tests/test_3dmath.h"

int main() {
  int failures = 0;
  failures += hs_test::math3d::run_3dmath_tests();
  return failures;
}
