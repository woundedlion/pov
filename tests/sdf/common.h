/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
/** @file
 * @brief Includes and basis helper shared by the SDF test sections.
 */

#include <algorithm>
#include <array>
#include <cstdio>

#include "core/render/sdf.h"
#include "core/render/sdf/ring_knots.h"
#include "core/mesh/solids.h"
#include "core/render/sdf/volume.h"
#include "core/render/sdf/face_class_bake.h"
#include "core/render/scan.h"
#include "core/math/geometry.h"
#include "tests/prepared_face_test_access.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <iterator>
#include <limits>
#include <numbers>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test::sdf_tests {

/**
 * @brief Builds the canonical equator-facing basis: v = +Y, u = +X, w = +Z.
 * @return A Basis oriented so its pole points along +Y.
 */
inline math::Basis equator_basis() {
  return math::Basis{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                     math::Vector(0, 0, 1)};
}

} // namespace hs_test::sdf_tests
