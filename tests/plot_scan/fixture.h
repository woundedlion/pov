/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for plotting, column scans, clipping and space-time helpers.
 */
#pragma once
/** @file
 * @brief Shared plot/scan test dependencies and compile-time contracts.
 */
#include <cstdio>
#include <array>
#include "tests/pixel_test_util.h"
#include "tests/temporal_oracle.h"

#include "core/animation/orientation.h"
#include "core/animation/transformer.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/animation/animation.h" // Segue::Dissolve
#include "core/render/filter.h"
#include "core/math/geometry.h"
#include "core/render/canvas.h"
#include "core/render/spacetime/crossing.h"
#include "core/render/spacetime/pose_table.h"
#include "core/render/spacetime/swept_bound.h"
#include "core/render/spacetime/worldlines.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/mesh_test_util.h"
#include "tests/vec_test_util.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>
#include <limits>
#include <type_traits>
#include <vector>

namespace hs_test {
namespace plot_scan_tests {

static_assert(std::is_constructible_v<Scan::ColumnBasisScope, Canvas &,
                                      Scan::ColumnBasis &>);
static_assert(std::is_constructible_v<Scan::ColumnBasisScope, Canvas &,
                                      const Scan::ColumnBasis &>);
static_assert(!std::is_constructible_v<Scan::ColumnBasisScope, Canvas &,
                                       Scan::ColumnBasis &&>);
static_assert(!std::is_constructible_v<Scan::ColumnBasisScope, Canvas &,
                                       const Scan::ColumnBasis &&>);
static_assert(std::is_constructible_v<Scan::DriftColumns<96>, Arena &,
                                      Scan::ColumnBasis &>);
static_assert(std::is_constructible_v<Scan::DriftColumns<96>, Arena &,
                                      const Scan::ColumnBasis &>);
static_assert(!std::is_constructible_v<Scan::DriftColumns<96>, Arena &,
                                       Scan::ColumnBasis &&>);
static_assert(!std::is_constructible_v<Scan::DriftColumns<96>, Arena &,
                                       const Scan::ColumnBasis &&>);

/** @brief Basis category accepted by the borrowing planar projection factory.
 * @tparam Chart Basis value or reference category. */
template <typename Chart>
concept BorrowableRasterBasis = requires(Chart &&basis) {
  Plot::RasterProjection::planar(std::forward<Chart>(basis));
};

static_assert(BorrowableRasterBasis<math::Basis &>);
static_assert(BorrowableRasterBasis<const math::Basis &>);
static_assert(!BorrowableRasterBasis<math::Basis>);
static_assert(!BorrowableRasterBasis<const math::Basis>);

static_assert(std::is_aggregate_v<Plot::RasterOptions>);
static_assert(std::is_trivially_copyable_v<Plot::RasterOptions>);
static_assert(sizeof(void *) != 8 || sizeof(Plot::RasterOptions) <= 88);
static_assert(
    !std::is_constructible_v<Plot::RasterLoop, bool, const Fragment *>);
static_assert(
    !std::is_constructible_v<Plot::RasterProjection, const math::Basis *,
                             std::span<const uint8_t>>);
static_assert(std::is_constructible_v<Plot::PointProjections,
                                      const float (&)[3], const float (&)[3]>);
static_assert(!std::is_constructible_v<Plot::PointProjections,
                                       const float (&)[3], const float (&)[2]>);
static_assert(!std::is_constructible_v<Plot::PointProjections, const float *,
                                       const float *, size_t>);
static_assert(!Plot::RasterLoop{}.is_closed());
static_assert(Plot::RasterLoop::closed().is_closed());

} // namespace plot_scan_tests
} // namespace hs_test
