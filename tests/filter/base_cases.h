/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief H_OFFSET precondition for the filter cases. */

#include "tests/filter/context.h"

namespace hs_test::filter_tests {

static_assert(hs::H_OFFSET == 0,
              "tests/filter/ cases assume row H-1 is the south pole; the "
              "offset-3 mapping is covered by tests/test_h_offset_renorm.h");
} // namespace hs_test::filter_tests
