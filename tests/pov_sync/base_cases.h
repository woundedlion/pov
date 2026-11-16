/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Pov sync module cases. */

#include "tests/pov_sync/context.h"

namespace hs_test::pov_sync_tests {

// ── Runner ──────────────────────────────────────────────────────────────────

/** @brief The shipping cadence rejoins within 25 revolutions. */
inline void test_shipping_rejoin_bound() {
  HS_EXPECT_LE(shipping_config().rejoin_bound_revs(), 25u);
}
} // namespace hs_test::pov_sync_tests
