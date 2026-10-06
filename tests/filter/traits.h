/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Trait structs — direct member values
// ============================================================================

/**
 * @brief Verifies each trait tag exposes the expected is_2d / has_history member
 *        constants.
 */
inline void test_trait_member_values() {
  HS_EXPECT_TRUE(Filter::Is2D::is_2d);
  HS_EXPECT_FALSE(Filter::Is2D::has_history);

  HS_EXPECT_FALSE(Filter::Is3D::is_2d);
  HS_EXPECT_FALSE(Filter::Is3D::has_history);

  HS_EXPECT_TRUE(Filter::Is2DWithHistory::is_2d);
  HS_EXPECT_TRUE(Filter::Is2DWithHistory::has_history);

  HS_EXPECT_FALSE(Filter::Is3DWithHistory::is_2d);
  HS_EXPECT_TRUE(Filter::Is3DWithHistory::has_history);
}
