/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// OrientationTrail
// ----------------------------------------------------------------------------
// Index 0 is the oldest snapshot, length()-1 the newest.
// ============================================================================

/**
 * @brief Verifies recorded snapshots are ordered oldest-first: index 0 is the
 * first recorded, the last index the newest.
 */
inline void test_orientation_trail_index_zero_is_oldest() {
  Animation::OrientationTrail<math::Orientation<4>, 8> trail;

  math::Orientation<4> a, b, c;
  a.set(math::make_rotation(math::Vector(0, 1, 0), 0.1f));
  b.set(math::make_rotation(math::Vector(0, 1, 0), 0.2f));
  c.set(math::make_rotation(math::Vector(0, 1, 0), 0.3f));

  trail.record(a);
  trail.record(b);
  trail.record(c);
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(3));

  // Index 0 == oldest (a); last index == newest (c).
  HS_EXPECT_NEAR(std::abs(math::dot(trail.get(0).get(), a.get())), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(std::abs(math::dot(trail.get(2).get(), c.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies expire() removes the oldest snapshot, leaving the survivor at
 * index 0.
 */
inline void test_orientation_trail_expire_drops_oldest() {
  Animation::OrientationTrail<math::Orientation<4>, 8> trail;
  math::Orientation<4> a, b;
  a.set(math::make_rotation(math::Vector(0, 1, 0), 0.1f));
  b.set(math::make_rotation(math::Vector(0, 1, 0), 0.2f));
  trail.record(a);
  trail.record(b);
  trail.expire(); // removes the oldest (a)
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(std::abs(math::dot(trail.get(0).get(), b.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies clear() empties the trail.
 */
inline void test_orientation_trail_clear() {
  Animation::OrientationTrail<math::Orientation<4>, 8> trail;
  math::Orientation<4> a;
  trail.record(a);
  trail.record(a);
  trail.clear();
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(0));
}
