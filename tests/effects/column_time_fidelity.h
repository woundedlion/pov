/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Column-time scan effects draw differently with temporal sampling on.
// ---------------------------------------------------------------------------

/**
 * @brief Checks that an effect's column-time scans change its frames when
 *        temporal sampling is on.
 * @tparam E Effect class template.
 * @param name Effect name for failure context.
 */
template <template <int, int> class E>
inline void check_column_time_scans_engage(const char *name) {
  HS_CONTEXT(name);
  constexpr int FRAMES = 12;
  HS_EXPECT_TRUE(
      (pinned_clock_frame_fold<E, SMALL_W, SMALL_H>(FRAMES, nullptr, true)
           .hash) !=
      (pinned_clock_frame_fold<E, SMALL_W, SMALL_H>(FRAMES, nullptr, false)
           .hash));
}

/**
 * @brief Column-time scan effects draw differently with temporal sampling on.
 */
inline void test_column_time_scans_engage_with_fidelity_on() {
  check_column_time_scans_engage<HankinSolids>("HankinSolids");
  check_column_time_scans_engage<IslamicStars>("IslamicStars");
  check_column_time_scans_engage<GnomonicStars>("GnomonicStars");
  check_column_time_scans_engage<ShapeShifter>("ShapeShifter");
  check_column_time_scans_engage<Comets>("Comets");
  check_column_time_scans_engage<Galaxies>("Galaxies");
  check_column_time_scans_engage<DisplacementField>("DisplacementField");
  check_column_time_scans_engage<Raymarch>("Raymarch");
  check_column_time_scans_engage<RingSpin>("RingSpin");
}
