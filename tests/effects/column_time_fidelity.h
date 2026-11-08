/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Column-time scan effects draw differently with temporal sampling on.
// ---------------------------------------------------------------------------

/** @brief Renders the selected fidelity mode through the shared fold. */
template <template <int, int> class E, int W, int H, int FRAMES,
          bool ENABLED = false>
inline PixelGolden column_time_fidelity_golden() {
  return pinned_clock_frame_golden<E, W, H>(FRAMES, nullptr, ENABLED);
}

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
      (column_time_fidelity_golden<E, 96, 20, FRAMES, true>().hash) !=
      (column_time_fidelity_golden<E, 96, 20, FRAMES, false>().hash));
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
