/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Fishbowl temporal oracle. */

/** @brief Warped Fishbowl endpoints, births and palette match frozen columns. */
inline void test_fishbowl_temporal_trail_oracle() {
  const ScopedTemporalSampling temporal_sampling(true);
  for (bool full : {false, true}) {
    auto render = [full](bool timed, float fraction, float displacement) {
      return FishbowlWhiteBox::render_temporal_trail(timed, fraction,
                                                     displacement, full);
    };
    hs_test::seam_compare_render<DEFAULT_W, DEFAULT_H>(render,
                                                       "Fishbowl trail");
  }
}
