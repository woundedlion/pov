/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief PetalFlow temporal oracle. */

/** @brief Production PetalFlow travel, twist, orientation and fade match frozen columns. */
inline void test_petalflow_temporal_source_oracle() {
  const bool enabled = !Spacetime::flat_surface();
  Spacetime::set_flat_surface(false);
  for (float rho : {0.0f, 1.25f, 2.75f}) {
    auto render = [rho](bool timed, float fraction, float displacement) {
      return PetalFlowWhiteBox::render_temporal_ring(timed, fraction,
                                                     displacement, rho);
    };
    hs_test::seam_compare_render<DEFAULT_W, DEFAULT_H>(
        render, rho == 0.0f
                    ? "PetalFlow equator"
                    : (rho == 1.25f ? "PetalFlow dense" : "PetalFlow fade"));
  }
  Spacetime::set_flat_surface(!enabled);
}
