/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ---------------------------------------------------------------------------
// Local arena for sampling.
// ---------------------------------------------------------------------------

/**
 * @brief Draws a cube-sampled, normalized random unit vector (not uniform).
 * @return A unit Vector; draws inside a 0.1-radius ball are rejected and
 *         redrawn, so the normalize is well conditioned.
 */
inline math::Vector rand_unit() {
  for (;;) {
    const float rx = hs::rand_f(-1, 1);
    const float ry = hs::rand_f(-1, 1);
    const float rz = hs::rand_f(-1, 1);
    math::Vector r(rx, ry, rz);
    if (r.length() > 0.1f)
      return r.normalized();
  }
}

/** @brief Backing storage for the module-local sampling arena. */
inline uint8_t plot_scan_arena_buf[256 * 1024];

/**
 * @brief Returns the module-local arena used to bind Fragments for sampling.
 * @return Reference to a function-static Arena over plot_scan_arena_buf.
 * @details Plot::*::sample takes a caller-bound Fragments, so this provides a
 *          dedicated arena rather than relying on the global scratch arena.
 *          The arena is not reset between tests: every caller must wrap its
 *          allocations in a ScratchScope, or allocations leak across tests and
 *          couple correctness to run order.
 */
inline Arena &plot_arena() {
  static Arena a(plot_scan_arena_buf, sizeof(plot_scan_arena_buf));
  return a;
}
