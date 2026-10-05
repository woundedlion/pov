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

/**
 * @brief The [0, π] paired trig path is bit-exact with the general functions.
 * @details The ~525k probed angles aggregate into a first-divergence capture
 *          and a sample counter.
 */
inline void test_geodesic_sincos_bit_parity() {
  float first_bad = -1.0f;
  int divergent = 0;
  int probed = 0;
  auto check = [&](float ang) {
    float s, c;
    math::fast_sincosf_0_pi(ang, s, c);
    ++probed;
    if (std::bit_cast<uint32_t>(s) !=
            std::bit_cast<uint32_t>(math::fast_sinf(ang)) ||
        std::bit_cast<uint32_t>(c) !=
            std::bit_cast<uint32_t>(math::fast_cosf(ang))) {
      if (divergent == 0)
        first_bad = ang;
      ++divergent;
    }
  };

  const uint32_t PI_BITS = std::bit_cast<uint32_t>(math::PI_F);
  const uint32_t HALF_PI_BITS = std::bit_cast<uint32_t>(math::PI_F * 0.5f);
  constexpr uint32_t BOUNDARY_ULPS = 65536;
  for (uint32_t bits = 0; bits <= BOUNDARY_ULPS; ++bits)
    check(std::bit_cast<float>(bits));
  for (uint32_t center : {HALF_PI_BITS, PI_BITS}) {
    for (uint32_t bits = center - BOUNDARY_ULPS;
         bits <= center + BOUNDARY_ULPS && bits <= PI_BITS; ++bits)
      check(std::bit_cast<float>(bits));
  }
  constexpr uint32_t STRIDE = 4093;
  for (uint64_t bits = 0; bits <= PI_BITS; bits += STRIDE)
    check(std::bit_cast<float>(static_cast<uint32_t>(bits)));
  check(math::PI_F);

  HS_EXPECT_EQ(divergent, 0);
  HS_EXPECT_EQ(first_bad, -1.0f);
  // A sweep that stopped generating angles would report zero divergences.
  HS_EXPECT_GT(probed, 500000);
}
