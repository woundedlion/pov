/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// Geodesic trigonometric parity.

/**
 * @brief The [0, π] paired trig path is bit-exact with the general functions.
 * @details Probes aggregate into a first-divergence capture and a sample
 *          counter.
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
