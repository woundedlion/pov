/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// sRGB <-> linear LUTs vs. float reference
// ============================================================================

/**
 * @brief Verifies sRGB 0 maps to linear 0 and sRGB 255 to max linear.
 */
inline void test_srgb_to_linear_endpoints() {
  HS_EXPECT_EQ(srgb_to_linear(0), 0);
  HS_EXPECT_EQ(srgb_to_linear(255), 65535);
}

/**
 * @brief Verifies the inverse LUT maps linear 0 to sRGB 0 and max to sRGB 255.
 */
inline void test_linear_to_srgb_endpoints() {
  HS_EXPECT_EQ(linear_to_srgb_lut[0], 0);
  HS_EXPECT_EQ(linear_to_srgb_lut[65535], 255);
}

// The ~1.5 KB split-decode must reproduce the 64 KB linear_to_srgb_lut for every
// one of the 65536 inputs.
inline void test_linear_to_srgb8_decode_matches_lut() {
  long mismatches = 0;
  for (int v = 0; v <= 65535; ++v)
    if (linear_to_srgb8((uint16_t)v) != linear_to_srgb_lut[v])
      ++mismatches;
  HS_EXPECT_EQ(mismatches, 0);
}

/**
 * @brief Pins MIN_ENCODABLE_ALPHA to the normalized first integer linear
 *        channel that encodes off zero.
 * @details Lower integer inputs encode to sRGB 0. This per-sample cutoff is
 *          distinct from the smaller fractional peak that rounds up to that
 *          first input during pixel accumulation.
 */
inline void test_min_encodable_alpha_is_the_encode_floor() {
  const int v = static_cast<int>(MIN_ENCODABLE_ALPHA * 65535.0f + 0.5f);
  HS_EXPECT_EQ(MIN_ENCODABLE_ALPHA, static_cast<float>(v) / 65535.0f);
  HS_EXPECT_GE(linear_to_srgb8(static_cast<uint16_t>(v)), 1);
  for (int u = 0; u < v; ++u)
    HS_EXPECT_EQ(linear_to_srgb8(static_cast<uint16_t>(u)), 0);
}

/**
 * @brief Verifies the 8-bit -> 16-bit linear LUT matches the float reference.
 * @details Matches the float reference (scaled to 16-bit) across all 256 entries
 *          within rounding tolerance.
 */
inline void test_srgb_linear_lut_vs_float_reference() {
  for (int s = 0; s <= 255; ++s) {
    float ref = srgb_to_linear_float(s / 255.0f) * 65535.0f;
    float lut = static_cast<float>(srgb_to_linear(static_cast<uint8_t>(s)));
    HS_EXPECT_NEAR(lut, ref, 0.5f);
  }
}

/**
 * @brief Verifies sRGB -> linear (LUT) -> sRGB (LUT) recovers the 8-bit value.
 */
inline void test_srgb_linear_roundtrip_lut() {
  for (int s = 0; s <= 255; ++s) {
    uint16_t lin = srgb_to_linear(static_cast<uint8_t>(s));
    uint8_t back = linear_to_srgb_lut[lin];
    HS_EXPECT_EQ(static_cast<int>(back), s);
  }
}

/**
 * @brief Verifies srgb_to_linear_interp recovers sub-pixel precision.
 * @details Interpolates between LUT entries so sub-8-bit fractions resolve to
 *          distinct, ordered, monotonic 16-bit values rather than collapsing to
 *          one bucket. Endpoints and 1/255 steps still match the integer LUT.
 */
inline void test_srgb_to_linear_interp_recovers_subpixel_precision() {
  // Endpoints match the integer LUT.
  HS_EXPECT_EQ(srgb_to_linear_interp(0.0f), 0);
  HS_EXPECT_EQ(srgb_to_linear_interp(1.0f), 65535);

  // At exact 1/255 steps the interpolated value equals the integer LUT entry
  // (within float rounding).
  for (int s = 0; s <= 255; ++s) {
    int interp = srgb_to_linear_interp(s / 255.0f);
    int lut = srgb_to_linear(static_cast<uint8_t>(s));
    HS_EXPECT_TRUE(std::abs(interp - lut) <= 2);
  }

  // Two sRGB values in the SAME 8-bit bucket (200) but at different fractions
  // must map to DIFFERENT, ordered 16-bit linear values.
  uint16_t a = srgb_to_linear_interp(200.2f / 255.0f);
  uint16_t b = srgb_to_linear_interp(200.8f / 255.0f);
  HS_EXPECT_TRUE(b > a);
  // ...and both lie between the two LUT entries being interpolated.
  uint16_t e200 = srgb_to_linear(200), e201 = srgb_to_linear(201);
  HS_EXPECT_TRUE(a >= e200 && a <= e201);
  HS_EXPECT_TRUE(b >= e200 && b <= e201);

  // Monotonic non-decreasing across the full range.
  uint16_t prev = 0;
  for (int k = 0; k <= 1000; ++k) {
    uint16_t v = srgb_to_linear_interp(k / 1000.0f);
    if (k > 0)
      HS_EXPECT_TRUE(v >= prev);
    prev = v;
  }
}

/**
 * @brief Verifies the float reference sRGB<->linear round-trip recovers values.
 * @details sRGB -> linear -> sRGB through the float reference functions recovers
 *          the original 8-bit value within half a level.
 */
inline void test_srgb_linear_roundtrip_float() {
  for (int s = 0; s <= 255; ++s) {
    float f = s / 255.0f;
    float lin = srgb_to_linear_float(f);
    float back = linear_to_srgb_float(lin);
    HS_EXPECT_NEAR(back * 255.0f, static_cast<float>(s), 0.5f);
  }
}
