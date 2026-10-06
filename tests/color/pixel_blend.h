/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// lerp16
// ============================================================================

/**
 * @brief Verifies frac=0 and frac=65535 recover the two endpoints exactly.
 */
inline void test_lerp16_endpoints() {
  Pixel a(1000, 2000, 3000);
  Pixel b(40000, 50000, 60000);

  Pixel at0 = a.lerp16(b, 0);
  HS_EXPECT_EQ(at0.r, a.r);
  HS_EXPECT_EQ(at0.g, a.g);
  HS_EXPECT_EQ(at0.b, a.b);

  Pixel at1 = a.lerp16(b, 65535);
  HS_EXPECT_EQ(at1.r, b.r);
  HS_EXPECT_EQ(at1.g, b.g);
  HS_EXPECT_EQ(at1.b, b.b);
}

/**
 * @brief Verifies frac~0.5 lands each channel on (a+b)/2 within rounding error.
 * @details Equal endpoints stay put.
 */
inline void test_lerp16_midpoint() {
  Pixel a(0, 100, 65535);
  Pixel b(65535, 300, 65535);
  Pixel mid = a.lerp16(b, 32768); // ~0.5

  HS_EXPECT_NEAR(static_cast<float>(mid.r), 32767.0f, 2.0f);
  HS_EXPECT_NEAR(static_cast<float>(mid.g), 200.0f, 2.0f);
  HS_EXPECT_EQ(mid.b, 65535); // both endpoints equal -> stays put
}

/**
 * @brief Verifies lerp16 rounds to nearest, not floor.
 * @details The reconstruction tail (x + (x>>16) + 32768) >> 16 adds half a
 *          quantum so the divide rounds; at frac = 49152 (~0.75) round-to-nearest
 *          and floor disagree on red and green.
 */
inline void test_lerp16_rounds_to_nearest() {
  Pixel a(0, 0, 0);
  Pixel b(1, 2, 4);
  Pixel m = a.lerp16(b, 49152); // 0.75
  // True values 0.75 / 1.5 / 3.0 -> round-to-nearest 1 / 2 / 3 (floor: 0 / 1 / 3).
  HS_EXPECT_EQ(m.r, 1);
  HS_EXPECT_EQ(m.g, 2);
  HS_EXPECT_EQ(m.b, 3);
}

/**
 * @brief Verifies every interpolated channel lies within the endpoint envelope.
 * @details One LSB of rounding slack below the minimum.
 */
inline void test_lerp16_bounded() {
  Pixel a(123, 45678, 60000);
  Pixel b(62000, 1000, 12345);
  for (uint32_t f = 0; f <= 65535; f += 4095) {
    Pixel m = a.lerp16(b, static_cast<uint16_t>(f));
    HS_EXPECT_LE(m.r, static_cast<uint16_t>(std::max(a.r, b.r)));
    HS_EXPECT_GE(m.r + 1, static_cast<uint16_t>(std::min(a.r, b.r)));
    HS_EXPECT_LE(m.g, static_cast<uint16_t>(std::max(a.g, b.g)));
    HS_EXPECT_GE(m.g + 1, static_cast<uint16_t>(std::min(a.g, b.g)));
    HS_EXPECT_LE(m.b, static_cast<uint16_t>(std::max(a.b, b.b)));
    HS_EXPECT_GE(m.b + 1, static_cast<uint16_t>(std::min(a.b, b.b)));
  }
}

inline void test_color4_lerp_straight_alpha() {
  const Color4 transparent_red(Pixel(65535, 0, 0), 0.0f);
  const Color4 opaque_blue(Pixel(0, 0, 65535), 1.0f);

  const Color4 midpoint = transparent_red.lerp(opaque_blue, 0.5f);
  HS_EXPECT_NEAR(static_cast<float>(midpoint.color.r), 32768.0f, 1.0f);
  HS_EXPECT_EQ(midpoint.color.g, 0);
  HS_EXPECT_NEAR(static_cast<float>(midpoint.color.b), 32768.0f, 1.0f);
  HS_EXPECT_NEAR(midpoint.alpha, 0.5f, 1e-6f);

  const Color4 before = transparent_red.lerp(opaque_blue, -1.0f);
  const Color4 after = transparent_red.lerp(opaque_blue, 2.0f);
  HS_EXPECT_EQ(before.color.r, transparent_red.color.r);
  HS_EXPECT_NEAR(before.alpha, transparent_red.alpha, 1e-6f);
  HS_EXPECT_EQ(after.color.b, opaque_blue.color.b);
  HS_EXPECT_NEAR(after.alpha, opaque_blue.alpha, 1e-6f);
}

inline void test_blend_outputs_denormal_alpha() {
#if defined(HS_TEST_FAST_MATH)
  hs_test::skip_case(
      __func__,
      "HS_TEST_FAST_MATH: the denormal-alpha rescale does not survive the flag "
      "pair");
#else
  const float alpha = std::numeric_limits<float>::denorm_min();
  const Color4 from(Pixel(1000, 2000, 3000), alpha);
  const Color4 to(Pixel(3000, 4000, 5000), alpha);
  const Color4 blended = blend_outputs(from, to, 0.5f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.r), 2000.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.g), 3000.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.b), 4000.0f, 1.0f);
  HS_EXPECT_EQ(blended.alpha, alpha);

  const Color4 quarter = blend_outputs(from, to, 0.25f);
  HS_EXPECT_NEAR(static_cast<float>(quarter.color.r), 1500.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(quarter.color.g), 2500.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(quarter.color.b), 3500.0f, 1.0f);
  HS_EXPECT_EQ(quarter.alpha, alpha);
#endif
}

/**
 * @brief Verifies mix 0 and mix 1 return the endpoints verbatim.
 * @details A transparent endpoint keeps its RGB, and out-of-range mixes clamp
 *          onto the endpoints.
 */
inline void test_blend_outputs_endpoints_verbatim() {
  const Color4 from(Pixel(1000, 2000, 3000), 0.0f);
  const Color4 to(Pixel(40000, 50000, 60000), 0.75f);
  for (float mix : {0.0f, -1.0f}) {
    const Color4 out = blend_outputs(from, to, mix);
    HS_EXPECT_EQ(out.color.r, from.color.r);
    HS_EXPECT_EQ(out.color.g, from.color.g);
    HS_EXPECT_EQ(out.color.b, from.color.b);
    HS_EXPECT_EQ(out.alpha, from.alpha);
  }
  for (float mix : {1.0f, 2.0f}) {
    const Color4 out = blend_outputs(from, to, mix);
    HS_EXPECT_EQ(out.color.r, to.color.r);
    HS_EXPECT_EQ(out.color.g, to.color.g);
    HS_EXPECT_EQ(out.color.b, to.color.b);
    HS_EXPECT_EQ(out.alpha, to.alpha);
  }
}

inline void test_blend_outputs_tiny_normal_alpha() {
  const float alpha = 4.0f * std::numeric_limits<float>::min();
  const Color4 from(Pixel(1000, 2000, 3000), alpha);
  const Color4 to(Pixel(3000, 4000, 5000), alpha);
  const Color4 blended = blend_outputs(from, to, 0.5f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.r), 2000.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.g), 3000.0f, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(blended.color.b), 4000.0f, 1.0f);
  HS_EXPECT_EQ(blended.alpha, alpha);
}

/**
 * @brief Round-to-nearest div-by-65535 lerp reference computed in double.
 * @param a First endpoint channel value in [0, 65535].
 * @param b Second endpoint channel value in [0, 65535].
 * @param frac Interpolation fraction in [0, 65535] mapping to [0, 1].
 * @return Interpolated channel value in [0, 65535], rounded to nearest.
 */
inline uint16_t lerp16_reference(uint16_t a, uint16_t b, uint16_t frac) {
  double t = static_cast<double>(frac) / 65535.0;
  return static_cast<uint16_t>(a * (1.0 - t) + b * t + 0.5);
}

/**
 * @brief Verifies lerp16 is correct across the full 0..65535 operand range.
 * @details Operands >= 32768 (a frac, an inverse-frac, or a bright channel)
 *          must not be read as negative by a signed multiply.
 */
inline void test_lerp16_full_range_correct() {
  // The midpoint between maximal and zero channels is half scale.
  Pixel hi(65535, 65535, 65535), lo(0, 0, 0);
  Pixel mid = hi.lerp16(lo, 32768);
  HS_EXPECT_NEAR(static_cast<float>(mid.r), 32768.0f, 2.0f);
  HS_EXPECT_NEAR(static_cast<float>(mid.g), 32768.0f, 2.0f);

  // Bright endpoints recovered exactly.
  Pixel a(65535, 49152, 40000), b(32768, 60000, 33000);
  HS_EXPECT_EQ(a.lerp16(b, 0), a);
  HS_EXPECT_EQ(a.lerp16(b, 65535), b);

  // Sweep high-operand pairs (all >= 32768) against the double reference — the
  // regime a signed multiply would corrupt.
  const uint16_t vals[] = {32768, 40000, 49152, 60000, 65535};
  for (uint16_t av : vals)
    for (uint16_t bv : vals)
      for (uint32_t f = 0; f <= 65535; f += 8191) {
        Pixel pa(av, 0, 0), pb(bv, 0, 0);
        uint16_t got = pa.lerp16(pb, static_cast<uint16_t>(f)).r;
        uint16_t ref = lerp16_reference(av, bv, static_cast<uint16_t>(f));
        HS_EXPECT_TRUE(
            std::abs(static_cast<int>(got) - static_cast<int>(ref)) <= 1);
      }
}

// ============================================================================
// Blend modes
// ============================================================================

/**
 * @brief Pins the device's packed (uqadd16) saturating-add lane layout.
 * @details The device path of operator+= packs g|b into one 32-bit
 *          uqadd16 lane and r into another, then unpacks.
 *          pixel_blend_add_packed shares that lane layout via the software
 *          uqadd16 and is checked against a per-channel saturating reference.
 */
inline void test_blend_add_packed_lane_layout() {
  auto ref = [](uint32_t x, uint32_t y) -> uint16_t {
    uint32_t s = x + y;
    return (uint16_t)(s > 65535 ? 65535 : s);
  };
  const Pixel cases[][2] = {
      {Pixel(60000, 1000, 40000), Pixel(10000, 200, 40000)}, // r,b sat; g not
      {Pixel(0, 65535, 0), Pixel(65535, 0, 65535)},          // each lane to max
      {Pixel(123, 45678, 9000), Pixel(40000, 30000, 50)},    // g sat only
      {Pixel(0, 0, 0), Pixel(0, 0, 0)},                      // zero
  };
  for (const auto &c : cases) {
    Pixel got = pixel_blend_add_packed(c[0], c[1]);
    HS_EXPECT_EQ(got.r, ref(c[0].r, c[1].r));
    HS_EXPECT_EQ(got.g, ref(c[0].g, c[1].g));
    HS_EXPECT_EQ(got.b, ref(c[0].b, c[1].b));
    // Host add operators must agree with the packed device layout.
    Pixel acc = c[0];
    acc += c[1];
    HS_EXPECT_EQ(acc, got);
  }
}

/**
 * @brief Verifies blend_alpha clamps its alpha to [0,1] before the float->int cast.
 * @details In-range values interpolate with round-to-nearest (+0.5f) weight
 *          quantization and out-of-range/overflowing/NaN alphas saturate to an
 *          endpoint instead of invoking cast UB.
 */
inline void test_blend_alpha_clamps_before_cast() {
  Pixel a(0, 0, 0);
  Pixel b(60000, 40000, 20000);

  HS_EXPECT_EQ(blend_alpha(0.0f)(a, b), a); // fully a
  HS_EXPECT_EQ(blend_alpha(1.0f)(a, b), b); // fully b

  // Alpha rounds to nearest (+0.5f): 0.5 -> weight 32768, not 32767.
  HS_EXPECT_EQ(blend_alpha(0.5f)(a, b), a.lerp16(b, 32768));

  // Out-of-range alpha saturates: a >= 1 -> full b; a <= 0 -> full a.
  HS_EXPECT_EQ(blend_alpha(1000.0f)(a, b), b);
  HS_EXPECT_EQ(blend_alpha(-5.0f)(a, b), a);
  // Large enough to overflow int in an unclamped (int)(a*65535).
  HS_EXPECT_EQ(blend_alpha(1e9f)(a, b), b);
  // NaN folds to the hi bound via hs::clamp.
  Pixel nan_res = blend_alpha(NAN)(a, b);
  HS_EXPECT_EQ(nan_res, b);
}

/**
 * @brief Verifies Pixel * float clamps each scaled channel into [0,65535]
 *        before the cast.
 * @details Overflowing scales saturate, negatives clamp to 0, and NaN folds to
 *          the hi bound rather than invoking cast UB.
 */
inline void test_pixel_scale_clamps_before_cast() {
  Pixel c(100, 2000, 30000);

  HS_EXPECT_EQ(c * 0.0f, Pixel(0, 0, 0));
  HS_EXPECT_EQ(c * 1.0f, c);
  HS_EXPECT_EQ(c * 2.0f, Pixel(200, 4000, 60000));

  // Half-LSB results round up: odd channels * 0.5 land on .5 and carry up.
  HS_EXPECT_EQ(Pixel(3, 5, 7) * 0.5f, Pixel(2, 3, 4));

  // Overflowing scale saturates at 65535.
  HS_EXPECT_EQ(c * 1e9f, Pixel(65535, 65535, 65535));
  // Negative scale clamps to zero.
  HS_EXPECT_EQ(c * -3.0f, Pixel(0, 0, 0));
  // NaN folds to the hi bound.
  HS_EXPECT_EQ(c * NAN, Pixel(65535, 65535, 65535));
}

/** @brief Quarter-scale accumulation rounds each sample before addition. */
inline void test_pixel_quarter_accumulation_rounds_per_sample() {
  uint32_t sum = 0;
  Pixel accumulated(0, 0, 0);
  for (uint32_t channel = 0; channel <= 65535u; ++channel) {
    Pixel sample(static_cast<uint16_t>(channel), 0, 0);
    accumulated = Pixel(0, 0, 0);
    sum = 0;
    for (int i = 0; i < 4; ++i) {
      accumulated += sample * 0.25f;
      sum += (channel + 2u) >> 2;
    }
    HS_EXPECT_EQ(accumulated.r,
                 static_cast<uint16_t>(sum > 65535u ? 65535u : sum));
  }
}
