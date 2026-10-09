/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// OKLab / OKLCH round-trips
// ============================================================================

/** Allowed absolute error of an OKLab/OKLCH round trip in 8-bit sRGB levels. */
inline constexpr float ROUNDTRIP_TOL255 = 1.0f;

/**
 * @brief Round-trips one sRGB color through OKLab and asserts recovery.
 * @param r Source red channel in [0, 255].
 * @param g Source green channel in [0, 255].
 * @param b Source blue channel in [0, 255].
 * @details Path: sRGB[0-255] -> linear float -> OKLab -> linear float ->
 *          sRGB[0-255].
 */
inline void roundtrip_oklab(uint8_t r, uint8_t g, uint8_t b) {
  float rf = srgb_to_linear_float(r / 255.0f);
  float gf = srgb_to_linear_float(g / 255.0f);
  float bf = srgb_to_linear_float(b / 255.0f);

  OKLab lab = linear_rgb_to_oklab(rf, gf, bf);
  float r2, g2, b2;
  oklab_to_linear_rgb(lab, r2, g2, b2);

  float rs = linear_to_srgb_float(hs::clamp(r2, 0.0f, 1.0f)) * 255.0f;
  float gs = linear_to_srgb_float(hs::clamp(g2, 0.0f, 1.0f)) * 255.0f;
  float bs = linear_to_srgb_float(hs::clamp(b2, 0.0f, 1.0f)) * 255.0f;

  HS_EXPECT_NEAR(rs, static_cast<float>(r), ROUNDTRIP_TOL255);
  HS_EXPECT_NEAR(gs, static_cast<float>(g), ROUNDTRIP_TOL255);
  HS_EXPECT_NEAR(bs, static_cast<float>(b), ROUNDTRIP_TOL255);
}

/**
 * @brief Verifies sRGB -> OKLab -> sRGB recovers a spread of colors.
 * @details Recovers grays, primaries/secondaries, and an arbitrary color to
 *          within ROUNDTRIP_TOL255 of the original 8-bit value.
 */
inline void test_oklab_roundtrip() {
  // Grays
  roundtrip_oklab(0, 0, 0);
  roundtrip_oklab(128, 128, 128);
  roundtrip_oklab(255, 255, 255);
  // Saturated primaries / secondaries
  roundtrip_oklab(255, 0, 0);
  roundtrip_oklab(0, 255, 0);
  roundtrip_oklab(0, 0, 255);
  roundtrip_oklab(255, 255, 0);
  roundtrip_oklab(0, 255, 255);
  roundtrip_oklab(255, 0, 255);
  // Arbitrary
  roundtrip_oklab(37, 142, 211);
}

/**
 * @brief Round-trips one sRGB color through OKLCH and asserts recovery.
 * @param r Source red channel in [0, 255].
 * @param g Source green channel in [0, 255].
 * @param b Source blue channel in [0, 255].
 * @details Path: sRGB[0-255] -> OKLCH -> OKLab -> linear -> sRGB[0-255].
 */
inline void roundtrip_oklch(uint8_t r, uint8_t g, uint8_t b) {
  OKLCH lch = srgb_to_oklch(r, g, b);
  float r2, g2, b2;
  oklab_to_linear_rgb(oklch_to_oklab(lch), r2, g2, b2);
  float rs = linear_to_srgb_float(hs::clamp(r2, 0.0f, 1.0f)) * 255.0f;
  float gs = linear_to_srgb_float(hs::clamp(g2, 0.0f, 1.0f)) * 255.0f;
  float bs = linear_to_srgb_float(hs::clamp(b2, 0.0f, 1.0f)) * 255.0f;
  HS_EXPECT_NEAR(rs, static_cast<float>(r), ROUNDTRIP_TOL255);
  HS_EXPECT_NEAR(gs, static_cast<float>(g), ROUNDTRIP_TOL255);
  HS_EXPECT_NEAR(bs, static_cast<float>(b), ROUNDTRIP_TOL255);
}

/**
 * @brief Verifies sRGB -> OKLCH -> sRGB recovers grays, primaries, and a sample.
 */
inline void test_oklch_roundtrip() {
  roundtrip_oklch(0, 0, 0);
  roundtrip_oklch(128, 128, 128);
  roundtrip_oklch(255, 255, 255);
  roundtrip_oklch(255, 0, 0);
  roundtrip_oklch(0, 255, 0);
  roundtrip_oklch(0, 0, 255);
  roundtrip_oklch(64, 180, 75);
}

/**
 * @brief Pins sRGB -> OKLab/OKLCH against published reference coordinates.
 * @details Canonical Ottosson sRGB references (white, pure red/green/blue) in
 *          the engine's units: L in [0,1], a/b Cartesian, h in radians.
 */
inline void test_oklab_reference_triples() {
  struct Ref {
    uint8_t r, g, b;
    float L, a, bb; // expected OKLab
  };
  const Ref refs[] = {
      {255, 255, 255, 1.0000f, 0.0000f, 0.0000f},
      {255, 0, 0, 0.6279f, 0.2249f, 0.1258f},
      {0, 255, 0, 0.8664f, -0.2339f, 0.1795f},
      {0, 0, 255, 0.4520f, -0.0324f, -0.3115f},
  };
  const float tol = 4e-3f;
  for (const Ref &x : refs) {
    float rf = srgb_to_linear_float(x.r / 255.0f);
    float gf = srgb_to_linear_float(x.g / 255.0f);
    float bf = srgb_to_linear_float(x.b / 255.0f);
    OKLab lab = linear_rgb_to_oklab(rf, gf, bf);
    HS_EXPECT_NEAR(lab.L, x.L, tol);
    HS_EXPECT_NEAR(lab.a, x.a, tol);
    HS_EXPECT_NEAR(lab.b, x.bb, tol);
  }

  // Pure red in OKLCH: L=0.6279, C=0.2577, h=29.23deg.
  OKLCH red = srgb_to_oklch(255, 0, 0);
  HS_EXPECT_NEAR(red.L, 0.6279f, tol);
  HS_EXPECT_NEAR(red.C, 0.2577f, tol);
  HS_EXPECT_NEAR(red.h, 29.23f * math::PI_F / 180.0f, 5e-3f);
}

/**
 * @brief Verifies pure grays have ~zero chroma in OKLCH.
 */
inline void test_oklch_gray_is_achromatic() {
  OKLCH g0 = srgb_to_oklch(0, 0, 0);
  OKLCH g1 = srgb_to_oklch(128, 128, 128);
  OKLCH g2 = srgb_to_oklch(255, 255, 255);
  HS_EXPECT_NEAR(g0.C, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(g1.C, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(g2.C, 0.0f, 1e-3f);
}

/**
 * @brief Verifies lerp_oklch hue handling for achromatic endpoints.
 * @details Two grays force hue to 0, and a gray/chromatic pair adopts the
 *          chromatic side's hue.
 */
inline void test_lerp_oklch_achromatic_hue() {
  // Both achromatic -> hue forced to 0, endpoints recovered in L.
  OKLCH a = srgb_to_oklch(0, 0, 0);
  OKLCH b = srgb_to_oklch(255, 255, 255);
  OKLCH mid = lerp_oklch(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.h, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(mid.C, 0.0f, 1e-3f);
  HS_EXPECT_GT(mid.L, a.L);
  HS_EXPECT_LT(mid.L, b.L);

  // One achromatic endpoint -> hue taken from the chromatic side.
  OKLCH red = srgb_to_oklch(255, 0, 0);
  OKLCH from_gray = lerp_oklch(a, red, 0.5f);
  HS_EXPECT_NEAR(from_gray.h, red.h, 1e-6f);
}

/**
 * @brief Verifies lerp_oklch interpolates hue along the short arc across the seam.
 * @details atan2f puts h in [-PI, PI], so two hues straddling that seam are only
 *          a small arc apart even though their numeric difference is ~2*PI; the
 *          lerp must pass through the seam at +/-PI, not through h=0.
 */
inline void test_lerp_oklch_shortest_arc_midpoint() {
  const float L = 0.6f, C = 0.15f;

  // Straddle the +/-PI seam: 2.8 and -2.8 rad are ~0.68 rad apart the short way
  // (through +/-PI), ~5.6 rad apart the long way (through 0).
  OKLCH a{L, C, 2.8f};
  OKLCH b{L, C, -2.8f};
  OKLCH mid = lerp_oklch(a, b, 0.5f);
  // Short arc midpoint sits at the seam (+/-PI), not at 0.
  HS_EXPECT_NEAR(std::fabs(mid.h), math::PI_F, 1e-4f);
  HS_EXPECT_NEAR(mid.L, L, 1e-5f);
  HS_EXPECT_NEAR(mid.C, C, 1e-5f);

  // Same seam, opposite winding.
  OKLCH c{L, C, -3.0f};
  OKLCH d{L, C, 3.0f};
  OKLCH mid2 = lerp_oklch(c, d, 0.5f);
  HS_EXPECT_NEAR(std::fabs(mid2.h), math::PI_F, 1e-4f);

  // A non-seam-crossing pair interpolates directly.
  OKLCH e{L, C, 0.5f};
  OKLCH f{L, C, 1.5f};
  HS_EXPECT_NEAR(lerp_oklch(e, f, 0.5f).h, 1.0f, 1e-4f);

  // Antipodal endpoints: the arc follows the sign of b.h - a.h, and swapping
  // the endpoints traverses the same arc.
  OKLCH g{L, C, 0.0f};
  OKLCH pos{L, C, math::PI_F};
  OKLCH neg{L, C, -math::PI_F};
  HS_EXPECT_EQ(lerp_oklch(g, pos, 0.5f).h, 0.5f * math::PI_F);
  HS_EXPECT_EQ(lerp_oklch(pos, g, 0.5f).h, 0.5f * math::PI_F);
  HS_EXPECT_EQ(lerp_oklch(g, neg, 0.5f).h, -0.5f * math::PI_F);
  HS_EXPECT_EQ(lerp_oklch(neg, g, 0.5f).h, -0.5f * math::PI_F);
}

/**
 * @brief Verifies lerp_oklch at amount 0/1 reproduces the endpoints' L and C.
 */
inline void test_lerp_oklch_endpoints() {
  OKLCH a = srgb_to_oklch(200, 30, 30);
  OKLCH b = srgb_to_oklch(30, 30, 200);
  OKLCH at0 = lerp_oklch(a, b, 0.0f);
  OKLCH at1 = lerp_oklch(a, b, 1.0f);
  HS_EXPECT_NEAR(at0.L, a.L, 1e-5f);
  HS_EXPECT_NEAR(at0.C, a.C, 1e-5f);
  HS_EXPECT_NEAR(at1.L, b.L, 1e-5f);
  HS_EXPECT_NEAR(at1.C, b.C, 1e-5f);
}

/**
 * @brief Verifies extrapolating amounts still yield a valid OKLCH.
 * @details Extrapolating amounts must clamp L to [0,1] and
 *          keep C non-negative, so an overshoot yields a valid L and cannot
 *          flip the hue 180deg through negative chroma.
 */
inline void test_lerp_oklch_extrapolation_clamped() {
  OKLCH dark{0.1f, 0.05f, 0.0f};
  OKLCH bright{0.9f, 0.20f, 1.0f};

  OKLCH under = lerp_oklch(dark, bright, -2.0f); // overshoots below dark
  HS_EXPECT_GE(under.L, 0.0f);
  HS_EXPECT_EQ(under.C, 0.0f);

  OKLCH over = lerp_oklch(bright, dark, 3.0f); // overshoots past dark toward 0
  HS_EXPECT_GE(over.L, 0.0f);
  HS_EXPECT_EQ(over.C, 0.0f);

  OKLCH high = lerp_oklch(dark, bright, 5.0f); // overshoots above bright
  HS_EXPECT_LE(high.L, 1.0f);
}

/**
 * @brief Verifies pixel conversion saturation and OKLCH gamut mapping.
 * @details Full lightness maps to white; an in-gamut neutral gray maps to
 *          linear L^3 on every channel, and an in-gamut chromatic colour keeps
 *          its linear sRGB value.
 */
inline void test_oklch_to_pixel_saturates_and_preserves_in_gamut() {
  HS_EXPECT_EQ(float_to_pixel16(1.5f), 65535);
  HS_EXPECT_EQ(float_to_pixel16(-0.25f), 0);
  HS_EXPECT_EQ(float_to_pixel16(NAN), 65535);
  OKLCH vivid{1.0f, 0.4f, 1.0f};
  Pixel hi = oklch_to_pixel(vivid);
  constexpr int WHITE_TOLERANCE = 4;
  HS_EXPECT_GE(hi.r, 65535 - WHITE_TOLERANCE);
  HS_EXPECT_GE(hi.g, 65535 - WHITE_TOLERANCE);
  HS_EXPECT_GE(hi.b, 65535 - WHITE_TOLERANCE);

  OKLCH gray{0.5f, 0.0f, 0.0f};
  Pixel mid = oklch_to_pixel(gray);
  HS_EXPECT_NEAR(mid.r, 0.125f * 65535.0f, 4.0f);
  HS_EXPECT_EQ(mid.r, mid.g);
  HS_EXPECT_EQ(mid.g, mid.b);

  constexpr float ROUND_TRIP_TOL = 16.0f;
  const Pixel green = oklch_to_pixel(srgb_to_oklch(64u, 180u, 75u));
  HS_EXPECT_NEAR(green.r, srgb_to_linear(64u), ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(green.g, srgb_to_linear(180u), ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(green.b, srgb_to_linear(75u), ROUND_TRIP_TOL);
}

// ============================================================================
// Chroma-reduction gamut mapping
// ============================================================================

/**
 * @brief Verifies the chroma-reduction map holds hue and lightness in-gamut.
 * @details A deeply out-of-gamut OKLCH (chroma far past the sRGB cusp at this L)
 *          must map back inside the cube by shrinking chroma: L and hue match
 *          the input within tolerance, and chroma is strictly smaller but still
 *          positive.
 */
inline void test_gamut_clip_preserves_hue() {
  const float L = 0.65f, C = 0.40f;
  for (float h : {0.3f, 1.2f, 2.5f, -1.0f, -2.6f}) {
    OKLab lab = oklch_to_oklab({L, C, h});

    // Precondition: this chroma is unreachable in sRGB at this L.
    float r0, g0, b0;
    oklab_to_linear_rgb(lab, r0, g0, b0);
    HS_EXPECT_FALSE(linear_rgb_in_gamut(r0, g0, b0));

    OKLab mapped = gamut_clip_preserve_chroma(lab);
    float r, g, b;
    oklab_to_linear_rgb(mapped, r, g, b);
    HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));

    OKLCH out = oklab_to_oklch(mapped);
    HS_EXPECT_NEAR(out.L, L, 1e-4f);
    HS_EXPECT_NEAR(wrap_hue_delta(out.h - h), 0.0f, 1e-3f);
    HS_EXPECT_LT(out.C, C);    // chroma reduced...
    HS_EXPECT_GT(out.C, 0.0f); // ...but not crushed
  }
}

/**
 * @brief Pins the shared OKLab matrices and their conversion column order.
 */
inline void test_gamut_refine_matrices_match_the_conversions() {
  HS_CONTEXT("shared OKLab conversion matrices");
  float l, m, s;

  oklab_to_lms_cbrt({1.0f, 0.0f, 0.0f}, l, m, s);
  HS_EXPECT_EQ(l, 1.0f);
  HS_EXPECT_EQ(m, 1.0f);
  HS_EXPECT_EQ(s, 1.0f);

  oklab_to_lms_cbrt({0.0f, 1.0f, 0.0f}, l, m, s);
  HS_EXPECT_EQ(l, 0.3963377774f);
  HS_EXPECT_EQ(m, -0.1055613458f);
  HS_EXPECT_EQ(s, -0.0894841775f);

  oklab_to_lms_cbrt({0.0f, 0.0f, 1.0f}, l, m, s);
  HS_EXPECT_EQ(l, 0.2158037573f);
  HS_EXPECT_EQ(m, -0.0638541728f);
  HS_EXPECT_EQ(s, -1.2914855480f);

  float r, g, b;
  lms_cbrt_to_linear_rgb(1.0f, 0.0f, 0.0f, r, g, b);
  HS_EXPECT_EQ(r, 4.0767416621f);
  HS_EXPECT_EQ(g, -1.2684380046f);
  HS_EXPECT_EQ(b, -0.0041960863f);

  lms_cbrt_to_linear_rgb(0.0f, 1.0f, 0.0f, r, g, b);
  HS_EXPECT_EQ(r, -3.3077115913f);
  HS_EXPECT_EQ(g, 2.6097574011f);
  HS_EXPECT_EQ(b, -0.7034186147f);

  lms_cbrt_to_linear_rgb(0.0f, 0.0f, 1.0f, r, g, b);
  HS_EXPECT_EQ(r, 0.2309699292f);
  HS_EXPECT_EQ(g, -0.3413193965f);
  HS_EXPECT_EQ(b, 1.7076147010f);

  // The cubic's constant term is L^3 in every channel, which holds only while
  // each RGB row sums to one.
  const float L = 0.61f;
  lms_cbrt_to_linear_rgb(L, L, L, r, g, b);
  HS_EXPECT_NEAR(r, L * L * L, 1e-7f);
  HS_EXPECT_NEAR(g, L * L * L, 1e-7f);
  HS_EXPECT_NEAR(b, L * L * L, 1e-7f);
}

/**
 * @brief Linear-RGB triple of an OKLab color in double precision.
 * @param L Lightness.
 * @param a OKLab a.
 * @param b OKLab b.
 * @param r Out: linear red.
 * @param g Out: linear green.
 * @param bl Out: linear blue.
 * @details Double-precision mirror of `oklab_to_linear_rgb` (core/color/color_space.h).
 */
inline void oklab_to_linear_rgb_ref(double L, double a, double b, double &r,
                                    double &g, double &bl) {
  double l_cbrt = L + 0.3963377774 * a + 0.2158037573 * b;
  double m_cbrt = L - 0.1055613458 * a - 0.0638541728 * b;
  double s_cbrt = L - 0.0894841775 * a - 1.2914855480 * b;
  double l = l_cbrt * l_cbrt * l_cbrt, m = m_cbrt * m_cbrt * m_cbrt,
         s = s_cbrt * s_cbrt * s_cbrt;
  r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
  g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
  bl = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s;
}

/**
 * @brief Estimates first-exit chroma from a sampled ray in double precision.
 * @param L Lightness held fixed.
 * @param ad Unit OKLab a of the hue direction.
 * @param bd Unit OKLab b of the hue direction.
 * @param cap Largest chroma considered; returned when all samples are in gamut.
 * @return The in-gamut bisection endpoint preceding the first outside sample.
 * @details Scans 512 intervals and bisects the first sampled exit. The gate's
 *          tolerance permits disconnected in-gamut intervals, so an outside
 *          interval narrower than the sampling step can be missed.
 */
inline double gamut_first_exit_ref(double L, double ad, double bd, double cap) {
  const double lo_b = -1e-4, hi_b = 1.0 + 1e-4;
  auto inside = [&](double c) {
    double r, g, b;
    oklab_to_linear_rgb_ref(L, ad * c, bd * c, r, g, b);
    return r >= lo_b && r <= hi_b && g >= lo_b && g <= hi_b && b >= lo_b &&
           b <= hi_b;
  };
  const int coarse = 512;
  int hit = -1;
  for (int i = 1; i <= coarse; ++i) {
    if (!inside(cap * i / coarse)) {
      hit = i;
      break;
    }
  }
  if (hit < 0)
    return cap;
  double lo = cap * (hit - 1) / coarse, hi = cap * hit / coarse;
  for (int i = 0; i < 50; ++i) {
    double mid = 0.5 * (lo + hi);
    if (inside(mid))
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

/**
 * @brief Verifies an out-of-gamut lower bound refines over the whole ray.
 * @details With @p lo outside the gamut the search covers [0, lo] and must
 *          return an in-gamut scale within lo/256 below the first exit.
 */
inline void test_gamut_bracket_refine_out_of_gamut_lower_bound() {
  const float LO = 0.45f, HI = 0.5f;
  for (int il = 3; il <= 8; ++il) {
    const float L = il / 10.0f;
    for (int ih = 0; ih < 32; ++ih) {
      HS_CONTEXT("L tenths, hue step", il, ih);
      const float h = math::TWO_PI_F * ih / 32.0f;
      const float a = cosf(h), b = sinf(h);
      float r, g, bl;
      oklab_to_linear_rgb({L, a * LO, b * LO}, r, g, bl);
      HS_EXPECT_FALSE(linear_rgb_in_gamut(r, g, bl));

      const float got = gamut_bracket_refine(L, a, b, LO, HI);
      oklab_to_linear_rgb({L, a * got, b * got}, r, g, bl);
      HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, bl));
      const float ref = static_cast<float>(gamut_first_exit_ref(L, a, b, LO));
      HS_EXPECT_LE(ref - got, LO / 256.0f + 1e-5f);
      HS_EXPECT_LE(got - ref, 1e-5f);
    }
  }
}

/**
 * @brief Bounds the angular gamut-boundary lookup against the sampled exit
 *        reference and checks agreement with the direction overload.
 * @details Sampled rays must stay within the overshoot and deficit bounds of the
 *          double-precision reference. Lightness runs to both ends of the
 *          sixteenth grid, where the bracket widens.
 */
inline void test_gamut_direction_lookup_matches_angle() {
  const float DEFICIT_BOUND = 5e-3f;
  const float OVERSAT_BOUND = 1e-6f;
  // Past the largest sRGB chroma in OKLab, so the ray always exits.
  const double CAP = 0.5;
  for (int il = 1; il < 16; ++il) {
    const float L = il / 16.0f;
    for (int ih = 0; ih < 64; ++ih) {
      const float h = math::TWO_PI_F * ih / 64.0f;
      const float got = gamut_max_chroma(L, h);
      HS_EXPECT_EQ(got, gamut_max_chroma(L, cosf(h), sinf(h)));
      const float ref =
          static_cast<float>(gamut_first_exit_ref(L, cosf(h), sinf(h), CAP));
      HS_EXPECT_LT(ref - got, DEFICIT_BOUND);
      HS_EXPECT_LE(got - ref, OVERSAT_BOUND);
    }
  }
}

/**
 * @brief Bounds clipping's deficit and residue against the sampled exit reference.
 * @details Reference overshoot beyond float rounding is limited in frequency
 *          and magnitude, including a ray with a disconnected in-gamut interval.
 * @param path Label naming which reduction path is armed, for failure context.
 */
inline void expect_clip_lands_on_first_exit(const char *path) {
  HS_CONTEXT(path);
  const float DEFICIT_BOUND = 5e-3f;
  // Below the LUT chroma quantum (1/GAMUT_LUT_SCALE).
  const float OVERSAT_BOUND = 1e-6f;
  const float RESIDUE_BOUND = 0.05f;
  const double CHROMA_IN[3] = {0.6, 0.35, 0.25};
  float worst_deficit = 0.0f, worst_oversat = -1.0f;
  size_t samples = 0, residues = 0;
  auto probe = [&](double lightness, double degrees, double chroma) {
    const double HUE = degrees * 3.14159265358979323846 / 180.0;
    const double AD = std::cos(HUE), BD = std::sin(HUE);
    const OKLab MAPPED = gamut_clip_preserve_chroma(
        {static_cast<float>(lightness), static_cast<float>(AD * chroma),
         static_cast<float>(BD * chroma)});
    float r, g, b;
    oklab_to_linear_rgb(MAPPED, r, g, b);
    HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));
    const float GOT = std::sqrt(MAPPED.a * MAPPED.a + MAPPED.b * MAPPED.b);
    const float REF =
        static_cast<float>(gamut_first_exit_ref(lightness, AD, BD, chroma));
    worst_deficit = fold_worst(worst_deficit, REF - GOT);
    worst_oversat = fold_worst(worst_oversat, GOT - REF);
    residues += GOT - REF > OVERSAT_BOUND;
    ++samples;
    return GOT - REF;
  };

  for (int il = 0; il <= 32; ++il) {
    const double L = 0.1 + 0.8 * il / 32.0;
    for (int ih = 0; ih < 360 + 96; ++ih) {
      // 360 even steps, then a fine fan across the blue vertex.
      const double deg = ih < 360 ? ih : 258.0 + 0.125 * (ih - 360);
      for (double chroma : CHROMA_IN)
        probe(L, deg, chroma);
    }
  }
  HS_EXPECT_GT(probe(0.165, 263.75, 0.6), OVERSAT_BOUND);
  HS_EXPECT_LT(worst_deficit, DEFICIT_BOUND);
  HS_EXPECT_LE(worst_oversat, RESIDUE_BOUND);
  HS_EXPECT_LE(residues, (samples + 49999) / 50000);
}

/**
 * @brief Bounds flash-master clipping against the sampled exit reference.
 * @details Runs the shared sweep with no arena LUT armed.
 */
inline void test_gamut_master_clip_lands_on_first_exit() {
  expect_clip_lands_on_first_exit("flash master");
}

inline void test_gamut_continuous_chroma_is_smooth_and_in_gamut() {
  for (int il = 1; il < 100; ++il) {
    const float L = il / 100.0f;
    for (int ih = 0; ih < 720; ++ih) {
      const float h = 2.0f * math::PI_F * ih / 720.0f;
      const float C = gamut_continuous_chroma(L, h);
      float r, g, b;
      oklab_to_linear_rgb(oklch_to_oklab({L, C, h}), r, g, b);
      HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));
    }
  }

  const float L = 0.43f;
  float previous = gamut_continuous_chroma(L, 4.58f);
  for (int i = 1; i <= 1000; ++i) {
    const float h = 4.58f + 0.06f * i / 1000.0f;
    const float current = gamut_continuous_chroma(L, h);
    HS_EXPECT_LT(fabsf(current - previous), 5e-4f);
    previous = current;
  }
}

// Coarsest supported grid.
inline constexpr int TEST_GAMUT_ANGLE_STEPS = GAMUT_LUT_MIN_ANGLE_STEPS;
inline constexpr int TEST_GAMUT_L_STEPS = GAMUT_LUT_MIN_L_STEPS;

/**
 * @brief Bounds coarsest-grid clipping against the sampled exit reference.
 * @details Uses the same in-gamut and chroma bounds as the flash-master sweep.
 */
inline void test_gamut_lut_clip_lands_on_first_exit() {
  alignas(uint16_t) static uint8_t
      lut_buf[gamut_lut_bytes(TEST_GAMUT_ANGLE_STEPS, TEST_GAMUT_L_STEPS)];
  Arena lut_arena(lut_buf, sizeof(lut_buf));
  init_gamut_lut(lut_arena, TEST_GAMUT_ANGLE_STEPS, TEST_GAMUT_L_STEPS);

  expect_clip_lands_on_first_exit("bracket LUT");

  // Outside the reported lightness band the bracket widens, but the in-gamut
  // guarantee is structural and must still hold at every L and hue.
  for (int il = 0; il <= 200; ++il) {
    const float L = il / 200.0f;
    for (int ih = 0; ih < 180; ++ih) {
      const float h = 6.28318531f * ih / 180.0f;
      for (float cin : {0.05f, 0.2f, 0.45f}) {
        OKLab mapped = gamut_clip_preserve_chroma(
            {L, cin * std::cos(h), cin * std::sin(h)});
        float r, g, b;
        oklab_to_linear_rgb(mapped, r, g, b);
        HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));
      }
    }
  }

  release_gamut_lut();
}

/**
 * @brief Verifies the downsample keeps every merged cell inside the coarse
 *        bracket.
 * @details Each coarse bracket contains every master-table bracket merged
 *          into that cell.
 */
inline void test_gamut_lut_downsample_preserves_bracket() {
  // Half the master on both axes, so cells actually merge.
  constexpr int A = GAMUT_LUT_ANGLE_STEPS / 2, NL = GAMUT_LUT_L_STEPS / 2;
  const int sa = GAMUT_LUT_ANGLE_STEPS / A, sl = GAMUT_LUT_L_STEPS / NL;
  alignas(uint16_t) static uint8_t lut_buf[gamut_lut_bytes(A, NL)];
  Arena lut_arena(lut_buf, sizeof(lut_buf));
  init_gamut_lut(lut_arena, A, NL);
  HS_EXPECT_TRUE(g_gamut_lut.angle_steps == A);
  HS_EXPECT_TRUE(g_gamut_lut.l_steps == NL);

  for (int l = 0; l < NL; ++l)
    for (int a = 0; a < A; ++a) {
      const uint16_t c_lo = g_gamut_lut.table[(l * A + a) * 2];
      const uint16_t c_hi = g_gamut_lut.table[(l * A + a) * 2 + 1];
      HS_EXPECT_LE(c_lo, c_hi);
      for (int dl = 0; dl < sl; ++dl)
        for (int da = 0; da < sa; ++da) {
          const int f =
              ((l * sl + dl) * GAMUT_LUT_ANGLE_STEPS + a * sa + da) * 2;
          HS_EXPECT_LE(c_lo, GAMUT_LUT[f]);
          HS_EXPECT_LE(GAMUT_LUT[f + 1], c_hi);
        }
    }

  release_gamut_lut();
}

inline void test_gamut_cell_nonfinite_coordinates() {
  const auto cell = gamut_cell(g_gamut_lut, NAN, NAN, NAN);
  HS_EXPECT_EQ(cell.angle_index, g_gamut_lut.angle_steps - 1);
  HS_EXPECT_EQ(cell.lightness_index, g_gamut_lut.l_steps - 1);
  for (float lightness : {-INFINITY, INFINITY}) {
    const auto boundary = gamut_cell(g_gamut_lut, lightness, 1.0f, 0.0f);
    HS_EXPECT_EQ(boundary.lightness_index,
                 lightness < 0 ? 0 : g_gamut_lut.l_steps - 1);
  }
}

/**
 * @brief Verifies an in-gamut color survives the arena-copy clip untouched and
 *        that releasing the copy leaves the clip working off the flash master.
 */
inline void test_gamut_lut_release_and_passthrough() {
  alignas(uint16_t) static uint8_t
      lut_buf[gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS)];
  Arena lut_arena(lut_buf, sizeof(lut_buf));
  init_gamut_lut(lut_arena, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);

  // Deep inside the cell minimum: returned unchanged, bit for bit.
  OKLab deep = oklch_to_oklab({0.5f, 0.02f, 1.0f});
  OKLab kept = gamut_clip_preserve_chroma(deep);
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.L),
               std::bit_cast<uint32_t>(deep.L));
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.a),
               std::bit_cast<uint32_t>(deep.a));
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.b),
               std::bit_cast<uint32_t>(deep.b));

  // Just inside the boundary but past the cell minimum: refined, not reduced.
  OKLab near_edge = oklch_to_oklab({0.6f, 0.144f, 1.0f});
  const auto cell =
      gamut_cell(g_gamut_lut, near_edge.L, near_edge.a, near_edge.b);
  const size_t cell_offset =
      (cell.lightness_index * g_gamut_lut.angle_steps + cell.angle_index) * 2;
  const float c_lo = g_gamut_lut.table[cell_offset] * GAMUT_LUT_INV_SCALE;
  HS_EXPECT_GT(0.144f, c_lo);
  kept = gamut_clip_preserve_chroma(near_edge);
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.L),
               std::bit_cast<uint32_t>(near_edge.L));
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.a),
               std::bit_cast<uint32_t>(near_edge.a));
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(kept.b),
               std::bit_cast<uint32_t>(near_edge.b));

  // An achromatic input must not divide by zero on the way through.
  OKLab gray = gamut_clip_preserve_chroma({0.5f, 0.0f, 0.0f});
  HS_EXPECT_NEAR(gray.a, 0.0f, 1e-9f);
  HS_EXPECT_NEAR(gray.b, 0.0f, 1e-9f);

  // Released: the pointer must leave the arena rather than dangle into it, and
  // land on the full-resolution flash master.
  release_gamut_lut();
  HS_EXPECT_TRUE(g_gamut_lut.table == GAMUT_LUT);
  HS_EXPECT_EQ(g_gamut_lut.angle_steps, GAMUT_LUT_ANGLE_STEPS);
  HS_EXPECT_EQ(g_gamut_lut.l_steps, GAMUT_LUT_L_STEPS);

  // Off the master the clip still maps a past-cusp color in.
  OKLab mapped =
      gamut_clip_preserve_chroma(oklch_to_oklab({0.65f, 0.4f, 1.2f}));
  float r, g, b;
  oklab_to_linear_rgb(mapped, r, g, b);
  HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));
}

/**
 * @brief Pins the single-step normalization used by the LUT gamut path.
 * @details The tolerance is relative because this module also runs under
 * -ffast-math, which reassociates the Newton step.
 */
inline void test_gamut_lut_boundary_scale_rounding() {
  const OKLab scaled = gamut_scale_to_boundary_lut({0.5f, 0.4f, 0.3f});
  HS_EXPECT_EQ(scaled.L, 0.5f);
  HS_EXPECT_NEAR_REL(scaled.a, 0.13113570f, 1e-6);
  HS_EXPECT_NEAR_REL(scaled.b, 0.09835178f, 1e-6);
  // One scale factor on (a, b): the input 0.4:0.3 ratio survives.
  HS_EXPECT_NEAR(scaled.a * 0.3f, scaled.b * 0.4f, 1e-7f);
  // The single Newton step is one-sided low, so the rescale lands in gamut.
  float r, g, b;
  oklab_to_linear_rgb(scaled, r, g, b);
  HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));
}

/**
 * @brief Verifies configure_arenas() drops an arena-resident copy.
 * @details The copy lives in the persistent arena, and configure_arenas() hands
 *          that storage out again.
 */
inline void test_configure_arenas_releases_gamut_lut() {
  init_gamut_lut(persistent_arena, TEST_GAMUT_ANGLE_STEPS, TEST_GAMUT_L_STEPS);
  HS_EXPECT_TRUE(g_gamut_lut.table != GAMUT_LUT);

  configure_arenas_default();
  HS_EXPECT_TRUE(g_gamut_lut.table == GAMUT_LUT);
}

/**
 * @brief Verifies oklch_to_pixel routes out-of-gamut colors through the
 *        chroma-reduction map, holding hue.
 * @details Realizes a past-cusp OKLCH as a Pixel, reads the realized color back
 *          through the exact forward transform, and checks the hue survived the
 *          16-bit quantization.
 */
inline void test_oklch_to_pixel_holds_hue_out_of_gamut() {
  const float L = 0.62f, C = 0.42f, h = 0.9f;
  Pixel p = oklch_to_pixel({L, C, h});

  OKLCH got = pixel_to_oklch(p);
  HS_EXPECT_NEAR(wrap_hue_delta(got.h - h), 0.0f, 2e-2f);
}

// ============================================================================
// Perceptual hue_rotate (OKLab)
// ============================================================================

/**
 * @brief Verifies hue_sincos tracks libm cosine and sine of a turn angle.
 * @details Covers turns outside [-0.5, 0.5) so the wrap is exercised; the
 *          approximation error stays under 2e-3.
 */
inline void test_hue_sincos_matches_libm() {
  for (int i = -192; i <= 384; ++i) {
    HS_CONTEXT("turn 192ths", i);
    const float turns = i / 192.0f + 0.0013f;
    float cosine, sine;
    hue_sincos(turns, cosine, sine);
    const double angle = 2.0 * 3.14159265358979323846 * turns;
    HS_EXPECT_NEAR(cosine, static_cast<float>(std::cos(angle)), 2e-3f);
    HS_EXPECT_NEAR(sine, static_cast<float>(std::sin(angle)), 2e-3f);
  }
}

/**
 * @brief Verifies a perceptual hue rotation leaves a gray unchanged.
 * @details A gray has zero chroma in OKLab, so the color must come back
 *          unchanged (within 1 LSB) for any amount, with alpha preserved.
 */
inline void test_hue_rotate_preserves_gray() {
  Color4 gray(128, 128, 128, 0.5f);
  for (float amt : {0.1f, 0.25f, 0.5f, 0.8f}) {
    Color4 out = hue_rotate(gray, amt);
    HS_EXPECT_NEAR(static_cast<float>(out.color.r),
                   static_cast<float>(gray.color.r), 1.0f);
    HS_EXPECT_NEAR(static_cast<float>(out.color.g),
                   static_cast<float>(gray.color.g), 1.0f);
    HS_EXPECT_NEAR(static_cast<float>(out.color.b),
                   static_cast<float>(gray.color.b), 1.0f);
    HS_EXPECT_NEAR(out.alpha, 0.5f, 1e-5f);
  }
}

/**
 * @brief Verifies a full-turn rotation returns to the original color.
 * @details The residual is the combined error of fast_cosf/fast_sinf at 2*PI
 *          and the fast_cbrt OKLab round-trip.
 */
inline void test_hue_rotate_full_turn_identity() {
  Color4 c(200, 60, 30, 1.0f);
  Color4 out = hue_rotate(c, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.r), static_cast<float>(c.color.r),
                 12.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.g), static_cast<float>(c.color.g),
                 12.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.b), static_cast<float>(c.color.b),
                 12.0f);
}

/**
 * @brief Verifies a full turn taken in N steps holds hue, chroma and lightness.
 * @details Applying a 1/32-turn rotation 32 times must land back on the input;
 *          per-step requantization and fast-math errors compound. Chroma may
 *          only shrink: a gamut-boundary color loses chroma where the clip
 *          pulls it in, while a base well inside the cusp holds it.
 */
inline void test_hue_rotate_full_turn_in_steps_holds_hue_and_chroma() {
  const int STEPS = 32;
  // Second entry is the interior base; the rest ride the gamut boundary.
  const Color4 cases[] = {Color4(200, 60, 30, 1.0f), Color4(90, 150, 200, 1.0f),
                          Color4(240, 250, 60, 1.0f), Color4(255, 0, 0, 1.0f)};
  for (int ci = 0; ci < 4; ++ci) {
    Color4 out = cases[ci];
    for (int i = 0; i < STEPS; ++i)
      out = hue_rotate(out, 1.0f / STEPS);

    const OKLCH before = pixel_to_oklch(cases[ci].color);
    const OKLCH after = pixel_to_oklch(out.color);
    HS_EXPECT_NEAR(wrap_hue_delta(after.h - before.h), 0.0f, 0.09f);
    HS_EXPECT_NEAR(after.L, before.L, 1e-3f);
    HS_EXPECT_LE(after.C, before.C + 1e-4f);
    if (ci == 1)
      HS_EXPECT_NEAR(after.C, before.C, 2e-4f);
  }
}

/**
 * @brief Verifies the precomputed-base hue rotation matches hue_rotate exactly.
 * @details Both paths run the same op sequence (forward transform hoisted vs
 *          inline), so the 16-bit channels must match bit-exactly across
 *          saturated, gamut-clipping, and achromatic bases.
 */
inline void test_hue_rotate_base_matches_direct() {
  const Color4 cases[] = {Color4(Pixel(52000, 9000, 3000), 0.8f),
                          Color4(Pixel(65535, 0, 40000), 1.0f),
                          Color4(128, 128, 128, 0.5f)};
  for (const Color4 &c : cases) {
    HueRotateBase hb = make_hue_rotate_base(c);
    for (float amt : {0.0f, 0.1f, 0.37f, 0.5f, 0.93f}) {
      Color4 direct = hue_rotate(c, amt);
      Color4 based = hue_rotate(hb, amt);
      HS_EXPECT_EQ(based.color.r, direct.color.r);
      HS_EXPECT_EQ(based.color.g, direct.color.g);
      HS_EXPECT_EQ(based.color.b, direct.color.b);
      HS_EXPECT_NEAR(based.alpha, direct.alpha, 1e-6f);
    }
  }
}
