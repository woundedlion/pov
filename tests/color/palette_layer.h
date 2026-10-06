/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Palette layer: source palettes, modifiers, compositions, and wrappers
// ============================================================================

/**
 * @brief Verifies ProceduralPalette evaluates its cosine color formula.
 * @details C(t) = a + b*cos(2*PI*(c*t + d)) in sRGB, then to linear. A
 *          {a=.5,b=.5,c=1,d=0} channel is cos-driven: t=0 -> 1.0 (full),
 *          t=0.5 -> 0.0.
 */
inline void test_procedural_palette_cosine() {
  ProceduralPalette pp({0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f},
                       {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f});
  Color4 c0 = pp.get(0.0f);
  Color4 chalf = pp.get(0.5f);
  HS_EXPECT_EQ(c0.color.r, 65535); // 0.5 + 0.5*cos(0) = 1.0
  HS_EXPECT_EQ(chalf.color.r, 0);  // 0.5 + 0.5*cos(PI) = 0.0
  HS_EXPECT_NEAR(c0.alpha, 1.0f, 1e-6f);
}

/**
 * @brief Verifies MutatingPalette blends between its two endpoint palettes.
 * @details mutate(0) reproduces palette #1, mutate(1) reproduces #2, and an
 *          interior amount lands between them.
 */
inline void test_mutating_palette_blends_endpoints() {
  // #1 -> white at t=0 ; #2 -> black everywhere.
  MutatingPalette m({0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {1.0f, 1.0f, 1.0f},
                    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
                    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f});
  m.mutate(0.0f);
  HS_EXPECT_EQ(m.get(0.0f).color.r, 65535);
  m.mutate(1.0f);
  HS_EXPECT_EQ(m.get(0.0f).color.r, 0);
  m.mutate(0.5f); // a=.25, b=.25 -> 0.25 + 0.25*cos(0) = 0.5 sRGB
  uint16_t mid = m.get(0.0f).color.r;
  HS_EXPECT_GT(mid, 0);
  HS_EXPECT_LT(mid, 65535);
}

/**
 * @brief Verifies each coordinate modifier's modify() in isolation.
 * @details Modifiers transform the palette coordinate; this exercises Scale,
 *          Cycle, Quantize, Fold, Breathe, Ripple, and Pinch one at a time.
 */
inline void test_palette_modifiers() {
  // Scale multiplies the coordinate.
  HS_EXPECT_NEAR(ScaleModifier(2.0f).modify(0.3f), 0.6f, 1e-5f);
  float dyn_scale = 3.0f;
  HS_EXPECT_NEAR(ScaleModifier(1.0f, &dyn_scale).modify(0.2f), 0.6f, 1e-5f);

  // Cycle adds the driver offset; a null driver passes through.
  float off = 0.25f;
  HS_EXPECT_NEAR(CycleModifier(&off).modify(0.5f), 0.75f, 1e-5f);
  HS_EXPECT_NEAR(CycleModifier(nullptr).modify(0.5f), 0.5f, 1e-5f);

  // Quantize snaps to the nearest step: round(0.3*4)/4 = 1/4 = 0.25.
  HS_EXPECT_NEAR(QuantizeModifier(4.0f).modify(0.3f), 0.25f, 1e-5f);
  // The top band reaches the endpoint exactly (bounded_output).
  HS_EXPECT_NEAR(QuantizeModifier(4.0f).modify(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(QuantizeModifier(4.0f).modify(0.9f), 1.0f, 1e-5f);

  // Fold is a triangle wave: 0->1, 0.25->0.5, 0.5->0.
  FoldModifier fold(2.0f);
  HS_EXPECT_NEAR(fold.modify(0.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(fold.modify(0.25f), 0.5f, 1e-5f);
  HS_EXPECT_NEAR(fold.modify(0.5f), 0.0f, 1e-5f);

  // Breathe is identity at phase zero and refreshes its memo when phase changes.
  float phase0 = 0.0f;
  BreatheModifier breathe(&phase0, .1f);
  HS_EXPECT_NEAR(breathe.modify(.5f), .5f, 1e-3f);
  phase0 = math::PI_F * .5f;
  HS_EXPECT_NEAR(breathe.modify(.5f), BreatheModifier(&phase0, .1f).modify(.5f),
                 1e-6f);

  // Ripple at t=0, phase=0 leaves the coordinate fixed (sin(0) = 0).
  float rphase = 0.0f;
  HS_EXPECT_NEAR(RippleModifier(&rphase, 3.0f, 0.1f).modify(0.0f), 0.0f, 1e-5f);

  // Pinch (positive tension) pulls an off-center coordinate toward 0.5.
  float tension = 0.5f;
  float pinched = PinchModifier(&tension).modify(0.25f);
  HS_EXPECT_GT(pinched, 0.25f);
  HS_EXPECT_LT(pinched, 0.5f);
  // A null tension driver passes through.
  HS_EXPECT_NEAR(PinchModifier(nullptr).modify(0.3f), 0.3f, 1e-5f);
  // The top endpoint survives the pinch exactly (bounded_output).
  HS_EXPECT_NEAR(PinchModifier(&tension).modify(1.0f), 1.0f, 1e-5f);

  // Fold's negative-reduction guard: a negative coordinate must stay in [0,1]
  // and match an independent triangle-wave reference (abs of the raw fmod).
  auto tri_ref = [](float x) { return fabsf(1.0f - fabsf(fmodf(x, 2.0f))); };
  HS_EXPECT_NEAR(fold.modify(-0.25f), tri_ref(-0.5f), 1e-5f);
  HS_EXPECT_NEAR(fold.modify(-0.9f), tri_ref(-1.8f), 1e-5f);
  HS_EXPECT_GE(fold.modify(-0.25f), 0.0f);
  HS_EXPECT_LE(fold.modify(-0.25f), 1.0f);

  // Negative phase driver also reduces into range.
  float neg_phase = -0.5f;
  FoldModifier fold_np(2.0f, &neg_phase);
  HS_EXPECT_NEAR(fold_np.modify(0.0f), tri_ref(-0.5f), 1e-5f);
  HS_EXPECT_GE(fold_np.modify(-0.3f), 0.0f);
  HS_EXPECT_LE(fold_np.modify(-0.3f), 1.0f);

  // Non-zero Ripple: sin distorts the coordinate off its input.
  float rphase_nz = 0.0f;
  RippleModifier ripple_nz(&rphase_nz, 1.0f, 0.1f);
  HS_EXPECT_NEAR(ripple_nz.modify(0.25f),
                 0.25f + math::fast_sinf(0.25f * math::PI_F * 2.0f) * 0.1f,
                 1e-5f);

  // Non-zero Breathe: quarter-turn phase shifts by amplitude.
  float bphase = math::PI_F * 0.5f;
  HS_EXPECT_NEAR(BreatheModifier(&bphase, 0.1f).modify(0.5f), 0.6f, 1e-3f);

  // Pinch with a negative coordinate re-anchors to t's own integer cell.
  float tension_p = 0.5f;
  HS_EXPECT_NEAR(PinchModifier(&tension_p).modify(-1e-8f), -1e-8f, 1e-7f);
  float pinched_neg = PinchModifier(&tension_p).modify(-0.75f);
  HS_EXPECT_GE(pinched_neg, -1.0f);
  HS_EXPECT_LE(pinched_neg, 0.0f);
}

/**
 * @brief Verifies NoiseWarpModifier's displacement contract.
 * @details The warp must match the value-noise reference, stay within
 *          +/-amplitude, evolve with the time driver, and vanish at zero
 *          amplitude.
 */
inline void test_noise_warp_modifier() {
  float time = 1.7f;
  NoiseWarpModifier warp(&time, 3.0f, 0.1f, 5u);

  for (int i = 0; i <= 20; ++i) {
    float t = i * 0.05f;
    float expected =
        t + (math::value_noise_2d(t * 3.0f, time, 5u) - 0.5f) * 2.0f * 0.1f;
    HS_EXPECT_NEAR(warp.modify(t), expected, 1e-6f);
    HS_EXPECT_LE(std::fabs(warp.modify(t) - t), 0.1f + 1e-6f);
  }

  // The warp field evolves as the time driver advances.
  float before = warp.modify(0.3f);
  time = 4.9f;
  HS_EXPECT_TRUE(warp.modify(0.3f) != before);

  NoiseWarpModifier still(&time, 3.0f, 0.0f);
  HS_EXPECT_NEAR(still.modify(0.42f), 0.42f, 1e-6f);
}

/**
 * @brief Verifies DriftModifier's per-frame memoized noise-walk offset.
 * @details Within one frame the offset is a t-independent constant matching the
 *          value-noise reference; advancing the time driver refreshes it.
 */
inline void test_drift_modifier() {
  float time = 2.3f;
  DriftModifier drift(&time, 0.5f, 0.2f, 11u);

  float expected =
      (math::value_noise_1d(time * 0.5f, 11u) - 0.5f) * 2.0f * 0.2f;
  HS_EXPECT_NEAR(drift.modify(0.0f), expected, 1e-6f);
  HS_EXPECT_LE(std::fabs(drift.modify(0.0f)), 0.2f);

  // Same frame: every coordinate shifts by the identical offset.
  HS_EXPECT_NEAR(drift.modify(0.7f) - 0.7f, drift.modify(0.1f) - 0.1f, 1e-6f);

  // Amplitude applies outside the memo: it takes effect without time moving.
  drift.amplitude = 0.4f;
  HS_EXPECT_NEAR(drift.modify(0.0f), expected * 2.0f, 1e-6f);
  drift.amplitude = 0.2f;

  // New frame: the memo refreshes to the new walk position.
  time = 9.8f;
  float refreshed =
      (math::value_noise_1d(time * 0.5f, 11u) - 0.5f) * 2.0f * 0.2f;
  HS_EXPECT_NEAR(drift.modify(0.4f), 0.4f + refreshed, 1e-6f);
  HS_EXPECT_TRUE(refreshed != expected);
}

/**
 * @brief Largest inter-channel split a gray may carry out of an OKLab shade.
 * @details A gray has a = b = 0, so only the OKLab -> linear-RGB rows and the
 * float -> uint16 rounding separate its channels. Measured worst over the
 * 16-bit gray ramp is 1 LSB, IEEE and -ffast-math alike.
 */
inline constexpr float ACHROMATIC_TOL = 4.0f;

/**
 * @brief Per-channel budget of an identity pass through OKLab.
 * @details linear RGB -> fast_cbrt LMS -> OKLab and back is not exact; with the
 * rotation angle or the chroma scale left at identity, that round trip plus the
 * float -> uint16 rounding is the whole residue. Measured worst over the RGB
 * cube is 7 LSB, IEEE and -ffast-math alike.
 */
inline constexpr float OKLAB_ROUND_TRIP_TOL = 16.0f;

/**
 * @brief Agreement budget between HueSpinShade's folded matrix and hue_rotate.
 * @details Both apply the same OKLab rotation — the spin folds it into a
 * cbrt-LMS 3x3 and cube-roots through fast_cbrt3's shared divide, hue_rotate
 * applies it in OKLab — so they differ only by float reassociation, amplified
 * on the way out through the LMS cube. Measured over the 17^3 channel grid at
 * 64 rotation amounts spanning [-1, 1] turns, IEEE and
 * -ffast-math alike: mean 0.083 LSB, worst single channel 171 LSB on the
 * cbrt-steep saturated corner (40959, 65535, 8191) at 0.206 turns.
 */
inline constexpr float HUE_SPIN_MEAN_HEADROOM = 1.5f;
inline constexpr float HUE_SPIN_WORST_HEADROOM = 1.5f;
inline constexpr float MEASURED_HUE_SPIN_MEAN_LSB = 0.083f;
inline constexpr float MEASURED_HUE_SPIN_WORST_LSB = 171.0f;
inline constexpr float HUE_SPIN_FOLD_MEAN_TOL =
    MEASURED_HUE_SPIN_MEAN_LSB * HUE_SPIN_MEAN_HEADROOM;
inline constexpr float HUE_SPIN_FOLD_TOL =
    MEASURED_HUE_SPIN_WORST_LSB * HUE_SPIN_WORST_HEADROOM;

/**
 * @brief Verifies HueSpinShade against the hue_rotate reference.
 * @details The folded-matrix path must agree with hue_rotate(c, amount) across
 *          the channel grid and every rotation amount, leave grays achromatic,
 *          preserve alpha, and refresh its memo when the driver moves.
 */
inline void test_hue_spin_shade() {
  // Whole-grid agreement: gate the mean and the worst channel.
  {
    double sum = 0.0;
    long long n = 0;
    float worst = 0.0f;
    for (int ai = 0; ai < 64; ++ai) {
      float sweep_amount = -1.0f + ai * (2.0f / 63.0f);
      HueSpinShade sweep(&sweep_amount);
      for (int r = 0; r <= 16; ++r)
        for (int g = 0; g <= 16; ++g)
          for (int b = 0; b <= 16; ++b) {
            Color4 c(Pixel(static_cast<uint16_t>(r * 65535 / 16),
                           static_cast<uint16_t>(g * 65535 / 16),
                           static_cast<uint16_t>(b * 65535 / 16)),
                     1.0f);
            Color4 got = sweep.shade(c, 0.5f);
            Color4 want = hue_rotate(c, sweep_amount);
            const float d =
                std::max({std::fabs(static_cast<float>(got.color.r) -
                                    static_cast<float>(want.color.r)),
                          std::fabs(static_cast<float>(got.color.g) -
                                    static_cast<float>(want.color.g)),
                          std::fabs(static_cast<float>(got.color.b) -
                                    static_cast<float>(want.color.b))});
            sum += d;
            worst = std::max(worst, d);
            ++n;
          }
    }
    const float mean = static_cast<float>(sum / static_cast<double>(n));
    std::printf(
        "  [hue spin fold] mean=%.4f/%.4f worst=%.1f/%.1f LSB\n",
        static_cast<double>(mean), static_cast<double>(HUE_SPIN_FOLD_MEAN_TOL),
        static_cast<double>(worst), static_cast<double>(HUE_SPIN_FOLD_TOL));
    HS_EXPECT_LE(mean, HUE_SPIN_FOLD_MEAN_TOL);
    HS_EXPECT_LE(worst, HUE_SPIN_FOLD_TOL);
  }

  Color4 vivid(Pixel(52000, 9000, 3000), 0.8f);

  float amount = 0.25f;
  HueSpinShade spin(&amount);
  Color4 spun = spin.shade(vivid, 0.5f);
  Color4 ref = hue_rotate(vivid, amount);
  HS_EXPECT_NEAR(static_cast<float>(spun.color.r),
                 static_cast<float>(ref.color.r), HUE_SPIN_FOLD_TOL);
  HS_EXPECT_NEAR(static_cast<float>(spun.color.g),
                 static_cast<float>(ref.color.g), HUE_SPIN_FOLD_TOL);
  HS_EXPECT_NEAR(static_cast<float>(spun.color.b),
                 static_cast<float>(ref.color.b), HUE_SPIN_FOLD_TOL);
  HS_EXPECT_NEAR(spun.alpha, 0.8f, 1e-6f);

  // A quarter turn visibly moves the color off its input.
  HS_EXPECT_GT(std::abs(static_cast<int>(spun.color.r) -
                        static_cast<int>(vivid.color.r)),
               2000);

  // Gray has no chroma to rotate.
  Color4 gray(Pixel(30000, 30000, 30000), 1.0f);
  Color4 spun_gray = spin.shade(gray, 0.0f);
  HS_EXPECT_NEAR(static_cast<float>(spun_gray.color.r),
                 static_cast<float>(spun_gray.color.g), ACHROMATIC_TOL);
  HS_EXPECT_NEAR(static_cast<float>(spun_gray.color.g),
                 static_cast<float>(spun_gray.color.b), ACHROMATIC_TOL);

  // Moving the driver refreshes the memoized matrix.
  amount = 0.5f;
  Color4 half = spin.shade(vivid, 0.5f);
  HS_EXPECT_GT(
      std::abs(static_cast<int>(half.color.r) - static_cast<int>(spun.color.r)),
      2000);
}

/**
 * @brief Verifies HueWobbleShade's t-dependent rotation contract.
 * @details Each sample must match hue_rotate at the wobble's local angle;
 *          opposite sine lobes rotate in opposite directions; zero depth is a
 *          pass-through of the rotation's near-identity.
 */
inline void test_hue_wobble_shade() {
  Color4 vivid(Pixel(52000, 9000, 3000), 0.6f);
  float phase = 0.3f;
  HueWobbleShade wobble(&phase, 1.0f, 0.2f);

  for (float t : {0.0f, 0.25f, 0.6f}) {
    Color4 got = wobble.shade(vivid, t);
    Color4 ref = hue_rotate(
        vivid, 0.2f * math::fast_sinf(t * math::PI_F * 2.0f + phase));
    HS_EXPECT_EQ(got.color.r, ref.color.r);
    HS_EXPECT_EQ(got.color.g, ref.color.g);
    HS_EXPECT_EQ(got.color.b, ref.color.b);
    HS_EXPECT_NEAR(got.alpha, 0.6f, 1e-6f);
  }

  // Opposite sine lobes (t=0.25 vs t=0.75 at phase 0) spin opposite ways.
  float phase0 = 0.0f;
  HueWobbleShade sym(&phase0, 1.0f, 0.2f);
  Color4 up = sym.shade(vivid, 0.25f);
  Color4 down = sym.shade(vivid, 0.75f);
  HS_EXPECT_GT(
      std::abs(static_cast<int>(up.color.g) - static_cast<int>(down.color.g)),
      2000);

  // Zero depth rotates by 0 turns, leaving only the OKLab round trip.
  HueWobbleShade flat(&phase0, 1.0f, 0.0f);
  Color4 same = flat.shade(vivid, 0.4f);
  HS_EXPECT_NEAR(static_cast<float>(same.color.r),
                 static_cast<float>(vivid.color.r), OKLAB_ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(static_cast<float>(same.color.g),
                 static_cast<float>(vivid.color.g), OKLAB_ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(static_cast<float>(same.color.b),
                 static_cast<float>(vivid.color.b), OKLAB_ROUND_TRIP_TOL);
}

/**
 * @brief Verifies SparkleShade whitens only over-threshold noise sites.
 * @details Sub-threshold samples pass through bit-exact; over-threshold samples
 *          match the reference lerp toward white and never darken; alpha is
 *          untouched. A 0.5 threshold guarantees both cases occur in a domain
 *          sweep.
 */
inline void test_sparkle_shade() {
  Color4 base(Pixel(20000, 8000, 30000), 0.9f);
  float time = 3.1f;
  SparkleShade sparkle(&time, 16.0f, 0.5f, 21u);

  int lit = 0, dark = 0;
  for (int i = 0; i <= 100; ++i) {
    float t = i * 0.01f;
    float n = math::value_noise_2d(t * 16.0f, time, 21u);
    Color4 got = sparkle.shade(base, t);
    HS_EXPECT_NEAR(got.alpha, 0.9f, 1e-6f);
    if (n <= 0.5f) {
      dark++;
      HS_EXPECT_EQ(got.color.r, base.color.r);
      HS_EXPECT_EQ(got.color.g, base.color.g);
      HS_EXPECT_EQ(got.color.b, base.color.b);
    } else {
      lit++;
      float w = (n - 0.5f) / 0.5f;
      Pixel ref = base.color.lerp16(Pixel(65535, 65535, 65535), frac_to_q16(w));
      HS_EXPECT_EQ(got.color.r, ref.r);
      HS_EXPECT_EQ(got.color.g, ref.g);
      HS_EXPECT_EQ(got.color.b, ref.b);
      HS_EXPECT_GE(got.color.r, base.color.r);
      HS_EXPECT_GE(got.color.g, base.color.g);
      HS_EXPECT_GE(got.color.b, base.color.b);
    }
  }
  HS_EXPECT_GT(lit, 0);
  HS_EXPECT_GT(dark, 0);
}

/**
 * @brief Verifies ChromaPulseShade scales chroma while holding L, hue, and alpha.
 * @details Positive sine lobes boost measured OKLCH chroma, negative lobes cut
 *          it, a zero lobe is a near-identity, and grays stay achromatic.
 */
inline void test_chroma_pulse_shade() {
  Color4 mid(Pixel(30000, 12000, 6000), 0.7f);
  OKLCH before = pixel_to_oklch(mid.color);

  // sin(pi/2) = 1: chroma scales up by 1 + depth.
  float phase = math::PI_F * 0.5f;
  ChromaPulseShade pulse(&phase, 0.3f);
  Color4 boosted = pulse.shade(mid, 0.2f);
  OKLCH after = pixel_to_oklch(boosted.color);
  HS_EXPECT_GT(after.C, before.C * 1.1f);
  HS_EXPECT_NEAR(wrap_hue_delta(after.h - before.h), 0.0f, .01f);
  HS_EXPECT_NEAR(after.L, before.L, 0.02f);
  HS_EXPECT_NEAR(boosted.alpha, 0.7f, 1e-6f);

  // sin(-pi/2) = -1: chroma scales down toward gray.
  float neg_phase = -math::PI_F * 0.5f;
  ChromaPulseShade cut(&neg_phase, 0.3f);
  const OKLCH CUT = pixel_to_oklch(cut.shade(mid, 0.2f).color);
  HS_EXPECT_LT(CUT.C, before.C * 0.9f);
  HS_EXPECT_NEAR(wrap_hue_delta(CUT.h - before.h), 0.0f, .01f);

  // sin(0) = 0: pass-through within the OKLab round-trip budget.
  float zero = 0.0f;
  ChromaPulseShade flat(&zero, 0.3f);
  Color4 same = flat.shade(mid, 0.2f);
  // ChromaPulseShade refreshes its memo when phase changes.
  zero = math::PI_F * .5f;
  HS_EXPECT_EQ(flat.shade(mid, .2f).color,
               ChromaPulseShade(&zero, .3f).shade(mid, .2f).color);
  HS_EXPECT_NEAR(static_cast<float>(same.color.r),
                 static_cast<float>(mid.color.r), OKLAB_ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(static_cast<float>(same.color.g),
                 static_cast<float>(mid.color.g), OKLAB_ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(static_cast<float>(same.color.b),
                 static_cast<float>(mid.color.b), OKLAB_ROUND_TRIP_TOL);

  // Gray has no chroma to scale.
  Color4 gray(Pixel(20000, 20000, 20000), 1.0f);
  Color4 pulsed_gray = pulse.shade(gray, 0.0f);
  HS_EXPECT_NEAR(static_cast<float>(pulsed_gray.color.r),
                 static_cast<float>(pulsed_gray.color.g), ACHROMATIC_TOL);
  HS_EXPECT_NEAR(static_cast<float>(pulsed_gray.color.g),
                 static_cast<float>(pulsed_gray.color.b), ACHROMATIC_TOL);
}

/**
 * @brief Verifies LightnessGrainShade's uniform-gain brightness grain.
 * @details Each sample must match the reference uniform gain, keep channel
 *          ratios (hue) within rounding, leave alpha alone, and be an identity
 *          at zero amplitude.
 */
inline void test_lightness_grain_shade() {
  Color4 base(Pixel(40000, 16000, 8000), 0.5f);
  float time = 6.4f;
  LightnessGrainShade grain(&time, 8.0f, 0.25f, 13u);

  for (int i = 0; i <= 20; ++i) {
    float t = i * 0.05f;
    float n = math::value_noise_2d(t * 8.0f, time, 13u);
    float gain = 1.0f + 0.25f * (2.0f * n - 1.0f);
    Color4 got = grain.shade(base, t);
    Pixel ref = base.color * gain;
    HS_EXPECT_EQ(got.color.r, ref.r);
    HS_EXPECT_EQ(got.color.g, ref.g);
    HS_EXPECT_EQ(got.color.b, ref.b);
    HS_EXPECT_NEAR(got.alpha, 0.5f, 1e-6f);
    // Uniform gain preserves the channel ratio (hue) up to rounding.
    HS_EXPECT_NEAR(static_cast<float>(got.color.r) / got.color.g,
                   static_cast<float>(base.color.r) / base.color.g, 0.01f);
  }

  LightnessGrainShade still(&time, 8.0f, 0.0f);
  Color4 same = still.shade(base, 0.3f);
  HS_EXPECT_EQ(same.color.r, base.color.r);
  HS_EXPECT_EQ(same.color.g, base.color.g);
  HS_EXPECT_EQ(same.color.b, base.color.b);
}

/**
 * @brief Verifies IridescentShade's additive cosine overlay.
 * @details On black the output must equal the weighted sheen reference with the
 *          per-channel thirds phase offsets; the add saturates at white; zero
 *          weight and alpha are pass-throughs.
 */
inline void test_iridescent_shade() {
  float phase = 0.7f;
  IridescentShade sheen(&phase, 2.0f, 0.4f);

  Color4 black(Pixel(0, 0, 0), 0.3f);
  for (float t : {0.0f, 0.3f, 0.85f}) {
    float arg = t * 2.0f * math::PI_F * 2.0f + phase;
    constexpr float THIRD = 2.0f * math::PI_F / 3.0f;
    Pixel ref =
        Pixel(srgb_to_linear_interp(0.5f + 0.5f * math::fast_cosf(arg)),
              srgb_to_linear_interp(0.5f + 0.5f * math::fast_cosf(arg + THIRD)),
              srgb_to_linear_interp(
                  0.5f + 0.5f * math::fast_cosf(arg + 2.0f * THIRD))) *
        0.4f;
    Color4 got = sheen.shade(black, t);
    HS_EXPECT_EQ(got.color.r, ref.r);
    HS_EXPECT_EQ(got.color.g, ref.g);
    HS_EXPECT_EQ(got.color.b, ref.b);
    HS_EXPECT_NEAR(got.alpha, 0.3f, 1e-6f);
  }

  // The overlay saturates instead of wrapping on a near-white input.
  Color4 bright(Pixel(65000, 65000, 65000), 1.0f);
  Color4 sat = sheen.shade(bright, 0.1f);
  HS_EXPECT_GE(sat.color.r, bright.color.r);
  HS_EXPECT_GE(sat.color.g, bright.color.g);
  HS_EXPECT_GE(sat.color.b, bright.color.b);

  // Zero weight adds nothing.
  IridescentShade off(&phase, 2.0f, 0.0f);
  Color4 mid(Pixel(12345, 23456, 34567), 1.0f);
  Color4 same = off.shade(mid, 0.5f);
  HS_EXPECT_EQ(same.color.r, mid.color.r);
  HS_EXPECT_EQ(same.color.g, mid.color.g);
  HS_EXPECT_EQ(same.color.b, mid.color.b);
}

/**
 * @brief Verifies StaticPalette folds its modifier chain then queries the source.
 * @details Applies the modifier chain in order before sampling the source.
 */
inline void test_static_palette_composition() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // Single modifier: scale=2 turns get(0.25) into source.get(0.5).
  ScaleModifier scale(2.0f);
  StaticPalette<Gradient, Coords<ScaleModifier>> sp;
  sp.bind(&grad, &scale);
  HS_EXPECT_EQ(sp.get(0.25f).color.r, grad.get(0.5f).color.r);

  // Out-of-range coord modifiers are flagged; bounded ones are not.
  static_assert(coord_requires_wrap<ScaleModifier>());
  static_assert(coord_requires_wrap<CycleModifier>());
  static_assert(coord_requires_wrap<BreatheModifier>());
  static_assert(coord_requires_wrap<RippleModifier>());
  static_assert(coord_requires_wrap<NoiseWarpModifier>());
  static_assert(coord_requires_wrap<DriftModifier>());
  static_assert(!coord_requires_wrap<InsetModifier>());
  static_assert(!coord_requires_wrap<MirrorModifier>());
  static_assert(!coord_requires_wrap<QuantizeModifier>());

  // Modifiers that confine [0,1] to [0,1] and reach exactly 1.0 carry
  // bounded_output, so wrap_t would fold their top endpoint to 0.
  static_assert(coord_bounded_output<FoldModifier>());
  static_assert(coord_bounded_output<ReverseModifier>());
  static_assert(coord_bounded_output<MirrorModifier>());
  static_assert(coord_bounded_output<InsetModifier>());
  static_assert(coord_bounded_output<PinchModifier>());
  static_assert(coord_bounded_output<QuantizeModifier>());
  static_assert(!coord_bounded_output<ScaleModifier>());
  static_assert(!coord_bounded_output<CycleModifier>());
  static_assert(!coord_bounded_output<BreatheModifier>());
  static_assert(!coord_bounded_output<RippleModifier>());
  static_assert(!coord_bounded_output<NoiseWarpModifier>());
  static_assert(!coord_bounded_output<DriftModifier>());

  // Fold's triangle wave and Inset's clamp confine arbitrary input, so they can
  // absorb an out-of-range predecessor. Reverse, Mirror and Quantize are
  // bounded only on [0,1] and pass an out-of-range coordinate through.
  static_assert(coord_rebounds_input<FoldModifier>());
  static_assert(coord_rebounds_input<InsetModifier>());
  static_assert(coord_rebounds_input<WrapModifier>());
  // WrapModifier's fold stops short of 1.0, so it is not a bounded tail.
  static_assert(!coord_bounded_output<WrapModifier>());
  static_assert(!coord_rebounds_input<ReverseModifier>());
  static_assert(!coord_rebounds_input<MirrorModifier>());
  static_assert(!coord_rebounds_input<QuantizeModifier>());
  static_assert(!coord_rebounds_input<ScaleModifier>());
  static_assert(!coord_rebounds_input<CycleModifier>());

  // Only a re-bounding entry clears an unbounded predecessor, and an unbounded
  // tail undoes it again.
  static_assert(!coord_chain_leaves_unit<>());
  static_assert(coord_chain_leaves_unit<ScaleModifier>());
  static_assert(!coord_chain_leaves_unit<ScaleModifier, FoldModifier>());
  static_assert(coord_chain_leaves_unit<ScaleModifier, ReverseModifier>());
  static_assert(coord_chain_leaves_unit<CycleModifier, MirrorModifier>());
  static_assert(
      !coord_chain_leaves_unit<CycleModifier, WrapModifier, MirrorModifier>());
  static_assert(coord_chain_leaves_unit<InsetModifier, CycleModifier>());
  static_assert(!coord_chain_bounded_tail<>());
  static_assert(!coord_chain_bounded_tail<MirrorModifier, CycleModifier>());
  static_assert(coord_chain_bounded_tail<CycleModifier, MirrorModifier>());

  // MirrorModifier's designed peak lands on 1.0; Wrap=false keeps it at the
  // source's last stop instead of folding it to the first.
  MirrorModifier mirror;
  StaticPalette<Gradient, Coords<MirrorModifier>, Colors<>, /*Wrap=*/false>
      mirrored;
  mirrored.bind(&grad, &mirror);
  HS_EXPECT_EQ(mirrored.get(0.5f).color.r, grad.get(1.0f).color.r);

  // Quantize's top band snaps to 1.0; Wrap=false renders it as the source's
  // last stop, distinct from the band at 0.
  QuantizeModifier quant(4.0f);
  StaticPalette<Gradient, Coords<QuantizeModifier>, Colors<>, /*Wrap=*/false>
      quantized;
  quantized.bind(&grad, &quant);
  HS_EXPECT_EQ(quantized.get(0.95f).color.r, grad.get(1.0f).color.r);
  HS_EXPECT_EQ(quantized.get(0.05f).color.r, grad.get(0.0f).color.r);

  // Two modifiers apply in tuple order (scale THEN cycle): 0.2 -> 0.4 -> 0.5.
  float off = 0.1f;
  CycleModifier cycle(&off);
  StaticPalette<Gradient, Coords<ScaleModifier, CycleModifier>> sp2;
  sp2.bind(&grad, &scale, &cycle);
  HS_EXPECT_EQ(sp2.get(0.2f).color.r, grad.get(0.5f).color.r);

  // Scroll then mirror: a WrapModifier between the two absorbs the cycle's
  // out-of-range coordinate, so Wrap=false can keep the mirror's 1.0 peak.
  float quarter = 0.25f;
  CycleModifier quarter_cycle(&quarter);
  WrapModifier fold_unit;
  HS_EXPECT_NEAR(fold_unit.modify(-0.25f), 0.75f, 1e-6f);
  HS_EXPECT_NEAR(fold_unit.modify(1.25f), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(fold_unit.modify(1.0f), 0.0f, 1e-6f);
  StaticPalette<Gradient, Coords<CycleModifier, WrapModifier, MirrorModifier>,
                Colors<>, /*Wrap=*/false>
      scrolled_mirror;
  scrolled_mirror.bind(&grad, &quarter_cycle, &fold_unit, &mirror);
  // 0.9375 -> 1.1875 -> 0.1875 -> 0.375.
  HS_EXPECT_EQ(scrolled_mirror.get(0.9375f).color.r, grad.get(0.375f).color.r);
  // 0.25 -> 0.5 -> 0.5 -> the mirror's peak, which reaches the last stop.
  HS_EXPECT_EQ(scrolled_mirror.get(0.25f).color.r, grad.get(1.0f).color.r);

  EdgeAlphaShade edge_alpha;
  StaticPalette<Gradient, Coords<CycleModifier>, Colors<EdgeAlphaShade>>
      cycled_alpha;
  cycled_alpha.bind(&grad, &cycle, &edge_alpha);
  HS_EXPECT_GT(cycled_alpha.get(1.2f).alpha, 0.9f);
}

/**
 * @brief Constant 0.5 falloff function for AlphaFalloffShade.
 * @param Unused normalized coordinate (the falloff is constant).
 * @return Constant falloff factor 0.5.
 * @details Must be a plain function pointer to bind to AlphaFalloffShade.
 */
inline float half_falloff(float) { return 0.5f; }

/**
 * @brief Verifies each modifier composed via StaticPalette produces its remap.
 * @details Exercises Reverse, Mirror, Inset+EdgeFade, Inset+EdgeAlpha,
 *          AlphaFalloff, PaletteFacade, and SolidColor. Bounded remaps use
 *          Wrap=false so a coordinate landing exactly on 1.0 reaches the source's
 *          last stop rather than wrapping to 0.
 */
inline void test_palette_wrappers() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // ReverseModifier: t -> 1-t, so t=0 samples the white end and t=1 the black.
  ReverseModifier rev;
  StaticPalette<Gradient, Coords<ReverseModifier>, Colors<>, /*Wrap=*/false>
      revp;
  revp.bind(&grad, &rev);
  HS_EXPECT_GT(revp.get(0.0f).color.r, 60000); // source's t=1 (white)
  HS_EXPECT_EQ(revp.get(1.0f).color.r, 0);     // source's t=0 (black)

  // MirrorModifier: [0,1] -> [0,1,0]; the midpoint reaches the far end.
  MirrorModifier mir;
  StaticPalette<Gradient, Coords<MirrorModifier>, Colors<>, /*Wrap=*/false>
      mirp;
  mirp.bind(&grad, &mir);
  HS_EXPECT_EQ(mirp.get(0.0f).color.r, 0);
  HS_EXPECT_GT(mirp.get(0.5f).color.r, 60000);
  HS_EXPECT_EQ(mirp.get(1.0f).color.r, 0);

  // Opaque vignette = InsetModifier + EdgeFadeShade: fades to black at the
  // edges, source color in the middle band.
  InsetModifier inset;
  EdgeFadeShade fade;
  StaticPalette<Gradient, Coords<InsetModifier>, Colors<EdgeFadeShade>,
                /*Wrap=*/false>
      vig;
  vig.bind(&grad, &inset, &fade);
  HS_EXPECT_LT(vig.get(0.0f).color.r, 1000); // edge -> ~black
  uint16_t vmid = vig.get(0.5f).color.r;     // middle -> source.get(0.5)
  HS_EXPECT_GT(vmid, 1000);
  HS_EXPECT_LT(vmid, 64000);

  // Transparent vignette = InsetModifier + EdgeAlphaShade: alpha (not color)
  // fades at the edges.
  EdgeAlphaShade alpha_fade;
  StaticPalette<Gradient, Coords<InsetModifier>, Colors<EdgeAlphaShade>,
                /*Wrap=*/false>
      tv;
  tv.bind(&grad, &inset, &alpha_fade);
  HS_EXPECT_NEAR(tv.get(0.0f).alpha, 0.0f, 1e-3f); // edge -> transparent
  HS_EXPECT_NEAR(tv.get(0.1f).alpha, 0.5f, 1e-2f); // quintic(0.5) = 0.5
  HS_EXPECT_NEAR(tv.get(0.5f).alpha, 1.0f, 1e-3f); // middle -> opaque

  // AlphaFalloffShade scales alpha by the falloff function.
  AlphaFalloffShade afs(half_falloff);
  StaticPalette<Gradient, Coords<>, Colors<AlphaFalloffShade>, /*Wrap=*/false>
      afp;
  afp.bind(&grad, &afs);
  HS_EXPECT_NEAR(afp.get(0.5f).alpha, 0.5f, 1e-5f);
  HS_EXPECT_EQ(afp.get(0.5f).color.r, grad.get(0.5f).color.r);

  // PaletteFacade exposes a composition through the polymorphic Palette API.
  PaletteFacade<decltype(afp)> facade(&afp);
  const Palette &as_palette = facade;
  HS_EXPECT_NEAR(as_palette.get(0.5f).alpha, 0.5f, 1e-5f);

  // Solid color ignores t.
  SolidColorPalette solid(Color4(Pixel(111, 222, 333), 0.7f));
  HS_EXPECT_EQ(solid.get(0.9f).color.g, 222);
  HS_EXPECT_NEAR(solid.get(0.1f).alpha, 0.7f, 1e-6f);
}

/**
 * @brief Identity falloff function: reports the coordinate it was shaded with.
 * @param t Normalized coordinate.
 * @return t unchanged.
 * @details Must be a plain function pointer to bind to AlphaFalloffShade.
 */
inline float identity_falloff(float t) { return t; }

/**
 * @brief Verifies ShadeCoord picks the color chain's coordinate independently
 *        of Wrap.
 * @details AlphaFalloffShade reports the coordinate it was handed through
 *          alpha. MATCH_WRAP selects the source lookup coordinate under
 *          Wrap=true and the raw input under Wrap=false. LOOKUP and RAW_INPUT
 *          pin their own coordinate either way.
 */
inline void test_palette_shade_coord_policy() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  AlphaFalloffShade report(identity_falloff);

  StaticPalette<Gradient, Coords<>, Colors<AlphaFalloffShade>> matched;
  matched.bind(&grad, &report);
  HS_EXPECT_NEAR(matched.get(1.25f).alpha, 0.25f, 1e-5f);

  StaticPalette<Gradient, Coords<>, Colors<AlphaFalloffShade>, /*Wrap=*/true,
                ShadeCoord::RAW_INPUT>
      raw;
  raw.bind(&grad, &report);
  HS_EXPECT_NEAR(raw.get(1.25f).alpha, 1.25f, 1e-5f);

  ReverseModifier rev;
  StaticPalette<Gradient, Coords<ReverseModifier>, Colors<AlphaFalloffShade>,
                /*Wrap=*/false>
      unwrapped;
  unwrapped.bind(&grad, &rev, &report);
  HS_EXPECT_NEAR(unwrapped.get(0.25f).alpha, 0.25f, 1e-5f);

  StaticPalette<Gradient, Coords<ReverseModifier>, Colors<AlphaFalloffShade>,
                /*Wrap=*/false, ShadeCoord::LOOKUP>
      unwrapped_lookup;
  unwrapped_lookup.bind(&grad, &rev, &report);
  HS_EXPECT_NEAR(unwrapped_lookup.get(0.25f).alpha, 0.75f, 1e-5f);
}

/**
 * @brief Verifies the reusable noise-hue palette resolves spatial hue shifts
 *        from its prepared LUTs and preserves the source palette's alpha.
 */
inline void test_noise_hue_palette() {
  static std::array<Pixel, HueRotationLutView::SIZE> hue_rotation;
  static std::array<int8_t, HueNoiseLutView::SIZE> hue_noise;
  Gradient source{{0.0f, CPixel(255u, 40u, 10u)},
                  {1.0f, CPixel(10u, 80u, 255u)}};
  prepare_hue_rotation_lut(
      std::span<Pixel, HueRotationLutView::SIZE>(hue_rotation), source);
  FastNoiseLite noise;
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(6047);
  noise.SetFrequency(1.0f);
  prepare_hue_noise_lut(std::span<int8_t, HueNoiseLutView::SIZE>(hue_noise),
                        noise, 2.0f, 0.0f);

  NoiseHuePalette<Gradient> palette(&source, hue_rotation.data(),
                                    hue_noise.data());
  const float value = 0.375f;
  const float shift =
      palette.hue_shift(math::Vector(1.0f, 2.0f, 3.0f).normalized(), 0.4f);
  const Color4 actual = palette.get(value, shift);
  const Pixel expected =
      sample_hue_rotation_lut({hue_rotation.data(), true}, value, shift);
  HS_EXPECT_EQ(actual.color.r, expected.r);
  HS_EXPECT_EQ(actual.color.g, expected.g);
  HS_EXPECT_EQ(actual.color.b, expected.b);
  HS_EXPECT_NEAR(actual.alpha, source.get(value).alpha, 1e-6f);
  HS_EXPECT_GT(fabsf(palette.hue_shift(math::X_AXIS, 1.0f) -
                     palette.hue_shift(math::Y_AXIS, 1.0f)),
               1e-3f);
  const float uv_a = palette.noise_uv(1.0f, 0.0f, 1.0f, 0.0f);
  const float uv_b = palette.noise_uv(1.0f, 0.0f, 0.0f, 1.0f);
  HS_EXPECT_GT(fabsf(uv_a - uv_b), 1e-3f);
  HS_EXPECT_NEAR(uv_a, palette.noise_uv(1.0f, -0.0f, 1.0f, -0.0f), 1e-6f);
}

inline void test_noise_hue_palette_direct() {
  std::array<int8_t, HueNoiseLutView::SIZE> noise;
  noise.fill(64);
  SolidColorPalette source(Color4(Pixel(40000, 12000, 3000), 0.37f));
  NoiseHuePalette<SolidColorPalette> palette(&source, noise.data());
  for (float shift : {-0.4f, 0.0f, 0.25f, 1.5f}) {
    const Color4 original = source.get(0.3f);
    const Color4 expected =
        shift == 0.0f ? original : hue_rotate_lut_gamut(original, shift);
    const Color4 actual = palette.get(0.3f, shift);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_EQ(actual.alpha, original.alpha);
  }
  HS_EXPECT_NEAR(palette.hue_shift(math::X_AXIS, 0.5f), 32.0f / 127.0f, 1e-6f);
  HS_EXPECT_EQ(palette.get(0.3f, math::X_AXIS, 0.5f).color,
               palette.get(0.3f, 32.0f / 127.0f).color);

  static std::array<Pixel, HueRotationLutView::SIZE> rotation;
  palette.bind(&source, rotation.data(), noise.data());
  palette.bind(&source, noise.data());
  HS_EXPECT_EQ(palette.get(0.3f, 0.0f).color, source.get(0.3f).color);
}

inline void test_noise_shimmer_palette() {
  std::array<int8_t, HueNoiseLutView::SIZE> noise;
  for (Pixel base : {Pixel(4000, 8000, 12000), Pixel(65535, 0, 0),
                     Pixel(0, 0, 0), Pixel(65535, 65535, 65535)}) {
    SolidColorPalette source(Color4(base, 0.37f));
    NoiseShimmerPalette<SolidColorPalette> palette(&source, noise.data());
    HS_EXPECT_EQ(palette.get(0.3f, 0.0f).color, base);
    HS_EXPECT_EQ(palette.get(0.3f, -1.0f).color, base);
    const auto to_lab = [](Pixel pixel) {
      const LinRGB rgb = pixel_to_linrgb(pixel);
      return linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
    };
    const OKLab original = to_lab(base);
    for (float lift : {0.1f, 0.4f, 1.0f}) {
      const Color4 raised = palette.get(0.3f, lift);
      const OKLab actual = to_lab(raised.color);
      HS_EXPECT_EQ(raised.alpha, 0.37f);
      HS_EXPECT_NEAR(actual.L, original.L + (1.0f - original.L) * lift, 0.005f);
      if (original.a * original.a + original.b * original.b > .0004f &&
          actual.a * actual.a + actual.b * actual.b > .0004f)
        HS_EXPECT_NEAR(wrap_hue_delta(atan2f(actual.b, actual.a) -
                                      atan2f(original.b, original.a)),
                       0.0f, .03f);
    }
    HS_EXPECT_EQ(palette.get(0.3f, 2.0f).color, palette.get(0.3f, 1.0f).color);
    noise.fill(-127);
    HS_EXPECT_EQ(palette.get(0.3f, math::X_AXIS, 0.5f).color, base);
    noise.fill(127);
    HS_EXPECT_EQ(palette.lightness_shift(math::Y_AXIS, 0.5f), 0.5f);
    HS_EXPECT_EQ(palette.get(0.3f, math::X_AXIS, 0.5f).color,
                 palette.get(0.3f, 0.5f).color);
    NoiseHuePalette<SolidColorPalette> hue(&source, noise.data());
    FastNoiseLite field;
    field.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    field.SetSeed(6047);
    field.SetFrequency(1.0f);
    prepare_hue_noise_lut(std::span<int8_t, HueNoiseLutView::SIZE>(noise),
                          field, 2.0f, 0.0f);
    HS_EXPECT_GT(
        fabsf(palette.noise(math::X_AXIS) - palette.noise(math::Y_AXIS)),
        1e-3f);
    for (const auto &direction : {math::X_AXIS, math::Y_AXIS}) {
      HS_EXPECT_EQ(palette.noise(direction),
                   sample_hue_noise_lut({noise.data(), true}, direction));
      HS_EXPECT_EQ(hue.noise(direction), palette.noise(direction));
    }
  }
}

/**
 * @brief Verifies the shared hue-noise bake cache rebuilds the table exactly
 *        when an input moved, and leaves it untouched otherwise.
 */
template <bool OPENSIMPLEX2_UNIT>
inline void test_hue_noise_bake_cache_variant() {
  static std::array<int8_t, HueNoiseLutView::SIZE> hue_noise;
  FastNoiseLite noise;
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(6047);
  noise.SetFrequency(1.0f);
  // Unreachable through the quantizer, whose range is [-127, 127].
  constexpr int8_t POISON = -128;
  constexpr size_t LAST = HueNoiseLutView::SIZE - 1;

  HueNoiseBakeCache cache;
  HS_EXPECT_EQ(cache.scale, 0.0f);
  HS_EXPECT_EQ(cache.phase, 0.0f);
  HS_EXPECT_TRUE(
      cache.refresh<OPENSIMPLEX2_UNIT>(hue_noise, noise, 2.0f, 0.25f));
  HS_EXPECT_EQ(cache.scale, 2.0f);
  HS_EXPECT_EQ(cache.phase, 0.25f);

  hue_noise[0] = POISON;
  hue_noise[LAST] = POISON;
  HS_EXPECT_TRUE(
      !cache.refresh<OPENSIMPLEX2_UNIT>(hue_noise, noise, 2.0f, 0.25f));
  HS_EXPECT_EQ(static_cast<int>(hue_noise[0]), static_cast<int>(POISON));
  HS_EXPECT_EQ(static_cast<int>(hue_noise[LAST]), static_cast<int>(POISON));

  HS_EXPECT_TRUE(
      cache.refresh<OPENSIMPLEX2_UNIT>(hue_noise, noise, 2.5f, 0.25f));
  HS_EXPECT_NE(static_cast<int>(hue_noise[0]), static_cast<int>(POISON));
  HS_EXPECT_NE(static_cast<int>(hue_noise[LAST]), static_cast<int>(POISON));
  HS_EXPECT_EQ(cache.scale, 2.5f);

  hue_noise[0] = POISON;
  hue_noise[LAST] = POISON;
  HS_EXPECT_TRUE(
      cache.refresh<OPENSIMPLEX2_UNIT>(hue_noise, noise, 2.5f, 0.5f));
  HS_EXPECT_NE(static_cast<int>(hue_noise[0]), static_cast<int>(POISON));
  HS_EXPECT_NE(static_cast<int>(hue_noise[LAST]), static_cast<int>(POISON));
  HS_EXPECT_EQ(cache.phase, 0.5f);

#if defined(HS_TEST_FAST_MATH)
  constexpr int N = HueNoiseLutView::FACE_STEPS;
  constexpr float STEP = 2.0f / (N - 1);
  constexpr float NOISE_TOLERANCE = 2e-5f;
  constexpr float ROUNDING_TOLERANCE = 0.5f + 127.0f * NOISE_TOLERANCE;
  const math::Vector OFFSET = math::noise_sphere_loop_offset(0.5f);
  for (int face = 0; face < HueNoiseLutView::FACE_COUNT; ++face)
    for (int y = 0; y < N; ++y)
      for (int x = 0; x < N; ++x) {
        const math::Vector D =
            hue_noise_face_direction(face, -1.0f + STEP * x, -1.0f + STEP * y);
        const math::Vector Q = 2.5f * D + OFFSET;
        const float SAMPLE =
            hs::clamp(noise.GetNoiseSingle(Q.x, Q.y, Q.z), -1.0f, 1.0f);
        HS_EXPECT_NEAR(static_cast<float>(hue_noise[face * N * N + y * N + x]),
                       SAMPLE * 127.0f, ROUNDING_TOLERANCE);
      }
#else
  static std::array<int8_t, HueNoiseLutView::SIZE> reference;
  prepare_hue_noise_lut<OPENSIMPLEX2_UNIT>(
      std::span<int8_t, HueNoiseLutView::SIZE>(reference), noise, 2.5f, 0.5f);
  HS_EXPECT_EQ(
      std::memcmp(hue_noise.data(), reference.data(), HueNoiseLutView::SIZE),
      0);
#endif
}

inline void test_hue_noise_bake_cache() {
  test_hue_noise_bake_cache_variant<false>();
  test_hue_noise_bake_cache_variant<true>();
}

/** @brief Compares paired cube-face bakes with independent scalar face walks. */
inline void test_hue_noise_paired_bakes_match_reference() {
  static std::array<int8_t, HueNoiseLutView::SIZE> actual, reference;
#if defined(HS_TEST_FAST_MATH)
  static std::array<float, HueNoiseLutView::SIZE> reference_samples;
#endif
  const auto check_bake = [&] {
#if defined(HS_TEST_FAST_MATH)
    // Fast-math reassociation can cross an int8 rounding boundary.
    constexpr float NOISE_TOLERANCE = 2e-5f;
    constexpr float ROUNDING_TOLERANCE = 0.5f + 127.0f * NOISE_TOLERANCE;
    for (size_t i = 0; i < actual.size(); ++i)
      HS_EXPECT_NEAR(static_cast<float>(actual[i]), reference_samples[i],
                     ROUNDING_TOLERANCE);
#else
    HS_EXPECT_EQ(std::memcmp(actual.data(), reference.data(), actual.size()),
                 0);
#endif
  };
  for (int variant = 0; variant < 4; ++variant) {
    FastNoiseLite noise;
    noise.SetSeed(6047 + variant);
    noise.SetNoiseType(variant == 2 ? FastNoiseLite::NoiseType_Perlin
                                    : FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetFrequency(variant == 1 ? 0.7f : 1.0f);
    if (variant == 3)
      noise.SetRotationType3D(FastNoiseLite::RotationType3D_ImproveXYPlanes);
    for (float scale : {0.0625f, 2.0f, 5.0f}) {
      for (float phase : {0.0f, 0.25f, 0.625f}) {
        const math::Vector offset = math::noise_sphere_loop_offset(phase);
        constexpr int N = HueNoiseLutView::FACE_STEPS;
        constexpr float STEP = 2.0f / (N - 1);
        for (int face = 0; face < HueNoiseLutView::FACE_COUNT; ++face) {
          for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
              const math::Vector d = hue_noise_face_direction(
                  face, -1.0f + STEP * x, -1.0f + STEP * y);
              const math::Vector q = scale * d + offset;
              const float sample =
                  hs::clamp(noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
              const int INDEX = face * N * N + y * N + x;
              reference[INDEX] = static_cast<int8_t>(
                  sample * 127.0f + (sample < 0.0f ? -0.5f : 0.5f));
#if defined(HS_TEST_FAST_MATH)
              reference_samples[INDEX] = sample * 127.0f;
#endif
            }
          }
        }
        prepare_hue_noise_lut(actual, noise, scale, phase);
        check_bake();
        if (variant == 0) {
          prepare_hue_noise_lut<true>(actual, noise, scale, phase);
          check_bake();
        }
      }
    }
  }
}

/**
 * @brief Verifies the cube-map hue-noise LUT agrees across a shared face edge
 *        and at a three-face corner, where the face tie-break picks different
 *        faces for neighbouring directions.
 */
inline void test_hue_noise_lut_seamless_across_faces() {
  static std::array<int8_t, HueNoiseLutView::SIZE> hue_noise;
  FastNoiseLite noise;
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(1301);
  noise.SetFrequency(1.0f);
  prepare_hue_noise_lut(std::span<int8_t, HueNoiseLutView::SIZE>(hue_noise),
                        noise, 3.0f, 0.0f);
  const HueNoiseLutView view{hue_noise.data(), true};

  constexpr float NUDGE = 1e-4f;
  float edge_low = 1.0f;
  float edge_high = -1.0f;
  for (int step = 0; step <= 16; ++step) {
    const float t = -1.0f + 0.125f * step;
    const float x_face = sample_hue_noise_lut(
        view, math::Vector(1.0f, 1.0f - NUDGE, t).normalized());
    const float y_face = sample_hue_noise_lut(
        view, math::Vector(1.0f - NUDGE, 1.0f, t).normalized());
    HS_EXPECT_NEAR(x_face, y_face, 2e-3f);
    edge_low = std::min(edge_low, x_face);
    edge_high = std::max(edge_high, x_face);
  }
  HS_EXPECT_GT(edge_high - edge_low, 1e-2f);

  const float corner_x = sample_hue_noise_lut(
      view, math::Vector(1.0f, 1.0f - NUDGE, 1.0f - NUDGE).normalized());
  const float corner_y = sample_hue_noise_lut(
      view, math::Vector(1.0f - NUDGE, 1.0f, 1.0f - NUDGE).normalized());
  const float corner_z = sample_hue_noise_lut(
      view, math::Vector(1.0f - NUDGE, 1.0f - NUDGE, 1.0f).normalized());
  HS_EXPECT_NEAR(corner_x, corner_y, 2e-3f);
  HS_EXPECT_NEAR(corner_x, corner_z, 2e-3f);
}

inline void test_hue_rotation_lut_clamps_out_of_range_value() {
  static std::array<Pixel, 3 * HueRotationLutView::SIZE> storage;
  storage.fill(Pixel(0, 0, 65535));
  Pixel *const lut = storage.data() + HueRotationLutView::SIZE;
  for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
    lut[hue] = Pixel(65535, 0, 0);
    lut[(HueRotationLutView::VALUE_STEPS - 1) * HueRotationLutView::HUE_STEPS +
        hue] = Pixel(0, 65535, 0);
  }

  const HueRotationLutView view{lut, true};
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, -0.5f, 0.25f), Pixel(65535, 0, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, 1.5f, 0.25f), Pixel(0, 65535, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, NAN, 0.25f), Pixel(0, 65535, 0));
  for (float amount : {-4.75f, -1e-7f, 3.5e8f, -3.5e8f, INFINITY, NAN})
    HS_EXPECT_EQ(sample_hue_rotation_lut(view, 0.0f, amount),
                 Pixel(65535, 0, 0));
}
