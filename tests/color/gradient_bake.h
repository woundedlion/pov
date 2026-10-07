/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Gradient::get
// ============================================================================

/**
 * @brief Verifies a black->white gradient yields black at t=0 and white at t~1.
 */
inline void test_gradient_endpoints() {
  // Black at t=0, white at t=1, with OKLCH interpolation between stops.
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  Color4 c0 = grad.get(0.0f);
  HS_EXPECT_EQ(c0.color.r, 0);
  HS_EXPECT_EQ(c0.color.g, 0);
  HS_EXPECT_EQ(c0.color.b, 0);

  Color4 c1 = grad.get(1.0f);
  HS_EXPECT_GT(c1.color.r, 60000);
  HS_EXPECT_GT(c1.color.g, 60000);
  HS_EXPECT_GT(c1.color.b, 60000);
  HS_EXPECT_EQ(c1.alpha, 1.0f);
}

/**
 * @brief Verifies the in-range ramp yields a non-decreasing red channel.
 */
inline void test_gradient_in_range_valid_and_monotone() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  uint16_t prev = 0;
  for (int i = 0; i <= 100; ++i) {
    float t = i / 100.0f;
    Color4 c = grad.get(t);
    if (i > 0)
      HS_EXPECT_GE(c.color.r, prev);
    prev = c.color.r;
  }
}

/**
 * @brief Verifies a gradient with identical stops returns one color for any t.
 */
inline void test_gradient_solid_color() {
  Gradient grad{{0.0f, CPixel(10u, 20u, 30u)}, {1.0f, CPixel(10u, 20u, 30u)}};
  const Pixel want = oklch_to_pixel(srgb_to_oklch(10u, 20u, 30u));
  constexpr float ROUND_TRIP_TOL = 16.0f;
  HS_EXPECT_NEAR(want.r, srgb_to_linear(10u), ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(want.g, srgb_to_linear(20u), ROUND_TRIP_TOL);
  HS_EXPECT_NEAR(want.b, srgb_to_linear(30u), ROUND_TRIP_TOL);
  for (float t : {0.0f, 0.25f, 0.5f, 1.0f}) {
    const Color4 got = grad.get(t);
    HS_EXPECT_EQ(got.color.r, want.r);
    HS_EXPECT_EQ(got.color.g, want.g);
    HS_EXPECT_EQ(got.color.b, want.b);
  }
}

/**
 * @brief Verifies Gradient::get interpolates between LUT entries.
 * @details Two t values inside the same cell (both truncate to index 200) yield
 *          distinct colors, each bracketed by its neighbouring entries.
 */
inline void test_gradient_interpolates_between_entries() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  Color4 lo = grad.get(200.25f / 255.0f);
  Color4 hi = grad.get(200.75f / 255.0f);
  HS_EXPECT_GT(hi.color.r, lo.color.r);

  Color4 e_lo = grad.get(200.0f / 255.0f);
  Color4 e_hi = grad.get(201.0f / 255.0f);
  HS_EXPECT_GE(lo.color.r, e_lo.color.r);
  HS_EXPECT_LE(hi.color.r, e_hi.color.r);
}

/**
 * @brief Verifies Gradient::get clamps t to [0,1] for out-of-range input.
 * @details Out-of-range input saturates to an endpoint and never indexes past
 *          the 256-entry table. NaN folds to the hi bound.
 */
inline void test_gradient_get_clamps_out_of_range() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // t < 0 saturates to the first entry (black).
  Color4 lo_end = grad.get(0.0f);
  Color4 below = grad.get(-0.5f);
  HS_EXPECT_EQ(below.color.r, lo_end.color.r);
  HS_EXPECT_EQ(below.color.g, lo_end.color.g);
  HS_EXPECT_EQ(below.color.b, lo_end.color.b);

  // t > 1 saturates to the last entry (white), same as the t==1 endpoint.
  Color4 hi_end = grad.get(1.0f);
  Color4 above = grad.get(1.5f);
  HS_EXPECT_EQ(above.color.r, hi_end.color.r);
  HS_EXPECT_EQ(above.color.g, hi_end.color.g);
  HS_EXPECT_EQ(above.color.b, hi_end.color.b);

  // NaN folds to the hi bound via hs::clamp -> last entry.
  Color4 nan_res = grad.get(NAN);
  HS_EXPECT_EQ(nan_res.color.r, hi_end.color.r);
}

/**
 * @brief Verifies a first stop at pos>0 flat-fills the LUT prefix with its color.
 * @details The constructor fills entries[0..first_stop] with the first stop's
 *          color, so a gradient whose first stop sits at 0.25 returns that color
 *          for all t in [0, 0.25]; only past the stop do the flanks interpolate.
 */
inline void test_gradient_first_stop_offset_flat_fills_prefix() {
  // First stop (pure red) at 0.25; second (pure blue) at 1.0.
  Gradient grad{{0.25f, CPixel(255u, 0u, 0u)}, {1.0f, CPixel(0u, 0u, 255u)}};

  Color4 at0 = grad.get(0.0f);
  Color4 flat = grad.get(0.1f); // still inside the [0,0.25] flat prefix
  // Pure red.
  HS_EXPECT_GT(at0.color.r, 60000);
  HS_EXPECT_EQ(at0.color.g, 0);
  HS_EXPECT_EQ(at0.color.b, 0);
  HS_EXPECT_EQ(flat.color.r, at0.color.r);
  HS_EXPECT_EQ(flat.color.b, at0.color.b);
  // Past the first stop the flank interpolates toward blue.
  Color4 ramp = grad.get(0.7f);
  HS_EXPECT_LT(ramp.color.r, at0.color.r);
  HS_EXPECT_GT(ramp.color.b, at0.color.b);
}

/**
 * @brief Verifies a >=3-stop gradient places the interior stop and interpolates flanks.
 * @details The interior color must appear near its position and each flanking
 *          segment must blend between its bracketing stops.
 */
inline void test_gradient_three_stops_interior_and_flanks() {
  // red -> green (interior, 0.5) -> blue.
  Gradient grad{{0.0f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(0u, 255u, 0u)},
                {1.0f, CPixel(0u, 0u, 255u)}};
  // Endpoints are the pure stops.
  Color4 a = grad.get(0.0f);
  Color4 b = grad.get(1.0f);
  HS_EXPECT_GT(a.color.r, 60000);
  HS_EXPECT_EQ(a.color.g, 0);
  HS_EXPECT_GT(b.color.b, 60000);
  HS_EXPECT_EQ(b.color.r, 0);
  // Interior green dominates at its stop.
  Color4 mid = grad.get(0.5f);
  HS_EXPECT_GT(mid.color.g, mid.color.r);
  HS_EXPECT_GT(mid.color.g, mid.color.b);
  // A flank between two saturated stops leaves the gamut, and the clip's
  // residual under-saturation pulls the result a few LSB off the cube face, so
  // the absent channel is near zero rather than zero.
  constexpr uint16_t ABSENT = 16;
  // First flank (red->green): both red and green present, blue absent.
  Color4 f1 = grad.get(0.25f);
  HS_EXPECT_GT(f1.color.r, 0);
  HS_EXPECT_GT(f1.color.g, 0);
  HS_EXPECT_LT(f1.color.b, ABSENT);
  // Second flank (green->blue): green and blue present, red absent.
  Color4 f2 = grad.get(0.75f);
  HS_EXPECT_GT(f2.color.g, 0);
  HS_EXPECT_GT(f2.color.b, 0);
  HS_EXPECT_LT(f2.color.r, ABSENT);
}

/**
 * @brief Verifies two stops at the same quantized index produce a hard stop.
 * @details Coincident stop positions leave end==start, so the segment is skipped
 *          and the LUT jumps abruptly: near-pure red below the boundary and
 *          near-pure blue above it.
 */
inline void test_gradient_hard_stop_is_abrupt() {
  Gradient grad{{0.0f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(0u, 0u, 255u)},
                {1.0f, CPixel(0u, 0u, 255u)}};
  // Entry immediately below the boundary is essentially pure red.
  Color4 lo = grad.get(127.0f / 255.0f);
  HS_EXPECT_GT(lo.color.r, 60000);
  HS_EXPECT_LT(lo.color.b, 100);
  // Entry at the boundary is essentially pure blue.
  Color4 hi = grad.get(128.0f / 255.0f);
  HS_EXPECT_GT(hi.color.b, 60000);
  HS_EXPECT_LT(hi.color.r, 100);
}

// ============================================================================
// BakedPalette::get  (requires an Arena)

template <typename T>
concept RebakeablePalette = requires(
    T &palette, const SolidColorPalette &source) { palette.rebake(source); };
static_assert(!RebakeablePalette<BakedPalette>);
static_assert(RebakeablePalette<BakedPaletteStorage>);
static_assert(std::is_copy_constructible_v<BakedPalette>);
static_assert(!std::is_copy_constructible_v<BakedPaletteStorage>);
static_assert(!std::is_copy_assignable_v<BakedPaletteStorage>);
static_assert(!std::is_convertible_v<BakedPaletteStorage &, BakedPalette &>);
static_assert(std::is_nothrow_move_constructible_v<BakedPaletteStorage>);
static_assert(std::is_trivially_destructible_v<BakedPaletteStorage>);
static_assert(sizeof(BakedPaletteStorage) == 2 * sizeof(void *));

inline void test_baked_palette_storage_and_views() {
  alignas(std::max_align_t)
      uint8_t buffer[3 * BakedPalette::required_arena_bytes()];
  Arena arena(buffer, sizeof(buffer));
  SolidColorPalette red(Color4(Pixel(65535, 0, 0), 1.0f));
  SolidColorPalette blue(Color4(Pixel(0, 0, 65535), 0.5f));
  BakedPaletteStorage writer;
  writer.bake(arena, red);
  const BakedPalette view = writer.view();
  BakedPaletteStorage clone;
  clone.clone_from(view, arena);
  const size_t before = arena.get_offset();
  const BakedPalette endpoint = bake_palette_blend(arena, view, clone, 0.0f);
  HS_EXPECT_EQ(arena.get_offset(), before);
  BakedPaletteStorage moved = std::move(writer);
  moved.rebake(blue);
  HS_EXPECT_EQ(view.get_color(0.3f).b, 65535);
  HS_EXPECT_EQ(endpoint.get_color(0.3f).b, 65535);
  HS_EXPECT_EQ(clone.get_color(0.3f).r, 65535);
  writer.bake(arena, red);
  HS_EXPECT_EQ(writer.get_color(0.3f).r, 65535);
  HS_EXPECT_EQ(moved.get_color(0.3f).b, 65535);
}
// ============================================================================

/**
 * @brief Verifies baked endpoint samples match a solid-color source.
 */
inline void test_baked_palette_matches_source_endpoints() {
  // Source: solid color palette so every entry is identical.
  Color4 target(Pixel(1000, 2000, 3000), 1.0f);
  SolidColorPalette src(target);

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPaletteStorage baked;
  baked.bake(arena, src);

  Color4 c0 = baked.get(0.0f);
  Color4 c1 = baked.get(1.0f);
  HS_EXPECT_EQ(c0.color.r, target.color.r);
  HS_EXPECT_EQ(c0.color.g, target.color.g);
  HS_EXPECT_EQ(c0.color.b, target.color.b);
  HS_EXPECT_EQ(c1.color.r, target.color.r);
  HS_EXPECT_NEAR(c0.alpha, 1.0f, 1e-6f);
}

/**
 * @brief Verifies baking a black->white gradient preserves the envelope.
 * @details Preserves the endpoints and keeps interior samples within the
 *          [black,white] envelope at full alpha.
 */
inline void test_baked_palette_in_range() {
  // Ramp source via Gradient (black->white), bake, then sample.
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPaletteStorage baked;
  baked.bake(arena, grad);

  // Endpoints: t=0 -> first entry (black), t=1 -> last entry (white).
  Color4 c0 = baked.get(0.0f);
  Color4 c1 = baked.get(1.0f);
  HS_EXPECT_EQ(c0.color.r, 0);
  HS_EXPECT_GT(c1.color.r, 60000);

  // Interior samples stay within the [c0, c1] envelope and alpha == 1.
  for (int i = 0; i <= 50; ++i) {
    float t = i / 50.0f;
    Color4 c = baked.get(t);
    HS_EXPECT_GE(c.color.r, c0.color.r);
    HS_EXPECT_LE(c.color.r, c1.color.r);
    HS_EXPECT_NEAR(c.alpha, 1.0f, 1e-6f);
  }
}

/**
 * @brief Verifies rebake samples the closed [0, 1] with divisor LUT_SIZE - 1.
 * @details The t = 1 endpoint pins the divisor. BakedPalette rejects Wrap=true
 *          sources at compile time.
 */
inline void test_baked_palette_rebake_samples_closed_interval() {
  struct Ramp {
    Color4 get(float t) const {
      return Color4(Pixel(static_cast<uint16_t>(t * 65535.0f + 0.5f), 0, 0),
                    1.0f);
    }
  } ramp;

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPaletteStorage baked;
  baked.bake(arena, ramp);
  HS_EXPECT_EQ(baked.get(0.0f).color.r, 0);
  HS_EXPECT_EQ(baked.get(1.0f).color.r, 65535);

  using Wrapped = StaticPalette<Gradient>;
  using Unwrapped = StaticPalette<Gradient, Coords<>, Colors<>, /*Wrap=*/false>;
  static_assert(Wrapped::WRAPS_COORDINATE);
  static_assert(!Unwrapped::WRAPS_COORDINATE);
  static_assert(!palette_wraps_coordinate<PaletteFacade<Unwrapped>>());
  static_assert(!palette_wraps_coordinate<Ramp>());
}

/**
 * @brief Produces a palette lookup coordinate.
 * @param i Sample number.
 * @return i / 65535 as a float.
 * @details noinline so every sampler reads one evaluation of the coordinate:
 * the samplers agree for a given index, not for a t each call site recomputes.
 */
__attribute__((noinline)) inline float palette_sample_coord(int i) {
  return static_cast<float>(i) / 65535.0f;
}

inline void test_baked_palette_color_sampler_matches_get() {
  struct Source {
    Color4 get(float t) const {
      return Color4(Pixel(static_cast<uint16_t>(123.0f + 59877.0f * t),
                          static_cast<uint16_t>(4567.0f + 17655.0f * t),
                          static_cast<uint16_t>(65535.0f - 65518.0f * t)),
                    0.1f + 0.8f * t);
    }
  } source;

  alignas(std::max_align_t) uint8_t buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage baked;
  baked.bake(arena, source);

  HS_EXPECT_NEAR(baked.get(0.0f).alpha, 0.1f, 2e-5f);
  HS_EXPECT_NEAR(baked.get(0.5f).alpha, 0.5f, 2e-5f);
  HS_EXPECT_NEAR(baked.get(1.0f).alpha, 0.9f, 2e-5f);

  for (int i = -128; i <= 65663; ++i) {
    const float t = palette_sample_coord(i);
    HS_EXPECT_EQ(baked.get_color(t), baked.get(t).color);
    if (i >= 0 && i <= 65535)
      HS_EXPECT_EQ(baked.get_color_unit(t), baked.get(t).color);
  }
  const float nan = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_EQ(baked.get_color(nan), baked.get(nan).color);
}

/**
 * @brief Verifies clone_from deep-copies the LUT so both palettes sample equal.
 * @details Bake a ramp, copy its LUT into a fresh arena allocation, and assert
 *          both palettes reproduce the same color at every sampled coordinate.
 */
inline void test_baked_palette_clone_from_matches_source() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  alignas(std::max_align_t) static uint8_t
      buf[2 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPaletteStorage src;
  src.bake(arena, grad);
  BakedPaletteStorage dst;
  dst.clone_from(src, arena);

  for (int i = 0; i <= 64; ++i) {
    float t = i / 64.0f;
    Color4 a = src.get(t);
    Color4 b = dst.get(t);
    HS_EXPECT_EQ(a.color.r, b.color.r);
    HS_EXPECT_EQ(a.color.g, b.color.g);
    HS_EXPECT_EQ(a.color.b, b.color.b);
    HS_EXPECT_NEAR(a.alpha, b.alpha, 1e-6f);
  }

  const Color4 CLONED = dst.get(0.5f);
  SolidColorPalette replacement(Color4(Pixel(123, 456, 789), 0.25f));
  src.rebake(replacement);
  HS_EXPECT_PIXEL(src.get(0.5f).color, 123, 456, 789);
  HS_EXPECT_PIXEL(dst.get(0.5f).color, CLONED.color.r, CLONED.color.g,
                  CLONED.color.b);
  HS_EXPECT_EQ(dst.get(0.5f).alpha, CLONED.alpha);
}

/**
 * @brief Verifies dot_key inverts DotKeyed's coordinate mapping.
 * @details A DotKeyed bake maps a LUT coordinate u to the cos value
 *          d = 1 - 2u and samples the source at acos(d)/PI; the fragment path
 *          keys that LUT by the raw dot product through dot_key(d). Sweeps the
 *          closed cos domain, pins the pole orientation, and checks that an
 *          out-of-domain dot product clamps onto a pole.
 */
inline void test_dot_key_inverts_dot_keyed_coordinate() {
  // Captures the angle coordinate a DotKeyed bake hands its source.
  struct CoordProbe {
    mutable float last_t = -1.0f;
    Color4 get(float t) const {
      last_t = t;
      return Color4(Pixel(0, 0, 0), 1.0f);
    }
  };
  const CoordProbe probe;
  const auto keyed = dot_keyed(probe);

  // d = +1 is the axis (angle 0), d = -1 the antipode (angle PI).
  HS_EXPECT_NEAR(dot_key(1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(dot_key(-1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(dot_key(0.0f), 0.5f, 1e-6f);

  constexpr int STEPS = 64;
  float previous_u = 2.0f;
  for (int i = 0; i <= STEPS; ++i) {
    const float d = -1.0f + 2.0f * (static_cast<float>(i) / STEPS);
    const float u = dot_key(d);
    HS_EXPECT_TRUE(u >= 0.0f && u <= 1.0f);
    if (i > 0)
      HS_EXPECT_LT(u, previous_u);
    previous_u = u;
    // The bake's u -> d leg recovers the dot product dot_key was handed.
    HS_EXPECT_NEAR(1.0f - 2.0f * u, d, 1e-6f);
    keyed.get(u);
    HS_EXPECT_NEAR(probe.last_t, math::fast_acos(d) / math::PI_F, 1e-6f);
  }

  HS_EXPECT_NEAR(dot_key(4.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(dot_key(-4.0f), 1.0f, 1e-6f);
}

/**
 * @brief Verifies a DotKeyed bake sampled at dot_key reproduces its source.
 * @details Bakes through dot_keyed(), then looks the LUT up by the raw dot
 *          product. A black-to-white ramp over t = angle/PI must land black on
 *          the axis, white on the antipode, and rise monotonically between
 *          them; interior samples match the source within LUT quantization.
 */
inline void test_dot_keyed_bake_round_trips_through_dot_key() {
  Gradient ramp{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage baked;
  baked.bake(arena, dot_keyed(ramp));

  HS_EXPECT_EQ(baked.get(dot_key(1.0f)).color.r, ramp.get(0.0f).color.r);
  HS_EXPECT_EQ(baked.get(dot_key(-1.0f)).color.r, ramp.get(1.0f).color.r);

  constexpr int STEPS = 64;
  int previous = -1;
  for (int i = STEPS; i >= 0; --i) {
    const float d = -1.0f + 2.0f * (static_cast<float>(i) / STEPS);
    const int got = baked.get(dot_key(d)).color.r;
    if (i < STEPS)
      HS_EXPECT_TRUE(got >= previous);
    previous = got;
  }

  // Away from the poles the LUT reproduces the source closely.
  for (int i = -3; i <= 3; ++i) {
    const float d = static_cast<float>(i) / 4.0f;
    const uint16_t got = baked.get(dot_key(d)).color.r;
    const uint16_t want = ramp.get(math::fast_acos(d) / math::PI_F).color.r;
    HS_EXPECT_NEAR(got, want, 32);
  }
}

/**
 * @brief Verifies a NaN blend weight selects the destination LUT.
 */
inline void test_bake_palette_blend_nan_weight_stays_finite() {
  SolidColorPalette black(Color4(Pixel(0, 0, 0), 0.25f));
  SolidColorPalette white(Color4(Pixel(65535, 65535, 65535), 1.0f));

  alignas(std::max_align_t) static uint8_t
      buf[3 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage from, to;
  from.bake(arena, black);
  to.bake(arena, white);

  BakedPalette dst;
  dst = bake_palette_blend(arena, from, to,
                           std::numeric_limits<float>::quiet_NaN());

  for (int i = 0; i <= 64; ++i) {
    Color4 c = dst.get(i / 64.0f);
    HS_EXPECT_TRUE(std::isfinite(c.alpha));
    HS_EXPECT_NEAR(c.alpha, 1.0f, 1e-6f);
    HS_EXPECT_EQ(c.color.r, to.get(i / 64.0f).color.r);
    HS_EXPECT_EQ(c.color.g, to.get(i / 64.0f).color.g);
    HS_EXPECT_EQ(c.color.b, to.get(i / 64.0f).color.b);
  }
}

/**
 * @brief Verifies step_wipe_rebake skips the arming frame, then decrements.
 * @details A ColorWipe is armed mid-step and first steps next frame, so the
 *          arming frame must be consumed without touching the frame counter;
 *          each later frame decrements it, and the count never underflows once
 *          exhausted.
 */
inline void test_step_wipe_rebake_skips_arming_then_decrements() {
  SolidColorPalette initial(Color4(Pixel(1000, 2000, 3000), 1.0f));
  SolidColorPalette src(Color4(Pixel(4000, 5000, 6000), 1.0f));
  SolidColorPalette after(Color4(Pixel(7000, 8000, 9000), 1.0f));

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage baked;
  baked.bake(arena, initial);

  bool wipe_pending = true;
  int frames = 2;

  step_wipe_rebake(wipe_pending, frames, baked, src);
  HS_EXPECT_FALSE(wipe_pending);
  HS_EXPECT_EQ(frames, 2);
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 1000, 2000, 3000);

  step_wipe_rebake(wipe_pending, frames, baked, src);
  HS_EXPECT_EQ(frames, 1);
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 4000, 5000, 6000);
  step_wipe_rebake(wipe_pending, frames, baked, src);
  HS_EXPECT_EQ(frames, 0);

  step_wipe_rebake(wipe_pending, frames, baked, after);
  HS_EXPECT_EQ(frames, 0);
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 4000, 5000, 6000);
}

/**
 * @brief Verifies PaletteWipe's arm/step cadence, zero-frame arm included.
 * @details arm() captures the endpoints and opens the rebake window; the
 *          arming frame is consumed without spending a rebake frame, the wipe
 *          reports in flight for exactly the armed count, and the counter
 *          floors at zero. A zero-frame arm is never in flight and its arming
 *          frame still clears, so the state is reusable without a guard.
 */
inline void test_palette_wipe_arm_step_cadence() {
  const GenerativePalette palette;
  SolidColorPalette initial(Color4(Pixel(1000, 2000, 3000), 1.0f));
  SolidColorPalette source(Color4(Pixel(4000, 5000, 6000), 1.0f));
  SolidColorPalette after(Color4(Pixel(7000, 8000, 9000), 1.0f));

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage baked;
  baked.bake(arena, initial);

  PaletteWipe wipe;
  HS_EXPECT_FALSE(wipe.in_flight());
  HS_EXPECT_FALSE(wipe.pending);

  constexpr int FRAMES = 3;
  wipe.arm(palette, palette.snapshot(), FRAMES);
  HS_EXPECT_TRUE(wipe.pending);
  HS_EXPECT_TRUE(wipe.in_flight());
  HS_EXPECT_EQ(wipe.frames_remaining, FRAMES);

  wipe.step(baked, source);
  HS_EXPECT_FALSE(wipe.pending);
  HS_EXPECT_EQ(wipe.frames_remaining, FRAMES);
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 1000, 2000, 3000);

  for (int left = FRAMES; left > 0; --left) {
    HS_EXPECT_TRUE(wipe.in_flight());
    wipe.step(baked, source);
    HS_EXPECT_EQ(wipe.frames_remaining, left - 1);
    HS_EXPECT_PIXEL(baked.get_color(0.5f), 4000, 5000, 6000);
  }
  HS_EXPECT_FALSE(wipe.in_flight());

  // Steps past the window neither underflow the counter nor re-arm.
  wipe.step(baked, after);
  HS_EXPECT_EQ(wipe.frames_remaining, 0);
  HS_EXPECT_FALSE(wipe.pending);
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 4000, 5000, 6000);

  wipe.arm(palette, palette.snapshot(), 0);
  HS_EXPECT_FALSE(wipe.in_flight());
  HS_EXPECT_TRUE(wipe.pending);
  wipe.step(baked, after);
  HS_EXPECT_FALSE(wipe.pending);
  HS_EXPECT_EQ(wipe.frames_remaining, 0);
  HS_EXPECT_FALSE(wipe.in_flight());
  HS_EXPECT_PIXEL(baked.get_color(0.5f), 4000, 5000, 6000);
}
