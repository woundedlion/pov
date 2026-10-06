/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Color death fixtures and guard cases.

inline void case_gamut_lut_scratch_a() {
  init_gamut_lut(scratch_arena_a, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
}

inline void case_gamut_lut_scratch_b() {
  init_gamut_lut(scratch_arena_b, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
}

inline void case_noise_hue_bake_invalid_scale() {
  static std::array<int8_t, HueNoiseLutView::SIZE> output{};
  FastNoiseLite noise;
  HueNoiseBakeCache cache;
  (void)cache.refresh(output, noise, opaque(0.0f), 0.0f);
}

/**
 * @brief Death case: a HueWobbleShade depth past the fast-trig argument range.
 */
inline void case_hue_wobble_depth_out_of_range() {
  static float phase = 0.0f;
  HueWobbleShade shade(&phase, 1.0f, opaque(1.0e6f));
  (void)shade;
}

/**
 * @brief Death case: a negative IridescentShade weight zeroes the overlay.
 */
inline void case_iridescent_weight_negative() {
  static float phase = 0.0f;
  IridescentShade shade(&phase, 3.0f, opaque(-0.5f));
  (void)shade;
}

/**
 * @brief Death case: AlphaFalloffShade requires a non-null callback.
 */
inline void case_alpha_falloff_null() {
  auto fn = opaque<AlphaFalloffShade::FalloffFunction>(nullptr);
  AlphaFalloffShade shade(fn);
  (void)shade;
}

inline void case_palette_cycler_mutated_policy() {
  static uint8_t storage[8192];
  Arena arena(storage, sizeof(storage));
  GenerativePalette first, second;
  const PaletteCycler::Entry entries[] = {first, second};
  PaletteCycler cycler;
  cycler.init(arena, entries, 2, 0, 4);
  PaletteRecipe changed;
  changed.chroma.headroom = 0.5f;
  second = GenerativePalette(changed);
  cycler.step();
  cycler.step();
}

inline void case_palette_cycler_restore_without_generated_init() {
  PaletteCycler cycler;
  GenerativePalette palette;
  cycler.restore_generated({}, palette, palette);
}

/**
 * @brief Death case: GeneratedPaletteBank::palette rejects a mode outside the
 *        harmony enum.
 */
inline void case_generated_palette_bank_unknown_mode() {
  enum class Mode : uint8_t { TRIADIC, COMPLEMENTARY, ANALOGOUS };
  GeneratedPaletteBank bank;
  (void)bank.palette(opaque(static_cast<Mode>(3)));
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null palette source.
 */
inline void case_noise_hue_palette_direct_null_source() {
  static int8_t noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), noise_lut);
}

inline void case_noise_hue_palette_direct_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const int8_t *>(nullptr));
}

inline void case_noise_shimmer_palette_null_source() {
  static int8_t noise_lut[1];
  NoiseShimmerPalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), noise_lut);
}

inline void case_noise_shimmer_palette_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  NoiseShimmerPalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const int8_t *>(nullptr));
}

inline void case_noise_hue_palette_null_source() {
  static Pixel hue_rotation_lut[1];
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), hue_rotation_lut,
               hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-rotation LUT.
 */
inline void case_noise_hue_palette_null_rotation_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const Pixel *>(nullptr), hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-noise LUT.
 */
inline void case_noise_hue_palette_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static Pixel hue_rotation_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, hue_rotation_lut, opaque<const int8_t *>(nullptr));
}

/**
 * @brief Death case: cloning a BakedPalette from itself must trap.
 * @details Color surface — clone_from allocates fresh storage into this handle
 *          before reading @c src, so a self-clone memcpys uninitialized arena
 *          onto itself and leaves the LUT filled with garbage.
 */
inline void case_baked_palette_clone_from_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPaletteStorage lut;
  lut.bake(a, src);
  const BakedPalette &self = opaque(&lut)->view();
  lut.clone_from(self, a); // -> HS_CHECK
}

/**
 * @brief Death case: blending a BakedPalette with itself as an endpoint must
 *        trap.
 * @details Color surface — bake_blend reallocates this handle before walking
 *          the endpoints, so an endpoint that is the output reads the fresh
 *          uninitialized arena instead of the baked LUT.
 */
inline void case_baked_palette_bake_blend_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPaletteStorage from, dst;
  from.bake(a, src);
  dst.bake(a, src);
  const BakedPalette &self = opaque(&dst)->view();
  dst.bake_blend(a, from, self, opaque(0.5f)); // -> HS_CHECK
}

/**
 * @brief Death case: a Gradient built from an empty stop list must trap.
 * @details The constructor requires at least one stop before reading the
 *          first stop's position and color.
 */
inline void case_gradient_no_stops() {
  Gradient grad({}); // empty stop list -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: a Gradient stop position outside [0,1] must trap.
 * @details Color surface — a stop position becomes a rounded LUT index via
 *          static_cast<int>(pos * 255 + 0.5f); sufficiently out-of-range
 *          positions can write beyond the table. The constructor traps the
 *          authoring error always-on at the cold literal-construction seam
 *          rather than corrupting memory.
 */
inline void case_gradient_stop_out_of_range() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)},
                {1.5f, CPixel(255u, 255u, 255u)}}; // pos > 1 -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: descending (unsorted) Gradient stops must trap.
 * @details Color surface — segments are only filled when end > start, so a
 *          transposed/unsorted pair would silently degenerate to wrong output.
 *          The constructor requires ascending positions and traps otherwise.
 */
inline void case_gradient_stops_unsorted() {
  Gradient grad{{0.6f, CPixel(0u, 0u, 0u)},
                {0.3f, CPixel(255u, 255u, 255u)}}; // descending -> HS_CHECK
  (void)grad;
}
