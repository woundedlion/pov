/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Color death cases.

/** @brief Traps on gamut lut scratch a. */
inline void case_gamut_lut_scratch_a() {
  init_gamut_lut(scratch_arena_a, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
}

/** @brief Traps on gamut lut scratch b. */
inline void case_gamut_lut_scratch_b() {
  init_gamut_lut(scratch_arena_b, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
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
 * @brief Death case: a negative IridescentShade weight must trap.
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

/** @brief Checks a valid interval before rejecting an invalid generated cadence. */
inline void reject_palette_display_interval(int interval, int fade_frames) {
  alignas(std::max_align_t) static uint8_t
      storage[2 * PaletteCycler::generated_arena_bytes()];
  Arena arena(storage, sizeof(storage));
  PaletteCycler cycler;
  const auto provider = [](void *, uint32_t, GenerativePalette &out) {
    out = GenerativePalette();
  };
  cycler.init_generated(arena, provider, nullptr, 0, fade_frames, nullptr,
                        nullptr, 2);
  cycler.init_generated(arena, provider, nullptr, 0, fade_frames, nullptr,
                        nullptr, opaque(interval));
}

/** @brief Generated display intervals must be nonzero. */
inline void case_palette_cadence_zero() {
  reject_palette_display_interval(0, 4);
}

/** @brief Generated display intervals cannot exceed their fade. */
inline void case_palette_cadence_past_fade() {
  reject_palette_display_interval(5, 4);
}

/** @brief Generated display intervals must fit their stored width. */
inline void case_palette_cadence_past_storage() {
  reject_palette_display_interval(256, 600);
}

/** @brief Traps on palette cycler mutated policy. */
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

/** @brief Traps on palette cycler restore without generated init. */
inline void case_palette_cycler_restore_without_generated_init() {
  PaletteCycler cycler;
  GenerativePalette palette;
  cycler.restore_generated({}, palette, palette);
}

/** @brief Restoring generated endpoints rejects incompatible policies. */
inline void case_palette_cycler_restore_incompatible() {
  static uint8_t storage[PaletteCycler::generated_arena_bytes()];
  Arena arena(storage, sizeof(storage));
  PaletteCycler cycler;
  cycler.init_generated(
      arena,
      [](void *, uint32_t, GenerativePalette &out) {
        out = GenerativePalette();
      },
      nullptr, 0, 4);
  const GenerativePalette from;
  cycler.restore_generated({.next_sequence = 2, .fade_active = true}, from,
                           from);
  PaletteRecipe recipe;
  recipe.domain = PaletteDomain::MIRROR;
  const GenerativePalette to(recipe);
  cycler.restore_generated({.next_sequence = 2, .fade_active = true}, from, to);
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

/** @brief Sphere field admission rejects an infinite direction. */
inline void case_modulation_infinite_direction() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(INFINITY), 0, 0));
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
 *          positions can write beyond the table.
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
 */
inline void case_gradient_stops_unsorted() {
  Gradient grad{{0.6f, CPixel(0u, 0u, 0u)},
                {0.3f, CPixel(255u, 255u, 255u)}}; // descending -> HS_CHECK
  (void)grad;
}

/** @brief BreatheModifier rejects a null driver. */
inline void case_composition_breathemodifier_null() {
  float driver = 0.0f;
  (void)BreatheModifier(&driver);
  (void)BreatheModifier(opaque<const float *>(nullptr));
}

/** @brief RippleModifier rejects a null driver. */
inline void case_composition_ripplemodifier_null() {
  float driver = 0.0f;
  (void)RippleModifier(&driver);
  (void)RippleModifier(opaque<const float *>(nullptr));
}

/** @brief NoiseWarpModifier rejects a null driver. */
inline void case_composition_noisewarpmodifier_null() {
  float driver = 0.0f;
  (void)NoiseWarpModifier(&driver);
  (void)NoiseWarpModifier(opaque<const float *>(nullptr));
}

/** @brief DriftModifier rejects a null driver. */
inline void case_composition_driftmodifier_null() {
  float driver = 0.0f;
  (void)DriftModifier(&driver);
  (void)DriftModifier(opaque<const float *>(nullptr));
}

/** @brief ScaleModifier rejects a null driver. */
inline void case_composition_scalemodifier_null() {
  float driver = 0.0f;
  (void)ScaleModifier(&driver);
  (void)ScaleModifier(opaque<const float *>(nullptr));
}

/** @brief HueSpinShade rejects a null driver. */
inline void case_composition_huespinshade_null() {
  float driver = 0.0f;
  (void)HueSpinShade(&driver);
  (void)HueSpinShade(opaque<const float *>(nullptr));
}

/** @brief HueWobbleShade rejects a null driver. */
inline void case_composition_huewobbleshade_null() {
  float driver = 0.0f;
  (void)HueWobbleShade(&driver);
  (void)HueWobbleShade(opaque<const float *>(nullptr));
}

/** @brief SparkleShade rejects a null driver. */
inline void case_composition_sparkleshade_null() {
  float driver = 0.0f;
  (void)SparkleShade(&driver);
  (void)SparkleShade(opaque<const float *>(nullptr));
}

/** @brief ChromaPulseShade rejects a null driver. */
inline void case_composition_chromapulseshade_null() {
  float driver = 0.0f;
  (void)ChromaPulseShade(&driver);
  (void)ChromaPulseShade(opaque<const float *>(nullptr));
}

/** @brief LightnessGrainShade rejects a null driver. */
inline void case_composition_lightnessgrainshade_null() {
  float driver = 0.0f;
  (void)LightnessGrainShade(&driver);
  (void)LightnessGrainShade(opaque<const float *>(nullptr));
}

/** @brief IridescentShade rejects a null driver. */
inline void case_composition_iridescentshade_null() {
  float driver = 0.0f;
  (void)IridescentShade(&driver);
  (void)IridescentShade(opaque<const float *>(nullptr));
}

/** @brief InsetModifier rejects an invalid range. */
inline void case_composition_inset_bounds() {
  (void)InsetModifier(0.0f, 1.0f);
  (void)InsetModifier(1.0f, 1.0f);
}

/** @brief SparkleShade rejects an invalid range. */
inline void case_composition_sparkle_threshold() {
  float driver = 0.0f;
  (void)SparkleShade(&driver, 1.0f, 0.5f);
  (void)SparkleShade(&driver, 1.0f, 1.0f);
}

/** @brief ChromaPulseShade rejects an invalid range. */
inline void case_composition_chroma_depth() {
  float driver = 0.0f;
  (void)ChromaPulseShade(&driver, 1.0f);
  (void)ChromaPulseShade(&driver, -1.0f);
}

/** @brief LightnessGrainShade rejects an invalid range. */
inline void case_composition_grain_amplitude() {
  float driver = 0.0f;
  (void)LightnessGrainShade(&driver, 1.0f, 1.0f);
  (void)LightnessGrainShade(&driver, 1.0f, 2.0f);
}

/** @brief EdgeFadeShade rejects an invalid range. */
inline void case_composition_edge_fade() {
  (void)EdgeFadeShade(0.5f);
  (void)EdgeFadeShade(0.0f);
}

/** @brief EdgeAlphaShade rejects an invalid range. */
inline void case_composition_edge_alpha() {
  (void)EdgeAlphaShade(0.5f);
  (void)EdgeAlphaShade(0.0f);
}

/** @brief StaticPalette rejects a missing source. */
inline void case_composition_null_source() {
  GenerativePalette source;
  CycleModifier coord;
  EdgeFadeShade color;
  StaticPalette<GenerativePalette, Coords<CycleModifier>, Colors<EdgeFadeShade>>
      palette;
  palette.bind(&source, &coord, &color);
  palette.bind(nullptr, &coord, &color);
}

/** @brief StaticPalette rejects a missing coord modifier. */
inline void case_composition_null_coord() {
  GenerativePalette source;
  CycleModifier coord;
  EdgeFadeShade color;
  StaticPalette<GenerativePalette, Coords<CycleModifier>, Colors<EdgeFadeShade>>
      palette;
  palette.bind(&source, &coord, &color);
  palette.bind(&source, nullptr, &color);
}

/** @brief StaticPalette rejects a missing color modifier. */
inline void case_composition_null_color() {
  GenerativePalette source;
  CycleModifier coord;
  EdgeFadeShade color;
  StaticPalette<GenerativePalette, Coords<CycleModifier>, Colors<EdgeFadeShade>>
      palette;
  palette.bind(&source, &coord, &color);
  palette.bind(&source, &coord, nullptr);
}

/** @brief PaletteFacade rejects a null composition. */
inline void case_composition_facade_construct() {
  using Composition =
      StaticPalette<GenerativePalette, Coords<>, Colors<>, false>;
  Composition palette;
  PaletteFacade<Composition> valid(&palette);
  (void)PaletteFacade<Composition>(opaque<const Composition *>(nullptr));
}

/** @brief PaletteFacade rejects a null composition. */
inline void case_composition_facade_bind() {
  using Composition =
      StaticPalette<GenerativePalette, Coords<>, Colors<>, false>;
  Composition palette;
  PaletteFacade<Composition> valid(&palette);
  valid.bind(nullptr);
}
