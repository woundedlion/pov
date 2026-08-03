/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_COLOR_INTERNAL
#error internal fragment of color.h; include "color/color.h" instead
#endif

// Palette-composition layer for color.h: the coordinate/color modifiers, the
// StaticPalette composition template, and baked-palette storage. Included by
// color.h after the palette core.

///////////////////////////////////////////////////////////////////////////////
// Palette Modifiers
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Linearly cycles the palette coordinate.
 *
 * Null offset driver is a deliberate "no cycling, static" mode (modify() passes
 * t through), not an error.
 */
struct CycleModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *offset;

  /**
   * @brief Constructs with an optional offset driver.
   * @param driver_offset Pointer to the per-frame offset, or null for static.
   */
  CycleModifier(const float *driver_offset = nullptr) : offset(driver_offset) {}

  /**
   * @brief Shifts the coordinate by the driver offset (pass-through if null).
   * @param t Input coordinate.
   * @return t plus the offset, or t unchanged when no driver is bound.
   * @note The result intentionally leaves [0,1] (t + a monotonically growing
   *       offset), relying on the consuming palette's `Wrap=true` (the
   *       `StaticPalette` default) to fold it back into range and produce the
   *       cycling. `requires_wrap` makes a `Wrap=false` composition a compile
   *       error.
   */
  float modify(float t) const { return offset ? t + *offset : t; }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *phase;
  float amplitude;
  /**
   * @brief Per-instance memo of fast_sinf(*phase).
   * @details *phase is frame-constant, so the sine is recomputed once per frame,
   * not per pixel. mutable so const modify() can update the memo.
   */
  mutable float cached_phase = 0.0f;
  mutable float cached_sin = 0.0f; /**< Memoized sine of cached_phase. */
  mutable bool primed = false;     /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory phase driver and amplitude.
   * @param driver_phase Pointer to the per-frame phase; must not be null.
   * @param amp Oscillation amplitude; defaults to 0.1.
   * @details Mandatory phase driver: a null one is trapped at construction so
   * per-pixel modify() can dereference unconditionally.
   */
  BreatheModifier(const float *driver_phase, float amp = 0.1f)
      : phase(driver_phase), amplitude(amp) {
    HS_CHECK(phase, "BreatheModifier: phase driver must not be null");
  }

  /**
   * @brief Oscillates the coordinate by amplitude * sin(phase).
   * @param t Input coordinate.
   * @return t plus the memoized oscillation term.
   */
  float modify(float t) const {
    if (!primed || *phase != cached_phase) {
      cached_phase = *phase;
      cached_sin = fast_sinf(*phase);
      primed = true;
    }
    return t + cached_sin * amplitude;
  }
};

/**
 * @brief Distorts the palette spatially with a sine wave, creating a liquid
 * ripple effect. Compresses and expands colors like waves on a spatial coord.
 */
struct RippleModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *phase;
  float frequency;
  float amplitude;

  /**
   * @brief Constructs with a mandatory phase driver, frequency, and amplitude.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param freq Spatial frequency of the ripple; defaults to 3.0.
   * @param amp Distortion amplitude; defaults to 0.1.
   * @details Mandatory phase driver (no default) — trap a null one at
   * construction rather than silently passing t through on every pixel.
   */
  RippleModifier(const float *phase, float freq = 3.0f, float amp = 0.1f)
      : phase(phase), frequency(freq), amplitude(amp) {
    HS_CHECK(phase, "RippleModifier: phase driver must not be null");
  }

  /**
   * @brief Distorts the coordinate with a sine wave of the given frequency.
   * @param t Input coordinate.
   * @return t plus the local sine distortion.
   */
  float modify(float t) const {
    return t + fast_sinf(t * frequency * PI_F * 2.0f + *phase) * amplitude;
  }
};

/**
 * @brief Warps the palette coordinate with smooth value noise — the organic,
 * aperiodic counterpart to RippleModifier's sine: colors wander and smear
 * instead of oscillating.
 */
struct NoiseWarpModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *time;
  float frequency;
  float amplitude;
  uint32_t seed;

  /**
   * @brief Constructs with a mandatory time driver, frequency, and amplitude.
   * @param time Pointer to the per-frame noise time axis; must not be null.
   * @param freq Spatial frequency of the noise over t; defaults to 3.0.
   * @param amp Peak displacement of the coordinate; defaults to 0.1.
   * @param seed Noise stream selector; defaults to 0.
   */
  NoiseWarpModifier(const float *time, float freq = 3.0f, float amp = 0.1f,
                    uint32_t seed = 0)
      : time(time), frequency(freq), amplitude(amp), seed(seed) {
    HS_CHECK(time, "NoiseWarpModifier: time driver must not be null");
  }

  /**
   * @brief Displaces the coordinate by centered 2D noise at (t*frequency, time).
   * @param t Input coordinate.
   * @return t plus a displacement in [-amplitude, amplitude].
   */
  float modify(float t) const {
    return t + (value_noise_2d(t * frequency, *time, seed) - 0.5f) * 2.0f *
                   amplitude;
  }
};

/**
 * @brief Meanders the whole palette along a smooth noise walk — unlike
 * CycleModifier's linear scroll, the offset wanders, hesitates, and reverses.
 */
struct DriftModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *time;
  float speed;
  float amplitude;
  uint32_t seed;
  /**
   * @brief Per-instance memo of the frame's centered walk sample.
   * @details *time is frame-constant, so the noise walk is sampled once per
   * frame, not per pixel. mutable so const modify() can update the memo.
   * Keyed on *time alone, so speed and seed must not change between frames.
   */
  mutable float cached_time = 0.0f;
  mutable float cached_walk = 0.0f; /**< Memoized walk in [-1, 1]. */
  mutable bool primed = false;      /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory time driver, walk speed, and amplitude.
   * @param time Pointer to the per-frame time; must not be null.
   * @param speed Walk rate in noise cells per time unit; defaults to 0.25.
   * @param amp Peak offset; defaults to 0.25.
   * @param seed Noise stream selector; defaults to 0.
   */
  DriftModifier(const float *time, float speed = 0.25f, float amp = 0.25f,
                uint32_t seed = 0)
      : time(time), speed(speed), amplitude(amp), seed(seed) {
    HS_CHECK(time, "DriftModifier: time driver must not be null");
  }

  /**
   * @brief Shifts the coordinate by the frame's noise-walk offset.
   * @param t Input coordinate.
   * @return t plus an offset in [-amplitude, amplitude].
   */
  float modify(float t) const {
    if (!primed || *time != cached_time) {
      cached_time = *time;
      cached_walk = (value_noise_1d(cached_time * speed, seed) - 0.5f) * 2.0f;
      primed = true;
    }
    return t + cached_walk * amplitude;
  }
};

/**
 * @brief Folds the palette back and forth like a kaleidoscope.
 * A folds value of 2.0 maps [0...1] to [1 -> 0 -> 1] (one full bounce);
 * each unit of folds adds another half-bounce.
 *
 * Null phase driver is the deliberate "no phase offset" mode (shift = 0).
 */
struct FoldModifier {
  /** @brief Output stays in [0,1] and hits 1; palette needs Wrap=false. */
  static constexpr bool bounded_output = true;
  /** @brief The triangle wave folds any input, in range or not, into [0,1]. */
  static constexpr bool rebounds_input = true;

  const float *phase;
  float folds;

  /**
   * @brief Constructs with a fold count and optional phase driver.
   * @param folds Number of bounces; defaults to 2.0 (one full bounce).
   * @param phase Pointer to an optional phase offset, or null for none.
   */
  FoldModifier(float folds = 2.0f, const float *phase = nullptr)
      : phase(phase), folds(folds) {}

  /**
   * @brief Folds the coordinate back and forth via a triangle wave.
   * @param t Input coordinate.
   * @return The folded coordinate in [0, 1].
   */
  float modify(float t) const {
    float shift = phase ? *phase : 0.0f;
    float scaled = (t * folds) + shift;

    // Triangle wave. fmodf keeps the dividend's sign, so reduce into [0, 2)
    // first — negative scaled would otherwise fold above 1.
    float m = fmodf(scaled, 2.0f);
    if (m < 0.0f)
      m += 2.0f;
    return fabsf(m - 1.0f);
  }
};

/**
 * @brief Pinches or expands the center of the palette.
 * positive tension pulls colors toward the center, negative pushes them to the
 * edges.
 *
 * Null tension driver is the deliberate "no pinch" pass-through mode.
 */
struct PinchModifier {
  /** @brief In-range input stays in [0,1] and hits 1; palette needs Wrap=false.
   */
  static constexpr bool bounded_output = true;

  const float
      *tension; /**< Pinch tension driver; expects roughly -0.9 to 0.9. */

  /**
   * @brief Constructs with an optional tension driver.
   * @param t Pointer to the tension value, or null for pass-through.
   */
  PinchModifier(const float *t = nullptr) : tension(t) {}

  /**
   * @brief Pinches or expands the coordinate around the domain center.
   * @param t Input coordinate.
   * @return The reshaped coordinate, or t unchanged when no driver is bound.
   */
  float modify(float t) const {
    if (!tension)
      return t;

    // Center the wrapped coordinate into [-1, 1].
    float wrapped_t = wrap_t(t);
    float centered = wrapped_t * 2.0f - 1.0f;
    float sign = centered < 0.0f ? -1.0f : 1.0f;

    float amount = hs::clamp(*tension, -0.99f, 0.99f);
    float power = (amount < 0.0f) ? (1.0f / (1.0f + std::abs(amount)))
                                  : (1.0f + amount * 3.0f);

    centered = sign * powf(std::abs(centered), power);

    // Re-anchor to t's own integer cell: floorf(t) pairs with wrap_t(t), correct
    // even for negative t.
    return floorf(t) + ((centered + 1.0f) * 0.5f);
  }
};

/**
 * @brief Snaps smooth gradients into harsh, distinct bands (Posterization).
 */
struct QuantizeModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *dynamic_steps;
  float base_steps;

  /**
   * @brief Constructs with a base step count and optional dynamic driver.
   * @param steps Base number of quantization steps.
   * @param d_steps Pointer to an animated step count, or null to use base.
   */
  QuantizeModifier(float steps, const float *d_steps = nullptr)
      : dynamic_steps(d_steps), base_steps(steps) {}

  /**
   * @brief Snaps the coordinate to the nearest multiple of 1/steps (steps+1
   * distinct levels over [0,1]).
   * @param t Input coordinate.
   * @return The quantized coordinate.
   */
  float modify(float t) const {
    float s = dynamic_steps ? *dynamic_steps : base_steps;
    if (s < 1.0f)
      s = 1.0f;

    // Round to nearest step in the infinite domain.
    return roundf(t * s) / s;
  }
};

/**
 * @brief Multiplies the palette coordinate, increasing the frequency
 * so the palette repeats multiple times across the domain.
 */
struct ScaleModifier {
  /** @brief Output leaves [0,1]; the consuming palette must have Wrap=true. */
  static constexpr bool requires_wrap = true;

  const float *dynamic_scale;
  float base_scale;

  /**
   * @brief Constructs with a base scale and optional dynamic driver.
   * @param s Base scale factor; defaults to 1.0.
   * @param d_scale Pointer to an animated scale, or null to use base.
   */
  ScaleModifier(float s = 1.0f, const float *d_scale = nullptr)
      : dynamic_scale(d_scale), base_scale(s) {}

  /**
   * @brief Constructs driven purely by an animated scale.
   * @param d_scale Pointer to an animated scale; must not be null.
   */
  ScaleModifier(const float *d_scale)
      : dynamic_scale(d_scale), base_scale(1.0f) {
    HS_CHECK(d_scale != nullptr,
             "ScaleModifier: dynamic scale must not be null");
  }

  /**
   * @brief Multiplies the coordinate by the active scale.
   * @param t Input coordinate.
   * @return The scaled coordinate.
   * @note For scale > 1 the result intentionally leaves [0,1], relying on the
   *       consuming palette's `Wrap=true` (the `StaticPalette` default) to fold
   *       it back into range — that fold IS the multiple-repeats effect.
   *       `requires_wrap` makes a `Wrap=false` composition a compile error.
   */
  float modify(float t) const {
    return t * (dynamic_scale ? *dynamic_scale : base_scale);
  }
};

/**
 * @brief Reverses the palette coordinate (t -> 1 - t).
 */
struct ReverseModifier {
  /** @brief In-range input stays in [0,1] and hits 1; palette needs Wrap=false.
   */
  static constexpr bool bounded_output = true;

  /**
   * @brief Reverses the coordinate.
   * @param t Input coordinate.
   * @return 1 - t.
   */
  float modify(float t) const { return 1.0f - t; }
};

/**
 * @brief Mirrors the coordinate so [0,1] maps to [0,1,0].
 * @details One symmetric bounce, for a seamless loop.
 */
struct MirrorModifier {
  /** @brief In-range input stays in [0,1] and hits 1; palette needs Wrap=false.
   */
  static constexpr bool bounded_output = true;

  /**
   * @brief Mirrors the coordinate into a symmetric bounce.
   * @param t Input coordinate.
   * @return The mirrored coordinate in [0, 1].
   */
  float modify(float t) const { return 1.0f - fabsf(2.0f * t - 1.0f); }
};

/**
 * @brief Compresses the source domain into an inset window [lo, hi] -> [0, 1].
 * @details Clamps outside so t below lo samples the first stop and t above hi
 * the last. Pairs with EdgeFadeShade / EdgeAlphaShade to build vignettes.
 */
struct InsetModifier {
  /** @brief Output stays in [0,1] and hits 1; palette needs Wrap=false. */
  static constexpr bool bounded_output = true;
  /** @brief The clamp confines any input, in range or not, to [0,1]. */
  static constexpr bool rebounds_input = true;

  float lo, hi;
  /**
   * @brief Constructs the inset window bounds.
   * @param lo Lower domain bound mapped to 0; defaults to 0.2.
   * @param hi Upper domain bound mapped to 1; defaults to 0.8.
   */
  InsetModifier(float lo = 0.2f, float hi = 0.8f) : lo(lo), hi(hi) {
    HS_CHECK(hi > lo,
             "InsetModifier: hi must be > lo (modify divides by hi - lo)");
  }
  /**
   * @brief Remaps the coordinate from [lo, hi] into [0, 1], clamping outside.
   * @param t Input coordinate.
   * @return The remapped coordinate in [0, 1].
   */
  float modify(float t) const {
    return hs::clamp((t - lo) / (hi - lo), 0.0f, 1.0f);
  }
};

///////////////////////////////////////////////////////////////////////////////
// Color Modifiers — reshape the sample after the source lookup.
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Rotates every sample's hue in OKLab by a driver amount, turning any
 * static palette into a continuously hue-cycling one.
 */
struct HueSpinShade {
  const float *amount; /**< Rotation driver in turns (0..1 = full turn). */
  /**
   * @brief Per-instance memo of the rotation folded into a cbrt-LMS 3x3.
   * @details *amount is frame-constant, so the matrix is rebuilt once per
   * frame; the per-sample cost is three fast_cbrt plus the folded transform.
   * mutable so const shade() can update the memo.
   */
  mutable float matrix[9] = {};
  mutable float cached_amount =
      0.0f;                    /**< Driver value the memo was built at. */
  mutable bool primed = false; /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory rotation driver.
   * @param amount Pointer to the per-frame rotation in turns; must not be null.
   */
  HueSpinShade(const float *amount) : amount(amount) {
    HS_CHECK(amount, "HueSpinShade: amount driver must not be null");
  }

  /**
   * @brief Rotates the sample's hue by the driver amount, preserving alpha.
   * @param c Sample color to reshape.
   * @param t Unused; the rotation is uniform over the domain.
   * @return The hue-rotated sample.
   */
  Color4 shade(Color4 c, float t) const {
    (void)t;
    if (!primed || *amount != cached_amount) {
      cached_amount = *amount;
      float ca, sa;
      turn_to_unit_cos_sin(cached_amount, ca, sa);
      hue_rotate_lms_matrix(ca, sa, matrix);
      primed = true;
    }
    LinRGB rgb = pixel_to_linrgb(c.color);
    LMS lms = linear_rgb_to_lms(rgb.r, rgb.g, rgb.b);
    lms_cbrt_transform_rgb(matrix, fast_cbrt(lms.l), fast_cbrt(lms.m),
                           fast_cbrt(lms.s), rgb.r, rgb.g, rgb.b);
    c.color = Pixel(float_to_pixel16(rgb.r), float_to_pixel16(rgb.g),
                    float_to_pixel16(rgb.b));
    return c;
  }
};

/**
 * @brief Rotates hue by an amount that varies along the palette domain, so
 * different parts of the gradient drift in opposite directions (iridescence).
 * @details Builds a rotation per sample; suited to bake-time sampling
 * (BakedPalette::rebake) rather than tight per-pixel loops.
 */
struct HueWobbleShade {
  const float *phase;
  float frequency;
  float depth;

  /**
   * @brief Constructs with a mandatory phase driver, frequency, and depth.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param freq Wobble frequency over the domain; defaults to 1.0.
   * @param depth Peak hue rotation in turns; defaults to 0.1.
   */
  HueWobbleShade(const float *phase, float freq = 1.0f, float depth = 0.1f)
      : phase(phase), frequency(freq), depth(depth) {
    HS_CHECK(phase, "HueWobbleShade: phase driver must not be null");
  }

  /**
   * @brief Rotates the sample's hue by depth * sin(t*frequency*2pi + phase)
   * turns, preserving alpha.
   * @param c Sample color to reshape.
   * @param t Coordinate driving the wobble.
   * @return The hue-rotated sample.
   */
  Color4 shade(Color4 c, float t) const {
    return hue_rotate(c,
                      depth * fast_sinf(t * frequency * PI_F * 2.0f + *phase));
  }
};

/**
 * @brief Ignites sparse traveling glints: where an evolving noise field over
 * the domain exceeds a threshold, the sample lerps toward white.
 */
struct SparkleShade {
  const float *time;
  float frequency;
  float threshold;
  uint32_t seed;

  /**
   * @brief Constructs with a mandatory time driver, density, and threshold.
   * @param time Pointer to the per-frame noise time axis; must not be null.
   * @param freq Glint density over the domain; defaults to 24.0.
   * @param threshold Noise level in [0, 1) above which a glint ignites; higher
   *        is sparser. Defaults to 0.75.
   * @param seed Noise stream selector; defaults to 0.
   */
  SparkleShade(const float *time, float freq = 24.0f, float threshold = 0.75f,
               uint32_t seed = 0)
      : time(time), frequency(freq), threshold(threshold), seed(seed) {
    HS_CHECK(time, "SparkleShade: time driver must not be null");
    HS_CHECK(threshold >= 0.0f && threshold < 1.0f,
             "SparkleShade: threshold must be in [0, 1)");
  }

  /**
   * @brief Whitens the sample where the noise field exceeds the threshold.
   * @param c Sample color to reshape.
   * @param t Coordinate locating the sample in the glint field.
   * @return The sample, lerped toward white by the over-threshold excess.
   */
  Color4 shade(Color4 c, float t) const {
    float n = value_noise_2d(t * frequency, *time, seed);
    if (n <= threshold)
      return c;
    float w = (n - threshold) / (1.0f - threshold);
    c.color = c.color.lerp16(Pixel(65535, 65535, 65535), frac_to_q16(w));
    return c;
  }
};

/**
 * @brief Breathes the palette's saturation: scales OKLab chroma by
 * 1 + depth * sin(phase), swinging every sample between pastel and vivid.
 */
struct ChromaPulseShade {
  const float *phase;
  float depth;
  /**
   * @brief Per-instance memo of the frame's chroma scale.
   * @details *phase is frame-constant, so the sine is recomputed once per
   * frame, not per sample. mutable so const shade() can update the memo.
   */
  mutable float cached_phase = 0.0f;
  mutable float cached_scale = 1.0f; /**< Memoized scale at cached_phase. */
  mutable bool primed = false; /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory phase driver and pulse depth.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param depth Pulse depth in [0, 1]: chroma swings over [1-depth, 1+depth].
   *        Defaults to 0.5.
   */
  ChromaPulseShade(const float *phase, float depth = 0.5f)
      : phase(phase), depth(depth) {
    HS_CHECK(phase, "ChromaPulseShade: phase driver must not be null");
    HS_CHECK(depth >= 0.0f && depth <= 1.0f,
             "ChromaPulseShade: depth must be in [0, 1]");
  }

  /**
   * @brief Scales the sample's OKLab chroma by the frame's pulse factor,
   * holding lightness and hue; over-gamut results chroma-clip.
   * @param c Sample color to reshape.
   * @param t Unused; the pulse is uniform over the domain.
   * @return The chroma-scaled sample, alpha untouched.
   */
  Color4 shade(Color4 c, float t) const {
    (void)t;
    if (!primed || *phase != cached_phase) {
      cached_phase = *phase;
      cached_scale = 1.0f + depth * fast_sinf(cached_phase);
      primed = true;
    }
    LinRGB rgb = pixel_to_linrgb(c.color);
    OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
    lab.a *= cached_scale;
    lab.b *= cached_scale;
    oklab_to_linear_rgb_gamut(lab, rgb.r, rgb.g, rgb.b);
    c.color = Pixel(float_to_pixel16(rgb.r), float_to_pixel16(rgb.g),
                    float_to_pixel16(rgb.b));
    return c;
  }
};

/**
 * @brief Grains the palette's brightness with an evolving noise field —
 * a subtler, hue-exact shimmer than SparkleShade's white glints.
 * @details Scales all three linear channels uniformly, so hue and saturation
 * ratios are exact below the saturation point; a gain above 1 clips bright
 * channels. No OKLab round-trip.
 */
struct LightnessGrainShade {
  const float *time;
  float frequency;
  float amplitude;
  uint32_t seed;

  /**
   * @brief Constructs with a mandatory time driver, grain density, and depth.
   * @param time Pointer to the per-frame noise time axis; must not be null.
   * @param freq Grain density over the domain; defaults to 12.0.
   * @param amp Gain swing in [0, 1]: brightness scales over [1-amp, 1+amp].
   *        Defaults to 0.25.
   * @param seed Noise stream selector; defaults to 0.
   */
  LightnessGrainShade(const float *time, float freq = 12.0f, float amp = 0.25f,
                      uint32_t seed = 0)
      : time(time), frequency(freq), amplitude(amp), seed(seed) {
    HS_CHECK(time, "LightnessGrainShade: time driver must not be null");
    HS_CHECK(amplitude >= 0.0f && amplitude <= 1.0f,
             "LightnessGrainShade: amplitude must be in [0, 1]");
  }

  /**
   * @brief Scales the sample's brightness by the local noise gain.
   * @param c Sample color to reshape.
   * @param t Coordinate locating the sample in the grain field.
   * @return The gain-scaled sample, alpha untouched.
   */
  Color4 shade(Color4 c, float t) const {
    float n = value_noise_2d(t * frequency, *time, seed);
    c.color = c.color * (1.0f + amplitude * (2.0f * n - 1.0f));
    return c;
  }
};

/**
 * @brief Adds a thin-film sheen: a phase-offset cosine overlay (the
 * ProceduralPalette waveform with per-channel thirds offsets) blended
 * additively over the sample, saturating at white.
 */
struct IridescentShade {
  const float *phase;
  float frequency;
  float weight;

  /**
   * @brief Constructs with a mandatory phase driver, frequency, and weight.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param freq Sheen frequency over the domain; defaults to 3.0.
   * @param weight Overlay strength; defaults to 0.25.
   */
  IridescentShade(const float *phase, float freq = 3.0f, float weight = 0.25f)
      : phase(phase), frequency(freq), weight(weight) {
    HS_CHECK(phase, "IridescentShade: phase driver must not be null");
  }

  /**
   * @brief Adds the weighted cosine sheen to the sample.
   * @param c Sample color to reshape.
   * @param t Coordinate locating the sample along the sheen.
   * @return The sample plus the overlay (per-channel saturating), alpha
   *         untouched.
   */
  Color4 shade(Color4 c, float t) const {
    float arg = t * frequency * PI_F * 2.0f + *phase;
    constexpr float THIRD = 2.0f * PI_F / 3.0f;
    Pixel sheen(
        srgb_to_linear_interp(0.5f + 0.5f * fast_cosf(arg)),
        srgb_to_linear_interp(0.5f + 0.5f * fast_cosf(arg + THIRD)),
        srgb_to_linear_interp(0.5f + 0.5f * fast_cosf(arg + 2.0f * THIRD)));
    c.color += sheen * weight;
    return c;
  }
};

/**
 * @brief Scales alpha by a caller-supplied falloff curve over the coordinate.
 */
struct AlphaFalloffShade {
  using FalloffFunction = float (*)(float);
  FalloffFunction fn;
  /**
   * @brief Constructs with the falloff function.
   * @param fn Non-null function mapping a coordinate to an alpha multiplier.
   */
  AlphaFalloffShade(FalloffFunction fn) : fn(fn) {
    HS_CHECK(fn != nullptr,
             "AlphaFalloffShade: falloff function must not be null");
  }
  /**
   * @brief Scales the sample's alpha by the falloff curve at t.
   * @param c Sample color to reshape.
   * @param t Coordinate passed to the falloff function.
   * @return The sample with alpha scaled.
   */
  Color4 shade(Color4 c, float t) const {
    c.alpha *= fn(t);
    return c;
  }
};

/**
 * @brief Fades the sample color to black near the coordinate edges.
 * @details Opaque vignette. Pair with InsetModifier so the edge bands resolve
 * to the source's first/last stop before fading.
 */
struct EdgeFadeShade {
  float edge;
  /**
   * @brief Constructs with the edge fade width.
   * @param edge Fraction of the domain over which each edge fades; default 0.2.
   */
  EdgeFadeShade(float edge = 0.2f) : edge(edge) {
    HS_CHECK(edge > 0.0f && edge <= 0.5f,
             "EdgeFadeShade: edge must be in (0, 0.5]");
  }
  /**
   * @brief Fades the sample color toward black within the edge bands.
   * @param c Sample color to reshape.
   * @param t Coordinate in [0, 1].
   * @return The sample with its color faded near the edges.
   */
  Color4 shade(Color4 c, float t) const {
    // 16-bit linear black, not CRGB: a CRGB temporary would route the blend
    // through an 8-bit sRGB lerp and band the fade.
    Pixel black(0, 0, 0);
    if (t < edge)
      return Color4(
          black.lerp16(c.color, frac_to_q16(quintic_kernel(t / edge))),
          c.alpha);
    if (t >= 1.0f - edge)
      return Color4(c.color.lerp16(black, frac_to_q16(quintic_kernel(
                                              (t - (1.0f - edge)) / edge))),
                    c.alpha);
    return c;
  }
};

/**
 * @brief Fades the sample alpha (not color) near the coordinate edges.
 * @details Transparent vignette. Pair with InsetModifier as with EdgeFadeShade.
 */
struct EdgeAlphaShade {
  float edge;
  /**
   * @brief Constructs with the edge fade width.
   * @param edge Fraction of the domain over which each edge fades; default 0.2.
   */
  EdgeAlphaShade(float edge = 0.2f) : edge(edge) {
    HS_CHECK(edge > 0.0f && edge <= 0.5f,
             "EdgeAlphaShade: edge must be in (0, 0.5]");
  }
  /**
   * @brief Fades the sample alpha within the edge bands.
   * @param c Sample color to reshape.
   * @param t Coordinate in [0, 1].
   * @return The sample with its alpha faded near the edges.
   */
  Color4 shade(Color4 c, float t) const {
    if (t < edge)
      c.alpha *= quintic_kernel(t / edge);
    else if (t >= 1.0f - edge)
      c.alpha *= quintic_kernel(1.0f - (t - (1.0f - edge)) / edge);
    return c;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Compile-Time Palette Composition
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Concept for a coordinate modifier.
 * @tparam T Type required to expose a const modify(float)->float method.
 * @details Remaps the lookup coordinate before the source is sampled.
 */
template <typename T>
concept CoordMod = requires(const T m, float t) {
  { m.modify(t) } -> std::convertible_to<float>;
};

/**
 * @brief Concept for a color modifier.
 * @tparam T Type required to expose a const shade(Color4, float)->Color4 method.
 * @details Reshapes the sample after the lookup, with the original coordinate
 * in hand.
 */
template <typename T>
concept ColorMod = requires(const T m, Color4 c, float t) {
  { m.shade(c, t) } -> std::convertible_to<Color4>;
};

/**
 * @brief Whether a coordinate modifier's output may leave [0,1].
 * @tparam M Coordinate modifier type.
 * @return M::requires_wrap when declared, false otherwise.
 */
template <typename M> constexpr bool coord_requires_wrap() {
  if constexpr (requires { M::requires_wrap; })
    return M::requires_wrap;
  else
    return false;
}

/**
 * @brief Whether a coordinate modifier maps [0,1] into [0,1] and reaches 1.
 * @tparam M Coordinate modifier type.
 * @return M::bounded_output when declared, false otherwise.
 */
template <typename M> constexpr bool coord_bounded_output() {
  if constexpr (requires { M::bounded_output; })
    return M::bounded_output;
  else
    return false;
}

/**
 * @brief Whether a coordinate modifier maps any input, in range or not, into
 * [0,1].
 * @tparam M Coordinate modifier type.
 * @return M::rebounds_input when declared, false otherwise.
 * @details Stronger than bounded_output, which only describes in-range input:
 * ReverseModifier and MirrorModifier are bounded on [0,1] but carry an
 * out-of-range coordinate straight through.
 */
template <typename M> constexpr bool coord_rebounds_input() {
  if constexpr (requires { M::rebounds_input; })
    return M::rebounds_input;
  else
    return false;
}

/**
 * @brief Whether the coordinate reaching the source may leave [0,1].
 * @tparam M Coordinate modifier types, in application order.
 * @return True when an unbounded modifier is not re-bounded by a later one.
 */
template <typename... M> constexpr bool coord_chain_leaves_unit() {
  bool leaves = false;
  ((leaves = coord_rebounds_input<M>() ? false
                                       : (leaves || coord_requires_wrap<M>())),
   ...);
  return leaves;
}

/**
 * @brief Whether the last modifier in the chain has bounded output.
 * @tparam M Coordinate modifier types, in application order.
 * @return coord_bounded_output of the final entry; false for an empty chain.
 */
template <typename... M> constexpr bool coord_chain_bounded_tail() {
  bool bounded = false;
  ((bounded = coord_bounded_output<M>()), ...);
  return bounded;
}

/**
 * @brief Type-list tag for the coordinate-modifier axis of a StaticPalette.
 * @tparam M Coordinate modifier types.
 */
template <typename... M> struct Coords {};
/**
 * @brief Type-list tag for the color-modifier axis of a StaticPalette.
 * @tparam M Color modifier types.
 */
template <typename... M> struct Colors {};

/**
 * @brief A compile-time composition of a Source palette, a coordinate-modifier
 * chain, and a color-modifier chain.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CoordList Coords<> type-list of coordinate modifiers.
 * @tparam ColorList Colors<> type-list of color modifiers.
 * @tparam Wrap Selects two linked behaviours: whether the coordinate is folded
 * into [0,1) before the source lookup, and which coordinate the color-modifier
 * chain receives — the folded post-coord-modifier value when true, the raw
 * pre-modifier input when false. Flipping it for the folding alone also changes
 * what every shade() sees.
 * @details Default construct, then bind() (ArenaVector idiom); both chains are
 * inlined by fold expression. get() applies the coord mods to t in order,
 * samples the source (wrapping the coordinate unless Wrap is false), then
 * applies the color mods with the wrapped coordinate when Wrap is true, or the
 * original coordinate when Wrap is false. Wrap=false suits inset/falloff
 * pipelines that must reach the source's exact endpoints
 * (wrap_t(1)==0 would otherwise fold the top edge). Both settings are checked
 * at compile time: `requires_wrap` on any unbounded modifier rejects Wrap=false
 * unless a later modifier declares `rebounds_input`, and `bounded_output` on the
 * final coord modifier rejects Wrap=true.
 */
template <typename Source, typename CoordList = Coords<>,
          typename ColorList = Colors<>, bool Wrap = true>
class StaticPalette;

/**
 * @brief Partial specialization splitting the two modifier type-lists.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CMods Coordinate modifier types.
 * @tparam XMods Color modifier types.
 * @tparam Wrap Folds the coordinate before the source lookup, and selects the
 * coordinate handed to the color chain (folded value, or raw input).
 */
template <typename Source, typename... CMods, typename... XMods, bool Wrap>
class StaticPalette<Source, Coords<CMods...>, Colors<XMods...>, Wrap> {
  static_assert((CoordMod<CMods> && ...), "Coords<> entries must be CoordMods");
  static_assert((ColorMod<XMods> && ...), "Colors<> entries must be ColorMods");
  static_assert(Wrap || !coord_chain_leaves_unit<CMods...>(),
                "Wrap=false composed with a coordinate modifier that leaves "
                "[0,1] (requires_wrap, e.g. CycleModifier/ScaleModifier): the "
                "source would be sampled out of range and the palette would "
                "freeze at its endpoint. Use Wrap=true, or follow it with a "
                "modifier that re-bounds arbitrary input (rebounds_input, e.g. "
                "FoldModifier/InsetModifier).");
  static_assert(!Wrap || !coord_chain_bounded_tail<CMods...>(),
                "Wrap=true with a bounded final coordinate modifier "
                "(bounded_output, e.g. MirrorModifier/InsetModifier): wrap_t "
                "folds its 1.0 output to 0.0, destroying the top endpoint. "
                "Use Wrap=false.");

public:
  static constexpr bool WRAPS_COORDINATE = Wrap;

  /**
   * @brief Default-constructs an unbound composition (bind() before use).
   */
  StaticPalette() = default;

  /**
   * @brief Binds the source and modifier chains by pointer.
   * @param src Source palette; must not be null.
   * @param cms Coordinate-modifier pointers, one per CMods entry; none null.
   * @param xms Color-modifier pointers, one per XMods entry; none null.
   * @details get()'s source assert is stripped on-device and a null read does
   * not fault on Teensy 4.x, so null binds are trapped here (always-on HS_CHECK,
   * empty packs fold to true) at the cold init seam.
   */
  void bind(const Source *src, const CMods *...cms, const XMods *...xms) {
    HS_CHECK(src != nullptr, "StaticPalette bound to null source");
    HS_CHECK(((cms != nullptr) && ...),
             "StaticPalette bound to null coord modifier");
    HS_CHECK(((xms != nullptr) && ...),
             "StaticPalette bound to null color modifier");
    source = src;
    coords = std::make_tuple(cms...);
    colors = std::make_tuple(xms...);
  }

  /**
   * @brief Applies the coord chain, samples the source, then the color chain.
   * @param t Lookup coordinate.
   * @return The fully modified color.
   * @details The coord mods remap t in order; the source is sampled (wrapping
   * the coordinate unless Wrap is false); then the color mods reshape the
   * sample with the wrapped coordinate, or the original coordinate when
   * Wrap=false.
   */
  Color4 get(float t) const {
    assert(source != nullptr && "StaticPalette used before bind()!");

    float ft = t;
    std::apply([&](const auto *...m) { ((ft = m->modify(ft)), ...); }, coords);

    float u = ft;
    if constexpr (Wrap)
      u = wrap_t(ft);
    Color4 c = source->get(u);

    const float shade_coordinate = Wrap ? u : t;
    std::apply(
        [&](const auto *...m) { ((c = m->shade(c, shade_coordinate)), ...); },
        colors);
    return c;
  }

private:
  const Source *source = nullptr;
  std::tuple<const CMods *...> coords{};
  std::tuple<const XMods *...> colors{};
};

/**
 * @brief Runtime Palette facade over a compile-time StaticPalette composition.
 * @tparam SP StaticPalette composition type exposing Color4 get(float) const.
 * @details Bridges a zero-overhead StaticPalette into the polymorphic
 * `const Palette*` world (preset tables, BakedPalette::bake). The virtual call
 * is paid only at bake time (cold), never on the per-pixel path.
 */
template <typename SP> class PaletteFacade : public Palette {
public:
  /**
   * @brief Default-constructs an unbound facade (bind() before use).
   */
  PaletteFacade() = default;
  /**
   * @brief Constructs a facade bound to a composition.
   * @param sp Composition to forward get() to.
   */
  explicit PaletteFacade(const SP *sp) : composition(sp) {}
  /**
   * @brief Binds the facade to a composition.
   * @param sp Composition to forward get() to; must not be null.
   */
  void bind(const SP *sp) {
    HS_CHECK(sp != nullptr, "PaletteFacade bound to null composition");
    composition = sp;
  }
  /**
   * @brief Forwards the lookup to the bound composition.
   * @param t Lookup coordinate.
   * @return The composition's color at t.
   */
  Color4 get(float t) const override {
    assert(composition != nullptr && "PaletteFacade used before bind()!");
    return composition->get(t);
  }

private:
  const SP *composition = nullptr;
};

/**
 * @brief Palette that returns one fixed color for every coordinate.
 */
class SolidColorPalette : public Palette {
public:
  /**
   * @brief Constructs with the fixed color.
   * @param color Color returned for every lookup.
   */
  SolidColorPalette(const Color4 &color) : color(color) {}
  /**
   * @brief Returns the fixed color.
   * @return The stored color, regardless of coordinate.
   */
  Color4 get(float) const override { return color; }

private:
  Color4 color;
};

/**
 * @brief Pre-baked 256-entry color/alpha LUT allocated in an arena.
 * @details Stores parallel Pixel and Q16-alpha tables with lerp
 * interpolation. Not a Palette subclass — call get(t) directly for
 * zero-overhead lookups.
 */
class BakedPalette {
public:
  static constexpr int LUT_SIZE = 256;

  /**
   * @brief Arena bytes bake() consumes, including worst-case alignment padding.
   */
  static constexpr size_t required_arena_bytes() {
    return LUT_SIZE * (sizeof(Pixel) + sizeof(uint16_t)) + alignof(Pixel) +
           alignof(uint16_t);
  }

  /**
   * @brief Default-constructs an unbaked palette (bake() before use).
   */
  BakedPalette() = default;

  /**
   * @brief Bakes any source into a 256-entry LUT in the given arena.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param arena Arena to allocate the LUT from.
   * @param source Source palette or composition to sample.
   * @details Works for a runtime Palette or a compile-time StaticPalette alike.
   */
  template <typename Source>
  HS_COLD_MEMBER void bake(Arena &arena, const Source &source) {
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    rebake(source);
  }

  /**
   * @brief Refills the existing LUT without allocating. Use for animated palettes.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param source Source palette or composition to sample.
   * @details Entry i samples t = i / (LUT_SIZE - 1), so the last entry lands on
   * t = 1 exactly. A composition with Wrap=true folds that sample back to 0 and
   * collapses its last entry onto its first — bake such sources with Wrap=false.
   */
  template <typename Source> HS_COLD_MEMBER void rebake(const Source &source) {
    if constexpr (requires { Source::WRAPS_COORDINATE; }) {
      static_assert(!Source::WRAPS_COORDINATE,
                    "BakedPalette cannot rebake a wrapping source");
    }
    HS_CHECK(colors != nullptr && alpha_q16 != nullptr,
             "BakedPalette::rebake before bake()");
    for (int i = 0; i < LUT_SIZE; ++i) {
      float t = static_cast<float>(i) / (LUT_SIZE - 1);
      const Color4 sample = source.get(t);
      colors[i] = sample.color;
      alpha_q16[i] = frac_to_q16(sample.alpha);
    }
  }

  /**
   * @brief Bakes this LUT as the w-blend of two baked palettes.
   * @param arena Arena the blended LUT is allocated from.
   * @param from The w = 0 endpoint; must be baked.
   * @param to The w = 1 endpoint; must be baked.
   * @param w Blend weight in (0, 1).
   * @details Walks the two source LUTs entry-wise with the fixed-point channel
   * lerp — no per-entry float resampling.
   */
  HS_COLD_MEMBER void bake_blend(Arena &arena, const BakedPalette &from,
                                 const BakedPalette &to, float w) {
    HS_CHECK(from.colors && from.alpha_q16 && to.colors && to.alpha_q16,
             "BakedPalette::bake_blend before bake()");
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    // Clamp before the cast: w < 0 or NaN is float->int UB, and a NaN weight
    // would otherwise reach every entry's alpha.
    const float wc = hs::clamp(w, 0.0f, 1.0f);
    const uint16_t weight = frac_to_q16(wc);
    for (int i = 0; i < LUT_SIZE; ++i) {
      colors[i] = from.colors[i].lerp16(to.colors[i], weight);
      alpha_q16[i] = lerp_q16(from.alpha_q16[i], to.alpha_q16[i], weight);
    }
  }

  /**
   * @brief Fast lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to the last entry).
   * @return The interpolated color.
   */
  Color4 get(float t) const {
    Color4 out;
    sample_into(t, out);
    return out;
  }

  /**
   * @brief Samples only the interpolated RGB channels.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The pixel get() gives for the same index, without interpolating
   * alpha.
   */
  __attribute__((always_inline)) Pixel get_color(float t) const {
    assert(colors != nullptr && "BakedPalette::get_color before bake()");
    float idx =
        hs::clamp(t * (LUT_SIZE - 1), 0.0f, static_cast<float>(LUT_SIZE - 1));
    return sample_color_index(idx);
  }

  /**
   * @brief Samples RGB for a coordinate already clamped to [0, 1].
   * @param t Finite lookup coordinate in [0, 1].
   * @return The pixel get() gives for the same index.
   */
  __attribute__((always_inline)) Pixel get_color_unit(float t) const {
    assert(colors != nullptr && "BakedPalette::get_color_unit before bake()");
    // Also traps NaN: both comparisons fail.
    assert(t >= 0.0f && t <= 1.0f);
    return sample_color_index(t * (LUT_SIZE - 1));
  }

  /**
   * @brief Deep-copies the LUT from another BakedPalette into the given arena.
   * @param src Source palette to copy; must already be baked.
   * @param arena Arena to allocate the new LUT from.
   * @details Used by Persist for arena compaction.
   */
  void clone_from(const BakedPalette &src, Arena &arena) {
    HS_CHECK(src.colors != nullptr && src.alpha_q16 != nullptr,
             "BakedPalette::clone_from before src bake()");
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    memcpy(colors, src.colors, LUT_SIZE * sizeof(Pixel));
    memcpy(alpha_q16, src.alpha_q16, LUT_SIZE * sizeof(uint16_t));
  }

private:
  static __attribute__((always_inline)) uint16_t lerp_q16(uint16_t a,
                                                          uint16_t b,
                                                          uint16_t weight) {
    const uint32_t inverse = 65535u - weight;
    const uint32_t x =
        static_cast<uint32_t>(a) * inverse + static_cast<uint32_t>(b) * weight;
    return static_cast<uint16_t>((x + (x >> 16) + 32768u) >> 16);
  }

  __attribute__((always_inline)) Pixel sample_color_index(float idx) const {
    if (idx <= 0.0f)
      return colors[0];
    return lut_sample_pixel(colors, LUT_SIZE, idx);
  }

  __attribute__((always_inline)) void sample_into(float t, Color4 &out) const {
    assert(colors != nullptr && alpha_q16 != nullptr &&
           "BakedPalette::get before bake()");
    // Clamp before the int cast: static_cast<int>(NaN) is UB. hs::clamp maps NaN
    // to the hi bound (last entry) and guarantees idx >= 0.
    float idx =
        hs::clamp(t * (LUT_SIZE - 1), 0.0f, static_cast<float>(LUT_SIZE - 1));
    // Split and weight through the same helpers lut_sample_pixel uses, so this
    // path and get_color/get_color_unit share one spelling of the arithmetic.
    const int lo = lut_index_lo(idx);
    if (lo >= LUT_SIZE - 1) {
      out = Color4(colors[LUT_SIZE - 1],
                   alpha_q16[LUT_SIZE - 1] * (1.0f / 65535.0f));
      return;
    }
    const uint16_t weight = lut_index_weight(idx, lo);
    out = Color4(colors[lo].lerp16(colors[lo + 1], weight),
                 lerp_q16(alpha_q16[lo], alpha_q16[lo + 1], weight) *
                     (1.0f / 65535.0f));
  }
  Pixel *colors = nullptr;
  uint16_t *alpha_q16 = nullptr;
};

/**
 * @brief Bake-time adapter mapping a LUT coordinate from the cos domain into a
 *        source palette's angle parameter.
 * @tparam Source Type exposing Color4 get(float) const over t = angle/PI.
 * @details Baking through this folds the d -> acos(d)/PI radial mapping into the
 * bake (256 acos per bake, not one per fragment): the fragment lookup keys the
 * LUT by the raw dot product via dot_key(d). dot_key inverts this mapping:
 * u -> d = 1 - 2u, get(u) returns the source at acos(d)/PI.
 */
template <typename Source> struct DotKeyed {
  const Source &source;
  Color4 get(float u) const {
    float d = hs::clamp(1.0f - 2.0f * u, -1.0f, 1.0f);
    return source.get(fast_acos(d) / PI_F);
  }
};

/**
 * @brief Wraps a palette source for a DotKeyed bake.
 * @param source Palette sampled over t = angle/PI; must outlive the bake call.
 */
template <typename Source>
inline DotKeyed<Source> dot_keyed(const Source &source) {
  return DotKeyed<Source>{source};
}

/// LUT coordinate for cos-value d = dot(axis, v); inverse of DotKeyed's mapping.
inline float dot_key(float d) {
  return (1.0f - hs::clamp(d, -1.0f, 1.0f)) * 0.5f;
}

/**
 * @brief Resolves a (from, to) baked-LUT pair at one crossfade weight.
 * @param dst Receives the resolved palette (a default-constructed unbaked
 * instance is fine).
 * @param arena Arena receiving the blended LUT when one is baked.
 * @param from The w = 0 endpoint.
 * @param to The w = 1 endpoint.
 * @param w Blend weight.
 * @details Weights at or beyond an endpoint alias that endpoint's LUT storage
 * (bitwise-exact, no allocation), so crossfade boundaries are exact by
 * construction; only 0 < w < 1 bakes a blended LUT into the arena.
 */
inline void bake_palette_blend(BakedPalette &dst, Arena &arena,
                               const BakedPalette &from, const BakedPalette &to,
                               float w) {
  if (w <= 0.0f)
    dst = from;
  else if (w >= 1.0f)
    dst = to;
  else
    dst.bake_blend(arena, from, to, w);
}

/**
 * @brief Steps a palette cross-fade rebake for one frame.
 * @tparam Source Palette type exposing Color4 get(float) const.
 * @param wipe_pending Set when a wipe was just armed; consumes the arming frame.
 * @param wipe_frames_remaining Frames left to rebake; decremented per step.
 * @param baked LUT rebaked from the wipe-mutated source while the wipe runs.
 * @param source Palette mutated in place by the in-flight ColorWipe.
 * @details A ColorWipe is armed mid-step and first steps next frame, so the
 * arming frame is skipped; each later frame rebakes the LUT the shader samples.
 */
template <typename Source>
inline void step_wipe_rebake(bool &wipe_pending, int &wipe_frames_remaining,
                             BakedPalette &baked, const Source &source) {
  if (wipe_pending) {
    wipe_pending = false;
  } else if (wipe_frames_remaining > 0) {
    baked.rebake(source);
    --wipe_frames_remaining;
  }
}

/**
 * @brief Bank of N baked palettes for bulk Persist/clone operations.
 */
struct BakedPaletteBank {
  static constexpr int N = 6;
  BakedPalette entries[N];

  /**
   * @brief Deep-copies all entries into a target arena.
   * @param src Source bank to copy from.
   * @param dst Destination bank to fill.
   * @param arena Arena to allocate the cloned LUTs from.
   * @details Required by Cloneable.
   */
  HS_COLD_MEMBER static void clone(const BakedPaletteBank &src,
                                   BakedPaletteBank &dst, Arena &arena) {
    for (int i = 0; i < N; ++i)
      dst.entries[i].clone_from(src.entries[i], arena);
  }
};
