/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Palette-composition layer for color.h: the coordinate/color modifiers, the
// StaticPalette composition template, and baked-palette storage. Included by
// color.h after the palette core; not a standalone header.

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
   *       cycling. Composing this modifier with a `Wrap=false` palette samples
   *       the source out of range — keep `Wrap=true` whenever a CycleModifier
   *       drives the coordinate.
   */
  float modify(float t) const { return offset ? t + *offset : t; }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier {
  const float *phase;
  float amplitude;
  /**
   * @brief Per-instance memo of fast_sinf(*phase).
   * @details *phase is frame-constant, so the sine is recomputed once per frame,
   * not per pixel. mutable so const modify() can update the memo.
   */
  mutable float cached_phase_ = 0.0f;
  mutable float cached_sin_ = 0.0f;   /**< Memoized sine of cached_phase_. */
  mutable bool primed_ = false;       /**< Whether the memo has been populated. */

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
    if (!primed_ || *phase != cached_phase_) {
      cached_phase_ = *phase;
      cached_sin_ = fast_sinf(*phase);
      primed_ = true;
    }
    return t + cached_sin_ * amplitude;
  }
};

/**
 * @brief Distorts the palette spatially with a sine wave, creating a liquid
 * ripple effect. Compresses and expands colors like waves on a spatial coord.
 */
struct RippleModifier {
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
    return t +
           (value_noise_2d(t * frequency, *time, seed) - 0.5f) * 2.0f * amplitude;
  }
};

/**
 * @brief Meanders the whole palette along a smooth noise walk — unlike
 * CycleModifier's linear scroll, the offset wanders, hesitates, and reverses.
 */
struct DriftModifier {
  const float *time;
  float speed;
  float amplitude;
  uint32_t seed;
  /**
   * @brief Per-instance memo of the frame's walk offset.
   * @details *time is frame-constant, so the noise walk is sampled once per
   * frame, not per pixel. mutable so const modify() can update the memo.
   */
  mutable float cached_time = 0.0f;
  mutable float cached_offset = 0.0f; /**< Memoized offset at cached_time. */
  mutable bool primed = false;        /**< Whether the memo has been populated. */

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
      cached_offset =
          (value_noise_1d(cached_time * speed, seed) - 0.5f) * 2.0f * amplitude;
      primed = true;
    }
    return t + cached_offset;
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
    if (m < 0.0f) m += 2.0f;
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
  const float *tension; /**< Pinch tension driver; expects roughly -0.9 to 0.9. */

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
   * @brief Multiplies the coordinate by the active scale.
   * @param t Input coordinate.
   * @return The scaled coordinate.
   * @note For scale > 1 the result intentionally leaves [0,1], relying on the
   *       consuming palette's `Wrap=true` (the `StaticPalette` default) to fold
   *       it back into range — that fold IS the multiple-repeats effect.
   *       Composing this modifier with a `Wrap=false` palette samples the source
   *       out of range — keep `Wrap=true` whenever a ScaleModifier with scale > 1
   *       drives the coordinate.
   */
  float modify(float t) const {
    return t * (dynamic_scale ? *dynamic_scale : base_scale);
  }
};

/**
 * @brief Reverses the palette coordinate (t -> 1 - t).
 */
struct ReverseModifier {
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
  float lo, hi;
  /**
   * @brief Constructs the inset window bounds.
   * @param lo Lower domain bound mapped to 0; defaults to 0.2.
   * @param hi Upper domain bound mapped to 1; defaults to 0.8.
   */
  InsetModifier(float lo = 0.2f, float hi = 0.8f) : lo(lo), hi(hi) {
    HS_CHECK(hi > lo, "InsetModifier: hi must be > lo (modify divides by hi - lo)");
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
  mutable float cached_amount = 0.0f; /**< Driver value the memo was built at. */
  mutable bool primed = false;        /**< Whether the memo has been populated. */

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
      float angle = cached_amount * (2.0f * PI_F);
      float ca = fast_cosf(angle);
      float sa = fast_sinf(angle);
      // renormalize fast trig so the rotation preserves chroma
      float inv = 1.0f / sqrtf(ca * ca + sa * sa);
      hue_rotate_lms_matrix(ca * inv, sa * inv, matrix);
      primed = true;
    }
    constexpr float INV16 = 1.0f / 65535.0f;
    float r = c.color.r * INV16, g = c.color.g * INV16, b = c.color.b * INV16;
    LMS lms = linear_rgb_to_lms(r, g, b);
    lms_cbrt_transform_rgb(matrix, fast_cbrt(lms.l), fast_cbrt(lms.m),
                           fast_cbrt(lms.s), r, g, b);
    c.color = Pixel(float_to_pixel16(r), float_to_pixel16(g),
                    float_to_pixel16(b));
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
  mutable bool primed = false;       /**< Whether the memo has been populated. */

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
    constexpr float INV16 = 1.0f / 65535.0f;
    float r = c.color.r * INV16, g = c.color.g * INV16, b = c.color.b * INV16;
    LMS lms = linear_rgb_to_lms(r, g, b);
    OKLab lab =
        lms_to_oklab(fast_cbrt(lms.l), fast_cbrt(lms.m), fast_cbrt(lms.s));
    lab.a *= cached_scale;
    lab.b *= cached_scale;
    oklab_to_linear_rgb_gamut(lab, r, g, b);
    c.color = Pixel(float_to_pixel16(r), float_to_pixel16(g),
                    float_to_pixel16(b));
    return c;
  }
};

/**
 * @brief Grains the palette's brightness with an evolving noise field —
 * a subtler, hue-exact shimmer than SparkleShade's white glints.
 * @details Scales all three linear channels uniformly, so hue and saturation
 * ratios are preserved exactly; no OKLab round-trip.
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
    Pixel sheen(srgb_to_linear_interp(0.5f + 0.5f * fast_cosf(arg)),
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
          black.lerp16(c.color, float_to_pixel16(quintic_kernel(t / edge))),
          c.alpha);
    if (t >= 1.0f - edge)
      return Color4(c.color.lerp16(black, float_to_pixel16(quintic_kernel(
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
 * @tparam Wrap Whether to wrap the coordinate before sampling the source.
 * @details Default construct, then bind() (ArenaVector idiom); both chains are
 * inlined by fold expression. get() applies the coord mods to t in order,
 * samples the source (wrapping the coordinate unless Wrap is false), then
 * applies the color mods with the *original* coordinate. Wrap=false suits
 * inset/falloff pipelines that must reach the source's exact endpoints
 * (wrap_t(1)==0 would otherwise fold the top edge).
 */
template <typename Source, typename CoordList = Coords<>,
          typename ColorList = Colors<>, bool Wrap = true>
class StaticPalette;

/**
 * @brief Partial specialization splitting the two modifier type-lists.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CMods Coordinate modifier types.
 * @tparam XMods Color modifier types.
 * @tparam Wrap Whether to wrap the coordinate before sampling the source.
 */
template <typename Source, typename... CMods, typename... XMods, bool Wrap>
class StaticPalette<Source, Coords<CMods...>, Colors<XMods...>, Wrap> {
  static_assert((CoordMod<CMods> && ...), "Coords<> entries must be CoordMods");
  static_assert((ColorMod<XMods> && ...), "Colors<> entries must be ColorMods");

public:
  /**
   * @brief Default-constructs an unbound composition (bind() before use).
   */
  StaticPalette() = default;

  /**
   * @brief Binds the source and modifier chains by pointer.
   * @param src Source palette; must not be null.
   * @param cms Coordinate-modifier pointers, one per CMods entry; none null.
   * @param xms Color-modifier pointers, one per XMods entry; none null.
   * @details get()'s source_ assert is stripped on-device and a null read does
   * not fault on Teensy 4.x, so null binds are trapped here (always-on HS_CHECK,
   * empty packs fold to true) at the cold init seam.
   */
  void bind(const Source *src, const CMods *...cms, const XMods *...xms) {
    HS_CHECK(src != nullptr, "StaticPalette bound to null source");
    HS_CHECK(((cms != nullptr) && ...), "StaticPalette bound to null coord modifier");
    HS_CHECK(((xms != nullptr) && ...), "StaticPalette bound to null color modifier");
    source_ = src;
    coords_ = std::make_tuple(cms...);
    colors_ = std::make_tuple(xms...);
  }

  /**
   * @brief Applies the coord chain, samples the source, then the color chain.
   * @param t Lookup coordinate.
   * @return The fully modified color.
   * @details The coord mods remap t in order; the source is sampled (wrapping
   * the coordinate unless Wrap is false); then the color mods reshape the
   * sample with the *original* coordinate.
   */
  Color4 get(float t) const {
    assert(source_ != nullptr && "StaticPalette used before bind()!");

    float ft = t;
    std::apply([&](const auto *...m) { ((ft = m->modify(ft)), ...); }, coords_);

    float u = ft;
    if constexpr (Wrap)
      u = wrap_t(ft);
    Color4 c = source_->get(u);

    std::apply([&](const auto *...m) { ((c = m->shade(c, t)), ...); }, colors_);
    return c;
  }

private:
  const Source *source_ = nullptr;
  std::tuple<const CMods *...> coords_{};
  std::tuple<const XMods *...> colors_{};
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
  explicit PaletteFacade(const SP *sp) : sp_(sp) {}
  /**
   * @brief Binds the facade to a composition.
   * @param sp Composition to forward get() to; must not be null.
   */
  void bind(const SP *sp) {
    HS_CHECK(sp != nullptr, "PaletteFacade bound to null composition");
    sp_ = sp;
  }
  /**
   * @brief Forwards the lookup to the bound composition.
   * @param t Lookup coordinate.
   * @return The composition's color at t.
   */
  Color4 get(float t) const override {
    assert(sp_ != nullptr && "PaletteFacade used before bind()!");
    return sp_->get(t);
  }

private:
  const SP *sp_ = nullptr;
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
  Color4 color;
};

/**
 * @brief Pre-baked 256-entry Pixel16 LUT allocated in an arena.
 * @details Converts any Palette into a fast table lookup with lerp
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
    return LUT_SIZE * sizeof(Color4) + alignof(Color4);
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
  template <typename Source> HS_COLD_MEMBER void bake(Arena &arena, const Source &source) {
    lut_ = arena.allocate_n<Color4>(LUT_SIZE);
    rebake(source);
  }

  /**
   * @brief Refills the existing LUT without allocating. Use for animated palettes.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param source Source palette or composition to sample.
   */
  template <typename Source> HS_COLD_MEMBER void rebake(const Source &source) {
    HS_CHECK(lut_ != nullptr, "BakedPalette::rebake before bake()");
    for (int i = 0; i < LUT_SIZE; ++i) {
      float t = static_cast<float>(i) / (LUT_SIZE - 1);
      lut_[i] = source.get(t);
    }
  }

  /**
   * @brief Fast lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to the last entry).
   * @return The interpolated color.
   */
  Color4 get(float t) const {
    assert(lut_ != nullptr && "BakedPalette::get before bake()");
    // Clamp before the int cast: static_cast<int>(NaN) is UB. hs::clamp maps NaN
    // to the hi bound (last entry) and guarantees idx >= 0.
    float idx = hs::clamp(t * (LUT_SIZE - 1), 0.0f,
                          static_cast<float>(LUT_SIZE - 1));
    int lo = static_cast<int>(idx);
    if (lo >= LUT_SIZE - 1) return lut_[LUT_SIZE - 1];
    float frac = idx - lo;
    const Color4 &a = lut_[lo];
    const Color4 &b = lut_[lo + 1];
    return Color4(a.color.lerp16(b.color, frac_to_q16(frac)),
                  hs::clamp(a.alpha + (b.alpha - a.alpha) * frac, 0.0f, 1.0f));
  }

  /**
   * @brief Deep-copies the LUT from another BakedPalette into the given arena.
   * @param src Source palette to copy; must already be baked.
   * @param arena Arena to allocate the new LUT from.
   * @details Used by Persist for arena compaction.
   */
  void clone_from(const BakedPalette &src, Arena &arena) {
    HS_CHECK(src.lut_ != nullptr, "BakedPalette::clone_from before src bake()");
    lut_ = arena.allocate_n<Color4>(LUT_SIZE);
    memcpy(lut_, src.lut_, LUT_SIZE * sizeof(Color4));
  }

private:
  Color4 *lut_ = nullptr;
};

/**
 * @brief Blend source for a palette crossfade: two baked LUTs lerped at a
 * fixed weight.
 */
struct RampBlend {
  const BakedPalette &from; /**< The w = 0 endpoint. */
  const BakedPalette &to;   /**< The w = 1 endpoint. */
  float w;                  /**< Blend weight in [0, 1]. */
  /**
   * @brief Samples both LUTs and lerps.
   * @param t Lookup coordinate.
   * @return The blended color at t.
   */
  Color4 get(float t) const { return from.get(t).lerp(to.get(t), w); }
};

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
    dst.bake(arena, RampBlend{from, to, w});
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
  static constexpr int N = 5;
  BakedPalette entries[N];

  /**
   * @brief Deep-copies all entries into a target arena.
   * @param src Source bank to copy from.
   * @param dst Destination bank to fill.
   * @param arena Arena to allocate the cloned LUTs from.
   * @details Required by Cloneable.
   */
  HS_COLD_MEMBER static void clone(const BakedPaletteBank &src, BakedPaletteBank &dst,
                    Arena &arena) {
    for (int i = 0; i < N; ++i)
      dst.entries[i].clone_from(src.entries[i], arena);
  }
};
