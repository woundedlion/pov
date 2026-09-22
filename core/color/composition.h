/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file composition.h
 * @brief Coordinate and color modifiers and static palette composition.
 */

#include <tuple>
#include "color/palette.h"
#include "color/palette_recipe.h"
#include "color/color_space.h"
#include "math/geometry.h"

///////////////////////////////////////////////////////////////////////////////
// Palette Modifiers
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Soft-limit on the |phase| a caller may hand to the fast-trig
 * modifiers below.
 * @details fast_sinf/fast_cosf range-reduce by `x - floor(x/2pi)*2pi`, whose
 * error is one ULP of the argument: past this bound the reduction error crosses
 * ~5e-4 rad and the oscillation quantizes and drifts. A free-running per-frame
 * accumulator — the natural driver for these — passes it in hours, so the
 * driver, not the modifier, owns the wrap: fold the accumulator back into
 * [0, 2pi) at the source. Folding here instead would cost a per-sample floorf
 * and could not recover the precision the driver already lost. Each consuming
 * body carries a stripped assert on the bound, so an unfolded driver trips in
 * the host build and costs nothing on device; a driver in turns asserts against
 * its radian-scaled value.
 */
inline constexpr float PALETTE_PHASE_ARG_LIMIT = 4096.0f;

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
   * @note The result intentionally leaves [0,1] (t plus the offset), relying on
   *       the consuming palette's `Wrap=true` (the `StaticPalette` default) to
   *       fold it back into range and produce the cycling. `requires_wrap`
   *       makes a `Wrap=false` composition a compile error. The driver owns
   *       the bound: keep the offset in [0,1), which `Animation::Driver`'s
   *       default `wrap` does. A free-running accumulator loses the fractional
   *       bits the palette coordinate is made of.
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
   * per-pixel modify() can dereference unconditionally. The driver must keep
   * |phase| under PALETTE_PHASE_ARG_LIMIT.
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
    assert(fabsf(*phase) < PALETTE_PHASE_ARG_LIMIT);
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
   * construction rather than silently passing t through on every pixel. The
   * caller must keep |t * freq * 2pi + phase| under PALETTE_PHASE_ARG_LIMIT.
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
    const float arg = t * frequency * PI_F * 2.0f + *phase;
    assert(fabsf(arg) < PALETTE_PHASE_ARG_LIMIT);
    return t + fast_sinf(arg) * amplitude;
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

    // Triangle wave over a [0, 2) reduction: a negative scaled would otherwise
    // fold above 1.
    return fabsf(wrap(scaled, 2.0f) - 1.0f);
  }
};

/**
 * @brief Pinches or expands the center of the palette.
 * positive tension pulls colors toward the center, negative pushes them to the
 * edges.
 *
 * Null tension driver is the deliberate "no pinch" pass-through mode.
 * @details A bound driver costs a powf per sample; suited to bake-time sampling
 * (BakedPaletteStorage::rebake) rather than tight per-pixel loops.
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
  /** @brief In-range input stays in [0,1] and hits 1; palette needs Wrap=false.
   */
  static constexpr bool bounded_output = true;

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
   * @return The quantized coordinate, capped at 1 to hold bounded_output for a
   *   fractional step count.
   */
  float modify(float t) const {
    float s = dynamic_steps ? *dynamic_steps : base_steps;
    if (s < 1.0f)
      s = 1.0f;

    // Round to nearest step in the infinite domain.
    return __builtin_fminf(roundf(t * s) / s, 1.0f);
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
  float modify(float t) const { return unit_bell(t); }
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

/**
 * @brief Folds the coordinate into [0,1) mid-chain.
 * @details Placed between a `requires_wrap` modifier and a `bounded_output`
 * tail, it absorbs the out-of-range coordinate the tail would otherwise carry
 * through, so the composition can use Wrap=false and keep the tail's 1.0
 * endpoint ("scroll the palette, then mirror it").
 */
struct WrapModifier {
  /** @brief The fold confines any input, in range or not, to [0,1). */
  static constexpr bool rebounds_input = true;

  /**
   * @brief Folds the coordinate into the unit domain.
   * @param t Input coordinate.
   * @return t folded into [0, 1).
   */
  float modify(float t) const { return wrap_t(t); }
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
   * frame; the per-sample cost is one fast_cbrt3 plus the folded transform.
   * mutable so const shade() can update the memo.
   */
  mutable float matrix[9] = {};
  mutable float cached_amount =
      0.0f;                    /**< Driver value the memo was built at. */
  mutable bool primed = false; /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory rotation driver.
   * @param amount Pointer to the per-frame rotation in turns; must not be null.
   * @details The driver must keep |amount * 2pi| under
   * PALETTE_PHASE_ARG_LIMIT.
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
    assert(fabsf(*amount) * (2.0f * PI_F) < PALETTE_PHASE_ARG_LIMIT);
    if (!primed || *amount != cached_amount) {
      cached_amount = *amount;
      float ca, sa;
      turn_to_unit_cos_sin(cached_amount, ca, sa);
      hue_rotate_lms_matrix(ca, sa, matrix);
      primed = true;
    }
    LinRGB rgb = pixel_to_linrgb(c.color);
    LMS lms = linear_rgb_to_lms(rgb.r, rgb.g, rgb.b);
    float cl, cm, cs;
    fast_cbrt3(lms.l, lms.m, lms.s, cl, cm, cs);
    lms_cbrt_transform_rgb(matrix, cl, cm, cs, rgb.r, rgb.g, rgb.b);
    c.color = linrgb_to_pixel(rgb);
    return c;
  }
};

/**
 * @brief Rotates hue by an amount that varies along the palette domain, so
 * different parts of the gradient drift in opposite directions (iridescence).
 * @details Builds a rotation per sample; suited to bake-time sampling
 * (BakedPaletteStorage::rebake) rather than tight per-pixel loops.
 */
struct HueWobbleShade {
  const float *phase;
  float frequency;
  float depth;

  /**
   * @brief Constructs with a mandatory phase driver, frequency, and depth.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param freq Wobble frequency over the domain; defaults to 1.0.
   * @param depth Peak hue rotation in turns; |depth| * 2pi must stay under
   *        PALETTE_PHASE_ARG_LIMIT. Defaults to 0.1.
   * @details The caller must keep |t * freq * 2pi + phase| under
   * PALETTE_PHASE_ARG_LIMIT.
   */
  HueWobbleShade(const float *phase, float freq = 1.0f, float depth = 0.1f)
      : phase(phase), frequency(freq), depth(depth) {
    HS_CHECK(phase, "HueWobbleShade: phase driver must not be null");
    HS_CHECK(fabsf(depth) * (2.0f * PI_F) < PALETTE_PHASE_ARG_LIMIT,
             "HueWobbleShade: depth must stay inside the fast-trig argument "
             "range");
  }

  /**
   * @brief Rotates the sample's hue by depth * sin(t*frequency*2pi + phase)
   * turns, preserving alpha.
   * @param c Sample color to reshape.
   * @param t Coordinate driving the wobble.
   * @return The hue-rotated sample.
   */
  Color4 shade(Color4 c, float t) const {
    const float arg = t * frequency * PI_F * 2.0f + *phase;
    assert(fabsf(arg) < PALETTE_PHASE_ARG_LIMIT);
    return hue_rotate(c, depth * fast_sinf(arg));
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
   * @brief Per-instance memo of fast_sinf(*phase).
   * @details *phase is frame-constant, so the sine is recomputed once per
   * frame, not per sample. depth is applied outside the memo, so a live depth
   * change lands on a frozen phase. mutable so const shade() can update the
   * memo.
   */
  mutable float cached_phase = 0.0f;
  mutable float cached_sin = 0.0f; /**< Memoized sine at cached_phase. */
  mutable bool primed = false;     /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory phase driver and pulse depth.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param depth Pulse depth in [0, 1]: chroma swings over [1-depth, 1+depth].
   *        Defaults to 0.5.
   * @details The driver must keep |phase| under PALETTE_PHASE_ARG_LIMIT.
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
    assert(fabsf(*phase) < PALETTE_PHASE_ARG_LIMIT);
    if (!primed || *phase != cached_phase) {
      cached_phase = *phase;
      cached_sin = fast_sinf(cached_phase);
      primed = true;
    }
    const float scale = 1.0f + depth * cached_sin;
    LinRGB rgb = pixel_to_linrgb(c.color);
    OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
    lab.a *= scale;
    lab.b *= scale;
    oklab_to_linear_rgb_gamut(lab, rgb.r, rgb.g, rgb.b);
    c.color = linrgb_to_pixel(rgb);
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
   * @param weight Overlay strength; must be non-negative. Defaults to 0.25.
   * @details The caller must keep |t * freq * 2pi + phase| under
   * PALETTE_PHASE_ARG_LIMIT.
   */
  IridescentShade(const float *phase, float freq = 3.0f, float weight = 0.25f)
      : phase(phase), frequency(freq), weight(weight) {
    HS_CHECK(phase, "IridescentShade: phase driver must not be null");
    HS_CHECK(weight >= 0.0f, "IridescentShade: weight must be non-negative");
  }

  /**
   * @brief Adds the weighted cosine sheen to the sample.
   * @param c Sample color to reshape.
   * @param t Coordinate locating the sample along the sheen.
   * @return The sample plus the overlay (per-channel saturating), alpha
   *         untouched.
   */
  Color4 shade(Color4 c, float t) const {
    const float arg = t * frequency * PI_F * 2.0f + *phase;
    assert(fabsf(arg) < PALETTE_PHASE_ARG_LIMIT);
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
      return Color4(
          black.lerp16(c.color, frac_to_q16(quintic_kernel((1.0f - t) / edge))),
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
      c.alpha *= quintic_kernel((1.0f - t) / edge);
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
 * @brief Which coordinate a StaticPalette hands to its color-modifier chain.
 */
enum class ShadeCoord : uint8_t {
  /** The coordinate the source was sampled at when Wrap, the raw input when
   *  not — the default. */
  MATCH_WRAP,
  /** The coordinate the source was sampled at, whatever Wrap is. */
  LOOKUP,
  /** The raw pre-modifier input, whatever Wrap is. */
  RAW_INPUT,
};

/**
 * @brief A compile-time composition of a Source palette, a coordinate-modifier
 * chain, and a color-modifier chain.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CoordList Coords<> type-list of coordinate modifiers.
 * @tparam ColorList Colors<> type-list of color modifiers.
 * @tparam Wrap Folds the coordinate into [0,1) before the source lookup.
 * @tparam Shade Which coordinate the color-modifier chain receives.
 * @details Default construct, then bind() (ArenaVector idiom); both chains are
 * inlined by fold expression. get() applies the coord mods to t in order,
 * samples the source (wrapping the coordinate unless Wrap is false), then
 * applies the color mods with the coordinate Shade selects. Wrap=false suits
 * inset/falloff pipelines that must reach the source's exact endpoints
 * (wrap_t(1)==0 would otherwise fold the top edge). Wrap is checked at compile
 * time: `requires_wrap` on any unbounded modifier rejects Wrap=false unless a
 * later modifier declares `rebounds_input`, and `bounded_output` on the final
 * coord modifier rejects Wrap=true — so a chain can force Wrap, which is why
 * the shading coordinate is a separate knob rather than a second meaning of it.
 */
template <typename Source, typename CoordList = Coords<>,
          typename ColorList = Colors<>, bool Wrap = true,
          ShadeCoord Shade = ShadeCoord::MATCH_WRAP>
class StaticPalette;

/**
 * @brief Partial specialization splitting the two modifier type-lists.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CMods Coordinate modifier types.
 * @tparam XMods Color modifier types.
 * @tparam Wrap Folds the coordinate before the source lookup.
 * @tparam Shade Which coordinate the color chain receives.
 */
template <typename Source, typename... CMods, typename... XMods, bool Wrap,
          ShadeCoord Shade>
class StaticPalette<Source, Coords<CMods...>, Colors<XMods...>, Wrap, Shade> {
  static_assert((CoordMod<CMods> && ...), "Coords<> entries must be CoordMods");
  static_assert((ColorMod<XMods> && ...), "Colors<> entries must be ColorMods");
  static_assert(Wrap || !coord_chain_leaves_unit<CMods...>(),
                "Wrap=false composed with a coordinate modifier that leaves "
                "[0,1] (requires_wrap, e.g. CycleModifier/ScaleModifier): the "
                "source would be sampled out of range and the palette would "
                "freeze at its endpoint. Use Wrap=true, or follow it with a "
                "modifier that re-bounds arbitrary input (rebounds_input, e.g. "
                "WrapModifier/FoldModifier/InsetModifier).");
  static_assert(!Wrap || !coord_chain_bounded_tail<CMods...>(),
                "Wrap=true with a bounded final coordinate modifier "
                "(bounded_output, e.g. MirrorModifier/InsetModifier): wrap_t "
                "folds its 1.0 output to 0.0, destroying the top endpoint. "
                "Use Wrap=false; if an earlier modifier leaves [0,1], insert a "
                "WrapModifier ahead of the bounded tail.");

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
   * sample with the coordinate Shade selects.
   */
  Color4 get(float t) const {
    assert(source != nullptr && "StaticPalette used before bind()!");

    float ft = t;
    std::apply([&](const auto *...m) { ((ft = m->modify(ft)), ...); }, coords);

    float u = ft;
    if constexpr (Wrap)
      u = wrap_t(ft);
    Color4 c = source->get(u);

    float shade_coordinate;
    if constexpr (Shade == ShadeCoord::LOOKUP)
      shade_coordinate = u;
    else if constexpr (Shade == ShadeCoord::RAW_INPUT)
      shade_coordinate = t;
    else
      shade_coordinate = Wrap ? u : t;
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
 * `const Palette*` world (preset tables, BakedPaletteStorage::bake). The virtual call
 * is paid only at bake time (cold), never on the per-pixel path.
 *
 * SP must not wrap its coordinate. A wrapping source folds t = 1 back to 0 and
 * collapses a bake's last entry onto its first; BakedPalette rejects one at
 * compile time, but the erasure to `const Palette&` this class performs is what
 * would hide it from that check.
 */
template <typename SP> class PaletteFacade : public Palette {
  static_assert(!palette_wraps_coordinate<SP>(),
                "PaletteFacade cannot erase a wrapping palette: bake and "
                "sample it through Wrap=false");

public:
  /**
   * @brief Default-constructs an unbound facade (bind() before use).
   */
  PaletteFacade() = default;
  /**
   * @brief Constructs a facade bound to a composition.
   * @param sp Composition to forward get() to; must not be null.
   */
  explicit PaletteFacade(const SP *sp) : composition(sp) {
    HS_CHECK(sp != nullptr, "PaletteFacade constructed with null composition");
  }
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
