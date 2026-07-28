/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "engine/transformers.h"
#include "color/color.h"

/**
 * @brief Feedback presets that bundle spatial/color transforms with scalar
 * parameters, plus the transform functions that consume them.
 * @details Style is POD-copyable — safe to store in Presets<> and lerp.
 * Typical usage:
 *   Feedback::Style style = Feedback::Style::Smoke();
 *   style.noise = &amp;my_noise_params;   // bind effect-owned state at init
 *   style.sync_noise();               // push scalars → NoiseParams each frame
 */
namespace Feedback {

struct Style;

/**
 * @brief Pointer type for a space (spatial warp) transform.
 * @details Receives the sample vector and the full Style for state access,
 * returning the warped sample direction.
 */
using SpaceFn = Vector (*)(const Vector &, const Style &);
/**
 * @brief Pointer type for a color (fade) transform.
 * @details Receives the pixel, the per-frame fade value, and the Style,
 * returning the faded/shifted pixel.
 */
using ColorFn = Pixel (*)(const Pixel &, float fade, const Style &);

// --- Transform implementations ------------------------------------------------

/**
 * @brief Noise-based spatial warping (default space transform).
 * @param v Sample direction on the unit sphere.
 * @param s Style supplying the bound NoiseParams.
 * @return Warped sample direction.
 */
inline Vector noise_warp(const Vector &v, const Style &s);

/**
 * @brief Downward melt space transform: slerps toward the north pole so the
 * image drips south, with added noise wobble.
 * @param v Sample direction on the unit sphere.
 * @param s Style supplying speed (drip rate) and the bound NoiseParams.
 * @return Warped sample direction.
 */
inline Vector melt_warp(const Vector &v, const Style &s);

/**
 * @brief Hue-rotating fade (default color transform).
 * @param p Source pixel color.
 * @param fade Per-frame scalar fade multiplier in [0, 1].
 * @param s Style supplying the precomputed hue rotation.
 * @return Faded and hue-rotated pixel.
 */
inline Pixel hue_fade(const Pixel &p, float fade, const Style &s);

// --- Style struct -------------------------------------------------------------

/**
 * @brief Named feedback preset: spatial/color transforms plus scalar params.
 * @details POD-copyable for Presets<> and lerp. Non-preset state (bound noise
 * pointer, per-frame hue cache) survives lerp() but a full-struct copy or
 * assignment (Presets::apply, `style = Style::Smoke()`) resets noise to nullptr
 * (noise_warp degrades to identity) and the hue cache to identity. Re-bind noise
 * and call sync_hue() after any such copy.
 */
struct Style {
  // --- Lerpable scalar params ---
  float fade = 0.95f;
  /** Hue rotation per e-fold decrease in feedback brightness, in turns. */
  float hue_shift = 0.0f;
  float amplitude = 0.5f;
  float frequency = 0.125f;
  float speed = 1.0f;
  float scale = 4.0f;

  // --- Function pointers (snap during lerp) ---
  SpaceFn space_fn = &noise_warp;
  ColorFn color_fn = &hue_fade;

  // --- Filter tuning (snap during lerp) ---
  /**
   * Coarse-grid downsample factor for the warp field. Higher = cheaper
   * (~DS^2 fewer space_fn / atan2 / acos calls), lower = more detail. One
   * flush holds the coarse x/y offset grid, the spherical control samples it
   * expands from, and one W-pixel row in the scratch arena at once; pole infill
   * puts the ring and sample counts above a flat (W/DS) x (H/DS) grid. At
   * 288x144, DS=4 ≈ 21KB, DS=2 ≈ 71KB.
   */
  int downsample = 4;

  // --- Bound state (set by effect at init, NOT part of presets) ---
  NoiseParams *noise = nullptr;

  // --- Per-frame derived cache (NOT a preset; refreshed by sync_hue) ---
  // cos/sin of the fade-scaled per-frame hue angle plus its cbrt-LMS rotation
  // matrix, read by hue_fade. Defaults to identity until the first sync_hue().
  float hue_ca = 1.0f;
  float hue_sa = 0.0f;
  float hue_k[9] = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};

  /**
   * @brief Interpolate this Style between two endpoints.
   * @param a Style at t = 0.
   * @param b Style at t = 1.
   * @param t Interpolation fraction in [0, 1].
   * @details Scalar params blend continuously; function pointers and discrete
   * tuning snap at t = 0.5. The bound noise pointer is left untouched (effect-
   * owned state, not preset data; pulling it from a preset would null it and
   * degrade noise_warp to identity).
   */
  void lerp(const Style &a, const Style &b, float t) {
    fade = hs::lerp(a.fade, b.fade, t);
    hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
    amplitude = hs::lerp(a.amplitude, b.amplitude, t);
    frequency = hs::lerp(a.frequency, b.frequency, t);
    speed = hs::lerp(a.speed, b.speed, t);
    scale = hs::lerp(a.scale, b.scale, t);
    space_fn = t < 0.5f ? a.space_fn : b.space_fn;
    color_fn = t < 0.5f ? a.color_fn : b.color_fn;
    downsample = t < 0.5f ? a.downsample : b.downsample;
  }

  /**
   * @brief Precompute the fade-scaled hue rotation into hue_ca/hue_sa and hue_k.
   * @details Each frame applies `hue_shift * -log(fade)` turns, making
   * accumulated hue depend only on remaining brightness. A zero fade uses the
   * identity rotation because no feedback remains visible.
   */
  void sync_hue() {
    float frame_shift = fade == 0.0f ? 0.0f : hue_shift * -logf(fade);
    turn_to_unit_cos_sin(frame_shift, hue_ca, hue_sa);
    hue_rotate_lms_matrix(hue_ca, hue_sa, hue_k);
  }

  /**
   * @brief Push this Style's scalar params into the bound NoiseParams.
   * @details No-op if no NoiseParams is bound. Copies amplitude, frequency,
   * speed, and scale, then calls NoiseParams::sync(). The feedback filter traps
   * on flush when a bound NoiseParams has drifted from these scalars, so any
   * write to them (slider, preset switch, lerp) must be followed by a call.
   */
  void sync_noise() const {
    if (!noise)
      return;
    noise->amplitude = amplitude;
    noise->frequency = frequency;
    noise->speed = speed;
    noise->scale = scale;
    noise->sync();
  }

  // --- Named presets ---
  // Params: fade, hue_shift, amplitude, frequency, speed, scale, space_fn, color_fn

  /**
   * @brief Branching, fast-moving distortion with pronounced hue rotation.
   * @return The ArcingLightning preset Style.
   */
  static constexpr Style ArcingLightning() {
    return {0.5f, 0.14426951f, 3.27f,       0.09f,
            1.5f, 50.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Broad, slowly evolving turbulence with gentle color drift.
   * @return The SlowFire preset Style.
   */
  static constexpr Style SlowFire() {
    return {0.8732f, 0.12316483f, 1.56f,       0.5297f,
            0.1f,    50.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Broad, quickly evolving turbulence with gentle color drift.
   * @return The EnergeticFire preset Style.
   */
  static constexpr Style EnergeticFire() {
    return {0.8732f, 0.12316483f, 1.56f,       0.22087f,
            0.9f,    50.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Gentle drifting haze with slow noise; classic smoke look.
   * @return The Smoke preset Style.
   */
  static constexpr Style Smoke() {
    return {0.9f,  0.09491219f, 0.51f,       0.42f,
            0.46f, 23.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Fine, slowly drifting turbulence with gentle color rotation.
   * @return The SlowDust preset Style.
   */
  static constexpr Style SlowDust() {
    return {0.83952f, 0.09546948f, 1.56f,       0.07237f,
            0.6f,     50.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Fine, rapidly moving distortion with pronounced color trails.
   * @return The WavyTrails preset Style.
   */
  static constexpr Style WavyTrails() {
    return {0.7257f, 0.22518973f, 1.95f,       0.01f,
            5.0f,    50.0f,       &noise_warp, &hue_fade};
  }

  /**
   * @brief Strong downward melt with slow drift and pronounced hue rotation.
   * @return The MeltingHi preset Style.
   */
  static constexpr Style MeltingHi() {
    return {0.59004f, 0.18955015f, 4.38f,      0.06346f,
            0.2f,     22.3554f,    &melt_warp, &hue_fade};
  }

  /**
   * @brief Gentle downward melt with slow drift and pronounced hue rotation.
   * @return The MeltingLo preset Style.
   */
  static constexpr Style MeltingLo() {
    return {0.59004f, 0.18955015f, 1.95f,      0.06346f,
            0.2f,     22.3554f,    &melt_warp, &hue_fade};
  }

  /**
   * @brief Drifting toxic haze: medium turbulence with slow drift and strong
   * per-frame hue cycling.
   * @return The Miasma preset Style.
   */
  static constexpr Style Miasma() {
    return {0.80586f, 0.234f,     2.61f,       0.05059f,
            0.725f,   26.297501f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Static high-amplitude twist over a medium scale; a loose swirling
   * tunnel with no temporal drift.
   * @return The LooseWormhole preset Style.
   */
  static constexpr Style LooseWormhole() {
    return {0.7257f, 0.22519f, 11.25f,      0.01f,
            0.0f,    10.1798f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Static high-amplitude twist over a tight scale; a tight swirling
   * tunnel with no temporal drift.
   * @return The TightWormhole preset Style.
   */
  static constexpr Style TightWormhole() {
    return {0.7257f, 0.22519f, 6.42f,       0.01f,
            0.0f,    7.8844f,  &noise_warp, &hue_fade};
  }

  /**
   * @brief Static twist over a broad scale; a wide wormhole with wandering
   * arms and no temporal drift.
   * @return The WigglingWormhole preset Style.
   */
  static constexpr Style WigglingWormhole() {
    return {0.7257f, 0.22519f, 7.11f,       0.01f,
            0.0f,    29.1917f, &noise_warp, &hue_fade};
  }
};

// Compile-time anchor pinning Miasma's resolved fields by name: a reorder of
// Style's same-typed scalar members would silently reassign the positional
// preset literals, which this catches at compile time.
static_assert(Style::Miasma().fade == 0.80586f &&
                  Style::Miasma().hue_shift == 0.234f &&
                  Style::Miasma().amplitude == 2.61f &&
                  Style::Miasma().frequency == 0.05059f &&
                  Style::Miasma().speed == 0.725f &&
                  Style::Miasma().scale == 26.297501f &&
                  Style::Miasma().space_fn == &noise_warp &&
                  Style::Miasma().color_fn == &hue_fade,
              "Style preset field order drifted from the positional brace-init "
              "in the *() presets; update the initializers or this anchor.");

// --- Deferred inline definitions (Style is now complete) ----------------------

inline Vector noise_warp(const Vector &v, const Style &s) {
  if (!s.noise)
    return v;
  return noise_transform(v, *s.noise);
}

inline Vector melt_warp(const Vector &v, const Style &s) {
  // Shift sample toward north pole → image appears to drip south. speed controls
  // drip rate; amplitude controls noise wobble.
  static constexpr Vector NORTH = {0.0f, 1.0f, 0.0f};
  // Slerp fraction toward the pole per frame at speed=1 (preset speeds scale it).
  static constexpr float MELT_STEP_PER_FRAME = 0.04f;
  // Amplitude floor below which the noise wobble is skipped.
  static constexpr float MELT_NOISE_AMP_FLOOR = 0.001f;
  float drip = s.speed * MELT_STEP_PER_FRAME;
  Vector drifted = slerp(v, NORTH, drip);

  if (s.noise && s.amplitude > MELT_NOISE_AMP_FLOOR) {
    return noise_transform(drifted, *s.noise);
  }
  return drifted;
}

/**
 * @brief Applies a cbrt-LMS rotation to linear-RGB channels.
 * @param k Rotation matrix already scaled by cbrt(fade/65535) (see hue_fade and
 *          the Feedback::flush fast path, which both prescale their fade here).
 * @param r Linear red channel (u16 magnitude, as float).
 * @param g Linear green channel.
 * @param b Linear blue channel.
 * @return The rotated, faded pixel.
 * @details Single source for the built-in hue-fade color math so hue_fade() and
 *          the flush() fast path that bypasses it cannot drift.
 */
HS_O3_FN inline Pixel hue_fade_apply(const float k[9], float r, float g,
                                     float b) {
  LMS lms = linear_rgb_to_lms(r, g, b);
  float cl, cm, cs;
  fast_cbrt3(lms.l, lms.m, lms.s, cl, cm, cs);
  float rr, gg, bb;
  lms_cbrt_transform_rgb(k, cl, cm, cs, rr, gg, bb);
  return Pixel(float_to_pixel16(rr), float_to_pixel16(gg),
               float_to_pixel16(bb));
}

/**
 * @brief Two-pixel hue_fade_apply sharing one code path.
 * @param k Rotation matrix already scaled by cbrt(fade/65535).
 * @param r0 Linear red of the first pixel (u16 magnitude, as float).
 * @param g0 Linear green of the first pixel.
 * @param b0 Linear blue of the first pixel.
 * @param r1 Linear red of the second pixel.
 * @param g1 Linear green of the second pixel.
 * @param b1 Linear blue of the second pixel.
 * @param p0 Out: the rotated, faded first pixel.
 * @param p1 Out: the rotated, faded second pixel.
 * @details Both pixels' cube roots go through one fast_cbrt6 call, so the six
 * seed/Halley chains schedule against each other and the pair costs one divide.
 * The shared reciprocal re-associates the arithmetic, moving results off the
 * scalar hue_fade_apply path by ~4e-7 relative — ~50x below the cube root's own
 * ~2.3e-5 error, and 0.002% of u16 channels by one LSB.
 */
HS_O3_FN inline void hue_fade_apply2(const float k[9], float r0, float g0,
                                     float b0, float r1, float g1, float b1,
                                     Pixel &p0, Pixel &p1) {
  LMS lms0 = linear_rgb_to_lms(r0, g0, b0);
  LMS lms1 = linear_rgb_to_lms(r1, g1, b1);
  const float lms6[6] = {lms0.l, lms0.m, lms0.s, lms1.l, lms1.m, lms1.s};
  float c6[6];
  fast_cbrt6(lms6, c6);
  float rr0, gg0, bb0, rr1, gg1, bb1;
  lms_cbrt_transform_rgb2(k, c6[0], c6[1], c6[2], c6[3], c6[4], c6[5], rr0, gg0,
                          bb0, rr1, gg1, bb1);
  p0 = Pixel(float_to_pixel16(rr0), float_to_pixel16(gg0),
             float_to_pixel16(bb0));
  p1 = Pixel(float_to_pixel16(rr1), float_to_pixel16(gg1),
             float_to_pixel16(bb1));
}

// The fade and the u16 normalization fold into the cbrt-LMS domain:
// cbrt(fade/65535 * LMS) = cbrt(fade/65535) * cbrt(LMS). Uses the rotation
// matrix precomputed once per frame by Style::sync_hue.
inline Pixel hue_fade(const Pixel &p, float fade, const Style &s) {
  const float sc = fast_cbrt(fade * (1.0f / 65535.0f));
  float k[9];
  for (int i = 0; i < 9; ++i)
    k[i] = s.hue_k[i] * sc;
  return hue_fade_apply(k, p.r, p.g, p.b);
}

} // namespace Feedback
