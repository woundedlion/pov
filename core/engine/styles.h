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
 * @brief Identity space transform: no spatial distortion.
 * @param v Sample direction on the unit sphere.
 * @return The unchanged sample direction.
 */
inline Vector identity_warp(const Vector &v, const Style &) { return v; }

/**
 * @brief Hue-rotating fade (default color transform).
 * @param p Source pixel color.
 * @param fade Per-frame scalar fade multiplier in [0, 1].
 * @param s Style supplying the precomputed hue rotation.
 * @return Faded and hue-rotated pixel.
 */
inline Pixel hue_fade(const Pixel &p, float fade, const Style &s);

/**
 * @brief Plain scalar fade color transform: no hue shift.
 * @param p Source pixel color.
 * @param fade Per-frame scalar fade multiplier in [0, 1].
 * @return Pixel scaled by fade.
 */
inline Pixel plain_fade(const Pixel &p, float fade, const Style &) {
  return p * fade;
}

// --- Style struct -------------------------------------------------------------

/**
 * @brief Named feedback preset: spatial/color transforms plus scalar params.
 * @details POD-copyable so it can live in Presets<> and be lerped. The lerpable
 * scalar fields blend continuously while function pointers and discrete tuning
 * snap at the lerp midpoint. Bound state (noise) and per-frame caches are not
 * preset data.
 */
struct Style {
  // --- Lerpable scalar params ---
  float fade      = 0.95f;
  float hue_shift = 0.0f;
  float amplitude = 0.5f;
  float frequency = 0.125f;
  float speed     = 1.0f;
  float scale     = 4.0f;

  // --- Function pointers (snap during lerp) ---
  SpaceFn space_fn = &noise_warp;
  ColorFn color_fn = &hue_fade;

  // --- Filter tuning (snap during lerp) ---
  /**
   * Coarse-grid downsample factor for the warp field. Higher = cheaper
   * (~DS^2 fewer space_fn / atan2 / acos calls) but loses detail in
   * high-frequency warps. Pick smaller values (2) for sharp noise,
   * larger (8) for slow swirling warps. Scratch arena must hold
   * (W/DS) * (H/DS) * 4 bytes — at 288x144, DS=4 ≈ 10KB, DS=2 ≈ 41KB.
   */
  int downsample = 4;

  // --- Bound state (set by effect at init, NOT part of presets) ---
  NoiseParams *noise = nullptr;

  // --- Per-frame derived cache (NOT a preset; refreshed by sync_hue) ---
  // cos/sin of hue_shift's turn angle, read by hue_fade. Defaults to the
  // identity rotation (angle 0) until the first sync_hue().
  float hue_ca = 1.0f;
  float hue_sa = 0.0f;

  /**
   * @brief Interpolate this Style between two endpoints.
   * @param a Style at t = 0.
   * @param b Style at t = 1.
   * @param t Interpolation fraction in [0, 1].
   * @details Scalar params blend continuously while function pointers and
   * discrete tuning snap at the midpoint (t = 0.5). The bound noise pointer is
   * left untouched: it is effect-owned state set at init, not preset data, so
   * pulling it from a preset would overwrite it with that preset's nullptr and
   * silently degrade noise_warp to identity.
   */
  void lerp(const Style &a, const Style &b, float t) {
    fade      = hs::lerp(a.fade,      b.fade,      t);
    hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
    amplitude = hs::lerp(a.amplitude, b.amplitude, t);
    frequency = hs::lerp(a.frequency, b.frequency, t);
    speed     = hs::lerp(a.speed,     b.speed,     t);
    scale     = hs::lerp(a.scale,     b.scale,     t);
    space_fn   = t < 0.5f ? a.space_fn   : b.space_fn;
    color_fn   = t < 0.5f ? a.color_fn   : b.color_fn;
    downsample = t < 0.5f ? a.downsample : b.downsample;
  }

  /**
   * @brief Precompute hue_shift's rotation (cos/sin) into hue_ca/hue_sa.
   * @details Caches the frame-constant rotation for the per-pixel hue_fade hot
   * path. Call once per frame before the feedback sampling loop.
   */
  void sync_hue() {
    float angle = hue_shift * (2.0f * PI_F);
    float ca = fast_cosf(angle);
    float sa = fast_sinf(angle);
    // fast trig is non-orthonormal; renormalize so hue_rotate preserves chroma
    // (else the scaling compounds per frame under feedback).
    float inv = 1.0f / sqrtf(ca * ca + sa * sa);
    hue_ca = ca * inv;
    hue_sa = sa * inv;
  }

  /**
   * @brief Push this Style's scalar params into the bound NoiseParams.
   * @details No-op if no NoiseParams is bound. Copies amplitude, frequency,
   * speed, and scale, then calls NoiseParams::sync().
   *
   * `const` while sibling sync_hue() is not: the constness reflects what each
   * writes, not a sync/no-sync distinction. sync_noise() leaves this Style's own
   * members untouched and writes only *through* the bound `noise` pointer (which
   * the pointed-to object's own constness, not ours, governs), so it is const.
   * sync_hue() mutates this Style's own hue_ca/hue_sa cache, so it cannot be.
   */
  void sync_noise() const {
    if (!noise) return;
    noise->amplitude = amplitude;
    noise->frequency = frequency;
    noise->speed     = speed;
    noise->scale     = scale;
    noise->sync();
  }

  // --- Named presets ---
  // Params: fade, hue_shift, amplitude, frequency, speed, scale, space_fn, color_fn

  /**
   * @brief Static fine-grained turbulence: high amplitude over a tight scale with
   * no temporal drift, giving a frozen, twisted distortion.
   * @return The SlowTwist preset Style.
   */
  static constexpr Style SlowTwist() {
    return {0.8158f, 0.041f, 6.36f, 0.21f, 0.0f, 2.1459f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Dense fine-grain turbulence with strong hue shift; tight scale, slow drift.
   * @return The Churn preset Style.
   */
  static constexpr Style Churn() {
    return {0.82f, 0.035f, 3.15f, 1.0f, 0.02f, 1.0f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Gentle drifting haze with slow noise; classic smoke look.
   * @return The Smoke preset Style.
   */
  static constexpr Style Smoke() {
    return {0.9f, 0.01f, 0.51f, 0.42f, 0.46f, 23.0f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Static frozen distortion — no temporal movement.
   * @return The Frozen preset Style.
   */
  static constexpr Style Frozen() {
    return {0.58f, 0.03f, 2.73f, 0.07f, 0.0f, 26.0f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Extreme static warping with fast decay; shattering glass look.
   * @return The Shatter preset Style.
   */
  static constexpr Style Shatter() {
    return {0.58f, 0.03f, 8.21f, 0.01f, 0.0f, 46.0f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Flowing medium-strength distortion; gentle liquid drift.
   * @return The Drift preset Style.
   */
  static constexpr Style Drift() {
    return {0.68f, 0.03f, 4.98f, 0.07f, 0.2f, 5.0f, &noise_warp, &hue_fade};
  }

  /**
   * @brief Image melts and drips downward off the sphere.
   * @return The Melting preset Style.
   */
  static constexpr Style Melting() {
    return {0.8158f, 0.0766f, 6.36f, 0.014f, 1.005f, 42.365f, &melt_warp, &hue_fade};
  }

  /**
   * @brief Fast downward swirl with strong distortion, no hue shift.
   * @return The Swirling preset Style.
   */
  static constexpr Style Swirling() {
    return {0.8158f, 0.0f, 6.36f, 0.014f, 1.465f, 42.365f, &melt_warp, &plain_fade};
  }
};

// Compile-time anchor pinning SlowTwist's resolved fields by name: a reorder of
// Style's same-typed scalar members would silently reassign the positional
// preset literals, which this catches at compile time.
static_assert(Style::SlowTwist().fade == 0.8158f &&
                  Style::SlowTwist().hue_shift == 0.041f &&
                  Style::SlowTwist().amplitude == 6.36f &&
                  Style::SlowTwist().frequency == 0.21f &&
                  Style::SlowTwist().speed == 0.0f &&
                  Style::SlowTwist().scale == 2.1459f &&
                  Style::SlowTwist().space_fn == &noise_warp &&
                  Style::SlowTwist().color_fn == &hue_fade,
              "Style preset field order drifted from the positional brace-init "
              "in the *() presets; update the initializers or this anchor.");

// --- Deferred inline definitions (Style is now complete) ----------------------

inline Vector noise_warp(const Vector &v, const Style &s) {
  if (!s.noise) return v;
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

// Uses the rotation precomputed once per frame by Style::sync_hue.
inline Pixel hue_fade(const Pixel &p, float fade, const Style &s) {
  return hue_rotate(Color4(p * fade, 1.0f), s.hue_ca, s.hue_sa).color;
}

} // namespace Feedback
