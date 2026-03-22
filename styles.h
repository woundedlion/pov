/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "transformers.h"
#include "color.h"
#include "filter.h"

/**
 * @brief Named feedback presets that bundle spatial/color transforms with
 * scalar parameters. POD-copyable — safe to store in Presets<> and lerp.
 *
 * Usage:
 *   Feedback::Style style = Feedback::Style::Smoke();
 *   style.noise = &my_noise_params;   // bind effect-owned state at init
 *   style.sync_noise();               // push scalars → NoiseParams each frame
 */
namespace Feedback {

struct Style;

/// Space transform: receives the vector and the full Style for state access.
using SpaceFn = Vector (*)(const Vector &, const Style &);
/// Color transform: receives the pixel, fade value, and Style.
using ColorFn = Pixel (*)(const Pixel &, float fade, const Style &);

// --- Transform implementations ------------------------------------------------

/// Noise-based spatial warping (default).
inline Vector noise_warp(const Vector &v, const Style &s);

/// Identity — no spatial distortion.
inline Vector identity_warp(const Vector &v, const Style &) { return v; }

/// Hue-rotating fade (default).
inline Pixel hue_fade(const Pixel &p, float fade, const Style &s);

/// Plain scalar fade — no color shift.
inline Pixel plain_fade(const Pixel &p, float fade, const Style &) {
  return p * fade;
}

// --- Style struct -------------------------------------------------------------

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

  // --- Bound state (set by effect at init, NOT part of presets) ---
  NoiseParams *noise = nullptr;

  void lerp(const Style &a, const Style &b, float t) {
    fade      = ::lerp(a.fade,      b.fade,      t);
    hue_shift = ::lerp(a.hue_shift, b.hue_shift, t);
    amplitude = ::lerp(a.amplitude, b.amplitude, t);
    frequency = ::lerp(a.frequency, b.frequency, t);
    speed     = ::lerp(a.speed,     b.speed,     t);
    scale     = ::lerp(a.scale,     b.scale,     t);
    // Function pointers snap at midpoint
    space_fn = t < 0.5f ? a.space_fn : b.space_fn;
    color_fn = t < 0.5f ? a.color_fn : b.color_fn;
    // Bound pointer preserved (not lerped)
    noise = a.noise;
  }

  /// Push scalar values into the bound NoiseParams.
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

  /// Dense fine-grain turbulence with strong hue shift. Tight scale, slow drift.
  static constexpr Style Churn() {
    return {0.82f, 0.035f, 3.15f, 1.0f, 0.02f, 1.0f, &noise_warp, &hue_fade};
  }

  /// Gentle drifting haze with slow noise. Classic smoke look.
  static constexpr Style Smoke() {
    return {0.9f, 0.01f, 0.51f, 0.42f, 0.46f, 23.0f, &noise_warp, &hue_fade};
  }

  /// Static frozen distortion — no temporal movement.
  static constexpr Style Frozen() {
    return {0.58f, 0.03f, 2.73f, 0.07f, 0.0f, 26.0f, &noise_warp, &hue_fade};
  }

  /// Extreme static warping with fast decay. Shattering glass look.
  static constexpr Style Shatter() {
    return {0.58f, 0.03f, 8.21f, 0.01f, 0.0f, 46.0f, &noise_warp, &hue_fade};
  }

  /// Flowing medium-strength distortion. Gentle liquid drift.
  static constexpr Style Drift() {
    return {0.68f, 0.03f, 4.98f, 0.07f, 0.2f, 5.0f, &noise_warp, &hue_fade};
  }
};

// --- Deferred inline definitions (Style is now complete) ----------------------

inline Vector noise_warp(const Vector &v, const Style &s) {
  if (!s.noise) return v;
  return noise_transform(v, *s.noise);
}

inline Pixel hue_fade(const Pixel &p, float fade, const Style &s) {
  return hue_rotate(Color4(p * fade, 1.0f), s.hue_shift).color;
}

// --- Style-aware Feedback Filter ----------------------------------------------

/// Internal adapters — bridge Style's function pointers to FunctionRef callables.
struct SpaceAdapter_ {
  const Style *s;
  Vector operator()(const Vector &v) const { return s->space_fn(v, *s); }
};

struct ColorAdapter_ {
  const Style *s;
  Pixel operator()(const Pixel &p, float fade) const {
    return s->color_fn(p, fade, *s);
  }
};

/**
 * @brief Style-based Feedback filter. Takes a Style& directly — no template
 * params for transform types, no adapter boilerplate in effects.
 *
 * Usage:
 *   Feedback::Style style = Feedback::Style::Smoke();
 *   Pipeline<W, H, ..., Feedback::Filter<W, H>> filters(..., Feedback::Filter<W, H>(style));
 */
template <int W, int H>
class Filter : public Is2DWithHistory {
  using Inner = ::Filter::Pixel::Feedback<W, H, SpaceAdapter_, ColorAdapter_>;

public:
  explicit Filter(Style &style)
      : style_(&style),
        inner_(SpaceAdapter_{&style}, style.fade, ColorAdapter_{&style}) {}

  /// Pass-through: current-frame pixels go straight to next filter.
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFn2D pass) {
    inner_.plot(x, y, color, age, alpha, pass);
  }

  /// Blend distorted previous frame into current frame via the Style's transforms.
  void flush(Canvas &cv, const ScreenTrailFn &trail, float alpha,
             PassFn2D pass) {
    if (!enabled_) return;
    // Sync fade and function pointers from the (possibly lerping) Style
    inner_.set_fade(style_->fade);
    inner_.set_transform(SpaceAdapter_{style_});
    inner_.set_color_fn(ColorAdapter_{style_});
    inner_.flush(cv, trail, alpha, pass);
  }

  /// Enable/disable feedback (disabled = skip flush entirely).
  void set_enabled(bool e) { enabled_ = e; }

  /// Access the bound Style.
  Style &style() { return *style_; }
  const Style &style() const { return *style_; }

private:
  Style *style_;
  Inner inner_;
  bool enabled_ = true;
};

} // namespace Feedback
