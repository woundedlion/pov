/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include <cmath>

namespace SHMath {
inline float factorial(int n) {
  if (n <= 1)
    return 1.0f;
  float result = 1.0f;
  for (int i = 2; i <= n; i++)
    result *= i;
  return result;
}

inline float associatedLegendre(int l, int m, float x) {
  float pmm = 1.0f;
  if (m > 0) {
    float somx2 = sqrtf(std::max(0.0f, (1.0f - x) * (1.0f + x)));
    float fact = 1.0f;
    for (int i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0f;
    }
  }
  if (l == m)
    return pmm;

  float pmmp1 = x * (2.0f * m + 1.0f) * pmm;
  if (l == m + 1)
    return pmmp1;

  float pll = 0;
  for (int ll = m + 2; ll <= l; ll++) {
    pll = ((2.0f * ll - 1.0f) * x * pmmp1 - (ll + m - 1.0f) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

/// Precompute normalization factor for given (l, m) — constant per shape.
inline float normalization(int l, int m) {
  int absM = std::abs(m);
  float N = sqrtf(((2.0f * l + 1.0f) / (4.0f * PI_F)) *
                  (factorial(l - absM) / factorial(l + absM)));
  return (m != 0) ? sqrtf(2.0f) * N : N;
}

/// Evaluate SH with precomputed norm and cos_phi (avoids cosf(acos(x)) roundtrip).
inline float sphericalHarmonic(int l, int m, float theta, float cos_phi, float N) {
  int absM = std::abs(m);
  float P = associatedLegendre(l, absM, cos_phi);
  return N * P * ((m > 0) ? fast_cosf(m * theta) : (m < 0) ? fast_sinf(absM * theta) : 1.0f);
}
} // namespace SHMath

template <int W, int H> class SphericalHarmonics : public Effect {
public:
  struct HarmonicBlob {
    int l1, m1;
    int l2, m2;
    float blend;
    float amplitude;
    Quaternion orientation;
    static constexpr bool is_solid = true;
    float N1, N2; // Precomputed normalization factors

    HarmonicBlob(int l1, int m1, int l2, int m2, float blend, float amp,
                 Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend), amplitude(amp),
          orientation(q),
          N1(SHMath::normalization(l1, m1)),
          N2(SHMath::normalization(l2, m2)) {}

    template <int H_> SDF::Bounds get_vertical_bounds() const {
      return {0, H_ - 1}; // Full Scan fallback
    }

    template <int W_scan, int H_, typename OutputIt>
    bool get_horizontal_intervals(int y, OutputIt out) const {
      return false; // Full Scan fallback
    }

    SDF::DistanceResult distance(const Vector &p) const {
      SDF::DistanceResult res;
      distance<true>(p, res);
      return res;
    }

    template <bool ComputeUVs = true>
    void distance(const Vector &p, SDF::DistanceResult &res) const {
      Vector local = rotate(p, orientation.conjugate());

      float theta = fast_atan2(local.z, local.x);
      if (theta < 0)
        theta += 2 * PI_F;
      // local.y IS cos(phi) — skip fast_acos + cosf roundtrip
      float cos_phi = hs::clamp(local.y, -1.0f, 1.0f);

      // Evaluate primary shape (N1 precomputed per-frame)
      float val1 = SHMath::sphericalHarmonic(l1, m1, theta, cos_phi, N1);
      float val = val1;

      // Only pay the math cost for the second shape if we are actively blending
      if (blend > 0.001f) {
        float val2 = SHMath::sphericalHarmonic(l2, m2, theta, cos_phi, N2);
        val = val1 + (val2 - val1) * blend;
      }

      res = SDF::DistanceResult(-1.0f, 0.0f, val, 0.0f, 1.0f);
    }
  };

  FLASHMEM SphericalHarmonics() : Effect(W, H), filters() {}

  void init() override {
    registerParam("Amplitude", &params.amplitude, 0.1f, 10.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Bake procedural palette into fast LUT
    baked_palette.bake(persistent_arena, Palettes::richSunset);

    // Initial shape
    current_idx = 6;
    next_idx = 6;

    // Spin
    Vector axis = Vector(0.5f, 1.0f, 0.2f).normalized();
    timeline.add(0, Animation::Rotation<W>(orientation, axis, 2 * PI_F * 100,
                                           10000, ease_mid, true));

    // Start morphing immediately
    start_morph();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Decode flat indices into (l, m) pairs
    int l1 = (int)sqrtf((float)current_idx);
    int m1 = current_idx - l1 * l1 - l1;

    int l2 = (int)sqrtf((float)next_idx);
    int m2 = next_idx - l2 * l2 - l2;

    HarmonicBlob blob(l1, m1, l2, m2, morph_alpha, params.amplitude,
                      orientation.get());

    auto shader = [&](const Vector &p, Fragment &frag) {
      float val = frag.v1;
      float abs_val = std::abs(val);

      // Positive palette (baked LUT — no cosf/powf per pixel)
      Color4 pos =
          baked_palette.get(std::min(1.0f, abs_val * params.amplitude));
      // Negative palette (channel-swapped)
      Color4 neg =
          Color4(Pixel(pos.color.b, static_cast<uint8_t>(pos.color.g * 0.8f),
                       pos.color.r),
                 pos.alpha);

      // Quintic-smoothed crossfade across the zero-crossing boundary
      constexpr float transition = 0.03f;
      float blend_t =
          quintic_kernel(hs::clamp(val / transition * 0.5f + 0.5f, 0.0f, 1.0f));
      Color4 base = pos.lerp(neg, 1.0f - blend_t);

      // Ambient Occlusion
      float shadow = hs::clamp((abs_val * params.amplitude) / 0.4f, 0.0f, 1.0f);
      float occlusion = 0.15f + 0.85f * shadow;
      base.color = base.color * occlusion;

      frag.color = base;
    };

    Scan::rasterize<W, H>(filters, canvas, blob, shader, params.debug_bb);
  }

private:
  void start_morph() {
    // Pick a random new harmonic (Modes 1 to 24 look the best)
    next_idx = static_cast<int>(hs::rand_int(1, 24));

    // Animate morph_alpha 0->1 over 90 frames using Sine Easing
    timeline.add(
        0, Animation::Transition(morph_alpha, 1.0f, 64, ease_mid, false, false)
               .then([this]() {
                 // Commit the morph and start next one immediately
                 current_idx = next_idx;
                 morph_alpha = 0.0f;
                 start_morph();
               }));
  }

  // Params
  struct Params {
    float amplitude = 3.2f;
    bool debug_bb = false;
  } params;

  Orientation<W> orientation;
  Timeline<W> timeline;
  Pipeline<W, H> filters;
  BakedPalette baked_palette;

  int current_idx;
  int next_idx;
  float morph_alpha = 0.0f;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(SphericalHarmonics)
