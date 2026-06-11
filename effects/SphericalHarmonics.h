/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include <cmath>
#include <utility>

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

// Flat harmonic index -> (l, m): idx = l*l + l + m, so l = floor(sqrt(idx)) and
// m = idx - l*l - l lands in [-l, l]. Seed l from a float sqrt, then snap it to
// the exact integer floor: sqrtf at (or just below) a perfect square can round
// to l-epsilon and truncate to l-1, which would push m to +l — outside the
// level's valid [-l, l] band and into the next level's order. The correction
// loops make l provably exact (l*l <= idx < (l+1)*(l+1)), so the returned order
// is always valid. Cold path (a few calls per frame).
inline std::pair<int, int> decode_lm(int idx) {
  int l = static_cast<int>(sqrtf(static_cast<float>(idx)));
  while ((l + 1) * (l + 1) <= idx)
    ++l;
  while (l > 0 && l * l > idx)
    --l;
  return {l, idx - l * l - l};
}
} // namespace SHMath

template <int W, int H> class SphericalHarmonics : public Effect {
public:
  // Field sampler for the visualizer. Intentionally NOT a deforming SDF: on the
  // unit-sphere display every pixel samples r=1, so a lobe-radius distance would
  // only modulate coverage near the harmonic's nodes. We want full coverage and
  // let the shader paint the harmonic value, so distance() reports a constant
  // "fully inside" and the field value rides out through raw_dist (frag.v1).
  struct HarmonicField {
    int l1, m1;
    int l2, m2;
    float blend;
    Quaternion orientation;
    float N1, N2; // Normalization factors, precomputed once per shape.
    static constexpr bool is_solid = true;

    HarmonicField(int l1, int m1, int l2, int m2, float blend, Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend), orientation(q),
          N1(SHMath::normalization(l1, m1)),
          N2(SHMath::normalization(l2, m2)) {}

    // Full-sphere scan: lobes can occupy any region, so no static bounding.
    template <int H_> SDF::Bounds get_vertical_bounds() const {
      return {0, H_ - 1};
    }
    template <int W_scan, int H_, typename OutputIt>
    bool get_horizontal_intervals(int, OutputIt) const {
      return false;
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
      // local.y IS cos(phi) — skip the fast_acos + cosf roundtrip.
      float cos_phi = hs::clamp(local.y, -1.0f, 1.0f);

      float val = SHMath::sphericalHarmonic(l1, m1, theta, cos_phi, N1);
      // Only pay for the second harmonic while actively blending.
      if (blend > 0.001f) {
        float val2 = SHMath::sphericalHarmonic(l2, m2, theta, cos_phi, N2);
        val += (val2 - val) * blend;
      }

      // Constant "inside" (see struct comment); field value -> frag.v1.
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

    auto [l1, m1] = SHMath::decode_lm(current_idx);
    auto [l2, m2] = SHMath::decode_lm(next_idx);
    HarmonicField field(l1, m1, l2, m2, morph_alpha, orientation.get());

    auto shader = [&](const Vector &p, Fragment &frag) {
      float val = frag.v1;
      float abs_val = std::abs(val);

      // Positive palette (baked LUT — no cosf/powf per pixel)
      Color4 pos =
          baked_palette.get(std::min(1.0f, abs_val * params.amplitude));
      // Negative palette (channel-swapped)
      Color4 neg =
          Color4(Pixel(pos.color.b, static_cast<uint16_t>(pos.color.g * 0.8f),
                       pos.color.r),
                 pos.alpha);

      // Quintic-smoothed crossfade across the zero-crossing boundary
      constexpr float transition = 0.03f;
      float blend_t =
          quintic_kernel(hs::clamp(val / transition * 0.5f + 0.5f, 0.0f, 1.0f));
      Color4 base = pos.lerp(neg, 1.0f - blend_t);

      // Ambient Occlusion: scale brightness from an ambient floor up to full as
      // the (amplitude-weighted) field magnitude saturates.
      float shadow =
          hs::clamp((abs_val * params.amplitude) / AO_FALLOFF, 0.0f, 1.0f);
      float occlusion = AO_AMBIENT + AO_RANGE * shadow;
      base.color = base.color * occlusion;

      frag.color = base;
    };

    Scan::rasterize<W, H>(filters, canvas, field, shader, params.debug_bb);
  }

private:
  // Ambient-occlusion shaping (see draw_frame shader): field magnitude at which
  // shading saturates, the ambient floor, and the range above it (AMBIENT +
  // RANGE = 1.0 → full brightness at saturation).
  static constexpr float AO_FALLOFF = 0.4f;
  static constexpr float AO_AMBIENT = 0.15f;
  static constexpr float AO_RANGE = 0.85f;

  void start_morph() {
    // Pick a random new harmonic (low modes 1..23 look the best)
    next_idx = static_cast<int>(hs::rand_int(1, 24));

    // Animate morph_alpha 0->1 over 64 frames with linear easing
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

  Orientation<> orientation;
  Timeline timeline;
  Pipeline<W, H> filters;
  BakedPalette baked_palette;

  int current_idx;
  int next_idx;
  float morph_alpha = 0.0f;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(SphericalHarmonics)
