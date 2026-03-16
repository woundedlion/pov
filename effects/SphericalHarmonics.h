/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
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

inline float sphericalHarmonic(int l, int m, float theta, float phi) {
  int absM = std::abs(m);
  float N = sqrtf(((2.0f * l + 1.0f) / (4.0f * PI_F)) *
                  (factorial(l - absM) / factorial(l + absM)));
  float P = associatedLegendre(l, absM, cosf(phi));

  if (m > 0) {
    return sqrtf(2.0f) * N * P * cosf(m * theta);
  } else if (m < 0) {
    return sqrtf(2.0f) * N * P * sinf(absM * theta);
  } else {
    return N * P;
  }
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

    HarmonicBlob(int l1, int m1, int l2, int m2, float blend, float amp,
                 Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend), amplitude(amp),
          orientation(q) {}

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

      float theta = atan2f(local.z, local.x);
      if (theta < 0)
        theta += 2 * PI_F;
      float phi = acosf(hs::clamp(local.y, -1.0f, 1.0f));

      // Evaluate primary shape
      float val1 = SHMath::sphericalHarmonic(l1, m1, theta, phi);
      float val = val1;

      // Only pay the math cost for the second shape if we are actively blending
      if (blend > 0.001f) {
        float val2 = SHMath::sphericalHarmonic(l2, m2, theta, phi);
        val = val1 + (val2 - val1) * blend;
      }

      res = SDF::DistanceResult(-1.0f, 0.0f, val, 0.0f, 1.0f);
    }
  };

  FLASHMEM SphericalHarmonics() : Effect(W, H), filters() {}

  void init() override {
    registerParam("Amplitude", &params.amplitude, 0.1f, 10.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Initial shape
    current_idx = 6;
    next_idx = 6;

    // Spin
    Vector axis = Vector(0.5f, 1.0f, 0.2f).normalize();
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
      float abs_val = std::abs(frag.v1);

      Color4 base;
      if (frag.v1 >= 0) {
        base = Palettes::richSunset.get(
            std::min(1.0f, abs_val * params.amplitude));
      } else {
        Color4 p = Palettes::richSunset.get(
            std::min(1.0f, abs_val * params.amplitude));
        base = Color4(
            Pixel(p.color.b, static_cast<uint8_t>(p.color.g * 0.8f), p.color.r),
            p.alpha);
      }

      // Ambient Occlusion
      float shadow =
          hs::clamp((abs_val * params.amplitude - 0.0f) / 0.4f, 0.0f, 1.0f);
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
        0, Animation::Transition(morph_alpha, 1.0f, 90, ease_mid, false, false)
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

  int current_idx;
  int next_idx;
  float morph_alpha = 0.0f;
};

#include "../effect_registry.h"
REGISTER_EFFECT(SphericalHarmonics)
