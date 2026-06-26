/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"
#include <cmath>
#include <utility>

namespace SHMath {
/**
 * @brief Factorial of n as a float.
 * @param n Non-negative integer whose factorial is computed; kept small.
 * @return n! as a float.
 * @details Float precision; large n would overflow float precision.
 */
inline float factorial(int n) {
  if (n <= 1)
    return 1.0f;
  float result = 1.0f;
  for (int i = 2; i <= n; i++)
    result *= i;
  return result;
}

/**
 * @brief Associated Legendre polynomial P_l^m(x).
 * @param l Degree (l >= 0).
 * @param m Order (0 <= m <= l).
 * @param x Argument, equal to cos(phi) with |x| <= 1.
 * @return Value of P_l^m(x).
 * @details Uses the standard upward recurrence: seed P_m^m, then P_{m+1}^m,
 * then recur in l.
 */
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

/**
 * @brief Precompute the spherical-harmonic normalization factor for (l, m).
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l].
 * @return Normalization factor N, constant per shape.
 */
inline float normalization(int l, int m) {
  int absM = std::abs(m);
  float N = sqrtf(((2.0f * l + 1.0f) / (4.0f * PI_F)) *
                  (factorial(l - absM) / factorial(l + absM)));
  return (m != 0) ? sqrtf(2.0f) * N : N;
}

/**
 * @brief Evaluate a real spherical harmonic with a precomputed norm.
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l]; sign selects cos/sin azimuthal factor.
 * @param theta Azimuthal angle in radians.
 * @param cos_phi Cosine of the polar angle, in [-1, 1].
 * @param N Precomputed normalization factor for (l, m).
 * @return The harmonic value.
 * @details Takes cos_phi directly to avoid a cosf(acos(x)) roundtrip.
 */
inline float sphericalHarmonic(int l, int m, float theta, float cos_phi, float N) {
  int absM = std::abs(m);
  float P = associatedLegendre(l, absM, cos_phi);
  return N * P * ((m > 0) ? fast_cosf(m * theta) : (m < 0) ? fast_sinf(absM * theta) : 1.0f);
}

/**
 * @brief Decode a flat harmonic index into its (l, m) pair.
 * @param idx Flat index where idx = l*l + l + m.
 * @return Pair {l, m} with l = floor(sqrt(idx)) and m in [-l, l].
 * @details Seeds l from a float sqrt, then snaps it to the exact integer floor:
 * sqrtf at (or just below) a perfect square can round to l-epsilon and truncate
 * to l-1, which would push m to +l — outside the level's valid [-l, l] band and
 * into the next level's order. The correction loops make l provably exact
 * (l*l <= idx < (l+1)*(l+1)), so the returned order is always valid. Cold path
 * (a few calls per frame).
 */
inline std::pair<int, int> decode_lm(int idx) {
  int l = static_cast<int>(sqrtf(static_cast<float>(idx)));
  while ((l + 1) * (l + 1) <= idx)
    ++l;
  while (l > 0 && l * l > idx)
    --l;
  return {l, idx - l * l - l};
}
} // namespace SHMath

/**
 * @brief Visualizer that paints real spherical harmonics on the unit sphere.
 * @tparam W Display width in pixels.
 * @tparam H Display height in pixels.
 * @details Continuously morphs between randomly chosen modes.
 */
template <int W, int H> class SphericalHarmonics : public Effect {
public:
  /**
   * @brief Field sampler that evaluates the (blended) harmonic at a world point.
   * @details Intentionally NOT a deforming SDF: on the unit-sphere display every
   * pixel samples r=1, so a lobe-radius distance would only modulate coverage
   * near the harmonic's nodes. We want full coverage and let the shader paint
   * the harmonic value, so distance() reports a constant "fully inside" and the
   * field value rides out through raw_dist (frag.v1).
   */
  struct HarmonicField {
    int l1, m1;
    int l2, m2;
    float blend;
    Quaternion orientation;
    float N1, N2; /**< Normalization factors, precomputed once per shape. */
    static constexpr bool is_solid = true;

    /**
     * @brief Construct a field for blending mode (l1, m1) into (l2, m2).
     * @param l1 Degree of the first harmonic.
     * @param m1 Order of the first harmonic.
     * @param l2 Degree of the second harmonic.
     * @param m2 Order of the second harmonic.
     * @param blend Morph fraction in [0, 1] from the first toward the second.
     * @param q Orientation quaternion of the shape.
     */
    HarmonicField(int l1, int m1, int l2, int m2, float blend, Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend), orientation(q),
          N1(SHMath::normalization(l1, m1)),
          N2(SHMath::normalization(l2, m2)) {}

    /**
     * @brief Vertical scan bounds for the field.
     * @tparam H_ Display height in pixels.
     * @return Full-height bounds; lobes can occupy any region, so no static
     * bounding.
     */
    template <int H_> SDF::Bounds get_vertical_bounds() const {
      return {0, H_ - 1};
    }
    /**
     * @brief Horizontal scan intervals for a given row.
     * @tparam W_scan Display width in pixels.
     * @tparam H_ Display height in pixels.
     * @tparam OutputIt Output iterator type for emitted intervals.
     * @return Always false: no horizontal interval narrowing (full-sphere scan).
     */
    template <int W_scan, int H_, typename OutputIt>
    bool get_horizontal_intervals(int, OutputIt) const {
      return false;
    }

    /**
     * @brief Convenience wrapper returning a fresh DistanceResult.
     * @param p World-space sample point.
     * @return The sampled DistanceResult for p.
     */
    SDF::DistanceResult distance(const Vector &p) const {
      SDF::DistanceResult res;
      distance<true>(p, res);
      return res;
    }

    /**
     * @brief Sample the (possibly blended) harmonic at world point p.
     * @tparam ComputeUVs Whether to compute UV coordinates (unused here).
     * @param p World-space sample point.
     * @param res Output DistanceResult; carries the field value in v1 with a
     * constant "inside" distance.
     * @details Rotates p into the shape's local frame, evaluates both modes in
     * spherical coords, and emits the field value through frag.v1.
     */
    template <bool ComputeUVs = true>
    void distance(const Vector &p, SDF::DistanceResult &res) const {
      Vector local = rotate(p, orientation.conjugate());

      float theta = fast_atan2(local.z, local.x);
      if (theta < 0)
        theta += 2 * PI_F;
      // local.y is cos(phi) in the LOCAL frame: the shape spins about an
      // arbitrary axis, so cos_phi varies across a screen row even though the
      // WORLD latitude is row-constant. Hoisting the associatedLegendre(cos_phi)
      // term to a per-row precompute is therefore invalid.
      float cos_phi = hs::clamp(local.y, -1.0f, 1.0f);

      // Per-pixel cost: associatedLegendre runs once here, and again for the
      // blend target, with no horizontal-interval narrowing
      // (get_horizontal_intervals returns false) and full-height vertical
      // bounds — the heaviest per-pixel evaluation in the harmonic group. Left
      // uncapped because morph indices are drawn from [1, 24] (start_morph), so
      // l = floor(sqrt(idx)) <= 4 keeps the Legendre recurrence short and the
      // cost bounded.
      float val = SHMath::sphericalHarmonic(l1, m1, theta, cos_phi, N1);
      if (blend > 0.001f) {
        float val2 = SHMath::sphericalHarmonic(l2, m2, theta, cos_phi, N2);
        val += (val2 - val) * blend;
      }

      // Constant "inside" (see struct comment); field value -> frag.v1.
      res = SDF::DistanceResult(-1.0f, 0.0f, val, 0.0f, 1.0f);
    }
  };

  /**
   * @brief Construct the visualizer with the display dimensions.
   */
  FLASHMEM SphericalHarmonics() : Effect(W, H), filters() {}

  /**
   * @brief One-time setup of params, palette, shape, spin, and first morph.
   * @details Registers params, bakes the palette, seeds the shape and the
   * continuous spin, and kicks off the first morph.
   */
  void init() override {
    registerParam("Amplitude", &params.amplitude, 0.1f, 10.0f);
    registerParam("Debug BB", &params.debug_bb);

    baked_palette.bake(persistent_arena, Palettes::richSunset);

    current_idx = SEED_MODE_IDX;

    Vector axis = Vector(0.5f, 1.0f, 0.2f).normalized();
    timeline.add(0, Animation::Rotation<W>(orientation, axis, 2 * PI_F * 100,
                                           10000, ease_linear, true));

    start_morph();
  }

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return true; the strip blanks to black immediately after each column is
   *         shown, so every column reads as a sharp slice with dark gaps
   *         between columns rather than persisting across the sweep.
   */
  bool strobe_columns() const override { return true; }

  /**
   * @brief Render one frame of the morphing harmonic.
   * @details Decodes the current and target modes, builds the field for this
   * frame's morph state, and rasterizes with the harmonic-coloring shader.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    auto [l1, m1] = SHMath::decode_lm(current_idx);
    auto [l2, m2] = SHMath::decode_lm(next_idx);
    HarmonicField field(l1, m1, l2, m2, morph_alpha, orientation.get());

    // Map the field value (carried in frag.v1) to color: diverging positive/
    // negative palettes crossfaded at the zero-crossing, then AO brightness.
    auto shader = [&](const Vector &, Fragment &frag) {
      float val = frag.v1;
      float abs_val = std::abs(val);

      Color4 pos =
          baked_palette.get(std::min(1.0f, abs_val * params.amplitude));
      // Negative palette: recolor pos by swapping R<->B and dimming green so
      // negative SH lobes read distinct. A stylized polarity cue, not OKLCH.
      constexpr float NEG_LOBE_GREEN_SCALE = 0.8f;
      Color4 neg = Color4(
          Pixel(pos.color.b,
                static_cast<uint16_t>(pos.color.g * NEG_LOBE_GREEN_SCALE),
                pos.color.r),
          pos.alpha);

      // Narrow quintic-smoothed seam at the zero-crossing; `transition` is the
      // half-width (field units) of the anti-aliasing band. blend_t is inverted
      // (1 - blend_t) below so a positive field maps to pos and a negative to neg.
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
  /**
   * @brief Ambient-occlusion shaping constants for the draw_frame shader.
   * @details AO_FALLOFF is the field magnitude at which shading saturates,
   * AO_AMBIENT is the ambient floor, and AO_RANGE is the range above it
   * (AMBIENT + RANGE = 1.0 → full brightness at saturation).
   */
  static constexpr float AO_FALLOFF = 0.4f;
  static constexpr float AO_AMBIENT = 0.15f;
  static constexpr float AO_RANGE = 0.85f;

  /**
   * @brief Choose the next harmonic and animate the morph toward it.
   * @details On completion commits it as current and recurses, yielding an
   * endless morph chain.
   */
  void start_morph() {
    // Re-roll on a match with current_idx: blending a basis function into
    // itself would freeze the sphere for the whole transition.
    do {
      next_idx = hs::rand_int(1, MAX_MODE_IDX + 1);
    } while (next_idx == current_idx);

    timeline.add(
        0, Animation::Transition(morph_alpha, 1.0f, 64, ease_linear, false, false)
               .then([this]() {
                 current_idx = next_idx;
                 morph_alpha = 0.0f;
                 start_morph();
               }));
  }

  /**
   * @brief User-tunable parameters for the visualizer.
   */
  struct Params {
    float amplitude = 3.2f; /**< Field-value gain applied before coloring. */
    bool debug_bb = false;  /**< Whether to draw bounding-box debug overlay. */
  } params;

  Orientation<> orientation;  /**< Current sphere orientation. */
  Timeline timeline;          /**< Drives spin and morph animations. */
  Pipeline<W, H> filters;     /**< Post-process filter pipeline. */
  BakedPalette baked_palette; /**< Precomputed color LUT for the shader. */

  // Highest harmonic degree the morph visits. Modes are flat-indexed
  // idx = l*l + l + m, so degrees [0, MAX_DEGREE] occupy idx [0, MAX_MODE_IDX].
  // The seed and the roll bound below derive from MAX_DEGREE so they cannot
  // drift from it (widen MAX_DEGREE alone to extend the visual range).
  static constexpr int MAX_DEGREE = 4;
  // Top flat index over those degrees: idx peaks at l = MAX_DEGREE, m = +MAX_DEGREE.
  static constexpr int MAX_MODE_IDX = (MAX_DEGREE + 1) * (MAX_DEGREE + 1) - 1; // 24
  // Initial mode (l=2, m=0); the constant mode (idx 0) is excluded from the roll.
  static constexpr int SEED_MODE_IDX = 6;
  static_assert(SEED_MODE_IDX > 0 && SEED_MODE_IDX <= MAX_MODE_IDX,
                "seed mode must be a valid, non-constant harmonic index");

  int current_idx = 0;        /**< Flat index of the currently displayed mode. */
  int next_idx = 0;           /**< Flat index of the mode being morphed toward. */
  float morph_alpha = 0.0f;   /**< Morph progress in [0, 1] from current to next. */
};

#include "core/effect_registry.h"
REGISTER_EFFECT(SphericalHarmonics)
