/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file SphericalHarmonics.h
 * @brief Real spherical harmonics painted on the unit sphere, continuously
 *        morphing between modes.
 */

#include "core/engine/engine.h"
#include "core/math/spherical_harmonics.h"
#include <cmath>
#include <utility>

// Unit-test accessor reaching the private morph-chain state and the baked
// palette the shader reads.
namespace hs_test {
namespace effects_tests {
struct SphericalHarmonicsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

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
   * @details Not a deforming SDF: distance() reports a constant "fully inside" so
   * the whole sphere is covered, and the field value rides out through raw_dist
   * (frag.v1) for the shader to paint.
   */
  struct HarmonicField {
    int l1, m1;
    int l2, m2;
    float blend;
    RotationMatrix orientation_conj; /**< World->local rotation (conjugate of
                                        the shape orientation). */
    float N1, N2; /**< Per-mode harmonic scales, precomputed once per shape. */
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
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend),
          orientation_conj(q.conjugate()), N1(SHMath::harmonic_scale(l1, m1)),
          N2(SHMath::harmonic_scale(l2, m2)) {}

    /**
     * @brief Vertical scan bounds for the field.
     * @tparam H_scan Display height in pixels.
     * @return Full-height bounds; lobes can occupy any region, so no static
     * bounding.
     */
    template <int H_scan> SDF::Bounds get_vertical_bounds() const {
      return {0, H_scan - 1};
    }
    /**
     * @brief Horizontal scan intervals for a given row.
     * @tparam W_scan Display width in pixels.
     * @tparam H_scan Display height in pixels.
     * @tparam OutputIt Output iterator type for emitted intervals.
     * @return Always false: no horizontal interval narrowing (full-sphere scan).
     */
    template <int W_scan, int H_scan, typename OutputIt>
    bool get_horizontal_intervals(int, OutputIt) const {
      return false;
    }

    /**
     * @brief Sample the (possibly blended) harmonic at world point p.
     * @tparam ComputeUVs Whether to compute UV coordinates (unused here).
     * @param p World-space sample point.
     * @param res Output DistanceResult; carries the field value in v1 with a
     * constant "inside" distance.
     * @details Rotates p into the shape's local frame, evaluates both modes
     * there, and emits the field value through frag.v1.
     */
    template <bool ComputeUVs = true>
    void distance(const Vector &p, SDF::DistanceResult &res) const {
      Vector local = orientation_conj.apply(p);

      // The shape spins about an arbitrary axis, so the local frame varies
      // across a screen row even though the WORLD latitude is row-constant.
      // Hoisting any part of the harmonic to a per-row precompute is invalid.
      float val = SHMath::spherical_harmonic(l1, m1, local, N1);
      if (blend > 0.001f) {
        float val2 = SHMath::spherical_harmonic(l2, m2, local, N2);
        val += (val2 - val) * blend;
      }

      // Constant "inside" (see struct comment); field value -> frag.v1.
      res = SDF::DistanceResult(-1.0f, 0.0f, val, 0.0f, 1.0f);
    }
  };

  /**
   * @brief Construct the visualizer with the display dimensions.
   */
  HS_COLD_MEMBER SphericalHarmonics()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})) {}

  /**
   * @brief One-time setup of params, palette, shape, spin, and first morph.
   * @details Registers params, bakes the palette, seeds the shape and the
   * continuous spin, and kicks off the first morph.
   */
  void init() override {
    register_param("Amplitude", &params.amplitude, 0.1f, 10.0f);
    register_param("Debug BB", &params.debug_bb);

    baked_palette.bake(persistent_arena, Palettes::RICH_SUNSET);

    current_idx = SEED_MODE_IDX;

    Vector axis = Vector(0.5f, 1.0f, 0.2f).normalized();
    timeline.add(0, Animation::Rotation<W>(orientation, axis, 2 * PI_F * 100,
                                           10000, ease_linear, true));

    start_morph();
  }

  /**
   * @brief Render one frame of the morphing harmonic.
   * @details Decodes the current and target modes, builds the field for this
   * frame's morph state, and rasterizes with the harmonic-coloring shader.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sh_timeline_step);
      timeline.step(canvas);
    }

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
      Color4 neg = Color4(Pixel(pos.color.b,
                                static_cast<uint16_t>(
                                    pos.color.g * NEG_LOBE_GREEN_SCALE + 0.5f),
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

    {
      HS_PROFILE(sh_rasterize);
      Scan::rasterize<W, H>(filters, canvas, field, shader, params.debug_bb);
    }
  }

private:
  friend struct ::hs_test::effects_tests::SphericalHarmonicsWhiteBox;

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
   * endless morph chain. Pausable: the chain is this effect's look
   * choreography, so "Pause Animation" holds the current blend.
   */
  HS_COLD_MEMBER void start_morph() {
#ifdef HS_PROFILE_ORDERED_CYCLE
    next_idx = (current_idx % MAX_MODE_IDX) + 1;
#else
    // Re-roll on a match with current_idx: blending a basis function into
    // itself would freeze the sphere for the whole transition.
    do {
      next_idx = hs::rand_int(1, MAX_MODE_IDX + 1);
    } while (next_idx == current_idx);
#endif
    hs::log("Mode: %d/%d", next_idx, MAX_MODE_IDX);

    timeline.add_pausable(
        0,
        Animation::Transition(morph_alpha, 1.0f, 64, ease_linear, false, false)
            .then([this]() {
              current_idx = next_idx;
              morph_alpha = 0.0f;
              start_morph();
            }),
        &anims_paused);
  }

  Orientation<> orientation;  /**< Current sphere orientation. */
  Timeline timeline;          /**< Drives spin and morph animations. */
  Pipeline<W, H> filters;     /**< Post-process filter pipeline. */
  BakedPalette baked_palette; /**< Precomputed color LUT for the shader. */

  // init() bakes one palette LUT into the persistent arena. Effect keeps the
  // default arena split, so the total must fit the device persistent partition.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes();
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "SphericalHarmonics persistent footprint exceeds the default partition");

  // Highest harmonic degree the morph visits. Modes are flat-indexed
  // idx = l*l + l + m, so degrees [0, MAX_DEGREE] occupy idx [0, MAX_MODE_IDX].
  static constexpr int MAX_DEGREE = 4;
  // normalization() divides factorials of up to 2*MAX_DEGREE; float holds exact
  // integers only up to 2^24, first exceeded at 11!.
  static_assert(2 * MAX_DEGREE <= 10,
                "MAX_DEGREE too large: factorial(2*MAX_DEGREE) exceeds float "
                "precision in normalization()");
  // Top flat index over those degrees: idx peaks at l = MAX_DEGREE, m = +MAX_DEGREE.
  static constexpr int MAX_MODE_IDX =
      (MAX_DEGREE + 1) * (MAX_DEGREE + 1) - 1; // 24
  // Initial mode (l=2, m=0); the constant mode (idx 0) is excluded from the roll.
  static constexpr int SEED_MODE_IDX = 6;
  static_assert(SEED_MODE_IDX > 0 && SEED_MODE_IDX <= MAX_MODE_IDX,
                "seed mode must be a valid, non-constant harmonic index");

  int current_idx = 0; /**< Flat index of the currently displayed mode. */
  int next_idx = 0;    /**< Flat index of the mode being morphed toward. */
  float morph_alpha =
      0.0f; /**< Morph progress in [0, 1] from current to next. */

  /**
   * @brief User-tunable parameters for the visualizer.
   */
  struct Params {
    float amplitude = 3.2f; /**< Field-value gain applied before coloring. */
    bool debug_bb = false;  /**< Whether to draw bounding-box debug overlay. */
  } params;
};

#include "core/control/registry.h"
REGISTER_EFFECT(SphericalHarmonics)
