/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Concentric rings warped by an animated sinusoidal radial displacement.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings are drawn in ring-space and warped by a sinusoidal radial
 * displacement whose amplitude is animated. Orientation random-walks over time
 * and each ring is shaded from a circular split-complementary palette.
 */
template <int W, int H> class DistortedRing : public Effect {
public:
  /**
   * @brief Builds the effect with palette, ring normal, and amplitude mutation.
   * @details Configures a circular split-complementary palette, an X-axis ring
   * normal, and an amplitude mutation that sweeps amplitude over
   * [-max_amplitude, +max_amplitude] as a unit sine wave (32 steps, eased,
   * looping).
   */
  FLASHMEM DistortedRing()
      : Effect(W, H),
        ringPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                    BrightnessProfile::FLAT),
        normal(X_AXIS),
        amplitude_mut(
            amplitude,
            [this](float t) {
              return sin_wave(-params.max_amplitude, params.max_amplitude, 1.0f,
                              0.0f)(t);
            },
            32, ease_linear, true) {}

  /**
   * @brief Registers params and builds the timeline.
   * @details Adds the ring sprite, orientation random walk, and amplitude
   * mutation animations to the timeline.
   */
  void init() override {
    // One pixel of azimuth in ring-space. Default stroke is 4 px, the Thickness
    // slider floor is 2 px and its cap is 12 px — all scaled by resolution so
    // the slider's pixel meaning (and the default's place within range) is the
    // same at any W. 12 px ≈ the prior fixed 0.75 rad cap at the canonical
    // Holosphere <96,20> device (12·2π/96 ≈ 0.785 rad), now consistent across
    // resolutions instead of a fixed angle whose pixel width grew with W.
    const float px = 2.0f * PI_F / W;
    params.thickness = 4.0f * px;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("MaxAmplitude", &params.max_amplitude, 0.0f, 2.0f);
    registerParam("Thickness", &params.thickness, 2.0f * px, 12.0f * px);
    registerParam("Rings", &params.num_rings, 1.0f, 10.0f);
    registerParam("Show Bounding", &params.debug_bb);

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1, 48, ease_linear, 0, ease_linear));

    timeline.add(0, Animation::RandomWalk<W>(orientation, normal, noise));

    timeline.add(0, amplitude_mut);
  }

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return false; each lit column persists in the POV sweep until the next
   *         column overwrites it, with no black strobe between columns.
   */
  bool strobe_columns() const override { return false; }

  /**
   * @brief Advances and renders the timeline for one frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /**
   * @brief Draws all rings for this frame.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details Ring radii are evenly spaced and the displacement wave amplitude
   * tracks the animated `amplitude` field.
   */
  void drawFn(Canvas &canvas, float opacity) {
    int nRings = static_cast<int>(params.num_rings);

    for (int i = 0; i < nRings; ++i) {
      float radius = 2.0f / (nRings + 1) * (i + 1);
      auto shiftFn = [this](float t) {
        return sin_wave(-amplitude, amplitude, 4.0f, 0.0f)(t);
      };

      Basis basis = make_basis(orientation.get(), normal);
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        // f.v0 is normalized azimuth (0..1)
        // f.v1 is distance from center line
        // f.size is thickness
        f.color = ringPalette.get(f.v0);

        float norm_dist = hs::clamp(f.v1 / f.size, 0.0f, 1.0f);
        // Smootherstep falloff: full alpha at the center line, fading to 0 at
        // the stroke edge.
        float falloff = quintic_kernel(1.0f - norm_dist);

        f.color.alpha = f.color.alpha * opacity * params.alpha * falloff;
      };

      Scan::DistortedRing::draw<W, H>(
          filters, canvas, basis, radius, params.thickness, shiftFn,
          std::abs(amplitude), fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;         /**< Overall ring opacity multiplier in [0, 1]. */
    float max_amplitude = 0.3f; /**< Peak radial displacement; the animated wave sweeps [-max_amplitude, +max_amplitude]. */
    float thickness = 1.0f;     /**< Ring stroke width, in radians of azimuth; init() reseeds it to 4 px scaled by resolution. */
    float num_rings = 1.0f;     /**< Number of evenly spaced concentric rings (truncated to int when drawn). */
    bool debug_bb = false;      /**< When true, draws each fragment's bounding box for debugging. */
  } params;

  float amplitude = 0;

  GenerativePalette ringPalette;
  Vector normal;
  Orientation<> orientation;
  // Built in the constructor (not inline in init() like the sprite/walk) so the
  // Mutation's reference_wrapper binding to `amplitude`/`params.max_amplitude`
  // is established once; init() then copies it into the timeline, and the copy
  // retains that binding. Held as a member only to be that copy source.
  Animation::Mutation amplitude_mut;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(DistortedRing)
