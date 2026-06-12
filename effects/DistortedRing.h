/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

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
            32, ease_mid, true) {}

  /**
   * @brief Registers params and builds the timeline.
   * @details Adds the ring sprite, orientation random walk, and amplitude
   * mutation animations to the timeline.
   */
  void init() override {
    // One pixel of azimuth in ring-space. Default stroke is 4 px and the
    // Thickness slider floor is 2 px, both scaled by resolution so the default
    // stays within slider range at any W.
    const float px = 2.0f * PI_F / W;
    params.thickness = 4.0f * px;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("MaxAmplitude", &params.max_amplitude, 0.0f, 2.0f);
    registerParam("Thickness", &params.thickness, 2.0f * px, 0.75f);
    registerParam("Rings", &params.numRings, 1.0f, 10.0f);
    registerParam("Show Bounding", &params.debugBB);

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1, 48, ease_mid, 0, ease_mid));

    timeline.add(0, Animation::RandomWalk<W>(orientation, normal, noise));

    timeline.add(0, amplitude_mut);
  }

  /**
   * @brief Reports whether the engine should clear the background each frame.
   * @return false; rings draw over the prior frame with no background clear.
   */
  bool show_bg() const override { return false; }

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
    int nRings = static_cast<int>(params.numRings);

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
          std::abs(amplitude), fragment_shader, 0.0f, params.debugBB);
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
    float alpha = 0.3f;
    float max_amplitude = 0.3f;
    float thickness = 1.0f;
    float numRings = 1.0f;
    bool debugBB = false;
  } params;

  float amplitude = 0;

  GenerativePalette ringPalette;
  Vector normal;
  Orientation<> orientation;
  Animation::Mutation amplitude_mut;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(DistortedRing)
