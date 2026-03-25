/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

template <int W, int H> class Test : public Effect {
public:
  FLASHMEM Test()
      : Effect(W, H),
        ringPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                    BrightnessProfile::FLAT),
        polyPalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                    BrightnessProfile::CUP),
        normal(X_AXIS),
        amplitude_mut(
            amplitude,
            [this](float t) {
              return sin_wave(-params.max_amplitude, params.max_amplitude, 1.0f,
                              0.0f)(t);
            },
            32, ease_mid, true) {}

  void init() override {
    params.thickness = 4.0f * (2.0f * PI_F / W);

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("MaxAmplitude", &params.max_amplitude, 0.0f, 2.0f);
    registerParam("Thickness", &params.thickness, 0.1f, 0.75f);
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

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Refresh Params
  }

  void drawFn(Canvas &canvas, float opacity) {
    int nRings = static_cast<int>(params.numRings);

    for (int i = 0; i < nRings; ++i) {
      float radius = 2.0f / (nRings + 1) * (i + 1);
      auto shiftFn = [this](float t) {
        return sin_wave(-amplitude, amplitude, 4.0f, 0.0f)(t);
      };

      Basis basis = make_basis(orientation.get(), normal);
      auto fragment_shader = [&](const Vector &p, Fragment &f) {
        // f.v0 is normalized azimuth (0..1)
        // f.v1 is distance from center line
        // f.size is thickness
        f.color = ringPalette.get(f.v0);

        float norm_dist = hs::clamp(f.v1 / f.size, 0.0f, 1.0f);
        float t = 1.0f - norm_dist;
        float falloff = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);

        f.color.alpha = f.color.alpha * opacity * params.alpha * falloff;
      };

      Scan::DistortedRing::draw<W, H>(
          filters, canvas, basis, radius, params.thickness, shiftFn,
          std::abs(amplitude), fragment_shader, 0.0f, params.debugBB);
    }
  }

private:
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H> filters;

  struct Params {
    float alpha = 0.3f;
    float max_amplitude = 0.3f;
    float thickness = 1.0f;
    float numRings = 1.0f;
    bool debugBB = false;
  } params;

  float amplitude = 0;

  GenerativePalette ringPalette;
  GenerativePalette polyPalette;
  Vector normal;
  Orientation<W> orientation;
  Animation::Mutation amplitude_mut;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Test)
