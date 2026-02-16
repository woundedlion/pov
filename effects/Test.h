/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../solids.h"

template <int W, int H>
class Test : public Effect {
public:
  Test() :
    Effect(W, H),
    ringPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    polyPalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS, BrightnessProfile::CUP),
    normal(X_AXIS),
    debugBB(false)
  {
    params.thickness = 4.0f * (2.0f * PI_F / W);
    
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Amplitude", &params.amplitude, 0.0f, 2.0f);
    registerParam("Thickness", &params.thickness, 0.1f, 10.0f);
    registerParam("Rings", &params.numRings, 1.0f, 10.0f);

    this->persist_pixels = false;
    timeline.add(0,
      Animation::Sprite([this](Canvas& canvas, float opacity) { this->drawFn(canvas, opacity); }, -1, 48, ease_mid, 0, ease_mid)
    );

    timeline.add(0,
      Animation::RandomWalk<W>(orientation, normal)
    );
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  void drawFn(Canvas& canvas, float opacity) {
    int nRings = static_cast<int>(params.numRings);
    
    for (int i = 0; i < nRings; ++i) {
      float radius = 2.0f / (nRings + 1) * (i + 1);
      auto shiftFn = [this](float t) {
        return sin_wave(this->params.amplitude, -this->params.amplitude, 4.0f, 0.0f)(t);
      };
      
      // Use actual amplitude for bounding box calculations
      float amp = this->params.amplitude;

      Basis basis = make_basis(orientation.get(), normal);
      auto fragment_shader = [&](const Vector& p, Fragment& f) {
          // f.v0 is normalized azimuth (0..1)
          // f.v1 is distance from center line
          // f.size is thickness
          
          f.color = ringPalette.get(f.v0);
          
          float norm_dist = std::clamp(f.v1 / f.size, 0.0f, 1.0f);
          float t = 1.0f - norm_dist;
          float falloff = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
          
          f.color.alpha = f.color.alpha * opacity * params.alpha * falloff;
      };

      Scan::DistortedRing::draw<W, H>(filters, canvas, basis, radius, params.thickness,
        shiftFn, amp,
        fragment_shader,
        debugBB
      );
    }
  }

private:
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  
  struct Params {
      float alpha = 0.3f;
      float amplitude = 0.0f;
      float thickness = 1.0f; 
      float numRings = 1.0f;
  } params;

  GenerativePalette ringPalette;
  GenerativePalette polyPalette;
  Vector normal;
  Orientation<W> orientation;
  
  bool debugBB;
};
