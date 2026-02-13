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
    alpha(0.3f),
    ringPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    polyPalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS, BrightnessProfile::CUP),
    normal(X_AXIS),
    amplitude(0.0f),
    amplitudeRange(0.3f),
    numRings(1),
    debugBB(false),
    thickness(4.0f * (2.0f * PI_F / W))
  {
    timeline.add(0,
      Animation::Sprite([this](Canvas& canvas, float opacity) { this->drawFn(canvas, opacity); }, -1, 48, ease_mid, 0, ease_mid)
    );

    timeline.add(0,
      Animation::RandomWalk<W>(orientation, normal)
    );

    timeline.add(0,
      Animation::Mutation(amplitude,
        sin_wave(-amplitudeRange, amplitudeRange, 1.0f, 0.0f), 32, ease_mid, true)
    );
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  void drawFn(Canvas& canvas, float opacity) {
    for (int i = 0; i < numRings; ++i) {
      float radius = 2.0f / (numRings + 1) * (i + 1);
      auto shiftFn = [this](float t) {
        return sin_wave(this->amplitude, -this->amplitude, 4.0f, 0.0f)(t);
      };
      float amp = this->amplitudeRange;

      Basis basis = make_basis(orientation.get(), normal);
      auto fragment_shader = [&](const Vector& p, const Fragment& f) -> Fragment {
          // f.v0 is normalized azimuth (0..1)
          Fragment out = f;
          out.color = ringPalette.get(f.v0);
          out.color.alpha = out.color.alpha * opacity * alpha;
          return out;
      };

      Scan::DistortedRing::draw<W, H>(filters, canvas, basis, radius, thickness,
        shiftFn, amp,
        fragment_shader,
        debugBB
      );
    }
  }

private:
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  
  float alpha;
  GenerativePalette ringPalette;
  GenerativePalette polyPalette;
  Vector normal;
  Orientation<W> orientation;
  
  float amplitude;
  float amplitudeRange;
  int numRings;
  bool debugBB;
  float thickness;
};
