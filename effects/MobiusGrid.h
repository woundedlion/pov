/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W>
class MobiusGrid : public Effect {
public:
  MobiusGrid() :
    Effect(W),
    palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    holeN(Z_AXIS),
    holeS(-Z_AXIS),
    filters(
      FilterHoleRef<W>(holeN, 1.2f),
      FilterHoleRef<W>(holeS, 1.2f),
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>()
    )
  {
    persist_pixels = false;

    timeline
      .add(0, MobiusWarp(params, 1.0f, 160, true))
      .add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400, ease_mid, true))
      .add(0, PeriodicTimer(120, [this](auto&) { wipe_palette(); }, true))
      .add(0, Mutation(num_rings, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true))
      .add(160, Mutation(num_lines, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float phase = fmodf(static_cast<float>(timeline.t), 120.0f) / 120.0f;

    // Calculate Stabilizing Counter-Rotation
    Vector n_in = Z_AXIS;
    Vector n_trans = inv_stereo(mobius(stereo(n_in), params));
    Vector s_in = -Z_AXIS;
    Vector s_trans = inv_stereo(mobius(stereo(s_in), params));
    Vector mid = (n_trans + s_trans);
    Quaternion q; // Default Identity
    if (mid.length() > 0.001f) {
      mid.normalize();
      q = make_rotation(mid, Z_AXIS);
    }

    // Update hole origins to match the rotated geometry
    holeN = rotate(n_trans, q).normalize();
    holeS = rotate(s_trans, q).normalize();

    draw_axis_rings(canvas, Z_AXIS, num_rings, phase, q);
    draw_longitudes(canvas, num_lines, phase, q);
  }

private:
  void wipe_palette() {
    next_palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT);
    timeline.add(0, ColorWipe(palette, next_palette, 60, ease_mid));
  }

  void draw_axis_rings(Canvas& canvas, const Vector& normal, float num, float phase, const Quaternion& q) {
    const float log_min = -2.5f;
    const float log_max = 2.5f;
    const float range = log_max - log_min;
    int count = static_cast<int>(std::ceil(num));

    for (int i = 0; i < count; ++i) {
      float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
      float log_r = log_min + t * range;
      float r_val = expf(log_r);
      float radius = (4.0f / PI_F) * atanf(1.0f / r_val);

      Points points;
      Basis basis = make_basis(Quaternion(), normal);
      Plot::Polygon::sample<W>(points, basis, radius, W / 4);

      Fragments fragments; 
      fragments.reserve(points.size());
      for (size_t k = 0; k < points.size(); ++k) {
        Vector transformed = inv_stereo(mobius(stereo(points[k]), params));
        Fragment f;
        f.pos = rotate(transformed, q).normalize();
        f.v0 = (float)k / (points.size() - 1); // t
        fragments.push_back(f);
      }

      float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
      
      auto fragment_shader = [&](const Vector&, const Fragment&) {
        Color4 c = palette.get(static_cast<float>(i) / num);
        c.alpha *= opacity;
        return c;
      };

      Plot::rasterize<W>(filters, canvas, fragments, fragment_shader, true);
    }
  }

  void draw_longitudes(Canvas& canvas, float num, float phase, const Quaternion& q) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      float theta = (static_cast<float>(i) / num) * PI_F;
      Vector normal(cosf(theta), sinf(theta), 0.0f);

      Points points;
      Basis basis = make_basis(Quaternion(), normal);
      Plot::Polygon::sample<W>(points, basis, 1.0f, W / 4);
      
      Fragments fragments;
      fragments.reserve(points.size());
      for (size_t k = 0; k < points.size(); ++k) {
        Vector transformed = inv_stereo(mobius(stereo(points[k]), params));
        Fragment f;
        f.pos = rotate(transformed, q).normalize();
        f.v0 = (float)k / points.size(); // t (0..1)
        fragments.push_back(f);
      }

      float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
      
      auto fragment_shader = [&](const Vector&, const Fragment& f_val) {
        float t_line = f_val.v0;
        // Approximate original Z to calculate log-gradient
        float original_idx = t_line * points.size();
        int idx1 = static_cast<int>(original_idx) % points.size();
        int idx2 = (idx1 + 1) % points.size();
        float f = original_idx - std::floor(original_idx);
        float z = points[idx1].k * (1.0f - f) + points[idx2].k * f;

        float R = sqrtf((1.0f + z) / (1.0f - z));
        float log_r = logf(R);
        const float log_min = -2.5f;
        const float log_max = 2.5f;
        float t = (log_r - log_min) / (log_max - log_min);

        Color4 c = palette.get(wrap(t - phase, 1.0f));
        c.alpha *= opacity;
        return c;
      };

      Plot::rasterize<W>(filters, canvas, fragments, fragment_shader, true);
    }
  }

  float alpha = 0.2f;
  float num_rings = 0;
  float num_lines = 0;
  GenerativePalette palette;
  GenerativePalette next_palette;
  MobiusParams params;
  Orientation orientation;
  Timeline timeline;


  Vector holeN;
  Vector holeS;

  Pipeline<W,
    FilterHoleRef<W>,
    FilterHoleRef<W>,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;
};