/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W, int H>
class MobiusGrid : public Effect {
  mutable Fragments m_points;
  mutable Fragments m_fragments; // Reusable buffers
public:
  MobiusGrid() :
    Effect(W, H),
    palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
    holeN(Z_AXIS),
    holeS(-Z_AXIS),
    filters(
      Filter::World::HoleRef<W>(holeN, 1.2f),
      Filter::World::HoleRef<W>(holeS, 1.2f),
      Filter::World::Orient<W>(orientation),
      Filter::Screen::AntiAlias<W, H>()
    )
  {
    registerParam("Rings", &params.num_rings, 0.0f, 20.0f);
    registerParam("Lines", &params.num_lines, 0.0f, 20.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    persist_pixels = false;

    timeline
      .add(0, Animation::MobiusWarpCircular(mobius_params, 1.0f, 160, true))
      .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400, ease_mid, true))
      .add(0, Animation::PeriodicTimer(120, [this](auto&) { wipe_palette(); }, true))
      .add(0, Animation::Mutation(params.num_rings, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true))
      .add(160, Animation::Mutation(params.num_lines, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float phase = fmodf(static_cast<float>(timeline.t), 120.0f) / 120.0f;

    // Calculate Stabilizing Counter-Rotation
    Vector n_in = Y_AXIS;
    Vector n_trans = inv_stereo(mobius(stereo(n_in), mobius_params));
    Vector s_in = -Y_AXIS;
    Vector s_trans = inv_stereo(mobius(stereo(s_in), mobius_params));
    Vector mid = (n_trans + s_trans);
    Quaternion q;
    if (mid.length() > 0.001f) {
      mid.normalize();
      q = make_rotation(mid, Z_AXIS);
    }

    // Update hole origins to match the rotated geometry
    holeN = rotate(n_trans, q).normalize();
    holeS = rotate(s_trans, q).normalize();

    draw_axis_rings(canvas, Y_AXIS, params.num_rings, phase, q);
    draw_longitudes(canvas, params.num_lines, phase, q);
  }

private:
  void wipe_palette() {
    next_palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT);
    timeline.add(0, Animation::ColorWipe(palette, next_palette, 60, ease_mid));
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

      m_points.clear();
      Basis basis = make_basis(Quaternion(), normal);
      Plot::SphericalPolygon::sample(m_points, basis, radius, W / 4);

      m_fragments.clear();
      m_fragments.reserve(m_points.size());
      for (size_t k = 0; k < m_points.size(); ++k) {
        Vector transformed = inv_stereo(mobius(stereo(m_points[k].pos), mobius_params));
        Fragment f;
        f.pos = rotate(transformed, q).normalize();
        f.v0 = (float)k / (m_points.size() - 1);
        m_fragments.push_back(f);
      }

      float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
      
      auto fragment_shader = [&](const Vector&, Fragment& f_val) {
        Color4 c = palette.get(static_cast<float>(i) / num);
        c.alpha *= opacity * params.alpha;
        f_val.color = c;
      };

      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader, true);
    }
  }

  void draw_longitudes(Canvas& canvas, float num, float phase, const Quaternion& q) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      float theta = (static_cast<float>(i) / num) * PI_F;
      Vector normal(cosf(theta), 0.0f, -sinf(theta));

      m_points.clear();
      // Explicit basis construction to match JS texture alignment
      // v = normal (in XZ plane)
      // w = Y_AXIS (Pole)
      // u = cross(v, w) (Tangent)
      Vector v = normal;
      Vector w = Y_AXIS;
      Vector u = cross(v, w).normalize();
      Basis basis = { u, v, w };

      Plot::SphericalPolygon::sample(m_points, basis, 1.0f, W / 4);
      
      m_fragments.clear();
      m_fragments.reserve(m_points.size());
      for (size_t k = 0; k < m_points.size(); ++k) {
        Vector transformed = inv_stereo(mobius(stereo(m_points[k].pos), mobius_params));
        Fragment f;
        f.pos = rotate(transformed, q).normalize();
        
        // Correct v0 range to [0, 1]
        f.v0 = (m_points.size() > 1) ? (float)k / (m_points.size() - 1) : 0.0f; 
        m_fragments.push_back(f);
      }

      float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
      
      auto fragment_shader = [&](const Vector&, Fragment& f_val) {
        float t_line = f_val.v0;
        float z = sinf(t_line * 2.0f * PI_F);

        float R = sqrtf((1.0f + z) / (1.0f - z));
        float log_r = logf(R);
        const float log_min = -2.5f;
        const float log_max = 2.5f;
        float t = (log_r - log_min) / (log_max - log_min);

        Color4 c = palette.get(wrap(t - phase, 1.0f));
        c.alpha *= opacity * params.alpha;
        f_val.color = c;
      };

      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader, true);
    }
  }

  GenerativePalette palette;
  GenerativePalette next_palette;
  MobiusParams mobius_params;
  
  struct Params {
      float num_rings = 0.0f;
      float num_lines = 0.0f;
      float alpha = 0.2f;
  } params;

  Orientation<W> orientation;
  Timeline<W> timeline;


  Vector holeN;
  Vector holeS;

  Pipeline<W, H,
    Filter::World::HoleRef<W>,
    Filter::World::HoleRef<W>,
    Filter::World::Orient<W>,
    Filter::Screen::AntiAlias<W, H>
  > filters;
};