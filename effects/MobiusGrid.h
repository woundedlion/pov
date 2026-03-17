/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include <functional>
#include <memory>

#include <map>
#include "../effects_engine.h"

template <int W, int H> class MobiusGrid : public Effect {

public:
  FLASHMEM MobiusGrid()
      : Effect(W, H),
        palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                BrightnessProfile::FLAT),
        next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::FLAT),
        mobius_gen(timeline), holeN(Z_AXIS), holeS(-Z_AXIS),
        filters(Filter::World::HoleRef<W>(holeN, 1.2f),
                Filter::World::HoleRef<W>(holeS, 1.2f),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  bool show_bg() const override { return false; }

  void init() override {
    // MobiusGrid requires very little persistent memory.
    // Give it 64KB for Scratch A and 64KB for Scratch B to comfortably handle rasterization arrays.
    configure_arenas(GLOBAL_ARENA_SIZE - 128 * 1024, 64 * 1024, 64 * 1024);

    registerParam("Rings", &params.num_rings, 0.0f, 20.0f);
    registerParam("Lines", &params.num_lines, 0.0f, 20.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Use Generator
    mobius_gen.spawn(0, 1.0f, 160, true);

    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400,
                                       ease_mid, true))
        .add(0, Animation::PeriodicTimer(
                    120, [this](auto &) { wipe_palette(); }, true))
        .add(0, Animation::Mutation(params.num_rings,
                                    sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320,
                                    ease_mid, true))
        .add(160, Animation::Mutation(params.num_lines,
                                      sin_wave(12.0f, 1.0f, 1.0f, PI_F / 2.0f),
                                      320, ease_mid, true));
  }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float phase = fmodf(static_cast<float>(timeline.t), 120.0f) / 120.0f;

    // Calculate Stabilizing Counter-Rotation
    Vector n_in = Y_AXIS;
    Vector n_trans = mobius_gen.transform(n_in);
    Vector s_in = -Y_AXIS;
    Vector s_trans = mobius_gen.transform(s_in);
    Vector mid = (n_trans + s_trans);
    Quaternion q;
    if (mid.length() > 0.001f) {
      mid.normalize();
      q = make_rotation(mid, Z_AXIS);
    }

    // Update hole origins to match the rotated geometry
    holeN = rotate(n_trans, q).normalized();
    holeS = rotate(s_trans, q).normalized();

    draw_axis_rings(canvas, Y_AXIS, params.num_rings, phase, q);
    draw_longitudes(canvas, params.num_lines, phase, q);
  }

private:
  void wipe_palette() {
    next_palette = GenerativePalette(GradientShape::CIRCULAR,
                                     HarmonyType::SPLIT_COMPLEMENTARY,
                                     BrightnessProfile::FLAT);
    timeline.add(0, Animation::ColorWipe(palette, next_palette, 60, ease_mid));
  }

  void draw_axis_rings(Canvas &canvas, const Vector &normal, float num,
                       float phase, const Quaternion &q) {
    const float log_min = -2.5f;
    const float log_max = 2.5f;
    const float range = log_max - log_min;
    int count = static_cast<int>(std::ceil(num));

    for (int i = 0; i < count; ++i) {
      ScratchScope _frag(scratch_arena_a);
      float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
      float log_r = log_min + t * range;
      float r_val = expf(log_r);
      float radius = (4.0f / PI_F) * atanf(1.0f / r_val);

      Fragments m_points;
      m_points.bind(scratch_arena_a, 144);
      Basis basis = make_basis(Quaternion(), normal);
      Plot::SphericalPolygon::sample(m_points, basis, radius, W / 4);

      Fragments m_fragments;
      m_fragments.bind(scratch_arena_a, 144);
      for (size_t k = 0; k < m_points.size(); ++k) {
        Vector transformed = mobius_gen.transform(m_points[k].pos);
        Fragment f;
        f.pos = rotate(transformed, q).normalized();
        f.v0 = (float)k / (m_points.size() - 1);
        m_fragments.push_back(f);
      }

      float opacity = hs::clamp(num - static_cast<float>(i), 0.0f, 1.0f);

      auto fragment_shader = [&](const Vector &, Fragment &f_val) {
        Color4 c = palette.get(static_cast<float>(i) / num);
        c.alpha *= opacity * params.alpha;
        f_val.color = c;
      };

      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader,
                            true);
    }
  }

  void draw_longitudes(Canvas &canvas, float num, float phase,
                       const Quaternion &q) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      ScratchScope _frag(scratch_arena_a);
      float theta = (static_cast<float>(i) / num) * PI_F;
      Vector normal(cosf(theta), 0.0f, -sinf(theta));

      Fragments m_points;
      m_points.bind(scratch_arena_a, 144);
      // Explicit basis construction to match JS texture alignment
      Vector v = normal;
      Vector w = Y_AXIS;
      Vector u = cross(v, w).normalized();
      Basis basis = {u, v, w};

      Plot::SphericalPolygon::sample(m_points, basis, 1.0f, W / 4);

      Fragments m_fragments;
      m_fragments.bind(scratch_arena_a, 144);
      for (size_t k = 0; k < m_points.size(); ++k) {
        Vector transformed = mobius_gen.transform(m_points[k].pos);
        Fragment f;
        f.pos = rotate(transformed, q).normalized();
        f.v0 = (m_points.size() > 1) ? (float)k / (m_points.size() - 1) : 0.0f;
        m_fragments.push_back(f);
      }

      float opacity = hs::clamp(num - static_cast<float>(i), 0.0f, 1.0f);

      auto fragment_shader = [&](const Vector &, Fragment &f_val) {
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

      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader,
                            true);
    }
  }

  GenerativePalette palette;
  GenerativePalette next_palette;
  MobiusWarpCircularTransformer<W, 1> mobius_gen;

  struct Params {
    float num_rings = 0.0f;
    float num_lines = 0.0f;
    float alpha = 0.2f;
  } params;

  Orientation<W> orientation;
  Timeline<W> timeline;

  Vector holeN;
  Vector holeS;

  Pipeline<W, H, Filter::World::HoleRef<W>, Filter::World::HoleRef<W>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(MobiusGrid)
