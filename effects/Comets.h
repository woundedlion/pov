/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "effects_engine.h"

// Configuration for Lissajous curves
struct LissajousConfig {
  float m1;
  float m2;
  float a;
  float domain;
};

template <int W, int H> class Comets : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115;

  struct Node {
    Orientation<W, 16> orientation;
    Animation::OrientationTrail<Orientation<W, 16>, TRAIL_LENGTH> trail;
    Vector v;

    Node() : v(Y_AXIS) {}
  };

  FLASHMEM Comets() : Effect(W, H), cur_function_idx(0) {}

  void init() override {
    node = static_cast<Node *>(
        persistent_arena.allocate(sizeof(Node), alignof(Node)));
    new (node) Node();

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.0f, 0.5f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Resolution", &params.resolution, 1.0f, 128.0f);
    registerParam("Debug BB", &params.debug_bb);

    update_path();
    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    timeline.add(0, Animation::Motion<W, 16>(node->orientation, path,
                                             (int)params.cycle_duration, true));
    timeline.add(0, Animation::PeriodicTimer(
                        2 * (int)params.cycle_duration,
                        [this](Canvas &c) {
                          cur_function_idx =
                              (cur_function_idx + 1) % functions.size();
                          update_path();
                          update_palette();
                        },
                        true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    node->trail.record(node->orientation);

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = palette.get(t);
        f.color.alpha *= quintic_kernel(t);
      };

      Vector v_local = rotate(node->v, q);
      Vector v_final = orientation.orient(v_local);
      Scan::Point::draw<W, H>(filters, canvas, v_final, params.thickness,
                              fragment_shader, params.debug_bb);
    });
  }

private:
  void update_path() {
    const LissajousConfig &config = functions[cur_function_idx];
    path.f = [&config](float t) {
      return lissajous(config.m1, config.m2, config.a, t * config.domain);
    };
  }

  void update_palette() {
    next_palette_ =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    timeline.add(0, Animation::ColorWipe(palette, next_palette_, 48, ease_mid));
  }

  FastNoiseLite noise;
  Timeline<W, 32> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ProceduralPath path;
  Orientation<W> orientation;
  GenerativePalette palette;
  std::array<LissajousConfig, 12> functions = {{{1.06f, 1.06f, 0, 5.909f},
                                                {6.06f, 1.0f, 0, 2 * PI_F},
                                                {6.02f, 4.01f, 0, 3.132f},
                                                {46.62f, 62.16f, 0, 0.404f},
                                                {46.26f, 69.39f, 0, 0.272f},
                                                {19.44f, 9.72f, 0, 0.646f},
                                                {8.51f, 17.01f, 0, 0.739f},
                                                {7.66f, 6.38f, 0, 4.924f},
                                                {8.75f, 5.0f, 0, 5.027f},
                                                {11.67f, 14.58f, 0, 2.154f},
                                                {11.67f, 8.75f, 0, 2.154f},
                                                {10.94f, 8.75f, 0, 2.872f}}};
  int cur_function_idx;
  Node *node = nullptr;
  GenerativePalette next_palette_;

  struct Params {
    float alpha = 1.0f;
    float thickness = 2.1f * 2 * PI_F / W;
    float cycle_duration = 80.0f;
    float resolution = 32.0f;
    bool debug_bb = false;
  } params;
};

#include "effect_registry.h"
REGISTER_EFFECT(Comets)
