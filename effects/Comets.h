/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <vector>
#include <array>
#include "../effects_engine.h"

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
    Orientation<W> orientation;
    Animation::OrientationTrail<Orientation<W>, TRAIL_LENGTH> trail;
    Vector v;

    // Default constructor needed for StaticCircularBuffer
    Node() : v(Y_AXIS) {}
  };

  FLASHMEM Comets() : Effect(W, H), cur_function_idx(0) {

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.0f, 0.5f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Resolution", &params.resolution, 1.0f, 128.0f);

    // Initialize Lissajous functions
    functions = {{1.06f, 1.06f, 0, 5.909f},   {6.06f, 1.0f, 0, 2 * PI_F},
                 {6.02f, 4.01f, 0, 3.132f},   {46.62f, 62.16f, 0, 0.404f},
                 {46.26f, 69.39f, 0, 0.272f}, {19.44f, 9.72f, 0, 0.646f},
                 {8.51f, 17.01f, 0, 0.739f},  {7.66f, 6.38f, 0, 4.924f},
                 {8.75f, 5.0f, 0, 5.027f},    {11.67f, 14.58f, 0, 2.154f},
                 {11.67f, 8.75f, 0, 2.154f},  {10.94f, 8.75f, 0, 2.872f}};

    update_path();
    timeline.add(0, Animation::RandomWalk<W>(orientation, random_vector()));
    timeline.add(0, Animation::Motion<W>(node.orientation, path,
                                         (int)params.cycle_duration, true));
    timeline.add(0, Animation::PeriodicTimer(
                        2 * (int)params.cycle_duration,
                        [this](Canvas &c) {
                          cur_function_idx = static_cast<int>(
                              hs::rand_int(0, functions.size()));
                          update_path();
                          update_palette();
                        },
                        true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    node.trail.record(node.orientation);

    deep_tween(node.trail, [&](const Quaternion &q, float t) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = palette.get(t);
        f.color.alpha *= quintic_kernel(t);
      };

      Vector v_local = rotate(node.v, q);
      Vector v_final = orientation.orient(v_local);
      Scan::Point::draw<W, H>(filters, canvas, v_final, params.thickness,
                              fragment_shader);
    });
  }

private:
  void update_path() {
    const auto &config = functions[cur_function_idx];
    // Capture config by value to bake it into the lambda
    path.f = [=](float t) {
      // Apply easing and domain scaling
      float t_eased = ease_mid(t);
      return lissajous(config.m1, config.m2, config.a, t_eased * config.domain);
    };
  }

  void update_palette() {
    static GenerativePalette next_palette(GradientShape::STRAIGHT,
                                          HarmonyType::TRIADIC,
                                          BrightnessProfile::ASCENDING);
    next_palette =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    timeline.add(0, Animation::ColorWipe(palette, next_palette, 48, ease_mid));
  }
  Timeline<W, 32> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ProceduralPath path;
  Orientation<W> orientation;
  GenerativePalette palette;
  std::vector<LissajousConfig> functions;
  int cur_function_idx;
  Node node;

  struct Params {
    float alpha = 1.0f;
    float thickness = 2.1f * 2 * PI_F / W;
    float cycle_duration = 80.0f;
    float resolution = 32.0f;
  } params;
};