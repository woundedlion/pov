/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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

template <int W>
class Comets : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 80;
  static constexpr int MAX_NODES = 4; // Capacity for nodes



  struct Node {
    Orientation orientation;
    OrientationTrail<TRAIL_LENGTH> trail;
    Vector v;

    // Default constructor needed for StaticCircularBuffer
    Node() : v(Y_AXIS) {}
  };

  Comets() :
    Effect(W),
    palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::DESCENDING),
    cur_function_idx(0),
    num_nodes(1),
    spacing(48),
    resolution(32),
    cycle_duration(80),
    alpha(0.5f),
    thickness(2.1f * 2 * PI_F / W)
  {
    persist_pixels = false;

    // Initialize Lissajous functions
    functions = {
        {1.06f, 1.06f, 0, 5.909f},
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
        {10.94f, 8.75f, 0, 2.872f}
    };

    update_path();

    // Spawn nodes
    for (int i = 0; i < num_nodes; ++i) {
      spawn_node();
    }

    timeline.add(0, PeriodicTimer(2 * cycle_duration, [this](Canvas& c) {
      cur_function_idx = static_cast<int>(hs::rand_int(0, functions.size()));
      update_path();
      update_palette();
      }, true));

    timeline.add(0, RandomWalk<W>(orientation, random_vector()));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    std::vector<Dot> render_points;

    for (auto& node : nodes) {
      node.trail.record(node.orientation);

      deep_tween(node.trail, [&](const Quaternion& q, float t_trail) {
        Color4 c = palette.get(t_trail);
        c.alpha = c.alpha * alpha * quintic_kernel(1.0f - t_trail);
        Vector v_local = rotate(node.v, q);
        Vector v_final = orientation.orient(v_local);
        render_points.emplace_back(Dot(v_final, c));
      });
    }

    for(const auto& dot : render_points) {
        Scan<W>::Point::draw(filters, canvas, dot.position, thickness, [&](const Vector&, float){ return dot.color; });
    }
  }

private:

  void update_path() {
    const auto& config = functions[cur_function_idx];
    float max_speed = sqrtf(config.m1 * config.m1 + config.m2 * config.m2);
    float length = config.domain * max_speed;
    int samples = std::max(128, static_cast<int>(ceilf(length * resolution)));

    Path<W> new_path;
    new_path.append_segment(
      [&](float t) { return lissajous(config.m1, config.m2, config.a, t); },
      config.domain,
      static_cast<float>(samples),
      ease_mid
    );
    path = new_path;
  }

  void update_palette() {
    static GenerativePalette next_palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);
    next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);

    timeline.add(0, ColorWipe(palette, next_palette, 48, ease_mid));
  }

  void spawn_node() {
    if (nodes.is_full()) return;

    int i = static_cast<int>(nodes.size());
    nodes.push_back(Node());
    Node& node = nodes.back();

    timeline.add(i * spacing, Motion<W>(node.orientation, path, cycle_duration, true));
  }

  Timeline timeline;
  Pipeline<W, FilterAntiAlias<W>> filters;
  Path<W> path;
  Orientation orientation;
  GenerativePalette palette;
  std::vector<LissajousConfig> functions;
  int cur_function_idx;
  StaticCircularBuffer<Node, MAX_NODES> nodes;

  int num_nodes;
  int spacing;
  int resolution;
  int cycle_duration;
  float alpha;
  float thickness;
};