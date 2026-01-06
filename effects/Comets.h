/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "../effects_engine.h"

template <int W>
class Comets : public Effect {
public:
  Comets() :
    Effect(W),
    palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING),
    next_palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING),
    trails(trail_length),
    filters(
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>())
  {
    persist_pixels = false;
    update_path();

    for (int i = 0; i < NUM_NODES; ++i) {
      spawn_node(i);
    }

    timeline.add(0,
      PeriodicTimer(2 * cycle_duration, [this](auto&) {
        increment_path();
        update_palette();
        }, true)
    );

    timeline.add(0, RandomWalk<W>(orientation, random_vector()));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    trails.render(canvas, filters,
      [this](const Vector& v, float t) {
        return palette.get(1.0f - t);
      }
    );
  }

private:

  struct Node {
    Orientation orientation;
    Vector v;
    const Path<W>& path;

    Node(const Path<W>& p) :
      v(Y_AXIS),
      path(p)
    {
    }
  };

  void increment_path() {
    cur_function = (cur_function + 1) % functions.size();
    update_path();
  }

  void update_path() {
    const auto& f = functions[cur_function];
    float max_speed = sqrtf(f.m1 * f.m1 + f.m2 * f.m2);
    float length = f.domain * max_speed;
    float samples = std::max(128.0f, ceilf(length * resolution));

    path.collapse();
    path.append_segment([f](float t) -> Vector {
      return lissajous(f.m1, f.m2, f.a, t);
      }, f.domain, samples, ease_mid);
  }

  void update_palette() {
    next_palette = GenerativePalette(
      GradientShape::STRAIGHT,
      HarmonyType::TRIADIC,
      BrightnessProfile::ASCENDING
    );

    timeline.add(0,
      ColorWipe(palette, next_palette, wipe_duration, ease_mid)
    );
  }

  void spawn_node(int i) {
    nodes.emplace_back(this->path);

    timeline.add(0,
      Sprite(
        [this, i](Canvas& canvas, float opacity) { draw_node(canvas, opacity, i); },
        -1,
        16, ease_mid,
        0, ease_mid
      )
    );

    timeline.add(i * spacing,
      Motion<W>(nodes[i].orientation, nodes[i].path, cycle_duration, true)
    );
  }

  void draw_node(Canvas& canvas, float opacity, int i) {
    Node& node = nodes[i];
    tween(node.orientation, [this, &canvas, opacity, &node](auto& q, auto t) {
      dots.clear();
      auto v = rotate(node.v, node.orientation.get()).normalize();
      ::draw_vector<W>(dots, v,
        [this](const auto& v, auto t) {
          return palette.get(t);
        });

      trails.record(dots, t, this->alpha * opacity);
      });
    node.orientation.collapse();
  }

  static constexpr int NUM_NODES = 1;
  float alpha = 0.5f;
  float resolution = 32.0f;
  size_t cycle_duration = 80;
  size_t trail_length = 80;
  size_t wipe_duration = 48;
  size_t spacing = 48;

  const std::array<LissajousParams, 12> functions = { {
    {1.06f, 1.06f, 0.0f, 5.909f},
    {6.06f, 1.0f, 0.0f, 2.0f * PI_F},
    {6.02f, 4.01f, 0.0f, 3.132f},
    {46.62f, 62.16f, 0.0f, 0.404f},
    {46.26f, 69.39f, 0.0f, 0.272f},
    {19.44f, 9.72f, 0.0f, 0.646f},
    {8.51f, 17.01f, 0.0f, 0.739f},
    {7.66f, 6.38f, 0.0f, 4.924f},
    {8.75f, 5.0f, 0.0f, 5.027f},
    {11.67f, 14.58f, 0.0f, 2.154f},
    {11.67f, 8.75f, 0.0f, 2.154f},
    {10.94f, 8.75f, 0.0f, 2.872f}
  } };

  Path<W> path;
  size_t cur_function = 0;
  Orientation orientation;
  GenerativePalette palette;
  GenerativePalette next_palette;

  DecayBuffer<W, 3000> trails;

  Pipeline<W,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;

  StaticCircularBuffer<Node, NUM_NODES> nodes;
  Timeline timeline;
  Dots dots;
};