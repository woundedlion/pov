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
class Dynamo : public Effect {
public:

  struct Node {
    Node() :
      x(0), y(0), v(0)
    {
    }

    int x;
    int y;
    int v;
  };

  Dynamo() :
    Effect(W),
    palettes({ GenerativePalette(
      GradientShape::VIGNETTE,
      HarmonyType::ANALOGOUS,
      BrightnessProfile::ASCENDING) }),
    palette_normal(Z_AXIS), // Z is UP in C++ (vs Y in JS)
    speed(2),
    gap(5),
    trail_length(8),

    filters(
      FilterWorldTrails<W, 10000>(trail_length),
      FilterReplicate<W>(3),
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>()
    )
  {
    persist_pixels = false;

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    timeline
      .add(0,
        RandomTimer(4, 64, [this](auto&) { reverse(); }, true))
      .add(0,
        RandomTimer(20, 64, [this](auto&) { color_wipe(); }, true))
      .add(0,
        RandomTimer(48, 160, [this](auto&) { rotate(); }, true));
  }

  bool show_bg() const override { return false; }

  void reverse() {
    speed *= -1;
  }

  void rotate() {
    timeline.add(0,
      Rotation<W>(
        orientation, random_vector(), PI_F, 40, ease_in_out_sin, false));
  }

  void color_wipe() {
    if (palettes.is_full()) {
      Serial.println("Palettes full, dropping color wipe!");
      return;
    }

    palettes.push_front(GenerativePalette(
      GradientShape::VIGNETTE,
      HarmonyType::ANALOGOUS,
      BrightnessProfile::ASCENDING));
    palette_boundaries.push_front(0);

    timeline.add(0,
      Transition(palette_boundaries.front(), PI_F, wipe_duration, ease_mid)
      .then([this]() {
        palette_boundaries.pop_back();
        palettes.pop_back();
        })
    );
  }

  Color4 color(const Vector& v, float t) {
    constexpr float blend_width = PI_F / 4;
    float a = angle_between(v, palette_normal);

    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      float boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return std::visit([t](auto& p) -> Color4 {
          using T = std::decay_t<decltype(p)>;
          if constexpr (std::is_same_v<T, std::monostate>) {
            return Color4(); // Should not happen if logic is correct
          } else {
            return p.get(t);
          }
        }, palettes[i]);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = std::clamp(blend_factor, 0.0f, 1.0f);

        Color4 c1 = std::visit([t](auto& p) -> Color4 {
          using T = std::decay_t<decltype(p)>;
          if constexpr (std::is_same_v<T, std::monostate>) { return Color4(); }
          else { return p.get(t); }
        }, palettes[i]);
        Color4 c2 = std::visit([t](auto& p) -> Color4 {
          using T = std::decay_t<decltype(p)>;
          if constexpr (std::is_same_v<T, std::monostate>) { return Color4(); }
          else { return p.get(t); }
        }, palettes[i + 1]);

        uint16_t fract = to_short(clamped_blend_factor);
        return Color4(
          c1.color.lerp16(c2.color, fract),
          lerp(c1.alpha, c2.alpha, clamped_blend_factor)
        );
      }

      auto next_boundary_lower_edge = (i + 1 < palette_boundaries.size()
        ? palette_boundaries[i + 1] - blend_width
        : 100.0f); // Infinity check

      if (a > upper_edge && a < next_boundary_lower_edge) {
        return std::visit([t](auto& p) -> Color4 {
          using T = std::decay_t<decltype(p)>;
          if constexpr (std::is_same_v<T, std::monostate>) { return Color4(); }
          else { return p.get(t); }
        }, palettes[i + 1]);
      }
    }

    return std::visit([t](auto& p) -> Color4 {
      using T = std::decay_t<decltype(p)>;
      if constexpr (std::is_same_v<T, std::monostate>) { return Color4(); }
      else { return p.get(t); }
    }, palettes[0]);
  }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    for (int i = std::abs(speed) - 1; i >= 0; --i) {
      pull(0);
      draw_nodes(canvas, static_cast<float>(i) / std::abs(speed));
    }
  }

private:

  void draw_nodes(Canvas& canvas, float age) {
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector<W>(nodes[i].x, nodes[i].y);
        Color4 c = color(from, 0);
        c.alpha *= 0.5f; 
        filters.plot(canvas, from, c.color, age, c.alpha);
      }
      else {
        auto from = pixel_to_vector<W>(nodes[i - 1].x, nodes[i - 1].y);
        auto to = pixel_to_vector<W>(nodes[i].x, nodes[i].y);
        Plot<W>::Line::draw(filters, canvas, from, to, [this](const Vector& v, float t) { 
            Color4 c = color(v, 0);
            c.alpha *= 0.5f;
            return c;
        }, 0, 1, false, true);
      }
    }
  }

  void pull(int leader) {
    nodes[leader].v = dir(speed);
    move(nodes[leader]);
    for (int i = leader - 1; i >= 0; --i) {
      drag(nodes[i + 1], nodes[i]);
    }
    for (size_t i = leader + 1; i < nodes.size(); ++i) {
      drag(nodes[i - 1], nodes[i]);
    }
  }

  void drag(Node& leader, Node& follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (shortest_distance(dest, leader.x, W) > gap) {
      follower.v = leader.v;
      while (shortest_distance(follower.x, leader.x, W) > gap) {
        move(follower);
      }
    }
    else {
      move(follower);
    }
  }

  void move(Node& node) {
    node.x = wrap(node.x + node.v, W);
  }

  int dir(int speed) const {
    return speed < 0 ? -1 : 1;
  }

  Timeline timeline;

  static constexpr size_t MAX_PALETTES = 16;
  static constexpr size_t NUM_NODES = H_VIRT;
  StaticCircularBuffer <PaletteVariant, MAX_PALETTES> palettes;
  StaticCircularBuffer <float, MAX_PALETTES - 1> palette_boundaries;

  Vector palette_normal;
  std::array<Node, NUM_NODES> nodes;
  int speed;
  int gap;
  uint32_t trail_length;
  Orientation orientation;

  int wipe_duration = 20;

  Pipeline<W,
    FilterWorldTrails<W, 10000>,
    FilterReplicate<W>,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;
};
