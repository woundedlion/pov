#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "effects_engine.h"

template <int W>
class RingSpin : public Effect {
public:

  struct Ring {
    Ring(const Vector& normal, Filter<W>& filters, const Palette& palette, uint8_t trail_length) :
      normal(normal),
      palette(palette),
      trails(trail_length)
    {    
      trails.chain(filters);
    }

    Vector normal;
    const Palette& palette;
    Orientation orientation;
    FilterDecay<W, 5000> trails;
  };
  
  RingSpin() :
    Effect(W),
    alpha(0.2),
    trail_length(5)
  {
    persist_pixels = false;
    rings.reserve(NUM_RINGS);
    for (int i = 0; i < NUM_RINGS; ++i) {
      spawn_ring(X_AXIS, *palettes[i]);
    }
  }

  bool show_bg() const { return false; }

  void spawn_ring(const Vector& normal, const Palette& palette) {
    auto ring_index = rings.size();
    rings.emplace_back(normal, filters, palette, trail_length);
    timeline.add(0,
      Sprite(
        [=](Canvas& canvas, double opacity) { draw_ring(canvas, opacity, ring_index); },
        -1,
        4, ease_mid,
        0, ease_mid
      ));
    timeline.add(0,
      RandomWalk<W>(rings[ring_index].orientation, rings[ring_index].normal));
  }

  void draw_ring(Canvas& canvas, double opacity, size_t ring_index) {
    auto& ring = rings[ring_index];
    size_t end = ring.orientation.length();
    size_t start = end == 1 ? 0 : 1;
    for (size_t i = start; i < end; ++i) {
      Dots dots;
      ::draw_ring<W>(dots, ring.orientation.orient(ring.normal, i), 1,
        [&](auto& v, auto t) { return vignette(ring.palette)(0); });
      plot_dots(dots, ring.trails, canvas, 
       (end - 1 - i) / end,
       alpha * opacity);
    }
    ring.orientation.collapse();

    ring.trails.trail(canvas,
      [&](double x, double y, double t) { return vignette(ring.palette)(1 - t); },
      alpha * opacity);
    
    ring.trails.decay();
    
  }

  void draw_frame() {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:

  static constexpr int NUM_RINGS = 6;
  std::vector<Ring> rings;
  double alpha;
  double trail_length;
  FilterAntiAlias<W> filters;
  std::array<const Palette*, 6> palettes = { &richSunset, &undersea, &mangoPeel, &lemonLime, &algae, &lateSunset };
  Timeline timeline;
};

template <int W>
class Dynamo : public Effect {
public:

  struct Node {
    Node() :
      x(0), y(0), v(0)
    {}

    int x;
    int y;
    int v;
  };

  Dynamo() :
    Effect(W),
    palettes({ GenerativePalette(GradientShape::VIGNETTE, HarmonyType::ANALOGOUS) }),
    palette_normal(X_AXIS),
    speed(2),
    gap(4),
    trail_length(5),
    filters(4),
    trails(trail_length),
    orient(orientation)
  {
    persist_pixels = false;

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    filters
      .chain(trails)
      .chain(orient)
      .chain(aa);

    timeline
      .add(0,
        RandomTimer(4, 64, [=](auto&) { reverse(); }, true))
      .add(0,
        RandomTimer(20, 64, [=](auto&) { color_wipe(); }, true))
      .add(0,
        RandomTimer(80, 150, [=](auto&) { rotate(); }, true));
  }

  bool show_bg() const { return false; }

  void reverse() {
    speed *= -1;
  }

  void rotate() {
    timeline.add(0,
      Rotation<W>(
        orientation, random_vector(), PI, 40, ease_in_out_sin, false));
  }

  void color_wipe() {
    palettes.push_back(GenerativePalette(GradientShape::VIGNETTE, HarmonyType::ANALOGOUS));
    palette_boundaries.push_back(0);
    timeline.add(0,
      Transition(palette_boundaries.back(), PI, 20, ease_mid)
      .then([this]() {
        palette_boundaries.pop();
        palettes.pop();
        })
    );
  }

  Pixel color(const Vector& v, double t) {
    constexpr double blend_width = PI / 8.0;
    double a = angle_between(v, palette_normal);

    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      double boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return std::visit([=](auto& p) { return p.get(t); },  palettes[i]);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = std::clamp(0.0, 1.0, blend_factor);
        auto c1 = std::visit([=](auto& p) { return p.get(t); }, palettes[i]);
        auto c2 = std::visit([=](auto& p) { return p.get(t); }, palettes[i + 1]);
        return c1.lerp16(c2, to_short(clamped_blend_factor));
      }

      auto next_boundary_lower_edge = (i + 1 < palette_boundaries.size()
        ? palette_boundaries[i + 1] - blend_width
        : PI);
      if (a > upper_edge && a < next_boundary_lower_edge) {
        return std::visit([=](auto& p) { return p.get(t); }, palettes[i + 1]);
      }
    }

    return std::visit([=](auto& p) { return p.get(t); }, palettes[0]);
  }

  void draw_frame() {
    Canvas canvas(*this);

    for (int i = std::abs(speed) - 1; i >= 0; --i) {
      pull(0);
      draw_nodes(canvas, static_cast<double>(i) / std::abs(speed));
    }
    trails.trail(
      canvas,
      [this](double x, double y, double t) { return color(pixel_to_vector<W>(x, y), t); },
      0.5
    );

    trails.decay();
    timeline.step(canvas);
  }

private:

  double node_y(const Node& node) const {
    return (node.y / (nodes.size() - 1)) * (H - 1);
  }

  void draw_nodes(Canvas& canvas, double age) {
    Dots dots;
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector<W>(nodes[i].x, node_y(nodes[i]));
        draw_vector(dots, from, [this](auto& v, auto t) { return color(v, 0); });
      } else {
        auto from = pixel_to_vector<W>(nodes[i - 1].x, node_y(nodes[i - 1]));
        auto to = pixel_to_vector<W>(nodes[i].x, node_y(nodes[i]));
        draw_line<W>(dots, from, to, [this](auto& v, auto t) { return color(v, 0); });
      }
    }
    plot_dots(dots, filters, canvas, age, 0.5);
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
    if (distance(dest, leader.x, W) > gap) {
      follower.v = leader.v;
      while (distance(follower.x, leader.x, W) > gap) {
        move(follower);
      }
    } else {
      move(follower);
    }
  }
  
  void move(Node& node) {
    node.x = wrap(node.x + node.v, W);
  }

  int dir(int speed) {
    return speed < 0 ? -1 : 1;
  }

  Timeline timeline;
 
  static constexpr size_t MAX_PALETTES = 16;
  static constexpr size_t NUM_NODES = H;
  StaticCircularBuffer <PaletteVariant, MAX_PALETTES> palettes;
  StaticCircularBuffer <double, MAX_PALETTES - 1> palette_boundaries; 

  Vector palette_normal;
  std::array<Node, NUM_NODES> nodes;
  int speed;
  int gap;
  uint32_t trail_length;
  Orientation orientation;
  
  FilterReplicate<W> filters;
  FilterDecay<W, 10000> trails;
  FilterOrient<W> orient;
  FilterAntiAlias<W> aa;
};

template <int W>
class RingShower2 : public Effect {
private:
  static constexpr size_t MAX_RINGS = 8;

  struct Ring {
    Vector normal;
    double speed;
    double radius;
    double duration;
    GenerativePalette palette;

    Ring() :
      normal(random_vector()),
      speed(0.0),
      radius(0.0),
      duration(0.0),
      palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS)
    {
    }

  };

public:
  RingShower2() :
    Effect(W)
  {
    persist_pixels = false;
    
    timeline.add(0,
      RandomTimer(1, 24,
        [this](auto&) { this->spawn_ring(); },
        true)
    );
    
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:

  void spawn_ring() {
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].duration <= 0) {
        Ring& ring = rings[i];
        ring.normal = random_vector();
        ring.duration = hs::rand_int(8, 72);
        ring.radius = 0;
        ring.palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS);

        timeline.add(0,
          Sprite(
            [this, i](Canvas& canvas, double opacity) { this->draw_ring(canvas, opacity, i); },
            ring.duration,
            4, ease_mid,
            0, ease_mid)
        );

        timeline.add(0,
          Transition(ring.radius, 2, ring.duration, ease_mid, false, false)
          .then([&ring]() { ring.duration = 0; })
        );

        return;
      }
    }
  }

  void draw_ring(Canvas& canvas, double opacity, size_t index) {
    Ring& ring = rings[index];
    Dots dots;
    ::draw_ring<W>(dots, orientation.orient(ring.normal), ring.radius,
      [&](auto& v, auto t) {
        Pixel color = ring.palette.get(t);
        return dim(color, opacity);
      },
      0);
    plot_dots<W>(dots, filters, canvas, 0, 1.0);
  }

  Ring rings[MAX_RINGS];
  FilterAntiAlias<W> filters;
  Orientation orientation;
  Timeline timeline;
};