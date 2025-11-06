#pragma once

#include "effects_infra.h"
#include <functional>
#include <memory>
#include <vector>
#include <map>

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
    FilterDecay<W, 3000> trails;
  };

  RingSpin() :
    Effect(W),
    alpha(0.2),
    trail_length(15)
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
    Dots dots;
    auto& ring = rings[ring_index];
    double s = ring.orientation.length();
    for (int i = 0; i < s; ++i) {
      ::draw_ring<W>(dots, ring.orientation.orient(ring.normal, i), 1,
        [&](auto& v, auto t) { return ring.palette.get(0); });
      plot_dots(dots, ring.trails, canvas, 
       (s - 1 - i) / s,
       alpha * opacity);
    }
    ring.trails.trail(canvas,
      [&](double x, double y, double t) { return vignette(ring.palette)(1 - t); },
      alpha * opacity);
    ring.trails.decay();
    ring.orientation.collapse();
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
    trail_length(10),
    filters(4),
    trails(trail_length),
    orient(orientation)
  {
    persist_pixels = false;

    for (int i = 0; i < NUM_NODES; ++i) {
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
    palettes.push_back(GenerativePalette(GradientShape::VIGNETTE));
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

    for (int i = 0; i < palette_boundaries.size(); ++i) {
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
    for (int i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector(nodes[i].x, node_y(nodes[i]));
        draw_vector(dots, from, [this](auto& v) { return color(v, 0); });
      } else {
        auto from = pixel_to_vector(nodes[i - 1].x, node_y(nodes[i - 1]));
        auto to = pixel_to_vector(nodes[i].x, node_y(nodes[i]));
        draw_line<W>(dots, from, to, [this](auto& v) { return color(v, 0); });
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
    for (int i = leader + 1; i < nodes.size(); ++i) {
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
  FilterDecay<W, 7000> trails;
  FilterOrient<W> orient;
  FilterAntiAlias<W> aa;
};

template <int W>
class RingShower2 : public Effect {
private:
  static constexpr size_t MAX_RINGS = 8; // Max simultaneous rings, fixed size for stability

  struct Ring {
    Vector normal;
    double speed;
    double radius;
    double duration;
    GenerativePalette palette;
    FilterDecay<W, 500> trails; // Use the V2 decay filter

    // The Sprite and Transition are not held here because they are managed by the Timeline,
    // but we need a Transition for the radius and a Sprite for the draw loop.

    Ring() :
      normal(random_vector()),
      speed(0.0),
      radius(0.0),
      duration(0.0),
      palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS), // Generate a random circular palette
      trails(10) // Small trail length
    {
      // Initialize the trail filter in the constructor
      trails.chain(filters);
    }
  };

public:
  RingShower2() :
    Effect(W)
  {
    persist_pixels = false; // We want a fresh canvas each frame

    // Initialize all rings and their persistent filters
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      rings[i].trails.chain(filters); // Chain the trail filter through the common filters
    }

    // T0: Start the continuous ring spawner
    timeline.add(0,
      RandomTimer(1, 24, // 1 to 24 frames (1/16s to 1.5s)
        [this](auto&) { this->spawn_ring(); },
        true) // Repeat forever
    );

    // Chain the final output filters once
    filters
      .chain(aa);
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  void spawn_ring() {
    // Find the first available (inactive) ring slot.
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].duration <= 0) {
        // Initialize new ring parameters
        Ring& ring = rings[i];
        ring.normal = random_vector();
        ring.duration = hs::rand_int(8, 72); // 8 to 72 frames
        ring.radius = 0.0;
        ring.palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS);

        // 1. Animation: Expand the radius from 0 to 2 (passing through 1)
        timeline.add(0,
          Transition(ring.radius, 2, ring.duration, ease_mid, false, false)
          .then([&ring]() { ring.duration = 0; }) // Mark ring slot as inactive when done
        );

        // 2. Sprite: The continuous drawing function
        timeline.add(0,
          Sprite(
            [this, i](Canvas& canvas, double opacity) { this->draw_ring(canvas, opacity, i); },
            ring.duration,
            4, ease_mid, // Fade in duration
            0, ease_mid)
        );
        return;
      }
    }
  }

  void draw_ring(Canvas& canvas, double opacity, size_t index) {
    Ring& ring = rings[index];
    Dots dots;

    // Use the Orientation for drawing, even if it's just the identity quaternion
    draw_ring<W>(dots, orientation.orient(ring.normal), ring.radius,
      [&](auto& v, auto t) {
        // Use the palette to get the color, dimming with overall opacity
        Pixel color = ring.palette.get(t);
        return dim(color, opacity);
      },
      0); // No phase shift

    // Plot the dots through the ring's private trail filter chain
    plot_dots<W>(dots, ring.trails, canvas, 0, 1.0);

    // Apply decay to the trails and trail it to the canvas
    // NOTE: The trail decay logic is handled implicitly by the Transition/Sprite duration.
    // We'll just decay the trails once per frame.
    ring.trails.decay();
    ring.trails.trail(canvas, vignette(ring.palette), 0.5); // Use the palette's vignette version for trails
  }

  // Static array of Ring structs to avoid dynamic allocation
  Ring rings[MAX_RINGS];

  // Common V2 infrastructure members
  FilterRaw<W> filters;
  FilterAntiAlias<W> aa;
  Orientation orientation; // Not used for rotation here, but required by some drawing fns
  Timeline timeline;
};