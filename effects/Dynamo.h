/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

template <int W, int H> class Dynamo : public Effect {
public:
  struct Node {
    Node() : x(0), y(0), v(0) {}

    int x;
    int y;
    int v;
  };

  FLASHMEM Dynamo()
      : Effect(W, H),
        palettes(
            {GenerativePalette(GradientShape::VIGNETTE, HarmonyType::ANALOGOUS,
                               BrightnessProfile::ASCENDING)}),
        palette_normal(Z_AXIS),
        filters(Filter::World::Trails<W, TRAIL_CAPACITY>(
                    (uint32_t)params.trail_length),
                Filter::World::Replicate<W>(3),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  bool show_bg() const override { return false; }

  void init() override {
    static_cast<Filter::World::Trails<W, TRAIL_CAPACITY> &>(filters).init_storage(
        persistent_arena);

    registerParam("Speed", &params.speed, -10.0f, 10.0f);
    registerParam("Gap", &params.gap, 1.0f, 20.0f);
    registerParam("Trail Len", &params.trail_length, 1.0f, 100.0f);
    registerParam("Wipe Dur", &params.wipe_duration, 1.0f, 100.0f);

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    timeline
        .add(0, Animation::RandomTimer(
                    4, 64, [this](auto &) { reverse(); }, true))
        .add(0, Animation::RandomTimer(
                    20, 64, [this](auto &) { color_wipe(); }, true))
        .add(0, Animation::RandomTimer(
                    48, 160, [this](auto &) { rotate(); }, true));
  }

  void reverse() { params.speed *= -1; }

  void rotate() {
    timeline.add(0, Animation::Rotation<W>(orientation, random_vector(), PI_F,
                                           40, ease_in_out_sin, false));
  }

  void color_wipe() {
    if (palettes.is_full()) {
      hs::log("Palettes full, dropping color wipe!");
      return;
    }

    palettes.push_front(GenerativePalette(GradientShape::VIGNETTE,
                                          HarmonyType::ANALOGOUS,
                                          BrightnessProfile::ASCENDING));
    palette_boundaries.push_front(0);

    timeline.add(0, Animation::Transition(palette_boundaries.front(), PI_F,
                                          (int)params.wipe_duration, ease_mid)
                        .then([this]() {
                          palette_boundaries.pop_back();
                          palettes.pop_back();
                        }));
  }

  Color4 color(const Vector &v, float t) {
    constexpr float blend_width = PI_F / 4;
    // +inf sentinel for "no next boundary": `a` is an angle in [0, PI], so any
    // value well above PI makes the `a < next_boundary_lower_edge` test pass.
    constexpr float NO_NEXT_BOUNDARY = 100.0f;
    float a = angle_between(v, palette_normal);

    // palettes[i+1] below is always in range. Invariant: palettes.size() ==
    // palette_boundaries.size() + 1, held rigidly because color_wipe() pushes
    // one of each together (guarded by palettes.is_full()) and the .then()
    // callback pops one of each together. Independently, palette_boundaries'
    // capacity is MAX_PALETTES-1, so i+1 <= MAX_PALETTES-1 < palettes capacity
    // even if that invariant were ever violated. (Do NOT "fix" the boundary
    // buffer to MAX_PALETTES — that removes this second guard and *introduces*
    // the OOB.)
    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      float boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return palettes[i].get(t);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = hs::clamp(blend_factor, 0.0f, 1.0f);

        Color4 c1 = palettes[i].get(t);
        Color4 c2 = palettes[i + 1].get(t);

        uint16_t fract = float_to_pixel16(clamped_blend_factor);
        return Color4(c1.color.lerp16(c2.color, fract),
                      hs::lerp(c1.alpha, c2.alpha, clamped_blend_factor));
      }

      auto next_boundary_lower_edge =
          (i + 1 < palette_boundaries.size()
               ? palette_boundaries[i + 1] - blend_width
               : NO_NEXT_BOUNDARY);

      if (a > upper_edge && a < next_boundary_lower_edge) {
        return palettes[i + 1].get(t);
      }
    }

    return palettes[0].get(t);
  }

  void draw_frame() override {
    Canvas canvas(*this);

    // Live "Trail Len" slider: push the (clamped) param into the Trails filter
    // each frame so dragging it retunes the trail length. Clamp into the
    // filter's [1,255] domain (the slider tops out at 100) so set_lifetime's
    // trap only ever fires on a genuine authoring error.
    static_cast<Filter::World::Trails<W, TRAIL_CAPACITY> &>(filters)
        .set_lifetime(hs::clamp((int)params.trail_length, 1, 255));

    timeline.step(canvas);

    // Hoist the divisor so the i/steps division is provably guarded: the loop
    // body only runs when steps >= 1, so speed == 0 produces zero iterations
    // (no divide-by-zero) and falls straight through to the flush below.
    const int steps = std::abs((int)params.speed);
    for (int i = steps - 1; i >= 0; --i) {
      pull(0);
      draw_nodes(canvas, static_cast<float>(i) / steps);
    }

    filters.flush(
        canvas, [this](const Vector &v, float t) { return color(v, t); }, 1.0f);
  }

private:
  void draw_nodes(Canvas &canvas, float age) {
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        Color4 c = color(from, 0);
        c.alpha *= 0.5f;
        filters.plot(canvas, from, c.color, age, c.alpha);
      } else {
        auto from = pixel_to_vector<W, H>(nodes[i - 1].x, nodes[i - 1].y);
        auto to = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        auto fragment_shader = [this](const Vector &v, Fragment &f) {
          f.color = color(v, 0);
          f.color.alpha *= 0.5f;
        };
        Fragment f_from; f_from.pos = from;
        Fragment f_to;   f_to.pos = to;
        Plot::Line::draw<W, H>(filters, canvas, f_from, f_to,
                               fragment_shader);
      }
    }
  }

  void pull(int leader) {
    nodes[leader].v = dir((int)params.speed);
    move(nodes[leader]);
    for (int i = leader - 1; i >= 0; --i) {
      drag(nodes[i + 1], nodes[i]);
    }
    for (size_t i = leader + 1; i < nodes.size(); ++i) {
      drag(nodes[i - 1], nodes[i]);
    }
  }

  void drag(Node &leader, Node &follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (shortest_distance(dest, leader.x, W) > (int)params.gap) {
      follower.v = leader.v;
      while (shortest_distance(follower.x, leader.x, W) > (int)params.gap) {
        move(follower);
      }
    } else {
      move(follower);
    }
  }

  void move(Node &node) { node.x = wrap(node.x + node.v, W); }

  int dir(int speed) const { return speed < 0 ? -1 : 1; }

  Timeline timeline;

  static constexpr size_t MAX_PALETTES = 16;
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static constexpr size_t NUM_NODES = H_VIRT;
  // Compile-time Trails storage capacity (max trail points). The active trail
  // length is the runtime "Trail Len" slider, pushed into the Trails filter
  // each frame via set_lifetime() in draw_frame(); this bound is only the
  // fixed upper limit on buffered points.
  static constexpr int TRAIL_CAPACITY = 10000;
  StaticCircularBuffer<GenerativePalette, MAX_PALETTES> palettes;
  StaticCircularBuffer<float, MAX_PALETTES - 1> palette_boundaries;

  Vector palette_normal;
  std::array<Node, NUM_NODES> nodes;

  struct Params {
    float speed = 2.0f;
    float gap = 5.0f;
    float trail_length = 8.0f;
    float wipe_duration = 20.0f;
  } params;

  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Trails<W, TRAIL_CAPACITY>,
           Filter::World::Replicate<W>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Dynamo)
