/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/effects_engine.h"

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

    // Allocate baked palette LUT (refilled each frame)
    baked_palette.bake(persistent_arena, palette);

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.0f, 0.5f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Debug BB", &params.debug_bb);

    update_path();
    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    // Motion + cycle timer are infinite and added before the only finite
    // animation (the periodic ColorWipe, always appended after), so the
    // timeline never moves them — these add_get() handles stay valid for live
    // Cycle Dur updates.
    motion_ = timeline.add_get(
        0, Animation::Motion<W, 16>(node->orientation, path,
                                    (int)params.cycle_duration, true));
    cycle_timer_ = timeline.add_get(
        0, Animation::PeriodicTimer(
               2 * (int)params.cycle_duration,
               [this](Canvas &c) {
                 cur_function_idx = (cur_function_idx + 1) % functions.size();
                 update_path();
                 update_palette();
               },
               true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Live-apply the Cycle Dur slider to the motion duration + cycle timer.
    int cd = (int)params.cycle_duration;
    if (cd != last_cycle_dur_) {
      last_cycle_dur_ = cd;
      if (motion_)
        motion_->set_duration(cd);
      if (cycle_timer_)
        cycle_timer_->set_period(2 * cd);
    }

    // Re-bake animated palette each frame (256 lookups vs ~1700)
    baked_palette.rebake(palette);

    node->trail.record(node->orientation);

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = baked_palette.get(t);
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
    LissajousParams config = functions[cur_function_idx];
    path.f = [config](float t) {
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
  Pipeline<W, H> filters;
  ProceduralPath path;
  Orientation<W> orientation;
  GenerativePalette palette;
  BakedPalette baked_palette;
  std::array<LissajousParams, 12> functions = {{{1.06f, 1.06f, 0, 5.909f},
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
  Animation::Motion<W, 16> *motion_ = nullptr;
  Animation::PeriodicTimer *cycle_timer_ = nullptr;
  int last_cycle_dur_ = -1;

  struct Params {
    float alpha = 1.0f;
    float thickness = 2.1f * 2 * PI_F / W;
    float cycle_duration = 80.0f;
    bool debug_bb = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Comets)
