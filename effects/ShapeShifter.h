/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

template <int W, int H> class ShapeShifter : public Effect {
public:
  enum class ShapeType { PlanarPolygon, SphericalPolygon, Flower, Star };
  enum class RenderMode { Plot, Scan };

  FLASHMEM ShapeShifter()
      : Effect(W, H), current_shape(ShapeType::PlanarPolygon) {}

  bool show_bg() const override { return false; }

  void init() override {
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Count", &params.num_shapes, 1.0f, 128.0f);
    registerParam("Radius", &params.radius, 0.1f, 5.0f);
    registerParam("Sides", &params.sides, 3.0f, 12.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);
    registerParam("Debug BB", &params.debug_bb);
    // Twist is driven by the Mutation in build(); flag it so the GUI
    // auto-pauses the animation when the user grabs the slider.
    markAnimated("Twist");

    build();
  }

  void draw_frame() override {
    Canvas canvas(*this);

    // Advance to the next shape type every 48 frames
    frame_count_++;
    if (frame_count_ % 48 == 0) {
      int next = (static_cast<int>(current_shape) + 1) % 4;
      current_shape = static_cast<ShapeType>(next);
    }

    timeline.step(canvas);
  }

  // The whole effect runs on a *fixed* five-slot timeline regardless of Count:
  // the camera tumble, the twist Mutation, two shared orientation tumbles (one
  // per render mode), and a single draw event that walks the live Count each
  // frame. The per-ring Sprite+Rotation pair the effect used to spawn was pure
  // overhead — the Sprite never faded (constant opacity 1) and every Plot ring
  // shared one identical rotation, as did every Scan ring — so it cost two
  // global-timeline slots per shape and overflowed the 64-slot array at
  // Count 16 (2 + 4*16 = 66), silently dropping the last ring's animations.
  // Collapsing to two shared orientations is visually identical (the per-ring
  // rotations were already in lockstep) and lets Count scale to the raster
  // budget, not the slot budget.
  FLASHMEM void build() {
    timeline.clear();

    // Single camera tumble.
    timeline.add(0, Animation::RandomWalk<W>(camera, X_AXIS, noise, {},
                                             hs::rand_int(0, 65535)));

    // Twist Mutation: sin wave.
    timeline.add(0, Animation::Mutation(
                        params.twist,
                        [](float t) { return (PI_F / 4.0f) * sinf(t * PI_F); },
                        480, ease_mid, true, &anims_paused_));

    // Two shared tumbles: every Plot ring rides plot_orient_ (+X), every Scan
    // ring rides scan_orient_ (-X). These also drive the once-per-frame
    // orientation collapse (motion blur) before the draw event runs below.
    timeline.add(0, Animation::Rotation<W>(plot_orient_, X_AXIS, 2 * PI_F, 160,
                                           ease_mid, true,
                                           Animation::Space::Local));
    timeline.add(0, Animation::Rotation<W>(scan_orient_, -X_AXIS, 2 * PI_F, 160,
                                           ease_mid, true,
                                           Animation::Space::Local));

    // One draw event for every shape. Reads Count live, so a GUI change needs no
    // rebuild. Added last so it steps after the rotations have collapsed.
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float) { drawAll(canvas); }, -1,
                        0, ease_mid, 0, ease_mid));
  }

  void drawAll(Canvas &canvas) {
    int count = static_cast<int>(params.num_shapes);
    if (count < 1)
      count = 1;

    // Basis is identical across all rings of a mode (shared orientation + fixed
    // normal), so compute each once per frame instead of once per ring.
    Quaternion cam_q = camera.get();
    Basis plot_basis = make_basis(cam_q * plot_orient_.get(), X_AXIS);
    Basis scan_basis = make_basis(cam_q * scan_orient_.get(), -X_AXIS);

    for (int i = count - 1; i >= 0; --i) {
      float t = static_cast<float>(i) / (count > 1 ? count - 1 : 1);
      Color4 color = Palettes::richSunset.get(t);
      drawRing(canvas, plot_basis, RenderMode::Plot, t, color, i);
      drawRing(canvas, scan_basis, RenderMode::Scan, t, color, i);
    }
  }

  void drawRing(Canvas &canvas, const Basis &basis, RenderMode mode, float scale,
                const Color4 &color, int layer_index) {
    auto fragment_shader = [&](const Vector &p, Fragment &f) {
      Color4 c = color;
      c.alpha = c.alpha * this->params.alpha;
      f.color = c;
    };

    float phase = layer_index * this->params.twist;
    float r = this->params.radius * scale;
    int sides_int = (int)params.sides;
    if (mode == RenderMode::Plot) {
      dispatchPlot(canvas, basis, r, sides_int, fragment_shader, phase);
    } else {
      dispatchScan(canvas, basis, r, sides_int, fragment_shader, phase);
    }
  }

private:
  template <typename F>
  __attribute__((noinline)) void
  dispatchPlot(Canvas &canvas, const Basis &basis, float r, int sides_int,
               const F &fragment_shader, float phase) {
    switch (current_shape) {
    case ShapeType::Flower:
      Plot::Flower::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                               fragment_shader, {}, phase);
      break;
    case ShapeType::Star:
      Plot::Star::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                             fragment_shader, phase);
      break;
    case ShapeType::PlanarPolygon:
      Plot::PlanarPolygon::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                                      fragment_shader, phase);
      break;
    default:
      Plot::SphericalPolygon::draw<W, H>(plot_filters, canvas, basis, r,
                                         sides_int, fragment_shader, phase);
      break;
    }
  }

  template <typename F>
  __attribute__((noinline)) void
  dispatchScan(Canvas &canvas, const Basis &basis, float r, int sides_int,
               const F &fragment_shader, float phase) {
    switch (current_shape) {
    case ShapeType::Flower:
      Scan::Flower::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                               fragment_shader, phase, params.debug_bb);
      break;
    case ShapeType::Star:
      Scan::Star::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                             fragment_shader, phase, params.debug_bb);
      break;
    case ShapeType::PlanarPolygon:
      Scan::PlanarPolygon::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                                      fragment_shader, phase, params.debug_bb);
      break;
    default: // SphericalPolygon
      Scan::SphericalPolygon::draw<W, H>(scan_filters, canvas, basis, r,
                                         sides_int, fragment_shader, phase,
                                         params.debug_bb);
      break;
    }
  }

  FastNoiseLite noise;
  Timeline timeline;
  Orientation<> camera;
  // Shared tumbles — one per render mode. Declared after `timeline` so they
  // outlive it: ~Timeline clears the Rotations that point here on teardown.
  Orientation<> plot_orient_;
  Orientation<> scan_orient_;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> plot_filters;
  Pipeline<W, H> scan_filters;

  struct Params {
    float alpha = 0.5f;
    float num_shapes = 7.0f;
    float radius = 1.0f;
    float sides = 5.0f;
    float twist = 0.0f;
    bool debug_bb = false;
  } params;

  ShapeType current_shape;
  int frame_count_ = 0;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(ShapeShifter)
