/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include <array>

template <int W, int H> class HopfFibration : public Effect {
public:
  static constexpr int MAX_TRAILS = 37500; // 300 KB at 8 bytes/item

  FLASHMEM HopfFibration()
      : Effect(W, H), filters(Filter::World::Trails<W, MAX_TRAILS>(40),
                              Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {}

  void init() override {
    registerParam("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    registerParam("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    registerParam("Folding", &params.folding, 0.0f, 2.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    static_cast<Filter::World::Trails<W, MAX_TRAILS> &>(filters).init_storage(
        persistent_arena);

    fibers = static_cast<Vector *>(persistent_arena.allocate(
        ACTUAL_FIBERS * sizeof(Vector), alignof(Vector)));
    prev_positions = static_cast<Vector *>(persistent_arena.allocate(
        ACTUAL_FIBERS * sizeof(Vector), alignof(Vector)));

    init_fibers();
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    advance_tumble();
    draw_fibers(canvas);
    first_frame_done = true;
    render_trails(canvas);
  }

  // Params
  struct Params {
    float flow_speed = 10.0f;
    float folding = 0.2f;
    float tumble_speed = 2.0f;
    float twist = 4.0f;
    float alpha = 0.4f;
  } params;

private:
  static constexpr int RINGS = 15;
  static constexpr int PER_RING = 14;
  static constexpr size_t ACTUAL_FIBERS = RINGS * PER_RING;

  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  bool first_frame_done = false;
  Vector *fibers = nullptr;
  Vector *prev_positions = nullptr;

  // Cached tumble rotation values (per-frame)
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;

  Orientation<W> orientation;
  Timeline<W> timeline;

  // Pipeline
  Pipeline<W, H, Filter::World::Trails<W, MAX_TRAILS>, Filter::World::Orient<W>,
           Filter::Screen::AntiAlias<W, H>>
      filters;

  void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float polar = PI_F * (i + 0.5f) / RINGS;
      for (int j = 0; j < PER_RING; ++j) {
        float azimuth = 2 * PI_F * j / PER_RING;
        fibers[idx++] = Vector(Spherical(azimuth, polar));
      }
    }
    first_frame_done = false;
  }

  void advance_tumble() {
    flow_offset += 0.02f * params.flow_speed * 0.2f;
    tumble_angle_x += 0.003f * params.tumble_speed;
    tumble_angle_y += 0.005f * params.tumble_speed;
    cx = cosf(tumble_angle_x);
    sx = sinf(tumble_angle_x);
    cy = cosf(tumble_angle_y);
    sy = sinf(tumble_angle_y);
    fold_base = sinf(tumble_angle_x * 0.5f) * 0.5f;
  }

  /// Project a base fiber through: folding → twist → S3 → tumble → stereo
  Vector hopf_project(size_t i) const {
    const Vector &base = fibers[i];
    Spherical sph(base);
    float theta = sph.phi;   // polar angle (co-latitude)
    float phi = sph.theta;   // azimuthal angle

    // Folding
    float eta = theta / 2.0f;
    eta += sinf(phi * 2.0f + tumble_angle_y + fold_base) * 0.1f *
           params.tumble_speed * params.folding;

    // Twist
    phi += eta * params.twist;

    // S3 point
    float phase = i * (PI_F / ACTUAL_FIBERS);
    float beta = flow_offset + phase;
    float q0 = cosf(eta) * cosf(phi + beta);
    float q1 = cosf(eta) * sinf(phi + beta);
    float q2 = sinf(eta) * cosf(beta);
    float q3 = sinf(eta) * sinf(beta);

    // Tumble (R_xw, R_yz)
    float q0_r = q0 * cx - q3 * sx;
    q3 = q0 * sx + q3 * cx;
    q0 = q0_r;
    float q1_r = q1 * cy - q2 * sy;
    q2 = q1 * sy + q2 * cy;
    q1 = q1_r;

    // Stereographic S3 → R3
    float factor = 1.0f / (1.001f - q3);
    return Vector(q0 * factor, q1 * factor, q2 * factor).normalize();
  }

  void draw_fibers(Canvas &canvas) {
    Color4 c = Palettes::richSunset.get(0.0f);
    c.alpha = params.alpha;

    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      Vector v = hopf_project(i);

      if (first_frame_done) {
        auto shader = [&](const Vector &, Fragment &f) { f.color = c; };
        Plot::Line::draw<W, H>(filters, canvas, Fragment(prev_positions[i]),
                               Fragment(v), shader);
      } else {
        filters.plot(canvas, v, c.color, 0, c.alpha);
      }
      prev_positions[i] = v;
    }
  }

  void render_trails(Canvas &canvas) {
    filters.flush(
        canvas,
        [](const Vector &, float t) {
          Color4 c = Palettes::richSunset.get(t);
          c.alpha *= (1.0f - t);
          return c;
        },
        1.0f);
  }
};

#include "../effect_registry.h"
REGISTER_EFFECT(HopfFibration)
