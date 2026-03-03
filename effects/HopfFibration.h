/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include <array>
#include <cmath>

template <int W, int H> class HopfFibration : public Effect {
public:
  static constexpr int MAX_TRAILS = 50000;

  FLASHMEM HopfFibration()
      : Effect(W, H), filters(Filter::World::Trails<W, MAX_TRAILS>(40),
                              Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {
    registerParam("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    registerParam("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    registerParam("Folding", &params.folding, 0.0f, 2.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    init_fibers();
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);

    timeline.step(canvas);

    // Update Params
    flow_offset += 0.02f * params.flow_speed * 0.2f;
    tumble_angle_x += 0.003f * params.tumble_speed;
    tumble_angle_y += 0.005f * params.tumble_speed;

    float cx = cosf(tumble_angle_x);
    float sx = sinf(tumble_angle_x);
    float cy = cosf(tumble_angle_y);
    float sy = sinf(tumble_angle_y);
    float fold_base = sinf(tumble_angle_x * 0.5f) * 0.5f;

    for (size_t i = 0; i < fibers.size(); ++i) {
      const Vector &base = fibers[i];

      // Hopf Fiber Params (S2 base)
      float theta = acosf(base.j); // y is up
      float phi = atan2f(base.k, base.i);

      // Folding
      float eta = theta / 2.0f;
      float folding_val = sinf(phi * 2.0f + tumble_angle_y + fold_base) * 0.1f *
                          params.tumble_speed * params.folding;
      eta += folding_val;

      // Twist
      phi += eta * params.twist;

      // Dot Generation
      float phase = i * (PI_F / fibers.size());
      float beta = flow_offset + phase;

      // 1. Construct point on S3
      // q = [cos(eta)cos(phi+beta), cos(eta)sin(phi+beta), sin(eta)cos(beta),
      // sin(eta)sin(beta)]
      float q0 = cosf(eta) * cosf(phi + beta);
      float q1 = cosf(eta) * sinf(phi + beta);
      float q2 = sinf(eta) * cosf(beta);
      float q3 = sinf(eta) * sinf(beta);

      // 2. Apply Tumble (R_xw, R_yz)
      // R_xw
      float q0_r = q0 * cx - q3 * sx;
      float q3_r = q0 * sx + q3 * cx;
      q0 = q0_r;
      q3 = q3_r;

      // R_yz
      float q1_r = q1 * cy - q2 * sy;
      float q2_r = q1 * sy + q2 * cy;
      q1 = q1_r;
      q2 = q2_r;

      // 3. Stereographic Projection S3 -> R3
      float div = 1.001f - q3;
      float factor = 1.0f / div;
      Vector v(q0 * factor, q1 * factor, q2 * factor);
      v = v.normalize();

      Color4 c = Palettes::richSunset.get(0.0f);
      c.alpha = params.alpha; // parameter alpha

      if (first_frame_done) {
        const Vector &prev = prev_positions[i];
        // Draw line segment
        auto fragment_shader = [&](const Vector &, Fragment &f) {
          f.color = c;
        };
        Plot::Line::draw<W, H>(filters, canvas, Fragment(prev), Fragment(v),
                               fragment_shader);

        prev_positions[i] = v;
      } else {
        // First frame
        filters.plot(canvas, v, c.color, 0, c.alpha);
        prev_positions[i] = v;
      }
    }

    first_frame_done = true;

    // Render Trails
    filters.flush(
        canvas,
        [](const Vector &v, float t) {
          Color4 c = Palettes::richSunset.get(t);
          c.alpha *= (1.0f - t);
          return c;
        },
        1.0f);
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
  static constexpr int NUM_FIBERS = 200;
  static constexpr int RINGS = 15;
  static constexpr int PER_RING = 14;
  static constexpr size_t ACTUAL_FIBERS = RINGS * PER_RING;

  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  bool first_frame_done = false;

  std::array<Vector, ACTUAL_FIBERS> fibers;
  std::array<Vector, ACTUAL_FIBERS> prev_positions;

  Orientation<W> orientation;
  Timeline<W> timeline;

  // Pipeline
  Pipeline<W, H, Filter::World::Trails<W, MAX_TRAILS>, Filter::World::Orient<W>,
           Filter::Screen::AntiAlias<W, H>>
      filters;

  void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float theta = PI_F * (i + 0.5f) / RINGS;
      float y = cosf(theta);
      float r = sinf(theta);

      for (int j = 0; j < PER_RING; ++j) {
        float phi = 2 * PI_F * j / PER_RING; // 0..2PI

        // Y-up
        float x = r * cosf(phi);
        float z = r * sinf(phi);

        fibers[idx++] = Vector(x, y, z);
      }
    }
    first_frame_done = false;
  }
};
