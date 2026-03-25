/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include <array>

template <int W, int H> class HopfFibration : public Effect {
public:
  static constexpr int TRAIL_LEN = 40;

  FLASHMEM HopfFibration()
      : Effect(W, H), fiber_pipeline(Filter::World::Orient<W>(orientation),
                                     Filter::Screen::AntiAlias<W, H>()),
        trail_pipeline(Filter::Screen::AntiAlias<W, H>()) {
    persist_pixels = false;
  }

  void init() override {
    registerParam("Flow Spd", &params.flow_speed, 0.0f, 20.0f);
    registerParam("Tumble Spd", &params.tumble_speed, 0.0f, 10.0f);
    registerParam("Folding", &params.folding, 0.0f, 2.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Bake palette for trail coloring
    baked_sunset.bake(persistent_arena, Palettes::richSunset);

    fibers = static_cast<Vector *>(persistent_arena.allocate(
        ACTUAL_FIBERS * sizeof(Vector), alignof(Vector)));

    // Allocate VectorTrails for each fiber
    trails = static_cast<Animation::VectorTrail<TRAIL_LEN> *>(
        persistent_arena.allocate(ACTUAL_FIBERS *
                                      sizeof(Animation::VectorTrail<TRAIL_LEN>),
                                  alignof(Animation::VectorTrail<TRAIL_LEN>)));
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      new (&trails[i]) Animation::VectorTrail<TRAIL_LEN>();
    }

    init_fibers();
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_mid, true));
    flow_driver_ = timeline.add_get(
        0, Animation::Driver(flow_offset, 0.02f * params.flow_speed * 0.2f,
                             false));
    tumble_x_driver_ = timeline.add_get(
        0,
        Animation::Driver(tumble_angle_x, 0.003f * params.tumble_speed, false));
    tumble_y_driver_ = timeline.add_get(
        0,
        Animation::Driver(tumble_angle_y, 0.005f * params.tumble_speed, false));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    advance_tumble();

    // Project fibers and record trail positions (oriented world space)
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      Vector v = hopf_project(i);
      // Store in world space (unoriented) — orient at render time
      trails[i].record(v);
    }

    // Render trails as polylines (already oriented, skip Orient filter)
    render_trails(canvas);
  }

  // Params
  struct Params {
    float flow_speed = 10.0f;
    float folding = 0.2f;
    float tumble_speed = 2.0f;
    float twist = 4.0f;
    float alpha = 1.0f;
  } params;

private:
  static constexpr int RINGS = 15;
  static constexpr int PER_RING = 14;
  static constexpr size_t ACTUAL_FIBERS = RINGS * PER_RING;

  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  Vector *fibers = nullptr;
  Animation::VectorTrail<TRAIL_LEN> *trails = nullptr;
  Animation::Driver *flow_driver_ = nullptr;
  Animation::Driver *tumble_x_driver_ = nullptr;
  Animation::Driver *tumble_y_driver_ = nullptr;

  // Cached tumble rotation values (per-frame)
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;

  Orientation<W> orientation;
  Timeline<W> timeline;
  BakedPalette baked_sunset;

  // Two pipelines: fibers go through Orient+AA, trails go through AA only
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      fiber_pipeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> trail_pipeline;

  void init_fibers() {
    int idx = 0;
    for (int i = 0; i < RINGS; ++i) {
      float polar = PI_F * (i + 0.5f) / RINGS;
      for (int j = 0; j < PER_RING; ++j) {
        float azimuth = 2 * PI_F * j / PER_RING;
        fibers[idx++] = Vector(Spherical(azimuth, polar));
      }
    }
  }

  void advance_tumble() {
    flow_driver_->set_speed(0.02f * params.flow_speed * 0.2f);
    tumble_x_driver_->set_speed(0.003f * params.tumble_speed);
    tumble_y_driver_->set_speed(0.005f * params.tumble_speed);
    cx = fast_cosf(tumble_angle_x);
    sx = fast_sinf(tumble_angle_x);
    cy = fast_cosf(tumble_angle_y);
    sy = fast_sinf(tumble_angle_y);
    fold_base = fast_sinf(tumble_angle_x * 0.5f) * 0.5f;
  }

  /// Project a base fiber through: folding → twist → S3 → tumble → stereo
  Vector hopf_project(size_t i) const {
    const Vector &base = fibers[i];
    Spherical sph(base);
    float theta = sph.phi; // polar angle (co-latitude)
    float phi = sph.theta; // azimuthal angle

    // Folding
    float eta = theta / 2.0f;
    eta += fast_sinf(phi * 2.0f + tumble_angle_y + fold_base) * 0.1f *
           params.tumble_speed * params.folding;

    // Twist
    phi += eta * params.twist;

    // S3 point
    float phase = i * (PI_F / ACTUAL_FIBERS);
    float beta = flow_offset + phase;
    float cos_eta = fast_cosf(eta), sin_eta = fast_sinf(eta);
    float q0 = cos_eta * fast_cosf(phi + beta);
    float q1 = cos_eta * fast_sinf(phi + beta);
    float q2 = sin_eta * fast_cosf(beta);
    float q3 = sin_eta * fast_sinf(beta);

    // Tumble (R_xw, R_yz)
    float q0_r = q0 * cx - q3 * sx;
    q3 = q0 * sx + q3 * cx;
    q0 = q0_r;
    float q1_r = q1 * cy - q2 * sy;
    q2 = q1 * sy + q2 * cy;
    q1 = q1_r;

    // Stereographic S3 → R3
    float factor = 1.0f / (1.001f - q3);
    return Vector(q0 * factor, q1 * factor, q2 * factor).normalized();
  }

  void render_trails(Canvas &canvas) {
    for (size_t i = 0; i < ACTUAL_FIBERS; ++i) {
      const auto &trail = trails[i];
      size_t len = trail.length();
      if (len < 2)
        continue;

      ScratchScope _sc(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, len);

      for (size_t j = 0; j < len; ++j) {
        Fragment f;
        // Orient from world space to current view at render time
        f.pos = orientation.orient(trail.get(j));
        f.v0 = (len > 1) ? static_cast<float>(j) / (len - 1) : 0.0f;
        f.age = 0;
        points.push_back(f);
      }

      auto shader = [this](const Vector &, Fragment &f) {
        float t = f.v0; // 0 = oldest, 1 = newest
        Color4 c = baked_sunset.get(1.0f - t);
        c.alpha *= t * params.alpha; // fade out tail
        f.color = c;
      };

      // Rasterize polyline through trail_pipeline (AA only, no Orient)
      Plot::rasterize<W, H>(trail_pipeline, canvas, points, shader);
    }
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(HopfFibration)
