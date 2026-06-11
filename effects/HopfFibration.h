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
      : Effect(W, H), trail_pipeline(Filter::Screen::AntiAlias<W, H>()) {
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
    // Ambient Y-axis spin is intentionally left ungated (keeps running while
    // paused, per the setAnimationsPaused contract). The flow/tumble Drivers
    // ARE gated so "Pause Animation" actually freezes the fiber motion.
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600,
                                           ease_mid, true));
    // Bound Drivers: each pulls its tuning slider (× per-unit rate) every step,
    // so the flow/tumble speeds stay live without a per-frame set_speed re-sync.
    timeline.add(0, Animation::Driver(flow_offset, &params.flow_speed, FLOW_RATE,
                                      false, &anims_paused_));
    timeline.add(0, Animation::Driver(tumble_angle_x, &params.tumble_speed,
                                      TUMBLE_X_RATE, false, &anims_paused_));
    timeline.add(0, Animation::Driver(tumble_angle_y, &params.tumble_speed,
                                      TUMBLE_Y_RATE, false, &anims_paused_));
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

  // Stereographic-projection guard: the denominator is (1 + eps) - q3, so the
  // projection pole (q3 == 1) maps to a finite 1/eps instead of dividing by
  // zero. normalized_or() then absorbs the degenerate direction.
  static constexpr float STEREO_POLE_EPSILON = 0.001f;

  // Per-unit driver speeds (multiplied by the corresponding tuning param via
  // the bound Drivers set up in init()).
  static constexpr float FLOW_RATE = 0.02f * 0.2f; // flow_offset / flow_speed
  static constexpr float TUMBLE_X_RATE = 0.003f;   // tumble_angle_x / tumble_speed
  static constexpr float TUMBLE_Y_RATE = 0.005f;   // tumble_angle_y / tumble_speed

  float flow_offset = 0.0f;
  float tumble_angle_x = 0.0f;
  float tumble_angle_y = 0.0f;
  Vector *fibers = nullptr;
  Animation::VectorTrail<TRAIL_LEN> *trails = nullptr;

  // Cached tumble rotation values (per-frame)
  float cx = 1.0f, sx = 0.0f, cy = 1.0f, sy = 0.0f, fold_base = 0.0f;

  Orientation<W> orientation;
  Timeline timeline;
  BakedPalette baked_sunset;

  // trail_pipeline applies AA only; trail points are oriented by hand in
  // render_trails before rasterizing, so no Orient filter is needed.
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

    // Stereographic S3 → R3. At a fiber pole (q0=q1=q2=0) the projected
    // direction is undefined; fall back to a stable axis instead of trapping.
    float factor = 1.0f / ((1.0f + STEREO_POLE_EPSILON) - q3);
    return normalized_or(Vector(q0 * factor, q1 * factor, q2 * factor),
                         Vector(1, 0, 0));
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
