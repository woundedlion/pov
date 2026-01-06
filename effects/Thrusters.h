/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include "../effects_engine.h"

template <int W>
class Thrusters : public Effect {
public:
  Thrusters() :
    Effect(W),
    palette({ 0.5f, 0.5f, 0.5f }, { 0.5f, 0.5f, 0.5f }, { 0.3f, 0.3f, 0.3f }, { 0.0f, 0.2f, 0.6f }),
    filters(FilterAntiAlias<W>()),
    ring_vec(0.5f, 0.5f, 0.5f),
    amplitude(0),
    radius(1.0f),
    warp_phase(0),
    t_global(0),
    warp_anim(amplitude, [](float) { return 0.0f; }, 0, ease_mid)
  {
    persist_pixels = false;
    ring_vec = ring_vec.normalize();

    timeline.add(0, Sprite(
      [this](Canvas& c, float opacity) { draw_ring(c, opacity); },
      -1, 16, ease_in_sin, 16, ease_out_sin
    ));

    timeline.add(0, RandomTimer(16, 48,
      [this](auto&) { on_fire_thruster(); }, true)
    );
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    if (!warp_anim.done()) {
      Canvas dummy(*this); // Mutation step doesn't need canvas but interface requires it
      warp_anim.step(dummy);
    }
    t_global++;
  }

private:

  struct ThrusterContext {
    Orientation orientation;
    Vector point;
    float radius;
    Transition motion;

    ThrusterContext() : radius(0), motion(radius, 0.3f, 8, ease_mid) {}

    ThrusterContext(const ThrusterContext& other)
      : orientation(other.orientation), point(other.point),
      radius(other.radius), motion(other.motion)
    {
      motion.rebind_mutant(radius);
    }

    ThrusterContext& operator=(const ThrusterContext& other) {
      if (this != &other) {
        orientation = other.orientation;
        point = other.point;
        radius = other.radius;
        motion = other.motion;
        motion.rebind_mutant(radius);
      }
      return *this;
    }

    void reset(const Orientation& o, const Vector& p) {
      orientation = o;
      point = p;
      radius = 0.0f;
      motion = Transition(radius, 0.3f, 8, ease_mid);
    }
  };

  StaticCircularBuffer<ThrusterContext, 16> thrusters;

  void on_fire_thruster() {
    warp_phase = hs::rand_f() * 2 * PI_F;

    auto r_fn = [this](float t) { return ring_fn(t); };
    Vector thrust_point = fn_point(r_fn, ring_vec, 1.0f, warp_phase);
    Vector thrust_opp = fn_point(r_fn, ring_vec, 1.0f, warp_phase + PI_F);

    // warp
    warp_anim = Mutation(
      amplitude,
      [](float t) { return 0.7f * expf(-2.0f * t); },
      32, ease_mid
    );

    // spin
    Vector thrust_axis = cross(orientation.orient(thrust_point), orientation.orient(ring_vec)).normalize();
    timeline.add(0, Rotation<MAX_W>(orientation, thrust_axis, 2 * PI_F, 8 * 16, ease_out_expo));

    // spawn
    spawn_thruster(thrust_point);
    spawn_thruster(thrust_opp);
  }

  void spawn_thruster(const Vector& point) {
    if (thrusters.is_full()) thrusters.pop();
    thrusters.push_back(ThrusterContext());
    ThrusterContext& ctx = thrusters.back();
    ctx.reset(orientation, point);

    timeline.add(0, Sprite(
      [this, &ctx](Canvas& c, float opacity) {
        ctx.motion.step(c);
        draw_thruster(c, ctx, opacity);
      },
      16, 0, ease_mid, 16, ease_out_expo
    ));
  }

  float ring_fn(float t) {
    return sin_wave(-1, 1, 2, warp_phase)(t)
      * sin_wave(-1, 1, 3, 0)(static_cast<float>(t_global % 32) / 32.0f)
      * amplitude;
  }

  void draw_thruster(Canvas& c, const ThrusterContext& ctx, float opacity) {
    Dots dots;
    ::draw_ring<W>(dots, ctx.orientation.get(), ctx.point, ctx.radius,
      [](const Vector&, float) { return Color4(CRGB::White); });
    plot_dots<W>(dots, filters, c, 0, opacity * 0.2f);
  }

  void draw_ring(Canvas& c, float opacity) {
    Dots dots;
    ::draw_fn<W>(dots, orientation.get(), ring_vec, radius,
      [this](float t) { return ring_fn(t); }, // Shift function
      [this](const Vector& v, float t) { // Color function
        Vector z_axis = orientation.orient(X_AXIS);
        float angle = angle_between(z_axis, orientation.orient(v));
        return palette.get(angle / PI_F);
      }
    );
    plot_dots<W>(dots, filters, c, 0, opacity * 0.2f);
  }

  ProceduralPalette palette;
  Pipeline<W, FilterAntiAlias<W>> filters;

  Vector ring_vec;
  float amplitude;
  float radius;
  float warp_phase;
  int t_global;

  Timeline timeline;
  Mutation warp_anim;
  Orientation orientation;
};