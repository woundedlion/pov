/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once
#include "../effects_engine.h"

template <int W, int H>
class Thrusters : public Effect {
public:
  Thrusters() :
    Effect(W, H),
    palette({ 0.5f, 0.5f, 0.5f }, { 0.5f, 0.5f, 0.5f }, { 0.3f, 0.3f, 0.3f }, { 0.0f, 0.2f, 0.6f }),
    filters(Filter::Screen::AntiAlias<W, H>()),
    ring_vec(0.5f, 0.5f, 0.5f),
    amplitude(0),
    radius(1.0f),
    warp_phase(0),
    t_global(0),
    warp_anim(amplitude, [](float) { return 0.0f; }, 0, ease_mid)
  {
    persist_pixels = false;
    ring_vec = ring_vec.normalize();

    timeline.add(0, Animation::Sprite(
      [this](Canvas& c, float opacity) { draw_ring(c, opacity); },
      -1, 16, ease_in_sin, 16, ease_out_sin
    ));

    timeline.add(0, Animation::RandomTimer(16, 48,
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
    Orientation<W> orientation;
    Vector point;
    float radius;
    Animation::Transition motion;

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

    void reset(const Orientation<W>& o, const Vector& p) {
      orientation = o;
      point = p;
      radius = 0.0f;
      motion = Animation::Transition(radius, 0.3f, 8, ease_mid);
    }
  };

  StaticCircularBuffer<ThrusterContext, 16> thrusters;

  void on_fire_thruster() {
    warp_phase = hs::rand_f() * 2 * PI_F;

    auto r_fn = [this](float t) { return ring_fn(t); };
    Basis basis = make_basis(Quaternion(), ring_vec);
    Vector thrust_point = Plot::DistortedRing::fn_point(r_fn, basis, 1.0f, warp_phase);
    Vector thrust_opp = Plot::DistortedRing::fn_point(r_fn, basis, 1.0f, warp_phase + PI_F);

    // warp
    warp_anim = Animation::Mutation(
      amplitude,
      [](float t) { return 0.7f * expf(-2.0f * t); },
      32, ease_mid
    );

    // spin
    Vector thrust_axis = cross(orientation.orient(thrust_point), orientation.orient(ring_vec)).normalize();
    timeline.add(0, Animation::Rotation<W>(orientation, thrust_axis, 2 * PI_F, 8 * 16, ease_out_expo));

    // spawn
    spawn_thruster(thrust_point);
    spawn_thruster(thrust_opp);
  }

  void spawn_thruster(const Vector& point) {
    if (thrusters.is_full()) thrusters.pop();
    thrusters.push_back(ThrusterContext());
    ThrusterContext& ctx = thrusters.back();
    ctx.reset(orientation, point);

    timeline.add(0, Animation::Sprite(
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
    Basis basis = make_basis(ctx.orientation.get(), ctx.point);
    auto fragment_shader = [](const Vector&, Fragment& f) { 
        f.color = Color4(CRGB(255, 255, 255));
    };
    Plot::Ring::draw<W, H>(filters, c, basis, ctx.radius, fragment_shader);
  }

  void draw_ring(Canvas& c, float opacity) {
    Basis basis = make_basis(orientation.get(), ring_vec);
    
    auto fragment_shader = [this](const Vector& v, Fragment& f) { // Color function
        Vector axis = orientation.orient(X_AXIS);
        float angle = angle_between(axis, orientation.orient(v));
        f.color = palette.get(angle / PI_F);
    };

    Plot::DistortedRing::draw<W, H>(filters, c, basis, radius,
      [this](float t) { return ring_fn(t); }, // Shift function
      fragment_shader
    );
  }

  ProceduralPalette palette;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  Vector ring_vec;
  float amplitude;
  float radius;
  float warp_phase;
  int t_global;

  Timeline<W> timeline;
  Animation::Mutation warp_anim;
  Orientation<W> orientation;
};