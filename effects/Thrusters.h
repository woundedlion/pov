/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/effects_engine.h"

template <int W, int H> class Thrusters : public Effect {
public:
  FLASHMEM Thrusters()
      // Inigo Quilez cosine palette: color(t) = bias + amp*cos(2π(freq*t + phase)).
      : Effect(W, H), palette(/*bias*/ {0.5f, 0.5f, 0.5f},
                              /*amp*/ {0.5f, 0.5f, 0.5f},
                              /*freq*/ {0.3f, 0.3f, 0.3f},
                              /*phase*/ {0.0f, 0.2f, 0.6f}),
        filters(Filter::Screen::AntiAlias<W, H>()), ring_vec(0.5f, 0.5f, 0.5f),
        amplitude(0), warp_phase(0), t_global(0),
        warp_anim(amplitude, [](float) { return 0.0f; }, 0, ease_mid) {}

  void init() override {
    registerParam("Radius", &params.radius, 0.1f, 2.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    ring_vec.normalize();

    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) { draw_ring(c, opacity); }, -1,
               16, ease_in_sin, 16, ease_out_sin));

    timeline.add(0, Animation::RandomTimer(
                        16, 48, [this](auto &) { on_fire_thruster(); }, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    // Stepped manually rather than via the timeline: each fire reassigns
    // warp_anim to restart the decay, and fires routinely overlap its life, so
    // a single restartable member is correct here where timeline.add() would
    // stack concurrent mutations fighting over `amplitude`.
    if (!warp_anim.done()) {
      warp_anim.step(canvas);
    }

    // Expire finished thrusters from the front (FIFO: oldest spawned first; all
    // share the same LIFE, so the front always expires first).
    while (!thrusters.is_empty() &&
           thrusters.front().age >= ThrusterContext::LIFE) {
      thrusters.pop();
    }

    // Advance and draw each live thruster directly from the ring. Lifetime,
    // radius growth and fade are driven here rather than by per-thruster
    // animations capturing &ctx: ctx is a recyclable ring slot, so an animation
    // outliving the slot's reuse would step/draw whatever thruster later lands
    // in it. Both radius and opacity are pure functions of `age`.
    for (size_t i = 0; i < thrusters.size(); ++i) {
      ThrusterContext &ctx = thrusters[i];
      float progress = static_cast<float>(ctx.age) / ThrusterContext::LIFE;
      float opacity = 1.0f - ease_out_expo(hs::clamp(progress, 0.0f, 1.0f));
      draw_thruster(canvas, ctx, ctx.radius_at(), opacity);
      ++ctx.age;
    }

    t_global++;
  }

private:
  struct ThrusterContext {
    // Visible lifetime in frames.
    static constexpr int LIFE = 16;
    // Ring radius grows from 0 to RADIUS_MAX over the first RADIUS_GROW_FRAMES
    // frames, then holds. Derived from `age` in draw_frame() (see radius_at()).
    static constexpr int RADIUS_GROW_FRAMES = 8;
    static constexpr float RADIUS_MAX = 0.3f;

    Orientation<W> orientation;
    Vector point;
    int age = 0;

    // Trivially copyable: the circular buffer relocates slots by plain
    // memberwise copy. No member holds a reference into this slot, so a copy
    // never has to be rebound.
    void reset(const Orientation<W> &o, const Vector &p) {
      orientation = o;
      point = p;
      age = 0;
    }

    // Eased ring radius for the current age. Mirrors the prior Transition's
    // phasing, which stepped once before the first draw (hence age + 1).
    float radius_at() const {
      float t = hs::clamp(static_cast<float>(age + 1) / RADIUS_GROW_FRAMES,
                          0.0f, 1.0f);
      return RADIUS_MAX * ease_mid(t);
    }
  };

  StaticCircularBuffer<ThrusterContext, 16> thrusters;

  void on_fire_thruster() {
    warp_phase = hs::rand_f() * 2 * PI_F;

    auto r_fn = [this](float t) { return ring_fn(t); };
    Basis basis = make_basis(Quaternion(), ring_vec);
    Vector thrust_point =
        Plot::DistortedRing::fn_point(r_fn, basis, 1.0f, warp_phase);
    Vector thrust_opp =
        Plot::DistortedRing::fn_point(r_fn, basis, 1.0f, warp_phase + PI_F);

    // warp
    warp_anim = Animation::Mutation(
        amplitude, [](float t) { return 0.7f * expf(-2.0f * t); }, 32,
        ease_mid);

    // spin
    Vector thrust_axis =
        cross(orientation.orient(thrust_point), orientation.orient(ring_vec))
            .normalized();
    timeline.add(0, Animation::Rotation<W>(orientation, thrust_axis, 2 * PI_F,
                                           8 * 16, ease_out_expo));

    // spawn
    spawn_thruster(thrust_point);
    spawn_thruster(thrust_opp);
  }

  void spawn_thruster(const Vector &point) {
    if (thrusters.is_full())
      thrusters.pop();
    thrusters.push_back(ThrusterContext());
    thrusters.back().reset(orientation, point);
    // Lifetime, radius and fade are advanced in draw_frame() by iterating the
    // ring rather than by a per-thruster animation capturing this slot
    // (see draw_frame()).
  }

  float ring_fn(float t) {
    // warp_phase is a radian-domain value (rand_f()*2*PI_F, shared with the
    // DistortedRing calls above); sin_wave's phase is in cycles, so convert.
    return sin_wave(-1, 1, 2, warp_phase / PI_F)(t) *
           sin_wave(-1, 1, 3, 0)(static_cast<float>(t_global % 32) / 32.0f) *
           amplitude;
  }

  void draw_thruster(Canvas &c, const ThrusterContext &ctx, float radius,
                     float opacity) {
    Basis basis = make_basis(ctx.orientation.get(), ctx.point);
    auto fragment_shader = [this, opacity](const Vector &, Fragment &f) {
      f.color = Color4(CRGB(255, 255, 255));
      f.color.color = f.color.color * opacity;
      f.color.alpha = opacity * params.alpha;
    };
    Plot::Ring::draw<W, H>(filters, c, basis, radius, fragment_shader);
  }

  void draw_ring(Canvas &c, float opacity) {
    Basis basis = make_basis(orientation.get(), ring_vec);

    auto fragment_shader = [this, opacity](const Vector &v, Fragment &f) {
      Vector axis = orientation.orient(X_AXIS);
      float angle = angle_between(axis, orientation.orient(v));
      f.color = palette.get(angle / PI_F);
      f.color.alpha *= params.alpha * opacity;
    };

    Plot::DistortedRing::draw<W, H>(
        filters, c, basis, params.radius,
        [this](float t) { return ring_fn(t); }, fragment_shader);
  }

  ProceduralPalette palette;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  Vector ring_vec;
  float amplitude;
  float warp_phase;
  int t_global;

  Timeline timeline;
  Animation::Mutation warp_anim;
  Orientation<W> orientation;

  struct Params {
    float radius = 1.0f;
    float alpha = 0.2f;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Thrusters)
