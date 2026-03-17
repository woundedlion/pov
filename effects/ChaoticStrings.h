/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "../effects_engine.h"

template <int W, int H> class ChaoticStrings : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115;

  struct Params {
    float alpha = 1.0f;
    float cycle_duration = 80.0f;
    float jitterAmp = 3.0f;
    float speed = 0.1f;
    float noiseFreq = 0.33f;
    float scaleFactor = 200.0f;
    float cycleSpeed = 0.1f;

    void lerp(const Params &a, const Params &b, float t) {
      alpha = ::lerp(a.alpha, b.alpha, t);
      cycle_duration = ::lerp(a.cycle_duration, b.cycle_duration, t);
      jitterAmp = ::lerp(a.jitterAmp, b.jitterAmp, t);
      speed = ::lerp(a.speed, b.speed, t);
      noiseFreq = ::lerp(a.noiseFreq, b.noiseFreq, t);
      scaleFactor = ::lerp(a.scaleFactor, b.scaleFactor, t);
      cycleSpeed = ::lerp(a.cycleSpeed, b.cycleSpeed, t);
    }
  } params;

  Presets<Params, 1> preset_manager = {{{{"gloopy",
                                          {/* alpha */ 1.0f,
                                           /* cycle duration */ 80.0f,
                                           /* jitter amp */ 1.7f,
                                           /* speed */ 0.04f,
                                           /* noise freq */ 0.32f}}}}};

  struct LissajousConfig {
    float m1;
    float m2;
    float a;
    float domain;
  };

  struct Node {
    Orientation<W, 16> orientation;
    Animation::OrientationTrail<Orientation<W, 16>, TRAIL_LENGTH> trail;
    Vector v;

    Node() : v(Y_AXIS) {}
  };

  FLASHMEM ChaoticStrings()
      : Effect(W, H), timeline(), filters(Filter::Screen::AntiAlias<W, H>()),
        path([this](float t) { return Vector(0, 1, 0); }), orientation(),
        palette_variant(),
        cycle_phase(0.0f), functions{{{12.0f, 5.0f, 0, 2 * PI_F}}},
        cur_function_idx(0), noise_xform(timeline) {}

  void init() override {

    configure_arenas(GLOBAL_ARENA_SIZE - 200 * 1024, 200 * 1024, 0);

    node = static_cast<Node *>(
        persistent_arena.allocate(sizeof(Node), alignof(Node)));
    new (node) Node();

    // Colors
    static_palette.bind(&palette_variant, &scale_mod, &cycle_mod);
    palette_variant = Palettes::fireAndIce;

    // Parameters
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Jitter Amp", &params.jitterAmp, 0.0f, 10.0f);
    registerParam("Noise Freq", &params.noiseFreq, 0.01f, 10.0f);
    registerParam("Scale Factor", &params.scaleFactor, 1.0f, 500.0f);
    registerParam("Cycle Speed", &params.cycleSpeed, 0.0f, 1.0f);

    // Initialize noise params
    noise_xform.params.amplitude = params.jitterAmp;
    noise_xform.params.frequency = params.noiseFreq;
    noise_xform.params.speed = params.speed;
    noise_xform.params.sync();

    // Initialize Lissajous functions

    update_path();
    noise_xform.spawn(0, -1);

    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    timeline.add(0, Animation::Motion<W, 16>(node->orientation, path,
                                             (int)params.cycle_duration, true));
    timeline.add(0, Animation::PeriodicTimer(
                        2 * (int)params.cycle_duration,
                        [this](Canvas &c) {
                          cur_function_idx =
                              (cur_function_idx + 1) % functions.size();
                          update_path();
                        },
                        true));
    timeline.add(0, Animation::PeriodicTimer(
                        4 * (int)params.cycle_duration,
                        [this](Canvas &c) {
                          cur_function_idx = static_cast<int>(
                              hs::rand_int(0, functions.size()));
                          next_preset();
                        },
                        true));

    driver_ =
        timeline.add_get(0, Animation::Driver(cycle_phase, params.cycleSpeed));

    params = preset_manager.get();
  }

  void next_preset() {
    // TODO: Re-enable with smaller Params if needed
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    ScratchScope _(scratch_arena_a);
    ArenaVector<Fragment> vertices(scratch_arena_a, 2000);
    timeline.step(canvas);

    noise_xform.params.frequency = params.noiseFreq;
    noise_xform.params.amplitude = params.jitterAmp;
    noise_xform.params.speed = params.speed;
    noise_xform.params.sync();

    // Update active noise entities to reflect live parameter changes
    for (auto &e : noise_xform.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
        e.params.sync();
      }
    }

    // Update cycle speed dynamically
    if (driver_ && &driver_->get_mutant() == &cycle_phase) {
      driver_->set_speed(params.cycleSpeed);
    }

    for (auto &e : noise_xform.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
        e.params.sync();
      }
    }

    // Record the current orientation snapshot
    node->trail.record(node->orientation);

    float current_t = timeline.t;
    constexpr float MAX_TRAIL = ChaoticStrings<W, H>::TRAIL_LENGTH;

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      // Re-apply transforms at draw time so the trail undulates with noise
      Vector pos =
          noise_xform.transform(orientation.orient(rotate(node->v, q)));
      Fragment f;
      f.pos = pos.normalized();
      f.age = current_t * t;
      f.v3 = t;
      vertices.push_back(f);
    });

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      float color_t = frag.age / MAX_TRAIL;
      frag.color = static_palette.get(color_t);
      frag.color.alpha *= quintic_kernel(frag.v3);
    };

    Plot::Multiline::draw<W, H>(filters, canvas, vertices, fragment_shader);
  }

private:
  void update_path() {
    const LissajousConfig &config = functions[cur_function_idx];
    path.f = [&config](float t) {
      return lissajous(config.m1, config.m2, config.a, t * config.domain);
    };
  }

  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ProceduralPath path;
  Orientation<W> orientation;
  ScaleModifier scale_mod{200.0f, &params.scaleFactor};
  CycleModifier cycle_mod{&cycle_phase};
  ProceduralPalette palette_variant;
  StaticPalette<ProceduralPalette, ScaleModifier, CycleModifier> static_palette;
  Animation::Driver *driver_ = nullptr;
  float cycle_phase = 0.0f;
  std::array<LissajousConfig, 1> functions;
  int cur_function_idx;
  Node *node = nullptr;
  NoiseTransformer<W, 1> noise_xform;
};

#include "../effect_registry.h"
REGISTER_EFFECT(ChaoticStrings)
