/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <vector>
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
    Orientation<W> orientation;
    Animation::OrientationTrail<Orientation<W>, TRAIL_LENGTH> trail;
    Vector v;

    Node() : v(Y_AXIS) {}
  };

  FLASHMEM ChaoticStrings()
      : Effect(W, H), timeline(), filters(Filter::Screen::AntiAlias<W, H>()),
        path([this](float t) { return Vector(0, 1, 0); }), orientation(),
        palette_variant(), animated_palette(&palette_variant),
        cycle_phase(0.0f), functions(), cur_function_idx(0), node(),
        noise(timeline), vertices() {

    // Colors
    animated_palette.add(ScaleModifier(200.0f, &params.scaleFactor))
        .add(CycleModifier(&cycle_phase));
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
    noise.params.amplitude = params.jitterAmp;
    noise.params.frequency = params.noiseFreq;
    noise.params.speed = params.speed;
    noise.params.sync();

    // Initialize Lissajous functions
    functions = {{12.0f, 5.0f, 0, 2 * PI_F}};

    update_path();
    noise.spawn(0, -1);

    timeline.add(0, Animation::RandomWalk<W>(orientation, random_vector()));
    timeline.add(0, Animation::Motion<W>(node.orientation, path,
                                         (int)params.cycle_duration, true));
    timeline.add(0, Animation::PeriodicTimer(
                        2 * (int)params.cycle_duration,
                        [this](Canvas &c) {
                          cur_function_idx = static_cast<int>(
                              hs::rand_int(0, functions.size()));
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

    timeline.add(0, Animation::Driver(cycle_phase, params.cycleSpeed));

    params = preset_manager.get();
  }

  void next_preset() {
    return;
    preset_manager.next();
    timeline.add(
        0, Animation::Lerp(params, preset_manager.get(), 160, ease_in_out_sin));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    noise.params.frequency = params.noiseFreq;
    noise.params.amplitude = params.jitterAmp;
    noise.params.speed = params.speed;
    noise.params.sync();

    // Update active noise entities to reflect live parameter changes
    for (auto &e : noise.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
        e.params.sync();
      }
    }

    // Update cycle speed dynamically
    for (int i = 0; i < timeline.num_events; ++i) {
      auto &e = timeline.events[i];
      if (auto *driver = std::get_if<Animation::Driver>(&e.animation)) {
        if (&driver->get_mutant() == &cycle_phase) {
          driver->set_speed(params.cycleSpeed);
        }
      }
    }
    for (auto &e : noise.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
        e.params.sync();
      }
    }

    node.trail.record(node.orientation);
    vertices.clear();

    float current_t = timeline.t;
    float trail_len = node.trail.length();
    constexpr float MAX_TRAIL = ChaoticStrings<W, H>::TRAIL_LENGTH;

    deep_tween(node.trail, [&](const Quaternion &q, float t) {
      Fragment f;
      f.pos = noise.transform(orientation.orient(rotate(node.v, q)));
      f.age = current_t - trail_len + 1.0f + (t * trail_len);
      f.v3 = t;
      vertices.push_back(f);
    });

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      float color_t = frag.age / MAX_TRAIL;
      frag.color = animated_palette.get(color_t);
      frag.color.alpha *= quintic_kernel(frag.v3);
    };

    Plot::Multiline::draw<W, H>(filters, canvas, vertices, fragment_shader);
  }

private:
  void update_path() {
    const auto &config = functions[cur_function_idx];
    path.f = [=](float t) {
      return lissajous(config.m1, config.m2, config.a,
                       ease_mid(t) * config.domain);
    };
  }

  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ProceduralPath path;
  Orientation<W> orientation;
  PaletteVariant palette_variant;
  AnimatedPalette animated_palette;
  float cycle_phase = 0.0f;
  std::vector<LissajousConfig> functions;
  int cur_function_idx;
  Node node;
  NoiseTransformer<W, 1> noise;
  StaticCircularBuffer<Fragment, 20000> vertices;
};