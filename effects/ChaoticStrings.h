/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <vector>
#include <array>
#include "../effects_engine.h"
#include "../waves.h"

template <int W, int H> class ChaoticStrings : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115;

  struct Params {
    float alpha = 1.0f;
    float cycle_duration = 80.0f;
    float jitterAmp = 3.0f;
    float speed = 0.1f;
    float noiseFreq = 0.33f;

    void lerp(const Params &a, const Params &b, float t) {
      alpha = ::lerp(a.alpha, b.alpha, t);
      cycle_duration = ::lerp(a.cycle_duration, b.cycle_duration, t);
      jitterAmp = ::lerp(a.jitterAmp, b.jitterAmp, t);
      speed = ::lerp(a.speed, b.speed, t);
      noiseFreq = ::lerp(a.noiseFreq, b.noiseFreq, t);
    }
  } params;

  Presets<Params> preset_manager = {{"gloopy",
                                     {/* alpha */ 1.0f,
                                      /* cycle duration */ 80.0f,
                                      /* jitter amp */ 1.7f,
                                      /* speed */ 0.04f,
                                      /* noise freq */ 0.32f}}};

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

  ChaoticStrings()
      : Effect(W, H), cur_function_idx(0),
        filters(Filter::Screen::AntiAlias<W, H>()),
        path([this](float t) { return Vector(0, 1, 0); }),
        animated_palette(&palette_variant), ripple_phase(0.0f),
        cycle_offset(0.0f), noise(timeline) {

    animated_palette.add(ScaleModifier(1.0f)).add(CycleModifier(&cycle_offset));
    palette_variant = Palettes::fireAndIce;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Jitter Amp", &params.jitterAmp, 0.0f, 10.0f);
    registerParam("Noise Freq", &params.noiseFreq, 0.01f, 10.0f);

    // Initialize noise params
    noise.params.amplitude = params.jitterAmp;
    noise.params.frequency = params.noiseFreq;
    noise.params.speed = params.speed;
    noise.params.sync();

    persist_pixels = false;

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

    timeline.add(0, Animation::Mutation(ripple_phase,
                                        sin_wave(0.0f, 10.0f, 0.05f, 0.0f), -1,
                                        ease_mid, true));
    //    timeline.add(0, Animation::Driver(cycle_offset, 0.001f));

    params = preset_manager.get();
  }

  void next_preset() {
    return;
    preset_manager.next();
    Params target = preset_manager.get();

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

    node.trail.record(node.orientation);
    vertices.clear();
    deep_tween(node.trail, [&](const Quaternion &q, float t) {
      Vector v_local = rotate(node.v, q);
      Vector v_final = orientation.orient(v_local);
      v_final = noise.transform(v_final);
      vertices.emplace_back(v_final);
    });

    auto fragment_shader = [&](const Vector &v, Fragment &frag) {
      frag.color = animated_palette.get(frag.v1);
      frag.color.alpha *= quintic_kernel(frag.v0);
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
  ProceduralPath<std::function<Vector(float)>> path;
  Orientation<W> orientation;
  PaletteVariant palette_variant;
  AnimatedPalette animated_palette;
  float ripple_phase;
  float cycle_offset;
  std::vector<LissajousConfig> functions;
  int cur_function_idx;
  Node node;
  NoiseTransformer<W, 1> noise;
  StaticCircularBuffer<Vector, 20000> vertices;
};