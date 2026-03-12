/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"

template <int W, int H> class MeshFeedback : public Effect {
public:
  struct Params {
    float fade = 0.93f;
    float amplitude = 0.9f;
    float frequency = 0.33f;
    float speed = 0.8f;
    float scale = 22.6f;

    void lerp(const Params &a, const Params &b, float t) {
      fade = ::lerp(a.fade, b.fade, t);
      amplitude = ::lerp(a.amplitude, b.amplitude, t);
      frequency = ::lerp(a.frequency, b.frequency, t);
      speed = ::lerp(a.speed, b.speed, t);
      scale = ::lerp(a.scale, b.scale, t);
    }
  } params;

  FLASHMEM MeshFeedback()
      : Effect(W, H), noise_params(), orientation(), timeline(), mesh(),
        palette(Palettes::richSunset),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H, TransformerFn>(
                    TransformerFn{this}, 0.95f)) {}

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 16384 - 2048, 16384, 2048);

    params = presets.get();
    apply_params();

    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.sync();

    filters.next.next.init_storage(Persistent(persistent_arena));

    mesh = generate_mesh<DodecahedronGenerator>(persistent_arena);

    registerParam("Fade", &params.fade, 0.5f, 0.99f);
    registerParam("Distort Amp", &params.amplitude, 0.0f, 30.0f);
    registerParam("Distort Freq", &params.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &params.scale, 0.1f, 50.0f);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));

    timeline.add(0, Animation::PeriodicTimer(
                        150,
                        [this](Canvas &) {
                          presets.next();
                          timeline.add(0, Animation::Lerp(params, presets.get(),
                                                          48, ease_mid));
                        },
                        true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    apply_params();
    noise_params.sync();
    filters.next.next.set_fade(params.fade);

    Plot::Mesh::draw<W, H>(filters, canvas, mesh,
                           [&](const Vector &v, Fragment &f) {
                             float t_val = (v.y + 1.0f) * 0.5f;
                             f.color = palette.get(t_val);
                           });

    filters.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
        1.0f);
  }

private:
  void apply_params() {
    noise_params.amplitude = params.amplitude;
    noise_params.frequency = params.frequency;
    noise_params.speed = params.speed;
    noise_params.scale = params.scale;
  }

  Presets<Params, 4> presets = {
      .entries = {{{"Flames", {0.94f, 0.51f, 0.42f, 0.46f, 23.0f}},
                   {"Organic", {0.93f, 0.9f, 0.33f, 0.8f, 22.6f}},
                   {"Intense", {0.96f, 5.0f, 0.15f, 1.5f, 10.0f}},
                   {"Subtle", {0.90f, 0.3f, 0.60f, 0.3f, 35.0f}}}},
      .current_idx = 0};
  NoiseParams noise_params;

  Orientation<W> orientation;
  Timeline<W> timeline;
  PolyMesh mesh;
  ProceduralPalette palette;

  struct TransformerFn {
    const MeshFeedback *self;
    Vector operator()(const Vector &v) const {
      return noise_transform(v, self->noise_params).normalize();
    }
  };

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H, TransformerFn>>
      filters;
};
