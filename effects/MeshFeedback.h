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
  FLASHMEM MeshFeedback()
      : Effect(W, H), params(), noise_params(), orientation(), timeline(),
        mesh(), palette(Palettes::richSunset),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H, TransformerFn>(
                    TransformerFn{this}, 0.95f)) {}

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 16384 - 2048, 16384, 2048);
    params.fade = 0.93f;
    noise_params.amplitude = 0.9f;
    noise_params.frequency = 0.33f;
    noise_params.speed = 1.0f;
    noise_params.scale = 22.6f;

    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.sync();

    filters.next.next.init_storage(Persistent(persistent_arena));

    // Initialize mesh (Dodecahedron matches Solids.get(3) or 'dodecahedron')
    mesh = generate_mesh<DodecahedronGenerator>(persistent_arena);

    registerParam("Fade", &params.fade, 0.5f, 0.99f);
    registerParam("Distort Amp", &noise_params.amplitude, 0.0f, 5.0f);
    registerParam("Distort Freq", &noise_params.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &noise_params.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &noise_params.scale, 0.1f, 50.0f);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

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

  struct Params {
    float fade;
  } params;

private:
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
