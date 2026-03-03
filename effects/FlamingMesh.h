/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"

template <int W, int H> class FlamingMesh : public Effect {
public:
public:
  struct DelayCalc {
    FlamingMesh *self;
    float operator()(float x, float y) const {
      return self->calculate_delay(x, y);
    }
  };

  FLASHMEM FlamingMesh()
      : Effect(W, H), params(), t(0.0f), orientation(), timeline(), noise(),
        mesh(), palette(Palettes::richSunset),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::Temporal<W, 100000, DelayCalc>(DelayCalc{this},
                                                               2.0f),
                Filter::Screen::AntiAlias<W, H>()) {

    // Initialize with JS defaults
    params.temporalEnabled = true;
    params.windowSize = 8;
    params.delayBase = 10;
    params.delayAmp = 10;
    params.speed = 1.0f;
    params.noiseFreq = 0.125f;

    // Initialize noise
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetFrequency(params.noiseFreq);

    // Initialize mesh (Dodecahedron matches Solids.get(3) or 'dodecahedron')
    ScratchContext ctx(scratch_arena_a, scratch_arena_b);
    mesh = DodecahedronGenerator().generate(geometry_arena, ctx);

    registerParam("Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Delay Base", &params.delayBase, 0.0f, 50.0f);
    registerParam("Delay Amp", &params.delayAmp, 1.0f, 50.0f);
    registerParam("Noise Freq", &params.noiseFreq, 0.01f, 1.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    noise.SetFrequency(params.noiseFreq);
    filters.next.set_window_size(params.windowSize);

    Plot::Mesh::draw<W, H>(filters, canvas, mesh,
                           [&](const Vector &v, Fragment &f) {
                             float t_val = (v.j + 1.0f) * 0.5f;
                             f.color = palette.get(t_val);
                           });

    filters.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
        1.0f);

    t += 1.0f;
  }

  // Parameters matching JS
  struct Params {
    bool temporalEnabled;
    int windowSize;
    float delayBase;
    float delayAmp;
    float speed;
    float noiseFreq;
  } params;

private:
  float t = 0.0f;

  Orientation<W> orientation;
  Timeline<W> timeline;
  FastNoiseLite noise;
  PolyMesh mesh;
  ProceduralPalette palette;

  float calculate_delay(float x, float y) {
    if (!params.temporalEnabled)
      return 0.0f;
    float scale = static_cast<float>(W) / (2.0f * PI_F);
    Vector v = pixel_to_vector<W, H>(x, y) * scale;
    float noiseVal = noise.GetNoise(v.i, v.j, v.k + (t * params.speed));
    return std::max(0.0f, params.delayBase +
                              (noiseVal * 0.5f + 0.5f) * params.delayAmp);
  }

  Pipeline<W, H, Filter::World::Orient<W>,
           Filter::Screen::Temporal<W, 100000, DelayCalc>,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};
