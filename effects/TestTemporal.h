/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <vector>
#include <functional>
#include "../effects_engine.h"

template <int W, int H> class TestTemporal : public Effect {
public:
  enum DelayMode {
    VerticalWave = 0,
    DiagonalSpiral = 1,
    LiquidTime = 2,
    QuantumTunnel = 3,
    Datamosh = 4
  };

public:
  struct DelayCalc {
    TestTemporal *self;
    float operator()(float x, float y) const {
      return self->calculate_delay(x, y);
    }
  };

  TestTemporal()
      : Effect(W, H), orientation(), noise(),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::Temporal<W, 100000, DelayCalc>(DelayCalc{this},
                                                               2.0f),
                Filter::Screen::AntiAlias<W, H>()),
        circular_source(Palettes::richSunset), palette(circular_source),
        modifier(0.02f) {
    registerParam("Mode", &params.mode, 0.0f, 4.0f);
    registerParam("Delay Base", &params.base, 0.0f, 50.0f);
    registerParam("Delay Amp", &params.amp, 0.0f, 50.0f);
    registerParam("Speed", &params.speed, 0.0f, 0.1f);
    registerParam("Param X", &params.paramX, 0.0f, 20.0f);
    registerParam("Param Y", &params.paramY, 0.0f, 20.0f);

    this->persist_pixels = false;

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    timeline.add(0, Animation::RandomWalk<W>(orientation, Vector(0, 1, 0)));
    palette.add(&modifier);
    timeline.add(0, Animation::PaletteAnimation(modifier));

    rebuild_mesh();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    t += 1.0f;

    noise.SetFrequency(0.03f);
    filters.next.set_window_size(2);

    Plot::Mesh::draw<W, H>(
        filters, canvas, mesh, [&](const Vector &v, Fragment &f) {
          float t_val = (v.j + 1.0f) * 0.5f;
          f.color = palette.get(t_val);

          // Lighting Logic (Edge Pulses)
          float phase = fmodf(t * 0.05f, 1.0f); // Fixed light speed
          if (phase < 0)
            phase += 1.0f;
          float dist = std::abs(f.v1 - phase);
          if (dist > 0.5f)
            dist = 1.0f - dist;

          const float width = 0.15f;
          if (dist < width) {
            float strength = powf(1.0f - (dist / width), 2.0f);

            // Lerp to white
            Pixel white(65535, 65535, 65535);
            float val = strength * 1.0f; // Fixed light alpha
            if (val > 1.0f)
              val = 1.0f;
            uint16_t frac = (uint16_t)(val * 65535.0f);
            f.color.color = f.color.color.lerp16(white, frac);
          }
        });

    filters.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
        1.0f);
  }

  // Generic Params for all modes
  struct Params {
    float mode = 0.0f; // VerticalWave, DiagonalSpiral, etc.
    float base = 8.0f;
    float amp = 4.0f;
    float speed = 0.01f;
    float paramX = 0.3f; // Generic param 1
    float paramY = 2.0f; // Generic param 2
  } params;

private:
  float t = 0;

  CircularPalette circular_source;
  AnimatedPalette palette;
  CycleModifier modifier;

  Orientation<W> orientation;
  Timeline<W> timeline;
  FastNoiseLite noise;
  Pipeline<W, H, Filter::World::Orient<W>,
           Filter::Screen::Temporal<W, 100000, DelayCalc>,
           Filter::Screen::AntiAlias<W, H>>
      filters;
  PolyMesh mesh;

  void rebuild_mesh() {
    ArenaMarker _(scratch_arena_a);
    ScratchContext ctx(scratch_arena_a, scratch_arena_b);
    mesh = Solids::finalize_solid(Solids::Platonic::icosahedron(ctx),
                                  geometry_arena);
  }

  float calculate_delay(float x, float y) {
    // Map mode index
    int modeIdx = (int)params.mode;

    switch (modeIdx) {
    case VerticalWave: {
      // X = Freq (0.3)
      float phase = y * params.paramX + t * params.speed;
      return std::max(0.0f, params.base + sinf(phase) * params.amp);
    }
    case DiagonalSpiral: {
      // X = xSpirals (2.0), Y = yFreq (0.3)
      float xPhase = (x / W) * PI_F * 2.0f * params.paramX;
      float yPhase = y * params.paramY;
      float phase = xPhase + yPhase + t * params.speed;
      return std::max(0.0f, params.base + sinf(phase) * params.amp);
    }
    case LiquidTime: {
      // X = TimeScale (2.0)
      // 3D Noise
      float noiseVal = noise.GetNoise(x, y, t * params.paramX);
      return params.base + (noiseVal + 1.0f) * 0.5f * params.amp;
    }
    case QuantumTunnel: {
      // X = Tightness (10.0), Y = Angle (5.0)
      // Polar
      float u = (x / W) * 2.0f - 1.0f;
      float v = (y / H) * 2.0f - 1.0f;
      float radius = sqrtf(u * u + v * v);
      float angle = atan2f(v, u);
      float spiral = sinf(radius * params.paramX - angle * params.paramY +
                          t * params.speed);
      return std::max(0.0f, params.base + (spiral + 1.0f) * params.amp);
    }
    case Datamosh: {
      // X = FlowSpeed (0.1), Y = GlitchScale (15.0)
      float flow = sinf(y * params.paramX + t * 0.05f);
      int blockSize = 8;
      int column = (int)(x / blockSize);
      float glitchOffset = sinf(column * 12.9898f) * params.paramY;
      float total = flow * 10.0f + glitchOffset;
      return std::max(0.0f, fmodf(params.base + std::abs(total), params.amp));
    }
    }
    return 0.0f;
  }
};
