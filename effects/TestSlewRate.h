/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H> class TestSlewRate : public Effect {
public:
  FLASHMEM TestSlewRate()
      : Effect(W, H), orientation(),
        // Pipeline: Orient -> Slew -> AntiAlias
        pipeline(Filter::World::Orient<W>(orientation),
                 Filter::Screen::Slew<W, 200000>(1.0f, 0.03f)),
        palette(&circular_source) {

    registerParam("Light Speed", &params.lightSpeed, 0.0f, 0.5f);
    registerParam("Light Alpha", &params.lightAlpha, 0.0f, 2.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, Vector(0, 1, 0)));
    palette.add(CycleModifier(&color_offset));
    timeline.add(0, Animation::Driver(color_offset, 0.02f));
    rebuild_mesh();
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    t += 1.0f;

    auto fragmentShader = [&](const Vector &v, Fragment &f) {
      Color4 baseColor = get_color(palette.get((v.j + 1.0f) * 0.5f), 1.0f);
      f.color = baseColor;

      // Lighting Logic (Edge Pulses)
      float phase = fmodf(t * params.lightSpeed, 1.0f);
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
        float val = strength * params.lightAlpha;
        if (val > 1.0f)
          val = 1.0f;
        uint16_t frac = (uint16_t)(val * 65535.0f);
        f.color.color = f.color.color.lerp16(white, frac);
      }
    };

    Plot::Mesh::draw<W, H>(pipeline, canvas, mesh, fragmentShader);
    pipeline.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
        1.0f);
  }

  // Params
  struct Params {
    float lightSpeed = 0.05f;
    float lightAlpha = 1.0f;
  } params;

  float t = 0;

private:
  Orientation<W> orientation;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::Slew<W, 200000>,
           Filter::Screen::AntiAlias<W, H>>
      pipeline;

  float color_offset = 0.0f;
  PaletteVariant source_variant{Palettes::richSunset};
  PaletteVariant circular_source{CircularPalette(&source_variant)};
  AnimatedPalette palette;
  Timeline<W> timeline;

  // Mesh
  PolyMesh mesh;

  void rebuild_mesh() {
    MemoryCtx ctx;
    ScopedScratch _(ctx.get_scratch_front());
    mesh = Solids::finalize_solid(Solids::Platonic::icosahedron(ctx),
                                  persistent_arena);
  }
};
