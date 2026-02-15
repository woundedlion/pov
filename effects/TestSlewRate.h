/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H>
class TestSlewRate : public Effect {
public:
  TestSlewRate()
    : Effect(W, H),
      orientation(),
      // Pipeline: Orient -> Slew -> AntiAlias
      pipeline(
        Filter::World::Orient<W>(orientation),
        Filter::Screen::Slew<W, 200000>(1.0f, 0.03f)
      ),
      palette(source_palette),
      modifier(0.02f)
  {
      this->persist_pixels = false;

      timeline.add(0, Animation::RandomWalk<W>(orientation, Vector(0, 1, 0)));
      palette.add(&modifier);
      timeline.add(0, Animation::PaletteAnimation(modifier));

      // Load Mesh (Icosahedron default)
      rebuild_mesh();
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    t += 1.0f; 
    
    // Filters need window size update if we fully port it, but for now just lights:
    
    auto fragmentShader =  [&] (const Vector& v, Fragment& f) {
        Color4 baseColor = palette.get((v.j + 1.0f) * 0.5f);
        f.color = baseColor;
        
        // Lighting Logic (Edge Pulses)
         float phase = fmodf(t * global.lightSpeed, 1.0f);
         if (phase < 0) phase += 1.0f;
         float dist = std::abs(f.v1 - phase);
         if (dist > 0.5f) dist = 1.0f - dist;
         
         const float width = 0.15f;
         if (dist < width) {
             float strength = powf(1.0f - (dist / width), 2.0f);
             
             // Lerp to white
             Pixel white(65535, 65535, 65535);
             float val = strength * global.lightAlpha;
             if (val > 1.0f) val = 1.0f;
             uint16_t frac = (uint16_t)(val * 65535.0f);
             f.color.color = f.color.color.lerp16(white, frac);
         }
    };

    Plot::Mesh::draw<W, H>(pipeline, canvas, mesh, fragmentShader); 
    pipeline.flush(canvas, [](float x, float y, float t) { return Color4(0,0,0,0); }, 1.0f); 
  }
  
  // Params
  struct GlobalParams {
      float lightSpeed = 0.05f;
      float lightAlpha = 1.0f;
  } global;
  
  float t = 0;

private:
  Orientation<W> orientation;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::Slew<W, 200000>, Filter::Screen::AntiAlias<W, H>> pipeline;
  
  CircularPalette source_palette = CircularPalette(Palettes::richSunset);
  AnimatedPalette palette;
  CycleModifier modifier;
  Timeline<W> timeline;
  
  // Mesh
  PolyMesh mesh; 
    
  void rebuild_mesh() {
        mesh = Solids::Platonic::icosahedron();
  }
};
