/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

template <int W>
class TestSlewRate : public Effect {
public:
  TestSlewRate()
    : Effect(W),
      orientation(),
      // Pipeline: Orient -> Slew -> AntiAlias
      slew(1.0f, 0.05f),
      pipeline(
        Filter::World::Orient<W>(orientation),
        slew,
        Filter::Screen::AntiAlias<W>()
      ),
      palette(CircularPalette(Palettes::richSunset)),
      modifier(0.02f)
  {
      timeline.add(0, Animation::RandomWalk<W>(orientation, Vector(0, 1, 0)));
      palette.add(&modifier);
      timeline.add(0, Animation::PaletteAnimation(modifier));

      // Load Mesh (Icosahedron default)
      rebuild_mesh();
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    timeline.step();
    t += 0.01f; 
    
    Plot::Mesh::draw<W>(pipeline, *this, mesh, [&](const Vector& v, const Fragment& f) -> Fragment {
        Color4 baseColor = palette.get((v.j + 1.0f) * 0.5f);
        
        float phase = fmodf(t * speed, 1.0f);
        if (phase < 0) phase += 1.0f;

        
        float dist = std::abs(f.v1 - phase);
        if (dist > 0.5f) dist = 1.0f - dist;
        
        float intensity = 0.0f;
        if (dist < lightSize) {
            float r = 1.0f - (dist / lightSize);
            intensity = r * r;
        }
        

        auto lerpColor = [](Color4 a, Color4 b, float t) {
             return Color4(
                static_cast<uint8_t>(a.color.r + (b.color.r - a.color.r) * t),
                static_cast<uint8_t>(a.color.g + (b.color.g - a.color.g) * t),
                static_cast<uint8_t>(a.color.b + (b.color.b - a.color.b) * t),
                a.alpha + (b.alpha - a.alpha) * t
             );
        };
        
        Fragment out = f;
        out.color = lerpColor(baseColor, Color4(255, 255, 255, 1.0f), intensity);
        return out;
    });

    pipeline.flush(*this, [](float x, float y, float t) { return Color4(0,0,0,0); }, 1.0f); 
  }
  
  // Params
  float t = 0;
  float speed = 0.05f;
  float lightSize = 0.15f;

private:
  Orientation orientation;
  Filter::Screen::Slew<W, 50000> slew; 
  
  Pipeline<W, Filter::World::Orient<W>, Filter::Screen::Slew<W, 50000>, Filter::Screen::AntiAlias<W>> pipeline;
  
  CircularPalette source_palette = CircularPalette(Palettes::richSunset);
  AnimatedPalette palette;
  CycleModifier modifier;
  Timeline timeline;
  
  // Mesh
  PolyMesh mesh; 
    
  void rebuild_mesh() {
        mesh = Solids::Platonic::icosahedron();
  }
};
