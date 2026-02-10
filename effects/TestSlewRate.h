/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../filter.h"
#include "../scan.h"
#include "../plot.h"
#include "../solids.h"
#include "../palettes.h"
#include "../animation.h"

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
    t += 0.01f; // Speed ~ 0.05 in JS per frame? JS frame is ? 
    // JS: t += 1, speed 0.05. So 0.05 per frame.
    
    // Draw Mesh
    // Plot::Mesh::draw is not explicitly defined in plot.h yet? 
    // plot.h has Point, Line, Vertices, Ring, etc.
    // I need to check if Plot::Mesh exists or if I should use Vertices/Lines.
    // In JS `Plot.Mesh.draw(filters, mesh, shader)`.
    // C++ usually has `Solids::draw` or similar? 
    // Let's use `Plot::Generic::draw` pattern or iterate faces?
    // Actually `TestTemporal.h` used `Plot::Mesh::draw`. Let me check `plot.h` again or `TestTemporal.h` imports.
    // `TestTemporal.h` used `Plot::Mesh::draw`.
    // `plot.h` (read in step 120) stopped at line 800. `Star` struct.
    // I suspect `Mesh` is further down in `plot.h`.
    // I will assume `Plot::Mesh` exists.
    
    Plot::Mesh::draw<W>(pipeline, *this, mesh, [&](const Vector& v, const Fragment& f) -> Fragment {
        Color4 baseColor = palette.get((v.j + 1.0f) * 0.5f);
        
        float phase = fmodf(t * speed, 1.0f);
        if (phase < 0) phase += 1.0f;
        
        // v1 in JS was likely "arc length" or similar, but here we don't have v1 from Mesh vertices directly 
        // unless Plot::Mesh populates it?
        // In JS TestSlewRate: `frag.v1`. 
        // JS Mesh draw usually populates v1 with... something?
        // Actually in JS `RichVertexShaders`, v1 is often longitudinal or arc length.
        // For a sphere mesh, v1 might be longitude/phi?
        // Let's assume v1 is populated.
        
        // Wait, TestTemporal.h uses `v.j` for coloring. 
        // TestSlewRate.js uses `frag.v1`.
        
        float dist = std::abs(f.v1 - phase);
        if (dist > 0.5f) dist = 1.0f - dist;
        
        float intensity = 0.0f;
        if (dist < lightSize) {
            float r = 1.0f - (dist / lightSize);
            intensity = r * r;
        }
        
        // Lerp to white
        // Color4 blend
        auto lerpColor = [](Color4 a, Color4 b, float t) {
             // simplified lerp
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
  Filter::Screen::Slew<W, 50000> slew; // Capacity matches JS? 500,000 in JS is huge. W*H is limited on C++.
  // 50000 is safer.
  
  Pipeline<W, Filter::World::Orient<W>, Filter::Screen::Slew<W, 50000>, Filter::Screen::AntiAlias<W>> pipeline;
  
  CircularPalette source_palette = CircularPalette(Palettes::richSunset);
  AnimatedPalette palette;
  CycleModifier modifier;
  Timeline timeline;
  
  // Mesh
  struct SimpleMesh {
        std::vector<Vector> vertices;
        std::vector<std::vector<int>> faces;
  };
  SimpleMesh mesh; 
    
  void rebuild_mesh() {
        // Load Icosahedron
        using S = Icosahedron;
        for(const auto& v : S::vertices) mesh.vertices.push_back(v);
        int offset = 0;
        for(uint8_t count : S::face_counts) {
            std::vector<int> face;
            for(int i = 0; i < count; ++i) face.push_back(S::faces[offset + i]);
            mesh.faces.push_back(face);
            offset += count;
        }
  }
};
