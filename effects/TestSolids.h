/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"
#include "../geometry.h"
#include "../solids.h"
#include "../scan.h"
#include "../plot.h"

template <int W>
class TestSolids : public Effect {
public:
  // Using MeshMorph::State directly to ensure compatibility
  using State = MeshMorph::State;

  TestSolids() : Effect(W), intensity(1.2f), opacity(1.0f), debug_bb(false), hankin(true), dual(false) {
      register_solids();
      rebuild();
  }

  void rebuild() {
      timeline = Timeline();
      orientation = Orientation();
      
      // Initial Solid
      current_solid_idx = 0;
      current_mesh = get_solid(current_solid_idx);
      
      // Start Rotation
      // Rotate around UP over 10 seconds (600 frames)
      timeline.add(0, Rotation<W>(orientation, UP, 2 * PI_F, 600, ease_mid, true)); // Loop
      
      // Start Morph Cycle
      start_morph_sequence();
      
      // Render Task
      // Sprite(draw_fn, duration, fade_in, fade_in_easing, fade_out, fade_out_easing)
      // duration -1 for infinite (using std::function constructor)
      timeline.add(0, Sprite(std::function<void(Canvas&, float)>([this](Canvas& canvas, float alpha) {
          this->draw_mesh(canvas, alpha);
      }), -1)); 
  }

  void draw_frame() override {
      Canvas canvas(*this);
      timeline.step(canvas);
  }
  
  bool show_bg() const override { return false; }

private:
  Timeline timeline;
  Orientation orientation;
  
  State current_mesh;
  int current_solid_idx;
  
  // Params
  float intensity;
  float opacity;
  bool debug_bb;
  bool hankin;
  bool dual;
  float hankin_angle = PI_F / 4.0f;

  std::vector<std::function<State()>> solid_generators;

  Pipeline<W> scan_pipeline;
  Pipeline<W, FilterAntiAlias<W>> plot_pipeline;

  template <typename T>
  void add_solid() {
      solid_generators.push_back([]() -> State {
          T t;
          State s;
          s.vertices = t.vertices;
          s.faces = t.faces;
          return s;
      });
  }

  void register_solids() {
      solid_generators.clear();
      add_solid<Tetrahedron>();               // 1
      add_solid<Cube>();                      // 2
      add_solid<Octahedron>();                // 3
      add_solid<Icosahedron>();               // 4
      add_solid<Dodecahedron>();              // 5
      add_solid<TruncatedTetrahedron>();      // 6
      add_solid<Cuboctahedron>();             // 7
      add_solid<TruncatedCube>();             // 8
      add_solid<TruncatedOctahedron>();       // 9
      add_solid<Rhombicuboctahedron>();       // 10
      add_solid<TruncatedCuboctahedron>();    // 11
      add_solid<SnubCube>();                  // 12
      add_solid<Icosidodecahedron>();         // 13
      add_solid<TruncatedDodecahedron>();     // 14
      add_solid<TruncatedIcosahedron>();      // 15
      add_solid<Rhombicosidodecahedron>();    // 16
      add_solid<SnubDodecahedron>();          // 17
  }

  State get_solid(int index) {
      if (solid_generators.empty()) return State();
      index = index % solid_generators.size();
      return solid_generators[index]();
  }

  void start_morph_sequence() {
      // Morph every 2 seconds (128 frames)
      // PeriodicTimer(period, callback, repeat)
      timeline.add(0, PeriodicTimer(128, [this](Canvas& c) {
          this->trigger_morph();
      }, true));
  }
  
  void trigger_morph() {
      // Next solid
      if (solid_generators.empty()) return;
      int next_idx = (current_solid_idx + 1) % solid_generators.size();
      State next_solid = get_solid(next_idx);
      
      // Add MeshMorph animation to timeline
      // MeshMorph(output_ptr, source, dest, duration, repeat, easing)
      State start_state = current_mesh; // Capture current state (vertices)
      
      timeline.add(0, MeshMorph(&current_mesh, start_state, next_solid, 64, false, ease_in_out_sin));
      
      current_solid_idx = next_idx;
  }

  void draw_mesh(Canvas& canvas, float alpha) {
      // Transform vertices
      VertexList transformed_verts = orientation.orient(current_mesh.vertices);
      
      // Temporary struct to pass to Scan::Mesh
      struct OrientedMesh {
          const std::vector<Vector>& vertices;
          const std::vector<std::vector<int>>& faces;
      } mesh_ref = { transformed_verts, current_mesh.faces };
      
      auto color_fn = [&](const Vector& p, float t, float d, int faceIdx) {
          // Select palette based on N (sides)
          const Palette* palette = &richSunset;
          
          size_t n = 0;
          if (faceIdx >= 0 && faceIdx < (int)mesh_ref.faces.size()) {
              n = mesh_ref.faces[faceIdx].size();
          }

          if (n == 3) palette = &lavenderLake;
          else if (n == 4) palette = &lavenderLake;
          else if (n == 5) palette = &emeraldForest;
          
          float distFromEdge = -d; // d is negative inside
          float val = std::clamp(distFromEdge * intensity, 0.0f, 1.0f);
          
          Color4 c = palette->get(val);
          c.alpha *= opacity;
          return c;
      };
      
      Scan<W>::Mesh::draw(scan_pipeline, canvas, mesh_ref, color_fn); 
  }
};
