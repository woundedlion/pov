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
      // duration -1 for infinite
      timeline.add(0, Sprite([this](Canvas& canvas, float alpha) {
          this->draw_mesh(canvas, alpha);
      }, -1)); // Infinite duration
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

  Pipeline<W> scan_pipeline;
  Pipeline<W, FilterAntiAlias<W>> plot_pipeline;

  State get_solid(int index) {
      // 0: Tetrahedron, 1: Cube, 2: Octahedron, 3: Dodecahedron, 4: Icosahedron
      index = index % 5;
      State s;
      if (index == 0) {
          Tetrahedron t; s.vertices = t.vertices; s.faces = t.faces;
      } else if (index == 1) {
          Cube c; s.vertices = c.vertices; s.faces = c.faces;
      } else if (index == 2) {
          Octahedron o; s.vertices = o.vertices; s.faces = o.faces;
      } else if (index == 3) {
          Dodecahedron d; s.vertices = d.vertices; s.faces = d.faces;
      } else {
          Icosahedron i; s.vertices = i.vertices; s.faces = i.faces;
      }
      return s;
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
      int next_idx = (current_solid_idx + 1) % 5;
      State next_solid = get_solid(next_idx);
      
      // Add MeshMorph animation to timeline
      // MeshMorph(output_ptr, source, dest, duration, repeat, easing)
      State start_state = current_mesh; // Capture current state (vertices)
      
      // Note: MeshMorph constructor writes *source* to *output* immediately in my modification.
      // This is good.
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
      
  }
};
