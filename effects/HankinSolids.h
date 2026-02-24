/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "../solids.h"

template <int W, int H> class HankinSolids : public Effect {
public:
  // Buffer for morphing operations
  Animation::MorphBuffer morph_buffer;

  // Active compiled states
  CompiledHankin primary_hankin;
  CompiledHankin secondary_hankin;
  MeshState primary_mesh;
  MeshState secondary_mesh;

  HankinSolids() : Effect(W, H), filters() {
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);

    persist_pixels = false;

    // Continuous Random Walk
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS,
                        Animation::RandomWalk<W>::Options::Languid()));

    // Init State
    solid_idx = 0;

    // Pre-allocate geometry buffers to guarantee purely zero-overhead swapping
    preallocate_buffers();

    // Initial Geometry Generation
    {
      ArenaMarker scratch_guard(scratch_arena_a);
      SolidGenerator gen(solid_idx);
      PolyMesh base = gen.generate(scratch_arena_a, scratch_arena_a);

      MeshOps::compile_hankin(base, primary_hankin, geometry_arena,
                              scratch_arena_a);
      MeshOps::update_hankin(primary_hankin, primary_mesh, geometry_arena,
                             params.hankin_angle);

      // Initial Topology & Color
      primary_topology =
          MeshOps::classify_faces_by_topology(primary_mesh, scratch_arena_a);
      pick_palettes(primary_palettes);
    }

    // Start the Loop
    start_hankin_cycle();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  Orientation<W> orientation;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;

  void preallocate_buffers() {
    // Exact maximums needed for the largest Archimedean solid
    // (Truncated Icosidodecahedron)
    constexpr int MAX_V = 150;
    constexpr int MAX_F = 100;
    constexpr int MAX_I = 400;

    // The output sizes after the Hankin math is applied
    constexpr int OUT_V = 600;
    constexpr int OUT_F = 250;
    constexpr int OUT_I = 1600;

    // Primary Mesh
    primary_mesh.vertices.initialize(geometry_arena, OUT_V);
    primary_mesh.face_counts.initialize(geometry_arena, OUT_F);
    primary_mesh.face_offsets.initialize(geometry_arena, OUT_F);
    primary_mesh.faces.initialize(geometry_arena, OUT_I);

    // Secondary Mesh
    secondary_mesh.vertices.initialize(geometry_arena, OUT_V);
    secondary_mesh.face_counts.initialize(geometry_arena, OUT_F);
    secondary_mesh.face_offsets.initialize(geometry_arena, OUT_F);
    secondary_mesh.faces.initialize(geometry_arena, OUT_I);

    // Primary Hankin Compiler
    primary_hankin.baseVertices.initialize(geometry_arena, MAX_V);
    primary_hankin.staticVertices.initialize(geometry_arena, MAX_I);
    primary_hankin.dynamicVertices.initialize(geometry_arena, MAX_I);
    primary_hankin.dynamicInstructions.initialize(geometry_arena, MAX_I);
    primary_hankin.face_counts.initialize(geometry_arena, OUT_F);
    primary_hankin.faces.initialize(geometry_arena, OUT_I);

    // Secondary Hankin Compiler
    secondary_hankin.baseVertices.initialize(geometry_arena, MAX_V);
    secondary_hankin.staticVertices.initialize(geometry_arena, MAX_I);
    secondary_hankin.dynamicVertices.initialize(geometry_arena, MAX_I);
    secondary_hankin.dynamicInstructions.initialize(geometry_arena, MAX_I);
    secondary_hankin.face_counts.initialize(geometry_arena, OUT_F);
    secondary_hankin.faces.initialize(geometry_arena, OUT_I);
  }

  int solid_idx = 0;

  // Topology & Color State
  std::vector<int> primary_topology;
  std::vector<int> secondary_topology;
  std::vector<const Palette *> primary_palettes;
  std::vector<const Palette *> secondary_palettes;

  // Palette Pool (matching IslamicStars)
  const std::vector<const Palette *> source_palettes_pool = {
      &Palettes::embers, &Palettes::richSunset, &Palettes::brightSunrise,
      &Palettes::bruisedMoss, &Palettes::lavenderLake};

  void pick_palettes(std::vector<const Palette *> &target) {
    target = source_palettes_pool;
    static std::mt19937 g(12345 + (int)timeline.t);
    std::shuffle(target.begin(), target.end(), g);
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;

    // 1. Mutation: Animate Hankin Angle
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() { this->start_morph_cycle(); }));

    // 2. Sprite: Draw the mesh
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          MeshOps::update_hankin(primary_hankin, primary_mesh,
                                                 geometry_arena,
                                                 params.hankin_angle);
                          draw_topology_mesh(c, primary_mesh, primary_topology,
                                             primary_palettes, opacity);
                        },
                        DURATION));
  }

  void start_morph_cycle() {
    constexpr int DURATION = 16;

    // Identify Next Solid
    int next_idx = (solid_idx + 1) % Solids::Collections::num_simple_solids;

    const auto &current_entry = Solids::Collections::simple_solids[solid_idx];
    const auto &next_entry = Solids::Collections::simple_solids[next_idx];

    std::cout << "Morphing: " << current_entry.name << " -> " << next_entry.name
              << std::endl;

    // Generate Target Mesh (Secondary)
    {
      ArenaMarker scratch_guard(scratch_arena_a);
      SolidGenerator next_gen(next_idx);
      PolyMesh next_base = next_gen.generate(scratch_arena_a, scratch_arena_a);

      MeshOps::compile_hankin(next_base, secondary_hankin, geometry_arena,
                              scratch_arena_a);
      MeshOps::update_hankin(secondary_hankin, secondary_mesh, geometry_arena,
                             params.hankin_angle);

      // Force primary_mesh to exactly match the final angle before morph
      // capturing it!
      MeshOps::update_hankin(primary_hankin, primary_mesh, geometry_arena,
                             params.hankin_angle);

      // Prepare Secondary State
      secondary_topology =
          MeshOps::classify_faces_by_topology(secondary_mesh, scratch_arena_a);
      pick_palettes(secondary_palettes);
    }

    // 1. Morph Animation
    timeline.add(0, Animation::MeshMorph(&primary_mesh, &secondary_mesh,
                                         &morph_buffer, primary_mesh,
                                         secondary_mesh, DURATION, false,
                                         ease_in_out_sin));

    // 2. Sprite: Outgoing
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          draw_topology_mesh(c, primary_mesh, primary_topology,
                                             primary_palettes, opacity);
                        },
                        DURATION, 0, ease_mid, DURATION, ease_mid));

    // 3. Sprite: Incoming
    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) {
                 draw_topology_mesh(c, secondary_mesh, secondary_topology,
                                    secondary_palettes, opacity);
               },
               DURATION, DURATION, ease_mid, 0, ease_mid)
               .then([this, next_idx]() {
                 // Commit State
                 this->solid_idx = next_idx;

                 // Move Secondary -> Primary via pointer swap (safe due to
                 // pre-allocation)
                 std::swap(this->primary_topology, this->secondary_topology);
                 std::swap(this->primary_palettes, this->secondary_palettes);
                 std::swap(this->primary_hankin, this->secondary_hankin);
                 std::swap(this->primary_mesh, this->secondary_mesh);

                 // Save scratch offset to prevent memory leaks from temporary
                 // geometry
                 ArenaMarker scratch_guard(scratch_arena_a);

                 // Re-compile Hankin for the NEW solid (Generated purely in
                 // scratch memory)
                 SolidGenerator new_gen(solid_idx);
                 PolyMesh new_base =
                     new_gen.generate(scratch_arena_a, scratch_arena_a);

                 MeshOps::compile_hankin(new_base, primary_hankin,
                                         geometry_arena, scratch_arena_a);
                 MeshOps::update_hankin(primary_hankin, primary_mesh,
                                        geometry_arena, params.hankin_angle);

                 // Loop
                 this->start_hankin_cycle();
               }));
  }

  void draw_topology_mesh(Canvas &canvas, const MeshState &mesh,
                          const std::vector<int> &topology,
                          const std::vector<const Palette *> &palettes,
                          float opacity) {
    if (mesh.vertices.empty())
      return;
    if (opacity < 0.01f)
      return;

    // Create a pristine geometry copy on the scratch matrix buffer for dynamic
    // traversal
    ArenaMarker scratch_guard(scratch_arena_a);
    MeshState rotated_mesh;
    MeshOps::transform(mesh, rotated_mesh, scratch_arena_a);

    Quaternion q = orientation.get();
    for (auto &v : rotated_mesh.vertices) {
      v = rotate(v, q);
    }

    auto shader = [&](const Vector &p, Fragment &f) {
      int faceIdx = (int)std::round(f.v2);
      int topoIdx = 0;
      if (faceIdx >= 0 && faceIdx < (int)topology.size()) {
        topoIdx = topology[faceIdx];
      }

      // Safety: Ensure palettes has content
      if (palettes.empty()) {
        f.color = Color4(1, 0, 1, 1); // Error pink
        return;
      }

      const Palette *pal = palettes[topoIdx % palettes.size()];

      // Edge Distance Intensity (v1 is -dist)
      float distFromEdge = -f.v1;
      float size = f.size;
      float normalizedDist = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
      float t = hs::clamp(normalizedDist * params.intensity, 0.0f, 1.0f);

      f.color = pal->get(t);
      f.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader);
  }

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
  } params;
};
