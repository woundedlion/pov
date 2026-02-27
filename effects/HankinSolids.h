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
  HankinSolids() : Effect(W, H), filters() {
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);

    persist_pixels = false;

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS,
                        Animation::RandomWalk<W>::Options::Languid()));

    solid_idx = 0;
    preallocate_buffers();

    {
      ArenaMarker _a(scratch_arena_a);
      ArenaMarker _b(scratch_arena_b);
      ScratchContext ctx(scratch_arena_a, scratch_arena_b);
      primary.load(solid_idx, params.hankin_angle, geometry_arena, ctx,
                   source_palettes_pool, timeline.t);
    }

    start_hankin_cycle();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  struct ShapeState {
    CompiledHankin hankin;
    MeshState mesh; // Post-Hankin (for static display)
    MeshState base; // Pre-Hankin (clean geometry for the morpher)
    std::vector<int> topology;
    std::vector<int> base_topology;
    std::vector<const Palette *> palettes;

    void preallocate(Arena &arena) {
      constexpr int MAX_V = 150, MAX_F = 100, MAX_I = 400;
      constexpr int OUT_V = 600, OUT_F = 250, OUT_I = 1600;

      mesh.vertices.initialize(arena, OUT_V);
      mesh.face_counts.initialize(arena, OUT_F);
      mesh.face_offsets.initialize(arena, OUT_F);
      mesh.faces.initialize(arena, OUT_I);

      base.vertices.initialize(arena, MAX_V);
      base.face_counts.initialize(arena, MAX_F);
      base.face_offsets.initialize(arena, MAX_F);
      base.faces.initialize(arena, MAX_I);

      hankin.baseVertices.initialize(arena, MAX_V);
      hankin.staticVertices.initialize(arena, MAX_I);
      hankin.dynamicVertices.initialize(arena, MAX_I);
      hankin.dynamicInstructions.initialize(arena, MAX_I);
      hankin.face_counts.initialize(arena, OUT_F);
      hankin.faces.initialize(arena, OUT_I);
    }

    void load(int idx, float angle, Arena &geom, ScratchContext &ctx,
              const std::vector<const Palette *> &pool, float time) {
      SolidGenerator gen(idx);
      hs::log("Transitioning to '%s'", Solids::get_entry(idx).name);

      PolyMesh temp_base = gen.generate(*(ctx.source), ctx);

      // Populate the 3 core mesh structures
      MeshOps::compile(temp_base, base, geom);
      MeshOps::compile_hankin(temp_base, hankin, ctx);
      MeshOps::update_hankin(hankin, mesh, geom, angle);

      topology = MeshOps::classify_faces_by_topology(mesh, *(ctx.source));
      base_topology = MeshOps::classify_faces_by_topology(base, *(ctx.source));

      // Shuffle palettes
      palettes = pool;
      std::mt19937 g(12345 + (int)time);
      std::shuffle(palettes.begin(), palettes.end(), g);
    }

    void swap(ShapeState &other) {
      std::swap(hankin, other.hankin);
      std::swap(mesh, other.mesh);
      std::swap(base, other.base);
      std::swap(topology, other.topology);
      std::swap(base_topology, other.base_topology);
      std::swap(palettes, other.palettes);
    }
  };

  void preallocate_buffers() {
    primary.preallocate(geometry_arena);
    secondary.preallocate(geometry_arena);

    constexpr int MAX_V = 150;
    constexpr int MAX_F = 100;
    constexpr int MAX_I = 400;

    active_base.vertices.initialize(geometry_arena, MAX_V);
    active_base.face_counts.initialize(geometry_arena, MAX_F);
    active_base.face_offsets.initialize(geometry_arena, MAX_F);
    active_base.faces.initialize(geometry_arena, MAX_I);
  }

  // Palette Pool (matching IslamicStars)
  const std::vector<const Palette *> source_palettes_pool = {
      &Palettes::embers, &Palettes::richSunset, &Palettes::brightSunrise,
      &Palettes::bruisedMoss, &Palettes::lavenderLake};

  // Helper to convert flat MeshState back to dynamic PolyMesh
  PolyMesh to_polymesh(const MeshState &ms, Arena &target) {
    PolyMesh pm;
    pm.vertices.initialize(target, ms.vertices.size());
    for (const auto &v : ms.vertices)
      pm.vertices.push_back(v);
    pm.face_counts.initialize(target, ms.face_counts.size());
    for (auto f : ms.face_counts)
      pm.face_counts.push_back(f);
    pm.faces.initialize(target, ms.faces.size());
    for (auto idx : ms.faces)
      pm.faces.push_back(idx);
    return pm;
  }

  void draw_dynamic_morph(Canvas &c, float opacity) {
    ArenaMarker _a(scratch_arena_a);
    ArenaMarker _b(scratch_arena_b);
    ScratchContext ctx(scratch_arena_a, scratch_arena_b);

    // 1. Convert active_base to PolyMesh for the compiler
    PolyMesh temp_pm = to_polymesh(active_base, scratch_arena_a);

    // 2. Compile and Apply Hankin dynamically to the stretching geometry
    CompiledHankin temp_hankin;
    MeshState temp_mesh;
    MeshOps::compile_hankin(temp_pm, temp_hankin, ctx);
    MeshOps::update_hankin(temp_hankin, temp_mesh, scratch_arena_a,
                           params.hankin_angle);

    // 3. Use the cached resting topology colors!
    // Do NOT recalculate hashes on the stretching faces.
    const std::vector<int> &active_topo = morph_buffer.using_dest_topology
                                              ? secondary.base_topology
                                              : primary.base_topology;

    draw_topology_mesh(c, temp_mesh, active_topo, primary.palettes,
                       secondary.palettes, morph_alpha, opacity);
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() { this->start_morph_cycle(); }));

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          MeshOps::update_hankin(primary.hankin, primary.mesh,
                                                 geometry_arena,
                                                 params.hankin_angle);
                          draw_topology_mesh(c, primary.mesh, primary.topology,
                                             primary.palettes, primary.palettes,
                                             0.0f, opacity);
                        },
                        DURATION));
  }

  void start_morph_cycle() {
    constexpr int DURATION = 16;
    int next_idx = (solid_idx + 1) % Solids::Collections::num_simple_solids;

    { // Load the next shape into the secondary buffer safely
      ArenaMarker _a(scratch_arena_a);
      ArenaMarker _b(scratch_arena_b);
      ScratchContext ctx(scratch_arena_a, scratch_arena_b);
      secondary.load(next_idx, params.hankin_angle, geometry_arena, ctx,
                     source_palettes_pool, timeline.t);
    }

    morph_alpha = 0.0f;

    // 1. Color Crossfade
    timeline.add(
        0, Animation::Transition(morph_alpha, 1.0f, DURATION, ease_in_out_sin));

    // 2. Slerp the Base Geometry
    timeline.add(0, Animation::MeshMorph(&active_base, &morph_buffer,
                                         &geometry_arena, primary.base,
                                         secondary.base, DURATION, false,
                                         ease_in_out_sin));

    // 3. Render the Dynamic Hankin Stars
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          draw_dynamic_morph(c, opacity);
                        },
                        DURATION)
                        .then([this, next_idx]() {
                          // Clean Swap and Restart!
                          this->solid_idx = next_idx;
                          primary.swap(secondary);
                          MeshOps::update_hankin(primary.hankin, primary.mesh,
                                                 geometry_arena,
                                                 params.hankin_angle);
                          this->start_hankin_cycle();
                        }));
  }

  void draw_topology_mesh(Canvas &canvas, const MeshState &mesh,
                          const std::vector<int> &topology,
                          const std::vector<const Palette *> &palettes_a,
                          const std::vector<const Palette *> &palettes_b,
                          float color_mix, float opacity,
                          bool is_morph = false) {
    if (mesh.vertices.empty())
      return;
    if (opacity < 0.01f)
      return;

    ArenaMarker _(scratch_arena_a);
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

      // Edge Distance Intensity (v1 is -dist)
      float distFromEdge = -f.v1;
      float size = f.size;
      float normalizedDist = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
      float t = hs::clamp(normalizedDist * params.intensity, 0.0f, 1.0f);

      Color4 c_a, c_b;

      if (is_morph) {
        if (morph_buffer.using_dest_topology) { // Growing (A -> B, but B has
                                                // more faces)
          c_a = palettes_a[topoIdx % palettes_a.size()]->get(t);
          c_b = palettes_b[topoIdx % palettes_b.size()]->get(t);
        } else { // Shrinking (A -> B, but A has more faces)
          c_a = palettes_a[topoIdx % palettes_a.size()]->get(t);
          c_b = palettes_b[topoIdx % palettes_b.size()]->get(t);
        }
      } else {
        c_a = palettes_a[topoIdx % palettes_a.size()]->get(t);
        c_b = palettes_b[topoIdx % palettes_b.size()]->get(t);
      }

      f.color = c_a.lerp(c_b, color_mix);
      f.color.alpha *= opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader);
  }

  // --- Encapsulated State ---
  ShapeState primary;
  ShapeState secondary;

  MeshState active_base;
  Animation::MorphBuffer morph_buffer;
  float morph_alpha = 0.0f;

  Orientation<W> orientation;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  int solid_idx = 0;

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
  } params;
};
