/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"

#include <algorithm>
#include <map>
#include <random>

#include "../solids.h"

template <int W, int H> class HankinSolids : public Effect {
public:
  FLASHMEM HankinSolids() : Effect(W, H), filters() {}

  void init() override {
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);
    registerParam("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    solid_idx = 0;

    primary.preallocate(Persistent(persistent_arena));
    secondary.preallocate(Persistent(persistent_arena));
    morph_buffer.preallocate(Persistent(persistent_arena), 544);

    {
      MemoryCtx ctx;
      ScopedScratch _a(ctx.get_scratch_front());
      ScopedScratch _b(ctx.get_scratch_back());
      primary.load(solid_idx, params.hankin_angle, Persistent(persistent_arena),
                   ctx, source_palettes_pool, timeline.t);
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
    MeshState mesh;        // The static resting geometry
    MeshState base;        // The simple base geometry
    MeshState active_mesh; // The dynamic animation puppet!
    std::array<ProceduralPalette, 5> palettes;

    void preallocate(Persistent arena) {
      // Measured maxima across all 18 simple solids:
      // Base: V=120 F=92 I=360  Output: V=540 F=182 I=1440
      constexpr int MAX_V = 128, MAX_F = 96, MAX_I = 368;
      constexpr int OUT_V = 544, OUT_F = 192, OUT_I = 1440;

      mesh.vertices.initialize(arena, OUT_V);
      mesh.face_counts.initialize(arena, OUT_F);
      if constexpr (requires { mesh.face_offsets; })
        mesh.face_offsets.initialize(arena, OUT_F);
      mesh.faces.initialize(arena, OUT_I);

      base.vertices.initialize(arena, MAX_V);
      base.face_counts.initialize(arena, MAX_F);
      if constexpr (requires { base.face_offsets; })
        base.face_offsets.initialize(arena, MAX_F);
      base.faces.initialize(arena, MAX_I);

      hankin.baseVertices.initialize(arena, MAX_V);
      hankin.staticVertices.initialize(arena, (MAX_I / 2) + 1);
      hankin.dynamicVertices.initialize(arena, MAX_I);
      hankin.dynamicInstructions.initialize(arena, MAX_I);
      hankin.face_counts.initialize(arena, OUT_F);
      hankin.faces.initialize(arena, OUT_I);

      active_mesh.vertices.initialize(arena, OUT_V);
      active_mesh.face_counts.initialize(arena, OUT_F);
      if constexpr (requires { active_mesh.face_offsets; })
        active_mesh.face_offsets.initialize(arena, OUT_F);
      active_mesh.faces.initialize(arena, OUT_I);
    }

    void load(int idx, float angle, Persistent geom, MemoryCtx &ctx,
              const std::array<ProceduralPalette, 5> &pool, float time) {
      auto solids = Solids::Collections::get_simple_solids();
      hs::log("Transitioning to '%s'", solids[idx].name);

      PolyMesh temp_base = Solids::finalize_solid(solids[idx].generate(ctx), ctx.get_scratch_front());
      MeshOps::compile(temp_base, base, geom);
      MeshOps::compile_hankin(temp_base, hankin, ctx.get_scratch_back(),
                              ctx.get_scratch_front());

      // Calculate initial state (angle = 0)
      MeshOps::update_hankin(hankin, mesh, geom, angle);
      ctx.swap_scratch();
      MeshOps::classify_faces_by_topology(mesh, ctx);

      palettes = pool;
      std::shuffle(palettes.begin(), palettes.end(), hs::random());
    }

    void swap(ShapeState &other) {
      std::swap(hankin, other.hankin);
      std::swap(mesh, other.mesh);
      std::swap(base, other.base);
      std::swap(active_mesh, other.active_mesh); // Swap the puppet too!
      std::swap(palettes, other.palettes);
    }
  };

  const std::array<ProceduralPalette, 5> source_palettes_pool = {
      Palettes::embers, Palettes::richSunset, Palettes::brightSunrise,
      Palettes::bruisedMoss, Palettes::lavenderLake};

  void draw_dynamic_morph(Canvas &c, float opacity,
                          const ArenaVector<int> &p_topo,
                          const std::array<ProceduralPalette, 5> &p_pal,
                          const ArenaVector<int> &s_topo,
                          const std::array<ProceduralPalette, 5> &s_pal) {
    MemoryCtx ctx;
    ScopedScratch _a(ctx.get_scratch_front());

    Quaternion q = orientation.get();

    float op_A = (1.0f - morph_alpha) * opacity;
    if (op_A > 0.01f) {
      draw_topology_mesh(c, primary.active_mesh, p_topo, p_pal, op_A, q);
    }

    float op_B = morph_alpha * opacity;
    if (op_B > 0.01f) {
      draw_topology_mesh(c, secondary.active_mesh, s_topo, s_pal, op_B, q);
    }
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() { this->start_morph_cycle(); }));

    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) {
                 MeshOps::update_hankin(primary.hankin, primary.mesh,
                                        persistent_arena, params.hankin_angle);
                 draw_topology_mesh(c, primary.mesh, primary.mesh.topology,
                                    primary.palettes, opacity);
               },
               DURATION));
  }

  FLASHMEM void start_morph_cycle() {
    constexpr int DURATION = 16;
    auto solids = Solids::Collections::get_simple_solids();
    int next_idx = (solid_idx + 1) % solids.size();

    {
      MemoryCtx ctx;
      ScopedScratch _a(ctx.get_scratch_front());
      ScopedScratch _b(ctx.get_scratch_back());
      // Load generates the fully compiled Hankin topology at angle = 0
      secondary.load(next_idx, params.hankin_angle,
                     Persistent(persistent_arena), ctx, source_palettes_pool,
                     timeline.t);
    }

    morph_alpha = 0.0f;
    timeline.add(
        0, Animation::Transition(morph_alpha, 1.0f, DURATION, ease_in_out_sin));

    // The Dual Elastic Stretch DIRECTLY on the encapsulated active meshes
    timeline.add(0, Animation::MeshMorph(
                        &primary.active_mesh, &secondary.active_mesh,
                        &morph_buffer, &persistent_arena, primary.mesh,
                        secondary.mesh, DURATION, false, ease_in_out_sin));

    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) {
                 draw_dynamic_morph(c, opacity, primary.mesh.topology,
                                    primary.palettes, secondary.mesh.topology,
                                    secondary.palettes);
               },
               DURATION)
               .then([this, next_idx]() {
                 this->solid_idx = next_idx;
                 primary.swap(secondary);

                 // Ensure resting state is perfect before resuming
                 // cycle
                 MeshOps::update_hankin(primary.hankin, primary.mesh,
                                        persistent_arena, params.hankin_angle);
                 this->start_hankin_cycle();
               }));
  }

  void draw_topology_mesh(Canvas &canvas, const MeshState &mesh,
                          const ArenaVector<int> &topology,
                          const std::array<ProceduralPalette, 5> &palettes,
                          float opacity, const Quaternion &q) {
    if (mesh.vertices.empty() || opacity < 0.01f)
      return;

    MemoryCtx ctx;
    ScopedScratch _(ctx.get_scratch_front());
    MeshState rotated_mesh;
    MeshOps::transform(mesh, rotated_mesh, ctx.get_scratch_front());

    for (auto &v : rotated_mesh.vertices)
      v = rotate(v, q);

    auto shader = [&](const Vector &p, Fragment &f) {
      int faceIdx = (int)std::round(f.v2);
      int topoIdx = (faceIdx >= 0 && faceIdx < (int)topology.size())
                        ? topology[faceIdx]
                        : 0;

      float distFromEdge = -f.v1;
      float normalizedDist =
          (f.size > 0.0001f) ? (distFromEdge / f.size) : 0.0f;
      float t = hs::clamp(normalizedDist * params.intensity, 0.0f, 1.0f);

      f.color = get_color(palettes[topoIdx % palettes.size()], t);
      f.color.alpha *= opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader,
                           params.debug_bb);
  }

  // Overload without pre-computed quaternion (for non-morph paths)
  void draw_topology_mesh(Canvas &canvas, const MeshState &mesh,
                          const ArenaVector<int> &topology,
                          const std::array<ProceduralPalette, 5> &palettes,
                          float opacity) {
    draw_topology_mesh(canvas, mesh, topology, palettes, opacity,
                       orientation.get());
  }

  ShapeState primary;
  ShapeState secondary;

  Animation::MorphBuffer morph_buffer;
  float morph_alpha = 0.0f;

  Orientation<W> orientation;
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  int solid_idx = 0;

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
    bool debug_bb = false;
  } params;
};
