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

    // Load initial shape into carousel's front slot
    load_shape(carousel.current(), compiled_hankin,
               palettes_slots[carousel.front_index()], solid_idx,
               params.hankin_angle);

    start_hankin_cycle();
  }

  bool show_bg() const override { return false; }
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  const std::array<ProceduralPalette, 5> source_palettes_pool = {
      Palettes::embers, Palettes::richSunset, Palettes::brightSunrise,
      Palettes::bruisedMoss, Palettes::lavenderLake};

  /// Load a shape: generate base, compile hankin, produce output mesh,
  /// classify topology, shuffle palettes.
  void load_shape(MeshState &out_mesh, CompiledHankin &out_hankin,
                  std::array<ProceduralPalette, 5> &out_palettes, int idx,
                  float angle) {
    auto solids = Solids::Collections::get_simple_solids();
    hs::log("Loading shape: '%s'", solids[idx].name);

    scratch_arena_a.reset();
    scratch_arena_b.reset();
    ScopedScratch _a(scratch_arena_a);
    ScopedScratch _b(scratch_arena_b);

    PolyMesh temp_base = Solids::finalize_solid(solids[idx].generate(scratch_arena_a, scratch_arena_b),
                                                scratch_arena_a);

    // Compile hankin instructions onto persistent arena
    out_hankin = CompiledHankin();
    MeshOps::compile_hankin(temp_base, out_hankin, persistent_arena,
                            scratch_arena_a);

    // Produce output mesh at the given angle
    out_mesh.clear();
    MeshOps::update_hankin(out_hankin, out_mesh, persistent_arena, angle);

    // Classify topology
    {
      scratch_arena_a.reset();
      scratch_arena_b.reset();
      MeshOps::classify_faces_by_topology(out_mesh, scratch_arena_a, scratch_arena_b, persistent_arena);
    }

    // Shuffle palettes
    out_palettes = source_palettes_pool;
    std::shuffle(out_palettes.begin(), out_palettes.end(), hs::random());
  }

  void draw_mesh(Canvas &canvas, const MeshState &mesh,
                 const ArenaVector<int> &topology,
                 const std::array<ProceduralPalette, 5> &palettes,
                 float opacity, const Quaternion &q) {
    if (mesh.vertices.empty() || opacity < 0.01f)
      return;

    ScopedScratch _(scratch_arena_a);
    MeshState rotated_mesh;
    MeshOps::transform(mesh, rotated_mesh, scratch_arena_a);

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

  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() { this->start_morph_cycle(); }));

    int front = carousel.front_index();
    timeline.add(
        0, Animation::Sprite(
               [this, front](Canvas &c, float opacity) {
                 MeshOps::update_hankin(compiled_hankin, carousel.slot(front),
                                        persistent_arena, params.hankin_angle);
                 draw_mesh(c, carousel.slot(front),
                           carousel.slot(front).topology, palettes_slots[front],
                           opacity, orientation.get());
               },
               DURATION));
  }

  FLASHMEM void start_morph_cycle() {
    constexpr int MORPH_FRAMES = 16;
    auto solids = Solids::Collections::get_simple_solids();
    int next_idx = (solid_idx + 1) % solids.size();

    int old_front = carousel.front_index();
    int new_slot = 1 - old_front;

    // Load incoming shape
    load_shape(carousel.slot(new_slot), compiled_hankin_staging,
               palettes_slots[new_slot], next_idx, params.hankin_angle);

    // MeshMorph needs SEPARATE active meshes (self-clone zeroes data)
    size_t max_v = std::max(carousel.slot(old_front).vertices.size(),
                            carousel.slot(new_slot).vertices.size());
    morph_buffer.preallocate(persistent_arena, max_v);

    morph_old_slot_ = old_front;
    morph_new_slot_ = new_slot;

    // Schedule self-rendering morph (MeshMorph draws both meshes with
    // crossfade)
    timeline.add(
        0,
        Animation::MeshMorph(
            &active_mesh_A, &active_mesh_B, &morph_buffer, &persistent_arena,
            carousel.slot(old_front), carousel.slot(new_slot),
            draw_morph_outgoing_fn_, draw_morph_incoming_fn_,
            MORPH_FRAMES, false, ease_in_out_sin)
            .then([this, next_idx, new_slot]() {
              solid_idx = next_idx;
              carousel.set_front(new_slot);

              // Promote staged hankin to primary
              compiled_hankin = std::move(compiled_hankin_staging);
              compiled_hankin_staging = CompiledHankin();

              // Drop temporary morph meshes and compact arena
              active_mesh_A = MeshState();
              active_mesh_B = MeshState();
              carousel.incoming() = MeshState();
              morph_buffer = Animation::MorphBuffer();
              {
                // Backup live data to BOTH scratch arenas (matching old split)
                ScratchScope sa(scratch_arena_a);
                ScratchScope sb(scratch_arena_b);
                CompiledHankin backup_hankin;
                MeshOps::clone(compiled_hankin, backup_hankin, sa.raw());
                MeshState backup_mesh;
                MeshOps::clone(carousel.current(), backup_mesh, sb.raw());
                // Wipe persistent and clone back
                persistent_arena.reset();
                carousel.current() = MeshState();
                compiled_hankin = CompiledHankin();
                MeshOps::clone(backup_mesh, carousel.current(), persistent_arena);
                MeshOps::clone(backup_hankin, compiled_hankin, persistent_arena);
              }

              // Ensure resting mesh state is correct
              MeshOps::update_hankin(compiled_hankin, carousel.current(),
                                     persistent_arena, params.hankin_angle);
              start_hankin_cycle();
            }));
  }

  MeshCarousel<W> carousel;
  CompiledHankin compiled_hankin;         // Active during hankin cycle
  CompiledHankin compiled_hankin_staging; // Built during morph cycle
  std::array<ProceduralPalette, 5> palettes_slots[2];

  // MeshMorph transient state
  MeshState active_mesh_A;
  MeshState active_mesh_B;
  Animation::MorphBuffer morph_buffer;
  int morph_old_slot_ = 0;
  int morph_new_slot_ = 1;

  /// Draw callbacks for morph — stored as members so FunctionRef doesn't dangle.
  void draw_morph_outgoing(Canvas &c, const MeshState &mesh, float opacity) {
    draw_mesh(c, mesh, carousel.slot(morph_old_slot_).topology,
              palettes_slots[morph_old_slot_], opacity, orientation.get());
  }
  void draw_morph_incoming(Canvas &c, const MeshState &mesh, float opacity) {
    draw_mesh(c, mesh, carousel.slot(morph_new_slot_).topology,
              palettes_slots[morph_new_slot_], opacity, orientation.get());
  }
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_outgoing_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_morph_outgoing(c, m, o);
      }};
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_incoming_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_morph_incoming(c, m, o);
      }};

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
