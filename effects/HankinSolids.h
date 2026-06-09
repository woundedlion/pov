/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

#include <algorithm>

template <int W, int H> class HankinSolids : public Effect {
public:
  FLASHMEM HankinSolids() : Effect(W, H), filters() {}

  void init() override {
    // scratch_b must hold CompiledHankin (~10KB) + BakedPaletteBank (~15KB)
    configure_arenas(GLOBAL_ARENA_SIZE - 16 * 1024 - 32 * 1024, 16 * 1024,
                     32 * 1024);
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);
    registerParam("Debug BB", &params.debug_bb);
    // Angle is swept by the Mutation in start_hankin_cycle(); flag it so the GUI
    // auto-pauses the sweep (and its morph cycle) when the user grabs the slider.
    markAnimated("Angle");

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    solid_idx = 0;

    palette_bank_.bake_all(persistent_arena);

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
  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;


  PolyMesh generate_base_solid(int idx, Arena &a, Arena &b) {
    auto solids = Solids::Collections::get_simple_solids();
    hs::log("Loading shape: '%s'", solids[idx].name);
    return Solids::finalize_solid(solids[idx].generate(a, b), a);
  }

  void classify_mesh_topology(MeshState &mesh) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
    MeshOps::classify_faces_by_topology(mesh, scratch_arena_a, scratch_arena_b,
                                        persistent_arena);
  }

  /// Full pipeline: generate → compile → evaluate → classify → palette indices.
  void load_shape(MeshState &out_mesh, CompiledHankin &out_hankin,
                  std::array<int, NUM_PALETTES> &out_palette_idx, int idx,
                  float angle) {
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = generate_base_solid(idx, a, b);

      out_hankin = CompiledHankin();
      MeshOps::compile_hankin(base, out_hankin, target, a);

      out_mesh.clear();
      MeshOps::update_hankin(out_hankin, out_mesh, target, angle);
    });

    classify_mesh_topology(out_mesh);
    MeshPaletteBank::shuffle_indices(out_palette_idx);
  }

  void draw_mesh(Canvas &canvas, const MeshState &mesh,
                 const ArenaVector<int> &topology,
                 const std::array<int, NUM_PALETTES> &palette_idx,
                 float opacity, const Quaternion &q) {
    if (mesh.vertices.empty() || opacity < 0.01f)
      return;

    ScratchScope _(scratch_arena_a);
    MeshState rotated_mesh;
    
    MeshOps::transform(mesh, rotated_mesh, scratch_arena_a);

    for (auto &v : rotated_mesh.vertices)
      v = rotate(v, q);

    auto fragment_shader = [&](const Vector &p, Fragment &f) {
      int faceIdx = (int)std::round(f.v2);
      int topoIdx = (faceIdx >= 0 && faceIdx < (int)topology.size())
                        ? topology[faceIdx]
                        : 0;

      float t = hs::clamp(fragment_edge_dist(f) * params.intensity, 0.0f, 1.0f);

      f.color = palette_bank_[palette_idx[topoIdx % NUM_PALETTES]].get(t);
      f.color.alpha *= opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  void promote_staged_hankin() {
    compiled_hankin = std::move(compiled_hankin_staging);
    compiled_hankin_staging = CompiledHankin();
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false, &anims_paused_)
                        .then([this]() {
                          this->start_morph_cycle();
                        }));

    int front = carousel.front_index();
    // Gate the render sprite on the same pause flag (holds the frame instead of
    // expiring) so a paused angle keeps rendering rather than going blank.
    timeline.add(
        0, Animation::Sprite(
               [this, front](Canvas &c, float opacity) {
                 MeshOps::update_hankin(compiled_hankin, carousel.slot(front),
                                        persistent_arena, params.hankin_angle);
                 draw_mesh(c, carousel.slot(front),
                           carousel.slot(front).topology, palettes_slots[front],
                           opacity, orientation.get());
               },
               DURATION, 0, ease_mid, 0, ease_mid, &anims_paused_));
  }

  FLASHMEM void start_morph_cycle() {
    constexpr int MORPH_FRAMES = 16;
    auto solids = Solids::Collections::get_simple_solids();
    int next_idx = (solid_idx + 1) % solids.size();

    int old_front = carousel.front_index();
    int new_slot = 1 - old_front;

    load_shape(carousel.slot(new_slot), compiled_hankin_staging,
               palettes_slots[new_slot], next_idx, params.hankin_angle);

    // Set slot indices before creating MeshMorph (draw callbacks reference
    // these)
    morph_old_slot_ = old_front;
    morph_new_slot_ = new_slot;

    timeline.add(
        0, Animation::MeshMorph(
               carousel.slot(old_front), carousel.slot(new_slot),
               persistent_arena, draw_morph_outgoing_fn_,
               draw_morph_incoming_fn_, MORPH_FRAMES, ease_in_out_sin)
               .then([this, next_idx, new_slot]() {
                 solid_idx = next_idx;
                 carousel.set_front(new_slot);
                 promote_staged_hankin();
                 // Manual compaction: preserve both carousel slots +
                 // compiled_hankin
                 {
                   Persist<CompiledHankin> ph(compiled_hankin, scratch_arena_b,
                                              persistent_arena);
                   Persist<MeshState> p0(carousel.slot(0), scratch_arena_a,
                                         persistent_arena);
                   Persist<MeshState> p1(carousel.slot(1), scratch_arena_a,
                                         persistent_arena);
                   Persist<MeshPaletteBank> pp(palette_bank_, scratch_arena_b,
                                               persistent_arena);
                   persistent_arena.reset();
                   hs::log("morph_cycle_then: finished arena compaction");
                 }

                 MeshOps::update_hankin(compiled_hankin, carousel.current(),
                                        persistent_arena, params.hankin_angle);
                 start_hankin_cycle();
               }));
  }

  MeshCarousel<W> carousel;
  CompiledHankin compiled_hankin;         // Active during hankin cycle
  CompiledHankin compiled_hankin_staging; // Built during morph cycle
  std::array<int, NUM_PALETTES> palettes_slots[2];

  // Slot indices for morph draw callbacks (set before each morph)
  int morph_old_slot_ = 0;
  int morph_new_slot_ = 1;

  /// Morph draw callbacks — members for stable FunctionRef lifetime
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_outgoing_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_old_slot_).topology,
                  palettes_slots[morph_old_slot_], o, orientation.get());
      }};
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_incoming_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_new_slot_).topology,
                  palettes_slots[morph_new_slot_], o, orientation.get());
      }};

  Orientation<W> orientation;
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H> filters;
  int solid_idx = 0;

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
    bool debug_bb = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(HankinSolids)
