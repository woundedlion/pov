/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include <algorithm>

template <int W, int H> class IslamicStars : public Effect {

public:
  FLASHMEM IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {}

  void init() override {
    // scratch arenas + room for BakedPaletteBank (~15KB) in persistent
    configure_arenas(GLOBAL_ARENA_SIZE - (120 + 120) * 1024, 120 * 1024,
                     120 * 1024);

    palette_bank_.bake_all(persistent_arena);

    // Set the ripple defaults BEFORE registering their pointers so the captured
    // registered defaults match the actual runtime values (registerParam snaps
    // *ptr as the default; setting them afterwards desynced the two).
    ripple_gen.params.amplitude = 0.4f;
    ripple_gen.params.thickness = 0.7f;
    ripple_gen.params.decay = 0.1f;

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Fade", &params.fade, 0.0f, 96.0f);
    registerParam("Burst", &params.burst_size, 1.0f, 8.0f);
    registerParam("Ripp Amp", &ripple_gen.params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.params.decay, 0.0f, 5.0f);
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 300.0f);
    registerParam("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    // Ripple now and schedule more ripples
    timeline.add(0, Animation::PeriodicTimer(
                        0, [this](auto &canvas) { ripple(canvas); }, false));
    timeline.add(0, Animation::PeriodicTimer(
                        96, [this](auto &canvas) { ripple(canvas); }, true));
    spawn_shape();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    ripple_gen.prepare_frame();
    timeline.step(canvas);
  }

private:
  Orientation<W> orientation;
  Timeline timeline;
  Pipeline<W, H> filters;
  RippleTransformer<8> ripple_gen;
  FastNoiseLite noise;
  float ripple_duration = 80.0f;
  int solid_idx = -1;
  MeshCarousel carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;
  std::array<int, NUM_PALETTES> palettes_slots[2];

  void ripple(Canvas &canvas) {
    Vector origin = random_vector();
    for (int i = 0; i < (int)params.burst_size; i++) {
      ripple_gen.spawn(i * 16, origin, PI_F / ripple_duration,
                       static_cast<int>(ripple_duration));
    }
  }

  void draw_shape(Canvas &canvas, float opacity, const MeshState &base_state,
                  const ArenaVector<int> &faceIndices,
                  const std::array<int, NUM_PALETTES> &palette_idx) {
    if (opacity <= 0.005f)
      return;
    ScratchScope _a(scratch_arena_a);
    MeshState transformed_state;
    OrientTransformer<W> camera(orientation);
    MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                       ripple_gen, camera);

    const int *raw_indices = faceIndices.data();

    const int num_faces = static_cast<int>(faceIndices.size());
    auto fragment_shader = [&](const Vector &p, Fragment &frag) {
      int faceIdx = static_cast<int>(frag.v2);
      // Defensive bounds guard matching HankinSolids: fall back to face 0 on a
      // malformed face index rather than reading raw_indices out of bounds.
      int topoIdx =
          (faceIdx >= 0 && faceIdx < num_faces) ? raw_indices[faceIdx] : 0;

      float intensity = hs::clamp(fragment_edge_dist(frag), 0.0f, 1.0f);

      frag.color =
          palette_bank_[palette_idx[topoIdx % NUM_PALETTES]].get(intensity);
      frag.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int capture_idx = 1 - carousel.front_index();
    MeshPaletteBank::shuffle_indices(palettes_slots[capture_idx]);

    int idx = solid_idx; // capture for generate lambda

    auto draw_fn = [this, capture_idx](Canvas &canvas, float opacity) {
      const MeshState &mesh = carousel.slot(capture_idx);
      this->draw_shape(canvas, opacity, mesh, mesh.topology,
                       palettes_slots[capture_idx]);
    };

    int back = 1 - carousel.front_index();

    // 1. Free the back slot and compact, rebaking palettes into the fresh
    // arena instead of tracking them through the evacuation (saves
    // fragmentation/OOM).
    carousel.compact_keep_front(
        [this](Arena &arena) { palette_bank_.bake_all(arena); });

    // 2. Generate new shape
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh mesh = solids[idx].generate(a, b);
      carousel.slot(back).clear();
      MeshOps::compile(mesh, carousel.slot(back), target);
    });

    scratch_arena_a.reset();
    scratch_arena_b.reset();
    MeshOps::classify_faces_by_topology(carousel.slot(back), scratch_arena_a,
                                        scratch_arena_b, persistent_arena);

    // 3. Flip front eagerly for the overlapping sprite
    carousel.set_front(back);

    // Live-read the Duration slider per shape cycle. Clamp fade to dur/2 so
    // the fade-in/out windows never overlap and the respawn delay (dur - fade)
    // stays >= dur/2 — a larger fade piles up sprites and per-frame respawns.
    int dur = static_cast<int>(params.duration);
    int fade = std::min(static_cast<int>(params.fade), dur / 2);
    timeline.add(0, Animation::Sprite(draw_fn, dur, fade, ease_mid, fade,
                                      ease_mid));

    // Topology indices are NOT pre-reduced here: draw_shape applies
    // `topoIdx % NUM_PALETTES` at render time, so reducing the stored topology
    // in place would be a redundant double-modulo (and a surprising mutation of
    // the carousel's mesh state). HankinSolids relies on the same render-time
    // modulo.

    // Log
    const auto &entry = solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)carousel.current().vertices.size(),
            (int)carousel.current().faces.size());

    // Schedule next overlapping
    int next_delay = dur - fade; // duration - fade_out
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  struct Params {
    float duration = 160.0f; // shape display period, in frames
    float fade = 32.0f;
    float burst_size = 4.0f; // ripples per burst; float-backed for registerParam
    bool debug_bb = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
