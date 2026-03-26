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
    configure_arenas(GLOBAL_ARENA_SIZE - (128 + 128) * 1024, 128 * 1024,
                     128 * 1024);

    for (int i = 0; i < NUM_PALETTES; ++i)
      baked_palette_bank_.entries[i].bake(persistent_arena,
                                          *source_palettes[i]);

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Ripp Amp", &ripple_gen.params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.params.decay, 0.0f, 5.0f);
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 300.0f);
    registerParam("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    ripple_gen.params.amplitude = 0.4f;
    ripple_gen.params.thickness = 0.7f;
    ripple_gen.params.decay = 0.1f;

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
  Timeline<W> timeline;
  Pipeline<W, H> filters;
  RippleTransformer<W, 8> ripple_gen;
  FastNoiseLite noise;
  float ripple_duration = 80.0f;
  int solid_idx = -1;
  MeshCarousel<W> carousel;

  static constexpr int NUM_PALETTES = 5;
  static constexpr std::array<const ProceduralPalette *, NUM_PALETTES>
      source_palettes = {&Palettes::embers, &Palettes::richSunset,
                         &Palettes::brightSunrise, &Palettes::bruisedMoss,
                         &Palettes::lavenderLake};
  BakedPaletteBank baked_palette_bank_;
  std::array<int, NUM_PALETTES> palettes_slots[2];

  void ripple(Canvas &canvas) {
    Vector origin = random_vector();
    for (int i = 0; i < params.burst_size; i++) {
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

    auto fragment_shader = [&](const Vector &p, Fragment &frag) {
      int faceIdx = static_cast<int>(frag.v2);
      int topoIdx = raw_indices[faceIdx];

      float size = frag.size;
      float intensity = (size > 0.0001f) ? (-frag.v1 / size) : 0.0f;
      intensity = hs::clamp(intensity, 0.0f, 1.0f);

      frag.color =
          baked_palette_bank_.entries[palette_idx[topoIdx % NUM_PALETTES]].get(
              intensity);
      frag.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  void shuffle_palette_indices(std::array<int, NUM_PALETTES> &out) {
    for (int i = 0; i < NUM_PALETTES; ++i)
      out[i] = i;
    std::shuffle(out.begin(), out.end(), hs::random());
  }

  void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int capture_idx = 1 - carousel.front_index();
    shuffle_palette_indices(palettes_slots[capture_idx]);

    int idx = solid_idx; // capture for generate lambda

    auto draw_fn = [this, capture_idx](Canvas &canvas, float opacity) {
      const MeshState &mesh = carousel.slot(capture_idx);
      this->draw_shape(canvas, opacity, mesh, mesh.topology,
                       palettes_slots[capture_idx]);
    };

    int back = 1 - carousel.front_index();

    // 1. Clear back slot and perform custom compaction
    {
      carousel.slot(back) = MeshState();
      Persist<MeshState> p(carousel.slot(carousel.front_index()),
                           scratch_arena_b, persistent_arena);
      Persist<BakedPaletteBank> pp(baked_palette_bank_, scratch_arena_b,
                                   persistent_arena);
      persistent_arena.reset();
      hs::log("IslamicStars: Finished arena compaction");
    }

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

    timeline.add(0, Animation::Sprite(draw_fn,
                                      160,          // duration
                                      32, ease_mid, // fade_in
                                      32, ease_mid  // fade_out
                                      ));

    // After transition(), front has flipped, so capture_idx is now front
    ArenaVector<int> &faceIndices = carousel.slot(capture_idx).topology;
    for (size_t i = 0; i < faceIndices.size(); ++i) {
      faceIndices[i] = faceIndices[i] % NUM_PALETTES;
    }

    // Log
    const auto &entry = solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)carousel.current().vertices.size(),
            (int)carousel.current().faces.size());

    // Schedule next overlapping
    int next_delay = 160 - 32; // duration - fade_out
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  struct Params {
    float duration = 96.0f;
    float fade = 32.0f;
    int burst_size = 4;
    bool debug_bb = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
