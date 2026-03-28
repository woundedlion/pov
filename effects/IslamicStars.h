/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include <algorithm>

// Uses HS_OS_CYCLES() from platform.h


template <int W, int H> class IslamicStars : public Effect {

public:
  FLASHMEM IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {}

  void init() override {
    // scratch arenas + room for BakedPaletteBank (~15KB) in persistent
    configure_arenas(GLOBAL_ARENA_SIZE - (120 + 120) * 1024, 120 * 1024,
                     120 * 1024);

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
    c_transform = 0;
    c_raster = 0;
    c_face_count = 0;
    hs::g_scan_metrics.reset();
    uint32_t t_start = HS_OS_CYCLES();

    Canvas canvas(*this);
    
    uint32_t t0 = HS_OS_CYCLES();
    ripple_gen.prepare_frame();
    c_ripple_prep = HS_OS_CYCLES() - t0;

    uint32_t t1 = HS_OS_CYCLES();
    timeline.step(canvas);
    c_timeline = HS_OS_CYCLES() - t1;

    static int frame_count = 0;
    if (frame_count++ % 60 == 0) {
      uint32_t t_total = HS_OS_CYCLES() - t_start;
      auto &m = hs::g_scan_metrics;
      uint32_t total_pixels = m.pixels_tested;
      uint32_t total_sdf = m.lut_hits + m.exact_hits;
      hs::log("=== IslamicStars Frame %d ===", frame_count);
      hs::log("  Total: %luus  Timeline: %luus", t_total/600, c_timeline/600);
      hs::log("  Phases: RipplePrep=%luus  Transform=%luus  Raster=%luus  Other=%luus",
              c_ripple_prep/600, c_transform/600, c_raster/600,
              (c_timeline - c_transform - c_raster)/600);
      hs::log("  Raster sub-phases (faces=%lu):", c_face_count);
      hs::log("    face_setup=%luus  scan_loop=%luus",
              m.face_setup/600, m.scan_loop/600);
      hs::log("    sdf_dist=%luus  frag_shader=%luus  plot=%luus",
              m.sdf_dist/600, m.frag_shader/600, m.plot/600);
      hs::log("  Pixels: tested=%lu culled=%lu (%.0f%%)",
              total_pixels, m.pixels_culled,
              total_pixels > 0 ? 100.0f * m.pixels_culled / total_pixels : 0.0f);
      hs::log("  SDF: lut_hits=%lu exact=%lu (lut%%=%.0f%%)",
              m.lut_hits, m.exact_hits,
              total_sdf > 0 ? 100.0f * m.lut_hits / total_sdf : 0.0f);
      if (c_face_count > 0) {
        hs::log("  Per-face avg: setup=%luus  scan=%luus  total=%luus",
                m.face_setup/600/c_face_count, m.scan_loop/600/c_face_count,
                c_raster/600/c_face_count);
      }
    }
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

  uint32_t c_transform = 0;
  uint32_t c_raster = 0;
  uint32_t c_ripple_prep = 0;
  uint32_t c_timeline = 0;
  uint32_t c_face_count = 0;

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
    uint32_t t0 = HS_OS_CYCLES();
    MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                       ripple_gen, camera);
    c_transform += (HS_OS_CYCLES() - t0);

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

    c_face_count += transformed_state.get_face_counts_size();
    uint32_t t1 = HS_OS_CYCLES();
    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           scratch_arena_a, params.debug_bb);
    c_raster += (HS_OS_CYCLES() - t1);
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
      persistent_arena.reset();

      // Rebake palettes into the fresh arena instead of tracking them 
      // during evacuation (saves fragmentation/OOM)
      for (int i = 0; i < NUM_PALETTES; ++i)
        baked_palette_bank_.entries[i].bake(persistent_arena,
                                            *source_palettes[i]);

      hs::log("IslamicStars: Finished arena compaction & rebake");
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
