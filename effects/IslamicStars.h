/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"
#include "../transformers.h"
#include <algorithm>
#include <map>
#include <random>

#include "../solids.h"

template <int W, int H> class IslamicStars : public Effect {

public:
  FLASHMEM IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {}

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - (128 + 128) * 1024, 128 * 1024,
                     128 * 1024);

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
  std::array<ProceduralPalette, 5> palettes_history[2];

  std::array<ProceduralPalette, 5> palettes = {
      Palettes::embers, Palettes::richSunset, Palettes::brightSunrise,
      Palettes::bruisedMoss, Palettes::lavenderLake};

  void ripple(Canvas &canvas) {
    Vector origin = random_vector();
    for (int i = 0; i < params.burst_size; i++) {
      ripple_gen.spawn(i * 16, origin, PI_F / ripple_duration,
                       static_cast<int>(ripple_duration));
    }
  }

  void draw_shape(Canvas &canvas, float opacity, const MeshState &base_state,
                  const ArenaVector<int> &faceIndices,
                  const std::array<ProceduralPalette, 5> &palettes) {
    // Early exit if barely visible to save an entire pass
    if (opacity <= 0.005f)
      return;

    ScopedScratch _a(scratch_arena_a);
    MeshState transformed_state;
    OrientTransformer<W> camera(orientation);
    MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                       ripple_gen, camera);

    const int *raw_indices = faceIndices.data();

    auto fragment_shader = [&](const Vector &p, Fragment &frag) {
      // Strip safety bounds checking; trust the rasterizer's face index output
      int faceIdx = static_cast<int>(frag.v2);
      const ProceduralPalette &pal = palettes[raw_indices[faceIdx]];

      float size = frag.size;
      // Pre-negate v1 and simplify
      float intensity = (size > 0.0001f) ? (-frag.v1 / size) : 0.0f;
      intensity = hs::clamp(intensity, 0.0f, 1.0f);

      frag.color = get_color(pal, intensity);
      frag.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           params.debug_bb);
  }

  void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();

    // Capture the slot index for this shape's draw lambda
    int capture_idx = 1 - carousel.front_index();

    // Prepare palettes for the incoming slot
    static std::mt19937 g(12345 + (int)timeline.t);
    palettes_history[capture_idx] = palettes;
    std::shuffle(palettes_history[capture_idx].begin(),
                 palettes_history[capture_idx].end(), g);

    int idx = solid_idx; // capture for generate lambda

    auto draw_fn = [this, capture_idx](Canvas &canvas, float opacity) {
      const MeshState &mesh = carousel.slot(capture_idx);
      // Map topology to palette indices
      this->draw_shape(canvas, opacity, mesh, mesh.topology,
                       palettes_history[capture_idx]);
    };

    carousel.transition(
        timeline,
        // generate_fn
        [idx, &solids](Arena &a, Arena &b) { return solids[idx].generate(a, b); },
        // draw_outgoing, draw_incoming
        draw_fn, draw_fn,
        160, // duration
        32,  // fade_in
        32   // fade_out
    );

    // Map topology -> palette index on the newly populated slot
    int num_colors = static_cast<int>(palettes_history[capture_idx].size());
    // After transition(), front has flipped, so capture_idx is now front
    ArenaVector<int> &faceIndices = carousel.slot(capture_idx).topology;
    for (size_t i = 0; i < faceIndices.size(); ++i) {
      faceIndices[i] = faceIndices[i] % num_colors;
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
