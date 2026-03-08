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
    PersistentTracker::register_mesh(&mesh_states[0]);
    PersistentTracker::register_mesh(&mesh_states[1]);

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Ripp Amp", &ripple_gen.params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.params.decay, 0.0f, 5.0f);
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 300.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP));

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
  float ripple_duration = 80.0f;
  int solid_idx = -1;
  MeshState mesh_states[2];
  std::array<ProceduralPalette, 5> palettes_history[2];
  int current_mesh_idx = 0;

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

    MemoryCtx ctx;
    ScopedScratch _a(ctx.get_scratch_front());
    MeshState transformed_state;
    OrientTransformer<W> camera(orientation);
    MeshOps::transform(base_state, transformed_state, ctx.get_scratch_front(),
                       ripple_gen, camera);

    // Hoist raw pointer out of the lambda to avoid arena overhead per-fragment
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

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader);
  }

  void spawn_shape() {
    MemoryCtx ctx;
    solid_idx = (solid_idx + 1) % Solids::Collections::num_islamic_solids;
    current_mesh_idx = (current_mesh_idx + 1) % 2;

    // Forced compaction
    mesh_states[current_mesh_idx] = MeshState();
    PersistentTracker::auto_compact(ctx.get_scratch_back());

    // Generate new shape
    {
      ScopedScratch _(ctx.get_scratch_front());
      SolidGenerator gen(solid_idx + Solids::Collections::num_simple_solids);
      PolyMesh local_mesh = gen.generate(ctx.get_scratch_front(), ctx);
      ctx.update_persistent(mesh_states[current_mesh_idx], local_mesh);
    }
    MeshOps::classify_faces_by_topology(mesh_states[current_mesh_idx], ctx);

    // Log transition
    const auto &entry = Solids::Collections::islamic_solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)mesh_states[current_mesh_idx].vertices.size(),
            (int)mesh_states[current_mesh_idx].faces.size());

    // Prepare Palettes
    static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
    palettes_history[current_mesh_idx] = palettes;
    std::shuffle(palettes_history[current_mesh_idx].begin(),
                 palettes_history[current_mesh_idx].end(), g);

    // Map Topology to palette
    int num_colors =
        static_cast<int>(palettes_history[current_mesh_idx].size());
    ArenaVector<int> &faceIndices = mesh_states[current_mesh_idx].topology;
    for (size_t i = 0; i < faceIndices.size(); ++i) {
      faceIndices[i] = faceIndices[i] % num_colors;
    }

    // Create Sprite
    int duration = 160;
    int fade_in = 32;
    int fade_out = 32;
    int capture_idx = current_mesh_idx;
    auto draw_fn = [this, capture_idx](Canvas &canvas, float opacity) {
      this->draw_shape(canvas, opacity, mesh_states[capture_idx],
                       mesh_states[capture_idx].topology,
                       palettes_history[capture_idx]);
    };

    timeline.add(0, Animation::Sprite(draw_fn, duration, fade_in, ease_mid,
                                      fade_out, ease_mid));

    // Schedule Next overlapping
    int next_delay = duration - fade_out;
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  struct Params {
    float duration = 96.0f;
    float fade = 32.0f;
    int burst_size = 4;
  } params;
};
