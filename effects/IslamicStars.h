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
#include <vector>
#include "../solids.h"

template <int W, int H> class IslamicStars : public Effect {

public:
  FLASHMEM IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {
    mesh_states[0] =
        new (persistent_arena.allocate(sizeof(MeshState))) MeshState();
    mesh_states[1] =
        new (persistent_arena.allocate(sizeof(MeshState))) MeshState();

    PersistentTracker::register_mesh(mesh_states[0]);
    PersistentTracker::register_mesh(mesh_states[1]);

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
    timeline.step(canvas);
  }

private:
  Orientation<W> orientation;
  Timeline<W> timeline;
  Pipeline<W, H> filters;
  RippleTransformer<W, 8> ripple_gen;
  float ripple_duration = 80.0f;

  int solid_idx = -1;
  MeshState *mesh_states[2];
  std::vector<int> topologies[2];
  int current_mesh_idx = 0;

  std::array<PaletteVariant, 5> palettes = {
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
                  const std::vector<int> &faceIndices,
                  const std::array<PaletteVariant, 5> &palettes) {
    MemoryCtx ctx;
    ScopedScratch _(ctx);
    MeshState transformed_state;
    OrientTransformer<W> camera(orientation);
    MeshOps::transform(base_state, transformed_state, ctx.get_scratch(),
                       ripple_gen, camera);

    auto fragment_shader = [&](const Vector &p, Fragment &frag) {
      int faceIdx = static_cast<int>(frag.v2);
      int color_idx = 0;
      if (faceIdx >= 0 && faceIdx < (int)faceIndices.size()) {
        color_idx = faceIndices[faceIdx];
      }
      const PaletteVariant &pal = palettes[color_idx];

      float distFromEdge = -frag.v1;
      float size = frag.size;
      float intensity = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
      intensity = hs::clamp(intensity, 0.0f, 1.0f);

      frag.color = get_color(pal, intensity);
      frag.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader);
  }

  void spawn_shape() {
    {
      MemoryCtx math_context;
      solid_idx = (solid_idx + 1) % Solids::Collections::num_islamic_solids;

      {
        ScopedScratch _(math_context);
        SolidGenerator gen(solid_idx + Solids::Collections::num_simple_solids);
        PolyMesh local_mesh =
            gen.generate(math_context.get_scratch(), math_context);

        current_mesh_idx = (current_mesh_idx + 1) % 2;

        math_context.update_persistent(*(mesh_states[current_mesh_idx]),
                                       local_mesh);

        math_context.swap_scratch();
        topologies[current_mesh_idx] = MeshOps::classify_faces_by_topology(
            *(mesh_states[current_mesh_idx]), math_context.get_scratch());
      }
    }

    const auto &entry = Solids::Collections::islamic_solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)mesh_states[current_mesh_idx]->vertices.size(),
            (int)mesh_states[current_mesh_idx]->faces.size());

    // Prepare Palettes
    static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
    std::shuffle(palettes.begin(), palettes.end(), g);

    // Map Topology to palette
    int num_colors = static_cast<int>(palettes.size());
    std::vector<int> &faceIndices = topologies[current_mesh_idx];
    for (size_t i = 0; i < faceIndices.size(); ++i) {
      faceIndices[i] = faceIndices[i] % num_colors;
    }

    // Create Sprite
    int duration = 160;
    int fade_in = 32;
    int fade_out = 32;
    int capture_idx = current_mesh_idx;
    auto pal_copy = palettes;

    auto draw_fn = [this, capture_idx, pal_copy](Canvas &canvas,
                                                 float opacity) {
      this->draw_shape(canvas, opacity, *(mesh_states[capture_idx]),
                       topologies[capture_idx], pal_copy);
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
