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
#include <string>
#include <vector>

template <int W, int H> class IslamicStars : public Effect {

public:
  IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {
    persist_pixels = false;

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Ripp Amp", &ripple_gen.params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.params.decay, 0.0f, 5.0f);

    // Ripple Duration (Speed is PI / duration)
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 300.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP));

    // Init Ripple Defaults
    ripple_gen.params.amplitude = 0.4f;
    ripple_gen.params.thickness = 0.7f;
    ripple_gen.params.decay = 0.1f;

    // Ripple now and schedule more
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

  void ripple(Canvas &canvas) {
    Vector origin = random_vector();
    for (int i = 0; i < params.burst_size; i++) {
      // Speed = PI / duration.
      // Generic spawn: delay, center, speed, duration
      ripple_gen.spawn(i * 16, origin, PI_F / ripple_duration,
                       static_cast<int>(ripple_duration));
    }
  }

  void spawn_shape() {
    solid_idx = (solid_idx + 1) % Solids::Collections::num_islamic_solids;

    // 1. GENERATE
    SolidGenerator gen(solid_idx + Solids::Collections::num_simple_solids);
    PolyMesh local_mesh;
    gen.generate(local_mesh);

    // Log Shape Name
    const auto &entry = Solids::Collections::islamic_solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)local_mesh.vertices.size(), (int)local_mesh.faces.size());

    // 2. COMPILE
    MeshState base_state = MeshOps::compile(local_mesh);

    auto faceIndices = MeshOps::classify_faces_by_topology(base_state);

    // Prepare Palettes
    std::vector<const Palette *> palettes = {
        &Palettes::embers, &Palettes::richSunset, &Palettes::brightSunrise,
        &Palettes::bruisedMoss, &Palettes::lavenderLake};
    static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
    std::shuffle(palettes.begin(), palettes.end(), g);

    // Map the topology to the effect's specific visual needs (the palette)
    int num_colors = static_cast<int>(palettes.size());
    for (size_t i = 0; i < faceIndices.size(); ++i) {
      faceIndices[i] = faceIndices[i] % num_colors;
    }

    // Create Sprite
    int duration = 160;
    int fade_in = 32;
    int fade_out = 32;
    auto draw_fn = [this, base_state, faceIndices, palettes](Canvas &canvas,
                                                             float opacity) {
      // 3. TRANSFORM
      MeshState transformed_state = base_state;
      OrientTransformer<W> camera(orientation);
      MeshOps::transform(base_state, transformed_state, ripple_gen, camera);

      // 4. SHADER
      auto fragment_shader = [&](const Vector &p, Fragment &frag) {
        int faceIdx = static_cast<int>(frag.v2);
        int color_idx = 0;
        if (faceIdx >= 0 && faceIdx < (int)faceIndices.size()) {
          color_idx = faceIndices[faceIdx];
        }
        const Palette *pal = palettes[color_idx];

        float distFromEdge = -frag.v1;
        float size = frag.size;
        float intensity = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
        intensity = hs::clamp(intensity, 0.0f, 1.0f);

        frag.color = pal->get(intensity);
        frag.color.alpha = opacity;
      };

      // 5. RASTERIZE
      Scan::Mesh::draw<W, H>(filters, canvas, transformed_state,
                             fragment_shader);
    };

    timeline.add(0, Animation::Sprite(draw_fn, duration, fade_in, ease_mid,
                                      fade_out, ease_mid));

    // Schedule Next
    // Overlap = fade. Next delay = duration - overlap.
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
