/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
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
    ripple_gen.set_amplitude(0.5f);
    ripple_gen.set_thickness(0.7f);
    ripple_gen.set_decay(0.1f);

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
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  RippleGenerator<W, 32> ripple_gen;
  float ripple_duration = 80.0f;

  MeshState poly_to_state(const PolyMesh &src) const {
    MeshState dst;
    dst.vertices = src.vertices;

    dst.faces.clear();
    dst.face_counts.clear();
    dst.face_counts.reserve(src.faces.size());

    for (const auto &f : src.faces) {
      dst.face_counts.push_back((uint8_t)f.size());
      for (int idx : f)
        dst.faces.push_back(idx);
    }

    build_bvh(dst);
    return dst;
  }

  int solid_idx = -1;

  void ripple(Canvas &canvas) {
    Vector origin = random_vector();
    for (int i = 0; i < params.burst_size; i++) {
      ripple_gen.ripple(origin, static_cast<int>(ripple_duration), i * 16);
    }
  }

  void spawn_shape() {
    solid_idx = (solid_idx + 1) % Solids::Collections::num_islamic_solids;
    const auto &entry = Solids::Collections::islamic_solids[solid_idx];
    PolyMesh mesh = entry.generate();

    auto faceIndices = MeshOps::classify_faces_by_topology(mesh);

    // Log Shape Name
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)mesh.vertices.size(), (int)mesh.faces.size());

    // Flatten for Rendering
    MeshState mesh_state = poly_to_state(mesh);

    // Prepare Palettes
    std::vector<const Palette *> palettes = {
        &Palettes::embers, &Palettes::richSunset, &Palettes::brightSunrise,
        &Palettes::bruisedMoss, &Palettes::lavenderLake};
    static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
    std::shuffle(palettes.begin(), palettes.end(), g);

    // Create Sprite
    int duration = 96;
    int fade_in = 32;
    int fade_out = 32;
    auto draw_fn = [this, mesh_state, faceIndices, palettes](Canvas &canvas,
                                                             float opacity) {
      std::vector<Vector> transformed_vertices =
          mesh_state.vertices; // Copy only vertices
      Quaternion q = orientation.get();
      for (auto &v : transformed_vertices) {
        v = ripple_gen.apply(v);
        v = rotate(v, q);
      }

      // 2. Fragment Shader
      auto fragment_shader = [&](const Vector &p, Fragment &frag) {
        int faceIdx = (int)std::round(frag.v2);
        int topoIdx = 0;
        if (faceIdx >= 0 && faceIdx < (int)faceIndices.size()) {
          topoIdx = faceIndices[faceIdx];
        }
        const Palette *pal = palettes[topoIdx % palettes.size()];

        float distFromEdge = -frag.v1;
        float size = frag.size;
        float intensity = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
        intensity = std::clamp(intensity, 0.0f, 1.0f);

        frag.color = pal->get(intensity);
        frag.color.alpha = opacity;
      };

      // 3. Draw using Scan::Mesh::drawRef
      struct MeshRef {
        const std::vector<Vector> &vertices;
        const std::vector<int> &faces;
        const std::vector<uint8_t> &face_counts;
        size_t num_faces;
        size_t num_vertices;
      } mesh_ref = {transformed_vertices, mesh_state.faces,
                    mesh_state.face_counts, mesh_state.face_counts.size(),
                    transformed_vertices.size()};

      Scan::Mesh::draw<W, H>(filters, canvas, mesh_ref, fragment_shader);
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
