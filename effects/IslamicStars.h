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
  enum class SolidType {
    Icosahedron_Hk59_Bitruncate033,
    Octahedron_Hk17_Ambo_Hk72,
    Icosahedron_Kis_Gyro,
    TruncatedIcosidodecahedron_Truncate05_Ambo_Dual,
    Icosidodecahedron_Truncate05_Ambo_Dual,
    SnubDodecahedron_Truncate05_Ambo_Dual,
    Octahedron_Hk34_Ambo_Hk72,
    Rhombicuboctahedron_Hk63_Ambo_Hk63,
    TruncatedIcosahedron_Hk54_Ambo_Hk72,
    Dodecahedron_Hk54_Ambo_Hk72,
    Dodecahedron_Hk72_Ambo_Dual_Hk20,
    TruncatedIcosahedron_Truncate05_Ambo_Dual,
    Last // Sentinel
  };

  IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {
    persist_pixels = false;

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Ripp Amp", &ripple_gen.params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.params.decay, 0.0f, 5.0f);

    // Ripple Duration (Speed is PI / duration)
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 300.0f);

    timeline.add(
        0, Animation::RandomWalk<W>(orientation, UP)); // Slow continuous spin

    // Init Ripple Defaults
    ripple_gen.set_amplitude(0.5f);
    ripple_gen.set_thickness(1.0f);
    ripple_gen.set_decay(0.1f);

    for (int i = 0; i < params.burst_size; i++) {
      timeline.add(i * 16,
                   Animation::PeriodicTimer(
                       96,
                       [this](Canvas &) {
                         ripple_gen.ripple(random_vector(),
                                           static_cast<int>(ripple_duration));
                       },
                       true));
    }
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

  void spawn_shape() {
    solid_idx = (solid_idx + 1) % (int)SolidType::Last;
    PolyMesh mesh = generate_solid((SolidType)solid_idx);
    auto faceIndices = MeshOps::classify_faces_by_topology(mesh);

    // Log Shape Name
    const char *names[] = {"Icosahedron_Hk59_Bitruncate033",
                           "Octahedron_Hk17_Ambo_Hk72",
                           "Icosahedron_Kis_Gyro",
                           "TruncatedIcosidodecahedron_Truncate05_Ambo_Dual",
                           "Icosidodecahedron_Truncate05_Ambo_Dual",
                           "SnubDodecahedron_Truncate05_Ambo_Dual",
                           "Octahedron_Hk34_Ambo_Hk72",
                           "Rhombicuboctahedron_Hk63_Ambo_Hk63",
                           "TruncatedIcosahedron_Hk54_Ambo_Hk72",
                           "Dodecahedron_Hk54_Ambo_Hk72",
                           "Dodecahedron_Hk72_Ambo_Dual_Hk20",
                           "TruncatedIcosahedron_Truncate05_Ambo_Dual"};
    if (solid_idx >= 0 && solid_idx < (int)SolidType::Last) {
      hs::log("Spawning Shape: %s (V=%d, F=%d)", names[solid_idx],
              (int)mesh.vertices.size(), (int)mesh.faces.size());
    }

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

  PolyMesh generate_solid(SolidType type) {
    switch (type) {
    case SolidType::Icosahedron_Hk59_Bitruncate033:
      return Solids::IslamicStarPatterns::icosahedron_hk59_bitruncate033();
    case SolidType::Octahedron_Hk17_Ambo_Hk72:
      return Solids::IslamicStarPatterns::octahedron_hk17_ambo_hk72();
    case SolidType::Icosahedron_Kis_Gyro:
      return Solids::IslamicStarPatterns::icosahedron_kis_gyro();
    case SolidType::TruncatedIcosidodecahedron_Truncate05_Ambo_Dual:
      return Solids::IslamicStarPatterns::
          truncatedIcosidodecahedron_truncate05_ambo_dual();
    case SolidType::Icosidodecahedron_Truncate05_Ambo_Dual:
      return Solids::IslamicStarPatterns::
          icosidodecahedron_truncate05_ambo_dual();
    case SolidType::SnubDodecahedron_Truncate05_Ambo_Dual:
      return Solids::IslamicStarPatterns::
          snubDodecahedron_truncate05_ambo_dual();
    case SolidType::Octahedron_Hk34_Ambo_Hk72:
      return Solids::IslamicStarPatterns::octahedron_hk34_ambo_hk72();
    case SolidType::Rhombicuboctahedron_Hk63_Ambo_Hk63:
      return Solids::IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63();
    case SolidType::TruncatedIcosahedron_Hk54_Ambo_Hk72:
      return Solids::IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72();
    case SolidType::Dodecahedron_Hk54_Ambo_Hk72:
      return Solids::IslamicStarPatterns::dodecahedron_hk54_ambo_hk72();
    case SolidType::Dodecahedron_Hk72_Ambo_Dual_Hk20:
      return Solids::IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20();
    case SolidType::TruncatedIcosahedron_Truncate05_Ambo_Dual:
      return Solids::IslamicStarPatterns::
          truncatedIcosahedron_truncate05_ambo_dual();
    default:
      return Solids::Platonic::dodecahedron();
    }
  }

  struct Params {
    float duration = 96.0f;
    float fade = 32.0f;
    int burst_size = 3;
  } params;
};
