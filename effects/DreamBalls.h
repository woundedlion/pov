/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

#include <array>

template <int W, int H> class DreamBalls : public Effect {
public:
  struct Params {
    const char *solid_name;
    float num_copies;
    float offset_radius;
    float offset_speed;
    float warp_scale;
    const Palette *palette;
    float alpha;
    bool enable_slice = false;
  };

  FLASHMEM DreamBalls()
      : Effect(W, H), filters(Filter::Screen::AntiAlias<W, H>()),
        mobius_gen(timeline) {}

  void init() override {
    setup_presets();

    params = preset_manager.get();
    baked_palette.bake(persistent_arena, *params.palette);

    registerParam("Copies", &params.num_copies, 1.0f, 20.0f);
    registerParam("Radius", &params.offset_radius, 0.0f, 1.0f);
    registerParam("Speed", &params.offset_speed, 0.0f, 5.0f);
    registerParam("Warp", &params.warp_scale, 0.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Start Sequence
    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &c) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    spawn_sprite(0);
    timeline.add(0, Animation::Driver(t, 0.01f, false));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  Params params;

private:
  float t = 0;

  struct Tangent {
    Vector u;
    Vector v;
  };

  struct PresetData {
    // Base Mesh Data
    MeshState mesh_state;
    ArenaVector<Tangent> tangents;
    ArenaVector<Plot::Mesh::Edge> edges;
  };

  std::array<PresetData, 4> loaded_presets;

  FastNoiseLite noise;
  Timeline<W> timeline;

  Orientation<W> global_orientation;

  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  MobiusWarpTransformer<W, 1> mobius_gen;
  BakedPalette baked_palette;

  ProceduralPalette bloodStreamPalette = Palettes::bloodStream;
  AlphaFalloffPalette bloodStreamFalloff{[](float t) { return 1.0f - t; },
                                         &bloodStreamPalette};

  Presets<Params, 4> preset_manager = {
      .entries = {{{{"rhombicuboctahedron", 18.0f, 0.3f, 0.4f, 0.3f,
                     &bloodStreamFalloff, 0.7f, false}},
                   {{"rhombicosidodecahedron", 6.0f, 0.05f, 1.0f, 1.8f,
                     &bloodStreamFalloff, 0.7f, false}},
                   {{"truncatedCuboctahedron", 6.0f, 0.16f, 1.0f, 2.0f,
                     &Palettes::richSunset, 0.3f, false}},
                   {{"icosidodecahedron", 10.0f, 0.16f, 1.0f, 0.5f,
                     &Palettes::lavenderLake, 0.3f, false}}}},
      .current_idx = 0};

  FLASHMEM void setup_presets() {
    const auto &entries = preset_manager.get_entries();

    int preset_idx = 0;
    for (const auto &entry : entries) {
      const auto &p = entry.params;
      auto &data = loaded_presets[preset_idx];

      PolyMesh m = generate(persistent_arena, Solids::get_by_name,
                            std::string_view(p.solid_name));

      // Store Verts (Deep Copy)
      data.mesh_state.vertices.bind(persistent_arena, m.vertices.size());
      for (const auto &v : m.vertices) {
        data.mesh_state.vertices.push_back(v);
      }

      // Store Faces
      data.mesh_state.faces.bind(persistent_arena, m.faces.size());
      data.mesh_state.face_counts.bind(persistent_arena, m.face_counts.size());

      int flat_idx = 0;
      for (size_t i = 0; i < m.face_counts.size(); ++i) {
        int count = m.face_counts[i];
        data.mesh_state.face_counts.push_back((uint8_t)count);
        for (int c = 0; c < count; ++c) {
          data.mesh_state.faces.push_back(m.faces[flat_idx++]);
        }
      }

      // Compute Tangents
      data.tangents.bind(persistent_arena, data.mesh_state.vertices.size());
      for (const auto &v : data.mesh_state.vertices) {
        Vector axis = (std::abs(v.y) > 0.99f) ? X_AXIS : Y_AXIS;
        Vector u = cross(v, axis).normalized();
        Vector frame_v = cross(v, u).normalized();
        data.tangents.push_back({u, frame_v});
      }

      // Precompute unique edge list (topology is static)
      size_t max_edges = data.mesh_state.faces.size(); // upper bound
      data.edges.bind(persistent_arena, max_edges);
      Plot::Mesh::extract_edges(data.mesh_state, data.edges);

      preset_idx++;
    }
  }

  void spawn_sprite(int idx) {
    auto entries = preset_manager.get_entries();
    int safe_idx = idx % entries.size();

    params = entries[safe_idx].params;
    int period = 288;
    mobius_gen.spawn(0, this->params.warp_scale, period, false);
    baked_palette.rebake(*params.palette);

    auto draw_fn = [this, safe_idx](Canvas &canvas, float opacity) {
      const auto &preset = loaded_presets[safe_idx];
      Params ip = preset_manager.get_entries()[safe_idx].params;
      ScratchScope _(scratch_arena_a);
      MeshState target_mesh;
      MeshOps::transform(preset.mesh_state, target_mesh, scratch_arena_a);

      this->draw_scene(canvas, ip, opacity, preset.mesh_state, target_mesh,
                       preset.tangents, preset.edges);
    };

    timeline
        .add(0, Animation::Sprite(draw_fn, 320, 32, ease_in_out_sin, 32,
                                  ease_in_out_sin))
        .add(period,
             Animation::PeriodicTimer(
                 0, [this, idx](Canvas &c) { this->spawn_sprite(idx + 1); },
                 false));
  }

  void update_displaced_mesh(const MeshState &base, MeshState &target,
                             const ArenaVector<Tangent> &tangents,
                             const Params &p, float angle_offset) {
    size_t count = base.vertices.size();
    float r = p.offset_radius;
    float speed = p.offset_speed;

    for (size_t i = 0; i < count; ++i) {
      const Vector &v = base.vertices[i];
      const auto &tan = tangents[i];

      float phase = i * 0.1f;
      float angle = t * speed * 2 * PI_F + phase + angle_offset;

      float cosA = fast_cosf(angle);
      float sinA = fast_sinf(angle);

      Vector disp = v + (tan.u * cosA + tan.v * sinA) * r;
      target.vertices[i] = disp.normalized();
    }
  }

  void draw_scene(Canvas &canvas, const Params &p, float opacity,
                  const MeshState &base, MeshState &target,
                  const ArenaVector<Tangent> &tangents,
                  const ArenaVector<Plot::Mesh::Edge> &edges) {

    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      Color4 c = baked_palette.get(f.v0);
      c.alpha *= p.alpha * opacity;
      f.color = c;
    };

    for (int i = 0; i < p.num_copies; ++i) {
      float offset = (static_cast<float>(i) / p.num_copies) * 2 * PI_F;
      update_displaced_mesh(base, target, tangents, p, offset);

      for (size_t vi = 0; vi < target.vertices.size(); ++vi) {
        target.vertices[vi] = mobius_gen.transform(target.vertices[vi]);
        target.vertices[vi] = global_orientation.orient(target.vertices[vi]);
      }

      Plot::Mesh::draw<W, H>(filters, canvas, target, edges,
                             fragment_shader);
    }
  }

  void spin_slices() {
    Vector axis = random_vector();
    timeline.add(0, Animation::Rotation<W>(global_orientation, axis, 2 * PI_F,
                                           80, ease_in_out_sin, false));
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(DreamBalls)
