/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../transformers.h"
#include "../generators.h"

#include "../effects_engine.h"
#include <vector>
#include <map>
#include <map>

template <int W, int H> class DreamBalls : public Effect {
public:
  struct Params {
    std::string solid_name;
    float num_copies;
    float offset_radius;
    float offset_speed;
    float warp_scale;
    const Palette *palette;
    float alpha;
    bool enable_slice = false;
  };

  DreamBalls()
      : Effect(W, H),
        filters(Filter::World::OrientSlice<W>(orientations, Y_AXIS),
                Filter::Screen::AntiAlias<W, H>()),

        slice_filter(filters), // Filters inherits Head (FilterOrientSlice)
        mobius_gen(timeline) {
    persist_pixels = false;

    // Initialize Orientations
    orientations.resize(2); // 2 slices

    // Initialize Presets
    setup_presets();

    // Initialize default params
    // if (!presets.empty()) params = presets[0];
    // replacing with preset_manager
    params = preset_manager.get();

    registerParam("Copies", &params.num_copies, 1.0f, 20.0f);
    registerParam("Radius", &params.offset_radius, 0.0f, 1.0f);
    registerParam("Speed", &params.offset_speed, 0.0f, 5.0f);
    registerParam("Warp", &params.warp_scale, 0.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Start Sequence
    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &c) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS,
                        Animation::RandomWalk<W>::Options::Languid()));

    spawn_sprite(0); // Start first preset
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);

    slice_filter.get().enabled = false;

    t += 0.01f;
    timeline.step(canvas);
  }

  // Exposed for consistency
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
    std::vector<Tangent> tangents;
  };

  std::vector<PresetData> loaded_presets;

  Timeline<W> timeline;

  std::vector<Orientation<W>> orientations;
  Orientation<W> global_orientation;

  Pipeline<W, H, Filter::World::OrientSlice<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
  std::reference_wrapper<Filter::World::OrientSlice<W>> slice_filter;
  std::vector<Params> presets; // DELETE
  MobiusWarpTransformer<W, 64> mobius_gen;

  AlphaFalloffPalette bloodStreamFalloff{[](float t) { return 1.0f - t; },
                                         Palettes::bloodStream};

  Presets<Params> preset_manager = {{"rhombicuboctahedron",
                                     {"rhombicuboctahedron", 18.0f, 0.3f, 0.4f,
                                      0.3f, &bloodStreamFalloff, 0.7f}},
                                    {"rhombicosidodecahedron",
                                     {"rhombicosidodecahedron", 6.0f, 0.05f,
                                      1.0f, 1.8f, &bloodStreamFalloff, 0.7f}},
                                    {"truncatedCuboctahedron",
                                     {"truncatedCuboctahedron", 6.0f, 0.16f,
                                      1.0f, 2.0f, &Palettes::richSunset, 0.3f}},
                                    {"icosidodecahedron",
                                     {"icosidodecahedron", 10.0f, 0.16f, 1.0f,
                                      0.5f, &Palettes::lavenderLake, 0.3f}}};

  void setup_presets() {
    // Pre-load Geometry
    const auto &entries = preset_manager.get_entries();
    loaded_presets.reserve(entries.size());

    for (const auto &entry : entries) {
      const auto &p = entry.params;
      loaded_presets.emplace_back();
      auto &data = loaded_presets.back();

      SolidNameGenerator gen(p.solid_name);
      PolyMesh m;
      gen.generate(m);

      // Store Verts
      data.mesh_state.vertices = m.vertices;

      // Store Faces
      data.mesh_state.faces.clear();
      data.mesh_state.face_counts.clear();
      for (const auto &f : m.faces) {
        data.mesh_state.face_counts.push_back((uint8_t)f.size());
        for (int idx : f)
          data.mesh_state.faces.push_back(idx);
      }

      // Compute Tangents
      for (const auto &v : data.mesh_state.vertices) {
        Vector axis = (std::abs(v.j) > 0.99f) ? X_AXIS : Y_AXIS;
        Vector u = cross(v, axis).normalize();
        Vector frame_v = cross(v, u).normalize();
        data.tangents.push_back({u, frame_v});
      }
    }
  }

  void spawn_sprite(int idx) {
    auto &entries = preset_manager.get_entries();
    int safe_idx = idx % entries.size();

    // safe_idx corresponds to loaded_presets index too since we iterated
    // entries
    params = entries[safe_idx].params;
    Params instance_params = entries[safe_idx].params;
    const PresetData *preset_ptr = &loaded_presets[safe_idx];
    int period = 288;
    mobius_gen.spawn(0, this->params.warp_scale, period, true);

    auto draw_fn = [this, preset_ptr, instance_params](Canvas &canvas,
                                                       float opacity) {
      MeshState target_mesh = preset_ptr->mesh_state;
      this->draw_scene(canvas, instance_params, opacity, preset_ptr->mesh_state,
                       target_mesh, preset_ptr->tangents);
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
                             const std::vector<Tangent> &tangents,
                             const Params &p, float angle_offset) {
    size_t count = base.vertices.size();
    float r = p.offset_radius;
    float speed = p.offset_speed;

    for (size_t i = 0; i < count; ++i) {
      const Vector &v = base.vertices[i];
      const auto &tan = tangents[i];

      float phase = i * 0.1f;
      float angle = t * speed * 2 * PI_F + phase + angle_offset;

      float cosA = cosf(angle);
      float sinA = sinf(angle);

      Vector disp = v + (tan.u * cosA + tan.v * sinA) * r;
      target.vertices[i] = disp.normalize();
    }
  }

  void draw_scene(Canvas &canvas, const Params &p, float opacity,
                  const MeshState &base, MeshState &target,
                  const std::vector<Tangent> &tangents) {

    auto vertex_shader = [&](Fragment &f) {
      f.pos = mobius_gen.transform(f.pos);
      f.pos = global_orientation.orient(f.pos);
    };

    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      Color4 c = p.palette->get(f.v0);
      c.alpha *= p.alpha * opacity;
      f.color = c;
    };

    for (int i = 0; i < p.num_copies; ++i) {
      float offset = (static_cast<float>(i) / p.num_copies) * 2 * PI_F;
      update_displaced_mesh(base, target, tangents, p, offset);
      Plot::Mesh::draw<W, H>(filters, canvas, target, fragment_shader,
                             vertex_shader);
    }
  }

  void spin_slices() {
    Vector axis = random_vector();
    slice_filter.get().axis = axis;

    for (size_t i = 0; i < orientations.size(); ++i) {
      float direction = (i % 2 == 0) ? 1.0f : -1.0f;
      timeline.add(0, Animation::Rotation<W>(orientations[i], axis,
                                             direction * 2 * PI_F, 80,
                                             ease_in_out_sin, false));
    }
  }
};
