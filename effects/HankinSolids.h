/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

template <int W, int H> class HankinSolids : public Effect {
public:
  // Buffer for morphing operations
  Animation::MorphBuffer morph_buffer;

  // Active compiled states
  CompiledHankin primary_hankin;
  CompiledHankin secondary_hankin;
  MeshState primary_mesh;
  MeshState secondary_mesh;

  HankinSolids() : Effect(W, H), filters() {
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);

    persist_pixels = false;

    // Continuous Random Walk
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS,
                        Animation::RandomWalk<W>::Options::Languid()));

    // Init State
    solid_idx = 0;

    // Initial Geometry Generation
    SolidGenerator gen(solid_idx);
    PolyMesh base;
    gen.generate(base);

    if (enable_hankin) {
      primary_hankin = MeshOps::compile_hankin(base);
      primary_mesh = MeshOps::update_hankin<MeshState>(primary_hankin,
                                                       params.hankin_angle);
    } else {
      primary_mesh = MeshOps::compile(base);
      primary_hankin.dynamicInstructions.clear();
    }

    // Initial Topology & Color
    primary_topology = MeshOps::classify_faces_by_topology(primary_mesh);
    pick_palettes(primary_palettes);

    // Start the Loop
    start_hankin_cycle();
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

  int solid_idx = 0;
  bool enable_dual = false;
  bool enable_hankin = true;

  // Topology & Color State
  std::vector<int> primary_topology;
  std::vector<int> secondary_topology;
  std::vector<const Palette *> primary_palettes;
  std::vector<const Palette *> secondary_palettes;

  // Palette Pool (matching IslamicStars)
  const std::vector<const Palette *> source_palettes_pool = {
      &Palettes::embers, &Palettes::richSunset, &Palettes::brightSunrise,
      &Palettes::bruisedMoss, &Palettes::lavenderLake};

  void pick_palettes(std::vector<const Palette *> &target) {
    target = source_palettes_pool;
    static std::mt19937 g(12345 + (int)timeline.t);
    std::shuffle(target.begin(), target.end(), g);
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;

    // 1. Mutation: Animate Hankin Angle
    // JS: this.timeline.add(0, new Animation.Mutation(this.params,
    // 'hankinAngle', ...)) Here we animate member variable params.hankin_angle
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() { this->start_morph_cycle(); }));

    // 2. Sprite: Draw the mesh
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          if (enable_hankin) {
                            primary_mesh = MeshOps::update_hankin<MeshState>(
                                primary_hankin, params.hankin_angle);
                          }
                          draw_topology_mesh(c, primary_mesh, primary_topology,
                                             primary_palettes, opacity);
                        },
                        DURATION));
  }

  void start_morph_cycle() {
    constexpr int DURATION = 16;

    // Identify Next Solid
    int next_idx = (solid_idx + 1) % Solids::Collections::num_simple_solids;

    const auto &current_entry = Solids::Collections::simple_solids[solid_idx];
    const auto &next_entry = Solids::Collections::simple_solids[next_idx];

    std::cout << "Morphing: " << current_entry.name << " -> " << next_entry.name
              << std::endl;

    // Generate Target Mesh (Secondary)
    SolidGenerator next_gen(next_idx);
    PolyMesh next_base;
    next_gen.generate(next_base);

    if (enable_dual)
      next_base = MeshOps::dual(next_base);

    if (enable_hankin) {
      secondary_hankin = MeshOps::compile_hankin(next_base);
      secondary_mesh = MeshOps::update_hankin<MeshState>(secondary_hankin,
                                                         params.hankin_angle);
    } else {
      secondary_mesh = MeshOps::compile(next_base);
      secondary_hankin.dynamicInstructions.clear();
    }

    // Prepare Secondary State
    secondary_topology = MeshOps::classify_faces_by_topology(secondary_mesh);
    pick_palettes(secondary_palettes);

    // 1. Morph Animation
    timeline.add(0, Animation::MeshMorph(&primary_mesh, &secondary_mesh,
                                         &morph_buffer, primary_mesh,
                                         secondary_mesh, DURATION, false,
                                         ease_in_out_sin)
                        .then([this, next_idx]() {
                          // Commit State
                          this->solid_idx = next_idx;

                          // Move Secondary -> Primary
                          this->primary_topology = this->secondary_topology;
                          this->primary_palettes = this->secondary_palettes;
                          this->primary_hankin = this->secondary_hankin;
                          this->primary_mesh = this->secondary_mesh;

                          // Re-compile Hankin for the NEW solid
                          SolidGenerator new_gen(solid_idx);
                          PolyMesh new_base;
                          new_gen.generate(new_base);

                          if (enable_dual)
                            new_base = MeshOps::dual(new_base);

                          if (enable_hankin) {
                            primary_hankin = MeshOps::compile_hankin(new_base);
                            primary_mesh = MeshOps::update_hankin<MeshState>(
                                primary_hankin, params.hankin_angle);
                          } else {
                            primary_mesh = MeshOps::compile(new_base);
                            primary_hankin.dynamicInstructions.clear();
                          }

                          // Loop
                          this->start_hankin_cycle();
                        }));

    // 2. Sprite: Outgoing
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          draw_topology_mesh(c, primary_mesh, primary_topology,
                                             primary_palettes, opacity);
                        },
                        DURATION, 0, ease_mid, DURATION, ease_mid));

    // 3. Sprite: Incoming
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &c, float opacity) {
                          draw_topology_mesh(c, secondary_mesh,
                                             secondary_topology,
                                             secondary_palettes, opacity);
                        },
                        DURATION, DURATION, ease_mid, 0, ease_mid));
  }

  void draw_topology_mesh(Canvas &canvas, const MeshState &mesh,
                          const std::vector<int> &topology,
                          const std::vector<const Palette *> &palettes,
                          float opacity) {
    if (mesh.vertices.empty())
      return;
    if (opacity < 0.01f)
      return;

    // Create rotated copy
    MeshState rotated_mesh = mesh; // Deep copy of vectors

    Quaternion q = orientation.get();
    for (auto &v : rotated_mesh.vertices) {
      v = rotate(v, q);
    }

    auto shader = [&](const Vector &p, Fragment &f) {
      int faceIdx = (int)std::round(f.v2);
      int topoIdx = 0;
      if (faceIdx >= 0 && faceIdx < (int)topology.size()) {
        topoIdx = topology[faceIdx];
      }

      // Safety: Ensure palettes has content
      if (palettes.empty()) {
        f.color = Color4(1, 0, 1, 1); // Error pink
        return;
      }

      const Palette *pal = palettes[topoIdx % palettes.size()];

      // Edge Distance Intensity (v1 is -dist)
      float distFromEdge = -f.v1;
      float size = f.size;
      float normalizedDist = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
      float t = hs::clamp(normalizedDist * params.intensity, 0.0f, 1.0f);

      f.color = pal->get(t);
      f.color.alpha = opacity;
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader);
  }

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
  } params;
};
