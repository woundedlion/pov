/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

#include <array>

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// reaches the private preset-cycle bookkeeping (spawn_sprite, active_bake_,
// last_preset_idx_). The smoke/determinism harness only renders ~120 frames,
// short of the 288-frame re-spawn period, so the preset advance / bake
// ping-pong / reseed-guard paths are driven directly through this seam.
namespace hs_test {
namespace effects_tests {
struct DreamBallsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Orbiting copies of a polyhedral wireframe, cycling through four
 *        solid presets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each vertex is displaced in its own tangent plane and
 *          Mobius-warped before being re-projected onto the sphere.
 */
template <int W, int H> class DreamBalls : public Effect {
public:
  /**
   * @brief Live, slider-bound render parameters; also the per-preset value set.
   */
  struct Params {
    const char *solid_name;
    float num_copies;
    float offset_radius;
    float offset_speed;
    float warp_scale;
    const Palette *palette;
    float alpha;
  };

  /**
   * @brief Constructs the effect with the anti-alias screen filter.
   * @details The Mobius generator hangs off the timeline so its warps animate
   *          in step.
   */
  FLASHMEM DreamBalls()
      : Effect(W, H), filters(Filter::Screen::AntiAlias<W, H>()),
        mobius_gen(timeline) {}

  /**
   * @brief Builds mesh data, bakes palette LUTs, registers live sliders, and
   *        starts the spawn/spin/orbit animation chain.
   */
  void init() override {
    setup_presets();

    // Wire the blood-stream composition before any preset bakes through it.
    bloodStreamComposition.bind(&bloodStreamPalette, &bloodStreamFade);

    params = preset_manager.get();
    // Bake both LUT slots up front so each owns its persistent-arena allocation;
    // spawn_sprite only ever rebakes (no further allocation) thereafter.
    baked_palettes_[0].bake(persistent_arena, *params.palette);
    baked_palettes_[1].bake(persistent_arena, *params.palette);

    registerParam("Copies", &params.num_copies, 1.0f, 20.0f);
    registerParam("Radius", &params.offset_radius, 0.0f, 1.0f);
    registerParam("Speed", &params.offset_speed, 0.0f, 5.0f);
    registerParam("Warp", &params.warp_scale, 0.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    // Each preset cycle reseeds these from the preset; flag them so the standard
    // "Pause Animation" toggle freezes the cycle and lets the user take over.
    for (const char *n : {"Copies", "Radius", "Speed", "Warp", "Alpha"})
      markAnimated(n);

    // Start Sequence
    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    spawn_sprite(0);
    // Accumulate orbit phase by integrating the live Speed slider each frame
    // (wrapped to [0,1)) so a Speed change alters the rate going forward only
    // rather than retroactively rescaling all elapsed phase. The wrap keeps fast
    // trig in precise range over multi-day runs.
    timeline.add(0, Animation::Driver(orbit_phase, &params.offset_speed, 0.01f,
                                      /*wrap=*/true));
  }

  /**
   * @brief Reports whether the engine should clear the background each frame.
   * @return Always false; sprites fade over the prior frame instead.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances the timeline, which drives all spawning and rendering.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /** @brief Live, slider-bound render parameters for the active preset. */
  Params params;

private:
  friend struct ::hs_test::effects_tests::DreamBallsWhiteBox;

  /** Orbit phase in turns, wrapped to [0,1) by the live-speed Driver below. */
  float orbit_phase = 0.0f;
  int last_preset_idx_ = -1; /**< Last preset whose values were copied into params. */

  /** Per-vertex phase increment (radians) for the orbit stagger: vertex i
       leads the next by this much so the surface ripples instead of pulsing in
       unison. */
  static constexpr float VERTEX_PHASE_STAGGER = 0.1f;

  /**
   * @brief Orthonormal tangent-plane basis (u, v) at a vertex.
   * @details The orbit displacement is spanned by this frame before
   *          re-projecting to the sphere.
   */
  struct Tangent {
    Vector u; /**< First tangent basis vector. */
    Vector v; /**< Second tangent basis vector, orthogonal to u. */
  };

  /**
   * @brief Precomputed, static per-preset geometry baked into the persistent
   *        arena.
   */
  struct PresetData {
    MeshState mesh_state;          /**< Baked vertices and faces. */
    ArenaVector<Tangent> tangents; /**< Per-vertex tangent frames. */
    ArenaVector<Plot::Mesh::Edge> edges; /**< Unique edge list (topology is static). */
  };

  std::array<PresetData, 4> loaded_presets;

  FastNoiseLite noise;
  Timeline timeline;

  Orientation<> global_orientation;

  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  MobiusWarpTransformer<1> mobius_gen;
  /**
   * @brief Two baked LUTs, ping-ponged per spawn.
   * @details On a preset change at most two sprites overlap (the outgoing
   *          fade-out and the incoming fade-in); each sprite captures its own
   *          slot so the outgoing one keeps the palette it was spawned with
   *          through its fade instead of hard-cutting to the freshly rebaked
   *          colors. The non-color live params stay shared on purpose
   *          (Copies/Radius/Speed/Alpha slider liveness).
   */
  BakedPalette baked_palettes_[2];
  int active_bake_ = 0; /**< Index of the slot the next spawn rebakes into. */

  ProceduralPalette bloodStreamPalette = Palettes::bloodStream;
  AlphaFalloffShade bloodStreamFade{[](float t) { return 1.0f - t; }};
  /**
   * @brief Composition of the bloodStream palette under the alpha falloff
   *        shade.
   * @details Wrap=false samples the source at the raw coordinate (the falloff
   *          curve owns the [0,1] domain, so no wrap). The facade lets the
   *          composition sit in the polymorphic preset table beside plain
   *          palettes; it is only ever baked.
   */
  StaticPalette<ProceduralPalette, Coords<>, Colors<AlphaFalloffShade>,
                /*Wrap=*/false>
      bloodStreamComposition;
  PaletteFacade<decltype(bloodStreamComposition)> bloodStreamFalloff{
      &bloodStreamComposition};

  Presets<Params, 4> preset_manager{std::array<PresetEntry<Params>, 4>{{
      {{"rhombicuboctahedron", 18.0f, 0.3f, 0.4f, 0.3f,
        &bloodStreamFalloff, 0.7f}},
      {{"rhombicosidodecahedron", 6.0f, 0.05f, 1.0f, 1.8f,
        &bloodStreamFalloff, 0.7f}},
      {{"truncatedCuboctahedron", 6.0f, 0.16f, 1.0f, 2.0f,
        &Palettes::richSunset, 0.3f}},
      {{"icosidodecahedron", 10.0f, 0.16f, 1.0f, 0.5f,
        &Palettes::lavenderLake, 0.3f}}}}};

  /**
   * @brief Generates each preset's solid and bakes its geometry into the
   *        persistent arena once at init.
   * @details Bakes vertices, faces, tangent frames, and the unique edge list
   *          for every entry in the preset table.
   */
  FLASHMEM void setup_presets() {
    const auto &entries = preset_manager.get_entries();

    int preset_idx = 0;
    for (const auto &entry : entries) {
      const auto &p = entry.params;
      auto &data = loaded_presets[preset_idx];

      PolyMesh m = generate(persistent_arena, Solids::get_by_name,
                            std::string_view(p.solid_name));

      // Deep-copy vertices into the persistent arena (the generator's scratch
      // mesh does not outlive this loop).
      data.mesh_state.vertices.bind(persistent_arena, m.vertices.size());
      for (const auto &v : m.vertices) {
        data.mesh_state.vertices.push_back(v);
      }

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

      // Build a tangent frame per vertex; pick X_AXIS as the seed near the poles
      // to avoid a degenerate cross product with the near-parallel Y_AXIS.
      data.tangents.bind(persistent_arena, data.mesh_state.vertices.size());
      for (const auto &v : data.mesh_state.vertices) {
        Vector axis = (std::abs(v.y) > 0.99f) ? X_AXIS : Y_AXIS;
        Vector u = cross(v, axis).normalized();
        Vector frame_v = cross(v, u).normalized();
        data.tangents.push_back({u, frame_v});
      }

      // Precompute unique edge list (topology is static). For a closed
      // 2-manifold every edge is shared by exactly two faces, so the flattened
      // corner count (faces.size() = Σ face degrees) is exactly 2·E — size to E.
      // These presets are all closed Solids:: polyhedra; a non-manifold/open
      // mesh would under-size and trap loudly in push_back (and extract_edges
      // already HS_CHECKs vertex indices), matching the engine's fail-fast
      // contract rather than carrying a standing 2× persistent over-allocation.
      size_t edge_count = data.mesh_state.faces.size() / 2;
      data.edges.bind(persistent_arena, edge_count);
      Plot::Mesh::extract_edges(data.mesh_state, data.edges);

      preset_idx++;
    }
  }

  /**
   * @brief Spawns one fading sprite for the given preset and schedules the next
   *        spawn one period later.
   * @param idx Preset index; taken modulo the preset count.
   * @details Rebakes the inactive palette slot, arms a Mobius warp, and reseeds
   *          the live params only when the preset actually changes.
   */
  void spawn_sprite(int idx) {
    auto entries = preset_manager.get_entries();
    int safe_idx = idx % entries.size();

    // Reseed the slider-bound params from the preset only when the preset
    // actually changes. While paused the chain re-spawns the SAME preset (see
    // the scheduler below), so the user's live edits persist instead of being
    // clobbered each cycle.
    if (safe_idx != last_preset_idx_) {
      params = entries[safe_idx].params;
      last_preset_idx_ = safe_idx;
    }
    int period = 288;
    // Bind the warp magnitude to the live "Warp" slider so dragging it takes
    // effect this frame instead of only at the next spawn (~5 s later). The
    // returned pointer is used transiently — never retained — so it is safe
    // against timeline compaction.
    if (auto *warp = mobius_gen.spawn(0, this->params.warp_scale, period, false))
      warp->bind_scale(this->params.warp_scale);
    // Rebake into the inactive slot so the previous sprite (still fading out)
    // keeps reading the palette it spawned with. The slot two spawns back is
    // long gone (sprites overlap at most pairwise), so this never disturbs a
    // live sprite.
    active_bake_ ^= 1;
    baked_palettes_[active_bake_].rebake(*params.palette);
    const int bake_slot = active_bake_;

    auto draw_fn = [this, safe_idx, bake_slot](Canvas &canvas, float opacity) {
      const auto &preset = loaded_presets[safe_idx];
      ScratchScope _(scratch_arena_a);
      MeshState target_mesh;
      MeshOps::transform(preset.mesh_state, target_mesh, scratch_arena_a);

      // Render from the live, slider-bound params (not a fresh preset copy) so
      // the registered Copies/Radius/Speed/Alpha sliders actually affect the
      // output, but from this sprite's own baked LUT so its color stays
      // continuous across a preset change.
      this->draw_scene(canvas, this->params, opacity, preset.mesh_state,
                       target_mesh, preset.tangents, preset.edges,
                       baked_palettes_[bake_slot]);
    };

    timeline
        .add(0, Animation::Sprite(draw_fn, 320, 32, ease_in_out_sin, 32,
                                  ease_in_out_sin))
        .add(period,
             Animation::PeriodicTimer(
                 0,
                 [this, idx](Canvas &) {
                   // Paused: re-spawn the same preset (params hold); otherwise
                   // advance to the next one.
                   this->spawn_sprite(animationsPaused() ? idx : idx + 1);
                 },
                 false));
  }

  /**
   * @brief Orbits each vertex in a small circle within its own tangent plane,
   *        then re-projects onto the unit sphere.
   * @param base Source mesh whose vertices are orbited.
   * @param target Destination mesh receiving the displaced vertices.
   * @param tangents Per-vertex (u, v) tangent frames spanning the orbit plane.
   * @param p Render params; p.offset_radius is the orbit radius.
   * @param angle_offset Per-copy phase offset in radians.
   * @details The per-vertex phase (i * VERTEX_PHASE_STAGGER) staggers the
   *          orbits so the surface ripples.
   */
  void update_displaced_mesh(const MeshState &base, MeshState &target,
                             const ArenaVector<Tangent> &tangents,
                             const Params &p, float angle_offset) {
    size_t count = base.vertices.size();
    float r = p.offset_radius;

    for (size_t i = 0; i < count; ++i) {
      const Vector &v = base.vertices[i];
      const auto &tan = tangents[i];

      // orbit_phase already integrates Speed (see the Driver in init()); it is a
      // fraction of a turn, so scale to radians here. Per-vertex stagger and the
      // copy's angle_offset spread the orbits in phase.
      float phase = i * VERTEX_PHASE_STAGGER;
      float angle = orbit_phase * 2 * PI_F + phase + angle_offset;

      float cosA = fast_cosf(angle);
      float sinA = fast_sinf(angle);

      Vector disp = v + (tan.u * cosA + tan.v * sinA) * r;
      target.vertices[i] = disp.normalized();
    }
  }

  /**
   * @brief Draws p.num_copies orbiting wireframe shells of the preset's solid.
   * @param canvas Render target.
   * @param p Live render params (copy count, radius, alpha, etc.).
   * @param opacity Sprite fade factor in [0,1]; multiplies p.alpha.
   * @param base Source mesh supplying the undisplaced vertices.
   * @param target Scratch mesh reused for each copy's displaced vertices.
   * @param tangents Per-vertex tangent frames for the displacement.
   * @param edges Unique edge list defining the wireframe topology.
   * @param baked Baked palette LUT supplying edge colors.
   * @details Each copy displaces vertices in their tangent frames (staggered in
   *          phase by an even 2*PI/num_copies offset), Mobius-warps and orients
   *          them, then plots the edges shaded from the baked palette at
   *          p.alpha * opacity.
   */
  void draw_scene(Canvas &canvas, const Params &p, float opacity,
                  const MeshState &base, MeshState &target,
                  const ArenaVector<Tangent> &tangents,
                  const ArenaVector<Plot::Mesh::Edge> &edges,
                  const BakedPalette &baked) {

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      Color4 c = baked.get(f.v0);
      c.alpha *= p.alpha * opacity;
      f.color = c;
    };

    // num_copies is a float param (slider/preset), but it is a count: hoist it
    // to an int bound so the loop compares int-to-int instead of promoting the
    // counter to float each iteration. Identical for the integer values it ever
    // holds (presets and the slider's whole-number stops). Clamp to >= 1 so the
    // "at least one shell" contract and the i/num_copies divisor hold even if the
    // registered min (currently 1.0) is ever lowered into [0, 1) — without it a
    // sub-1 value would truncate the bound to 0 and silently draw nothing.
    const int num_copies_raw = static_cast<int>(p.num_copies);
    const int num_copies = num_copies_raw < 1 ? 1 : num_copies_raw;
    for (int i = 0; i < num_copies; ++i) {
      float offset = (static_cast<float>(i) / num_copies) * 2 * PI_F;
      update_displaced_mesh(base, target, tangents, p, offset);

      for (size_t vi = 0; vi < target.vertices.size(); ++vi) {
        target.vertices[vi] = mobius_gen.transform(target.vertices[vi]);
        target.vertices[vi] = global_orientation.orient(target.vertices[vi]);
      }

      Plot::Mesh::draw<W, H>(filters, canvas, target, edges,
                             fragment_shader);
    }
  }

  /**
   * @brief Kicks off a full-turn rotation of the global orientation about a
   *        fresh random axis.
   * @details Scheduled periodically to keep the whole cluster slowly tumbling.
   */
  void spin_slices() {
    Vector axis = random_vector();
    timeline.add(0, Animation::Rotation<W>(global_orientation, axis, 2 * PI_F,
                                           80, ease_in_out_sin, false));
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(DreamBalls)
