/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

#include <array>

// Unit-test accessor reaching the private preset-cycle bookkeeping; the smoke
// harness renders ~120 frames, short of the 288-frame re-spawn period, so those
// paths are driven directly through this seam.
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
      : Effect(W, H, {.strobe = true}), filters(Filter::Screen::AntiAlias<W, H>()),
        mobius_gen(timeline) {}

  /**
   * @brief Builds mesh data, bakes palette LUTs, registers live sliders, and
   *        starts the spawn/spin/orbit animation chain.
   */
  void init() override {
    setup_presets();

    // Bind before any preset bakes through it.
    bloodStreamComposition.bind(&bloodStreamPalette, &bloodStreamFade);

    params = preset_manager.get();
    baked_palettes_[0].bake(persistent_arena, *params.palette);
    baked_palettes_[1].bake(persistent_arena, *params.palette);

    registerParam("Copies", &params.num_copies, 1.0f, 20.0f);
    registerParam("Radius", &params.offset_radius, 0.0f, 1.0f);
    registerParam("Speed", &params.offset_speed, 0.0f, 5.0f);
    registerParam("Warp", &params.warp_scale, 0.0f, 5.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    // Flag the preset-reseeded params so "Pause Animation" lets the user take
    // over.
    for (const char *n : {"Copies", "Radius", "Speed", "Warp", "Alpha"})
      markAnimated(n);

    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    spawn_sprite(0);
    // Wrap the integrated phase to [0,1) so the orbit trig stays in precise range.
    timeline.add(0, Animation::Driver(orbit_phase, &params.offset_speed, 0.01f,
                                      /*wrap=*/true));
  }

  /**
   * @brief Advances the timeline, which drives all spawning and rendering.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    // Mirror live slider edits into the active sprite's snapshot so the incoming
    // shape tracks the sliders while a still-fading outgoing sprite keeps the
    // frozen snapshot it was spawned with.
    param_slots_[active_bake_] = params;
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
  static constexpr int SPRITE_LIFE = 320;  /**< Visible frames per sprite. */
  static constexpr int SPAWN_PERIOD = 288; /**< Frames between spawns. */
  // The two-slot ping-pong is safe only while at most two sprites overlap, i.e.
  // a sprite finishes before the spawn two periods later reuses its slot.
  static_assert(SPRITE_LIFE < 2 * SPAWN_PERIOD,
                "DreamBalls ping-pong needs at most two overlapping sprites");

  BakedPalette baked_palettes_[2];
  int active_bake_ = 0; /**< Index of the slot the next spawn rebakes into. */
  /**
   * @brief Per-sprite render-param snapshots, ping-ponged with baked_palettes_.
   * @details Each sprite renders from the slot it was spawned with so an
   *          outgoing fade keeps its own copy count / radius / warp / alpha
   *          instead of jumping to the incoming preset's values. draw_frame()
   *          mirrors the live slider-bound params into the active slot every
   *          frame, so the incoming sprite stays live-editable while the
   *          outgoing slot is frozen.
   */
  Params param_slots_[2];

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

      // Deep-copy into the persistent arena; the generator's scratch mesh does
      // not outlive this loop.
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

      // Per-vertex tangent frame; seed with X_AXIS near the poles to avoid a
      // degenerate cross product with the near-parallel Y_AXIS.
      data.tangents.bind(persistent_arena, data.mesh_state.vertices.size());
      for (const auto &v : data.mesh_state.vertices) {
        Vector axis = (std::abs(v.y) > 0.99f) ? X_AXIS : Y_AXIS;
        Vector u = cross(v, axis).normalized();
        Vector frame_v = cross(v, u).normalized();
        data.tangents.push_back({u, frame_v});
      }

      // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
      size_t edge_count = data.mesh_state.faces.size() / 2;
      data.edges.bind(persistent_arena, edge_count);
      Plot::Mesh::extract_edges(data.mesh_state, data.edges);
      HS_CHECK(data.edges.size() == edge_count,
               "DreamBalls edge extraction over/under-filled the edge bind");

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

    // Reseed the slider-bound params only when the preset actually changes, so a
    // paused re-spawn of the same preset keeps the user's live edits.
    if (safe_idx != last_preset_idx_) {
      params = entries[safe_idx].params;
      last_preset_idx_ = safe_idx;
    }
    int period = SPAWN_PERIOD;
    // Ping-pong to the inactive slot so the still-fading previous sprite keeps
    // its palette and params (sprites overlap at most pairwise, so the slot two
    // spawns back is already gone). Seed this sprite's param snapshot from the
    // current live params; draw_frame() keeps the active slot tracking sliders.
    active_bake_ ^= 1;
    baked_palettes_[active_bake_].rebake(*params.palette);
    const int bake_slot = active_bake_;
    param_slots_[bake_slot] = params;

    // Bind the warp magnitude to this spawn's scale so dragging "Warp" takes
    // effect this frame. The transformer has a single slot, so a crossfade
    // shares one warp across both sprites; the outgoing warp has relaxed to
    // identity by the period boundary where the next sprite spawns, so the
    // shared magnitude is unnoticeable.
    if (auto *warp = mobius_gen.spawn(0, param_slots_[bake_slot].warp_scale,
                                      period, false))
      warp->bind_scale(param_slots_[bake_slot].warp_scale);

    auto draw_fn = [this, safe_idx, bake_slot](Canvas &canvas, float opacity) {
      const auto &preset = loaded_presets[safe_idx];
      ScratchScope scratch_a_guard(scratch_arena_a);
      MeshState target_mesh;
      MeshOps::transform(preset.mesh_state, target_mesh, scratch_arena_a);

      // This sprite's own param + palette snapshot, so geometry and color both
      // stay continuous across a preset change instead of snapping to the
      // incoming preset.
      this->draw_scene(canvas, param_slots_[bake_slot], opacity,
                       preset.mesh_state, target_mesh, preset.tangents,
                       preset.edges, baked_palettes_[bake_slot]);
    };

    timeline
        .add(0, Animation::Sprite(draw_fn, SPRITE_LIFE, 32, ease_in_out_sin, 32,
                                  ease_in_out_sin))
        .add(period,
             Animation::PeriodicTimer(
                 0,
                 [this, safe_idx](Canvas &) {
                   // Paused: re-spawn the same preset (params hold); otherwise
                   // advance to the next one. Pass the wrapped index so the
                   // re-passed value stays bounded instead of incrementing an int
                   // without limit.
                   this->spawn_sprite(animationsPaused() ? safe_idx
                                                         : safe_idx + 1);
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

    // MeshOps::transform pre-sizes target.vertices to match base; the indexed
    // writes below rely on it.
    HS_CHECK(target.vertices.size() == base.vertices.size(),
             "DreamBalls: displaced-mesh target not pre-sized to base");

    for (size_t i = 0; i < count; ++i) {
      const Vector &v = base.vertices[i];
      const auto &tan = tangents[i];

      // orbit_phase is a fraction of a turn; scale to radians.
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

    // Floor at 1 so the i/num_copies divisor below can't hit zero.
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

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(DreamBalls)
