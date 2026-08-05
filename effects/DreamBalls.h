/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file DreamBalls.h
 * @brief Orbiting copies of a polyhedral wireframe, cycling through solid
 *        presets.
 */

#include "core/engine/engine.h"

#include <array>

// Unit-test accessor reaching the private preset-cycle bookkeeping; the smoke
// harness renders ~120 frames, short of the 320-frame re-spawn period, so those
// paths are driven directly through this seam.
namespace hs_test {
namespace effects_tests {
struct DreamBallsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Orbiting copies of a polyhedral wireframe, cycling through five
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
    const char *solid_name = nullptr;
    float num_copies = COPIES_MIN;
    float offset_radius = RADIUS_MIN;
    float offset_speed = SPEED_MIN;
    float warp_scale = WARP_MIN;
    const Palette *palette = nullptr;
    float alpha = ALPHA_MIN;
  };

  /**
   * @brief Constructs the effect with the anti-alias screen filter.
   * @details The Mobius generator hangs off the timeline so its warps animate
   *          in step.
   */
  HS_COLD_MEMBER DreamBalls()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        filters(Filter::Screen::AntiAlias<W, H>()), mobius_gen(timeline) {}

  /**
   * @brief Builds mesh data, bakes palette LUTs, registers live sliders, and
   *        starts the spawn/spin/orbit animation chain.
   */
  void init() override {
    mobius_gen.init_storage(persistent_arena);
    setup_presets();

    std::span<PresetEntry<Params>> entries = preset_manager.get_entries();
    for (auto &entry : entries) {
      entry.params.palette = &palette;
    }

    params = preset_manager.get();
    baked_palettes[0].bake(persistent_arena, *params.palette);
    baked_palettes[1].bake(persistent_arena, *params.palette);

    register_animated_param("Copies", &params.num_copies, COPIES_MIN,
                            COPIES_MAX);
    register_animated_param("Radius", &params.offset_radius, RADIUS_MIN,
                            RADIUS_MAX);
    register_animated_param("Speed", &params.offset_speed, SPEED_MIN,
                            SPEED_MAX);
    register_animated_param("Warp", &params.warp_scale, WARP_MIN, WARP_MAX);
    register_animated_param("Alpha", &params.alpha, ALPHA_MIN, ALPHA_MAX);

    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    crossfade.overlap = CROSSFADE_OVERLAP;
    spawn_sprite();
    // Wrap the integrated phase to [0,1) so the orbit trig stays in precise range.
    timeline.add(0, Animation::Driver(orbit_phase, &params.offset_speed, 0.01f,
                                      /*wrap=*/true));
  }

  /**
   * @brief Advances the timeline, which drives all spawning and rendering.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    // Mirror live slider edits into the active sprite's snapshot; the previous
    // sprite's slot stays at the values it was spawned with.
    param_slots[active_bake] = params;
    {
      HS_PROFILE(db_timeline_step);
      timeline.step(canvas);
    }
  }

private:
  friend struct ::hs_test::effects_tests::DreamBallsWhiteBox;

  static constexpr float COPIES_MIN = 1.0f, COPIES_MAX = 20.0f;
  static constexpr float RADIUS_MIN = 0.0f, RADIUS_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float WARP_MIN = 0.0f, WARP_MAX = 5.0f;
  static constexpr float ALPHA_MIN = 0.0f, ALPHA_MAX = 1.0f;
  static constexpr size_t PRESET_COUNT = 5;

  /** Orbit phase in turns, wrapped to [0,1) by the live-speed Driver below. */
  float orbit_phase = 0.0f;
  int last_preset_idx =
      -1; /**< Last preset whose values were copied into params. */

  /** Per-vertex phase increment (radians) for the orbit stagger, so the surface
       ripples instead of pulsing in unison. */
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
    ArenaVector<Plot::Mesh::Edge>
        edges; /**< Unique edge list (topology is static). */
  };

  std::array<PresetData, PRESET_COUNT> loaded_presets;

  /**
   * @brief One preset's baked geometry element counts.
   */
  struct GeometryCounts {
    size_t vertices;   /**< Vertex count. */
    size_t face_slots; /**< Flattened face-index slots, i.e. 2*E. */
    size_t faces;      /**< Face count. */
  };
  /** Counts of the five Archimedean solids named in the preset table, in that
      table's order; setup_presets() re-checks every bake against them. */
  static constexpr std::array<GeometryCounts, PRESET_COUNT> PRESET_GEOMETRY = {
      {{24, 96, 26},
       {60, 240, 62},
       {48, 144, 26},
       {30, 120, 32},
       {24, 120, 38}}};
  static constexpr size_t PRESET_VERTICES =
      PRESET_GEOMETRY[0].vertices + PRESET_GEOMETRY[1].vertices +
      PRESET_GEOMETRY[2].vertices + PRESET_GEOMETRY[3].vertices +
      PRESET_GEOMETRY[4].vertices;
  static constexpr size_t PRESET_FACE_SLOTS =
      PRESET_GEOMETRY[0].face_slots + PRESET_GEOMETRY[1].face_slots +
      PRESET_GEOMETRY[2].face_slots + PRESET_GEOMETRY[3].face_slots +
      PRESET_GEOMETRY[4].face_slots;
  static constexpr size_t PRESET_FACES =
      PRESET_GEOMETRY[0].faces + PRESET_GEOMETRY[1].faces +
      PRESET_GEOMETRY[2].faces + PRESET_GEOMETRY[3].faces +
      PRESET_GEOMETRY[4].faces;

  FastNoiseLite noise;
  Timeline timeline;

  Orientation<> global_orientation;

  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  MobiusWarpTransformer<1> mobius_gen;
  static constexpr int SPRITE_LIFE = 320; /**< Visible frames per sprite. */
  static constexpr int FADE_WINDOW = 32;  /**< Fade-in/out length in frames. */
  /** Frames consecutive sprites coexist; 0 keeps every frame at a single mesh
      render. */
  static constexpr int CROSSFADE_OVERLAP = 0;
  // The two-slot ping-pong is safe only while at most two sprites overlap, i.e.
  // a sprite finishes before the spawn two periods later reuses its slot.
  static_assert(SPRITE_LIFE < 2 * (SPRITE_LIFE - CROSSFADE_OVERLAP),
                "DreamBalls ping-pong needs at most two overlapping sprites");

  /**
   * @brief Two baked LUTs, ping-ponged per spawn.
   * @details A spawn rebakes the inactive slot, so a rebake never lands on the
   *          slot a sprite is drawing from.
   */
  BakedPalette baked_palettes[2];
  // Persistent allocations: two ping-ponged palette LUTs plus every preset's
  // baked vertices, faces, face counts, tangent frames, and unique edge list.
  // All five presets are baked at init and live for the effect's whole life.
  static constexpr size_t FOOTPRINT_BYTES =
      2 * BakedPalette::required_arena_bytes() +
      PRESET_VERTICES * (sizeof(Vector) + sizeof(Tangent)) +
      PRESET_FACE_SLOTS * sizeof(uint16_t) + PRESET_FACES * sizeof(uint8_t) +
      (PRESET_FACE_SLOTS / 2) * sizeof(Plot::Mesh::Edge);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "DreamBalls persistent footprint exceeds the default partition");
  int active_bake = 0; /**< Slot of the current (most-recently baked) sprite;
                             the next spawn flips this before baking. */
  /**
   * @brief Per-sprite render-param snapshots, ping-ponged with baked_palettes.
   * @details Each sprite renders from its spawn-time slot; draw_frame() mirrors
   *          live sliders into the active slot only, so edits reach the drawing
   *          sprite without touching the previous snapshot.
   */
  Params param_slots[2];
  /** @brief Sprite hand-off crossfade; overlap is CROSSFADE_OVERLAP, set at
   *         init(). */
  Segue::Crossfade crossfade;

  GenerativePalette palette{GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                            BrightnessProfile::CUP};

  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{"rhombicuboctahedron", 18.0f, 0.3f, 0.4f, 0.3f, nullptr, 0.7f}},
      {{"rhombicosidodecahedron", 6.0f, 0.05f, 1.0f, 1.8f, nullptr, 0.7f}},
      {{"truncatedCuboctahedron", 6.0f, 0.16f, 1.0f, 2.0f, nullptr, 0.3f}},
      {{"icosidodecahedron", 10.0f, 0.16f, 1.0f, 0.5f, nullptr, 0.3f}},
      {{"snubCube", 4.534f, 0.153f, 2.025f, 0.0f, nullptr, 0.3f}},
  }};

  static constexpr bool preset_in_ranges(const Params &p) {
    return p.num_copies >= COPIES_MIN && p.num_copies <= COPIES_MAX &&
           p.offset_radius >= RADIUS_MIN && p.offset_radius <= RADIUS_MAX &&
           p.offset_speed >= SPEED_MIN && p.offset_speed <= SPEED_MAX &&
           p.warp_scale >= WARP_MIN && p.warp_scale <= WARP_MAX &&
           p.alpha >= ALPHA_MIN && p.alpha <= ALPHA_MAX;
  }

  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "a DreamBalls preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset");

  Presets<Params, PRESET_COUNT> preset_manager{PRESETS};

  /**
   * @brief Generates each preset's solid and bakes its geometry into the
   *        persistent arena once at init.
   * @details Bakes vertices, faces, tangent frames, and the unique edge list
   *          for every entry in the preset table.
   */
  HS_COLD_MEMBER void setup_presets() {
    const auto &entries = preset_manager.get_entries();

    int preset_idx = 0;
    for (const auto &entry : entries) {
      const auto &p = entry.params;
      auto &data = loaded_presets[preset_idx];

      hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        // Build the raw solid into scratch; only the deep copy below, bound to
        // the persistent target, outlives this generate() call.
        PolyMesh m =
            Solids::get_by_name(a, a, b, std::string_view(p.solid_name));

        data.mesh_state.vertices.bind(target, m.vertices.size());
        for (const auto &v : m.vertices) {
          data.mesh_state.vertices.push_back(v);
        }

        data.mesh_state.faces.bind(target, m.faces.size());
        data.mesh_state.face_counts.bind(target, m.face_counts.size());

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
        data.tangents.bind(target, data.mesh_state.vertices.size());
        for (const auto &v : data.mesh_state.vertices) {
          Vector axis = (std::abs(v.y) > 0.99f) ? X_AXIS : Y_AXIS;
          Vector u = cross(v, axis).normalized();
          Vector frame_v = cross(v, u).normalized();
          data.tangents.push_back({u, frame_v});
        }

        // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
        size_t edge_count = data.mesh_state.faces.size() / 2;
        data.edges.bind(target, edge_count);
        Plot::Mesh::extract_edges(data.mesh_state, data.edges);
        HS_CHECK(data.edges.size() == edge_count,
                 "DreamBalls edge extraction over/under-filled the edge bind");
      });

      // Pins the footprint table to what the generators actually produced.
      const GeometryCounts &counts = PRESET_GEOMETRY[preset_idx];
      HS_CHECK(data.mesh_state.vertices.size() == counts.vertices &&
                   data.mesh_state.faces.size() == counts.face_slots &&
                   data.mesh_state.face_counts.size() == counts.faces,
               "DreamBalls preset %d baked %d/%d/%d verts/face-slots/faces, "
               "footprint table says %d/%d/%d",
               preset_idx, (int)data.mesh_state.vertices.size(),
               (int)data.mesh_state.faces.size(),
               (int)data.mesh_state.face_counts.size(), (int)counts.vertices,
               (int)counts.face_slots, (int)counts.faces);

      preset_idx++;
    }
  }

  /**
   * @brief Spawns one fading sprite for the current preset and schedules the
   *        next spawn one period later.
   * @details Rebakes the inactive palette slot, arms a Mobius warp, and reseeds
   *          the live params only when the preset actually changes.
   */
  HS_COLD_MEMBER void spawn_sprite() {
    auto entries = preset_manager.get_entries();
    int safe_idx = static_cast<int>(preset_manager.current_index());

    // Reseed the slider-bound params only when the preset actually changes, so a
    // paused re-spawn of the same preset keeps the user's live edits.
    if (safe_idx != last_preset_idx) {
      params = entries[safe_idx].params;
      last_preset_idx = safe_idx;
    }
    // Ping-pong to the inactive slot so the rebake never lands on the previous
    // sprite's palette and params. draw_frame() keeps the active slot tracking
    // sliders.
    active_bake ^= 1;
    baked_palettes[active_bake].rebake(*params.palette);
    const int bake_slot = active_bake;
    param_slots[bake_slot] = params;

    auto draw_fn = [this, safe_idx, bake_slot](Canvas &canvas, float opacity) {
      HS_PROFILE(db_draw);
      if (!crossfade.visible(opacity))
        return;
      const auto &preset = loaded_presets[safe_idx];
      ScratchScope scratch_a_guard(scratch_arena_a);
      MeshState target_mesh;
      {
        HS_PROFILE(db_mesh_copy);
        MeshOps::transform(preset.mesh_state, target_mesh, scratch_arena_a);
      }

      // Slot captured at spawn, not active_bake: this sprite renders from its
      // own param + palette snapshot.
      this->draw_scene(canvas, param_slots[bake_slot],
                       crossfade.opacity(opacity), preset.mesh_state,
                       target_mesh, preset.tangents, preset.edges,
                       baked_palettes[bake_slot]);
    };

    const int period = crossfade.schedule(timeline, draw_fn, SPRITE_LIFE,
                                          FADE_WINDOW, &anims_paused);

    // Single-slot pool: free here only because the previous warp runs exactly
    // `period` frames and so completes earlier in the same step() that fires the
    // re-spawn timer below (events step in insertion order). The magnitude binds
    // to the live params rather than this spawn's frozen slot, so a dropped
    // spawn leaves the running warp tracking "Warp" instead of going inert.
    if (auto *warp = mobius_gen.spawn(0, params.warp_scale, period, false))
      warp->bind_scale(params.warp_scale);

    timeline.add_pausable(period,
                          Animation::PeriodicTimer(
                              0,
                              [this](Canvas &) {
                                preset_manager.next();
                                this->spawn_sprite();
                              },
                              false),
                          &anims_paused);
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

      float cos_a = fast_cosf(angle);
      float sin_a = fast_sinf(angle);

      Vector disp = v + (tan.u * cos_a + tan.v * sin_a) * r;
      target.vertices[i] = disp.normalized();
    }
  }

  /**
   * @brief Draws p.num_copies orbiting wireframe shells of the preset's solid.
   * @param canvas Render target.
   * @param p Live render params (copy count, radius, alpha, etc.).
   * @param opacity This sprite's crossfade opacity scaling edge alpha.
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
    HS_PROFILE(db_draw_scene);

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
      {
        HS_PROFILE(db_displace);
        update_displaced_mesh(base, target, tangents, p, offset);
      }

      {
        HS_PROFILE(db_warp_orient);
        for (size_t vi = 0; vi < target.vertices.size(); ++vi) {
          target.vertices[vi] = mobius_gen.transform(target.vertices[vi]);
          target.vertices[vi] = global_orientation.orient(target.vertices[vi]);
        }
      }

      {
        HS_PROFILE(db_mesh_plot);
        Plot::Mesh::draw<W, H>(filters, canvas, target, edges, fragment_shader);
      }
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

  /** @brief Live, slider-bound render parameters for the active preset. */
  Params params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(DreamBalls)
