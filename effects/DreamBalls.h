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

#include <algorithm>
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
  using BaseMesh = Solids::BaseMesh;

  /** @brief Selects automatic, source, or forced-medial weave topology. */
  enum class WeaveTopology : uint8_t {
    AUTOMATIC,
    ORIGINAL_WITH_DEFECTS,
    MEDIAL,
  };

  /**
   * @brief Live, slider-bound render parameters; also the per-preset value set.
   */
  struct Params {
    BaseMesh base_mesh = BaseMesh::TETRAHEDRON;
    WeaveTopology weave_topology = WeaveTopology::AUTOMATIC;
    float weave_gap = WEAVE_GAP_DEFAULT;
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
    loaded_solids = persistent_arena.allocate_n<SolidData>(SOLID_COUNT);
    for (size_t i = 0; i < SOLID_COUNT; ++i)
      new (&loaded_solids[i]) SolidData();
    setup_solids();

    blood_stream_composition.bind(&blood_stream_palette, &blood_stream_fade);
    std::span<PresetEntry<Params>> entries = preset_manager.get_entries();
    entries[0].params.palette = &blood_stream_falloff;
    entries[1].params.palette = &blood_stream_falloff;
    entries[2].params.palette = &Palettes::RICH_SUNSET;
    entries[3].params.palette = &Palettes::LAVENDER_LAKE;
    entries[4].params.palette = &Palettes::CORAL_BLUE;
    for (const auto &entry : entries) {
      HS_CHECK(entry.params.palette,
               "DreamBalls preset table left a null palette unpatched");
    }

    params = preset_manager.get();
    baked_palettes[0].bake(persistent_arena, *params.palette);
    baked_palettes[1].bake(persistent_arena, *params.palette);

    register_animated_param("Base Mesh", &params.base_mesh,
                            Solids::BASE_MESH_OPTIONS,
                            Solids::BASE_MESH_EXPORT_OPTIONS, SOLID_COUNT);
    register_animated_param(
        "Weave Topology", &params.weave_topology, WEAVE_TOPOLOGY_OPTIONS,
        WEAVE_TOPOLOGY_EXPORT_OPTIONS, std::size(WEAVE_TOPOLOGY_OPTIONS));
    register_animated_param("Weave Gap", &params.weave_gap, WEAVE_GAP_MIN,
                            WEAVE_GAP_MAX);
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
  static constexpr float WEAVE_GAP_MIN = 0.02f, WEAVE_GAP_MAX = 0.45f;
  static constexpr float WEAVE_GAP_DEFAULT = 0.18f;
  static constexpr size_t PRESET_COUNT = 5;
  static constexpr size_t SOLID_COUNT = Solids::BASE_MESH_COUNT;
  static constexpr const char *WEAVE_TOPOLOGY_OPTIONS[] = {
      "Automatic", "Original with defects", "Medial"};
  static constexpr const char *WEAVE_TOPOLOGY_EXPORT_OPTIONS[] = {
      "WeaveTopology::AUTOMATIC", "WeaveTopology::ORIGINAL_WITH_DEFECTS",
      "WeaveTopology::MEDIAL"};
  static_assert(SOLID_COUNT == std::size(Solids::simple_registry) +
                                   std::size(Solids::catalan_registry));
  static_assert(std::size(WEAVE_TOPOLOGY_OPTIONS) ==
                std::size(WEAVE_TOPOLOGY_EXPORT_OPTIONS));

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

  /** @brief Precomputed solid geometry baked into the persistent arena. */
  struct SolidData {
    MeshState mesh_state;          /**< Baked vertices and faces. */
    ArenaVector<Tangent> tangents; /**< Per-vertex tangent frames. */
    ArenaVector<Plot::Mesh::Edge>
        original_edges; /**< Unique edges of the source mesh. */
    ArenaVector<Plot::Mesh::Edge>
        automatic_edges; /**< Oriented source or medial edges. */
    ArenaVector<Plot::Mesh::Edge>
        medial_edges; /**< Forced-medial edges for four-regular solids. */
    bool four_regular = false; /**< Every source vertex has degree four. */
  };

  SolidData *loaded_solids = nullptr;
  static constexpr size_t MAX_SOLID_VERTICES = 120;
  static constexpr size_t MAX_SOLID_FACE_SLOTS = 360;
  static constexpr size_t MAX_SOLID_FACES = 120;
  static constexpr size_t MAX_SOLID_EDGES = MAX_SOLID_FACE_SLOTS / 2;
  static constexpr size_t FOUR_REGULAR_SOLID_COUNT = 5;

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
  static constexpr size_t FOOTPRINT_BYTES =
      2 * BakedPalette::required_arena_bytes() +
      SOLID_COUNT * sizeof(SolidData) +
      SOLID_COUNT * (MAX_SOLID_VERTICES * (sizeof(Vector) + sizeof(Tangent)) +
                     MAX_SOLID_FACE_SLOTS * sizeof(uint16_t) +
                     MAX_SOLID_FACES * sizeof(uint8_t) +
                     3 * MAX_SOLID_EDGES * sizeof(Plot::Mesh::Edge)) +
      FOUR_REGULAR_SOLID_COUNT * 2 * MAX_SOLID_EDGES * sizeof(Plot::Mesh::Edge);
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

  ProceduralPalette blood_stream_palette = Palettes::BLOOD_STREAM;
  AlphaFalloffShade blood_stream_fade{[](float t) { return 1.0f - t; }};
  StaticPalette<ProceduralPalette, Coords<>, Colors<AlphaFalloffShade>,
                /*Wrap=*/false>
      blood_stream_composition;
  PaletteFacade<decltype(blood_stream_composition)> blood_stream_falloff{
      &blood_stream_composition};

  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{BaseMesh::RHOMBICUBOCTAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 18.0f, 0.3f, 0.4f, 0.3f, nullptr, 0.7f}},
      {{BaseMesh::RHOMBICOSIDODECAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 6.0f, 0.05f, 1.0f, 1.8f, nullptr, 0.7f}},
      {{BaseMesh::TRUNCATED_CUBOCTAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 6.0f, 0.16f, 1.0f, 2.0f, nullptr, 0.3f}},
      {{BaseMesh::ICOSIDODECAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 10.0f, 0.16f, 1.0f, 0.5f, nullptr, 0.3f}},
      {{BaseMesh::SNUB_CUBE, WeaveTopology::AUTOMATIC, WEAVE_GAP_DEFAULT,
        4.534f, 0.153f, 2.025f, 0.0f, nullptr, 0.3f}},
  }};

  static constexpr bool preset_in_ranges(const Params &p) {
    return p.weave_gap >= WEAVE_GAP_MIN && p.weave_gap <= WEAVE_GAP_MAX &&
           p.num_copies >= COPIES_MIN && p.num_copies <= COPIES_MAX &&
           p.offset_radius >= RADIUS_MIN && p.offset_radius <= RADIUS_MAX &&
           p.offset_speed >= SPEED_MIN && p.offset_speed <= SPEED_MAX &&
           p.warp_scale >= WARP_MIN && p.warp_scale <= WARP_MAX &&
           p.alpha >= ALPHA_MIN && p.alpha <= ALPHA_MAX;
  }

  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "a DreamBalls preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset");

  Presets<Params, PRESET_COUNT> preset_manager{PRESETS};

  HS_COLD_MEMBER static bool source_is_four_regular(const MeshState &mesh,
                                                    Arena &scratch) {
    ScratchScope guard(scratch);
    uint8_t *degrees = scratch.allocate_n<uint8_t>(mesh.vertices.size());
    std::fill_n(degrees, mesh.vertices.size(), static_cast<uint8_t>(0));

    for (uint16_t vertex : mesh.faces) {
      HS_CHECK(vertex < mesh.vertices.size());
      degrees[vertex]++;
    }
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
      if (degrees[i] != 4)
        return false;
    return true;
  }

  HS_COLD_MEMBER static uint16_t
  find_edge_index(const ArenaVector<Plot::Mesh::Edge> &edges, uint16_t u,
                  uint16_t v) {
    for (size_t i = 0; i < edges.size(); ++i) {
      const auto &edge = edges[i];
      if ((edge.u == u && edge.v == v) || (edge.u == v && edge.v == u))
        return static_cast<uint16_t>(i);
    }
    HS_CHECK(false, "DreamBalls face edge missing from extracted topology");
    return 0;
  }

  HS_COLD_MEMBER static void
  build_four_regular_edges(const MeshState &mesh,
                           ArenaVector<Plot::Mesh::Edge> &edges,
                           Arena &scratch) {
    ScratchScope guard(scratch);
    HalfEdgeMesh half_edges(scratch, mesh);
    const size_t face_count = mesh.get_face_counts_size();
    uint8_t *colors = scratch.allocate_n<uint8_t>(face_count);
    uint16_t *queue = scratch.allocate_n<uint16_t>(face_count);
    std::fill_n(colors, face_count, static_cast<uint8_t>(UINT8_MAX));

    for (size_t seed = 0; seed < face_count; ++seed) {
      if (colors[seed] != UINT8_MAX)
        continue;
      size_t head = 0;
      size_t tail = 0;
      colors[seed] = 0;
      queue[tail++] = static_cast<uint16_t>(seed);

      while (head < tail) {
        const uint16_t face = queue[head++];
        const uint16_t start = half_edges.faces[face].half_edge;
        uint16_t current = start;
        do {
          const HalfEdge &he = half_edges.half_edges[current];
          HS_CHECK(he.pair != HE_NONE,
                   "DreamBalls weave requires a closed mesh");
          const uint16_t adjacent = half_edges.half_edges[he.pair].face;
          const uint8_t adjacent_color = colors[face] ^ 1;
          if (colors[adjacent] == UINT8_MAX) {
            colors[adjacent] = adjacent_color;
            queue[tail++] = adjacent;
          } else {
            HS_CHECK(colors[adjacent] == adjacent_color,
                     "four-regular mesh dual is not bipartite");
          }
          current = he.next;
        } while (current != start);
      }
    }

    const uint8_t *counts = mesh.get_face_counts_data();
    const uint16_t *faces = mesh.get_faces_data();
    size_t offset = 0;
    for (size_t face = 0; face < face_count; ++face) {
      const int count = counts[face];
      if (colors[face] == 0) {
        for (int i = 0; i < count; ++i)
          edges.push_back({faces[offset + i], faces[offset + (i + 1) % count]});
      }
      offset += count;
    }
  }

  HS_COLD_MEMBER static void
  build_medial_edges(const MeshState &mesh,
                     const ArenaVector<Plot::Mesh::Edge> &original_edges,
                     ArenaVector<Plot::Mesh::Edge> &medial_edges) {
    const uint8_t *counts = mesh.get_face_counts_data();
    const uint16_t *faces = mesh.get_faces_data();
    size_t offset = 0;
    for (size_t face = 0; face < mesh.get_face_counts_size(); ++face) {
      const int count = counts[face];
      for (int i = 0; i < count; ++i) {
        const uint16_t a = faces[offset + i];
        const uint16_t b = faces[offset + (i + 1) % count];
        const uint16_t c = faces[offset + (i + 2) % count];
        medial_edges.push_back({find_edge_index(original_edges, a, b),
                                find_edge_index(original_edges, b, c)});
      }
      offset += count;
    }
  }

  static float under_gap_alpha(float edge_t, float gap) {
    const float x = hs::clamp((1.0f - edge_t) / gap, 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
  }

  HS_COLD_MEMBER static Tangent tangent_frame(const Vector &normal) {
    const Vector axis = std::abs(normal.y) > 0.99f ? X_AXIS : Y_AXIS;
    const Vector u = cross(normal, axis).normalized();
    return {u, cross(normal, u).normalized()};
  }

  static Vector parallel_transport(const Vector &from, const Vector &to,
                                   const Vector &tangent) {
    const float denominator = 1.0f + dot(from, to);
    HS_CHECK(denominator > math::TOLERANCE,
             "DreamBalls weave edge has antipodal endpoints");
    return tangent - (from + to) * (dot(tangent, to) / denominator);
  }

  HS_COLD_MEMBER static Vector woven_vertex(const SolidData &solid, bool medial,
                                            size_t vertex) {
    if (!medial)
      return solid.mesh_state.vertices[vertex];
    const auto &edge = solid.original_edges[vertex];
    return normalized_or(solid.mesh_state.vertices[edge.u] +
                             solid.mesh_state.vertices[edge.v],
                         solid.mesh_state.vertices[edge.u]);
  }

  HS_COLD_MEMBER static void
  prepare_woven_buffers(const SolidData &solid, bool medial,
                        const ArenaVector<Plot::Mesh::Edge> &edges,
                        ArenaVector<Vector> &base_vertices,
                        ArenaVector<Vector> &frame_u,
                        ArenaVector<Vector> &offsets, MeshState &framed_mesh,
                        ArenaVector<Plot::Mesh::Edge> &framed_edges) {
    const size_t vertex_count =
        medial ? solid.original_edges.size() : solid.mesh_state.vertices.size();
    const size_t edge_count = edges.size();

    base_vertices.bind(scratch_arena_a, vertex_count);
    frame_u.bind(scratch_arena_a, vertex_count);
    for (size_t vertex = 0; vertex < vertex_count; ++vertex) {
      const Vector base = woven_vertex(solid, medial, vertex);
      base_vertices.push_back(base);
      frame_u.push_back(tangent_frame(base).u);
    }

    offsets.bind(scratch_arena_a, vertex_count);
    for (size_t vertex = 0; vertex < vertex_count; ++vertex)
      offsets.push_back(Vector());

    framed_mesh.vertices.bind(scratch_arena_a, vertex_count + edge_count);
    for (size_t vertex = 0; vertex < vertex_count + edge_count; ++vertex)
      framed_mesh.vertices.push_back(Vector());

    framed_edges.bind(scratch_arena_b, edge_count);
    for (size_t edge = 0; edge < edge_count; ++edge)
      framed_edges.push_back(
          {edges[edge].u, static_cast<uint16_t>(vertex_count + edge)});
  }

  /**
   * @brief Generates each selectable solid and bakes its geometry into the
   *        persistent arena once at init.
   * @details Bakes vertices, faces, tangent frames, and the unique edge list.
   */
  HS_COLD_MEMBER void setup_solids() {
    for (size_t solid_idx = 0; solid_idx < SOLID_COUNT; ++solid_idx) {
      const Solids::Entry &entry = Solids::get_entry(solid_idx);
      auto &data = loaded_solids[solid_idx];

      hs::generate(persistent_arena, [&](Arena &target, Arena &a,
                                         Arena &b) HS_COLD_MEMBER {
        PolyMesh m = Solids::finalize_solid(entry.generate(a, b), a);

        HS_CHECK(m.vertices.size() <= MAX_SOLID_VERTICES &&
                     m.faces.size() <= MAX_SOLID_FACE_SLOTS &&
                     m.face_counts.size() <= MAX_SOLID_FACES,
                 "DreamBalls selectable solid exceeds footprint bounds");

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

        data.tangents.bind(target, data.mesh_state.vertices.size());
        for (const auto &v : data.mesh_state.vertices)
          data.tangents.push_back(tangent_frame(v));

        // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
        size_t edge_count = data.mesh_state.faces.size() / 2;
        data.original_edges.bind(target, edge_count);
        Plot::Mesh::extract_edges(data.mesh_state, data.original_edges);
        HS_CHECK(data.original_edges.size() == edge_count,
                 "DreamBalls edge extraction over/under-filled the edge bind");

        data.four_regular = source_is_four_regular(data.mesh_state, b);
        const size_t automatic_edge_count =
            data.four_regular ? edge_count : 2 * edge_count;
        data.automatic_edges.bind(target, automatic_edge_count);
        if (data.four_regular) {
          build_four_regular_edges(data.mesh_state, data.automatic_edges, b);
          data.medial_edges.bind(target, 2 * edge_count);
          build_medial_edges(data.mesh_state, data.original_edges,
                             data.medial_edges);
          HS_CHECK(data.medial_edges.size() == 2 * edge_count,
                   "DreamBalls medial topology edge count mismatch");
        } else {
          build_medial_edges(data.mesh_state, data.original_edges,
                             data.automatic_edges);
        }
        HS_CHECK(data.automatic_edges.size() == automatic_edge_count,
                 "DreamBalls automatic topology edge count mismatch");
      });
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

    auto draw_fn = [this, bake_slot](Canvas &canvas, float opacity) {
      HS_PROFILE(db_draw);
      if (!crossfade.visible(opacity))
        return;
      const Params &sprite_params = param_slots[bake_slot];
      const auto &solid =
          loaded_solids[static_cast<size_t>(sprite_params.base_mesh)];
      ScratchScope scratch_a_guard(scratch_arena_a);

      // Slot captured at spawn, not active_bake: this sprite renders from its
      // own param + palette snapshot.
      this->draw_scene(canvas, sprite_params, crossfade.opacity(opacity), solid,
                       baked_palettes[bake_slot]);
    };

    const int period = crossfade.schedule(timeline, draw_fn, SPRITE_LIFE,
                                          FADE_WINDOW, &anims_paused);

    // Single-slot pool: free here only because the previous warp runs exactly
    // `period` frames and so completes earlier in the same step() that fires the
    // re-spawn timer below (events step in insertion order).
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
  HS_COLD_MEMBER void
  update_displaced_mesh(const MeshState &base, MeshState &target,
                        const ArenaVector<Tangent> &tangents, const Params &p,
                        float angle_offset) {
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

  void draw_woven_scene(Canvas &canvas, const Params &p, const SolidData &solid,
                        bool medial, FragmentShaderFn fragment_shader) {
    const auto &edges = medial && solid.four_regular ? solid.medial_edges
                                                     : solid.automatic_edges;
    const size_t vertex_count =
        medial ? solid.original_edges.size() : solid.mesh_state.vertices.size();
    const size_t edge_count = edges.size();

    ArenaVector<Vector> base_vertices;
    ArenaVector<Vector> frame_u;
    ArenaVector<Vector> offsets;
    MeshState framed_mesh;
    ScratchScope scratch_b_guard(scratch_arena_b);
    ArenaVector<Plot::Mesh::Edge> framed_edges;
    prepare_woven_buffers(solid, medial, edges, base_vertices, frame_u, offsets,
                          framed_mesh, framed_edges);

    const int num_copies_raw = static_cast<int>(p.num_copies);
    const int num_copies = num_copies_raw < 1 ? 1 : num_copies_raw;
    for (int copy = 0; copy < num_copies; ++copy) {
      const float copy_phase =
          orbit_phase * 2 * PI_F +
          (static_cast<float>(copy) / num_copies) * 2 * PI_F;
      {
        HS_PROFILE(db_displace);
        for (size_t vertex = 0; vertex < vertex_count; ++vertex) {
          const Vector &base = base_vertices[vertex];
          const Vector &u = frame_u[vertex];
          const Vector v = cross(base, u);
          const float phase =
              copy_phase + static_cast<float>(vertex) * VERTEX_PHASE_STAGGER;
          offsets[vertex] = u * fast_cosf(phase) + v * fast_sinf(phase);
          framed_mesh.vertices[vertex] =
              normalized_or(base + offsets[vertex] * p.offset_radius, base);
        }

        for (size_t edge_index = 0; edge_index < edge_count; ++edge_index) {
          const auto &edge = edges[edge_index];
          const Vector &from = base_vertices[edge.u];
          const Vector &to = base_vertices[edge.v];
          const Vector head_offset =
              parallel_transport(from, to, offsets[edge.u]);
          framed_mesh.vertices[vertex_count + edge_index] =
              normalized_or(to + head_offset * p.offset_radius, to);
        }
      }

      {
        HS_PROFILE(db_warp_orient);
        for (auto &vertex : framed_mesh.vertices) {
          vertex = mobius_gen.transform(vertex);
          vertex = global_orientation.orient(vertex);
        }
      }

      {
        HS_PROFILE(db_mesh_plot);
        Plot::Mesh::draw<W, H>(filters, canvas, framed_mesh, framed_edges,
                               fragment_shader);
      }
    }
  }

  HS_COLD_MEMBER void draw_original_scene(Canvas &canvas, const Params &p,
                                          const SolidData &solid,
                                          FragmentShaderFn fragment_shader) {
    MeshState target;
    {
      HS_PROFILE(db_mesh_copy);
      MeshOps::transform(solid.mesh_state, target, scratch_arena_a);
    }

    const int num_copies_raw = static_cast<int>(p.num_copies);
    const int num_copies = num_copies_raw < 1 ? 1 : num_copies_raw;
    for (int copy = 0; copy < num_copies; ++copy) {
      const float offset = (static_cast<float>(copy) / num_copies) * 2 * PI_F;
      {
        HS_PROFILE(db_displace);
        update_displaced_mesh(solid.mesh_state, target, solid.tangents, p,
                              offset);
      }

      {
        HS_PROFILE(db_warp_orient);
        for (auto &vertex : target.vertices) {
          vertex = mobius_gen.transform(vertex);
          vertex = global_orientation.orient(vertex);
        }
      }

      {
        HS_PROFILE(db_mesh_plot);
        Plot::Mesh::draw<W, H>(filters, canvas, target, solid.original_edges,
                               fragment_shader);
      }
    }
  }

  /**
   * @brief Draws p.num_copies orbiting wireframe shells of the preset's solid.
   * @param canvas Render target.
   * @param p Live render params (copy count, radius, alpha, etc.).
   * @param opacity This sprite's crossfade opacity scaling edge alpha.
   * @param solid Source geometry and precomputed weave topology.
   * @param baked Baked palette LUT supplying edge colors.
   * @details Automatic mode parallel-transports each crossing's outgoing frame
   *          to the hidden end of every incoming edge. Defect mode retains the
   *          source mesh's shared vertex framing.
   */
  void draw_scene(Canvas &canvas, const Params &p, float opacity,
                  const SolidData &solid, const BakedPalette &baked) {
    HS_PROFILE(db_draw_scene);

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      Color4 c = baked.get(f.v0);
      c.alpha *= p.alpha * opacity * under_gap_alpha(f.v0, p.weave_gap);
      f.color = c;
    };

    if (p.weave_topology == WeaveTopology::ORIGINAL_WITH_DEFECTS) {
      draw_original_scene(canvas, p, solid, fragment_shader);
    } else {
      const bool medial =
          p.weave_topology == WeaveTopology::MEDIAL || !solid.four_regular;
      draw_woven_scene(canvas, p, solid, medial, fragment_shader);
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
