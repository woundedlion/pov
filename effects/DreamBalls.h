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

#include "core/control/choreography.h"
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
 * @brief DreamBalls' live, slider-bound render parameters; also the per-preset
 *        value set. The per-preset palette lives beside the effect's preset
 *        table, not here, so parameter snapshots stay pointer-free.
 */
struct DreamBallsParams {
  using BaseMesh = Solids::BaseMesh;

  /** @brief Selects automatic, source, or forced-medial weave topology. */
  enum class WeaveTopology : uint8_t {
    AUTOMATIC,
    ORIGINAL_WITH_DEFECTS,
    MEDIAL,
  };

  BaseMesh base_mesh = BaseMesh::TETRAHEDRON;
  WeaveTopology weave_topology = WeaveTopology::AUTOMATIC;
  float weave_gap = 0.18f;
  float num_copies = 1.0f;
  float offset_radius = 0.0f;
  float offset_speed = 0.0f;
  float warp_scale = 0.0f;
  float alpha = 0.0f;
};

/**
 * @brief Orbiting copies of a polyhedral wireframe, cycling through ten
 *        solid presets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each vertex is displaced in its own tangent plane and
 *          Mobius-warped before being re-projected onto the sphere.
 */
template <int W, int H>
class DreamBalls
    : public ChoreographedEffect<DreamBalls<W, H>, DreamBallsParams> {
  using Choreography = ChoreographedEffect<DreamBalls<W, H>, DreamBallsParams>;
  friend Choreography;

public:
  using Params = DreamBallsParams;
  using BaseMesh = Solids::BaseMesh;
  using WeaveTopology = Params::WeaveTopology;

  /** Preset policy: every origin snaps. The sprite chain owns the automatic
      cadence — a preset advances at each sprite hand-off, so the dwell
      countdown never runs and step_choreography() is never called. */
  static constexpr Segue::Preset::Snap PRESET_SEGUE{};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** Bookkeeping only (see PRESET_SEGUE); mirrors the sprite hand-off period. */
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;

  static bool valid_params(const Params &p) { return preset_in_ranges(p); }

  /**
   * @brief Constructs the effect with the anti-alias screen filter.
   * @details The Mobius generator hangs off the timeline so its warp animates
   *          in step.
   */
  HS_COLD_MEMBER DreamBalls()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        filters(Filter::Screen::AntiAlias<W, H>()), mobius_gen(this->timeline) {
  }

  /**
   * @brief Builds mesh data, bakes palette LUTs, registers live sliders, and
   *        starts the spawn/spin/orbit animation chain.
   */
  void init() override {
    begin_choreography();
    mobius_gen.init_storage(persistent_arena);
    loaded_solids = persistent_arena.make_n<SolidData>(SOLID_COUNT);
    setup_solids();

    blood_stream_composition.bind(&blood_stream_palette, &blood_stream_fade);
    preset_palettes = {&blood_stream_falloff,  &blood_stream_falloff,
                       &Palettes::RICH_SUNSET, &Palettes::LAVENDER_LAKE,
                       &Palettes::CORAL_BLUE,  &Palettes::CORAL_BLUE,
                       &Palettes::CORAL_BLUE,  &Palettes::CORAL_BLUE,
                       &Palettes::CORAL_BLUE,  &Palettes::CORAL_BLUE};
    live_palette = preset_palettes[0];
    baked_palettes[0].bake(persistent_arena, *live_palette);
    baked_palettes[1].bake(persistent_arena, *live_palette);

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

    // Ahead of every sprite, so a frame draws the warp and orbit phase it just
    // advanced rather than the previous frame's.
    arm_warp();
    timeline.add(0, Animation::PeriodicTimer(
                        160, [this](Canvas &) { this->spin_slices(); }, true));
    timeline.add(9, Animation::RandomWalk<W>(
                        global_orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
    // Wrap the integrated phase to [0,1) so the orbit trig stays in precise range.
    timeline.add(0, Animation::Driver(orbit_phase, &params.offset_speed, 0.01f,
                                      /*wrap=*/true));

    crossfade.overlap = CROSSFADE_OVERLAP;
    spawn_sprite();
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
  using Choreography::anims_paused;
  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::register_animated_param;
  using Choreography::timeline;

  /**
   * @brief Tracks the committed preset's palette; a non-automatic change also
   *        retargets the live sprite.
   * @details An automatic advance leaves the rebake to the spawn_sprite() call
   * that follows it in the hand-off timer, so the new palette lands on the
   * fresh sprite's slot rather than the finished one's.
   */
  HS_FLASH_MEMBER void
  preset_changed(const Effect::PresetChange &change) override {
    live_palette = preset_palettes[change.to];
    if (change.origin != Effect::PresetChangeOrigin::AUTOMATIC) {
      baked_palettes[active_bake].rebake(*live_palette);
      param_slots[active_bake] = params;
    }
#ifdef HS_PROFILE_ENABLE
    hs::log("Preset: %u/%u", static_cast<unsigned>(change.to + 1),
            static_cast<unsigned>(this->getPresetCount()));
#endif
  }

  friend struct ::hs_test::effects_tests::DreamBallsWhiteBox;

  static constexpr float COPIES_MIN = 1.0f, COPIES_MAX = 20.0f;
  static constexpr float RADIUS_MIN = 0.0f, RADIUS_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float WARP_MIN = 0.0f, WARP_MAX = 5.0f;
  static constexpr float ALPHA_MIN = 0.0f, ALPHA_MAX = 1.0f;
  static constexpr float WEAVE_GAP_MIN = 0.02f, WEAVE_GAP_MAX = 0.45f;
  static constexpr float WEAVE_GAP_DEFAULT = 0.18f;
  static constexpr size_t PRESET_COUNT = 10;
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

  /** Per-vertex phase increment (radians) for the orbit stagger, so the surface
       ripples instead of pulsing in unison. */
  static constexpr float VERTEX_PHASE_STAGGER = 0.1f;

  /** @brief Precomputed solid geometry baked into the persistent arena. */
  struct SolidData {
    MeshState mesh_state; /**< Baked vertices and faces. */
    /** First tangent-basis vector per vertex; the second is cross(vertex, u). */
    ArenaVector<Vector> tangent_u;
    ArenaVector<Plot::Mesh::Edge>
        original_edges; /**< Unique edges of the source mesh. */
    ArenaVector<Plot::Mesh::Edge>
        automatic_edges; /**< Oriented source or medial edges. */
    ArenaVector<Plot::Mesh::Edge>
        medial_edges; /**< Forced-medial edges for four-regular solids. */
    bool four_regular = false; /**< Every source vertex has degree four. */
  };

  SolidData *loaded_solids = nullptr;
  /** @brief Registry solids the medial-edge term of FOOTPRINT_BYTES covers. */
  static constexpr size_t FOUR_REGULAR_SOLID_COUNT = 5;
  static_assert(Solids::MAX_SOLID_VERTICES <=
                    static_cast<size_t>(Plot::Mesh::DEDUP_CAPACITY),
                "a baked solid must fit the wireframe edge-dedup bitset");

  FastNoiseLite noise;

  Orientation<> global_orientation;

  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  MobiusWarpTransformer<1> mobius_gen;
  static constexpr int SPRITE_LIFE = 320; /**< Visible frames per sprite. */
  static constexpr int FADE_WINDOW = 32;  /**< Fade-in/out length in frames. */
  /** Frames consecutive sprites coexist; 0 keeps every frame at a single mesh
      render. */
  static constexpr int CROSSFADE_OVERLAP = 0;
  /** Warp cycle length: the pinned warp repeats in lockstep with the sprite
      hand-off, so each sprite spans exactly one warp. */
  static constexpr int WARP_PERIOD = SPRITE_LIFE - CROSSFADE_OVERLAP;
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
      SOLID_COUNT * (Solids::MAX_SOLID_VERTICES * 2 * sizeof(Vector) +
                     Solids::MAX_SOLID_FACE_SLOTS * sizeof(uint16_t) +
                     Solids::MAX_SOLID_FACES * sizeof(uint8_t) +
                     3 * Solids::MAX_SOLID_EDGES * sizeof(Plot::Mesh::Edge)) +
      FOUR_REGULAR_SOLID_COUNT * 2 * Solids::MAX_SOLID_EDGES *
          sizeof(Plot::Mesh::Edge);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "DreamBalls persistent footprint exceeds the default partition");
  // Medial topology turns every source edge into a woven vertex, so the vertex
  // bound is MAX_SOLID_EDGES rather than MAX_SOLID_VERTICES.
  static constexpr size_t WOVEN_VERTEX_BOUND =
      Solids::MAX_SOLID_EDGES > Solids::MAX_SOLID_VERTICES
          ? Solids::MAX_SOLID_EDGES
          : Solids::MAX_SOLID_VERTICES;
  // draw_woven_scene stages three per-vertex buffers plus the framed
  // vertex + edge-head mesh in scratch_a, all live across Plot::Mesh::draw's
  // per-edge fragment scope and rasterize's sub-step cache.
  static constexpr size_t SCRATCH_A_PEAK_BYTES =
      (4 * WOVEN_VERTEX_BOUND + 2 * Solids::MAX_SOLID_EDGES) * sizeof(Vector) +
      Plot::Mesh::EDGE_MAX_POINTS * sizeof(Fragment) +
      Plot::rasterize_scratch_a_bytes<W>();
  static_assert(
      SCRATCH_A_PEAK_BYTES <= DEFAULT_SCRATCH_A_SIZE,
      "DreamBalls woven staging exceeds the default scratch_a budget; "
      "retune the solid bounds or carve a larger scratch arena");
  // The framed edge list and the per-vertex weave-start owners live in
  // scratch_b across the same mesh draw.
  static constexpr size_t SCRATCH_B_PEAK_BYTES =
      2 * Solids::MAX_SOLID_EDGES * sizeof(Plot::Mesh::Edge) +
      WOVEN_VERTEX_BOUND * sizeof(uint16_t);
  static_assert(
      SCRATCH_B_PEAK_BYTES <= DEFAULT_SCRATCH_B_SIZE,
      "DreamBalls woven staging exceeds the default scratch_b budget; "
      "retune the solid bounds or carve a larger scratch arena");
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
        WEAVE_GAP_DEFAULT, 18.0f, 0.3f, 0.4f, 0.3f, 0.7f}},
      {{BaseMesh::RHOMBICOSIDODECAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 6.0f, 0.05f, 1.0f, 1.8f, 0.7f}},
      {{BaseMesh::TRUNCATED_CUBOCTAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 6.0f, 0.16f, 1.0f, 2.0f, 0.3f}},
      {{BaseMesh::ICOSIDODECAHEDRON, WeaveTopology::AUTOMATIC,
        WEAVE_GAP_DEFAULT, 10.0f, 0.16f, 1.0f, 0.5f, 0.3f}},
      {{BaseMesh::SNUB_CUBE, WeaveTopology::AUTOMATIC, WEAVE_GAP_DEFAULT,
        4.534f, 0.153f, 2.025f, 0.0f, 0.3f}},
      {{BaseMesh::TRUNCATED_DODECAHEDRON, WeaveTopology::AUTOMATIC, 0.18f,
        4.515f, 0.179f, 1.89f, 1.535f, 0.7f}},
      {{BaseMesh::TRIAKIS_ICOSAHEDRON, WeaveTopology::AUTOMATIC, 0.18f, 4.515f,
        0.131f, 1.89f, 1.535f, 0.7f}},
      {{BaseMesh::TRIAKIS_ICOSAHEDRON, WeaveTopology::AUTOMATIC, 0.18f, 6.0f,
        0.078f, 1.0f, 0.0f, 0.3f}},
      {{BaseMesh::DISDYAKIS_TRIACONTAHEDRON, WeaveTopology::AUTOMATIC, 0.18f,
        6.0f, 0.03f, 1.0f, 1.795f, 0.3f}},
      {{BaseMesh::TRIAKIS_ICOSAHEDRON, WeaveTopology::AUTOMATIC, 0.18f, 6.0f,
        0.03f, 1.0f, 1.795f, 0.3f}},
  }};

  /** @brief Per-preset palette, patched at init(); kept beside PRESETS rather
   *  than inside Params so parameter snapshots carry no pointer. */
  std::array<const Palette *, PRESET_COUNT> preset_palettes = {};
  /** @brief The committed preset's palette; spawns bake from this. */
  const Palette *live_palette = nullptr;

  static constexpr bool preset_in_ranges(const Params &p) {
    return static_cast<size_t>(p.base_mesh) < SOLID_COUNT &&
           static_cast<size_t>(p.weave_topology) <
               std::size(WEAVE_TOPOLOGY_OPTIONS) &&
           p.weave_gap >= WEAVE_GAP_MIN && p.weave_gap <= WEAVE_GAP_MAX &&
           p.num_copies >= COPIES_MIN && p.num_copies <= COPIES_MAX &&
           p.offset_radius >= RADIUS_MIN && p.offset_radius <= RADIUS_MAX &&
           p.offset_speed >= SPEED_MIN && p.offset_speed <= SPEED_MAX &&
           p.warp_scale >= WARP_MIN && p.warp_scale <= WARP_MAX &&
           p.alpha >= ALPHA_MIN && p.alpha <= ALPHA_MAX;
  }

  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "a DreamBalls preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset");

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

  static float under_gap_alpha(float edge_t, float gap) {
    return cubic_kernel((1.0f - edge_t) / gap);
  }

  HS_FLASH_MEMBER static void
  assign_woven_start_owners(const ArenaVector<Plot::Mesh::Edge> &edges,
                            uint16_t *owners, size_t vertex_count) {
    std::fill_n(owners, vertex_count, UINT16_MAX);
    for (size_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
      const uint16_t vertex = edges[edge_index].u;
      HS_CHECK(vertex < vertex_count);
      if (owners[vertex] == UINT16_MAX)
        owners[vertex] = static_cast<uint16_t>(edge_index);
    }
  }

  static bool
  owns_woven_start_sample(const ArenaVector<Plot::Mesh::Edge> &edges,
                          const uint16_t *owners, size_t edge_index) {
    HS_CHECK(edge_index < edges.size());
    return owners[edges[edge_index].u] == edge_index;
  }

  HS_FLASH_MEMBER static Vector woven_vertex(const SolidData &solid,
                                             bool medial, size_t vertex) {
    if (!medial)
      return solid.mesh_state.vertices[vertex];
    const auto &edge = solid.original_edges[vertex];
    return normalized_or(solid.mesh_state.vertices[edge.u] +
                             solid.mesh_state.vertices[edge.v],
                         solid.mesh_state.vertices[edge.u]);
  }

  HS_FLASH_MEMBER static void
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
      // Medial vertices are edge midpoints, absent from the baked table.
      frame_u.push_back(medial ? tangent_axis(base) : solid.tangent_u[vertex]);
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
    size_t four_regular_solids = 0;
    for (size_t solid_idx = 0; solid_idx < SOLID_COUNT; ++solid_idx) {
      const Solids::Entry &entry = Solids::get_entry(solid_idx);
      auto &data = loaded_solids[solid_idx];

      hs::generate(persistent_arena, [&](Arena &target, Arena &a,
                                         Arena &b) HS_COLD_MEMBER {
        PolyMesh m = Solids::finalize_solid(entry.generate(a, b), a);

        HS_CHECK(m.vertices.size() <= Solids::MAX_SOLID_VERTICES &&
                     m.faces.size() <= Solids::MAX_SOLID_FACE_SLOTS &&
                     m.face_counts.size() <= Solids::MAX_SOLID_FACES,
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

        data.tangent_u.bind(target, data.mesh_state.vertices.size());
        for (const auto &v : data.mesh_state.vertices)
          data.tangent_u.push_back(tangent_axis(v));

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
          ++four_regular_solids;
          HS_CHECK(four_regular_solids <= FOUR_REGULAR_SOLID_COUNT,
                   "DreamBalls registry holds more four-regular solids than "
                   "the persistent footprint budgets");
          Plot::Mesh::extract_four_regular_edges(data.mesh_state,
                                                 data.automatic_edges, b);
          data.medial_edges.bind(target, 2 * edge_count);
          Plot::Mesh::extract_medial_edges(data.mesh_state, data.original_edges,
                                           data.medial_edges);
          HS_CHECK(data.medial_edges.size() == 2 * edge_count,
                   "DreamBalls medial topology edge count mismatch");
        } else {
          Plot::Mesh::extract_medial_edges(data.mesh_state, data.original_edges,
                                           data.automatic_edges);
        }
        HS_CHECK(data.automatic_edges.size() == automatic_edge_count,
                 "DreamBalls automatic topology edge count mismatch");
      });
    }
  }

  /**
   * @brief Arms the one repeating Mobius warp every sprite renders through.
   * @details Pinned: it never completes, so the pool slot is never contended and
   *          the warp leads every sprite in the event order.
   */
  HS_COLD_MEMBER void arm_warp() {
    auto *warp =
        mobius_gen.spawn_pinned(0, params.warp_scale, WARP_PERIOD, true);
    HS_CHECK(warp, "DreamBalls: pinned warp spawn must succeed");
    warp->bind_scale(params.warp_scale);
  }

  /**
   * @brief Spawns one fading sprite for the current preset and schedules the
   *        next spawn one period later.
   * @details The preset's params were adopted when it was committed; the spawn
   *          rebakes the inactive palette slot for the fresh sprite.
   */
  HS_COLD_MEMBER void spawn_sprite() {
    // Ping-pong to the inactive slot so the rebake never lands on the previous
    // sprite's palette and params. draw_frame() keeps the active slot tracking
    // sliders.
    active_bake ^= 1;
    baked_palettes[active_bake].rebake(*live_palette);
    const int bake_slot = active_bake;
    param_slots[bake_slot] = params;

    auto draw_fn = [this, bake_slot](Canvas &canvas, float opacity) {
      HS_PROFILE(db_draw);
      if (!crossfade.visible(opacity))
        return;
      const Params &sprite_params = param_slots[bake_slot];
      // Alpha below one slider LSB: skip the whole weave build and plot.
      if (sprite_params.alpha < MIN_VISIBLE_ALPHA)
        return;
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

    HS_CHECK(period == WARP_PERIOD,
             "DreamBalls: sprite hand-off drifted off the pinned warp's cycle");

    timeline.add_pausable(period,
                          Animation::PeriodicTimer(
                              0,
                              [this](Canvas &) {
                                const bool advanced = this->advancePreset();
                                HS_CHECK(advanced);
                                this->spawn_sprite();
                              },
                              false),
                          &anims_paused);
  }

  /**
   * @brief Unit orbit offset in the tangent plane at a vertex.
   * @param base Unit vertex direction whose tangent plane carries the orbit.
   * @param u First tangent-basis vector; cross(base, u) spans the second.
   * @param phase Orbit angle in radians.
   * @return The unit offset direction at that phase.
   */
  __attribute__((always_inline)) static Vector
  tangent_orbit_offset(const Vector &base, const Vector &u, float phase) {
    return u * fast_cosf(phase) + cross(base, u) * fast_sinf(phase);
  }

  /**
   * @brief Orbits each vertex in a small circle within its own tangent plane,
   *        then re-projects onto the unit sphere.
   * @param base Source mesh whose vertices are orbited.
   * @param target Destination mesh receiving the displaced vertices.
   * @param tangent_u Per-vertex first tangent-basis vector; the second spanning
   *        vector is cross(vertex, u).
   * @param p Render params; p.offset_radius is the orbit radius.
   * @param angle_offset Per-copy phase offset in radians.
   * @details The per-vertex phase (i * VERTEX_PHASE_STAGGER) staggers the
   *          orbits so the surface ripples.
   */
  HS_COLD_MEMBER void
  update_displaced_mesh(const MeshState &base, MeshState &target,
                        const ArenaVector<Vector> &tangent_u, const Params &p,
                        float angle_offset) {
    size_t count = base.vertices.size();
    float r = p.offset_radius;

    // MeshOps::transform pre-sizes target.vertices to match base; the indexed
    // writes below rely on it.
    HS_CHECK(target.vertices.size() == base.vertices.size(),
             "DreamBalls: displaced-mesh target not pre-sized to base");

    for (size_t i = 0; i < count; ++i) {
      const Vector &v = base.vertices[i];

      // orbit_phase is a fraction of a turn; scale to radians.
      float phase = i * VERTEX_PHASE_STAGGER;
      float angle = orbit_phase * 2 * PI_F + phase + angle_offset;

      const Vector disp = v + tangent_orbit_offset(v, tangent_u[i], angle) * r;
      target.vertices[i] = normalized_or(disp, v);
    }
  }

  void draw_woven_scene(Canvas &canvas, const Params &p, const SolidData &solid,
                        bool medial, float opacity, const BakedPalette &baked) {
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

    uint16_t *start_owners = scratch_arena_b.allocate_n<uint16_t>(vertex_count);
    assign_woven_start_owners(edges, start_owners, vertex_count);

    auto woven_shader = [&](const Vector &, Fragment &f) {
      Color4 c = baked.get(f.v0);
      c.alpha *= p.alpha * opacity * under_gap_alpha(f.v0, p.weave_gap);
      if (f.v0 == 0.0f) {
        const size_t edge_index = static_cast<size_t>(f.v2);
        if (!owns_woven_start_sample(edges, start_owners, edge_index))
          c.alpha = 0.0f;
      }
      f.color = c;
    };

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
          const float phase =
              copy_phase + static_cast<float>(vertex) * VERTEX_PHASE_STAGGER;
          offsets[vertex] = tangent_orbit_offset(base, frame_u[vertex], phase);
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
                               woven_shader);
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
        update_displaced_mesh(solid.mesh_state, target, solid.tangent_u, p,
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

    if (p.weave_topology == WeaveTopology::ORIGINAL_WITH_DEFECTS) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        Color4 c = baked.get(f.v0);
        c.alpha *= p.alpha * opacity * under_gap_alpha(f.v0, p.weave_gap);
        f.color = c;
      };
      draw_original_scene(canvas, p, solid, fragment_shader);
    } else {
      const bool medial =
          p.weave_topology == WeaveTopology::MEDIAL || !solid.four_regular;
      draw_woven_scene(canvas, p, solid, medial, opacity, baked);
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

#include "core/control/registry.h"
REGISTER_EFFECT(DreamBalls)
