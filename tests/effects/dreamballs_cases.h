/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// DreamBalls geometry and lifecycle coverage.

/** @brief White-box accessor for DreamBalls' private preset-cycle bookkeeping. */
struct DreamBallsWhiteBox {
  using DB = DreamBalls<SMALL_W, SMALL_H>;
  static constexpr int PRESETS = static_cast<int>(DB::PRESETS.size());
  static constexpr size_t SOLID_COUNT = DB::SOLID_COUNT;
  static constexpr size_t MAX_SOLID_EDGES = Solids::MAX_SOLID_EDGES;
  static constexpr size_t SCRATCH_A_PEAK_BYTES = DB::SCRATCH_A_PEAK_BYTES;
  static constexpr float WEAVE_GAP_MIN = DB::WEAVE_GAP_MIN;

  static int active_bake(const DB &db) { return db.active_bake; }
  // Advance the choreography, then re-spawn.
  static void advance(DB &db) {
    HS_EXPECT_TRUE(db.advance_preset());
    db.spawn_sprite();
  }
  // Re-spawn of the current preset (no advance).
  static void respawn(DB &db) { db.spawn_sprite(); }
  static DB::BaseMesh preset_mesh(int idx) {
    return DB::PRESETS[idx].params.base_mesh;
  }
  static const DB::Params &preset_params(int idx) {
    return DB::PRESETS[idx].params;
  }
  static const Palette *preset_palette(const DB &db, int idx) {
    return db.preset_palettes[idx];
  }
  static const Palette *blood_stream_falloff(const DB &db) {
    return &db.blood_stream_falloff;
  }
  static bool solid_loaded(const DB &db, size_t idx) {
    return !db.loaded_solids[idx].mesh_state.vertices.is_empty();
  }
  static bool four_regular(const DB &db, size_t idx) {
    return db.loaded_solids[idx].four_regular;
  }
  static const auto &original_edges(const DB &db, size_t idx) {
    return db.loaded_solids[idx].original_edges;
  }
  static const auto &automatic_edges(const DB &db, size_t idx) {
    return db.loaded_solids[idx].automatic_edges;
  }
  static const auto &medial_edges(const DB &db, size_t idx) {
    const auto &solid = db.loaded_solids[idx];
    return solid.four_regular ? solid.medial_edges : solid.automatic_edges;
  }
  static std::vector<uint16_t> woven_start_owners(const DB &db, size_t idx,
                                                  bool medial) {
    const auto &solid = db.loaded_solids[idx];
    const auto &edges = medial && solid.four_regular ? solid.medial_edges
                                                     : solid.automatic_edges;
    const size_t vertex_count =
        medial ? solid.original_edges.size() : solid.mesh_state.vertices.size();
    std::vector<uint16_t> owners(vertex_count);
    DB::assign_woven_start_owners(edges, owners.data(), owners.size());
    return owners;
  }
  static bool owns_woven_start(const ArenaVector<Plot::Mesh::Edge> &edges,
                               const std::vector<uint16_t> &owners,
                               size_t edge_index) {
    return DB::owns_woven_start_sample(edges, owners.data(), edge_index);
  }
  static void
  assign_woven_start_owners(const ArenaVector<Plot::Mesh::Edge> &edges,
                            uint16_t *owners, size_t vertex_count) {
    DB::assign_woven_start_owners(edges, owners, vertex_count);
  }
  static size_t source_vertex_count(const DB &db, size_t idx) {
    return db.loaded_solids[idx].mesh_state.vertices.size();
  }
  static float under_gap_alpha(float edge_t, float gap) {
    return DB::under_gap_alpha(edge_t, gap);
  }
  static math::Vector woven_vertex(const DB &db, size_t solid, bool medial,
                                   size_t vertex) {
    return DB::woven_vertex(db.loaded_solids[solid], medial, vertex);
  }
  static DB::BaseMesh live_mesh(const DB &db) { return db.params.base_mesh; }
  static DB::WeaveTopology live_weave_topology(const DB &db) {
    return db.params.weave_topology;
  }
  static float live_weave_gap(const DB &db) { return db.params.weave_gap; }
  static float &num_copies(DB &db) { return db.params.num_copies; }
};

/**
 * @brief Drives spawn_sprite across a full preset cycle and asserts the bake-slot
 *        ping-pong, the modulo preset advance, and the reseed-on-change guard.
 * @details Each spawn flips the bake slot so a fading-out sprite keeps its own
 *          LUT. A re-spawn of the same preset holds params so a live slider
 *          edit survives.
 */
inline void test_dreamballs_preset_cycle_bookkeeping() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init(); // runs spawn_sprite() at preset 0

  // init() spawned preset 0 and flipped the bake slot once (0 -> 1).
  HS_EXPECT_EQ(db.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::active_bake(db), 1);
  HS_EXPECT_EQ(WB::live_mesh(db), WB::preset_mesh(0));

  // Pins each row's structural selections; magnitudes are range-checked
  // during the cycle.
  struct Row {
    WB::DB::BaseMesh base_mesh;
    WB::DB::WeaveTopology weave_topology;
    const Palette *palette;
  };
  using BaseMesh = WB::DB::BaseMesh;
  constexpr auto AUTOMATIC = WB::DB::WeaveTopology::AUTOMATIC;
  const Row rows[] = {
      {BaseMesh::RHOMBICUBOCTAHEDRON, AUTOMATIC, WB::blood_stream_falloff(db)},
      {BaseMesh::RHOMBICOSIDODECAHEDRON, AUTOMATIC,
       WB::blood_stream_falloff(db)},
      {BaseMesh::TRUNCATED_CUBOCTAHEDRON, AUTOMATIC, &Palettes::RICH_SUNSET},
      {BaseMesh::ICOSIDODECAHEDRON, AUTOMATIC, &Palettes::LAVENDER_LAKE},
      {BaseMesh::SNUB_CUBE, AUTOMATIC, &Palettes::MAUVE_FADE},
      {BaseMesh::TRUNCATED_DODECAHEDRON, AUTOMATIC, &Palettes::CORAL_BLUE},
      {BaseMesh::TRIAKIS_ICOSAHEDRON, AUTOMATIC, &Palettes::BRUISED_MOSS},
      {BaseMesh::TRIAKIS_ICOSAHEDRON, AUTOMATIC, &Palettes::LAVENDER_LAKE},
      {BaseMesh::DISDYAKIS_TRIACONTAHEDRON, AUTOMATIC, &Palettes::PLUM_SUNRISE},
      {BaseMesh::TRIAKIS_ICOSAHEDRON, AUTOMATIC, &Palettes::BRUISED_MANGO},
  };
  HS_EXPECT_EQ(std::size(rows), static_cast<size_t>(WB::PRESETS));
  for (size_t i = 0;
       i < std::min(std::size(rows), static_cast<size_t>(WB::PRESETS)); ++i) {
    HS_CONTEXT("preset", i);
    const auto &params = WB::preset_params(i);
    HS_EXPECT_EQ(params.base_mesh, rows[i].base_mesh);
    HS_EXPECT_EQ(params.weave_topology, rows[i].weave_topology);
    HS_EXPECT_TRUE(WB::preset_palette(db, i) == rows[i].palette);
  }

  // Preset rows are assigned straight into params, bypassing register_param's
  // range check.
  auto expect_in_range = [&]() {
    for (const auto &def : db.getParameters()) {
      HS_CONTEXT(def.name);
      const float v = def.get();
      HS_EXPECT_TRUE(std::isfinite(v));
      HS_EXPECT_GE(v, def.min);
      HS_EXPECT_LE(v, def.max);
      if (def.option_count > 0) {
        HS_EXPECT_EQ(v, std::floor(v));
        HS_EXPECT_LT(v, static_cast<float>(def.option_count));
      }
    }
  };

  // Two full cycles: the bake slot ping-pongs and params reseed each step.
  int expect_bake = WB::active_bake(db); // 1
  std::vector<std::vector<float>> live_rows;
  for (int step = 1; step <= 2 * WB::PRESETS; ++step) {
    WB::advance(db);
    expect_bake ^= 1;
    const int safe = step % WB::PRESETS;
    HS_CONTEXT("preset", safe);
    HS_EXPECT_EQ(WB::active_bake(db), expect_bake);
    HS_EXPECT_EQ(db.getPresetIndex(), static_cast<size_t>(safe));
    HS_EXPECT_EQ(WB::live_mesh(db), WB::preset_mesh(safe));
    expect_in_range();
    if (step <= WB::PRESETS) {
      std::vector<float> live_row;
      for (const auto &def : db.getParameters())
        live_row.push_back(def.get());
      live_rows.push_back(live_row);
    }
  }

  // Identical parameter vectors would be one preset visited twice.
  for (size_t i = 0; i < live_rows.size(); ++i)
    for (size_t j = i + 1; j < live_rows.size(); ++j) {
      HS_CONTEXT("preset pair", static_cast<int>(i), static_cast<int>(j));
      HS_EXPECT(live_rows[i] != live_rows[j],
                "each DreamBalls preset must be distinct");
    }

  // A same-preset re-spawn keeps a live slider edit but still flips the bake
  // slot.
  const float sentinel = WB::num_copies(db) + 5.0f;
  WB::num_copies(db) = sentinel;
  const size_t held_idx = db.getPresetIndex();
  expect_bake ^= 1;
  WB::respawn(db);
  HS_EXPECT_EQ(db.getPresetIndex(), held_idx);    // preset unchanged
  HS_EXPECT_EQ(WB::num_copies(db), sentinel);     // live edit preserved
  HS_EXPECT_EQ(WB::active_bake(db), expect_bake); // bake slot still flipped
}

/** @brief Verifies the Base Mesh dropdown covers every simple solid family. */
inline void test_dreamballs_base_mesh_selector() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  const auto *base_mesh = db.getParameters().find("Base Mesh");
  HS_EXPECT_TRUE(base_mesh != nullptr);
  if (!base_mesh)
    return;

  constexpr int EXPECTED_SOLIDS =
      static_cast<int>(Solids::PLATONIC_COUNT + Solids::ARCHIMEDEAN_COUNT +
                       Solids::CATALAN_COUNT);
  HS_EXPECT_EQ(base_mesh->option_count, EXPECTED_SOLIDS);
  HS_EXPECT_TRUE(base_mesh->animated);
  for (int i = 0; i < EXPECTED_SOLIDS; ++i)
    HS_EXPECT_TRUE(WB::solid_loaded(db, static_cast<size_t>(i)));
  HS_EXPECT_EQ(std::string_view(base_mesh->options[0]), "Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[17]), "Snub Dodecahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[18]), "Triakis Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->options[EXPECTED_SOLIDS - 1]),
               "Pentagonal Hexecontahedron");
  HS_EXPECT_EQ(std::string_view(base_mesh->export_options[EXPECTED_SOLIDS - 1]),
               "BaseMesh::PENTAGONAL_HEXECONTAHEDRON");

  HS_EXPECT_EQ(
      db.updateParameter("Base Mesh", static_cast<float>(EXPECTED_SOLIDS - 1)),
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_mesh(db), WB::DB::BaseMesh::PENTAGONAL_HEXECONTAHEDRON);
  HS_EXPECT_TRUE(db.animations_paused());

  db.setAnimationsPaused(false);
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }
  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(db);
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief Renders the widest solid the dropdown reaches under medial topology and
 *        checks SCRATCH_A_PEAK_BYTES against the real frame peak.
 * @details Forced medial topology stages one vertex per source edge and
 * one framed edge per medial edge at the static memory bound.
 */
inline void test_dreamballs_max_edge_solid_render() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  size_t widest = 0;
  size_t widest_edges = 0;
  for (size_t i = 0; i < WB::SOLID_COUNT; ++i) {
    const size_t edges = WB::original_edges(db, i).size();
    if (edges > widest_edges) {
      widest_edges = edges;
      widest = i;
    }
  }
  HS_EXPECT_EQ(widest_edges, WB::MAX_SOLID_EDGES);

  HS_EXPECT_EQ(db.updateParameter("Base Mesh", static_cast<float>(widest)),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(db.updateParameter("Weave Topology", 2.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db), WB::DB::WeaveTopology::MEDIAL);
  HS_EXPECT_EQ(WB::medial_edges(db, widest).size(), 2 * widest_edges);

  db.setAnimationsPaused(false);
  scratch_arena_a.reset_high_water_mark();
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }

  // Three per-vertex buffers plus the framed vertex + edge-head mesh, at the
  // medial bound.
  const size_t staged_bytes = 6 * widest_edges * sizeof(math::Vector);
  HS_EXPECT_GE(scratch_arena_a.get_high_water_mark(), staged_bytes);
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), WB::SCRATCH_A_PEAK_BYTES);

  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(db);
  HS_EXPECT_GT(energy, 0u);
}

/** @brief Verifies automatic, forced-medial, and defect weave topology. */
inline void test_dreamballs_weave_topology() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();

  WB::DB db;
  db.init();

  const auto *topology = db.getParameters().find("Weave Topology");
  const auto *gap = db.getParameters().find("Weave Gap");
  HS_EXPECT_TRUE(topology != nullptr);
  HS_EXPECT_TRUE(gap != nullptr);
  if (!topology || !gap)
    return;

  HS_EXPECT_EQ(topology->option_count, 3);
  HS_EXPECT_EQ(std::string_view(topology->options[0]), "Automatic");
  HS_EXPECT_EQ(std::string_view(topology->options[1]), "Original with defects");
  HS_EXPECT_EQ(std::string_view(topology->options[2]), "Medial");
  HS_EXPECT_EQ(std::string_view(topology->export_options[1]),
               "WeaveTopology::ORIGINAL_WITH_DEFECTS");
  HS_EXPECT_EQ(std::string_view(topology->export_options[2]),
               "WeaveTopology::MEDIAL");
  HS_EXPECT_TRUE(topology->animated);
  HS_EXPECT_TRUE(gap->animated);

  constexpr int EXPECTED_SOLIDS =
      static_cast<int>(Solids::PLATONIC_COUNT + Solids::ARCHIMEDEAN_COUNT +
                       Solids::CATALAN_COUNT);
  int four_regular_count = 0;
  for (int i = 0; i < EXPECTED_SOLIDS; ++i) {
    const auto &original = WB::original_edges(db, i);
    const auto &automatic = WB::automatic_edges(db, i);
    const auto &medial = WB::medial_edges(db, i);
    const bool four_regular = WB::four_regular(db, i);
    const bool expected_four_regular =
        i == static_cast<int>(WB::DB::BaseMesh::OCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::CUBOCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::RHOMBICUBOCTAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::ICOSIDODECAHEDRON) ||
        i == static_cast<int>(WB::DB::BaseMesh::RHOMBICOSIDODECAHEDRON);
    HS_EXPECT_EQ(four_regular, expected_four_regular);
    four_regular_count += four_regular ? 1 : 0;

    const size_t automatic_vertex_count =
        four_regular ? WB::source_vertex_count(db, i) : original.size();
    const auto automatic_start_owners =
        WB::woven_start_owners(db, i, !four_regular);
    HS_EXPECT_EQ(automatic.size(),
                 four_regular ? original.size() : 2 * original.size());
    HS_EXPECT_EQ(medial.size(), 2 * original.size());

    std::vector<int> incoming(automatic_vertex_count, 0);
    std::vector<int> outgoing(automatic_vertex_count, 0);
    std::vector<int> owned_starts(automatic_vertex_count, 0);
    for (size_t edge_index = 0; edge_index < automatic.size(); ++edge_index) {
      const auto &edge = automatic[edge_index];
      HS_EXPECT_LT(edge.u, automatic_vertex_count);
      HS_EXPECT_LT(edge.v, automatic_vertex_count);
      if (edge.u < automatic_vertex_count && edge.v < automatic_vertex_count) {
        outgoing[edge.u]++;
        incoming[edge.v]++;
        owned_starts[edge.u] +=
            WB::owns_woven_start(automatic, automatic_start_owners, edge_index)
                ? 1
                : 0;

        const math::Vector from =
            WB::woven_vertex(db, i, !four_regular, edge.u);
        const math::Vector to = WB::woven_vertex(db, i, !four_regular, edge.v);
        const math::Vector frame_u = math::tangent_axis(from);
        const math::Vector offset =
            frame_u * 0.6f + math::cross(from, frame_u) * 0.8f;
        const math::Vector transported =
            math::parallel_transport(from, to, offset);
        HS_EXPECT_NEAR(math::dot(transported, to), 0.0f, 2e-5f);
        HS_EXPECT_NEAR(math::dot(transported, transported), 1.0f, 2e-5f);
        HS_EXPECT_VEC(math::parallel_transport(to, from, transported), offset,
                      2e-5f);
      }
    }
    for (size_t vertex = 0; vertex < automatic_vertex_count; ++vertex) {
      HS_EXPECT_EQ(incoming[vertex], 2);
      HS_EXPECT_EQ(outgoing[vertex], 2);
      HS_EXPECT_EQ(owned_starts[vertex], 1);
      HS_EXPECT_LT(automatic_start_owners[vertex], automatic.size());
      if (automatic_start_owners[vertex] < automatic.size())
        HS_EXPECT_EQ(automatic[automatic_start_owners[vertex]].u, vertex);
    }

    const auto medial_start_owners = WB::woven_start_owners(db, i, true);
    std::fill(incoming.begin(), incoming.end(), 0);
    std::fill(outgoing.begin(), outgoing.end(), 0);
    incoming.resize(original.size(), 0);
    outgoing.resize(original.size(), 0);
    owned_starts.assign(original.size(), 0);
    for (size_t edge_index = 0; edge_index < medial.size(); ++edge_index) {
      const auto &edge = medial[edge_index];
      HS_EXPECT_LT(edge.u, original.size());
      HS_EXPECT_LT(edge.v, original.size());
      if (edge.u < original.size() && edge.v < original.size()) {
        outgoing[edge.u]++;
        incoming[edge.v]++;
        owned_starts[edge.u] +=
            WB::owns_woven_start(medial, medial_start_owners, edge_index) ? 1
                                                                          : 0;
      }
    }
    for (size_t vertex = 0; vertex < original.size(); ++vertex) {
      HS_EXPECT_EQ(incoming[vertex], 2);
      HS_EXPECT_EQ(outgoing[vertex], 2);
      HS_EXPECT_EQ(owned_starts[vertex], 1);
      HS_EXPECT_LT(medial_start_owners[vertex], medial.size());
      if (medial_start_owners[vertex] < medial.size())
        HS_EXPECT_EQ(medial[medial_start_owners[vertex]].u, vertex);
    }
  }
  HS_EXPECT_EQ(four_regular_count, 5);

  HS_EXPECT_NEAR(WB::under_gap_alpha(0.0f, 0.2f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(0.8f, 0.2f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(0.9f, 0.2f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(WB::under_gap_alpha(1.0f, 0.2f), 0.0f, 1e-6f);

  HS_EXPECT_EQ(db.updateParameter("Base Mesh", 0.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(db.updateParameter("Weave Topology", 1.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(db.updateParameter("Weave Gap", 0.25f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db),
               WB::DB::WeaveTopology::ORIGINAL_WITH_DEFECTS);
  HS_EXPECT_NEAR(WB::live_weave_gap(db), 0.25f, 1e-6f);

  HS_EXPECT_EQ(
      db.updateParameter("Base Mesh",
                         static_cast<float>(WB::DB::BaseMesh::OCTAHEDRON)),
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(db.updateParameter("Weave Topology", 2.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::live_weave_topology(db), WB::DB::WeaveTopology::MEDIAL);

  db.setAnimationsPaused(false);
  for (int frame = 0; frame < 20; ++frame) {
    db.draw_frame();
    db.advance_display();
  }
  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(db);
  HS_EXPECT_GT(energy, 0u);
}

/** @brief Renders the defect weave topology and its gap fade. */
inline void test_dreamballs_defect_weave_renders() {
  using WB = DreamBallsWhiteBox;
  auto render_defects = [](float weave_gap) {
    reset_effect_globals();
    WB::DB defects;
    defects.init();
    HS_EXPECT_EQ(defects.updateParameter("Base Mesh", 0.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(defects.updateParameter("Weave Topology", 1.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(defects.updateParameter("Weave Gap", weave_gap),
                 ParamSetResult::APPLIED);
    defects.setAnimationsPaused(false);
    for (int frame = 0; frame < 20; ++frame) {
      defects.draw_frame();
      defects.advance_display();
    }
    HS_EXPECT_GT((frame_energy<SMALL_W, SMALL_H>(defects)), 0u);
    std::vector<Pixel> frame;
    capture_frame<SMALL_W, SMALL_H>(defects, frame);
    return frame;
  };
  HS_EXPECT(render_defects(WB::WEAVE_GAP_MIN) != render_defects(0.25f),
            "the defect weave gap must reach the rendered frame");
}

/**
 * @brief Verifies pause freezes DreamBalls' sprite and next-preset clock.
 */
inline void test_dreamballs_respawn_fires_and_honors_pause() {
  using WB = DreamBallsWhiteBox;
  reset_effect_globals();
  WB::DB db;
  db.init();
  const int held_bake = WB::active_bake(db);
  db.setAnimationsPaused(true);
  for (int f = 0; f < 340; ++f) {
    db.draw_frame();
    db.advance_display();
  }
  HS_EXPECT_EQ(db.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::active_bake(db), held_bake);

  db.setAnimationsPaused(false);
  for (int f = 0; f < 340; ++f) {
    db.draw_frame();
    db.advance_display();
  }
  HS_EXPECT_EQ(db.getPresetIndex(), 1u);
  HS_EXPECT_EQ(WB::active_bake(db), held_bake ^ 1);
}
