/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// MeshFeedback: draw ordering, base-mesh selection, export arity and storage
// ---------------------------------------------------------------------------

/** @brief White-box accessor for MeshFeedback's style / noise / preset state. */
struct MeshFeedbackWhiteBox {
  using MF = MeshFeedback<SMALL_W, SMALL_H>;

  static const Feedback::Style &style(const MF &fx) { return fx.params.style; }
  static MF::BaseMesh base_mesh(const MF &fx) { return fx.params.base_mesh; }
  static MF::BaseMesh active_base_mesh(const MF &fx) {
    return fx.active_base_mesh;
  }
  static MF::BaseMesh preset_base_mesh(MF &fx, size_t index) {
    return fx.preset_params(index).base_mesh;
  }
  static size_t mesh_vertices(const MF &fx) { return fx.mesh.vertices.size(); }
  static size_t mesh_edges(const MF &fx) { return fx.edges.size(); }
  static const Animation::NoiseParams &noise(const MF &fx) {
    return fx.noise_params;
  }
  static size_t preset_index(const MF &fx) { return fx.getPresetIndex(); }
  static size_t mesh_storage_mark(const MF &fx) { return fx.mesh_storage_mark; }
  static void rebuild_mesh(MF &fx, MF::BaseMesh base_mesh) {
    fx.rebuild_mesh(base_mesh);
  }
  /** @brief True when @p target lies inside the effect's Params storage. */
  static bool in_params(const MF &fx, const void *target) {
    const auto *begin = reinterpret_cast<const unsigned char *>(&fx.params);
    const auto *at = static_cast<const unsigned char *>(target);
    return at >= begin && at < begin + sizeof(fx.params);
  }
};

/**
 * @brief Renders MeshFeedback for `frames` frames and copies out the last frame.
 * @param out Receives the displayed frame, SMALL_W * SMALL_H pixels row-major.
 * @param frames Number of frames to render.
 * @param feedback Value written to the Feedback toggle before the first frame.
 */
inline void meshfeedback_capture(std::vector<Pixel> &out, int frames,
                                 bool feedback) {
  reset_effect_globals();
  hs::set_mock_time(0, 0);

  MeshFeedbackWhiteBox::MF fx;
  fx.init();
  HS_EXPECT_EQ(fx.updateParameter("Feedback", feedback ? 1.0f : 0.0f),
               ParamSetResult::APPLIED);
  for (int f = 0; f < frames; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();
  }
  hs::clear_mock_time();

  capture_frame<SMALL_W, SMALL_H>(fx, out);
}

/**
 * @brief Verifies the feedback flush never decays the same frame's wireframe.
 * @details The flush overwrites the draw buffer with the warped previous
 *          frame. With the flush leading, every channel of a feedback-on render
 *          is at least the feedback-off render's, and strictly greater somewhere.
 */
inline void test_meshfeedback_flush_precedes_mesh_draw() {
  // Well past the empty-history early-out, short of the preset rotation.
  constexpr int FRAMES = 48;
  std::vector<Pixel> lit, bare;
  meshfeedback_capture(lit, FRAMES, true);
  meshfeedback_capture(bare, FRAMES, false);

  int decayed = 0, accumulated = 0;
  const uint64_t bare_energy = frame_energy(bare);
  for (size_t i = 0; i < lit.size(); ++i) {
    const Pixel &a = lit[i], &b = bare[i];
    if (a.r < b.r || a.g < b.g || a.b < b.b)
      ++decayed;
    if (a.r > b.r || a.g > b.g || a.b > b.b)
      ++accumulated;
  }

  std::printf("  MeshFeedback flush order: decayed=%d accumulated=%d "
              "wireframe energy=%llu\n",
              decayed, accumulated,
              static_cast<unsigned long long>(bare_energy));
  // A pair of black frames would agree trivially.
  HS_EXPECT_GT(bare_energy, 0u);
  HS_EXPECT_EQ(decayed, 0);
  HS_EXPECT_GT(accumulated, 0);
}

/**
 * @brief Drives the preset rotation and pins the switch-frame noise sync.
 * @details Crosses two dwell boundaries. The bound NoiseParams must already
 *          carry the incoming preset's scalars when the switch frame ends.
 */
inline void test_meshfeedback_preset_rotation_syncs_noise() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();
  hs::set_mock_time(0, 0);

  MF fx;
  fx.init();
  HS_EXPECT_EQ(WB::preset_index(fx), 0u);

  const auto in_sync = [&fx]() {
    const Feedback::Style &s = WB::style(fx);
    const Animation::NoiseParams &n = WB::noise(fx);
    return n.amplitude == s.amplitude && n.frequency == s.frequency &&
           n.speed == s.speed && n.scale == s.scale;
  };
  const auto matches_preset = [&fx](size_t idx) {
    const Feedback::Style &s = WB::style(fx);
    const Feedback::Style &p = MF::PRESETS[idx].params.style;
    return s.fade == p.fade && s.hue_shift == p.hue_shift &&
           s.scale == p.scale && s.amplitude == p.amplitude &&
           s.frequency == p.frequency && s.speed == p.speed &&
           s.space_fn == p.space_fn && s.color_fn == p.color_fn;
  };

  int switches = 0, desynced = 0, wrong_preset = 0;
  for (int f = 1; f <= 2 * MF::PRESET_DWELL_FRAMES + 1; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    fx.draw_frame();
    fx.advance_display();

    const size_t expect_idx =
        static_cast<size_t>(f / MF::PRESET_DWELL_FRAMES) % MF::PRESETS.size();
    if (WB::preset_index(fx) != expect_idx)
      ++wrong_preset;
    if (f % MF::PRESET_DWELL_FRAMES == 0) {
      ++switches;
      // The live style is the incoming preset, and the noise the flush reads
      // already followed it within this same frame.
      HS_EXPECT_TRUE(matches_preset(expect_idx));
      HS_EXPECT_TRUE(in_sync());
    }
    if (!in_sync())
      ++desynced;
  }
  hs::clear_mock_time();

  std::printf("  MeshFeedback presets: %d switches, %d desynced frames, "
              "%d wrong index\n",
              switches, desynced, wrong_preset);
  HS_EXPECT_EQ(switches, 2);
  HS_EXPECT_EQ(wrong_preset, 0);
  HS_EXPECT_EQ(desynced, 0);
}

/** @brief Verifies MeshFeedback presets and GUI carry the selected solid. */
inline void test_meshfeedback_base_mesh_selector() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  const auto *base_mesh = effect.getParameters().find("Base Mesh");
  HS_EXPECT_TRUE(base_mesh != nullptr);
  if (!base_mesh)
    return;

  HS_EXPECT_EQ(base_mesh->option_count,
               static_cast<int>(Solids::BASE_MESH_COUNT));
  HS_EXPECT_TRUE(base_mesh->animated);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 0), MF::BaseMesh::ICOSAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 1), MF::BaseMesh::DODECAHEDRON);
  HS_EXPECT_EQ(WB::preset_base_mesh(effect, 11),
               MF::BaseMesh::PENTAGONAL_HEXECONTAHEDRON);

  const auto selected = MF::BaseMesh::PENTAGONAL_HEXECONTAHEDRON;
  HS_EXPECT_EQ(
      effect.updateParameter("Base Mesh", static_cast<float>(selected)),
      ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(WB::base_mesh(effect), selected);
  HS_EXPECT_EQ(WB::active_base_mesh(effect), selected);
  HS_EXPECT_GT(WB::mesh_vertices(effect), 0u);
  HS_EXPECT_GT(WB::mesh_edges(effect), 0u);
}

/**
 * @brief Verifies MeshFeedback's preset export covers exactly its Params.
 * @details A preset export emits one value per preset-flagged param into a
 *   PresetEntry<Params> aggregate, so a param backed by a member outside Params
 *   must be marked global or the pasted row will not compile.
 */
inline void test_meshfeedback_preset_export_arity() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  int preset_params = 0;
  for (const auto &def : effect.getParameters()) {
    HS_CONTEXT(def.name);
    HS_EXPECT_EQ(def.preset, WB::in_params(effect, def.target));
    preset_params += def.preset ? 1 : 0;
  }
  HS_EXPECT_GT(preset_params, 0);
}

/** @brief Verifies repeated MeshFeedback mesh changes reuse arena storage. */
inline void test_meshfeedback_mesh_rebuild_reuses_storage() {
  using WB = MeshFeedbackWhiteBox;
  using MF = WB::MF;
  reset_effect_globals();

  MF effect;
  effect.init();

  const size_t storage_mark = WB::mesh_storage_mark(effect);
  std::array<size_t, Solids::BASE_MESH_COUNT> first_cycle_offsets{};
  for (size_t i = 0; i < Solids::BASE_MESH_COUNT; ++i) {
    WB::rebuild_mesh(effect, static_cast<MF::BaseMesh>(i));
    first_cycle_offsets[i] = persistent_arena.get_offset();
    HS_EXPECT_GT(first_cycle_offsets[i], storage_mark);
  }

  for (size_t i = 0; i < Solids::BASE_MESH_COUNT; ++i) {
    WB::rebuild_mesh(effect, static_cast<MF::BaseMesh>(i));
    HS_EXPECT_EQ(persistent_arena.get_offset(), first_cycle_offsets[i]);
  }
}
