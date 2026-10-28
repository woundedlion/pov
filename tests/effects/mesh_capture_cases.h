/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Mesh effect capture and profile cases. */

/**
 * @brief Measures capture and scan work of the large scanned meshes under the
 *        device's segment clip: the top half of the rows and alternating
 *        half-width column runs.
 */
inline void measure_mesh_capture_work() {
  constexpr int W = 288, H = 144;
  const auto measure = [](auto &effect, const char *name, int frames,
                          auto draw) {
    Scan::MeshCaptureCounters total, peak;
    size_t pixels = 0, shaded = 0, pixel_peak = 0;
    int measured = 0;
    for (int frame = 0; frame < frames; ++frame) {
      const int first = frame % 2 ? W / 2 : 0;
      effect.set_clip(0, H / 2, first, first + W / 2);
      pin_frame_clock(frame);
      Scan::mesh_capture_counters = {};
      HS_SCAN_METRIC(hs::g_scan_metrics.reset());
      if (!draw(frame))
        break;
      ++measured;
      const auto &work = Scan::mesh_capture_counters;
      const auto add = [](size_t &sum, size_t &top, size_t value) {
        sum += value;
        top = std::max(top, value);
      };
      add(total.vertices, peak.vertices, work.vertices);
      add(total.evaluations, peak.evaluations, work.evaluations);
      add(total.captures, peak.captures, work.captures);
      add(total.shared_faces, peak.shared_faces, work.shared_faces);
      add(total.seam_faces, peak.seam_faces, work.seam_faces);
      add(total.cross_faces, peak.cross_faces, work.cross_faces);
      add(total.polar_faces, peak.polar_faces, work.polar_faces);
      add(total.swept_faces, peak.swept_faces, work.swept_faces);
      add(total.face_copies, peak.face_copies, work.face_copies);
      add(total.split_vertices, peak.split_vertices, work.split_vertices);
      add(total.scans, peak.scans, work.scans);
      add(total.face_setups, peak.face_setups, work.face_setups);
      HS_SCAN_METRIC(add(pixels, pixel_peak, hs::g_scan_metrics.pixels_tested));
      HS_SCAN_METRIC(shaded += hs::g_scan_metrics.shade_candidates);
    }
    std::printf(
        "  [capture work %s] frames=%d, total/peak per frame:\n"
        "    moving vertices %zu/%zu evaluations %zu/%zu captures %zu/%zu\n"
        "    shared faces %zu/%zu seam copies %zu/%zu cross %zu/%zu polar "
        "%zu/%zu (swept %zu/%zu) face copies %zu/%zu split vertices "
        "%zu/%zu\n"
        "    scans %zu/%zu face setups %zu/%zu pixels tested %zu/%zu "
        "shaded %zu\n",
        name, measured, total.vertices, peak.vertices, total.evaluations,
        peak.evaluations, total.captures, peak.captures, total.shared_faces,
        peak.shared_faces, total.seam_faces, peak.seam_faces, total.cross_faces,
        peak.cross_faces, total.polar_faces, peak.polar_faces,
        total.swept_faces, peak.swept_faces, total.face_copies,
        peak.face_copies, total.split_vertices, peak.split_vertices,
        total.scans, peak.scans, total.face_setups, peak.face_setups, pixels,
        pixel_peak, shaded);
  };
  for (const uint8_t node :
       {ConwayGraph::TRUNCATED_ICOSIDODECAHEDRON,
        ConwayGraph::TRUNCATED_ICOSAHEDRON, ConwayGraph::CUBOCTAHEDRON,
        ConwayGraph::TETRAHEDRON}) {
    reset_effect_globals();
    pin_frame_clock(0);
    HankinSolids<W, H> effect;
    effect.init();
    HankinSolidsWhiteBox::load_node(effect, node);
    HankinSolidsWhiteBox::rebuild(effect, 0.4f);
    const math::Vector axis = math::Vector(0.3f, 1.0f, 0.2f).normalized();
    char name[64];
    std::snprintf(name, sizeof(name), "Hankin %s",
                  Solids::simple_registry[node].name);
    measure(effect, name, 24, [&](int frame) {
      const auto &vertices = HankinSolidsWhiteBox::mesh_vertices(effect);
      std::vector<math::Vector> previous(vertices.begin(), vertices.end());
      HankinSolidsWhiteBox::rebuild(effect, 0.4f + 0.04f * (frame + 1));
      FaceCopies copies;
      HankinSolidsWhiteBox::draw_copies(
          effect, previous, math::make_rotation(axis, 0.02f * frame),
          math::make_rotation(axis, 0.02f * (frame + 1)), copies);
      return true;
    });
  }
  for (int shape : {6, 8}) {
    reset_effect_globals();
    pin_frame_clock(0);
    IslamicStars<W, H> effect;
    IslamicBuildProbe::set_trans_speed(effect, 4.0f);
    effect.init();
    IslamicBuildProbe::set_burst_size(effect, shape == 6 ? 8 : 1);
    IslamicBuildProbe::spawn_entry(effect,
                                   *IslamicBuildProbe::presets()[shape].params);
    measure(effect, shape == 6 ? "Islamic shape 6" : "Islamic shape 8", 120,
            [&](int) {
              effect.draw_frame();
              effect.advance_display();
              return IslamicBuildProbe::selected_entry(effect) ==
                     IslamicBuildProbe::presets()[shape].params;
            });
  }
  hs::clear_mock_time();
}

/// Frames of the device profile image's IslamicStars run.
inline constexpr int PROFILE_REPLAY_FRAMES = 280;

/**
 * @brief Replays the device profile image's IslamicStars run on the host and
 *        prints its deformation work per frame and per profile window.
 * @param preset Preset the image pins (HS_PROFILE_PRESET).
 * @param trans_speed The image's Trans Speed (HS_PROFILE_TRANS_SPEED).
 * @param frames Frames to draw, the first the full-width prime.
 * @param report Whether to print the work.
 * @details Mirrors the segmented driver on segment 0: the effect stream is
 * seeded from the effect's identity, the effect is built with its Trans
 * Speed, initialized, pinned and paused at @p preset, primed full width over
 * the segment's rows, then drawn one alternating half of them per frame.
 * Frame and window numbers match the image's `f` lines and windows. Run with
 * HS_MESH_CAPTURE_ONLY=3 (HS_REPLAY_PRESET, HS_REPLAY_TRANS_SPEED and
 * HS_REPLAY_FRAMES override 6, 4 and PROFILE_REPLAY_FRAMES).
 * HS_FACE_AUDIT_PIX names a binary output of native uint32_t records:
 * face, x, y, float distance bits, convex flag; UINT32_MAX, frame, 0, 0, 0
 * terminates each frame. HS_FACE_AUDIT_REFERENCE enables reference face
 * visibility when present. HS_FACE_AUDIT_IMAGES names a binary output for
 * every frame: native uint32_t frame index followed by W * H raw Pixels.
 */
inline void replay_islamicstars_profile(size_t preset, float trans_speed,
                                        int frames, bool report = true) {
  constexpr int W = 288, H = 144, WINDOW = 16;
  reset_effect_globals();
  pin_frame_clock(0);
  hs::seed_effect_stream<IslamicStars<W, H>>();
  IslamicStars<W, H> effect;
  // The image compiles its Trans Speed into the constructor.
  IslamicBuildProbe::set_trans_speed(effect, trans_speed);
  effect.init();
  effect.profile_select_preset(preset);
  HS_EXPECT_TRUE(effect.updateParameter("Trans Speed", trans_speed) ==
                 ParamSetResult::APPLIED);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *audit_path = report ? std::getenv("HS_FACE_AUDIT_PIX") : nullptr;
  Scan::face_audit_reference =
      report && std::getenv("HS_FACE_AUDIT_REFERENCE") != nullptr;
  const char *image_path =
      report ? std::getenv("HS_FACE_AUDIT_IMAGES") : nullptr;
  FILE *images = image_path ? std::fopen(image_path, "wb") : nullptr;
  FILE *audit = audit_path ? std::fopen(audit_path, "wb") : nullptr;
#pragma clang diagnostic pop
  static FILE *audit_sink = nullptr;
  audit_sink = audit;
  if (audit)
    Scan::face_audit_candidate = [](size_t face, int x, int y, float d,
                                    bool convex) {
      const uint32_t record[] = {
          static_cast<uint32_t>(face), static_cast<uint32_t>(x),
          static_cast<uint32_t>(y), std::bit_cast<uint32_t>(d),
          static_cast<uint32_t>(convex)};
      std::fwrite(record, sizeof(record), 1, audit_sink);
    };
  Scan::MeshCaptureCounters window;
  [[maybe_unused]] size_t probes = 0;
  for (int frame = 1; frame <= frames; ++frame) {
    if (frame == 1)
      effect.set_clip(0, H / 2, 0, W);
    else
      effect.set_clip(0, H / 2, frame % 2 ? W / 2 : 0, frame % 2 ? W : W / 2);
    pin_frame_clock(frame);
    Scan::mesh_capture_counters = {};
    std::fill_n(Scan::visibility_counts, 4, size_t(0));
    HS_SCAN_METRIC(hs::g_scan_metrics.reset());
    const size_t passes = IslamicBuildProbe::vertex_passes<W, H>();
    effect.draw_frame();
    effect.advance_display();
    if (audit) {
      const uint32_t record[] = {UINT32_MAX, static_cast<uint32_t>(frame), 0, 0,
                                 0};
      std::fwrite(record, sizeof(record), 1, audit);
    }
    if (images) {
      const uint32_t index = static_cast<uint32_t>(frame);
      std::fwrite(&index, sizeof(index), 1, images);
      std::fwrite(effect.display_buffer(), sizeof(Pixel), W * H, images);
    }
    const auto &w = Scan::mesh_capture_counters;
    HS_SCAN_METRIC(probes += hs::g_scan_metrics.pixels_tested);
    if (report)
      std::printf(
          "  visibility %d faces %zu/%zu vertices %zu/%zu reach %.6f setups %zu\n",
          frame, Scan::visibility_counts[1], Scan::visibility_counts[0],
          Scan::visibility_counts[3], Scan::visibility_counts[2],
          Scan::visibility_audit_reach, w.face_setups);

    if (report && w.deformations)
      std::printf(
          "  f %d deform %zu moving %zu path evaluations %zu (scan %zu) "
          "captures %zu polar %zu swept %zu copies %zu\n",
          frame, w.deformations, w.vertices,
          IslamicBuildProbe::vertex_passes<W, H>() - passes, w.evaluations,
          w.captures, w.polar_faces, w.swept_faces, w.face_copies);
    window.deformations += w.deformations;
    window.vertices += w.vertices;
    window.evaluations += w.evaluations;
    window.polar_faces += w.polar_faces;
    window.swept_faces += w.swept_faces;
    window.face_copies += w.face_copies;
    if (report && frame % WINDOW == 0)
      std::printf("  window %d-%d: deform %zu polar %zu swept %zu copies %zu "
                  "moving %zu evaluations %zu probes %zu\n",
                  frame - WINDOW + 1, frame, window.deformations,
                  window.polar_faces, window.swept_faces, window.face_copies,
                  window.vertices, window.evaluations, probes);
    if (frame % WINDOW == 0) {
      window = {};
      probes = 0;
    }
  }
  if (report)
    std::printf("  scratch peak a %zu of %zu, b %zu of %zu\n",
                scratch_arena_a.get_lifetime_high_water_mark(),
                scratch_arena_a.get_capacity(),
                scratch_arena_b.get_lifetime_high_water_mark(),
                scratch_arena_b.get_capacity());
  Scan::face_audit_candidate = nullptr;
  Scan::face_audit_reference = false;
  if (images)
    std::fclose(images);
  if (audit)
    std::fclose(audit);
  hs::clear_mock_time();
}

/**
 * @brief Checks that one preset's profile run deforms within the effect's
 *        scratch arenas, which trap on overflow.
 * @param preset Preset the run pins.
 * @param frames Frames the run draws.
 */
inline void check_islamicstars_motion_fits_scratch(size_t preset, int frames) {
  HS_CONTEXT("preset", static_cast<int>(preset));
  scratch_arena_a.reset_peak_tracking();
  scratch_arena_b.reset_peak_tracking();
  replay_islamicstars_profile(preset, 4.0f, frames, false);
  std::printf("  [islamic scratch, preset %zu] peak a %zu, b %zu\n", preset,
              scratch_arena_a.get_lifetime_high_water_mark(),
              scratch_arena_b.get_lifetime_high_water_mark());
}

/**
 * @brief Preset 6, the scratch arenas' tightest, fits through its next
 *        shape's spawn at frame 151, where scratch A peaks.
 */
inline void test_islamicstars_peak_fits_scratch() {
  check_islamicstars_motion_fits_scratch(6, 160);
}

/**
 * @brief Every preset's whole profile run, its builds and its ripples,
 *        deforms within the effect's scratch arenas.
 */
inline void test_islamicstars_motion_fits_scratch() {
  for (size_t preset = 0; preset < std::size(IslamicBuildProbe::presets());
       ++preset)
    check_islamicstars_motion_fits_scratch(preset, PROFILE_REPLAY_FRAMES);
}

/** @brief Runs the focused IslamicStars white-box cases. */
inline void run_islamicstars_cases() {
  HS_RUN(run_case, test_islamicstars_seed_sprite_fade_in);
  HS_RUN(run_case, test_islamicstars_opacity_shades_at_column_time);
  HS_RUN(run_case, test_islamicstars_ripple_motion_is_seam_free);
  HS_RUN(run_case, test_islamicstars_seam_copies_are_one_polygon);
  HS_RUN(run_case, test_islamicstars_playback_motion_is_seam_free);
  HS_RUN(run_case, test_islamicstars_build_legs_hand_their_path);
  HS_RUN(run_case, test_islamicstars_opacity_history);
  HS_RUN(run_case, test_islamicstars_burst_size_is_snapshotted_per_spawn);
  HS_RUN(run_case, test_islamicstars_smooth_recipe_completion);
  HS_RUN(run_case, test_islamicstars_default_cycle_frames);
  HS_RUN(run_case, test_islamicstars_presets_cover_roster);
  HS_RUN(run_case, test_islamicstars_manual_select_cuts_build);
  HS_RUN(run_case, test_islamicstars_pause_holds_shape);

  HS_RUN(run_case, test_islamicstars_prepared_build);
  HS_RUN(run_case, test_islamicstars_prepared_ripples);
  HS_RUN(run_case, test_islamicstars_endpoint_groups);
  HS_RUN(run_case, test_islamicstars_peak_fits_scratch);
}
