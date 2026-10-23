/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// HankinSolids white-box checks
// ---------------------------------------------------------------------------

/** @brief Test access to HankinSolids' morph-chain state and draw path. */
struct HankinSolidsWhiteBox {
  using EffectT = HankinSolids<SMALL_W, SMALL_H>;

  static bool morph_active(const EffectT &effect) {
    return effect.pending_landing != nullptr;
  }
  static uint8_t node(const EffectT &effect) { return effect.node; }
  static void draw_without_topology(EffectT &effect) {
    effect.hankin_mesh.topology.clear();
    effect.hankin_mesh.topology_key = 0;
    const BakedPalette *palettes[EffectT::NUM_PALETTES]{};
    Canvas canvas(effect);
    effect.draw_mesh(canvas, effect.hankin_mesh, palettes, palettes);
  }

  using FrameHandoff = HankinSolidsDetail::FrameHandoff;
  using ShapeWeights = HankinSolidsDetail::ShapeWeights;
  static constexpr int PALETTES = EffectT::NUM_PALETTES;

  static int cycle_frame(const EffectT &effect) {
    return effect.hankin_cycle_frame;
  }
  static int drawn_cycle_frame(const EffectT &effect) {
    return effect.drawn_cycle_frame;
  }
  static int record_cycle_draw(EffectT &effect, int frame) {
    return effect.record_cycle_draw(frame);
  }
  static FrameHandoff handoff(int frame, int previous) {
    return EffectT::frame_handoff(frame, previous);
  }
  static void set_angle(EffectT &effect, float angle) {
    MeshOps::update_hankin(effect.compiled_hankin, effect.hankin_mesh,
                           persistent_arena, angle);
  }
  /** @brief Arms every slot's strap crossfade from the next bank palette. */
  static void arm_strap_crossfade(EffectT &effect) {
    for (int s = 0; s < PALETTES; ++s)
      effect.strap_from[s] = (effect.palette_idx[s] + 1) % PALETTES;
    effect.strap_blend_mask = static_cast<uint8_t>((1u << PALETTES) - 1);
  }
  /**
   * @brief Draws the hankin mesh at the current camera.
   * @param effect Effect whose mesh and palettes are drawn.
   * @param handoff Handoff weights at both endpoints.
   * @param strap_weight Crossfade weight the strap LUTs are pre-blended at.
   * @return Displayed frame.
   */
  static std::vector<Pixel> draw(EffectT &effect, const FrameHandoff &handoff,
                                 float strap_weight) {
    ScratchScope blends(scratch_arena_b);
    const auto &bank = effect.palette_bank.bank.entries;
    BakedPalette blended[PALETTES];
    const BakedPalette *star[PALETTES];
    const BakedPalette *strap[PALETTES];
    for (int s = 0; s < PALETTES; ++s) {
      star[s] = &bank[effect.palette_idx[s]].view();
      blended[s] =
          bake_palette_blend(scratch_arena_b, bank[effect.strap_from[s]],
                             bank[effect.palette_idx[s]], strap_weight);
      strap[s] = &blended[s];
    }
    {
      Canvas canvas(effect);
      effect.draw_mesh(canvas, effect.hankin_mesh, star, strap, handoff);
    }
    effect.advance_display();
    std::vector<Pixel> out;
    capture_frame<SMALL_W, SMALL_H>(effect, out);
    return out;
  }

  /** @brief Returns the current Hankin mesh's model-space vertices. */
  template <typename E>
  static const ArenaVector<math::Vector> &mesh_vertices(const E &effect) {
    return effect.hankin_mesh.vertices;
  }
  static size_t faces(const EffectT &effect) {
    return effect.hankin_mesh.num_faces();
  }
  static const BakedPaletteBank &palette_bank(const EffectT &effect) {
    return effect.palette_bank.bank;
  }
  /**
   * @brief Replaces an effect's hankin mesh with one built on a graph node.
   * @param effect Effect at any size.
   * @param node Simple-registry index of the base solid.
   */
  template <typename E> static void load_node(E &effect, uint8_t node) {
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = Solids::simple_registry[node].generate(a, b);
      effect.compiled_hankin = CompiledHankin();
      MeshOps::compile_hankin(base, effect.compiled_hankin, target, a);
      effect.hankin_mesh.clear();
      MeshOps::update_hankin(effect.compiled_hankin, effect.hankin_mesh, target,
                             0.0f);
    });
    effect.classify_mesh_topology(effect.hankin_mesh);
  }
  /**
   * @brief Camera of an effect at its current endpoint.
   * @param effect Effect at any size.
   * @return The orientation's latest quaternion.
   */
  template <typename E> static math::Quaternion camera(const E &effect) {
    return effect.orientation.get();
  }
  /**
   * @brief Shades a hankin face from its class's star palette as draw_mesh()
   *        does with every handoff settled.
   * @param effect Effect whose palettes and intensity shade the face.
   * @param face Face index.
   * @param size The face's inradius.
   * @param shader Receives the face's palette and gradient scale.
   */
  template <int W, int H>
  static void shade_face(const HankinSolids<W, H> &effect, size_t face,
                         float size, FacePaletteShader &shader) {
    const int slot =
        MeshPaletteBank::slot_of(effect.hankin_mesh.topology[face]);
    shader.set_palette(
        &effect.palette_bank.bank.entries[effect.palette_idx[slot]].view());
    shader.scale =
        size > math::TOLERANCE ? effect.params.intensity / size : 0.0f;
  }
  /**
   * @brief Draws the moving hankin mesh through a moving camera, recording the
   *        faces scanned against each arm boundary.
   * @param effect Effect whose current mesh is drawn; its camera is replaced.
   * @param previous The mesh's model vertices at the previous endpoint.
   * @param from Camera at the previous endpoint.
   * @param to Camera at the current endpoint.
   * @param copies Receives the arm sides' copies.
   */
  template <int W, int H>
  static void draw_copies(HankinSolids<W, H> &effect,
                          const std::vector<math::Vector> &previous,
                          const math::Quaternion &from,
                          const math::Quaternion &to, FaceCopies &copies) {
    effect.orientation.set(from).push(to);
    const auto &mesh = effect.hankin_mesh;
    std::vector<math::Vector> before = previous;
    {
      Canvas canvas(effect);
      FacePaletteShader fragment;
      int window = -1;
      effect.scan_mesh(
          canvas, mesh, before.data(), false, [&](const MeshState &rotated) {
            window =
                record_face_copies(canvas, rotated, effect.orientation, copies);
            auto setup = [&](size_t f, float size) {
              shade_face(effect, f, size, fragment);
              if (window >= 0)
                copies[window][f].size = size;
            };
            Scan::Mesh::draw_specialized<W, H>(effect.filters, canvas, rotated,
                                               fragment, scratch_arena_a,
                                               nullptr, setup);
          });
    }
    effect.advance_display();
  }
  /**
   * @brief Interlace angle of an effect's hankin mesh.
   * @param effect Effect at any size.
   * @return The angle the last hankin draw used.
   */
  template <typename E> static float angle(const E &effect) {
    return effect.params.hankin_angle;
  }
  /**
   * @brief Rebuilds an effect's hankin mesh at an interlace angle.
   * @param effect Effect at any size.
   * @param angle Interlace angle (radians).
   */
  template <typename E> static void rebuild(E &effect, float angle) {
    MeshOps::update_hankin(effect.compiled_hankin, effect.hankin_mesh,
                           persistent_arena, angle);
  }
  /**
   * @brief Draws an effect's hankin mesh with every handoff settled, the mesh
   *        held over the frame.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @param effect Effect whose mesh, palettes and camera are drawn.
   * @return Displayed frame.
   */
  template <int W, int H>
  static std::vector<Pixel> draw_held(HankinSolids<W, H> &effect) {
    using E = HankinSolids<W, H>;
    const auto &bank = effect.palette_bank.bank.entries;
    const BakedPalette *star[E::NUM_PALETTES];
    for (int s = 0; s < E::NUM_PALETTES; ++s)
      star[s] = &bank[effect.palette_idx[s]].view();
    {
      Canvas canvas(effect);
      effect.draw_mesh(canvas, effect.hankin_mesh, star, star);
    }
    effect.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(effect, out);
    return out;
  }

  /** @brief One morph draw as the draw callback received it. */
  struct ConwayDraw {
    std::vector<math::Vector> vertices; ///< The drawn mesh's vertices.
    std::vector<math::Vector> previous; ///< Its previous vertices, if any.
    bool has_previous = false;          ///< Whether the leg supplied them.
  };
  /** @brief Morph draws recorded by record_conway_draws(). */
  struct ConwayRecorder {
    EffectT *effect = nullptr;     ///< Effect whose draws are forwarded.
    std::vector<ConwayDraw> draws; ///< Recorded draws, oldest first.
  };
  /**
   * @brief Routes an effect's morph draws through a recorder.
   * @param effect Effect whose draw callback is replaced.
   * @param recorder Receives every draw; must outlive the effect's legs.
   */
  static void record_conway_draws(EffectT &effect, ConwayRecorder &recorder) {
    recorder.effect = &effect;
    effect.draw_conway_fn = decltype(effect.draw_conway_fn)(
        [r = &recorder](Canvas &c, MeshState &m,
                        const Animation::OpLeg::Shading &sh) {
          ConwayDraw draw;
          draw.vertices.assign(m.vertices.begin(), m.vertices.end());
          draw.has_previous = sh.previous_vertices != nullptr;
          if (draw.has_previous)
            draw.previous.assign(sh.previous_vertices,
                                 sh.previous_vertices + m.vertices.size());
          r->draws.push_back(std::move(draw));
          r->effect->draw_conway_mesh(c, m, sh);
        });
  }

  static std::vector<Pixel> draw_conway(EffectT &effect,
                                        const Animation::OpLeg::Shading &sh) {
    {
      Canvas canvas(effect);
      effect.draw_conway_mesh(canvas, effect.hankin_mesh, sh);
    }
    effect.advance_display();
    std::vector<Pixel> out;
    capture_frame<SMALL_W, SMALL_H>(effect, out);
    return out;
  }
};

/** @brief Channel bound of a hankin column against its held oracle. */
inline constexpr int HANKIN_HANDOFF_BOUND = 2;

/** @brief Verifies a manual Angle edit freezes an in-flight Conway morph. */
inline void test_hankinsolids_manual_pause_holds_morph() {
  reset_effect_globals();
  HankinSolidsWhiteBox::EffectT effect;
  effect.init();

  int frames = 0;
  while (!HankinSolidsWhiteBox::morph_active(effect) && frames < 100) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
  }
  HS_EXPECT_TRUE(HankinSolidsWhiteBox::morph_active(effect));

  hs_test::run_frames(effect, 5);
  const uint8_t held_node = HankinSolidsWhiteBox::node(effect);
  HS_EXPECT_EQ(effect.updateParameter("Angle", 0.7f), ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.animations_paused());

  hs_test::run_frames(effect, 300);
  HS_EXPECT_TRUE(HankinSolidsWhiteBox::morph_active(effect));
  HS_EXPECT_EQ(HankinSolidsWhiteBox::node(effect), held_node);
  const auto *angle = effect.getParameters().find("Angle");
  HS_EXPECT_TRUE(angle != nullptr);
  if (!angle)
    return;
  HS_EXPECT_NEAR(angle->get(), 0.7f, 1e-6f);

  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(effect);
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief Bounds whole-solid generation, classification and rendering scratch
 * against HankinSolids' exported budgets at the device height.
 */
inline void test_hankinsolids_arena_budget_covers_every_solid() {
  constexpr int W = 288, H = 144;
  constexpr size_t SCRATCH_A = HankinSolids<W, H>::SCRATCH_A_BYTES;
  constexpr size_t SCRATCH_B = HankinSolids<W, H>::SCRATCH_B_BYTES;
  constexpr size_t MEASURE = 1024 * 1024; // headroom so a peak never traps here
  constexpr float ANGLE = math::PI_F / 4.0f;

  auto solids = Solids::Collections::get_simple_solids();
  for (size_t idx = 0; idx < solids.size(); ++idx) {
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * MEASURE, MEASURE, MEASURE);

    MeshPaletteBank palette_bank;
    palette_bank.bake_all(persistent_arena);

    // The effect's held graph-walk seed; dodecahedron is the largest Platonic.
    PolyMesh seed;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed =
          Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), target);
    });

    MeshState mesh;
    CompiledHankin hankin;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = Solids::finalize_solid(solids[idx].generate(a, b), a);
      hankin = CompiledHankin();
      MeshOps::compile_hankin(base, hankin, target, a);
      mesh.clear();
      MeshOps::update_hankin(hankin, mesh, target, ANGLE);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(mesh, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Render peak: transform into scratch_a, then Scan::Mesh::draw stacks a
    // FaceScratchBuffer on top.
    {
      ScratchScope a_guard(scratch_arena_a);
      math::Orientation<> orientation;
      OrientTransformer camera(orientation);
      MeshState rotated;
      MeshOps::transform(mesh, rotated, scratch_arena_a, camera);
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      Pipeline<W, H> filters;
      auto frag = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(1000, 1000, 1000), 1.0f);
      };
      Scan::Mesh::draw<W, H>(filters, canvas, rotated, frag, scratch_arena_a);
    }

    // Morph compaction peak: the CompiledHankin + palette bank survive into
    // scratch_b, the mesh + walk seed into scratch_a, then persistent is
    // reset.
    {
      Persist<CompiledHankin> ph(hankin, scratch_arena_b, persistent_arena);
      Persist<MeshState> pf(mesh, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank, scratch_arena_b,
                                  persistent_arena);
      Persist<PolyMesh> ps(seed, scratch_arena_a, persistent_arena);
      persistent_arena.reset();
    }

    const size_t a_peak = scratch_arena_a.get_high_water_mark();
    const size_t b_peak = scratch_arena_b.get_high_water_mark();
    if (a_peak > SCRATCH_A || b_peak > SCRATCH_B)
      std::printf("  HankinSolids arena OVER BUDGET solid[%zu] '%s': "
                  "scratchA=%zu/%zu scratchB=%zu/%zu\n",
                  idx, solids[idx].name, a_peak, SCRATCH_A, b_peak, SCRATCH_B);
    HS_EXPECT_TRUE(a_peak <= SCRATCH_A);
    HS_EXPECT_TRUE(b_peak <= SCRATCH_B);
  }
  configure_arenas_default();
}

/**
 * @brief Lerps every shape weight between two endpoints.
 * @param from Weights at the previous endpoint.
 * @param to Weights at the current endpoint.
 * @param u Frame fraction.
 * @return The weights at @p u.
 */
inline HankinSolidsDetail::ShapeWeights
lerp_shape_weights(const HankinSolidsDetail::ShapeWeights &from,
                   const HankinSolidsDetail::ShapeWeights &to, float u) {
  const auto lerp = [u](float a, float b) { return a + (b - a) * u; };
  return {lerp(from.strap_open, to.strap_open),
          lerp(from.strap_close, to.strap_close),
          lerp(from.strap_terminal, to.strap_terminal),
          lerp(from.star_close, to.star_close)};
}

/**
 * @brief Each column of a hankin draw shades the strap opening, closing and
 *        crossfade and the star rim blend at its sweep time, on both arms.
 */
inline void test_hankinsolids_handoffs_shade_at_column_time() {
  using WB = HankinSolidsWhiteBox;
  reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  WB::set_angle(effect, 0.6f);
  WB::arm_strap_crossfade(effect);
  for (const int frame : {3, 40, 78}) {
    HS_CONTEXT("cycle frame", frame);
    const WB::FrameHandoff handoff = WB::handoff(frame, frame - 1);
    HS_EXPECT_TRUE(handoff.moves);
    const std::vector<Pixel> drawn =
        WB::draw(effect, handoff, handoff.strap_blend);
    const std::vector<Pixel> expected =
        column_time_oracle<SMALL_W, SMALL_H>([&](float u) {
          const auto weights =
              lerp_shape_weights(handoff.previous, handoff.current, u);
          const float strap =
              handoff.previous_strap_blend +
              (handoff.strap_blend - handoff.previous_strap_blend) * u;
          return WB::draw(effect, {weights, weights, strap, strap, false},
                          strap);
        });
    HS_EXPECT_TRUE(expect_frames_within("hankin handoff", drawn, expected,
                                        HANKIN_HANDOFF_BOUND) > 100);
  }
}

/**
 * @brief A hankin draw's previous endpoint is the cycle frame drawn last,
 *        held while paused and reset at a cycle start and with temporal
 *        sampling off, whose frames keep the history.
 */
inline void test_hankinsolids_handoff_history() {
  using WB = HankinSolidsWhiteBox;
  reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  HS_EXPECT_EQ(WB::drawn_cycle_frame(effect), -1);
  for (int f = 0; f < 4; ++f) {
    const int frame = WB::cycle_frame(effect);
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_EQ(WB::drawn_cycle_frame(effect), frame);
    HS_EXPECT_EQ(WB::cycle_frame(effect), frame + 1);
  }
  HS_EXPECT_EQ(WB::record_cycle_draw(effect, 4), 3);
  HS_EXPECT_EQ(WB::record_cycle_draw(effect, 4), 4);
  Spacetime::set_flat_surface(true);
  HS_EXPECT_EQ(WB::record_cycle_draw(effect, 5), 5);
  Spacetime::set_flat_surface(false);
  HS_EXPECT_EQ(WB::record_cycle_draw(effect, 6), 5);
  HS_EXPECT_FALSE(WB::handoff(6, 6).moves);
}

/**
 * @brief Each column of a Conway-leg draw shades the leg's palette crossfade
 *        at its sweep time.
 */
inline void test_hankinsolids_conway_crossfade_shades_at_column_time() {
  using WB = HankinSolidsWhiteBox;
  using Shading = Animation::OpLeg::Shading;
  reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  WB::set_angle(effect, 0.6f);
  constexpr int RAMPS = 3;
  constexpr float FROM_WEIGHT = 0.3f, TO_WEIGHT = 0.6f;
  const uint8_t ramp_from[RAMPS] = {0, 1, 2};
  const uint8_t ramp_to[RAMPS] = {3, 4, 2};
  const size_t faces = WB::faces(effect);
  std::vector<uint8_t> face_ramp(faces);
  for (size_t f = 0; f < faces; ++f)
    face_ramp[f] = static_cast<uint8_t>(f % RAMPS);
  const auto draw = [&](float weight, float previous_weight, bool moving) {
    ScratchScope blends(scratch_arena_b);
    BakedPalette ramps[RAMPS];
    for (int r = 0; r < RAMPS; ++r)
      ramps[r] = bake_palette_blend(
          scratch_arena_b, WB::palette_bank(effect).entries[ramp_from[r]],
          WB::palette_bank(effect).entries[ramp_to[r]], weight);
    Shading shading{ramps, face_ramp.data(), faces};
    if (moving) {
      shading.bank = &WB::palette_bank(effect);
      shading.ramp_from = ramp_from;
      shading.ramp_to = ramp_to;
    }
    shading.weight = weight;
    shading.previous_weight = previous_weight;
    return WB::draw_conway(effect, shading);
  };
  const std::vector<Pixel> drawn = draw(TO_WEIGHT, FROM_WEIGHT, true);
  const std::vector<Pixel> expected =
      column_time_oracle<SMALL_W, SMALL_H>([&](float u) {
        const float weight = FROM_WEIGHT + (TO_WEIGHT - FROM_WEIGHT) * u;
        return draw(weight, weight, false);
      });
  HS_EXPECT_TRUE(expect_frames_within("conway crossfade", drawn, expected,
                                      HANKIN_HANDOFF_BOUND) > 100);
}

/** @brief Hankin sweep frame the seam test draws: handoffs settled at both
 *  endpoints, the interlace angle near its fastest. */
inline constexpr int HANKIN_SEAM_FRAME = 25;

/**
 * @brief The starting tetrahedron has no excess colour step at capture seams.
 * @details Nearly every face of the tetrahedron's hankin reaches a pole, and
 * its large faces sweep column pair by column pair.
 */
inline void test_hankinsolids_interlace_motion_is_seam_free() {
  constexpr int W = 288, H = 144;
  reset_effect_globals();
  pin_frame_clock(0);
  HankinSolids<W, H> effect;
  effect.init();
  std::vector<Pixel> previous, current;
  MeshSeamSteps steps;
  for (int f = 0; f < 120; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    capture_frame<W, H>(effect, current);
    if (f > 0)
      steps.add(previous, current, W, H);
    previous.swap(current);
  }
  expect_mesh_seam_continuous(steps, "hankin interlace");
  hs::clear_mock_time();
}

/**
 * @brief A hankin face across the arm boundary scans as one polygon, shaded
 *        alike, in one sweep frame's late copy and the next frame's early
 *        copy.
 * @details The interlace angles and cameras are the effect's own over three
 * sweep frames, on the starting tetrahedron's hankin.
 */
inline void test_hankinsolids_seam_copies_are_one_polygon() {
  using WB = HankinSolidsWhiteBox;
  constexpr int W = 288, H = 144;
  reset_effect_globals();
  pin_frame_clock(0);
  HankinSolids<W, H> effect;
  effect.init();
  std::array<float, 3> angles{};
  std::array<math::Quaternion, 3> cameras{};
  for (int f = 0; f <= HANKIN_SEAM_FRAME + 1; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    if (f >= HANKIN_SEAM_FRAME - 1) {
      angles[f - HANKIN_SEAM_FRAME + 1] = WB::angle(effect);
      cameras[f - HANKIN_SEAM_FRAME + 1] = WB::camera(effect);
    }
  }
  hs::clear_mock_time();
  WB::load_node(effect, ConwayGraph::TETRAHEDRON);
  std::array<std::vector<math::Vector>, 3> model;
  for (int j = 0; j < 3; ++j) {
    WB::rebuild(effect, angles[j]);
    model[j].assign(WB::mesh_vertices(effect).begin(),
                    WB::mesh_vertices(effect).end());
  }
  std::array<FaceCopies, 2> frames;
  float arc = 0.0f, turn = 0.0f;
  for (int j = 0; j < 2; ++j) {
    WB::rebuild(effect, angles[j + 1]);
    WB::draw_copies(effect, model[j], cameras[j], cameras[j + 1], frames[j]);
    arc = fmaxf(arc, Scan::max_vertex_arc(model[j].data(), model[j + 1].data(),
                                          model[j].size()));
    turn = fmaxf(turn, quaternion_angle(cameras[j], cameras[j + 1]));
  }
  HS_EXPECT_TRUE(Scan::mesh_motion_active<W>(arc));
  const auto path = [&](int frame, uint16_t i, float u) {
    return (model[frame][i] + (model[frame + 1][i] - model[frame][i]) * u)
        .normalized();
  };
  const size_t faces = expect_seams_match<W, H>(
      effect, "hankin seam copies", frames, cameras[1],
      [&](uint16_t i) {
        return seam_vertex_tolerance<W>([&](float u) { return path(0, i, u); },
                                        [&](float u) { return path(1, i, u); },
                                        cameras[1], turn);
      },
      [&](size_t f, float size, FacePaletteShader &shader) {
        WB::shade_face(effect, f, size, shader);
      });
  HS_EXPECT_GE(faces, size_t{4});
}

/**
 * @brief A morph leg hands each draw the vertices of the frame it drew last,
 *        none on its first draw or on a redraw of the same frame.
 */
inline void test_hankinsolids_morph_draws_carry_previous_vertices() {
  using WB = HankinSolidsWhiteBox;
  reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  WB::ConwayRecorder recorder;
  WB::record_conway_draws(effect, recorder);
  int frames = 0;
  while (recorder.draws.size() < 8 && frames < 400) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
  }
  HS_EXPECT_GE(recorder.draws.size(), size_t{8});
  if (recorder.draws.size() < 8)
    return;
  HS_EXPECT_FALSE(recorder.draws[0].has_previous);
  bool moved = false;
  for (size_t d = 1; d < recorder.draws.size(); ++d) {
    HS_CONTEXT("draw", static_cast<long long>(d));
    const WB::ConwayDraw &draw = recorder.draws[d];
    HS_EXPECT_TRUE(draw.has_previous);
    HS_EXPECT_TRUE(draw.previous == recorder.draws[d - 1].vertices);
    moved |= draw.previous != draw.vertices;
  }
  HS_EXPECT_TRUE(moved);
  HS_EXPECT_EQ(effect.updateParameter("Angle", 0.3f), ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_FALSE(recorder.draws.back().has_previous);
}

/**
 * @brief A morph draw moves its mesh from the previous vertices the leg
 *        hands it to the current ones across each arm, free of the seam the
 *        held mesh shows.
 */
inline void test_hankinsolids_morph_motion_is_seam_free() {
  using WB = HankinSolidsWhiteBox;
  using Shading = Animation::OpLeg::Shading;
  reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  constexpr float FROM = 0.5f, TO = 0.7f;
  WB::rebuild(effect, FROM);
  const size_t count = WB::mesh_vertices(effect).size();
  std::vector<math::Vector> from(WB::mesh_vertices(effect).begin(),
                                 WB::mesh_vertices(effect).end());
  WB::rebuild(effect, TO);
  const size_t faces = WB::faces(effect);
  const BakedPalette ramp = WB::palette_bank(effect).entries[0];
  std::vector<uint8_t> face_ramp(faces, 0);
  const auto draw = [&](math::Vector *previous) {
    Shading shading{&ramp, face_ramp.data(), faces};
    shading.previous_vertices = previous;
    return WB::draw_conway(effect, shading);
  };
  std::vector<math::Vector> previous = from;
  const std::vector<Pixel> drawn = draw(previous.data());
  const std::vector<Pixel> held = draw(nullptr);
  HS_EXPECT_EQ(count, from.size());
  expect_seam_resolved<SMALL_W, SMALL_H>(
      "hankin morph", drawn, held, [&](float u) {
        WB::rebuild(effect, FROM + (TO - FROM) * u);
        return draw(nullptr);
      });
}
