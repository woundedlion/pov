/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Hankin-sweep displacement and face-normal bounds.
// ---------------------------------------------------------------------------

using MeshOps::HankinSolve;

/** Solves each dynamic vertex and retains its branch diagnostics. */
inline void hankin_solve(const CompiledHankin &compiled, float angle,
                         std::vector<HankinSolve> &out) {
  const bool is_flat = std::abs(angle) < math::TOLERANCE;
  const float cos_ha = cosf(angle * 0.5f);
  const float sin_ha = sinf(angle * 0.5f);
  out.clear();
  out.reserve(compiled.dynamic_instructions.size());
  for (const auto &instr : compiled.dynamic_instructions)
    out.push_back(
        MeshOps::solve_hankin_vertex(compiled, instr, is_flat, cos_ha, sin_ha));
}

/**
 * @brief Newell normals of every compiled hankin face for a solved star set.
 * @param compiled Baked hankin topology.
 * @param dyn Solved dynamic vertices.
 * @param out Filled with one unnormalized Newell normal per face.
 */
inline void hankin_face_normals(const CompiledHankin &compiled,
                                const std::vector<HankinSolve> &dyn,
                                std::vector<math::Vector> &out) {
  auto vertex_at = [&](uint16_t idx) {
    return idx < compiled.static_offset ? compiled.static_vertices[idx]
                                        : dyn[idx - compiled.static_offset].pos;
  };
  out.assign(compiled.face_counts.size(), math::Vector());
  size_t base = 0;
  for (size_t f = 0; f < compiled.face_counts.size(); ++f) {
    const size_t n = compiled.face_counts[f];
    out[f] = newell_normal(static_cast<int>(n), [&](int k) {
      return vertex_at(compiled.faces[base + k]);
    });
    base += n;
  }
}

/** Newell magnitude (twice the face area) below which a face's normal
 * direction is numerical noise and its sign is not evidence of a fold. */
inline constexpr float HANKIN_FLAT_FACE = 1e-4f;

/** @brief Per-step stability metrics of one sweep sample. */
struct HankinStepStats {
  float max_disp = 0;  /**< Largest single-vertex chord vs the previous step. */
  float mean_disp = 0; /**< Mean dynamic-vertex chord vs the previous step. */
  int normal_flips =
      0;              /**< Non-degenerate faces whose Newell normal reversed. */
  int flat_faces = 0; /**< Faces below HANKIN_FLAT_FACE this step. */
  float max_corner_chord = 0; /**< Largest chord(star point, its corner). */
};

/**
 * @brief Fills @p stats from consecutive solved states.
 * @param compiled Baked hankin topology.
 * @param prev Previous step's solve (empty for the first step).
 * @param prev_normals Previous step's face normals.
 * @param curr Current step's solve.
 * @param curr_normals Current step's face normals.
 * @param stats Metrics to fill.
 */
inline void hankin_step_stats(const CompiledHankin &compiled,
                              const std::vector<HankinSolve> &prev,
                              const std::vector<math::Vector> &prev_normals,
                              const std::vector<HankinSolve> &curr,
                              const std::vector<math::Vector> &curr_normals,
                              HankinStepStats &stats) {
  for (size_t i = 0; i < curr.size(); ++i) {
    const math::Vector cn = math::normalized_or(
        compiled.corner(compiled.dynamic_instructions[i].v_corner),
        curr[i].pos);
    stats.max_corner_chord = hs_test::fold_worst(
        stats.max_corner_chord, (curr[i].pos - cn).magnitude());
  }
  for (const math::Vector &n : curr_normals)
    if (n.magnitude() < HANKIN_FLAT_FACE)
      ++stats.flat_faces;
  if (prev.empty())
    return;
  double sum = 0;
  for (size_t i = 0; i < curr.size(); ++i) {
    const float d = (curr[i].pos - prev[i].pos).magnitude();
    sum += d;
    stats.max_disp = hs_test::fold_worst(stats.max_disp, d);
  }
  stats.mean_disp = static_cast<float>(sum / curr.size());
  for (size_t f = 0; f < curr_normals.size(); ++f) {
    if (prev_normals[f].magnitude() < HANKIN_FLAT_FACE ||
        curr_normals[f].magnitude() < HANKIN_FLAT_FACE)
      continue;
    if (math::dot(prev_normals[f], curr_normals[f]) < 0)
      ++stats.normal_flips;
  }
}

/** @brief Sweep-wide roll-up of the per-step metrics. */
struct HankinSweepSummary {
  int total_normal_flips = 0;
  float worst_max_disp = 0;
  float worst_mean_disp = 0;
  float spike_ratio = 0;    /**< worst_max_disp / mean_disp at that step. */
  int worst_flat_faces = 0; /**< Most sub-HANKIN_FLAT_FACE faces in a step. */
  float worst_corner_chord =
      0; /**< Largest chord(star point, corner) over the sweep. */
};

/**
 * @brief Rolls the per-step table into a summary.
 * @param table Per-step metrics, index 0 being the first (baseline) sample.
 * @return The roll-up.
 */
inline HankinSweepSummary
hankin_summarize(const std::vector<HankinStepStats> &table) {
  HankinSweepSummary sum;
  for (const HankinStepStats &r : table) {
    sum.worst_flat_faces = std::max(sum.worst_flat_faces, r.flat_faces);
    sum.worst_corner_chord =
        hs_test::fold_worst(sum.worst_corner_chord, r.max_corner_chord);
  }
  for (size_t s = 1; s < table.size(); ++s) {
    sum.total_normal_flips += table[s].normal_flips;
    sum.worst_mean_disp =
        hs_test::fold_worst(sum.worst_mean_disp, table[s].mean_disp);
    if (std::isnan(table[s].max_disp) ||
        table[s].max_disp > sum.worst_max_disp) {
      sum.worst_max_disp = table[s].max_disp;
      sum.spike_ratio = table[s].mean_disp > 0
                            ? table[s].max_disp / table[s].mean_disp
                            : 0.0f;
    }
  }
  return sum;
}

/**
 * @brief Measures the drawn sweep's stability over the sampled dodecahedron,
 *        dodecahedron_hk62_ambo, octahedron and octahedron_hk17_ambo prefixes.
 * @details The gates bound displacement and face-normal reversals of the
 * frames the leg draws, which land on the snorm16 arrival vertices.
 */
inline void test_hankin_sweep_vertex_stability() {
  const auto split = islamic_split(IslamicStars<288, 144>::GENERATED_BUDGET);
  constexpr int SAMPLES = 32;

  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    HS_CONTEXT(site.name);
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    const PolyMesh seed = finalize_site_seed(site, persist);

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    CompiledHankin compiled;
    MeshOps::compile_hankin(seed, compiled, a, b);

    std::vector<HankinSolve> arrival;
    hankin_solve(compiled, site.theta_star, arrival);

    Arena leg_arena(morph_aux_buf, sizeof(morph_aux_buf));
    Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
    MeshPaletteBank bank;
    bank.bake_all(bank_arena);
    std::vector<uint8_t> palettes(seed.face_counts.size(), 0);
    Animation::OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                             .prev_face_palette =
                                                 palettes.data(),
                                             .prev_faces = palettes.size()};
    LegDrawProbe draw_probe;
    std::vector<std::vector<math::Vector>> drawn_frames;
    auto draw = [&](Canvas &, const MeshState &mesh,
                    const Animation::OpLeg::Shading &shading) {
      draw_probe.observe(mesh, shading);
      drawn_frames.push_back(draw_probe.prev_v);
    };
    Animation::OpLeg leg(
        seed,
        Animation::OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                          .theta_end = site.theta_star,
                                          .sweep_frames = SAMPLES - 1},
        leg_arena, draw, handoff);
    draw_probe.ramp_count = leg.landing().blend_pairs;
    hs_test::StubEffect effect(288, 144);
    for (int sample = 0; sample < SAMPLES; ++sample) {
      {
        Canvas canvas(effect);
        if (sample == 0)
          leg.step_paused(canvas);
        else
          leg.step(canvas);
      }
      effect.advance_display();
    }
    HS_EXPECT_SIZE_OR_RETURN(drawn_frames, static_cast<size_t>(SAMPLES));

    std::vector<HankinStepStats> table;
    std::vector<HankinSolve> prev, curr;
    std::vector<math::Vector> prev_normals, curr_normals;
    for (int sample = 0; sample < SAMPLES; ++sample) {
      HS_EXPECT_SIZE_OR_RETURN(drawn_frames[sample],
                               compiled.static_vertices.size() +
                                   arrival.size());
      curr.assign(arrival.size(), HankinSolve{});
      for (size_t i = 0; i < arrival.size(); ++i)
        curr[i] = {drawn_frames[sample][compiled.static_vertices.size() + i],
                   arrival[i].branch, 0.0f};
      hankin_face_normals(compiled, curr, curr_normals);
      HankinStepStats row;
      hankin_step_stats(compiled, prev, prev_normals, curr, curr_normals, row);
      table.push_back(row);
      prev = curr;
      prev_normals = curr_normals;
    }
    const HankinSweepSummary sum = hankin_summarize(table);
    HS_EXPECT_EQ(sum.total_normal_flips, 0);
    HS_EXPECT_LE(sum.worst_max_disp, 0.05f);
    HS_EXPECT_LE(sum.worst_mean_disp, 0.03f);
    HS_EXPECT_LE(sum.spike_ratio, 1.5f);
    HS_EXPECT_EQ(sum.worst_flat_faces, 0);
    HS_EXPECT_LE(sum.worst_corner_chord, 0.50f);

    PolyMesh bookend;
    Animation::OpLeg::arrival_mesh(leg.landing(), bookend, b);
    HS_EXPECT_EQ(drawn_frames.back().size(), bookend.vertices.size());
    if (drawn_frames.back().size() == bookend.vertices.size())
      for (size_t vertex = 0; vertex < bookend.vertices.size(); ++vertex)
        HS_EXPECT_VEC(drawn_frames.back()[vertex], bookend.vertices[vertex],
                      1e-6f);
    for (size_t i = 0; i < arrival.size(); ++i)
      HS_EXPECT_LE((drawn_frames.back()[compiled.static_vertices.size() + i] -
                    arrival[i].pos)
                       .magnitude(),
                   3e-5f);
  }
}

/**
 * @brief Smoke-tests an OpLeg HANKIN_SWEEP on every sweep seed: construction
 *        populates the landing (star faces first, in base-face order), and every
 *        step hands the draw callback a compiled
 *        mesh with the constant hankin face count and in-range ramp indices.
 */
inline void test_opleg_hankin_sweep_smoke() {
  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    const auto split = islamic_split(IslamicStars<288, 144>::GENERATED_BUDGET);
    hs::random().seed(2026u);

    Arena leg(morph_target_buf, sizeof(morph_target_buf));
    Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));

    MeshPaletteBank bank;
    bank.bake_all(bank_arena);

    Arena seed_scratch(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh seed = site.seed(leg, seed_scratch);
    std::vector<uint8_t> pal(seed.face_counts.size());
    for (size_t f = 0; f < seed.face_counts.size(); ++f)
      pal[f] = static_cast<uint8_t>(f % Animation::OpLeg::PALETTES);

    Animation::OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                             .prev_face_palette = pal.data(),
                                             .prev_faces =
                                                 seed.face_counts.size()};

    constexpr float MAX_STEP_CHORD = 0.15f;
    LegDrawProbe probe;
    auto cb = [&](Canvas &, const MeshState &m,
                  const Animation::OpLeg::Shading &sh) {
      probe.observe(m, sh);
    };

    constexpr int SWEEP = 8;
    Animation::OpLeg anim(
        seed,
        Animation::OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                          .theta_end = site.theta_star,
                                          .sweep_frames = SWEEP},
        leg, cb, handoff);

    const Animation::OpLeg::Landing &landing = anim.landing();
    probe.ramp_count = landing.blend_pairs;
    HS_EXPECT_EQ(landing.primary_faces, seed.face_counts.size());
    HS_EXPECT_EQ(landing.faces, seed.face_counts.size() + seed.vertices.size());
    HS_EXPECT_TRUE(landing.topology != nullptr);

    hs_test::StubEffect fx(288, 144);
    hs_test::drive_leg(anim, fx, SWEEP);
    HS_EXPECT_EQ(probe.drawn, (size_t)SWEEP);
    HS_EXPECT_EQ(probe.faces, landing.faces);
    HS_EXPECT_LT(probe.worst_step, MAX_STEP_CHORD);
    std::printf("  [opleg hankin] worst per-frame vertex step %.4f chord "
                "(bound %.2f)\n",
                (double)probe.worst_step, (double)MAX_STEP_CHORD);

    Arena arrival_arena(morph_aux_buf, sizeof(morph_aux_buf));
    PolyMesh arrival;
    Animation::OpLeg::arrival_mesh(landing, arrival, arrival_arena);
    std::vector<math::Vector> seed_centers(seed.face_counts.size());
    for (size_t f = 0, off = 0; f < seed.face_counts.size();
         off += seed.face_counts[f], ++f)
      seed_centers[f] = face_centroid_unit(seed, off, seed.face_counts[f]);
    HS_EXPECT_GE(arrival.face_counts.size(), landing.primary_faces);
    for (size_t f = 0, off = 0;
         f < landing.primary_faces && f < arrival.face_counts.size();
         off += arrival.face_counts[f], ++f) {
      const math::Vector c =
          face_centroid_unit(arrival, off, arrival.face_counts[f]);
      size_t best = 0;
      float best_d = INFINITY, runner_d = INFINITY;
      for (size_t g = 0; g < seed_centers.size(); ++g) {
        const float d = (seed_centers[g] - c).magnitude();
        if (d < best_d) {
          runner_d = best_d;
          best_d = d;
          best = g;
        } else if (d < runner_d) {
          runner_d = d;
        }
      }
      HS_CONTEXT(site.name, f);
      HS_EXPECT_EQ(best, f);
      HS_EXPECT_LT(best_d, 0.5f * runner_d);
    }
  }
}
