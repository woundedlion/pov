/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Hankin-sweep stability probe: per-step branch, displacement and face-normal
// diagnostics for the hankin legs.
// ---------------------------------------------------------------------------

/** @brief Branch update_hankin took for one dynamic vertex. */
enum class HankinBranch : uint8_t {
  INTERSECT, /**< Contact-plane intersection accepted. */
  COLLAPSED, /**< is_flat or zero-length edge: snapped to the corner. */
  FALLBACK, /**< Degenerate or far near-parallel intersection: edge-midpoint mean. */
  BLENDED, /**< Partial fallback: 0 < fallback_blend < 1. */
};

/** @brief One dynamic vertex's solved position and the branch that made it. */
struct HankinSolve {
  math::Vector pos;
  HankinBranch branch;
  /** dist^2(star, corner) / max(dist^2(m, corner)); full fallback requires
   * far_ratio >= STAR_FAR_RATIO_SQ and near-parallel contact planes
   * (parallel_gate = 1). Zero on degenerate and collapsed branches. */
  float far_ratio = 0;
};

/** @brief Mirror of update_hankin's local blend_above (core/mesh/hankin.h). */
inline float hankin_blend_above(float value, float start, float end) {
  const float t =
      std::max(0.0f, std::min(1.0f, (value - start) / (end - start)));
  return math::quintic_kernel(t);
}

/**
 * @brief Mirrors MeshOps::update_hankin's per-vertex solve, exposing the
 *        branch each dynamic vertex takes.
 * @param compiled Baked hankin topology.
 * @param angle Contact angle in radians.
 * @param out Filled with one entry per dynamic instruction.
 * @details Reproduces the regularization anchor, the plane_cross_sq parallel
 *   gate and the fallback blend; positions must track update_hankin exactly.
 */
inline void hankin_solve(const CompiledHankin &compiled, float angle,
                         std::vector<HankinSolve> &out) {
  const bool is_flat = std::abs(angle) < math::TOLERANCE;
  const float cos_ha = cosf(angle * 0.5f);
  const float sin_ha = sinf(angle * 0.5f);

  out.assign(compiled.dynamic_instructions.size(), HankinSolve{});
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const HankinInstruction &instr = compiled.dynamic_instructions[i];
    const math::Vector p_corner = compiled.corner(instr.v_corner);
    const math::Vector cn = math::normalized_or(p_corner, p_corner);

    if (is_flat) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const math::Vector m1 = compiled.static_vertices[instr.idx_m1];
    const math::Vector m2 = compiled.static_vertices[instr.idx_m2];
    const math::Vector cross1 =
        math::cross(compiled.corner(instr.v_prev), p_corner);
    const math::Vector cross2 =
        math::cross(p_corner, compiled.corner(instr.v_next));
    if (math::dot(cross1, cross1) < math::EPS_CROSS_SQ ||
        math::dot(cross2, cross2) < math::EPS_CROSS_SQ) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const math::Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y,
                              sin_ha * m1.z);
    const math::Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y,
                              -sin_ha * m2.z);
    math::Vector intersect = math::cross(math::rotate(cross1.normalized(), q1),
                                         math::rotate(cross2.normalized(), q2));
    const float plane_cross_sq = math::dot(intersect, intersect);

    math::Vector fallback = math::normalized_or(m1 + m2, cn);
    if (math::dot(fallback, p_corner) < 0.0f)
      fallback = -fallback;
    const math::Vector oriented_intersect =
        intersect * math::dot(intersect, p_corner);
    const math::Vector raw_star =
        math::normalized_or(oriented_intersect, fallback);
    const float local_sq = std::max(math::distance_squared(m1, cn),
                                    math::distance_squared(m2, cn));
    if (!(local_sq > math::EPS_LEN_SQ)) {
      out[i] = {fallback, HankinBranch::FALLBACK, 0.0f};
      continue;
    }

    const float raw_ratio_sq = math::distance_squared(raw_star, cn) / local_sq;
    const float conditioned = hankin_blend_above(
        raw_ratio_sq, MeshOps::HANKIN_CONDITIONED_NEAR_RATIO_SQ,
        MeshOps::HANKIN_CONDITIONED_FAR_RATIO_SQ);
    const float anchor =
        std::max(0.0f,
                 MeshOps::HANKIN_PARALLEL_REGULARIZATION_SQ - plane_cross_sq) +
        conditioned * std::max(0.0f, MeshOps::HANKIN_CONDITIONED_CLEAR_SQ -
                                         plane_cross_sq);
    intersect =
        math::normalized_or(oriented_intersect + fallback * anchor, fallback);

    const float far_ratio = math::distance_squared(intersect, cn) / local_sq;
    const float parallel_gate = hankin_blend_above(
        -plane_cross_sq, -MeshOps::HANKIN_PARALLEL_GATE_HI_SQ,
        -MeshOps::HANKIN_PARALLEL_GATE_LO_SQ);
    const float fallback_blend =
        hankin_blend_above(far_ratio, MeshOps::STAR_FAR_BLEND_START_RATIO_SQ,
                           MeshOps::STAR_FAR_RATIO_SQ) *
        parallel_gate;
    if (fallback_blend >= 1.0f) {
      out[i] = {fallback, HankinBranch::FALLBACK, far_ratio};
      continue;
    }
    if (fallback_blend > 0.0f) {
      out[i] = {math::normalized_or(intersect * (1.0f - fallback_blend) +
                                        fallback * fallback_blend,
                                    fallback),
                HankinBranch::BLENDED, far_ratio};
      continue;
    }
    out[i] = {intersect, HankinBranch::INTERSECT, far_ratio};
  }
}

/** Chord tolerance between hankin_solve and MeshOps::update_hankin; the two
 * evaluate identical expressions under independently contracted float
 * arithmetic, so bitwise agreement is not guaranteed. */
inline constexpr float HANKIN_MIRROR_TOL = 1e-5f;

/**
 * @brief Pins hankin_solve's star points to MeshOps::update_hankin's.
 * @param compiled Baked hankin topology.
 * @param angle Contact angle in radians.
 * @param arena Scratch arena backing the reference mesh; rewound on return.
 * @return Largest chord between a mirrored and a shipping star point.
 */
inline float hankin_check_mirror(CompiledHankin &compiled, float angle,
                                 Arena &arena) {
  ScratchScope guard(arena);
  std::vector<HankinSolve> mirror;
  hankin_solve(compiled, angle, mirror);
  PolyMesh ref;
  MeshOps::update_hankin(compiled, ref, arena, angle);
  const size_t base = compiled.static_vertices.size();
  HS_EXPECT_EQ(ref.vertices.size(), base + mirror.size());
  float max_chord = 0;
  for (size_t i = 0; i < mirror.size() && base + i < ref.vertices.size(); ++i)
    max_chord = hs_test::fold_worst(
        max_chord, (ref.vertices[base + i] - mirror[i].pos).magnitude());
  HS_EXPECT_LT(max_chord, HANKIN_MIRROR_TOL);
  return max_chord;
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
  float theta = 0; /**< Contact angle of this sample, radians. */
  int branch_flips =
      0; /**< Vertices whose branch changed vs the previous step. */
  float max_disp = 0;  /**< Largest single-vertex chord vs the previous step. */
  float mean_disp = 0; /**< Mean dynamic-vertex chord vs the previous step. */
  int normal_flips =
      0;              /**< Non-degenerate faces whose Newell normal reversed. */
  int flat_faces = 0; /**< Faces below HANKIN_FLAT_FACE this step. */
  float max_far_ratio = 0;    /**< Largest far_ratio this step. */
  float max_corner_chord = 0; /**< Largest chord(star point, its corner). */
};

/**
 * @brief Fills @p stats from consecutive solved states.
 * @param compiled Baked hankin topology.
 * @param prev Previous step's solve (empty for the first step).
 * @param prev_normals Previous step's face normals.
 * @param curr Current step's solve.
 * @param curr_normals Current step's face normals.
 * @param stats Metrics to fill; theta must already be set.
 */
inline void hankin_step_stats(const CompiledHankin &compiled,
                              const std::vector<HankinSolve> &prev,
                              const std::vector<math::Vector> &prev_normals,
                              const std::vector<HankinSolve> &curr,
                              const std::vector<math::Vector> &curr_normals,
                              HankinStepStats &stats) {
  for (size_t i = 0; i < curr.size(); ++i) {
    stats.max_far_ratio =
        hs_test::fold_worst(stats.max_far_ratio, curr[i].far_ratio);
    const math::Vector cn = math::normalized_or(
        compiled.base_vertices[compiled.dynamic_instructions[i].v_corner],
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
    if (curr[i].branch != prev[i].branch)
      ++stats.branch_flips;
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
  int total_branch_flips = 0;
  int total_normal_flips = 0;
  int steps_with_normal_flips = 0;
  float worst_max_disp = 0;
  float worst_mean_disp = 0;
  int worst_step = 0;        /**< Step index owning worst_max_disp. */
  float spike_ratio = 0;     /**< worst_max_disp / mean_disp at that step. */
  int worst_flat_faces = 0;  /**< Most sub-HANKIN_FLAT_FACE faces in a step. */
  float worst_far_ratio = 0; /**< Largest far_ratio over the sweep. */
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
    sum.worst_far_ratio =
        hs_test::fold_worst(sum.worst_far_ratio, r.max_far_ratio);
    sum.worst_corner_chord =
        hs_test::fold_worst(sum.worst_corner_chord, r.max_corner_chord);
  }
  for (size_t s = 1; s < table.size(); ++s) {
    sum.total_branch_flips += table[s].branch_flips;
    sum.total_normal_flips += table[s].normal_flips;
    if (table[s].normal_flips > 0)
      ++sum.steps_with_normal_flips;
    sum.worst_mean_disp =
        hs_test::fold_worst(sum.worst_mean_disp, table[s].mean_disp);
    if (std::isnan(table[s].max_disp) ||
        table[s].max_disp > sum.worst_max_disp) {
      sum.worst_max_disp = table[s].max_disp;
      sum.worst_step = static_cast<int>(s);
      sum.spike_ratio = table[s].mean_disp > 0
                            ? table[s].max_disp / table[s].mean_disp
                            : 0.0f;
    }
  }
  return sum;
}

/**
 * @brief Measures sweep stability over the sampled dodecahedron,
 *        dodecahedron_hk62_ambo, octahedron and octahedron_hk17_ambo prefixes
 *        under the shipping slerp-from-corner parameterization.
 * @details Re-solve modes are diagnostic comparisons. The slerp gates bound
 * displacement and face-normal reversals using snorm16 arrival vertices and
 * the production K_EPS opening fraction.
 */
inline void test_hankin_sweep_vertex_stability() {
  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  constexpr int SAMPLES = 32;
  constexpr float THETA_EPS = Animation::OpLeg::THETA_EPS;

  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    CompiledHankin compiled;
    MeshOps::compile_hankin(seed, compiled, a, b);

    std::vector<HankinSolve> collapsed, arrival;
    hankin_solve(compiled, 0.0f, collapsed);
    hankin_solve(compiled, site.theta_star, arrival);
    std::vector<math::Vector> packed_arrival;
    for (const auto &point : arrival)
      packed_arrival.push_back(math::Snorm3::encode(point.pos).decode());

    // Every metric reads hankin_solve; pin it to the shipping solver first.
    float mirror_chord =
        std::max(hankin_check_mirror(compiled, 0.0f, b),
                 hankin_check_mirror(compiled, site.theta_star, b));
    for (int s = 0; s < SAMPLES; ++s) {
      const float u = static_cast<float>(s) / (SAMPLES - 1);
      const float span = site.theta_star - THETA_EPS;
      mirror_chord = std::max(
          {mirror_chord, hankin_check_mirror(compiled, THETA_EPS + span * u, b),
           hankin_check_mirror(
               compiled, THETA_EPS + span * math::ease_in_out_sin(u), b)});
    }

    // Guard headroom: the far-star guard is scaled by the corner's local edge
    // scale, so on coarse seeds STAR_FAR_RATIO_SQ * local_sq can exceed the
    // 4.0 maximum squared chord, making the guard unreachable for that vertex.
    int unreachable = 0;
    float max_local_sq = 0;
    for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
      const HankinInstruction &instr = compiled.dynamic_instructions[i];
      const math::Vector cn =
          math::normalized_or(compiled.base_vertices[instr.v_corner],
                              compiled.base_vertices[instr.v_corner]);
      const float local_sq = std::max(
          math::distance_squared(compiled.static_vertices[instr.idx_m1], cn),
          math::distance_squared(compiled.static_vertices[instr.idx_m2], cn));
      max_local_sq = std::max(max_local_sq, local_sq);
      if (MeshOps::STAR_FAR_RATIO_SQ * local_sq >= 4.0f)
        ++unreachable;
    }

    std::printf("  [hankin-stability] %s: theta* = %.1f deg, base F=%zu, "
                "dyn V=%zu, hankin F=%zu, guard-unreachable %d/%zu "
                "(max local_sq %.4f)\n",
                site.name, site.theta_star * 180.0f / math::PI_F,
                seed.face_counts.size(), collapsed.size(),
                compiled.face_counts.size(), unreachable, collapsed.size(),
                max_local_sq);
    std::printf(
        "      mirror vs update_hankin: max chord %.3e over %d angles\n",
        mirror_chord, 2 * SAMPLES + 2);

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

    // Three parameterizations over the same sample grid.
    std::vector<HankinStepStats> tables[3];
    std::vector<std::vector<math::Vector>> resolve_pos(SAMPLES);
    float path_dev = 0;
    const char *mode_name[3] = {"resolve-uniform", "resolve-eased",
                                "slerp-eased    "};
    for (int mode = 0; mode < 3; ++mode) {
      std::vector<HankinSolve> prev, curr;
      std::vector<math::Vector> prev_normals, curr_normals;
      for (int s = 0; s < SAMPLES; ++s) {
        const float u = static_cast<float>(s) / (SAMPLES - 1);
        const float k = mode == 0 ? u : math::ease_in_out_sin(u);
        HankinStepStats row;
        if (mode == 2) {
          const float shipping_k =
              Animation::OpLeg::K_EPS + (1.0f - Animation::OpLeg::K_EPS) * k;
          row.theta = site.theta_star;
          curr.assign(arrival.size(), HankinSolve{});
          for (size_t i = 0; i < arrival.size(); ++i) {
            curr[i] = {shipping_k >= 1.0f
                           ? packed_arrival[i].normalized()
                           : math::slerp(collapsed[i].pos.normalized(),
                                         packed_arrival[i], shipping_k),
                       arrival[i].branch, 0.0f};
            HS_EXPECT_VEC(drawn_frames[s][compiled.static_vertices.size() + i],
                          curr[i].pos, 1e-6f);
            path_dev = std::max(path_dev,
                                (curr[i].pos - resolve_pos[s][i]).magnitude());
          }
        } else {
          row.theta = THETA_EPS + (site.theta_star - THETA_EPS) * k;
          hankin_solve(compiled, row.theta, curr);
          if (mode == 1) {
            resolve_pos[s].resize(curr.size());
            for (size_t i = 0; i < curr.size(); ++i)
              resolve_pos[s][i] = curr[i].pos;
          }
        }
        hankin_face_normals(compiled, curr, curr_normals);
        hankin_step_stats(compiled, prev, prev_normals, curr, curr_normals,
                          row);
        tables[mode].push_back(row);
        prev = curr;
        prev_normals = curr_normals;
      }
    }

    for (int mode = 0; mode < 3; ++mode) {
      const HankinSweepSummary sum = hankin_summarize(tables[mode]);
      std::printf("      %s  branch_flips=%d normal_flips=%d (in %d steps) "
                  "worst_max=%.5f @s%d spike=%.1fx worst_mean=%.5f "
                  "flat<=%d far<=%.1f corner<=%.3f\n",
                  mode_name[mode], sum.total_branch_flips,
                  sum.total_normal_flips, sum.steps_with_normal_flips,
                  sum.worst_max_disp, sum.worst_step, sum.spike_ratio,
                  sum.worst_mean_disp, sum.worst_flat_faces,
                  sum.worst_far_ratio, sum.worst_corner_chord);
      if (mode == 2) {
        HS_CONTEXT(site.name, mode);
        HS_EXPECT_EQ(sum.total_normal_flips, 0);
        HS_EXPECT_LE(sum.worst_max_disp, 0.05f);
        HS_EXPECT_LE(sum.worst_mean_disp, 0.03f);
        HS_EXPECT_LE(sum.spike_ratio, 1.5f);
        HS_EXPECT_EQ(sum.worst_flat_faces, 0);
        HS_EXPECT_LE(sum.worst_corner_chord, 0.50f);
      }
    }

    std::printf("      path deviation slerp vs resolve at equal k: max=%.5f\n",
                path_dev);

    PolyMesh bookend;
    Animation::OpLeg::arrival_mesh(leg.landing(), bookend, b);
    HS_EXPECT_EQ(drawn_frames.back().size(), bookend.vertices.size());
    if (drawn_frames.back().size() == bookend.vertices.size())
      for (size_t vertex = 0; vertex < bookend.vertices.size(); ++vertex)
        HS_EXPECT_VEC(drawn_frames.back()[vertex], bookend.vertices[vertex],
                      1e-6f);
    for (size_t i = 0; i < arrival.size(); ++i)
      HS_EXPECT_LE(
          (packed_arrival[i].normalized() - arrival[i].pos).magnitude(), 3e-5f);

    // Opening bookend: chord between the collapsed form and the leg's first
    // drawn angle, in sphere radii (sub-pixel iff below ~1/display radius).
    float eps_chord = 0;
    std::vector<HankinSolve> at_eps;
    hankin_solve(compiled, THETA_EPS, at_eps);
    for (size_t i = 0; i < at_eps.size(); ++i)
      eps_chord =
          std::max(eps_chord, (at_eps[i].pos - collapsed[i].pos).magnitude());
    std::printf("      theta_eps=%.3f opening chord max=%.6f radii "
                "(%.3f px at r=64)\n",
                THETA_EPS, eps_chord, eps_chord * 64.0f);

    std::vector<HankinSolve> probe(arrival.size());
    for (size_t i = 0; i < arrival.size(); ++i)
      probe[i] = {math::slerp(collapsed[i].pos, arrival[i].pos,
                              Animation::OpLeg::K_EPS),
                  arrival[i].branch, 0.0f};
    std::vector<math::Vector> probe_normals;
    hankin_face_normals(compiled, probe, probe_normals);
    for (const math::Vector &normal : probe_normals)
      HS_EXPECT_GE(normal.magnitude(), HANKIN_FLAT_FACE);
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
    reset_globals();
    const ScopedArenaSplit split(
        IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
        IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
        IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
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
    Animation::OpLeg anim(seed,
                          Animation::OpLeg::HankinSweepSpec{
                              .theta_start = Animation::OpLeg::THETA_EPS,
                              .theta_end = site.theta_star,
                              .sweep_frames = SWEEP},
                          leg, cb, handoff);

    const Animation::OpLeg::Landing &landing = anim.landing();
    probe.ramp_count = landing.blend_pairs;
    HS_EXPECT_EQ(landing.primary_faces, seed.face_counts.size());
    HS_EXPECT_EQ(landing.faces, seed.face_counts.size() + seed.vertices.size());
    HS_EXPECT_TRUE(landing.topology != nullptr);

    hs_test::StubEffect fx(288, 144);
    for (int f = 0; f < SWEEP; ++f) {
      {
        Canvas c(fx);
        anim.step(c);
      }
      fx.advance_display();
    }
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
