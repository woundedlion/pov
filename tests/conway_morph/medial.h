/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Medial (Conway-dual) bridge: MeshOps::medial yields one rectified
// connectivity with endpoint vertex sets a_e (== ambo(P)) and b_e
// (== ambo(dual(P))); the medial leg slerps a_e -> b_e.
// ---------------------------------------------------------------------------

inline PolyMesh probe_icosa_kis_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .kis()
      .snub()
      .build();
}
inline PolyMesh probe_toct_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Archimedean::truncatedOctahedron(a, b), a,
                              b)
      .snub()
      .build();
}
inline PolyMesh probe_dodeca_hk72_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(72.0f * D2R)
      .ambo()
      .build();
}
inline PolyMesh probe_icosidodeca_trunc5_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::TRUNCATE_T_NEAR;
  return Solids::SolidBuilder(Solids::Archimedean::icosidodecahedron(a, b), a,
                              b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .build();
}
/** @brief Smooth-dual bridge seeds, including the dt macro's truncate prefix. */
inline constexpr StepLegSite DUAL_LEG_SITES[] = {
    {"icosahedron_kis_gyro", probe_icosa_kis_snub, 0.0f},
    {"truncatedOctahedron_gyro", probe_toct_snub, 0.0f, nullptr, true},
    {"dodecahedron_hk72_ambo_dual", probe_dodeca_hk72_ambo, 0.0f},
    {"icosidodecahedron_truncate5d_ambo_dual", probe_icosidodeca_trunc5_ambo,
     0.0f},
    {"truncatedIcosahedron_ambo_relax100_hk54_needle",
     build_ticosa_ambo_relax100_hk54, 0.0f, nullptr, true, 5e-4},
    {"dodecahedron_bevel2_relax_gyro",
     recipe_step_seed<Solids::DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE,
                      Solids::Op::DUAL>,
     0.0f},
    {"snubDodecahedron_truncate5d_ambo_dual",
     recipe_step_seed<Solids::SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE,
                      Solids::Op::DUAL>,
     0.0f, nullptr, false, 8.5e-4},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual",
     recipe_step_seed<
         Solids::TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE,
         Solids::Op::DUAL>,
     0.0f, nullptr, false, 1e-3, .65f},
    {"truncatedIcosahedron_truncate50d_ambo_dual",
     recipe_step_seed<
         Solids::TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE,
         Solids::Op::DUAL>,
     0.0f},
};

/** @brief Max nearest-vertex distance from every vertex of @p x to @p y. */
inline float medial_vertex_set_dist(const PolyMesh &x, const PolyMesh &y) {
  float worst = 0.0f;
  for (const auto &vx : x.vertices) {
    float best = 1e9f;
    for (const auto &vy : y.vertices)
      best = std::min(best, math::distance_between(vx, vy));
    worst = fold_worst(worst, best);
  }
  return worst;
}
inline float medial_vertex_set_dist(const ArenaVector<math::Vector> &x,
                                    const PolyMesh &y) {
  float worst = 0.0f;
  for (const auto &vx : x) {
    float best = 1e9f;
    for (const auto &vy : y.vertices)
      best = std::min(best, math::distance_between(vx, vy));
    worst = fold_worst(worst, best);
  }
  return worst;
}

/** @brief Signed total solid angle of a mesh (per-face fan from its centroid). */
inline double medial_total_solid_angle(const PolyMesh &m) {
  double total = 0.0;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    math::Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    for (int k = 0; k < n; ++k) {
      const math::Vector a = m.vertices[m.faces[off + k]];
      const math::Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      const double num = math::dot(c, math::cross(a, b));
      const double den =
          1.0 + math::dot(c, a) + math::dot(a, b) + math::dot(b, c);
      total += 2.0 * std::atan2(num, den);
    }
    off += n;
  }
  return total;
}

/** @brief Smallest face area, summed from its planar triangle fan. */
inline double medial_min_face_area(const PolyMesh &m) {
  double mn = 1e9;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    math::Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    double area = 0.0;
    for (int k = 0; k < n; ++k) {
      const math::Vector a = m.vertices[m.faces[off + k]];
      const math::Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      area += 0.5 * math::cross(a - c, b - c).length();
    }
    mn = std::min(mn, area);
    off += n;
  }
  return mn;
}

/** @brief Faces whose area-weighted normal points inward (a fold/inversion). */
inline int medial_inverted_faces(const PolyMesh &m) {
  int inv = 0;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    math::Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    math::Vector nrm(0, 0, 0);
    for (int k = 0; k < n; ++k) {
      const math::Vector a = m.vertices[m.faces[off + k]];
      const math::Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      nrm = nrm + math::cross(a, b);
    }
    if (math::dot(nrm, c) < 0.0f)
      ++inv;
    off += n;
  }
  return inv;
}

/**
 * @brief Gates the medial dual bridge on every DUAL-leg seed: endpoint match
 *        plus slerp well-formedness (task-validated failure modes).
 */
inline void test_medial_dual_bridge_wellformed() {
  constexpr int SAMPLES = 33;
  // Well clear of an antipodal/coincident slerp singularity.
  constexpr float MIN_ENDPOINT_DOT = 0.9f;
  constexpr float MAX_INTER_SAMPLE_STEP = 0.05f;
  constexpr float ENDPOINT_TOL = 1e-4f;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh P = build_step_leg_seed(site, persist);

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

    PolyMesh med_a;
    ArenaVector<math::Vector> med_b;
    MeshOps::medial(P, med_a, med_b, a, b);

    // s = 0 is ambo(P), s = 1 is ambo(dual(P)): one rectified connectivity, so
    // the face count is constant. Vertex identity is checked by position: a
    // hankin seed's dual is lossy, so ambo(dual(P)) merges coincident edge
    // midpoints.
    PolyMesh ambo_p = MeshOps::ambo(P, b, aux);
    PolyMesh dual_p = MeshOps::dual(P, b, aux);
    PolyMesh ambo_dual_p = MeshOps::ambo(dual_p, aux, b);
    HS_EXPECT_EQ(med_a.vertices.size(), ambo_p.vertices.size());
    HS_EXPECT_EQ(med_a.face_counts.size(), ambo_p.face_counts.size());
    HS_EXPECT_EQ(med_a.face_counts.size(), ambo_dual_p.face_counts.size());
    HS_EXPECT_TRUE(
        std::equal(med_a.face_counts.begin(), med_a.face_counts.end(),
                   ambo_p.face_counts.begin(), ambo_p.face_counts.end()));
    HS_EXPECT_TRUE(std::equal(med_a.faces.begin(), med_a.faces.end(),
                              ambo_p.faces.begin(), ambo_p.faces.end()));

    // out_a is bit-identical to ambo(P); pin it before snorm16 packing can
    // absorb a drift.
    const size_t pair_n =
        std::min(med_a.vertices.size(), ambo_p.vertices.size());
    size_t bit_diff = 0;
    float worst_bit = 0.0f;
    for (size_t v = 0; v < pair_n; ++v) {
      const math::Vector &mv = med_a.vertices[v];
      const math::Vector &av = ambo_p.vertices[v];
      if (std::bit_cast<uint32_t>(mv.x) != std::bit_cast<uint32_t>(av.x) ||
          std::bit_cast<uint32_t>(mv.y) != std::bit_cast<uint32_t>(av.y) ||
          std::bit_cast<uint32_t>(mv.z) != std::bit_cast<uint32_t>(av.z)) {
        ++bit_diff;
        worst_bit = fold_worst(worst_bit, math::distance_between(mv, av));
      }
    }
    HS_EXPECT_EQ(bit_diff, size_t(0));
    if (bit_diff != 0)
      std::printf("    [medial] %s out_a != ambo(P) exactly: %zu/%zu vertices, "
                  "worst |d|=%.3e\n",
                  site.name, bit_diff, pair_n, (double)worst_bit);

    // Gate the snorm16-decoded positions the MEDIAL_SLERP leg slerps.
    for (auto &v : med_a.vertices)
      v = math::Snorm3::encode(v).decode().normalized();
    for (auto &v : med_b)
      v = math::Snorm3::encode(v).decode().normalized();

    // Every medial vertex sits on an ambo(P) vertex at s=0 and an ambo(dual(P))
    // vertex at s=1 (set containment; a lossy dual makes s=1 many-to-one).
    const float end_a = medial_vertex_set_dist(med_a, ambo_p);
    const float end_b = medial_vertex_set_dist(med_b, ambo_dual_p);
    HS_EXPECT_LT(end_a, ENDPOINT_TOL);
    HS_EXPECT_LT(end_b, ENDPOINT_TOL);

    // No antipodal/coincident slerp inputs.
    float min_dot = 2.0f;
    for (size_t v = 0; v < med_a.vertices.size(); ++v)
      min_dot = std::min(min_dot, math::dot(med_a.vertices[v].normalized(),
                                            med_b[v].normalized()));
    HS_EXPECT_GT(min_dot, MIN_ENDPOINT_DOT);

    // Slerp sweep: fixed connectivity, per-vertex slerp; well-formed throughout.
    int total_inv = 0;
    double worst_4pi = 0.0, min_area = 1e9;
    float max_step = 0.0f;
    SweepFingerprint first;
    std::vector<math::Vector> prev(med_a.vertices.size());
    for (int s = 0; s < SAMPLES; ++s) {
      const float k = static_cast<float>(s) / (SAMPLES - 1);
      ScratchScope frame_a(aux);
      PolyMesh frame;
      frame.vertices.bind(aux, med_a.vertices.size());
      for (size_t v = 0; v < med_a.vertices.size(); ++v)
        frame.vertices.push_back(math::slerp(med_a.vertices[v], med_b[v], k));
      frame.face_counts.bind(aux, med_a.face_counts.size());
      frame.face_counts.append_bulk(med_a.face_counts.data(),
                                    med_a.face_counts.size());
      frame.faces.bind(aux, med_a.faces.size());
      frame.faces.append_bulk(med_a.faces.data(), med_a.faces.size());

      {
        ScratchScope ca(a);
        ScratchScope cb(b);
        const SweepFingerprint fp = check_sweep_sample(frame, a, b);
        if (s == 0)
          first = fp;
        else
          expect_same_fingerprint(fp, first); // fixed emission order/count
      }
      total_inv += medial_inverted_faces(frame);
      worst_4pi =
          fold_worst(worst_4pi, std::abs(medial_total_solid_angle(frame) -
                                         4.0 * 3.14159265358979323846));
      min_area = std::min(min_area, medial_min_face_area(frame));
      if (s > 0)
        for (size_t v = 0; v < frame.vertices.size(); ++v)
          max_step = std::max(
              max_step, math::distance_between(frame.vertices[v], prev[v]));
      prev.assign(frame.vertices.begin(), frame.vertices.end());
    }

    HS_EXPECT_EQ(total_inv, 0);
    HS_EXPECT_LT(worst_4pi, 1e-5);
    HS_EXPECT_GT(min_area, site.min_face_area);
    HS_EXPECT_LT(max_step, MAX_INTER_SAMPLE_STEP);

    if (hs_test::stats().failed != failed_before)
      std::printf("    [medial] %s FAILED (endA=%.2e endB=%.2e min a.b=%.4f "
                  "inv=%d 4pi-err=%.2e minA=%.3e step=%.4f)\n",
                  site.name, (double)end_a, (double)end_b, (double)min_dot,
                  total_inv, worst_4pi, min_area, (double)max_step);
    else
      std::printf(
          "  [medial] %s: F=%zu endpoints<%.0e min a.b=%.4f 4pi-err=%.1e "
          "minA=%.2e step=%.4f\n",
          site.name, med_a.face_counts.size(), (double)ENDPOINT_TOL,
          (double)min_dot, worst_4pi, min_area, (double)max_step);
  }
}

/**
 * @brief Drives a MEDIAL_SLERP leg to completion through OpLeg on every DUAL-leg
 *        seed: constant compiled face count, bounded per-frame motion, a total
 *        palette mapping, and a landing describing the whole medial face list.
 * @details The leg departs from ambo(P) (one face per primal face + one per
 * primal vertex), so the handoff is class-keyed on that mesh with its face
 * centroids for the geometric provenance mapping.
 */
inline void test_opleg_medial_leg_smoke() {
  using Animation::OpLeg;
  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  constexpr int SWEEP = 24;
  constexpr float MAX_STEP_CHORD = 0.15f;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh P = build_step_leg_seed(site, persist);

    Arena leg(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

    // Departed mesh: ambo(P), class-keyed palettes plus face centroids.
    PolyMesh ambo_p = Solids::finalize_solid(MeshOps::ambo(P, leg, temp), leg);
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(ambo_p, scratch_arena_a,
                                          scratch_arena_b, leg);
    }
    const size_t prev_faces = ambo_p.face_counts.size();
    std::vector<uint8_t> pal(prev_faces);
    std::vector<math::Vector> centroid(prev_faces);
    size_t off = 0;
    for (size_t f = 0; f < prev_faces; ++f) {
      pal[f] = static_cast<uint8_t>(
          math::wrap(static_cast<int>(ambo_p.topology[f]), OpLeg::PALETTES));
      const int n = ambo_p.face_counts[f];
      centroid[f] = face_centroid_unit(ambo_p, off, n);
      off += n;
    }

    OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                  .prev_face_palette = pal.data(),
                                  .prev_faces = prev_faces,
                                  .prev_face_centroid = centroid.data()};

    LegDrawProbe probe;
    auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
      probe.observe(m, sh);
    };

    OpLeg anim(P, OpLeg::MedialSpec{.sweep_frames = SWEEP}, leg, cb, handoff);
    const OpLeg::Landing &landing = anim.landing();
    probe.ramp_count = landing.blend_pairs;
    // The medial connectivity is ambo(P), so the leg's face list is the whole
    // rectified polyhedron and every face lives the whole slerp.
    HS_EXPECT_EQ(landing.faces, prev_faces);
    HS_EXPECT_EQ(landing.primary_faces, prev_faces);

    hs_test::StubEffect fx(288, 144);
    for (int f = 0; f < SWEEP; ++f) {
      {
        Canvas c(fx);
        anim.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(probe.drawn, (size_t)SWEEP);
    HS_EXPECT_EQ(landing.faces, probe.faces);
    HS_EXPECT_LT(probe.worst_step, MAX_STEP_CHORD);
    HS_EXPECT_TRUE(landing.from_palette != nullptr);

    std::printf("  [opleg medial] %s: F=%zu across %d frames, worst step "
                "%.4f%s\n",
                site.name, landing.faces, SWEEP, (double)probe.worst_step,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/**
 * @brief Drives the dual bridge's medial -> closing-leg handoff on every
 *        DUAL-leg seed and pins the face correspondence across the seam.
 * @details The closing leg's face list is block-transposed against the
 * medial's ([D-faces][D-vertex orbits] vs [P-faces][P-vertex orbits]). The
 * permutation is derived by exact centroid matching at the ambo point, and the
 * closing leg's from-palettes must follow it. A rendered A/B (leg-2 last frame
 * vs leg-3 first frame) must stay near one in-leg step. The needle site
 * (truncate(X, 1/3), unequal blocks) runs on the bridge arena split.
 */
inline void test_opleg_dual_bridge_seam_correspondence() {
  using Animation::OpLeg;
  constexpr int SWEEP = 24;
  constexpr int RW = 288, RH = 144;
  constexpr float SEAM_MATCH_TOL = 0.02f;
  // Largest measured frame-wide sum of absolute channel deltas across sites: 917321.
  constexpr long long SEAM_SUMABS_MAX = 1100000ll;

  static Pipeline<RW, RH> filters;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    reset_globals();
    const ScopedArenaSplit split(
        IslamicStars<288, 144>::BRIDGE_BUDGET.persistent(GLOBAL_ARENA_SIZE),
        IslamicStars<288, 144>::BRIDGE_BUDGET.scratch_a,
        IslamicStars<288, 144>::BRIDGE_BUDGET.scratch_b);
    hs::random().seed(2026u);

    Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
    MeshPaletteBank bank;
    bank.bake_all(bank_arena);

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    Arena leg(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

    PolyMesh P = build_step_leg_seed(site, persist);
    const size_t PF = P.face_counts.size();

    // Departed mesh ambo(P): class-keyed palettes in medial face order.
    PolyMesh ambo_p;
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      ambo_p = Solids::finalize_solid(MeshOps::ambo(P, aux, temp), leg);
    }
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(ambo_p, scratch_arena_a,
                                          scratch_arena_b, leg);
    }
    const size_t nf = ambo_p.face_counts.size();
    std::vector<uint8_t> pal2(nf);
    for (size_t f = 0; f < nf; ++f)
      pal2[f] = static_cast<uint8_t>(
          math::wrap(static_cast<int>(ambo_p.topology[f]), OpLeg::PALETTES));
    hs_test::StubEffect fx(RW, RH);
    std::vector<Pixel> snaps[3]; // leg-2 last, leg-3 first, leg-3 second
    int drawn = 0, rasterize_at = -1;
    auto cb = [&](Canvas &c, const MeshState &m, const OpLeg::Shading &sh) {
      ++drawn;
      if (drawn != rasterize_at)
        return;
      auto shader = [&](const math::Vector &, Fragment &frag) {
        int fi = static_cast<int>(frag.v2);
        int ramp =
            (fi >= 0 && fi < static_cast<int>(sh.faces)) ? sh.face_ramp[fi] : 0;
        float t = hs::clamp(fragment_edge_dist(frag), 0.0f, 1.0f);
        frag.color = sh.ramps[ramp].get(t);
        frag.color.alpha = 1.0f;
      };
      Scan::Mesh::draw<RW, RH>(filters, c, m, shader, scratch_arena_b);
    };
    auto snap = [&](std::vector<Pixel> &out) {
      capture_frame<RW, RH>(fx, out);
    };

    // Leg 2: the medial slerp, departed from ambo(P), with crossfading
    // handoffs; the seeded RNG keeps the target shuffles deterministic.
    OpLeg::PaletteHandoff handoff2{.bank = &bank.bank,
                                   .prev_face_palette = pal2.data(),
                                   .prev_faces = nf,
                                   .correspondence =
                                       OpLeg::FaceCorrespondence::IDENTITY};
    OpLeg leg2(P, OpLeg::MedialSpec{.sweep_frames = SWEEP}, leg, cb, handoff2);
    const OpLeg::Landing &landing2 = leg2.landing();
    HS_EXPECT_EQ(landing2.faces, nf);
    HS_EXPECT_TRUE(landing2.arrival_topology != nullptr);
    HS_EXPECT_TRUE(landing2.arrival_point != nullptr);
    HS_EXPECT_TRUE(landing2.from_palette != nullptr);
    if (landing2.faces != nf || !landing2.arrival_topology ||
        !landing2.arrival_point || !landing2.from_palette)
      continue;
    HS_EXPECT_EQ(landing2.arrival_topology->face_counts.size(), nf);
    if (landing2.arrival_topology->face_counts.size() != nf)
      continue;
    for (size_t f = 0; f < nf; ++f)
      HS_EXPECT_EQ(landing2.from_palette[f], pal2[f]);
    rasterize_at = SWEEP;
    for (int f = 0; f < SWEEP; ++f) {
      {
        Canvas c(fx);
        leg2.step(c);
      }
      fx.advance_display();
    }
    snap(snaps[0]);

    // Leg-3 handoff, exactly as schedule_dual_untruncate builds it: departed
    // centroids and palettes cached by leg 2.
    std::vector<uint8_t> pal3(nf);
    std::vector<math::Vector> cen3(nf);
    const PolyMesh &arrival = *landing2.arrival_topology;
    size_t off = 0;
    for (size_t f = 0; f < nf; ++f) {
      math::Vector c(0.0f, 0.0f, 0.0f);
      for (int j = 0; j < arrival.face_counts[f]; ++j) {
        const size_t v = arrival.faces[off + j];
        c = c + landing2.arrival_point[v].decode().normalized();
      }
      cen3[f] = c.normalized();
      pal3[f] = landing2.landed_palette(f);
      off += arrival.face_counts[f];
    }

    PolyMesh D;
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      D = Solids::finalize_solid(MeshOps::dual(P, aux, temp), leg);
    }
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(D, scratch_arena_a, scratch_arena_b,
                                          leg);
    }
    // dual drops sub-3 orbits, so D's faces are exactly the medial's P-vertex
    // block: the two blocks partition the shared face count.
    const size_t DF = D.face_counts.size();
    HS_EXPECT_EQ(DF, nf - PF);

    OpLeg::PaletteHandoff handoff3{.bank = &bank.bank,
                                   .prev_face_palette = pal3.data(),
                                   .prev_faces = nf,
                                   .correspondence =
                                       OpLeg::FaceCorrespondence::DUAL_CLOSING};
    OpLeg::BookendClasses bookend3{.topology = D.topology.data(), .faces = DF};
    OpLeg leg3(OpLeg::SweepSeed::borrow(D),
               OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                     .t_start = 0.5f,
                                     .t_end = 0.0f,
                                     .sweep_frames = SWEEP,
                                     .bridge_provenance = true},
               leg, cb, handoff3, bookend3);
    const OpLeg::Landing &landing3 = leg3.landing();
    HS_EXPECT_EQ(landing3.faces, nf);
    HS_EXPECT_TRUE(landing3.from_palette != nullptr);

    // Derive the true seam permutation from the exact geometry: ambo(D) has
    // the closing leg's face order and the seam's exact positions.
    std::vector<int> perm(nf, -1);
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh ambo_d = MeshOps::ambo(D, aux, temp);
      HS_EXPECT_EQ(ambo_d.face_counts.size(), nf);
      if (ambo_d.face_counts.size() != nf)
        continue;
      std::vector<char> used(nf, 0);
      int bad_match = 0, non_bijective = 0;
      size_t off = 0;
      for (size_t l = 0; l < nf; ++l) {
        math::Vector c(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < ambo_d.face_counts[l]; ++j)
          c = c + ambo_d.vertices[ambo_d.faces[off + j]];
        c = c.normalized();
        off += ambo_d.face_counts[l];
        size_t best = 0;
        float bd = 1e9f;
        for (size_t m = 0; m < nf; ++m) {
          const float d = math::distance_between(c, cen3[m]);
          if (d < bd) {
            bd = d;
            best = m;
          }
        }
        if (bd > SEAM_MATCH_TOL)
          ++bad_match;
        if (used[best])
          ++non_bijective;
        used[best] = 1;
        perm[l] = static_cast<int>(best);
      }
      HS_EXPECT_EQ(bad_match, 0);
      HS_EXPECT_EQ(non_bijective, 0);
    }

    // Block structure: D-faces are the medial's P-vertex orbits in emission
    // order; D-vertex orbit faces land in the medial's P-face block.
    int block1_viol = 0, block2_viol = 0, from_mismatch = 0;
    for (size_t l = 0; l < nf; ++l) {
      if (l < DF && perm[l] != static_cast<int>(PF + l))
        ++block1_viol;
      if (l >= DF && perm[l] >= static_cast<int>(PF))
        ++block2_viol;
      if (landing3.from_palette[l] != pal3[perm[l]])
        ++from_mismatch;
    }
    HS_EXPECT_EQ(block1_viol, 0);
    HS_EXPECT_EQ(block2_viol, 0);
    HS_EXPECT_EQ(from_mismatch, 0);

    // Rendered seam A/B plus a one-step in-leg control.
    drawn = 0;
    rasterize_at = 1;
    {
      Canvas c(fx);
      leg3.step(c);
    }
    fx.advance_display();
    snap(snaps[1]);
    rasterize_at = 2;
    {
      Canvas c(fx);
      leg3.step(c);
    }
    fx.advance_display();
    snap(snaps[2]);
    HS_EXPECT_EQ(drawn, 2);
    for (const auto &frame : snaps) {
      size_t lit = 0;
      for (const Pixel &pixel : frame)
        lit += !is_black(pixel);
      HS_EXPECT_GT(lit, static_cast<size_t>(RW * RH / 100));
    }

    auto diff = [&](const std::vector<Pixel> &a, const std::vector<Pixel> &b,
                    long long &sumabs, int &changed) {
      sumabs = 0;
      changed = 0;
      for (size_t i = 0; i < a.size(); ++i) {
        const int dr = std::abs((a[i].r >> 8) - (b[i].r >> 8));
        const int dg = std::abs((a[i].g >> 8) - (b[i].g >> 8));
        const int db = std::abs((a[i].b >> 8) - (b[i].b >> 8));
        sumabs += dr + dg + db;
        if (dr | dg | db)
          ++changed;
      }
    };
    long long seam_sum = 0, ctrl_sum = 0;
    int seam_px = 0, ctrl_px = 0;
    diff(snaps[0], snaps[1], seam_sum, seam_px);
    diff(snaps[1], snaps[2], ctrl_sum, ctrl_px);

    const int seam_quiet = RW * RH - seam_px;
    const int ctrl_quiet = RW * RH - ctrl_px;
    HS_EXPECT_GT(ctrl_quiet, 0);
    if (ctrl_quiet == 0)
      continue;
    const float quiet_ratio =
        static_cast<float>(seam_quiet) / static_cast<float>(ctrl_quiet);
    HS_EXPECT_GT(quiet_ratio, site.seam_quiet_ratio);
    HS_EXPECT_LT(seam_sum, SEAM_SUMABS_MAX);

    std::printf("  [opleg seam] %s: F=%zu blocks %zu/%zu, seam diff "
                "sumabs=%lld px=%d (control %lld/%d) quiet=%.3f/%.3f%s\n",
                site.name, nf, DF, nf - DF, seam_sum, seam_px, ctrl_sum,
                ctrl_px, static_cast<double>(quiet_ratio),
                static_cast<double>(site.seam_quiet_ratio),
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/** @brief Recipe-step leg kind driven by the smoke test. */
enum class StepLegKind { TRUNCATE, SNUB, RELAX };

/** Paused redraws check_step_leg_smoke issues at its pause frame. */
inline constexpr int PAUSED_REDRAWS = 2;

/**
 * @brief Drives one recipe-step leg to completion through OpLeg, gating
 *        compiled-face-count constancy, ramp indices and per-frame motion.
 * @param kind Leg kind under test.
 * @param site Seed site and arrival parameter.
 * @param frames Leg length in frames.
 * @param max_step_chord Per-frame max-vertex-motion bound.
 * @param easing Easing driving the sweep parameter.
 * @param frames_out Optional sink receiving each drawn frame's vertices.
 * @param pause_after Leg frame after which two step_paused() redraws run; 0
 *        drives the leg unpaused.
 * @details The departed palettes come from the seed's own classification and
 * the bookend from the step's clean endpoint.
 */
inline void check_step_leg_smoke(
    StepLegKind kind, const StepLegSite &site, int frames, float max_step_chord,
    EasingFn easing = math::ease_in_out_sin,
    std::vector<std::vector<math::Vector>> *frames_out = nullptr,
    int pause_after = 0) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;

  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    PolyMesh raw;
    switch (kind) {
    case StepLegKind::TRUNCATE:
      raw = MeshOps::truncate(seed, ea, eb, site.param);
      break;
    case StepLegKind::SNUB:
      raw = MeshOps::snub(seed, ea, eb, site.param, 0.0f);
      break;
    case StepLegKind::RELAX:
      raw = MeshOps::relax_baked(seed, ea, *site.bake);
      break;
    }
    endpoint = Solids::finalize_solid(raw, leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  math::Vector *prev_centroid = leg_arena.allocate_n<math::Vector>(prev_faces);
  size_t off = 0;
  for (size_t f = 0; f < prev_faces; ++f) {
    prev_pal[f] = static_cast<uint8_t>(
        slots[math::wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);
    const int n = seed.face_counts[f];
    prev_centroid[f] = face_centroid_unit(seed, off, n);
    off += n;
  }

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces,
                                .prev_face_centroid = prev_centroid};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  LegDrawProbe probe;
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    probe.observe(m, sh);
    if (frames_out)
      frames_out->push_back(probe.prev_v);
  };

  hs_test::StubEffect fx(288, 144);

  auto run = [&](OpLeg &&leg) {
    const OpLeg::Landing &landing = leg.landing();
    probe.ramp_count = landing.blend_pairs;
    HS_EXPECT_EQ(landing.primary_faces, seed.face_counts.size());
    HS_EXPECT_TRUE(landing.topology != nullptr);
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
      if (f + 1 != pause_after)
        continue;
      for (int p = 0; p < PAUSED_REDRAWS; ++p) {
        {
          Canvas c(fx);
          leg.step_paused(c);
        }
        fx.advance_display();
      }
    }
    HS_EXPECT_EQ(probe.drawn,
                 (size_t)(frames + (pause_after > 0 ? PAUSED_REDRAWS : 0)));
    HS_EXPECT_EQ(probe.faces, landing.faces);
  };

  switch (kind) {
  case StepLegKind::TRUNCATE:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::SNUB:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::RELAX:
    run(OpLeg(seed, OpLeg::RelaxSpec{.bake = site.bake, .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  }

  HS_EXPECT_LT(probe.worst_step, max_step_chord);
  const char *label = kind == StepLegKind::TRUNCATE ? "truncate"
                      : kind == StepLegKind::SNUB   ? "snub"
                                                    : "relax";
  std::printf("  [opleg %s] %s: %zu faces, worst per-frame vertex step %.4f "
              "chord (bound %.2f)%s\n",
              label, site.name, probe.faces, (double)probe.worst_step,
              (double)max_step_chord,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/** Seed fixtures for kis and dual gated-swap legs. */
inline constexpr StepLegSite GATE_LEG_SITES[] = {
    {"icosahedron", probe_icosahedron, 0.0f},
    {"icosahedron_snub", probe_icosa_snub, 0.0f},
};

/**
 * @brief Drives one gated-swap leg to completion through OpLeg, gating the
 *        per-side compiled-face-count constancy, the single swap, the ramp
 *        indices.
 * @param op Partition operator the leg swaps to.
 * @param site Seed site.
 * @param gate Half-gate length in frames.
 * @details Derives this test-local handoff and bookend: the
 * departed palettes from the seed's own classification, the bookend from the
 * op's clean endpoint.
 */
inline void check_gated_leg_smoke(Animation::OpLeg::SwapOp op,
                                  const StepLegSite &site, int gate) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;
  const bool is_kis = op == OpLeg::SwapOp::KIS;

  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    endpoint = Solids::finalize_solid(is_kis ? MeshOps::kis(seed, ea, eb)
                                             : MeshOps::dual(seed, ea, eb),
                                      leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  for (size_t f = 0; f < prev_faces; ++f)
    prev_pal[f] = static_cast<uint8_t>(
        slots[math::wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  const int frames = 2 * gate + 1;
  int drawn = 0, swaps = 0, swap_frame = -1;
  size_t side_faces = 0;
  int target_ramps = 0;
  std::array<bool, OpLeg::PALETTES> seed_palettes{};
  for (size_t f = 0; f < prev_faces; ++f)
    seed_palettes[prev_pal[f]] = true;
  const int SEED_RAMPS = static_cast<int>(
      std::count(seed_palettes.begin(), seed_palettes.end(), true));
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                   drawn < gate ? SEED_RAMPS : target_ramps);

    if (drawn == 0) {
      side_faces = sh.faces;
    } else if (sh.faces != side_faces) {
      ++swaps;
      swap_frame = drawn;
      side_faces = sh.faces;
    }
    ++drawn;
  };

  hs_test::StubEffect fx(288, 144);

  OpLeg leg(seed, OpLeg::GatedSwapSpec{.op = op, .gate_frames = gate},
            leg_arena, cb, handoff, bookend);
  const OpLeg::Landing &landing = leg.landing();
  target_ramps = landing.blend_pairs;
  for (int f = 0; f < frames; ++f) {
    {
      Canvas c(fx);
      leg.step(c);
    }
    fx.advance_display();
  }

  HS_EXPECT_EQ(drawn, frames);
  // Topology changes exactly once, at the swap frame, and the leg lands on the
  // arrival face count.
  HS_EXPECT_EQ(swaps, 1);
  HS_EXPECT_EQ(swap_frame, gate);
  HS_EXPECT_EQ(side_faces, landing.faces);
  HS_EXPECT_EQ(landing.faces, endpoint.face_counts.size());
  HS_EXPECT_EQ(landing.primary_faces, prev_faces);
  HS_EXPECT_TRUE(landing.from_palette != nullptr);

  std::printf("  [opleg %s] %s: F %zu -> %zu, swap at frame %d of %d%s\n",
              is_kis ? "kis" : "dual", site.name, prev_faces, landing.faces,
              swap_frame, frames,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/**
 * @brief Smoke-tests a gated-swap leg for both partition ops on every gate
 *        site, with a test-local 6-frame half-gate.
 */
inline void test_opleg_gated_swap_smoke() {
  using Animation::OpLeg;
  for (const StepLegSite &site : GATE_LEG_SITES) {
    check_gated_leg_smoke(OpLeg::SwapOp::KIS, site, 6);
    check_gated_leg_smoke(OpLeg::SwapOp::DUAL, site, 6);
  }
}

/**
 * @brief Smoke-tests one leg of each recipe-step kind end to end.
 */
inline void test_opleg_step_leg_smoke() {
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES, 0.15f);
  check_step_leg_smoke(StepLegKind::SNUB, SNUB_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES, 0.15f);
  for (const StepLegSite &site : RELAX_LEG_SITES)
    check_step_leg_smoke(StepLegKind::RELAX, site,
                         RecipeLegLengths::RELAX_LEG_FRAMES, 0.15f);
}

/**
 * @brief Pins a paused leg's redraw: the held frame, drawn again, with the
 *        sweep clock untouched.
 * @details Both paused draws must reproduce the held frame's vertices bitwise,
 * and the resumed frames must match the unpaused run frame for frame.
 */
inline void test_opleg_step_paused_holds_frame() {
  constexpr int FRAMES = 8;
  constexpr int PAUSE_AFTER = 4;
  // Loose chord bound: motion is not under test.
  constexpr float CHORD_MAX = 1.0f;
  std::vector<std::vector<math::Vector>> unpaused, held;
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &unpaused);
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &held, PAUSE_AFTER);
  HS_EXPECT_SIZE_OR_RETURN(unpaused, (size_t)FRAMES);
  HS_EXPECT_EQ(held.size(), (size_t)(FRAMES + PAUSED_REDRAWS));
  if (held.size() != (size_t)(FRAMES + PAUSED_REDRAWS))
    return;

  auto identical = [](const std::vector<math::Vector> &a,
                      const std::vector<math::Vector> &b) {
    if (a.size() != b.size() || a.empty())
      return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z)
        return false;
    return true;
  };

  // Every paused draw repeats the frame the leg was holding.
  for (int p = 0; p < PAUSED_REDRAWS; ++p)
    HS_EXPECT_TRUE(identical(held[PAUSE_AFTER - 1], held[PAUSE_AFTER + p]));

  // The clock never moved: the resumed frames are the unpaused run's.
  size_t drifted = 0;
  for (int f = 0; f < FRAMES; ++f) {
    const size_t k = f < PAUSE_AFTER ? (size_t)f : (size_t)(f + PAUSED_REDRAWS);
    if (!identical(unpaused[f], held[k]))
      ++drifted;
  }
  HS_EXPECT_EQ(drifted, (size_t)0);
  std::printf("  [opleg paused] truncate: %d frames, %d paused redraws at "
              "frame %d, %zu resumed frames drifted\n",
              FRAMES, PAUSED_REDRAWS, PAUSE_AFTER, drifted);
}

/**
 * @brief Drives swept legs under an easing whose range leaves [0, 1], gating
 *        the sweep parameter against extrapolation past the arrival.
 * @details ease_out_elastic exceeds 1 over x in [0.075, 0.225], peaking near
 * 1.35. Clamped, frames 2-5 of 24 hold the arrival exactly (the
 * frame-4-vs-last comparison); frame 1 is still mid-sweep. The chord bound is
 * loose because elastic covers half the sweep in one frame.
 */
inline void test_opleg_step_leg_overshooting_easing() {
  constexpr StepLegSite NEAR_AMBO{"icosahedron_ambo_truncate049",
                                  probe_icosa_ambo, 0.49f};
  constexpr int FRAMES = 24;
  for (int k = 0; k < 2; ++k) {
    std::vector<std::vector<math::Vector>> drawn;
    const bool truncate = k == 0;
    check_step_leg_smoke(truncate ? StepLegKind::TRUNCATE : StepLegKind::SNUB,
                         truncate ? NEAR_AMBO : SNUB_LEG_SITES[0], FRAMES, 2.0f,
                         math::ease_out_elastic, &drawn);
    HS_EXPECT_SIZE_OR_RETURN(drawn, FRAMES);
    const std::vector<math::Vector> &arrival = drawn.back();
    const std::vector<math::Vector> &peak = drawn[3];
    const std::vector<math::Vector> &opening = drawn[0];
    HS_EXPECT_SIZE_OR_RETURN(peak, arrival.size());
    HS_EXPECT_SIZE_OR_RETURN(opening, arrival.size());
    float worst_peak = 0.0f, worst_opening = 0.0f;
    for (size_t i = 0; i < arrival.size(); ++i) {
      worst_peak =
          fold_worst(worst_peak, math::distance_between(peak[i], arrival[i]));
      worst_opening = fold_worst(
          worst_opening, math::distance_between(opening[i], arrival[i]));
    }
    HS_EXPECT_LT(worst_peak, 1e-5f);
    HS_EXPECT_GT(worst_opening, 1e-3f);
    std::printf("  [opleg elastic] %s: overshoot frame %.2e from arrival, "
                "opening frame %.4f\n",
                truncate ? "truncate" : "snub", (double)worst_peak,
                (double)worst_opening);
  }
}

/**
 * @brief Pins the crossfade colour model on a ConwayGraph edge leg under the
 *        classic_blend default: exact `from` at frame 1, moved off `from` by
 *        mid-leg, exact `to` at the arrival frame.
 */
inline void test_opleg_edge_leg_crossfade() {
  using Animation::OpLeg;
  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  hs_test::StubEffect fx(288, 144);

  // LUT grid-aligned sample coordinates for exact ramp-color comparisons.
  constexpr float SAMPLES[] = {0.0f, 0.5f, 1.0f};
  const OpLeg::Landing *lp = nullptr;
  std::vector<char> all_from, all_to; // per drawn frame: exact from/to palettes
  auto matches = [&](const OpLeg::Shading &sh, size_t f, uint8_t pal) {
    for (float t : SAMPLES) {
      const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
      const Color4 exp = bank.bank.entries[pal].get(t);
      if (got.color.r != exp.color.r || got.color.g != exp.color.g ||
          got.color.b != exp.color.b)
        return false;
    }
    return true;
  };
  auto cb = [&](Canvas &, const MeshState &, const OpLeg::Shading &sh) {
    bool from_ok = true, to_ok = true;
    for (size_t f = 0; f < sh.faces; ++f) {
      const uint8_t from = lp->from_palette[f];
      const uint8_t to = lp->to_palette[math::wrap(
          static_cast<int>(lp->topology[f]), OpLeg::PALETTES)];
      from_ok = from_ok && matches(sh, f, from);
      to_ok = to_ok && matches(sh, f, to);
    }
    all_from.push_back(from_ok);
    all_to.push_back(to_ok);
  };
  auto divergent_faces = [&]() {
    int n = 0;
    for (size_t f = 0; f < lp->faces; ++f)
      if (lp->from_palette[f] !=
          lp->to_palette[math::wrap(static_cast<int>(lp->topology[f]),
                                    OpLeg::PALETTES)])
        ++n;
    return n;
  };
  auto run_frames = [&](OpLeg &leg, int frames) {
    all_from.clear();
    all_to.clear();
    lp = &leg.landing();
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(all_from.size(), (size_t)frames);
  };
  {
    const int failed_before = hs_test::stats().failed;
    Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, leg_arena);
    uint8_t pal[16];
    for (size_t f = 0; f < cube.face_counts.size(); ++f)
      pal[f] = static_cast<uint8_t>(f % OpLeg::PALETTES);
    const int edge = find_directed_edge(ConwayGraph::EDGES, ConwayGraph::CUBE,
                                        ConwayGraph::TRUNCATED_CUBE);
    HS_EXPECT_GE(edge, 0);
    if (edge < 0)
      return;
    OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                  .prev_face_palette = pal,
                                  .prev_faces = cube.face_counts.size()};
    constexpr int EDGE_FRAMES = 24;
    OpLeg leg(cube,
              OpLeg::EdgeSweepSpec{.edge = &ConwayGraph::EDGES[edge],
                                   .sweep_frames = EDGE_FRAMES},
              leg_arena, cb, handoff);
    run_frames(leg, EDGE_FRAMES);
    HS_EXPECT_SIZE_OR_RETURN(all_from, EDGE_FRAMES);
    HS_EXPECT_SIZE_OR_RETURN(all_to, EDGE_FRAMES);
    HS_EXPECT_GT(divergent_faces(), 0);
    HS_EXPECT_TRUE(all_from[0]);
    HS_EXPECT_TRUE(!all_from[EDGE_FRAMES / 2 - 1]);
    HS_EXPECT_TRUE(all_to[EDGE_FRAMES - 1]);
    std::printf("  [crossfade] edge leg: F=%zu mid-leg blend intact%s\n",
                lp->faces,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/**
 * @brief Pins recipe-step acceptance at the supported sweep boundaries.
 */
inline void test_unsweepable_recipe_steps_are_gated() {
  using Solids::Op;
  using Solids::IslamicStarPatterns::D2R;
  using Solids::IslamicStarPatterns::TRUNCATE_T_FAR;

  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.33f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::SNUB, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 8.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::HANKIN, 62.0f * D2R}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::AMBO}));
  // 0.01 is below T_EPS but sweepable: the leg births at the derived
  // per-arrival floor.
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.01f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::TRUNCATE, 0.001f}));
  // At t == 1 the cut faces collapse.
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, TRUNCATE_T_FAR}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::TRUNCATE, 1.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::CHAMFER, 0.63f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::KIS}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::DUAL}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::CHAMFER, 0.001f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::CHAMFER, 0.9f}));
  // Zero snub inset, zero hankin angle and bake-less relax below one iteration
  // are not morphable.
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::SNUB, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::HANKIN, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::RELAX, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::RELAX, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 1.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step(
      {.op = Op::RELAX,
       .bake = &Solids::RelaxBakes::dodecahedron_bevel20_converged}));
}
