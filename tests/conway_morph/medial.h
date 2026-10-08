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
  return recipe_step_seed<Solids::ICOSAHEDRON_KIS_GYRO_RECIPE, Solids::Op::DUAL,
                          0>(a, b);
}
inline PolyMesh probe_toct_snub(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE,
                          Solids::Op::DUAL, 0>(a, b);
}
inline PolyMesh probe_dodeca_hk72_ambo(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE,
                          Solids::Op::DUAL, 0>(a, b);
}
inline PolyMesh probe_icosidodeca_trunc5_ambo(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE,
                          Solids::Op::DUAL, 0>(a, b);
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
     0.0f, nullptr, false, 1e-3, .55f},
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
 *        plus slerp well-formedness.
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
    HS_EXPECT_EQ(med_b.size(), med_a.vertices.size());
    if (med_b.size() != med_a.vertices.size())
      continue;

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
 * medial's ([D-faces][D-vertex orbits] vs [P-faces][P-vertex orbits]); the
 * closing leg's from-palettes must follow the permutation, and the rendered
 * seam must stay near one in-leg step.
 */
inline void test_opleg_dual_bridge_seam_correspondence() {
  using Animation::OpLeg;
  constexpr int SWEEP = 24;
  constexpr int RW = 288, RH = 144;
  constexpr float SEAM_MATCH_TOL = 0.02f;
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
    if (!landing3.from_palette)
      continue;

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
