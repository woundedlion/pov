/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// Ambo-on-hankin probe (docs/specs/opchain_morph_spec.md, "Validation
// contract"): HANKIN_AMBO_SITES samples ambo legs on hankin meshes as a
// truncate sweep whose compiled face count must not move within the leg and
// whose sampled meshes retain two-face edge incidence and Euler characteristic 2.
// ---------------------------------------------------------------------------

/** @brief One ambo-on-hankin sweep seed from the Islamic registry chains. */
struct HankinAmboSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the ambo. */
};

inline PolyMesh probe_dodeca_hk62(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(62.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk72(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(72.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk35(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk35_ambo_hk62(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .build();
}
inline PolyMesh probe_octa_hk17(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::octahedron(a, b), a, b)
      .hankin(17.0f * D2R)
      .build();
}
inline PolyMesh probe_octa_hk34(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::octahedron(a, b), a, b)
      .hankin(34.0f * D2R)
      .build();
}
inline PolyMesh probe_rhombicubocta_hk63(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::rhombicuboctahedron(a, b), a,
                              b)
      .hankin(63.0f * D2R)
      .build();
}
inline PolyMesh probe_ticosa_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .hankin(54.0f * D2R)
      .build();
}

inline constexpr HankinAmboSite HANKIN_AMBO_SITES[] = {
    {"dodecahedron_hk62", probe_dodeca_hk62},
    {"dodecahedron_hk72", probe_dodeca_hk72},
    {"dodecahedron_hk35", probe_dodeca_hk35},
    {"dodecahedron_hk35_ambo_hk62", probe_dodeca_hk35_ambo_hk62},
    {"dodecahedron_hk54", probe_dodeca_hk54},
    {"octahedron_hk17", probe_octa_hk17},
    {"octahedron_hk34", probe_octa_hk34},
    {"rhombicuboctahedron_hk63", probe_rhombicubocta_hk63},
    {"truncatedIcosahedron_hk54", probe_ticosa_hk54},
};

/**
 * @brief Steps a truncate (ambo-leg) sweep on HANKIN_AMBO_SITES, asserting constant raw and compiled face counts, two-face
 *        edge incidence, and Euler characteristic 2 at each sampled parameter.
 */
inline void test_ambo_leg_on_hankin_seed_holds_topology() {
  constexpr int SAMPLES = 33;
  constexpr float T_HI = 0.5f - ConwayGraph::T_EPS_AMBO;

  for (const HankinAmboSite &site : HANKIN_AMBO_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t =
          T_EPS + (T_HI - T_EPS) * (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = MeshOps::truncate(seed, a, b, t);
      MeshState compiled;
      MeshOps::compile(swept, compiled, a, b);
      if (s == 0) {
        v0 = swept.vertices.size();
        f0 = swept.face_counts.size();
        i0 = swept.faces.size();
        compiled0 = compiled.face_counts.size();
        expect_op_counts(v0, f0, i0,
                         morph_op_counts(ConwayGraph::MorphOp::TRUNCATE, seed));
      } else {
        HS_EXPECT_EQ(swept.vertices.size(), v0);
        HS_EXPECT_EQ(swept.face_counts.size(), f0);
        HS_EXPECT_EQ(swept.faces.size(), i0);
        HS_EXPECT_EQ(compiled.face_counts.size(), compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_characteristic_two(swept);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [hankin-ambo] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [hankin-ambo] %s: F=%zu compiled=%zu across %d samples\n",
                  site.name, f0, compiled0, SAMPLES);
  }
}
