/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Ambo-on-hankin probe: ambo legs on hankin meshes keep their face counts and
// manifold topology across the sweep.
// ---------------------------------------------------------------------------

inline constexpr float HANKIN_AMBO_T = 0.5f - ConwayGraph::T_EPS_AMBO;
inline constexpr StepLegSite HANKIN_AMBO_SITES[] = {
    {"dodecahedron_hk62",
     recipe_step_seed<Solids::DODECAHEDRON_HK62_AMBO_HK62_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"dodecahedron_hk72",
     recipe_step_seed<Solids::DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"dodecahedron_hk35",
     recipe_step_seed<
         Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE,
         Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"dodecahedron_hk35_ambo_hk62",
     recipe_step_seed<
         Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE,
         Solids::Op::AMBO, 1>,
     HANKIN_AMBO_T},
    {"dodecahedron_hk54",
     recipe_step_seed<Solids::DODECAHEDRON_HK54_AMBO_HK72_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"octahedron_hk17",
     recipe_step_seed<Solids::OCTAHEDRON_HK17_AMBO_HK73_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"octahedron_hk34",
     recipe_step_seed<Solids::OCTAHEDRON_HK34_AMBO_HK72_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"rhombicuboctahedron_hk63",
     recipe_step_seed<Solids::RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"truncatedIcosahedron_hk54",
     recipe_step_seed<Solids::TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE,
                      Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
    {"icosahedron_truncate25_hk60",
     recipe_step_seed<
         Solids::ICOSAHEDRON_TRUNCATE25_HK60_AMBO_RELAX100_HK60_RECIPE,
         Solids::Op::AMBO, 0>,
     HANKIN_AMBO_T},
};

/** @brief Sweeps truncate legs on the ambo-on-hankin recipe seeds. */
inline void test_ambo_leg_on_hankin_seed_holds_topology() {
  check_step_leg_sweep(StepLegKind::TRUNCATE, HANKIN_AMBO_SITES);
}
