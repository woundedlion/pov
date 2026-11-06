/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Ambo-on-hankin probe: ambo legs on hankin meshes keep their face counts and
// manifold topology across the sweep.
// ---------------------------------------------------------------------------

inline PolyMesh probe_dodeca_hk62(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::DODECAHEDRON_HK62_AMBO_HK62_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_dodeca_hk72(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_dodeca_hk35(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE,
      Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_dodeca_hk35_ambo_hk62(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE,
      Solids::Op::AMBO, 1>(a, b);
}
inline PolyMesh probe_dodeca_hk54(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::DODECAHEDRON_HK54_AMBO_HK72_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_octa_hk17(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::OCTAHEDRON_HK17_AMBO_HK73_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_octa_hk34(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::OCTAHEDRON_HK34_AMBO_HK72_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_rhombicubocta_hk63(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_ticosa_hk54(Arena &a, Arena &b) {
  return recipe_step_seed<Solids::TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE,
                          Solids::Op::AMBO, 0>(a, b);
}
inline PolyMesh probe_icosa_truncate25_hk60(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::ICOSAHEDRON_TRUNCATE25_HK60_AMBO_RELAX100_HK60_RECIPE,
      Solids::Op::AMBO, 0>(a, b);
}

inline constexpr float HANKIN_AMBO_T = 0.5f - ConwayGraph::T_EPS_AMBO;
inline constexpr StepLegSite HANKIN_AMBO_SITES[] = {
    {"dodecahedron_hk62", probe_dodeca_hk62, HANKIN_AMBO_T},
    {"dodecahedron_hk72", probe_dodeca_hk72, HANKIN_AMBO_T},
    {"dodecahedron_hk35", probe_dodeca_hk35, HANKIN_AMBO_T},
    {"dodecahedron_hk35_ambo_hk62", probe_dodeca_hk35_ambo_hk62, HANKIN_AMBO_T},
    {"dodecahedron_hk54", probe_dodeca_hk54, HANKIN_AMBO_T},
    {"octahedron_hk17", probe_octa_hk17, HANKIN_AMBO_T},
    {"octahedron_hk34", probe_octa_hk34, HANKIN_AMBO_T},
    {"rhombicuboctahedron_hk63", probe_rhombicubocta_hk63, HANKIN_AMBO_T},
    {"truncatedIcosahedron_hk54", probe_ticosa_hk54, HANKIN_AMBO_T},
    {"icosahedron_truncate25_hk60", probe_icosa_truncate25_hk60, HANKIN_AMBO_T},
};

/** @brief Sweeps truncate legs on the ambo-on-hankin recipe seeds. */
inline void test_ambo_leg_on_hankin_seed_holds_topology() {
  check_step_leg_sweep(StepLegKind::TRUNCATE, HANKIN_AMBO_SITES);
}
