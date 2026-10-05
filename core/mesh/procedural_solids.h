/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/mesh/solid_generators.h.

// ==========================================================================================
// PROCEDURAL GENERATORS
// ==========================================================================================

namespace Platonic {
/**
 * @brief Builds a tetrahedron (V=4, F=4, I=12).
 * @param a Arena that backs the returned mesh.
 * @return The tetrahedron mesh.
 */
FLASHMEM static PolyMesh tetrahedron(Arena &a, Arena &) {
  return to_polymesh<Tetrahedron>(a);
}
/**
 * @brief Builds a cube (V=8, F=6, I=24).
 * @param a Arena that backs the returned mesh.
 * @return The cube mesh.
 */
FLASHMEM static PolyMesh cube(Arena &a, Arena &) {
  return to_polymesh<Cube>(a);
}
/**
 * @brief Builds an octahedron (V=6, F=8, I=24).
 * @param a Arena that backs the returned mesh.
 * @return The octahedron mesh.
 */
FLASHMEM static PolyMesh octahedron(Arena &a, Arena &) {
  return to_polymesh<Octahedron>(a);
}
/**
 * @brief Builds a dodecahedron (V=20, F=12, I=60).
 * @param a Arena that backs the returned mesh.
 * @return The dodecahedron mesh.
 */
FLASHMEM static PolyMesh dodecahedron(Arena &a, Arena &) {
  return to_polymesh<Dodecahedron>(a);
}
/**
 * @brief Builds an icosahedron (V=12, F=20, I=60).
 * @param a Arena that backs the returned mesh.
 * @return The icosahedron mesh.
 */
FLASHMEM static PolyMesh icosahedron(Arena &a, Arena &) {
  return to_polymesh<Icosahedron>(a);
}
} // namespace Platonic

namespace Archimedean {
using namespace Platonic;

/**
 * @brief Builds a truncated tetrahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated tetrahedron mesh.
 */
FLASHMEM static PolyMesh truncatedTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(tetrahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a cuboctahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The cuboctahedron mesh.
 */
FLASHMEM static PolyMesh cuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).ambo().build();
}
/**
 * @brief Builds a truncated cube.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated cube mesh.
 */
FLASHMEM static PolyMesh truncatedCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).truncate(T_TRUNC_CUBE).build();
}
/**
 * @brief Builds a truncated octahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated octahedron mesh.
 */
FLASHMEM static PolyMesh truncatedOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(octahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a rhombicuboctahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The rhombicuboctahedron mesh.
 */
FLASHMEM static PolyMesh rhombicuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b).expand().build();
}
/**
 * @brief Builds a truncated cuboctahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated cuboctahedron mesh.
 */
FLASHMEM static PolyMesh truncatedCuboctahedron(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b)
      .bevel(T_TRUNC_CUBE)
      .relax_baked(RelaxBakes::truncated_cuboctahedron_converged)
      .build();
}
/**
 * @brief Builds a snub cube.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The snub cube mesh.
 */
FLASHMEM static PolyMesh snubCube(Arena &a, Arena &b) {
  return SolidBuilder(cube(a, b), a, b)
      .snub(T_SNUB_CUBE, SNUB_CUBE_TWIST)
      .relax_baked(RelaxBakes::snub_cube_converged)
      .build();
}
/**
 * @brief Builds an icosidodecahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The icosidodecahedron mesh.
 */
FLASHMEM static PolyMesh icosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).ambo().build();
}
/**
 * @brief Builds a truncated dodecahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated dodecahedron mesh.
 */
FLASHMEM static PolyMesh truncatedDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b).truncate(T_TRUNC_ICOS).build();
}
/**
 * @brief Builds a truncated icosahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated icosahedron mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosahedron(a, b), a, b).truncate(T_TRUNC_THIRD).build();
}
/**
 * @brief Builds a rhombicosidodecahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The rhombicosidodecahedron mesh.
 */
FLASHMEM static PolyMesh rhombicosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .expand()
      .relax_baked(RelaxBakes::rhombicosidodecahedron_converged)
      .build();
}
/**
 * @brief Builds a truncated icosidodecahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The truncated icosidodecahedron mesh.
 */
FLASHMEM static PolyMesh truncatedIcosidodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .bevel(T_TRUNC_ICOS)
      .relax_baked(RelaxBakes::truncated_icosidodecahedron_converged)
      .build();
}
/**
 * @brief Builds a snub dodecahedron.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The snub dodecahedron mesh.
 */
FLASHMEM static PolyMesh snubDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(dodecahedron(a, b), a, b)
      .snub(0.5f)
      .relax_baked(RelaxBakes::snub_dodecahedron_converged)
      .build();
}
} // namespace Archimedean

namespace Catalan {
using namespace Archimedean;

/**
 * @brief Builds a triakis tetrahedron (dual of the truncated tetrahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The triakis tetrahedron mesh.
 */
FLASHMEM static PolyMesh triakisTetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedTetrahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a rhombic dodecahedron (dual of the cuboctahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The rhombic dodecahedron mesh.
 */
FLASHMEM static PolyMesh rhombicDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(cuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a triakis octahedron (dual of the truncated cube).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The triakis octahedron mesh.
 */
FLASHMEM static PolyMesh triakisOctahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCube(a, b), a, b).dual().build();
}
/**
 * @brief Builds a tetrakis hexahedron (dual of the truncated octahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The tetrakis hexahedron mesh.
 */
FLASHMEM static PolyMesh tetrakisHexahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedOctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a deltoidal icositetrahedron (dual of the rhombicuboctahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The deltoidal icositetrahedron mesh.
 */
FLASHMEM static PolyMesh deltoidalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a disdyakis dodecahedron (dual of the truncated cuboctahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The disdyakis dodecahedron mesh.
 */
FLASHMEM static PolyMesh disdyakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedCuboctahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentagonal icositetrahedron (dual of the snub cube).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The pentagonal icositetrahedron mesh.
 */
FLASHMEM static PolyMesh pentagonalIcositetrahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubCube(a, b), a, b).dual().build();
}
/**
 * @brief Builds a rhombic triacontahedron (dual of the icosidodecahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The rhombic triacontahedron mesh.
 */
FLASHMEM static PolyMesh rhombicTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(icosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a triakis icosahedron (dual of the truncated dodecahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The triakis icosahedron mesh.
 */
FLASHMEM static PolyMesh triakisIcosahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedDodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentakis dodecahedron (dual of the truncated icosahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The pentakis dodecahedron mesh.
 */
FLASHMEM static PolyMesh pentakisDodecahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a deltoidal hexecontahedron (dual of the
 * rhombicosidodecahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The deltoidal hexecontahedron mesh.
 */
FLASHMEM static PolyMesh deltoidalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(rhombicosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a disdyakis triacontahedron (dual of the truncated
 * icosidodecahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The disdyakis triacontahedron mesh.
 */
FLASHMEM static PolyMesh disdyakisTriacontahedron(Arena &a, Arena &b) {
  return SolidBuilder(truncatedIcosidodecahedron(a, b), a, b).dual().build();
}
/**
 * @brief Builds a pentagonal hexecontahedron (dual of the snub dodecahedron).
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The pentagonal hexecontahedron mesh.
 */
FLASHMEM static PolyMesh pentagonalHexecontahedron(Arena &a, Arena &b) {
  return SolidBuilder(snubDodecahedron(a, b), a, b).dual().build();
}
} // namespace Catalan

namespace IslamicStarPatterns {

/** Degrees-to-radians conversion factor. */
inline constexpr float D2R = math::PI_F / 180.0f;

/** Truncation depth of the `*_truncate5d_*` recipes, bit-exactly 5.0f * D2R
 * and named for it, consumed by truncate as a dimensionless edge fraction
 * short of the ambo pinch at t = 0.5. */
inline constexpr float TRUNCATE_T_NEAR = 0.0872664601f;
/** Truncation depth of the `*_truncate50d_*` recipes, bit-exactly 50.0f * D2R
 * and named for it, consumed by truncate as a dimensionless edge fraction past
 * the ambo pinch, where the cut faces self-intersect by design. */
inline constexpr float TRUNCATE_T_FAR = 0.87266463f;
static_assert(TRUNCATE_T_NEAR == 5.0f * D2R);
static_assert(TRUNCATE_T_FAR == 50.0f * D2R);

/**
 * @brief Builds the truncatedIcosahedron_hk58_chamfer63 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_hk58_chamfer63(Arena &a,
                                                             Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .hankin(58.0f * D2R)
      .chamfer(0.63f)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk62_ambo_hk62 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk62_ambo_hk62(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(62.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .build();
}
/**
 * @brief Builds the octahedron_hk17_ambo_hk73 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh octahedron_hk17_ambo_hk73(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::octahedron(a, b), a, b)
      .hankin(17.0f * D2R)
      .ambo()
      .hankin(73.0f * D2R)
      .build();
}
/**
 * @brief Builds the icosahedron_kis_gyro star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_kis_gyro(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b).kis().gyro().build();
}
/**
 * @brief Builds the truncatedIcosidodecahedron_truncate50d_ambo_dual star
 * pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosidodecahedron_truncate50d_ambo_dual(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosidodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_FAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the icosidodecahedron_truncate5d_ambo_dual star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosidodecahedron_truncate5d_ambo_dual(Arena &a,
                                                                Arena &b) {
  return SolidBuilder(Archimedean::icosidodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the snubDodecahedron_truncate5d_ambo_dual star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh snubDodecahedron_truncate5d_ambo_dual(Arena &a,
                                                               Arena &b) {
  return SolidBuilder(Archimedean::snubDodecahedron(a, b), a, b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the octahedron_hk34_ambo_hk72 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh octahedron_hk34_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::octahedron(a, b), a, b)
      .hankin(34.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the rhombicuboctahedron_hk63_ambo_hk63 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh rhombicuboctahedron_hk63_ambo_hk63(Arena &a,
                                                            Arena &b) {
  return SolidBuilder(Archimedean::rhombicuboctahedron(a, b), a, b)
      .hankin(63.0f * D2R)
      .ambo()
      .hankin(63.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_hk54_ambo_hk72 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_hk54_ambo_hk72(Arena &a,
                                                             Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk54_ambo_hk72 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk54_ambo_hk72(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .ambo()
      .hankin(72.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk72_ambo_dual_hk20 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_hk72_ambo_dual_hk20(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(72.0f * D2R)
      .ambo()
      .dual()
      .hankin(20.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_truncate50d_ambo_dual star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedIcosahedron_truncate50d_ambo_dual(Arena &a,
                                                                    Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .truncate(TRUNCATE_T_FAR)
      .ambo()
      .dual()
      .build();
}
/**
 * @brief Builds the icosahedron_snub_relax_truncate033_hankin62 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_snub_relax_truncate033_hankin62(Arena &a,
                                                                     Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .snub()
      .relax_baked(RelaxBakes::icosahedron_snub_converged)
      .truncate(0.33f)
      .hankin(62.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 * @details The final contact angle sits clear of the ~43-degree resonance
 * where one corner class's contact planes go near-parallel.
 */
FLASHMEM static PolyMesh dodecahedron_hk35_ambo_hk62_ambo_relax_hk42(Arena &a,
                                                                     Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .ambo()
      .relax_baked(RelaxBakes::dodecahedron_hankin_ambo_hankin_ambo_converged)
      .hankin(42.0f * D2R)
      .build();
}
/**
 * @brief Builds the icosahedron_ambo_truncate033_hankin59 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh icosahedron_ambo_truncate033_hankin59(Arena &a,
                                                               Arena &b) {
  return SolidBuilder(Platonic::icosahedron(a, b), a, b)
      .ambo()
      .truncate(0.33f)
      .hankin(59.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin59 star
 * pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin59(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.01f)
      .hankin(59.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate001_hankin73 star
 * pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate001_hankin73(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.01f)
      .hankin(73.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedOctahedron_gyro_kis_hk17 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh truncatedOctahedron_gyro_kis_hk17(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedOctahedron(a, b), a, b)
      .gyro()
      .kis()
      .hankin(17.0f * D2R)
      .build();
}
/**
 * @brief Builds the truncatedIcosidodecahedron_bevel5_relax_hk77 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosidodecahedron_bevel5_relax_hk77(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosidodecahedron(a, b), a, b)
      .bevel(0.5f)
      .relax_baked(RelaxBakes::truncated_icosidodecahedron_bevel50_relax100)
      .hankin(77.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_bevel2_relax_gyro star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_bevel2_relax_gyro(Arena &a, Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .bevel(0.2f)
      .relax_baked(RelaxBakes::dodecahedron_bevel20_converged)
      .gyro()
      .build();
}
/**
 * @brief Builds the truncatedIcosahedron_ambo_relax_truncate33_hk64 star
 * pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh
truncatedIcosahedron_ambo_relax_truncate33_hk64(Arena &a, Arena &b) {
  return SolidBuilder(Archimedean::truncatedIcosahedron(a, b), a, b)
      .ambo()
      .relax_baked(RelaxBakes::truncated_icosahedron_ambo_converged)
      .truncate(0.33f)
      .hankin(64.0f * D2R)
      .build();
}
/**
 * @brief Builds the dodecahedron_ambo_bevel33_relax_hk66 star pattern.
 * @param a First arena in the alternating construction pair.
 * @param b Second arena; the result may borrow storage from either arena.
 * @return The resulting star-pattern mesh.
 */
FLASHMEM static PolyMesh dodecahedron_ambo_bevel33_relax_hk66(Arena &a,
                                                              Arena &b) {
  return SolidBuilder(Platonic::dodecahedron(a, b), a, b)
      .ambo()
      .bevel(0.33f)
      .relax_baked(RelaxBakes::dodecahedron_ambo_bevel33_converged)
      .hankin(66.0f * D2R)
      .build();
}
} // namespace IslamicStarPatterns
