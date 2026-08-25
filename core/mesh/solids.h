/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file solids.h
 * @brief The Platonic, Archimedean, Catalan and Islamic star-pattern solid
 *        registries: the authored Recipe mirrors, the registry tables and
 *        the name/index lookups over them.
 */

#include "mesh/solid_generators.h"

#include <array>
#include "math/geometry.h"
#include "mesh/mesh.h" // For MeshOps
#include "mesh/relax_bakes_generated.h"
#include <cmath>
#include <string_view>
#include <span>

namespace Solids {

HS_O3_BEGIN

/**
 * @brief Coarse generator-cost hint surfaced to the picker UI.
 * @details Display label only (Complex = Islamic star-pattern registry, Simple
 * = everything else); no runtime path gates on it.
 */
enum class Category { Simple, Complex };

/**
 * @brief Authored Conway operator, including the composite ops.
 * @details expand_to_primitives() (mesh/recipe.h) lowers composites to the
 * primitives plus AMBO; authored recipes keep the composite.
 */
enum class Op : uint8_t {
  TRUNCATE,
  EXPAND,
  SNUB,
  CHAMFER,
  HANKIN,
  RELAX,
  KIS,
  DUAL,
  AMBO,
  BEVEL,
  GYRO,
  META,
  NEEDLE,
  ZIP
};

/**
 * @brief One authored step in a recipe's op chain.
 */
struct OpStep {
  Op op; /**< Operator applied at this step. */
  /**
   * @brief t / contact angle (radians) / RELAX iterations.
   * @details Unread on a RELAX step carrying a `bake`; such steps leave it at
   * zero rather than naming a count the replay never runs. A bake-less RELAX
   * step must name at least one iteration: the zero default would replay as a
   * normalize-only pass-through, so apply_step traps on it.
   */
  float param = 0.0f;
  float twist = 0.0f; /**< SNUB face rotation, radians. */
  /**
   * @brief RELAX bake this step lands on, mirroring the generator's
   * relax_baked() call; null replays `param` live iterations.
   * @details When set, every replay of the step (bitwise gate, on-screen build
   * leg, and eager clean endpoint) resolves through the baked converged mesh
   * instead of `param` smoothing steps, so the recipe lands on the shipped
   * geometry rather than a mid-convergence freeze. Left null for a mid-chain
   * relax on a mesh with no bake, which keeps iterating.
   */
  const MeshOps::RelaxBake *bake = nullptr;
};

/**
 * @brief Declarative op chain mirroring a registry generator.
 * @details Parallel declaration of a generator's chain, proven bitwise-equal
 * to it in tests/test_solids.h. The generators stay the source of truth for
 * shipping geometry.
 */
struct Recipe {
  uint8_t seed;        /**< simple_registry index of the base solid. */
  const OpStep *steps; /**< Authored op chain, applied left to right. */
  uint8_t count;       /**< Number of steps. */
};

/**
 * @brief One named solid in a registry.
 */
struct Entry {
  const char *name; /**< Registry key / display name. */
  PolyMesh (*generate)(
      Arena &a, Arena &b); /**< Generator building into arena pair (a, b). */
  Category category;       /**< Simple or Complex cost class. */
  const Recipe *recipe = nullptr; /**< Declarative chain mirror; null = none. */
};

/**
 * @brief Registry of the Platonic and the 13 Archimedean solids.
 * @details Order is load-bearing:
 * Collections::get_platonic/archimedean_solids() slice this array by fixed
 * offsets (Platonic 0-4, Archimedean 5-17).
 */
inline constexpr Entry simple_registry[] = {

    // Platonic
    {"tetrahedron", Platonic::tetrahedron, Category::Simple},
    {"cube", Platonic::cube, Category::Simple},
    {"octahedron", Platonic::octahedron, Category::Simple},
    {"dodecahedron", Platonic::dodecahedron, Category::Simple},
    {"icosahedron", Platonic::icosahedron, Category::Simple},

    // Archimedean
    {"truncatedTetrahedron", Archimedean::truncatedTetrahedron,
     Category::Simple},
    {"cuboctahedron", Archimedean::cuboctahedron, Category::Simple},
    {"truncatedCube", Archimedean::truncatedCube, Category::Simple},
    {"truncatedOctahedron", Archimedean::truncatedOctahedron, Category::Simple},
    {"rhombicuboctahedron", Archimedean::rhombicuboctahedron, Category::Simple},
    {"truncatedCuboctahedron", Archimedean::truncatedCuboctahedron,
     Category::Simple},
    {"snubCube", Archimedean::snubCube, Category::Simple},
    {"icosidodecahedron", Archimedean::icosidodecahedron, Category::Simple},
    {"truncatedDodecahedron", Archimedean::truncatedDodecahedron,
     Category::Simple},
    {"truncatedIcosahedron", Archimedean::truncatedIcosahedron,
     Category::Simple},
    {"rhombicosidodecahedron", Archimedean::rhombicosidodecahedron,
     Category::Simple},
    {"truncatedIcosidodecahedron", Archimedean::truncatedIcosidodecahedron,
     Category::Simple},
    {"snubDodecahedron", Archimedean::snubDodecahedron, Category::Simple}};

/**
 * @brief Registry of Catalan solids (duals of the Archimedean solids).
 */
inline constexpr Entry catalan_registry[] = {
    // Catalan
    {"triakisTetrahedron", Catalan::triakisTetrahedron, Category::Simple},
    {"rhombicDodecahedron", Catalan::rhombicDodecahedron, Category::Simple},
    {"triakisOctahedron", Catalan::triakisOctahedron, Category::Simple},
    {"tetrakisHexahedron", Catalan::tetrakisHexahedron, Category::Simple},
    {"deltoidalIcositetrahedron", Catalan::deltoidalIcositetrahedron,
     Category::Simple},
    {"disdyakisDodecahedron", Catalan::disdyakisDodecahedron, Category::Simple},
    {"pentagonalIcositetrahedron", Catalan::pentagonalIcositetrahedron,
     Category::Simple},
    {"rhombicTriacontahedron", Catalan::rhombicTriacontahedron,
     Category::Simple},
    {"triakisIcosahedron", Catalan::triakisIcosahedron, Category::Simple},
    {"pentakisDodecahedron", Catalan::pentakisDodecahedron, Category::Simple},
    {"deltoidalHexecontahedron", Catalan::deltoidalHexecontahedron,
     Category::Simple},
    {"disdyakisTriacontahedron", Catalan::disdyakisTriacontahedron,
     Category::Simple},
    {"pentagonalHexecontahedron", Catalan::pentagonalHexecontahedron,
     Category::Simple}};

// Recipe seed indices into simple_registry, pinned to the registry order so a
// reorder fails to compile.
inline constexpr uint8_t SEED_OCTAHEDRON = 2;
inline constexpr uint8_t SEED_DODECAHEDRON = 3;
inline constexpr uint8_t SEED_ICOSAHEDRON = 4;
inline constexpr uint8_t SEED_TRUNCATED_OCTAHEDRON = 8;
inline constexpr uint8_t SEED_RHOMBICUBOCTAHEDRON = 9;
inline constexpr uint8_t SEED_ICOSIDODECAHEDRON = 12;
inline constexpr uint8_t SEED_TRUNCATED_ICOSAHEDRON = 14;
inline constexpr uint8_t SEED_TRUNCATED_ICOSIDODECAHEDRON = 16;
inline constexpr uint8_t SEED_SNUB_DODECAHEDRON = 17;

static_assert(std::string_view(simple_registry[SEED_OCTAHEDRON].name) ==
              "octahedron");
static_assert(std::string_view(simple_registry[SEED_DODECAHEDRON].name) ==
              "dodecahedron");
static_assert(
    std::string_view(simple_registry[SEED_RHOMBICUBOCTAHEDRON].name) ==
    "rhombicuboctahedron");
static_assert(std::string_view(simple_registry[SEED_ICOSIDODECAHEDRON].name) ==
              "icosidodecahedron");
static_assert(std::string_view(simple_registry[SEED_SNUB_DODECAHEDRON].name) ==
              "snubDodecahedron");
static_assert(std::string_view(simple_registry[SEED_ICOSAHEDRON].name) ==
              "icosahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_OCTAHEDRON].name) ==
    "truncatedOctahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_ICOSAHEDRON].name) ==
    "truncatedIcosahedron");
static_assert(
    std::string_view(simple_registry[SEED_TRUNCATED_ICOSIDODECAHEDRON].name) ==
    "truncatedIcosidodecahedron");

/** Step table for dodecahedron_hk62_ambo_hk62. */
inline constexpr OpStep DODECAHEDRON_HK62_AMBO_HK62_STEPS[] = {
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk62_ambo_hk62. */
inline constexpr Recipe DODECAHEDRON_HK62_AMBO_HK62_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK62_AMBO_HK62_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK62_AMBO_HK62_STEPS))};

/** Step table for octahedron_hk17_ambo_hk73. */
inline constexpr OpStep OCTAHEDRON_HK17_AMBO_HK73_STEPS[] = {
    {Op::HANKIN, 17.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 73.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::octahedron_hk17_ambo_hk73. */
inline constexpr Recipe OCTAHEDRON_HK17_AMBO_HK73_RECIPE = {
    SEED_OCTAHEDRON, OCTAHEDRON_HK17_AMBO_HK73_STEPS,
    static_cast<uint8_t>(std::size(OCTAHEDRON_HK17_AMBO_HK73_STEPS))};

/** Step table for octahedron_hk34_ambo_hk72. */
inline constexpr OpStep OCTAHEDRON_HK34_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 34.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::octahedron_hk34_ambo_hk72. */
inline constexpr Recipe OCTAHEDRON_HK34_AMBO_HK72_RECIPE = {
    SEED_OCTAHEDRON, OCTAHEDRON_HK34_AMBO_HK72_STEPS,
    static_cast<uint8_t>(std::size(OCTAHEDRON_HK34_AMBO_HK72_STEPS))};

/** Step table for dodecahedron_hk54_ambo_hk72. */
inline constexpr OpStep DODECAHEDRON_HK54_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 54.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk54_ambo_hk72. */
inline constexpr Recipe DODECAHEDRON_HK54_AMBO_HK72_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK54_AMBO_HK72_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK54_AMBO_HK72_STEPS))};

/** Step table for icosahedron_ambo_truncate033_hankin59. */
inline constexpr OpStep ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS[] = {
    {Op::AMBO},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 59.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59.
 */
inline constexpr Recipe ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_STEPS))};

/** Step table for rhombicuboctahedron_hk63_ambo_hk63. */
inline constexpr OpStep RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS[] = {
    {Op::HANKIN, 63.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 63.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63. */
inline constexpr Recipe RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE = {
    SEED_RHOMBICUBOCTAHEDRON, RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS,
    static_cast<uint8_t>(std::size(RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_STEPS))};

/** Step table for icosahedron_kis_gyro. */
inline constexpr OpStep ICOSAHEDRON_KIS_GYRO_STEPS[] = {{Op::KIS}, {Op::GYRO}};
/** Recipe mirror of IslamicStarPatterns::icosahedron_kis_gyro. */
inline constexpr Recipe ICOSAHEDRON_KIS_GYRO_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_KIS_GYRO_STEPS,
    static_cast<uint8_t>(std::size(ICOSAHEDRON_KIS_GYRO_STEPS))};

/** Step table for dodecahedron_hk72_ambo_dual_hk20. */
inline constexpr OpStep DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS[] = {
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::DUAL},
    {Op::HANKIN, 20.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20. */
inline constexpr Recipe DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_HK72_AMBO_DUAL_HK20_STEPS))};

/** Step table for icosidodecahedron_truncate5d_ambo_dual. */
inline constexpr OpStep ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_NEAR},
    {Op::AMBO},
    {Op::DUAL}};
/** Recipe mirror of IslamicStarPatterns::icosidodecahedron_truncate5d_ambo_dual. */
inline constexpr Recipe ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE = {
    SEED_ICOSIDODECAHEDRON, ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS))};

/** Step table for snubDodecahedron_truncate5d_ambo_dual. */
inline constexpr OpStep SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_NEAR},
    {Op::AMBO},
    {Op::DUAL}};
/** Recipe mirror of IslamicStarPatterns::snubDodecahedron_truncate5d_ambo_dual. */
inline constexpr Recipe SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE = {
    SEED_SNUB_DODECAHEDRON, SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_STEPS))};

/** Step table for dodecahedron_ambo_bevel33_relax_hk66. */
inline constexpr OpStep DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS[] = {
    {Op::AMBO},
    {Op::BEVEL, 0.33f},
    {.op = Op::RELAX, .bake = &RelaxBakes::dodecahedron_ambo_bevel33_converged},
    {Op::HANKIN, 66.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_ambo_bevel33_relax_hk66. */
inline constexpr Recipe DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS,
    static_cast<uint8_t>(
        std::size(DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_STEPS))};

/** Step table for truncatedIcosahedron_hk58_chamfer63. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS[] = {
    {Op::HANKIN, 58.0f * IslamicStarPatterns::D2R}, {Op::CHAMFER, 0.63f}};
/** Recipe mirror of IslamicStarPatterns::truncatedIcosahedron_hk58_chamfer63. */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON, TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate33_hk64. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.33f},
        {Op::HANKIN, 64.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate33_hk64.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS,
        static_cast<uint8_t>(
            std::size(TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_STEPS))};

/** Step table for dodecahedron_bevel2_relax_gyro. */
inline constexpr OpStep DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS[] = {
    {Op::BEVEL, 0.2f},
    {.op = Op::RELAX, .bake = &RelaxBakes::dodecahedron_bevel20_converged},
    {Op::GYRO}};
/** Recipe mirror of IslamicStarPatterns::dodecahedron_bevel2_relax_gyro. */
inline constexpr Recipe DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS,
    static_cast<uint8_t>(std::size(DODECAHEDRON_BEVEL2_RELAX_GYRO_STEPS))};

/** Step table for truncatedIcosidodecahedron_bevel5_relax_hk77. */
inline constexpr OpStep TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS[] =
    {{Op::BEVEL, 0.5f},
     {.op = Op::RELAX,
      .bake = &RelaxBakes::truncated_icosidodecahedron_bevel50_relax100},
     {Op::HANKIN, 77.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosidodecahedron_bevel5_relax_hk77.
 */
inline constexpr Recipe TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_RECIPE = {
    SEED_TRUNCATED_ICOSIDODECAHEDRON,
    TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_STEPS))};

/** Step table for truncatedOctahedron_gyro_kis_hk17. */
inline constexpr OpStep TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS[] = {
    {Op::GYRO}, {Op::KIS}, {Op::HANKIN, 17.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17. */
inline constexpr Recipe TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE = {
    SEED_TRUNCATED_OCTAHEDRON, TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS,
    static_cast<uint8_t>(std::size(TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate001_hankin59. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.01f},
        {Op::HANKIN, 59.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin59.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_STEPS))};

/** Step table for truncatedIcosahedron_ambo_relax_truncate001_hankin73. */
inline constexpr OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS[] = {
        {Op::AMBO},
        {.op = Op::RELAX,
         .bake = &RelaxBakes::truncated_icosahedron_ambo_converged},
        {Op::TRUNCATE, 0.01f},
        {Op::HANKIN, 73.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin73.
 */
inline constexpr Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_RECIPE = {
        SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_STEPS))};

/** Step table for dodecahedron_hk35_ambo_hk62_ambo_relax_hk42. */
inline constexpr OpStep DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS[] = {
    {Op::HANKIN, 35.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {.op = Op::RELAX,
     .bake = &RelaxBakes::dodecahedron_hankin_ambo_hankin_ambo_converged},
    {Op::HANKIN, 42.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_relax_hk42.
 */
inline constexpr Recipe DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE = {
    SEED_DODECAHEDRON, DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS,
    static_cast<uint8_t>(
        std::size(DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS))};

/** Step table for truncatedIcosidodecahedron_truncate50d_ambo_dual. */
inline constexpr OpStep
    TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
        {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_FAR},
        {Op::AMBO},
        {Op::DUAL}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosidodecahedron_truncate50d_ambo_dual.
 */
inline constexpr Recipe
    TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE = {
        SEED_TRUNCATED_ICOSIDODECAHEDRON,
        TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS,
        static_cast<uint8_t>(std::size(
            TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS))};

/** Step table for truncatedIcosahedron_hk54_ambo_hk72. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 54.0f * IslamicStarPatterns::D2R},
    {Op::AMBO},
    {Op::HANKIN, 72.0f * IslamicStarPatterns::D2R}};
/** Recipe mirror of IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72. */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON, TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_STEPS))};

/** Step table for truncatedIcosahedron_truncate50d_ambo_dual. */
inline constexpr OpStep TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, IslamicStarPatterns::TRUNCATE_T_FAR},
    {Op::AMBO},
    {Op::DUAL}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::truncatedIcosahedron_truncate50d_ambo_dual.
 */
inline constexpr Recipe TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE = {
    SEED_TRUNCATED_ICOSAHEDRON,
    TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS,
    static_cast<uint8_t>(
        std::size(TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_STEPS))};

/** Step table for icosahedron_snub_relax_truncate033_hankin62. */
inline constexpr OpStep ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS[] = {
    {Op::SNUB, 0.5f, 0.0f},
    {.op = Op::RELAX, .bake = &RelaxBakes::icosahedron_snub_converged},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 62.0f * IslamicStarPatterns::D2R}};
/**
 * Recipe mirror of
 * IslamicStarPatterns::icosahedron_snub_relax_truncate033_hankin62.
 */
inline constexpr Recipe ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE = {
    SEED_ICOSAHEDRON, ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS,
    static_cast<uint8_t>(
        std::size(ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_STEPS))};

/**
 * @brief Registry of Islamic star-pattern solids.
 */
inline constexpr Entry islamic_registry[] = {
    {"dodecahedron_hk62_ambo_hk62",
     IslamicStarPatterns::dodecahedron_hk62_ambo_hk62, Category::Complex,
     &DODECAHEDRON_HK62_AMBO_HK62_RECIPE},
    {"truncatedIcosahedron_hk58_chamfer63",
     IslamicStarPatterns::truncatedIcosahedron_hk58_chamfer63,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_HK58_CHAMFER63_RECIPE},
    {"dodecahedron_ambo_bevel33_relax_hk66",
     IslamicStarPatterns::dodecahedron_ambo_bevel33_relax_hk66,
     Category::Complex, &DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate33_hk64",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate33_hk64,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE},
    {"dodecahedron_bevel2_relax_gyro",
     IslamicStarPatterns::dodecahedron_bevel2_relax_gyro, Category::Complex,
     &DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE},
    {"truncatedIcosidodecahedron_bevel5_relax_hk77",
     IslamicStarPatterns::truncatedIcosidodecahedron_bevel5_relax_hk77,
     Category::Complex, &TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_RECIPE},
    {"truncatedOctahedron_gyro_kis_hk17",
     IslamicStarPatterns::truncatedOctahedron_gyro_kis_hk17, Category::Complex,
     &TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin59",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin59,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin73",
     IslamicStarPatterns::truncatedIcosahedron_ambo_relax_truncate001_hankin73,
     Category::Complex,
     &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN73_RECIPE},
    {"icosahedron_ambo_truncate033_hankin59",
     IslamicStarPatterns::icosahedron_ambo_truncate033_hankin59,
     Category::Complex, &ICOSAHEDRON_AMBO_TRUNCATE033_HANKIN59_RECIPE},
    {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42",
     IslamicStarPatterns::dodecahedron_hk35_ambo_hk62_ambo_relax_hk42,
     Category::Complex, &DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE},
    {"octahedron_hk17_ambo_hk73",
     IslamicStarPatterns::octahedron_hk17_ambo_hk73, Category::Complex,
     &OCTAHEDRON_HK17_AMBO_HK73_RECIPE},
    {"icosahedron_kis_gyro", IslamicStarPatterns::icosahedron_kis_gyro,
     Category::Complex, &ICOSAHEDRON_KIS_GYRO_RECIPE},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosidodecahedron_truncate50d_ambo_dual,
     Category::Complex,
     &TRUNCATED_ICOSIDODECAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE},
    {"icosidodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::icosidodecahedron_truncate5d_ambo_dual,
     Category::Complex, &ICOSIDODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE},
    {"snubDodecahedron_truncate5d_ambo_dual",
     IslamicStarPatterns::snubDodecahedron_truncate5d_ambo_dual,
     Category::Complex, &SNUB_DODECAHEDRON_TRUNCATE5D_AMBO_DUAL_RECIPE},
    {"octahedron_hk34_ambo_hk72",
     IslamicStarPatterns::octahedron_hk34_ambo_hk72, Category::Complex,
     &OCTAHEDRON_HK34_AMBO_HK72_RECIPE},
    {"rhombicuboctahedron_hk63_ambo_hk63",
     IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63, Category::Complex,
     &RHOMBICUBOCTAHEDRON_HK63_AMBO_HK63_RECIPE},
    {"truncatedIcosahedron_hk54_ambo_hk72",
     IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_HK54_AMBO_HK72_RECIPE},
    {"dodecahedron_hk54_ambo_hk72",
     IslamicStarPatterns::dodecahedron_hk54_ambo_hk72, Category::Complex,
     &DODECAHEDRON_HK54_AMBO_HK72_RECIPE},
    {"dodecahedron_hk72_ambo_dual_hk20",
     IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20, Category::Complex,
     &DODECAHEDRON_HK72_AMBO_DUAL_HK20_RECIPE},
    {"truncatedIcosahedron_truncate50d_ambo_dual",
     IslamicStarPatterns::truncatedIcosahedron_truncate50d_ambo_dual,
     Category::Complex, &TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL_RECIPE},
    {"icosahedron_snub_relax_truncate033_hankin62",
     IslamicStarPatterns::icosahedron_snub_relax_truncate033_hankin62,
     Category::Complex, &ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE}};

/** Total number of solids across all three registries. */
inline constexpr int NUM_ENTRIES =
    sizeof(simple_registry) / sizeof(simple_registry[0]) +
    sizeof(catalan_registry) / sizeof(catalan_registry[0]) +
    sizeof(islamic_registry) / sizeof(islamic_registry[0]);

// simple_registry is [Platonic | Archimedean]; the static_asserts check the two
// counts exactly tile it and name the entries either side of the slice
// boundary, so a boundary move can't silently mis-slice.
inline constexpr size_t PLATONIC_COUNT = 5;
inline constexpr size_t ARCHIMEDEAN_COUNT = 13;
static_assert(PLATONIC_COUNT + ARCHIMEDEAN_COUNT == std::size(simple_registry),
              "PLATONIC_COUNT + ARCHIMEDEAN_COUNT must equal simple_registry "
              "size; update the counts if the registry layout changes");
static_assert(std::string_view(simple_registry[PLATONIC_COUNT - 1].name) ==
                  "icosahedron",
              "PLATONIC_COUNT must end the Platonic run");
static_assert(std::string_view(simple_registry[PLATONIC_COUNT].name) ==
                  "truncatedTetrahedron",
              "PLATONIC_COUNT must start the Archimedean run");

inline constexpr size_t CATALAN_COUNT = 13;
inline constexpr size_t ISLAMIC_COUNT = 23;
static_assert(CATALAN_COUNT == std::size(catalan_registry),
              "catalan_registry size changed; update CATALAN_COUNT and the "
              "README registry table");
static_assert(ISLAMIC_COUNT == std::size(islamic_registry),
              "islamic_registry size changed; update ISLAMIC_COUNT and the "
              "README registry table");
static_assert(PLATONIC_COUNT + ARCHIMEDEAN_COUNT + CATALAN_COUNT +
                      ISLAMIC_COUNT ==
                  static_cast<size_t>(NUM_ENTRIES),
              "NUM_ENTRIES must equal the sum of the per-registry counts");

// A BaseMesh ordinal is a get_entry() global index: simple_registry then
// catalan_registry. MeshFeedback and DreamBalls resolve solids through it and
// MindSplatter slices its Platonic head, so the enum, the label arrays and the
// two registries share one order. The static_asserts name the entries either
// side of both slice boundaries and at both ends, so inserting, removing or
// reordering a registry entry or an enumerator fails to compile.
static_assert(PLATONIC_BASE_MESH_COUNT == PLATONIC_COUNT,
              "PLATONIC_BASE_MESH_COUNT must equal the Platonic run of "
              "simple_registry");
static_assert(BASE_MESH_COUNT ==
                  PLATONIC_COUNT + ARCHIMEDEAN_COUNT + CATALAN_COUNT,
              "BaseMesh must enumerate simple_registry then catalan_registry");
static_assert(
    std::string_view(simple_registry[static_cast<size_t>(BaseMesh::TETRAHEDRON)]
                         .name) == "tetrahedron",
    "BaseMesh::TETRAHEDRON must open simple_registry");
static_assert(
    std::string_view(simple_registry[static_cast<size_t>(BaseMesh::ICOSAHEDRON)]
                         .name) == "icosahedron",
    "BaseMesh::ICOSAHEDRON must end the Platonic run");
static_assert(
    std::string_view(
        simple_registry[static_cast<size_t>(BaseMesh::TRUNCATED_TETRAHEDRON)]
            .name) == "truncatedTetrahedron",
    "BaseMesh::TRUNCATED_TETRAHEDRON must start the Archimedean run");
static_assert(
    std::string_view(
        simple_registry[static_cast<size_t>(BaseMesh::SNUB_DODECAHEDRON)]
            .name) == "snubDodecahedron",
    "BaseMesh::SNUB_DODECAHEDRON must end simple_registry");
static_assert(
    std::string_view(
        catalan_registry[static_cast<size_t>(BaseMesh::TRIAKIS_TETRAHEDRON) -
                         std::size(simple_registry)]
            .name) == "triakisTetrahedron",
    "BaseMesh::TRIAKIS_TETRAHEDRON must open catalan_registry");
static_assert(
    std::string_view(
        catalan_registry[static_cast<size_t>(
                             BaseMesh::PENTAGONAL_HEXECONTAHEDRON) -
                         std::size(simple_registry)]
            .name) == "pentagonalHexecontahedron",
    "BaseMesh::PENTAGONAL_HEXECONTAHEDRON must close catalan_registry");

namespace Collections {
/**
 * @brief Returns the five Platonic solids.
 * @return Span over the Platonic entries (offset 0, count 5) of
 * simple_registry.
 */
inline std::span<const Entry> get_platonic_solids() {
  return std::span<const Entry>(simple_registry, PLATONIC_COUNT);
}
/**
 * @brief Returns the 13 Archimedean solids.
 * @return Span over the Archimedean entries (offset 5, count 13) of
 * simple_registry.
 */
inline std::span<const Entry> get_archimedean_solids() {
  return std::span<const Entry>(simple_registry + PLATONIC_COUNT,
                                ARCHIMEDEAN_COUNT);
}
/**
 * @brief Returns all simple (Platonic and Archimedean) solids.
 * @return Span over the entire simple_registry.
 */
inline std::span<const Entry> get_simple_solids() {
  return std::span<const Entry>(simple_registry);
}
/**
 * @brief Returns all Catalan solids.
 * @return Span over the entire catalan_registry.
 */
inline std::span<const Entry> get_catalan_solids() {
  return std::span<const Entry>(catalan_registry);
}
/**
 * @brief Returns all Islamic star-pattern solids.
 * @return Span over the entire islamic_registry.
 */
inline std::span<const Entry> get_islamic_solids() {
  return std::span<const Entry>(islamic_registry);
}
} // namespace Collections

/**
 * @brief The three solid registries in flat global-index order.
 * @return Spans over the simple, then Catalan, then Islamic registries.
 * @details Single source of truth for the registry enumeration order; the
 * index- and name-based lookups below all derive from it.
 */
inline constexpr std::array<std::span<const Entry>, 3> all_registries() {
  return {std::span<const Entry>(simple_registry),
          std::span<const Entry>(catalan_registry),
          std::span<const Entry>(islamic_registry)};
}

/**
 * @brief Tests whether every registered recipe's seed indexes simple_registry.
 * @return True when no entry carries an out-of-range Recipe::seed.
 */
inline constexpr bool all_recipe_seeds_in_range() {
  for (std::span<const Entry> reg : all_registries())
    for (const Entry &entry : reg)
      if (entry.recipe && entry.recipe->seed >= std::size(simple_registry))
        return false;
  return true;
}
static_assert(all_recipe_seeds_in_range(),
              "Recipe::seed must index simple_registry; the replay sites index "
              "it without a runtime bound");

/**
 * @brief Tests whether every entry's Category matches the registry it sits in.
 * @return True when islamic_registry is uniformly Complex and the other two are
 * uniformly Simple.
 */
inline constexpr bool all_categories_match_registry() {
  for (const Entry &entry : simple_registry)
    if (entry.category != Category::Simple)
      return false;
  for (const Entry &entry : catalan_registry)
    if (entry.category != Category::Simple)
      return false;
  for (const Entry &entry : islamic_registry)
    if (entry.category != Category::Complex)
      return false;
  return true;
}
static_assert(all_categories_match_registry(),
              "Category must be Complex on islamic_registry and Simple "
              "everywhere else");

/**
 * @brief Tests whether every Islamic star-pattern entry carries a Recipe mirror.
 * @return True when no islamic_registry entry has a null recipe.
 */
inline constexpr bool all_islamic_entries_have_recipes() {
  for (const Entry &entry : islamic_registry)
    if (!entry.recipe)
      return false;
  return true;
}
static_assert(all_islamic_entries_have_recipes(),
              "every islamic_registry entry needs its Recipe mirror; the "
              "lowering and morph-feasibility sweeps skip recipe-less entries");

/**
 * @brief Finds a registry entry by name across all registries.
 * @param name Candidate solid name.
 * @return Pointer to the matching entry, or nullptr when no name matches.
 */
inline const Entry *find_entry(std::string_view name) {
  for (std::span<const Entry> reg : all_registries())
    for (const Entry &entry : reg)
      if (name == entry.name)
        return &entry;
  return nullptr;
}

/**
 * @brief Looks up a registry entry by global index across all three registries.
 * @param index Zero-based index in [0, NUM_ENTRIES); traps if out of range.
 * @return Reference to the entry at that index.
 */
inline const Entry &get_entry(size_t index) {
  HS_CHECK(index < static_cast<size_t>(NUM_ENTRIES),
           "Solids::get_entry: index out of range");

  for (std::span<const Entry> reg : all_registries()) {
    if (index < reg.size())
      return reg[index];
    index -= reg.size();
  }
  __builtin_unreachable();
}

/**
 * @brief Builds the solid with the given name into the geometry arena.
 * @param geom Long-lived arena that backs the returned mesh.
 * @param a Output arena for even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @param name Registry name of the solid to build; traps if unknown.
 * @return The finalized solid mesh owned by geom.
 * @details For trusted (firmware) callers; an unknown name fails fast.
 */
[[maybe_unused]] FLASHMEM static PolyMesh
get_by_name(Arena &geom, Arena &a, Arena &b, std::string_view name) {
  const Entry *entry = find_entry(name);
  HS_CHECK(entry, "Solids::get_by_name: unknown solid name");
  return finalize_solid(entry->generate(a, b), geom);
}

/**
 * @brief Builds a registry solid's unit vertex directions plus per-vertex
 *        orientation quaternions and nearest-neighbour gaps.
 * @param scratch Arena backing the intermediate mesh; nothing is retained
 *        after return.
 * @param temp Alternate scratch arena for odd pipeline stages.
 * @param name Registry name of the solid; traps if unknown.
 * @param max_points Capacity of the output arrays; traps if exceeded.
 * @param points Out: vertex directions projected onto the unit sphere (Catalan
 *        vertices sit at multiple radii).
 * @param quats Out: per-vertex Y-axis-to-direction rotations.
 * @param nn_angle Out: per-vertex nearest-neighbour angle (radians), for
 *        sizing per-vertex geometry to its local gap.
 * @return The vertex count written.
 * @details Shared by the volume-scatter effects (Raymarch). HS_COLD:
 * setup-only, keeps the build loops out of ITCM.
 */
[[maybe_unused]] HS_COLD static int
build_vertex_directions(Arena &scratch, Arena &temp, std::string_view name,
                        int max_points, Vector *points, Quaternion *quats,
                        float *nn_angle) {
  const Entry *entry = find_entry(name);
  HS_CHECK(entry, "build_vertex_directions: unknown solid name");
  // Read straight out of the generator's arena pair; nothing outlives the call,
  // so finalize_solid's long-lived copy would buy nothing.
  PolyMesh mesh = entry->generate(scratch, temp);
  int count = static_cast<int>(mesh.vertices.size());
  HS_CHECK(count <= max_points,
           "build_vertex_directions: vertex count exceeds capacity");
  for (int i = 0; i < count; ++i) {
    points[i] = mesh.vertices[i].normalized();
    quats[i] = make_rotation(Y_AXIS, points[i]);
  }
  // O(n^2) nearest-neighbor scan is intentional: cold setup path, vertex counts
  // are small, so the KD-tree's build overhead is not worth it here.
  for (int i = 0; i < count; ++i) {
    float max_dot = -1.0f;
    for (int j = 0; j < count; ++j)
      if (j != i)
        max_dot = std::max(max_dot, dot(points[i], points[j]));
    nn_angle[i] = acosf(hs::clamp(max_dot, -1.0f, 1.0f));
  }
  return count;
}

HS_O3_END

} // namespace Solids
