/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effects.h
 * @brief Effect roster: pulls in every effect header plus the HS_EFFECT_LIST
 *        X-macro.
 *
 * Include this only from a build target (run_tests, the firmware entry point).
 * An effect must NOT include it — that would pull in all the other effects and
 * recurse. Effects include core/engine/engine.h (the API umbrella) instead.
 */

#include "effects/BZReactionDiffusion.h"
#include "effects/ChaoticStrings.h"
#include "effects/Comets.h"
#include "effects/CurlLattice.h"
#include "effects/DisplacementField.h"
#include "effects/DreamBalls.h"
#include "effects/Dynamo.h"
#include "effects/FacetGrid.h"
#include "effects/GnomonicStars.h"
#include "effects/GSReactionDiffusion.h"
#include "effects/HankinSolids.h"
#include "effects/HopfFibration.h"
#include "effects/IslamicStars.h"
#include "effects/MeshFeedback.h"
#include "effects/MindSplatter.h"
#include "effects/MobiusRings.h"
#include "effects/PetalFlow.h"
#include "effects/PromotedShaderLooks.h"
#include "effects/Raymarch.h"
#include "effects/RingShower.h"
#include "effects/RingSpin.h"
#include "effects/ShapeShifter.h"
#include "effects/SphericalHarmonics.h"
#include "effects/Thrusters.h"
#include "effects/Voronoi.h"

/**
 * @brief Single source of truth for the registered effect roster, as an X-macro.
 * @param X Function-like macro applied to each effect type name in the roster.
 * @details The `#include` list above and this X-macro list must stay in lock-step.
 *     * WASM:   the self-registering EffectRegistry size is checked against
 *               HS_EFFECT_COUNT at engine startup
 *               (targets/wasm/engine_bindings.h), so a registered-but-unlisted
 *               (or listed-but-unregistered) effect traps.
 *     * Native: the effect smoke suite iterates this X-macro list, so its coverage
 *               is derived from the list rather than hand-maintained, and it runs
 *               the same registry-count oracle unconditionally
 *               (tests/test_effects.h), so the same drift fails the suite.
 *   Adding an effect therefore means: add the `#include` above, the
 *   REGISTER_EFFECT in its header, and one X() row here.
 */
#define HS_EFFECT_LIST(X)                                                      \
  X(BZReactionDiffusion)                                                       \
  X(ChaoticStrings)                                                            \
  X(Comets)                                                                    \
  X(ContourLattice)                                                            \
  X(CurlLattice)                                                               \
  X(DisplacementField)                                                         \
  X(DreamBalls)                                                                \
  X(Dynamo)                                                                    \
  X(EquatorGrid)                                                               \
  X(FacetGrid)                                                                 \
  X(FacetWave)                                                                 \
  X(GnomonicStars)                                                             \
  X(GSReactionDiffusion)                                                       \
  X(GlitchGrid)                                                                \
  X(HankinSolids)                                                              \
  X(HopfFibration)                                                             \
  X(IslamicStars)                                                              \
  X(HexWave)                                                                   \
  X(AlienOcean)                                                                \
  X(KaleidoWave)                                                               \
  X(MeshFeedback)                                                              \
  X(MindSplatter)                                                              \
  X(MobiusGrid)                                                                \
  X(MobiusRings)                                                               \
  X(PetalFlow)                                                                 \
  X(PrismLattice)                                                              \
  X(Raymarch)                                                                  \
  X(RingShower)                                                                \
  X(RingSpin)                                                                  \
  X(Shader)                                                                    \
  X(ShapeShifter)                                                              \
  X(SignalWeave)                                                               \
  X(SphericalHarmonics)                                                        \
  X(CosmicEyeball)                                                             \
  X(Thrusters)                                                                 \
  X(VectorFacets)                                                              \
  X(Voronoi)

/**
 * @brief Phantasm playlist: HS_EFFECT_LIST minus Shader and the low-res-only
 *        effects Dynamo, MobiusRings, and Thrusters.
 * @param X Function-like macro applied to each effect type name and its show
 *          duration in seconds.
 * @details Same order as HS_EFFECT_LIST. Only the Phantasm firmware target
 *   consumes this; the registry, tests, and gallery stay on the full roster.
 *   The static_assert below HS_PHANTASM_EFFECT_COUNT forces this list to be
 *   revisited whenever HS_EFFECT_LIST gains or loses an entry.
 */
#define HS_PHANTASM_EFFECT_LIST(X)                                             \
  X(BZReactionDiffusion, 120)                                                  \
  X(ChaoticStrings, 120)                                                       \
  X(Comets, 120)                                                               \
  X(ContourLattice, 6)                                                         \
  X(CurlLattice, 13)                                                           \
  X(DisplacementField, 120)                                                    \
  X(DreamBalls, 120)                                                           \
  X(EquatorGrid, 19)                                                           \
  X(FacetGrid, 19)                                                             \
  X(FacetWave, 6)                                                              \
  X(GnomonicStars, 120)                                                        \
  X(GSReactionDiffusion, 120)                                                  \
  X(GlitchGrid, 6)                                                             \
  X(HankinSolids, 120)                                                         \
  X(HopfFibration, 120)                                                        \
  X(HexWave, 6)                                                                \
  X(IslamicStars, 120)                                                         \
  X(AlienOcean, 7)                                                             \
  X(KaleidoWave, 7)                                                            \
  X(MeshFeedback, 181)                                                         \
  X(MindSplatter, 120)                                                         \
  X(MobiusGrid, 120)                                                           \
  X(PetalFlow, 120)                                                            \
  X(PrismLattice, 6)                                                           \
  X(Raymarch, 120)                                                             \
  X(RingShower, 120)                                                           \
  X(RingSpin, 120)                                                             \
  X(ShapeShifter, 135)                                                         \
  X(SignalWeave, 7)                                                            \
  X(SphericalHarmonics, 120)                                                   \
  X(CosmicEyeball, 6)                                                          \
  X(VectorFacets, 6)                                                           \
  X(Voronoi, 120)

/** Shader promotion product group in gallery and device order. */
#define HS_SHADER_PRODUCT_GROUP(X)                                             \
  X(SignalWeave, 7)                                                            \
  X(KaleidoWave, 7)                                                            \
  X(AlienOcean, 7)                                                             \
  X(GlitchGrid, 6)                                                             \
  X(FacetWave, 6)                                                              \
  X(ContourLattice, 6)                                                         \
  X(CurlLattice, 13)                                                           \
  X(PrismLattice, 6)                                                           \
  X(VectorFacets, 6)                                                           \
  X(FacetGrid, 19)                                                             \
  X(HexWave, 6)                                                                \
  X(EquatorGrid, 19)                                                           \
  X(CosmicEyeball, 6)                                                          \
  X(MobiusGrid, 6)

#define HS_SHADER_GROUP_SECONDS(name, seconds) +seconds
constexpr int HS_SHADER_PRODUCT_GROUP_SECONDS =
    0 HS_SHADER_PRODUCT_GROUP(HS_SHADER_GROUP_SECONDS);
#undef HS_SHADER_GROUP_SECONDS
static_assert(HS_SHADER_PRODUCT_GROUP_SECONDS == 120,
              "the promoted Shader group must retain ShaderBall's airtime");

#define HS_SHADER_GROUP_REACHABLE(name, seconds)                               \
  static_assert(name<96, 20>::PRESET_IDS.size() > 0,                           \
                #name " must expose a reachable preset");
HS_SHADER_PRODUCT_GROUP(HS_SHADER_GROUP_REACHABLE)
#undef HS_SHADER_GROUP_REACHABLE

/**
 * @brief Expands to +1 so HS_EFFECT_LIST can be summed into an entry count.
 * @param name Effect type name supplied by HS_EFFECT_LIST (unused).
 */
#define HS_EFFECT_COUNT_ADD(name) +1
/**
 * @brief Number of entries in HS_EFFECT_LIST, derived rather than hand-counted.
 */
constexpr int HS_EFFECT_COUNT = 0 HS_EFFECT_LIST(HS_EFFECT_COUNT_ADD);
#define HS_PHANTASM_EFFECT_COUNT_ADD(name, duration_seconds) +1
/**
 * @brief Number of entries in HS_PHANTASM_EFFECT_LIST, derived rather than
 *        hand-counted.
 */
constexpr int HS_PHANTASM_EFFECT_COUNT =
    0 HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_EFFECT_COUNT_ADD);
#undef HS_PHANTASM_EFFECT_COUNT_ADD
#undef HS_EFFECT_COUNT_ADD

/**
 * @brief Compile-time equality of two NUL-terminated names.
 * @param a First name.
 * @param b Second name.
 * @return True when both spell the same string.
 */
constexpr bool hs_effect_name_eq(const char *a, const char *b) {
  while (*a && *a == *b) {
    ++a;
    ++b;
  }
  return *a == *b;
}

/**
 * @brief True when `name` is one of the HS_EFFECT_LIST class names.
 * @param name Effect class name to look up.
 */
constexpr bool hs_in_effect_list(const char *name) {
#define HS_EFFECT_NAME_MATCH(cls) || hs_effect_name_eq(name, #cls)
  return false HS_EFFECT_LIST(HS_EFFECT_NAME_MATCH);
#undef HS_EFFECT_NAME_MATCH
}

/**
 * @brief True when `name` is one of the HS_PHANTASM_EFFECT_LIST class names.
 * @param name Effect class name to look up.
 */
constexpr bool hs_in_phantasm_effect_list(const char *name) {
#define HS_PHANTASM_NAME_MATCH(cls, duration_seconds)                          \
  || hs_effect_name_eq(name, #cls)
  return false HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_MATCH);
#undef HS_PHANTASM_NAME_MATCH
}

/**
 * @brief Occurrences of `name` in HS_PHANTASM_EFFECT_LIST.
 * @param name Effect class name to count.
 */
constexpr int hs_phantasm_effect_list_count(const char *name) {
#define HS_PHANTASM_NAME_COUNT(cls, duration_seconds)                          \
  +(hs_effect_name_eq(name, #cls) ? 1 : 0)
  return 0 HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_COUNT);
#undef HS_PHANTASM_NAME_COUNT
}

/** @brief True when no HS_PHANTASM_EFFECT_LIST name appears twice. */
constexpr bool hs_phantasm_effect_list_is_distinct() {
#define HS_PHANTASM_NAME_ONCE(cls, duration_seconds)                           \
  &&(hs_phantasm_effect_list_count(#cls) == 1)
  return true HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_ONCE);
#undef HS_PHANTASM_NAME_ONCE
}

/** @brief True when every HS_PHANTASM_EFFECT_LIST name is in HS_EFFECT_LIST. */
constexpr bool hs_phantasm_effect_list_is_subset() {
#define HS_PHANTASM_NAME_ON_ROSTER(cls, duration_seconds)                      \
  &&hs_in_effect_list(#cls)
  return true HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_ON_ROSTER);
#undef HS_PHANTASM_NAME_ON_ROSTER
}

// Drift guard: an effect added to (or removed from) HS_EFFECT_LIST must also be
// deliberately added to or excluded from the Phantasm playlist above. The count
// pins the cardinality; the name scans pin which entries are missing, so
// swapping one exclusion for another cannot ride green. Distinctness closes the
// last hole: a duplicated entry paired with an omission holds the count and
// passes both exclusion scans while an effect silently drops off the playlist.
// Containment rejects a playlist entry that names no roster effect at all.
static_assert(hs_phantasm_effect_list_is_distinct(),
              "HS_PHANTASM_EFFECT_LIST names an effect twice — the duplicate "
              "is masking an omitted effect");
static_assert(hs_phantasm_effect_list_is_subset(),
              "HS_PHANTASM_EFFECT_LIST names an effect that is not in "
              "HS_EFFECT_LIST — a rename or typo left the playlist off-roster");
static_assert(HS_PHANTASM_EFFECT_COUNT == HS_EFFECT_COUNT - 4,
              "HS_PHANTASM_EFFECT_LIST out of sync with HS_EFFECT_LIST "
              "(full roster minus Shader, Dynamo, MobiusRings and Thrusters)");
static_assert(hs_in_effect_list("Shader") && hs_in_effect_list("Dynamo") &&
                  hs_in_effect_list("MobiusRings") &&
                  hs_in_effect_list("Thrusters"),
              "Phantasm exclusion names a non-roster effect — a rename left "
              "the exclusion guard below vacuous");
static_assert(!hs_in_phantasm_effect_list("Shader") &&
                  !hs_in_phantasm_effect_list("Dynamo") &&
                  !hs_in_phantasm_effect_list("MobiusRings") &&
                  !hs_in_phantasm_effect_list("Thrusters"),
              "HS_PHANTASM_EFFECT_LIST must exclude Shader, Dynamo, "
              "MobiusRings and Thrusters");
