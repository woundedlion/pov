/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "engine/effects.h"

/**
 * @file phantasm_playlist.h
 * @brief Phantasm device playlist: the HS_PHANTASM_EFFECT_LIST X-macro with
 *        per-entry show durations, plus its roster drift guards.
 *
 * Only the Phantasm firmware target and the roster cross-checks consume this;
 * the registry, tests, and gallery stay on HS_EFFECT_LIST's full roster.
 * tools/docs_check.py, scripts/effect_roster.mjs, and tools/profile_sweep.sh
 * parse the macro body out of this file's source text.
 */

/**
 * @brief Phantasm playlist: HS_EFFECT_LIST minus Shader, ShaderChain, and the
 *        low-res-only effects Dynamo, MobiusRings, and Thrusters.
 * @param X Function-like macro applied to each effect type name and its show
 *          duration in seconds.
 * @details Entry order is the device show order, chosen independently of
 *   HS_EFFECT_LIST. The static_assert below HS_PHANTASM_EFFECT_COUNT forces
 *   this list to be revisited whenever HS_EFFECT_LIST gains or loses an entry.
 */
#define HS_PHANTASM_EFFECT_LIST(X)                                             \
  X(BZReactionDiffusion, 120)                                                  \
  X(Fishbowl, 120)                                                             \
  X(Comets, 120)                                                               \
  X(GridSpace, hs_preset_window_seconds<GridSpace>())                          \
  X(HyperLattice, hs_preset_window_seconds<HyperLattice>())                    \
  X(AshCloud, hs_preset_window_seconds<AshCloud>())                            \
  X(LatticeMelt, hs_preset_window_seconds<LatticeMelt>())                      \
  X(DisplacementField, 120)                                                    \
  X(DreamBalls, 120)                                                           \
  X(KaleidoscopeFlowers, hs_preset_window_seconds<KaleidoscopeFlowers>())      \
  X(KaleidoscopeSmooth, hs_preset_window_seconds<KaleidoscopeSmooth>())        \
  X(KaleidoscopeMandala, hs_preset_window_seconds<KaleidoscopeMandala>())      \
  X(GnomonicStars, 120)                                                        \
  X(GSReactionDiffusion, 120)                                                  \
  X(AlienCore, hs_preset_window_seconds<AlienCore>())                          \
  X(HankinSolids, 120)                                                         \
  X(HopfFibration, 120)                                                        \
  X(KaleidoscopeHexBright, hs_preset_window_seconds<KaleidoscopeHexBright>())  \
  X(IslamicStars, 120)                                                         \
  X(AlienOcean, hs_preset_window_seconds<AlienOcean>())                        \
  X(KaleidoscopeHexSoft, hs_preset_window_seconds<KaleidoscopeHexSoft>())      \
  X(MeshFeedback, 181)                                                         \
  X(MindSplatter, 120)                                                         \
  X(MobiusGrid, hs_preset_window_seconds<MobiusGrid>())                        \
  X(PetalFlow, 120)                                                            \
  X(KaleidoscopePentBright,                                                    \
    hs_preset_window_seconds<KaleidoscopePentBright>())                        \
  X(KaleidoscopeHexOil, hs_preset_window_seconds<KaleidoscopeHexOil>())        \
  X(Raymarch, 120)                                                             \
  X(RingShower, 120)                                                           \
  X(RingSpin, 120)                                                             \
  X(ShapeShifter, 135)                                                         \
  X(AlienBrain, hs_preset_window_seconds<AlienBrain>())                        \
  X(SphericalHarmonics, 120)                                                   \
  X(CosmicEyeball, hs_preset_window_seconds<CosmicEyeball>())                  \
  X(KaleidoscopeStainedGlass,                                                  \
    hs_preset_window_seconds<KaleidoscopeStainedGlass>())                      \
  X(Voronoi, 120)

#define HS_PHANTASM_EFFECT_COUNT_ADD(name, duration_seconds) +1
/**
 * @brief Number of entries in HS_PHANTASM_EFFECT_LIST, derived rather than
 *        hand-counted.
 */
constexpr int HS_PHANTASM_EFFECT_COUNT =
    0 HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_EFFECT_COUNT_ADD);
#undef HS_PHANTASM_EFFECT_COUNT_ADD

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

/**
 * @brief Show duration HS_PHANTASM_EFFECT_LIST assigns `name`, 0 when absent.
 * @param name Effect class name to look up.
 */
constexpr int hs_phantasm_duration_seconds(const char *name) {
#define HS_PHANTASM_NAME_DURATION(cls, duration_seconds)                       \
  +(hs_effect_name_eq(name, #cls) ? (duration_seconds) : 0)
  return 0 HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_DURATION);
#undef HS_PHANTASM_NAME_DURATION
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
static_assert(HS_PHANTASM_EFFECT_COUNT == HS_EFFECT_COUNT - 3 -
                                              HS_ENABLE_SHADER_WORKBENCH -
                                              HS_ENABLE_CHAIN_INTERPRETER,
              "HS_PHANTASM_EFFECT_LIST out of sync with HS_EFFECT_LIST "
              "(full roster minus Shader, ShaderChain, Dynamo, MobiusRings "
              "and Thrusters)");
static_assert((!HS_ENABLE_SHADER_WORKBENCH || hs_in_effect_list("Shader")) &&
                  (!HS_ENABLE_CHAIN_INTERPRETER ||
                   hs_in_effect_list("ShaderChain")) &&
                  hs_in_effect_list("Dynamo") &&
                  hs_in_effect_list("MobiusRings") &&
                  hs_in_effect_list("Thrusters"),
              "Phantasm exclusion names a non-roster effect — a rename left "
              "the exclusion guard below vacuous");
static_assert(!hs_in_phantasm_effect_list("Shader") &&
                  !hs_in_phantasm_effect_list("ShaderChain") &&
                  !hs_in_phantasm_effect_list("Dynamo") &&
                  !hs_in_phantasm_effect_list("MobiusRings") &&
                  !hs_in_phantasm_effect_list("Thrusters"),
              "HS_PHANTASM_EFFECT_LIST must exclude Shader, ShaderChain, "
              "Dynamo, MobiusRings and Thrusters");

// Product-group durations mirror the Phantasm playlist.
#define HS_SHADER_GROUP_DURATION_MATCHES(cls, duration_seconds)                \
  static_assert(hs_phantasm_duration_seconds(#cls) == (duration_seconds), #cls \
                " duration disagrees between HS_SHADER_PRODUCT_GROUP and "     \
                "HS_PHANTASM_EFFECT_LIST");
HS_SHADER_PRODUCT_GROUP(HS_SHADER_GROUP_DURATION_MATCHES)
#undef HS_SHADER_GROUP_DURATION_MATCHES
