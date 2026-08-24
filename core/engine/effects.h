/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

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
#include "effects/AlienOcean.h"
#include "effects/Fishbowl.h"
#include "effects/Comets.h"
#include "effects/AshCloud.h"
#include "effects/LatticeMelt.h"
#include "effects/DisplacementField.h"
#include "effects/DreamBalls.h"
#include "effects/Dynamo.h"
#include "effects/GridSpace.h"
#include "effects/HyperLattice.h"
#include "effects/CosmicEyeball.h"
#include "effects/KaleidoscopeFlowers.h"
#include "effects/KaleidoscopeSmooth.h"
#include "effects/KaleidoscopeMandala.h"
#include "effects/AlienCore.h"
#include "effects/KaleidoscopeHexBright.h"
#include "effects/KaleidoscopeHexSoft.h"
#include "effects/MobiusGrid.h"
#include "effects/KaleidoscopePentBright.h"
#include "effects/KaleidoscopeHexOil.h"
#include "effects/AlienBrain.h"
#include "effects/KaleidoscopeStainedGlass.h"
#include "effects/GnomonicStars.h"
#include "effects/GSReactionDiffusion.h"
#include "effects/HankinSolids.h"
#include "effects/HopfFibration.h"
#include "effects/IslamicStars.h"
#include "effects/MeshFeedback.h"
#include "effects/MindSplatter.h"
#include "effects/MobiusRings.h"
#include "effects/PetalFlow.h"
#if HS_ENABLE_SHADER_WORKBENCH
#include "workbench/ShaderWorkbench.h"
#define HS_SHADER_WORKBENCH_EFFECT(X) X(Shader)
#else
#define HS_SHADER_WORKBENCH_EFFECT(X)
#endif
#if HS_ENABLE_CHAIN_INTERPRETER
#include "workbench/ShaderChain.h"
#define HS_CHAIN_INTERPRETER_EFFECT(X) X(ShaderChain)
#else
#define HS_CHAIN_INTERPRETER_EFFECT(X)
#endif
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
  X(Fishbowl)                                                                  \
  X(Comets)                                                                    \
  X(GridSpace)                                                                 \
  X(HyperLattice)                                                              \
  X(AshCloud)                                                                  \
  X(LatticeMelt)                                                               \
  X(DisplacementField)                                                         \
  X(DreamBalls)                                                                \
  X(Dynamo)                                                                    \
  X(KaleidoscopeFlowers)                                                       \
  X(KaleidoscopeSmooth)                                                        \
  X(KaleidoscopeMandala)                                                       \
  X(GnomonicStars)                                                             \
  X(GSReactionDiffusion)                                                       \
  X(AlienCore)                                                                 \
  X(HankinSolids)                                                              \
  X(HopfFibration)                                                             \
  X(IslamicStars)                                                              \
  X(KaleidoscopeHexBright)                                                     \
  X(AlienOcean)                                                                \
  X(KaleidoscopeHexSoft)                                                       \
  X(MeshFeedback)                                                              \
  X(MindSplatter)                                                              \
  X(MobiusGrid)                                                                \
  X(MobiusRings)                                                               \
  X(PetalFlow)                                                                 \
  X(KaleidoscopePentBright)                                                    \
  X(KaleidoscopeHexOil)                                                        \
  X(Raymarch)                                                                  \
  X(RingShower)                                                                \
  X(RingSpin)                                                                  \
  HS_SHADER_WORKBENCH_EFFECT(X)                                                \
  HS_CHAIN_INTERPRETER_EFFECT(X)                                               \
  X(ShapeShifter)                                                              \
  X(AlienBrain)                                                                \
  X(SphericalHarmonics)                                                        \
  X(CosmicEyeball)                                                             \
  X(Thrusters)                                                                 \
  X(KaleidoscopeStainedGlass)                                                  \
  X(Voronoi)

/// Phantasm renders one frame per half-revolution, so 480 RPM is 16 fps.
constexpr int HS_SHOW_FRAMES_PER_SECOND = 16;

/**
 * @brief Show duration that gives every preset one full dwell.
 * @tparam EffectT Effect class template whose preset cadence sizes the window.
 * @return Whole seconds covering every dwell and the intervening segues.
 */
template <template <int, int> class EffectT>
constexpr int hs_preset_window_seconds() {
  using Effect = EffectT<96, 20>;
  constexpr size_t PRESET_COUNT = Effect::PRESET_IDS.size();
  static_assert(PRESET_COUNT > 0);
  constexpr size_t FRAMES = PRESET_COUNT * Effect::PRESET_DWELL_FRAMES +
                            (PRESET_COUNT - 1) * Effect::PRESET_SEGUE.frames;
  return static_cast<int>((FRAMES + HS_SHOW_FRAMES_PER_SECOND - 1) /
                          HS_SHOW_FRAMES_PER_SECOND);
}

/**
 * @brief Phantasm playlist: HS_EFFECT_LIST minus Shader, ShaderChain, and the
 *        low-res-only effects Dynamo, MobiusRings, and Thrusters.
 * @param X Function-like macro applied to each effect type name and its show
 *          duration in seconds.
 * @details Entry order is the device show order, chosen independently of
 *   HS_EFFECT_LIST. Only the Phantasm firmware target consumes this; the
 *   registry, tests, and gallery stay on the full roster.
 *   The static_assert below HS_PHANTASM_EFFECT_COUNT forces this list to be
 *   revisited whenever HS_EFFECT_LIST gains or loses an entry.
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

/** Shader promotion product group in gallery and fixed-pipeline roster
 * order; the device show order is HS_PHANTASM_EFFECT_LIST's. */
#define HS_SHADER_PRODUCT_GROUP(X)                                             \
  X(AlienBrain, hs_preset_window_seconds<AlienBrain>())                        \
  X(KaleidoscopeHexSoft, hs_preset_window_seconds<KaleidoscopeHexSoft>())      \
  X(AlienOcean, hs_preset_window_seconds<AlienOcean>())                        \
  X(AlienCore, hs_preset_window_seconds<AlienCore>())                          \
  X(KaleidoscopeMandala, hs_preset_window_seconds<KaleidoscopeMandala>())      \
  X(GridSpace, hs_preset_window_seconds<GridSpace>())                          \
  X(LatticeMelt, hs_preset_window_seconds<LatticeMelt>())                      \
  X(AshCloud, hs_preset_window_seconds<AshCloud>())                            \
  X(KaleidoscopePentBright,                                                    \
    hs_preset_window_seconds<KaleidoscopePentBright>())                        \
  X(KaleidoscopeHexOil, hs_preset_window_seconds<KaleidoscopeHexOil>())        \
  X(KaleidoscopeStainedGlass,                                                  \
    hs_preset_window_seconds<KaleidoscopeStainedGlass>())                      \
  X(KaleidoscopeSmooth, hs_preset_window_seconds<KaleidoscopeSmooth>())        \
  X(KaleidoscopeHexBright, hs_preset_window_seconds<KaleidoscopeHexBright>())  \
  X(KaleidoscopeFlowers, hs_preset_window_seconds<KaleidoscopeFlowers>())      \
  X(CosmicEyeball, hs_preset_window_seconds<CosmicEyeball>())                  \
  X(MobiusGrid, hs_preset_window_seconds<MobiusGrid>())

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
