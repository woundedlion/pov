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
#include "effects/DisplacementField.h"
#include "effects/DreamBalls.h"
#include "effects/Dynamo.h"
#include "effects/Flyby.h"
#include "effects/GnomonicStars.h"
#include "effects/GSReactionDiffusion.h"
#include "effects/HankinSolids.h"
#include "effects/HopfFibration.h"
#include "effects/IslamicStars.h"
#include "effects/Liquid2D.h"
#include "effects/MeshFeedback.h"
#include "effects/MindSplatter.h"
#include "effects/MobiusGrid.h"
#include "effects/PetalFlow.h"
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
 *               HS_EFFECT_COUNT at engine startup (targets/wasm/wasm.cpp), so a
 *               registered-but-unlisted (or listed-but-unregistered) effect traps.
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
  X(DisplacementField)                                                         \
  X(DreamBalls)                                                                \
  X(Dynamo)                                                                    \
  X(Flyby)                                                                     \
  X(GnomonicStars)                                                             \
  X(GSReactionDiffusion)                                                       \
  X(HankinSolids)                                                              \
  X(HopfFibration)                                                             \
  X(IslamicStars)                                                              \
  X(Liquid2D)                                                                  \
  X(MeshFeedback)                                                              \
  X(MindSplatter)                                                              \
  X(MobiusGrid)                                                                \
  X(PetalFlow)                                                                 \
  X(Raymarch)                                                                  \
  X(RingShower)                                                                \
  X(RingSpin)                                                                  \
  X(ShapeShifter)                                                              \
  X(SphericalHarmonics)                                                        \
  X(Thrusters)                                                                 \
  X(Voronoi)

/**
 * @brief Phantasm (288x144) playlist: HS_EFFECT_LIST minus the low-res-only
 *        effects (Dynamo, Thrusters — Holosphere 96x20 only).
 * @param X Function-like macro applied to each effect type name in the playlist.
 * @details Same order as HS_EFFECT_LIST. Only the Phantasm firmware target
 *   consumes this; the registry, tests, and gallery stay on the full roster.
 *   The static_assert below HS_PHANTASM_EFFECT_COUNT forces this list to be
 *   revisited whenever HS_EFFECT_LIST gains or loses an entry.
 */
#define HS_PHANTASM_EFFECT_LIST(X)                                             \
  X(BZReactionDiffusion)                                                       \
  X(ChaoticStrings)                                                            \
  X(Comets)                                                                    \
  X(DisplacementField)                                                         \
  X(DreamBalls)                                                                \
  X(Flyby)                                                                     \
  X(GnomonicStars)                                                             \
  X(GSReactionDiffusion)                                                       \
  X(HankinSolids)                                                              \
  X(HopfFibration)                                                             \
  X(IslamicStars)                                                              \
  X(Liquid2D)                                                                  \
  X(MeshFeedback)                                                              \
  X(MindSplatter)                                                              \
  X(MobiusGrid)                                                                \
  X(PetalFlow)                                                                 \
  X(Raymarch)                                                                  \
  X(RingShower)                                                                \
  X(RingSpin)                                                                  \
  X(ShapeShifter)                                                              \
  X(SphericalHarmonics)                                                        \
  X(Voronoi)

/**
 * @brief Expands to +1 so HS_EFFECT_LIST can be summed into an entry count.
 * @param name Effect type name supplied by HS_EFFECT_LIST (unused).
 */
#define HS_EFFECT_COUNT_ADD(name) +1
/**
 * @brief Number of entries in HS_EFFECT_LIST, derived rather than hand-counted.
 */
constexpr int HS_EFFECT_COUNT = 0 HS_EFFECT_LIST(HS_EFFECT_COUNT_ADD);
/**
 * @brief Number of entries in HS_PHANTASM_EFFECT_LIST, derived rather than
 *        hand-counted.
 */
constexpr int HS_PHANTASM_EFFECT_COUNT =
    0 HS_PHANTASM_EFFECT_LIST(HS_EFFECT_COUNT_ADD);
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
#define HS_PHANTASM_NAME_MATCH(cls) || hs_effect_name_eq(name, #cls)
  return false HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_MATCH);
#undef HS_PHANTASM_NAME_MATCH
}

/**
 * @brief Occurrences of `name` in HS_PHANTASM_EFFECT_LIST.
 * @param name Effect class name to count.
 */
constexpr int hs_phantasm_effect_list_count(const char *name) {
#define HS_PHANTASM_NAME_COUNT(cls) +(hs_effect_name_eq(name, #cls) ? 1 : 0)
  return 0 HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_COUNT);
#undef HS_PHANTASM_NAME_COUNT
}

/** @brief True when no HS_PHANTASM_EFFECT_LIST name appears twice. */
constexpr bool hs_phantasm_effect_list_is_distinct() {
#define HS_PHANTASM_NAME_ONCE(cls) &&(hs_phantasm_effect_list_count(#cls) == 1)
  return true HS_PHANTASM_EFFECT_LIST(HS_PHANTASM_NAME_ONCE);
#undef HS_PHANTASM_NAME_ONCE
}

// Drift guard: an effect added to (or removed from) HS_EFFECT_LIST must also be
// deliberately added to or excluded from the Phantasm playlist above. The count
// pins the cardinality; the name scans pin *which* two are missing, so swapping
// one exclusion for another cannot ride green. Distinctness closes the last
// hole: a duplicated entry paired with an omission holds the count and passes
// both exclusion scans while an effect silently drops off the playlist.
static_assert(hs_phantasm_effect_list_is_distinct(),
              "HS_PHANTASM_EFFECT_LIST names an effect twice — the duplicate "
              "is masking an omitted effect");
static_assert(HS_PHANTASM_EFFECT_COUNT == HS_EFFECT_COUNT - 2,
              "HS_PHANTASM_EFFECT_LIST out of sync with HS_EFFECT_LIST "
              "(full roster minus Dynamo and Thrusters)");
static_assert(hs_in_effect_list("Dynamo") && hs_in_effect_list("Thrusters"),
              "Phantasm exclusion names a non-roster effect — a rename left "
              "the exclusion guard below vacuous");
static_assert(!hs_in_phantasm_effect_list("Dynamo") &&
                  !hs_in_phantasm_effect_list("Thrusters"),
              "HS_PHANTASM_EFFECT_LIST must exclude Dynamo and Thrusters "
              "(Holosphere 96x20 only)");
