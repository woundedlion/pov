/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Effect roster: pulls in every effect header plus the HS_EFFECT_LIST X-macro.
// Include this only from a build target (run_tests, the firmware entry point).
// An effect must NOT include it — that would pull in all the other effects and
// recurse. Effects include core/engine.h (the API umbrella) instead.

#include "effects/BZReactionDiffusion.h"
#include "effects/ChaoticStrings.h"
#include "effects/Comets.h"
#include "effects/DistortedRing.h"
#include "effects/DreamBalls.h"
#include "effects/Dynamo.h"
#include "effects/FlowField.h"
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
#include "effects/Moire.h"
#include "effects/PetalFlow.h"
#include "effects/Raymarch.h"
#include "effects/RingShower.h"
#include "effects/RingSpin.h"
#include "effects/ShapeShifter.h"
#include "effects/SphericalHarmonics.h"
#include "effects/SplineFlow.h"
#include "effects/Thrusters.h"
#include "effects/Voronoi.h"

/**
 * @brief Single source of truth for the registered effect roster, as an X-macro.
 * @param X Function-like macro applied to each effect type name in the roster.
 * @details The `#include` list above and this X-macro list must stay in lock-step.
 *     * WASM:   the self-registering EffectRegistry size is checked against
 *               HS_EFFECT_COUNT at engine startup (targets/wasm/wasm.cpp), so a
 *               registered-but-unlisted (or listed-but-unregistered) effect traps.
 *               This is the load-bearing anti-drift check.
 *     * Native: the effect smoke suite iterates this X-macro list, so its coverage
 *               is derived from the list rather than hand-maintained
 *               (tests/test_effects.h). Because the suite is driven BY the list, a
 *               forgotten X() row silently drops that effect from native coverage
 *               rather than failing the build.
 *   Adding an effect therefore means: add the `#include` above, the
 *   REGISTER_EFFECT in its header, and one X() row here.
 */
#define HS_EFFECT_LIST(X)                                                       \
  X(BZReactionDiffusion)                                                        \
  X(ChaoticStrings)                                                             \
  X(Comets)                                                                     \
  X(DistortedRing)                                                              \
  X(DreamBalls)                                                                 \
  X(Dynamo)                                                                     \
  X(FlowField)                                                                  \
  X(Flyby)                                                                      \
  X(GnomonicStars)                                                              \
  X(GSReactionDiffusion)                                                        \
  X(HankinSolids)                                                               \
  X(HopfFibration)                                                              \
  X(IslamicStars)                                                               \
  X(Liquid2D)                                                                   \
  X(MeshFeedback)                                                               \
  X(MindSplatter)                                                               \
  X(MobiusGrid)                                                                 \
  X(Moire)                                                                      \
  X(PetalFlow)                                                                  \
  X(Raymarch)                                                                   \
  X(RingShower)                                                                 \
  X(RingSpin)                                                                   \
  X(ShapeShifter)                                                               \
  X(SphericalHarmonics)                                                         \
  X(SplineFlow)                                                                 \
  X(Thrusters)                                                                  \
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
#undef HS_EFFECT_COUNT_ADD

