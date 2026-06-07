/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

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
#include "effects/Metaballs.h"
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

// ---------------------------------------------------------------------------
// Single source of truth for the registered effect roster.
//
// The #include list above and this X-macro list must stay in lock-step. The
// coupling is enforced from both build modes so the roster can never silently
// drift from what actually ships or gets tested:
//   * WASM:   the self-registering EffectRegistry size is checked against
//             HS_EFFECT_COUNT at engine startup (targets/wasm/wasm.cpp).
//   * Native: the effect smoke suite generates exactly one case per entry, so
//             its coverage is derived from this list, not hand-maintained
//             (tests/test_effects.h).
// Adding an effect therefore means: add the #include above, the REGISTER_EFFECT
// in its header, and one X() row here — miss any and a build/startup check fires.
// ---------------------------------------------------------------------------
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
  X(Metaballs)                                                                  \
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

// Count of entries in HS_EFFECT_LIST (derived, never hand-counted).
#define HS_EFFECT_COUNT_ADD(name) +1
constexpr int HS_EFFECT_COUNT = 0 HS_EFFECT_LIST(HS_EFFECT_COUNT_ADD);
#undef HS_EFFECT_COUNT_ADD

