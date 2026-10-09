/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"
#include "core/platform/constants.h" // MAX_W/MAX_H

/**
 * @file effects.h
 * @brief Effect roster: pulls in every effect header plus the HS_EFFECT_LIST
 *        X-macro.
 *
 * Roster effects include core/engine/engine.h instead; including this header
 * from one recurses.
 */

#include "targets/composed_effect_includes.h"
#include "effects/BZReactionDiffusion.h"
#include "effects/Fishbowl.h"
#include "effects/Galaxies.h"
#include "effects/Comets.h"
#include "effects/DisplacementField.h"
#include "effects/DreamBalls.h"
#include "effects/Dynamo.h"
#include "effects/HyperLattice.h"
#include "effects/GnomonicStars.h"
#include "effects/GSReactionDiffusion.h"
#include "effects/HankinSolids.h"
#include "effects/HopfFibration.h"
#include "effects/IslamicStars.h"
#include "effects/MeshFeedback.h"
#include "effects/MindSplatter.h"
#include "effects/MobiusRings.h"
#include "effects/PetalFlow.h"
#if HS_ENABLE_CHAIN_INTERPRETER
#include "workbench/shader/chain_host.h"
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
 * @details Composed includes are generated; other effect includes and all X()
 * rows are authored here.
 */
#define HS_EFFECT_LIST(X)                                                      \
  X(BZReactionDiffusion)                                                       \
  X(Fishbowl)                                                                  \
  X(Comets)                                                                    \
  X(GridSpace)                                                                 \
  X(HyperLattice)                                                              \
  X(AshCloud)                                                                  \
  X(LatticeMelt)                                                               \
  X(ChromaticLichen)                                                           \
  X(MermaidSkin)                                                               \
  X(JewelMelt)                                                                 \
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
  HS_CHAIN_INTERPRETER_EFFECT(X)                                               \
  X(ShapeShifter)                                                              \
  X(AlienBrain)                                                                \
  X(SphericalHarmonics)                                                        \
  X(CosmicEyeball)                                                             \
  X(Thrusters)                                                                 \
  X(KaleidoscopeStainedGlass)                                                  \
  X(Voronoi)                                                                   \
  X(Galaxies)

/// Phantasm renders one frame per half-revolution, so 480 RPM is 16 fps.
constexpr int HS_SHOW_FRAMES_PER_SECOND = 16;

/** Canvas the preset cadence is read at, so targets with different
 * CANVAS_W/CANVAS_H agree on a show length. */
inline constexpr int HS_PRESET_WINDOW_W = 96;
/** Height half of the cadence-reading canvas (see HS_PRESET_WINDOW_W). */
inline constexpr int HS_PRESET_WINDOW_H = 20;

/**
 * @brief Frames covering every preset dwell and the intervening segues.
 * @tparam EffectT Effect class template whose preset cadence sizes the window.
 * @tparam W Canvas width the cadence constants are read at.
 * @tparam H Canvas height the cadence constants are read at.
 * @return Whole frames.
 */
template <template <int, int> class EffectT, int W, int H>
constexpr size_t hs_preset_window_frames() {
  using EffectType = EffectT<W, H>;
  static_assert(EffectType::authored_preset_count() > 0);
  return EffectType::cycle_frames();
}

/**
 * @brief Show duration that gives every preset one full dwell.
 * @tparam EffectT Effect class template whose preset cadence sizes the window.
 * @return Whole seconds covering every dwell and the intervening segues.
 */
template <template <int, int> class EffectT>
constexpr int hs_preset_window_seconds() {
  constexpr size_t FRAMES = hs_preset_window_frames<EffectT, HS_PRESET_WINDOW_W,
                                                    HS_PRESET_WINDOW_H>();
  static_assert(FRAMES == (hs_preset_window_frames<EffectT, MAX_W, MAX_H>()),
                "preset cadence varies with canvas size, so a show duration "
                "read at one resolution no longer covers the rendered one");
  return static_cast<int>((FRAMES + HS_SHOW_FRAMES_PER_SECOND - 1) /
                          HS_SHOW_FRAMES_PER_SECOND);
}

/** Shader promotion product group in gallery and composed-effect roster
 * order; the device show order is HS_PHANTASM_EFFECT_LIST's. */
#define HS_SHADER_PRODUCT_GROUP(X)                                             \
  X(AlienBrain, hs_preset_window_seconds<AlienBrain>())                        \
  X(KaleidoscopeHexSoft, hs_preset_window_seconds<KaleidoscopeHexSoft>())      \
  X(AlienOcean, hs_preset_window_seconds<AlienOcean>())                        \
  X(AlienCore, hs_preset_window_seconds<AlienCore>())                          \
  X(KaleidoscopeMandala, hs_preset_window_seconds<KaleidoscopeMandala>())      \
  X(GridSpace, hs_preset_window_seconds<GridSpace>())                          \
  X(LatticeMelt, hs_preset_window_seconds<LatticeMelt>())                      \
  X(ChromaticLichen, hs_preset_window_seconds<ChromaticLichen>())              \
  X(MermaidSkin, hs_preset_window_seconds<MermaidSkin>())                      \
  X(JewelMelt, hs_preset_window_seconds<JewelMelt>())                          \
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
 * @brief Number of entries in HS_EFFECT_LIST.
 */
constexpr int HS_EFFECT_COUNT = 0 HS_EFFECT_LIST(HS_EFFECT_COUNT_ADD);
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

/** @brief Occurrences of `name` in HS_SHADER_PRODUCT_GROUP. */
constexpr int hs_shader_product_group_count(const char *name) {
#define HS_SHADER_PRODUCT_NAME_COUNT(cls, duration_seconds)                    \
  +(hs_effect_name_eq(name, #cls) ? 1 : 0)
  return 0 HS_SHADER_PRODUCT_GROUP(HS_SHADER_PRODUCT_NAME_COUNT);
#undef HS_SHADER_PRODUCT_NAME_COUNT
}

/** @brief True when no HS_SHADER_PRODUCT_GROUP name appears twice. */
constexpr bool hs_shader_product_group_is_distinct() {
#define HS_SHADER_PRODUCT_NAME_ONCE(cls, duration_seconds)                     \
  &&(hs_shader_product_group_count(#cls) == 1)
  return true HS_SHADER_PRODUCT_GROUP(HS_SHADER_PRODUCT_NAME_ONCE);
#undef HS_SHADER_PRODUCT_NAME_ONCE
}

/** @brief True when every HS_SHADER_PRODUCT_GROUP name is in HS_EFFECT_LIST. */
constexpr bool hs_shader_product_group_is_subset() {
#define HS_SHADER_PRODUCT_NAME_ON_ROSTER(cls, duration_seconds)                \
  &&hs_in_effect_list(#cls)
  return true HS_SHADER_PRODUCT_GROUP(HS_SHADER_PRODUCT_NAME_ON_ROSTER);
#undef HS_SHADER_PRODUCT_NAME_ON_ROSTER
}

static_assert(hs_shader_product_group_is_distinct(),
              "HS_SHADER_PRODUCT_GROUP names an effect twice");
static_assert(hs_shader_product_group_is_subset(),
              "HS_SHADER_PRODUCT_GROUP names an effect that is not in "
              "HS_EFFECT_LIST — a rename or typo left the group off-roster");
