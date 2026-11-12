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
/**
 * @def HS_CHAIN_INTERPRETER_EFFECT(X)
 * @brief Applies @p X to ShaderChain when the chain interpreter is built;
 *        empty otherwise.
 * @param X Function-like macro applied to the effect type name.
 */
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
  X(Millefiori)                                                                \
  X(GridSpace)                                                                 \
  X(Involution)                                                                \
  X(GridWorms)                                                                 \
  X(GridWarp)                                                                  \
  X(LiquidWindows)                                                             \
  X(Iridescence)                                                               \
  X(HyperLattice)                                                              \
  X(AshCloud)                                                                  \
  X(KaleidoscopeMelt)                                                          \
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
  X(MeshFeedback)                                                              \
  X(MindSplatter)                                                              \
  X(MobiusGrid)                                                                \
  X(MobiusRings)                                                               \
  X(PetalFlow)                                                                 \
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
 * @brief Rounds a frame count up to whole show seconds.
 * @param frames Duration in frames.
 * @return Whole seconds covering the frames.
 */
constexpr int hs_show_seconds(size_t frames) {
  return static_cast<int>((frames + HS_SHOW_FRAMES_PER_SECOND - 1) /
                          HS_SHOW_FRAMES_PER_SECOND);
}

/**
 * @brief Show duration covering an effect's full cycle.
 * @tparam EffectT Effect class template exposing FULL_CYCLE_FRAMES.
 * @return Whole seconds covering the cycle.
 */
template <template <int, int> class EffectT>
constexpr int hs_full_cycle_seconds() {
  return hs_show_seconds(
      EffectT<HS_PRESET_WINDOW_W, HS_PRESET_WINDOW_H>::FULL_CYCLE_FRAMES);
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
  return hs_show_seconds(FRAMES);
}

/** Shader promotion product group; the device show order is
 * HS_PHANTASM_EFFECT_LIST's. */
#define HS_SHADER_PRODUCT_GROUP(X)                                             \
  X(AlienBrain)                                                                \
  X(AlienOcean)                                                                \
  X(AlienCore)                                                                 \
  X(KaleidoscopeMandala)                                                       \
  X(Millefiori)                                                                \
  X(GridSpace)                                                                 \
  X(Involution)                                                                \
  X(GridWorms)                                                                 \
  X(GridWarp)                                                                  \
  X(LiquidWindows)                                                             \
  X(Iridescence)                                                               \
  X(KaleidoscopeMelt)                                                          \
  X(ChromaticLichen)                                                           \
  X(MermaidSkin)                                                               \
  X(JewelMelt)                                                                 \
  X(AshCloud)                                                                  \
  X(KaleidoscopeHexOil)                                                        \
  X(KaleidoscopeStainedGlass)                                                  \
  X(KaleidoscopeSmooth)                                                        \
  X(KaleidoscopeHexBright)                                                     \
  X(KaleidoscopeFlowers)                                                       \
  X(CosmicEyeball)                                                             \
  X(MobiusGrid)

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

#define HS_EFFECT_NAME_ENTRY(cls) #cls,
/// Engine roster class names.
inline constexpr const char *const HS_EFFECT_NAMES[] = {
    HS_EFFECT_LIST(HS_EFFECT_NAME_ENTRY)};
/// Shader product class names.
inline constexpr const char *const HS_SHADER_PRODUCT_NAMES[] = {
    HS_SHADER_PRODUCT_GROUP(HS_EFFECT_NAME_ENTRY)};
#undef HS_EFFECT_NAME_ENTRY

/**
 * @brief Counts an effect name in a name array.
 * @tparam N Array length.
 * @param names Names to search.
 * @param name Effect class name.
 * @return Number of matching entries.
 */
template <size_t N>
constexpr int hs_name_count(const char *const (&names)[N], const char *name) {
  int count = 0;
  for (const char *entry : names)
    count += hs_effect_name_eq(entry, name);
  return count;
}

/**
 * @brief Tests whether all names are distinct.
 * @tparam N Array length.
 * @param names Names to test.
 * @return True when no name appears twice.
 */
template <size_t N>
constexpr bool hs_names_distinct(const char *const (&names)[N]) {
  for (const char *name : names)
    if (hs_name_count(names, name) != 1)
      return false;
  return true;
}

/**
 * @brief Tests whether every name belongs to another array.
 * @tparam N Candidate array length.
 * @tparam M Superset array length.
 * @param names Candidate subset.
 * @param superset Names admitting the subset.
 * @return True when every candidate is admitted.
 */
template <size_t N, size_t M>
constexpr bool hs_names_subset(const char *const (&names)[N],
                               const char *const (&superset)[M]) {
  for (const char *name : names)
    if (hs_name_count(superset, name) == 0)
      return false;
  return true;
}

static_assert(hs_names_distinct(HS_SHADER_PRODUCT_NAMES),
              "HS_SHADER_PRODUCT_GROUP names an effect twice");
static_assert(hs_names_subset(HS_SHADER_PRODUCT_NAMES, HS_EFFECT_NAMES),
              "HS_SHADER_PRODUCT_GROUP names an effect that is not in "
              "HS_EFFECT_LIST — a rename or typo left the group off-roster");
