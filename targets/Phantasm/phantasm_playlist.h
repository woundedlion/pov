/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <iterator>

#include "targets/effects.h"

/**
 * @file phantasm_playlist.h
 * @brief Phantasm device playlist tables derived from phantasm_playlist.def,
 *        plus its roster drift guards.
 */

/** @brief Device playlist class names, in show order. */
inline constexpr const char *const HS_PHANTASM_NAMES[] = {
#define HS_PHANTASM_ROW(name, duration_seconds) #name,
#include "targets/Phantasm/phantasm_playlist_rows.h"
};

/** @brief Number of show rows in the playlist. */
constexpr int HS_PHANTASM_EFFECT_COUNT =
    static_cast<int>(std::size(HS_PHANTASM_NAMES));

/** @brief Device playlist seeds in show order. */
template <int W, int H>
inline constexpr uint64_t HS_PHANTASM_EFFECT_SEEDS[] = {
#define HS_PHANTASM_ROW(name, duration_seconds)                                \
  hs::stable_effect_seed(hs::effect_seed_key<name<W, H>>(#name)),
#include "targets/Phantasm/phantasm_playlist_rows.h"
};

/**
 * @brief Show duration the playlist assigns `name`, 0 when absent.
 * @param name Effect class name to look up.
 * @return Duration in seconds.
 */
constexpr int hs_phantasm_duration_seconds(const char *name) {
  int seconds = 0;
#define HS_PHANTASM_ROW(cls, duration_seconds)                                 \
  if (hs_effect_name_eq(name, #cls))                                           \
    seconds = (duration_seconds);
#include "targets/Phantasm/phantasm_playlist_rows.h"
  return seconds;
}

/** @brief Roster names the playlist excludes. Dynamo, MobiusRings and
 *  Thrusters are low-resolution only; other exclusions are product curation. */
inline constexpr const char *const HS_PHANTASM_EXCLUDED_NAMES[] = {
#define HS_PHANTASM_EXCLUDED_ROW(name) #name,
#include "targets/Phantasm/phantasm_playlist_rows.h"
};

/** @brief Number of exclude rows in the playlist. */
constexpr int HS_PHANTASM_EXCLUDED_COUNT =
    static_cast<int>(std::size(HS_PHANTASM_EXCLUDED_NAMES));

/**
 * @brief True when no excluded name is also shown.
 * @return True when every exclusion is absent from the show rows.
 */
constexpr bool hs_phantasm_exclusions_are_omitted() {
  for (const char *name : HS_PHANTASM_EXCLUDED_NAMES)
    if (hs_name_count(HS_PHANTASM_NAMES, name) != 0)
      return false;
  return true;
}

// Segmented playlist roster admission.
static_assert(hs_names_distinct(HS_PHANTASM_NAMES),
              "phantasm_playlist.def shows an effect twice — the duplicate "
              "is masking an omitted effect");
static_assert(hs_names_subset(HS_PHANTASM_NAMES, HS_EFFECT_NAMES),
              "phantasm_playlist.def shows an effect that is not in "
              "HS_EFFECT_LIST — a rename or typo left the playlist off-roster");
static_assert(hs_names_distinct(HS_PHANTASM_EXCLUDED_NAMES),
              "phantasm_playlist.def excludes an effect twice");
static_assert(HS_PHANTASM_EFFECT_COUNT == HS_EFFECT_COUNT -
                                              HS_PHANTASM_EXCLUDED_COUNT -
                                              HS_ENABLE_CHAIN_INTERPRETER,
              "phantasm_playlist.def out of sync with HS_EFFECT_LIST: every "
              "roster effect but the build-gated ShaderChain needs a show or "
              "an exclude row");
static_assert(!HS_ENABLE_CHAIN_INTERPRETER ||
                  hs_name_count(HS_EFFECT_NAMES, "ShaderChain") != 0,
              "HS_EFFECT_LIST must include HS_CHAIN_INTERPRETER_EFFECT when "
              "the chain interpreter is enabled");
static_assert(hs_names_subset(HS_PHANTASM_EXCLUDED_NAMES, HS_EFFECT_NAMES),
              "phantasm_playlist.def excludes a non-roster effect — a rename "
              "left the exclusion guard below vacuous");
static_assert(hs_name_count(HS_PHANTASM_NAMES, "ShaderChain") == 0 &&
                  hs_phantasm_exclusions_are_omitted(),
              "phantasm_playlist.def must not show ShaderChain or an excluded "
              "effect");
