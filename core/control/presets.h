/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file presets.h
 * @brief PresetEntry and all_presets_in_ranges: the preset-table row type and
 *        its compile-time range check. ChoreographedEffect
 *        (control/choreography.h) owns the tables built from these.
 */

#include <array>
#include <cstddef>

/**
 * @brief Standalone entry wrapping a single Params preset.
 * @tparam Params The preset parameter type stored in each entry.
 * @details A dedicated entry type lets the CTAD deduction guides avoid
 * dependent-name issues when deducing the preset count.
 */
template <typename Params> struct PresetEntry {
  Params params; /**< The stored preset parameters. */
};

/**
 * @brief True iff every entry of @p entries satisfies @p in_ranges.
 * @tparam Params The preset parameter type stored in each entry.
 * @tparam N The number of entries.
 * @tparam Predicate Callable accepting one const Params reference.
 * @param entries The preset table to check.
 * @param in_ranges Predicate testing one preset against its slider ranges.
 * @return True when every entry passes.
 * @details For a `static_assert` over a preset table: the loop covers appended
 * entries, which an unrolled conjunction over literal indices does not.
 */
template <typename Params, size_t N, typename Predicate>
constexpr bool
all_presets_in_ranges(const std::array<PresetEntry<Params>, N> &entries,
                      Predicate &&in_ranges) {
  for (const auto &e : entries)
    if (!in_ranges(e.params))
      return false;
  return true;
}
