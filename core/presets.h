/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <span>
#include <utility>

#include "platform.h" // HS_CHECK

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
 * @brief Fixed-size, cyclic selector over a set of Params presets.
 * @tparam Params The preset parameter type held by each entry.
 * @tparam Size The number of entries; must be greater than zero.
 * @details Tracks the current entry plus the one active before the last move,
 * so callers can crossfade between the outgoing and incoming presets.
 */
template <typename Params, size_t Size> class Presets {
  static_assert(Size > 0, "Presets requires at least one entry");

public:
  using Entry = PresetEntry<Params>;

  /**
   * @brief Builds a selector over the given entries, starting at index 0.
   * @param e The preset entries; copied into the container.
   * @details Non-explicit so brace-init (`Presets p = {{...}}`) and the CTAD
   * guide below keep working. The indices are deliberately NOT constructor
   * inputs — they start at 0 and only next()/prev() ever move them, which is
   * what keeps them in [0, Size) (see the private section).
   */
  constexpr Presets(std::array<Entry, Size> e) : entries(std::move(e)) {}

  /**
   * @brief Returns the params of the currently selected entry.
   * @return Const reference to the current entry's params.
   */
  const Params &get() const {
    HS_CHECK(current_idx >= 0 && static_cast<size_t>(current_idx) < Size,
             "Presets::get: current_idx out of range");
    return entries[current_idx].params;
  }

  /**
   * @brief Advances to the next entry, wrapping past the end.
   * @details Records the outgoing index in prev_idx before advancing.
   */
  void next() {
    prev_idx = current_idx;
    current_idx = (current_idx + 1) % Size;
  }

  /**
   * @brief Steps back to the previous entry, wrapping past the front.
   * @details Records the outgoing index in prev_idx before stepping back.
   */
  void prev() {
    prev_idx = current_idx;
    current_idx = (current_idx - 1 + Size) % Size;
  }

  /**
   * @brief Copies the current entry's params into target.
   * @param target Destination params, overwritten with the current entry.
   */
  void apply(Params &target) const { target = get(); }

  /**
   * @brief Returns a read-only view over all entries.
   * @return Span covering all Size entries in order.
   */
  std::span<const Entry> get_entries() const {
    return std::span<const Entry>(entries);
  }

  /**
   * @brief Returns a mutable view over all entries, for in-place rebinding.
   * @return Span covering all Size entries in order.
   */
  std::span<Entry> get_entries() { return std::span<Entry>(entries); }

  /** @brief Index of the currently selected entry; always in [0, Size). */
  int current_index() const { return current_idx; }
  /** @brief Index active before the last next()/prev(); always in [0, Size). */
  int prev_index() const { return prev_idx; }

  /**
   * @brief Returns the params of the entry active before the last move.
   * @return Const reference to the previous entry's params; for crossfades.
   */
  const Params &prev_get() const {
    HS_CHECK(prev_idx >= 0 && static_cast<size_t>(prev_idx) < Size,
             "Presets::prev_get: prev_idx out of range");
    return entries[prev_idx].params;
  }

private:
  std::array<Entry, Size> entries; /**< The backing store of preset entries. */
  int current_idx = 0; /**< Index of the currently selected entry. */
  int prev_idx = 0;    /**< Index active before the last next()/prev(); for crossfades. */
};

/**
 * @brief CTAD guide deducing Size from a brace-initialized entry array.
 * @tparam Params The preset parameter type held by each entry.
 * @tparam N The deduced number of entries.
 * @details Lets `Presets presets = {{...}}` deduce Size = N from the array.
 */
template <typename Params, size_t N>
Presets(std::array<PresetEntry<Params>, N>) -> Presets<Params, N>;
