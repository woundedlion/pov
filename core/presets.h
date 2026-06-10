/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <span>

// Standalone entry type so CTAD deduction guides avoid dependent-name issues.
template <typename Params> struct PresetEntry {
  Params params;
};

template <typename Params, size_t Size> class Presets {
  // An empty preset set is meaningless for a selector and would make get()/
  // apply()/prev_get() index a zero-length array (UB). Reject it up front so
  // every accessor has at least one valid entry.
  static_assert(Size > 0, "Presets requires at least one entry");

public:
  using Entry = PresetEntry<Params>;

  const Params &get() const { return entries[current_idx].params; }

  void next() {
    prev_idx = current_idx;
    current_idx = (current_idx + 1) % Size;
  }

  void prev() {
    prev_idx = current_idx;
    current_idx = (current_idx - 1 + Size) % Size;
  }

  void apply(Params &target) const { target = get(); }

  std::span<const Entry> get_entries() const {
    return std::span<const Entry>(entries);
  }

  std::array<Entry, Size> entries;
  int current_idx = 0;
  int prev_idx = 0; // index active before the last next()/prev(); for crossfades.

  const Params &prev_get() const { return entries[prev_idx].params; }
};

// CTAD: deduce Size from the number of entries in the initializer.
// Presets presets = {{...}}  →  Size = N, deduced from std::array.
template <typename Params, size_t N>
Presets(std::array<PresetEntry<Params>, N>) -> Presets<Params, N>;

template <typename Params, size_t N>
Presets(std::array<PresetEntry<Params>, N>, int) -> Presets<Params, N>;

template <typename Params, size_t N>
Presets(std::array<PresetEntry<Params>, N>, int, int) -> Presets<Params, N>;
