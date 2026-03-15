/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <string_view>

#include <array>
#include <span>

// Standalone entry type so CTAD deduction guides avoid dependent-name issues.
template <typename Params> struct PresetEntry {
  const char *name;
  Params params;
};

template <typename Params, size_t Size> class Presets {
public:
  using Entry = PresetEntry<Params>;

  const Params &get(const char *name) const {
    std::string_view target(name);
    for (const auto &entry : entries) {
      if (std::string_view(entry.name) == target) {
        return entry.params;
      }
    }
    return entries[0].params;
  }

  const Params &get() const { return entries[current_idx].params; }

  void next() {
    if constexpr (Size == 0)
      return;
    prev_idx = current_idx;
    current_idx = (current_idx + 1) % Size;
  }

  void prev() {
    if constexpr (Size == 0)
      return;
    current_idx = (current_idx - 1 + Size) % Size;
  }

  const char *get_current_name() const {
    if constexpr (Size == 0)
      return "";
    return entries[current_idx].name;
  }

  void apply(Params &target) const { target = get(); }

  std::span<const Entry> get_entries() const {
    return std::span<const Entry>(entries);
  }

  std::array<Entry, Size> entries;
  int current_idx = 0;
  int prev_idx = 0;

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
