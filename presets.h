/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <initializer_list>
#include <string_view>

#include <array>
#include <span>
#include <utility>
#include <algorithm>

template <typename Params, size_t Size> class Presets {
public:
  struct Entry {
    const char *name;
    Params params;
  };
  // Aggregate Initialization
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
    if (Size == 0)
      return;
    current_idx = (current_idx + 1) % Size;
  }

  void prev() {
    if (Size == 0)
      return;
    current_idx = (current_idx - 1 + Size) % Size;
  }

  const char *get_current_name() const {
    if (Size == 0)
      return "";
    return entries[current_idx].name;
  }

  void apply(Params &target) const { target = get(); }

  std::span<const Entry> get_entries() const {
    return std::span<const Entry>(entries);
  }

  std::array<Entry, Size> entries;
  int current_idx = 0;
};
