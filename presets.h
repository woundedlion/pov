/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <initializer_list>
#include <string_view>
#include <vector>
#include <utility>
#include <algorithm>

template <typename Params> class Presets {
public:
  struct Entry {
    const char *name;
    Params params;
  };

  Presets(std::initializer_list<std::pair<const char *, Params>> init)
      : current_idx(0) {
    for (const auto &item : init) {
      entries.push_back({item.first, item.second});
    }
  }

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
    if (entries.empty())
      return;
    current_idx = (current_idx + 1) % entries.size();
  }

  void prev() {
    if (entries.empty())
      return;
    current_idx = (current_idx - 1 + entries.size()) % entries.size();
  }

  const char *get_current_name() const {
    if (entries.empty())
      return "";
    return entries[current_idx].name;
  }

  void apply(Params &target) const { target = get(); }

  const std::vector<Entry> &get_entries() const { return entries; }

private:
  std::vector<Entry> entries;
  int current_idx;
};
