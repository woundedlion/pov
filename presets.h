/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <initializer_list>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <algorithm>

template <typename Params> class Presets {
public:
  struct Entry {
    std::string name;
    Params params;
  };

  Presets(std::initializer_list<std::pair<const char *, Params>> init)
      : current_idx(0) {
    for (const auto &item : init) {
      entries.push_back({item.first, item.second});
    }
  }

  const Params &get(const char *name) const {
    for (const auto &entry : entries) {
      if (entry.name == name) {
        return entry.params;
      }
    }
    if (!entries.empty())
      return entries[0].params;
    throw std::runtime_error("No presets available");
  }

  const Params &get() const {
    if (entries.empty())
      throw std::runtime_error("No presets available");
    return entries[current_idx].params;
  }

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
    return entries[current_idx].name.c_str();
  }

  void apply(Params &target) const { target = get(); }

  const std::vector<Entry> &get_entries() const { return entries; }

private:
  std::vector<Entry> entries;
  int current_idx;
};
