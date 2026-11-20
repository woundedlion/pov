/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file trail_store.h
 * @brief Shared storage and decay for frame-based trails. */

#include <algorithm>
#include <cstddef>
#include <span>
#include "core/memory.h"

namespace Filter::Detail {
/**
 * @brief Arena-backed trail slots with a positive remaining lifetime.
 * @tparam Item Sample record with a ttl field.
 * @tparam Capacity Maximum number of live samples.
 */
template <typename Item, int Capacity> class TrailStore {
public:
  static constexpr size_t STORAGE_BYTES =
      Capacity * sizeof(Item); ///< Arena bytes.

  /** @brief Allocates empty storage.
   * @param arena Arena providing Capacity slots.
   * @param owner Trail domain used in initialization diagnostics. */
  void init_storage(Arena &arena, [[maybe_unused]] const char *owner) {
#ifndef NDEBUG
    HS_CHECK(!slots || !stamp.block_alive(slots, STORAGE_BYTES),
             "%s filter: storage already initialized", owner);
#endif
    slots = arena.allocate_n<Item>(Capacity);
    count = 0;
#ifndef NDEBUG
    stamp.record(arena);
#endif
  }

  /** @brief Reports whether the store has been allocated. */
  explicit operator bool() const { return slots != nullptr; }

  /** @brief Checks allocation liveness under structural audits.
   * @param owner Trail type used in liveness diagnostics. */
  void check_storage_alive([[maybe_unused]] const char *owner) const {
#ifndef NDEBUG
    HS_AUDIT_CHECK(stamp.block_alive(slots, STORAGE_BYTES),
                   "%s use-after-free!", owner);
#endif
  }

  /** @brief Accesses live slots.
   * @return Live samples in slot order. Compaction can reorder them. */
  std::span<Item> live() { return {slots, count}; }

  /** @brief Counts live slots.
   * @return Number of buffered samples. */
  size_t size() const { return count; }

  /** @brief Appends a positive-TTL item, replacing the last at capacity.
   * @param item Trail sample. */
  void append(const Item &item) {
    if (count == Capacity)
      --count;
    slots[count++] = item;
  }

  /** @brief Ages each sample one frame and removes expired slots by swapping. */
  void age_and_compact() {
    for (size_t i = 0; i < count;) {
      Item &item = slots[i];
      if (--item.ttl <= 0)
        item = slots[--count];
      else
        ++i;
    }
  }

private:
  Item *slots = nullptr;
  size_t count = 0;
#ifndef NDEBUG
  ArenaBlockStamp stamp;
#endif
};

/** @brief Computes fade progress from a sample's remaining lifetime.
 * @param ttl Remaining frames.
 * @param lifetime Fade divisor in frames.
 * @return Progress in [0, 1]. */
HS_HOT_INLINE inline float trail_progress(float ttl, float lifetime) {
  return hs::clamp(1.0f - ttl / lifetime, 0.0f, 1.0f);
}

/** @brief Computes nonnegative age in the trail's frame unit.
 * @tparam T Trail frame register type.
 * @param ttl Remaining frames.
 * @param lifetime Fade divisor in frames.
 * @return Sample age in frames. */
template <typename T> HS_HOT_INLINE inline T trail_age(T ttl, T lifetime) {
  return std::max(T{0}, lifetime - ttl);
}
} // namespace Filter::Detail
