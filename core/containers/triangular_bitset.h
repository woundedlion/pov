/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file triangular_bitset.h
 * @brief Compact upper-triangular bitset for unordered pair membership.
 */

#include "platform/platform.h"
#include <climits>
#include <cstdint>
#include <cstring>

/**
 * @brief Upper-triangular bitset for O(1) pair deduplication.
 * @tparam MAX_V Maximum vertex/element index (exclusive).
 * @details Stores one bit per unique unordered pair (a, b) where a < b < MAX_V.
 * Total storage: MAX_V * (MAX_V - 1) / 2 bits.
 */
template <int MAX_V> struct TriangularBitset {
  // BITS (here) and index() below form an intermediate product ~MAX_V^2 in `int`;
  // for MAX_V >= ~46341 that overflows int32 and corrupts the bit layout. The
  // static_assert pins the ceiling so a future large-mesh instantiation fails at
  // compile time, not at runtime.
  static_assert(
      static_cast<long long>(MAX_V) * MAX_V <= INT_MAX,
      "TriangularBitset: MAX_V too large; index() product overflows int");
  // Below 2 there is no unordered pair, so BYTES is 0 and data[] is a
  // zero-length array (a GNU extension, not ISO C++).
  static_assert(MAX_V >= 2, "TriangularBitset: MAX_V must be at least 2");
  static constexpr int BITS = MAX_V * (MAX_V - 1) / 2;
  static constexpr int BYTES = (BITS + 7) / 8;
  uint8_t data[BYTES] = {}; /**< Packed bit storage; zero-initialized so a pair
                                 read before clear() reads "unset" rather than UB. */

  /**
   * @brief Clears every pair bit to zero.
   */
  void clear() { memset(data, 0, BYTES); }

  /**
   * @brief Bit index for pair (small, large) where small < large.
   * @param small Lower index of the pair, in [0, large).
   * @param large Higher index of the pair, in (small, MAX_V).
   * @return Bit index into the packed storage.
   */
  static int index(int small, int large) {
    // The triangular layout is only valid for an ordered, in-range pair: a
    // swapped pair aliases the wrong bit (dedup corruption) and an out-of-range
    // one writes adjacent memory. HS_CHECK (survives NDEBUG) fails fast; this runs
    // on the per-edge mesh-dedup setup path (plot.h draw()), not a per-pixel loop.
    HS_CHECK(small >= 0 && small < large && large < MAX_V,
             "TriangularBitset::index: pair (%d, %d) violates "
             "0 <= small < large < %d",
             small, large, MAX_V);
    return small * (2 * MAX_V - small - 1) / 2 + (large - small - 1);
  }

  /**
   * @brief Tests whether pair (a, b) is set.
   * @param a Lower index of the pair; requires a < b.
   * @param b Higher index of the pair (see index()).
   * @return True iff the pair's bit is set.
   */
  bool test(int a, int b) const {
    int bit = index(a, b);
    return (data[bit >> 3] >> (bit & 7)) & 1;
  }

  /**
   * @brief Tests and sets the bit for pair (a, b).
   * @param a Lower index of the pair; requires a < b.
   * @param b Higher index of the pair (see index()).
   * @return True if already set (hit), false if newly inserted (miss).
   */
  bool test_and_set(int a, int b) {
    int bit = index(a, b);
    uint8_t &byte = data[bit >> 3];
    uint8_t mask = 1 << (bit & 7);
    if (byte & mask)
      return true;
    byte |= mask;
    return false;
  }
};
