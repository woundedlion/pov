/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file srgb_decode.h
 * @brief Split-region linear16 to sRGB8 decode over the generated bucket
 *        tables.
 */

#include "color/srgb_decode_lut.h"
#include "engine/platform.h"
#include <array>
#include <cstdint>
#include <iterator>

// Bucket geometry of the split-decode, shared with the generator of record
// (scripts/generate_srgb_decode.cpp, which includes this header): the low
// region is 1<<LOW_SHIFT wide, the high region 1<<HIGH_SHIFT.
inline constexpr int SRGB_DECODE_LOW_SHIFT = 4;
inline constexpr int SRGB_DECODE_HIGH_SHIFT = 7;
inline constexpr int SRGB_DECODE_LOW_MASK = (1 << SRGB_DECODE_LOW_SHIFT) - 1;
inline constexpr int SRGB_DECODE_HIGH_MASK = (1 << SRGB_DECODE_HIGH_SHIFT) - 1;
inline constexpr int SRGB_DECODE_LOW_N =
    SRGB_DECODE_VSPLIT >> SRGB_DECODE_LOW_SHIFT;
inline constexpr int SRGB_DECODE_HIGH_N =
    (65536 - SRGB_DECODE_VSPLIT) >> SRGB_DECODE_HIGH_SHIFT;

// A committed table generated under different shifts would otherwise be copied
// out of bounds below.
static_assert(std::size(srgb_decode_low_src) == SRGB_DECODE_LOW_N,
              "srgb_decode_low_src length disagrees with the low-region shift");
static_assert(
    std::size(srgb_decode_high_src) == SRGB_DECODE_HIGH_N,
    "srgb_decode_high_src length disagrees with the high-region shift");

// DTCM copies (zero-wait, bypass the L1 D-cache) of the two-region decode
// tables. Residing in DTCM is what stops the concurrent render from evicting
// them, unlike a cacheable table.
// constinit is load-bearing: an inline variable's dynamic init is unordered
// against other translation units' static initializers, so a runtime fill would
// let linear_to_srgb8 read a zeroed table and return 0 for every input.
inline constinit std::array<uint16_t, SRGB_DECODE_LOW_N> srgb_decode_low =
    []() constexpr {
      std::array<uint16_t, SRGB_DECODE_LOW_N> t{};
      for (int i = 0; i < SRGB_DECODE_LOW_N; ++i)
        t[i] = srgb_decode_low_src[i];
      return t;
    }();
inline constinit std::array<uint16_t, SRGB_DECODE_HIGH_N> srgb_decode_high =
    []() constexpr {
      std::array<uint16_t, SRGB_DECODE_HIGH_N> t{};
      for (int i = 0; i < SRGB_DECODE_HIGH_N; ++i)
        t[i] = srgb_decode_high_src[i];
      return t;
    }();

/**
 * @brief Bit-exact linear-16 -> sRGB-8 encode via a two-region split-decode.
 * @param v Linear 16-bit channel value.
 * @return sRGB 8-bit output, identical to linear_to_srgb_lut[v] for all v.
 * @details Replaces the 64 KB linear_to_srgb_lut (which exceeds the 32 KB L1
 * D-cache and thrashes against the framebuffer read on the pack hot path) with
 * ~1.5 KB of DTCM tables. A fine 16-wide low region and a coarse 128-wide high
 * region each hold at most one output step, so each side is a single branchless
 * compare: base + (frac >= step). Equivalence over all 65536 inputs is checked
 * by unit_color's test_linear_to_srgb8_decode_matches_lut.
 */
inline __attribute__((always_inline)) uint8_t linear_to_srgb8(uint16_t v) {
  if (v < SRGB_DECODE_VSPLIT) {
    uint16_t e = srgb_decode_low[v >> SRGB_DECODE_LOW_SHIFT];
    return (uint8_t)((e & 0xFF) +
                     ((v & SRGB_DECODE_LOW_MASK) >= (e >> 8) ? 1 : 0));
  }
  uint32_t d = (uint32_t)v - SRGB_DECODE_VSPLIT;
  uint16_t e = srgb_decode_high[d >> SRGB_DECODE_HIGH_SHIFT];
  return (uint8_t)((e & 0xFF) +
                   ((d & SRGB_DECODE_HIGH_MASK) >= (e >> 8) ? 1 : 0));
}
