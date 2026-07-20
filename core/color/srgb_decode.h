#pragma once
#include "color/srgb_decode_lut.h"
#include <cstdint>

/**
 * @brief Bit-exact linear-16 -> sRGB-8 encode via an 8 KB split-decode.
 * @param v Linear 16-bit channel value.
 * @return sRGB 8-bit output, identical to linear_to_srgb_lut[v] for all v.
 * @details Replaces the 64 KB linear_to_srgb_lut (which exceeds the 32 KB L1
 * D-cache and thrashes against the framebuffer read on the pack hot path) with
 * two 4 KB fields packed into srgb_decode_lut. 16-wide buckets each span at
 * most one output step, so the decode is branchless: base + (frac >= step).
 * Equivalence over all 65536 inputs is checked by unit_srgb_decode.
 */
#ifdef HS_PACK_DECODE_DTCM
// Diagnostic: mirror the decode table into DTCM (zero-wait, bypasses L1) to test
// whether the residual pack cost is the concurrent render evicting the
// flash-resident 8 KB table. g_srgb_decode_dtcm points into the unused top of
// the arena (DTCM); the mirror is copied lazily on first use.
extern uint16_t *g_srgb_decode_dtcm; // filled at static init in memory.cpp
inline __attribute__((always_inline)) uint8_t linear_to_srgb8(uint16_t v) {
  uint16_t e = g_srgb_decode_dtcm[v >> 4];
  return (uint8_t)((e & 0xFF) + ((v & 0xF) >= (e >> 8) ? 1 : 0));
}
#else
inline __attribute__((always_inline)) uint8_t linear_to_srgb8(uint16_t v) {
  uint16_t e = srgb_decode_lut[v >> 4];
  return (uint8_t)((e & 0xFF) + ((v & 0xF) >= (e >> 8) ? 1 : 0));
}
#endif
