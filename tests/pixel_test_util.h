/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared framebuffer predicates for the suites that assert on a StubEffect's
 * pixels — canvas, filter, scan, and mesh raster.
 *
 * The two lit-pixel counters differ in the region they scan and are not
 * interchangeable: count_lit_canvas() takes the effect's own reported
 * dimensions, count_lit_region() a compile-time window that a test states
 * independently of the effect it built.
 */
#pragma once

#include "core/render/canvas.h"
#include "tests/test_fixture.h"

#include <cstddef>
#include <cstdint>

namespace hs_test {

/**
 * @brief Tests whether a pixel is fully black (cleared-frame / unwritten).
 * @param p Pixel to inspect.
 * @return True when all of the pixel's RGB channels are zero.
 */
inline bool is_black(const Pixel &p) {
  return p.r == 0 && p.g == 0 && p.b == 0;
}

/**
 * @brief Tests whether a pixel has exactly the given channel values.
 * @param p Pixel to inspect.
 * @param r Expected red channel value.
 * @param g Expected green channel value.
 * @param b Expected blue channel value.
 * @return True if all three channels match exactly.
 */
inline bool pix_eq(const Pixel &p, uint16_t r, uint16_t g, uint16_t b) {
  return p.r == r && p.g == g && p.b == b;
}

/**
 * @brief Counts the non-black pixels across the effect's reported canvas.
 * @param fx Effect whose framebuffer is scanned, fx.width() by fx.height().
 * @return Number of lit (non-black) pixels.
 */
inline size_t count_lit_canvas(const StubEffect &fx) {
  size_t n = 0;
  for (int y = 0; y < fx.height(); ++y)
    for (int x = 0; x < fx.width(); ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++n;
  return n;
}

/**
 * @brief Counts the non-black pixels across a fixed W by H canvas region.
 * @tparam W Region width in pixels.
 * @tparam H Region height in pixels.
 * @param fx Effect whose framebuffer is scanned.
 * @return Number of lit (non-black) pixels in the region.
 */
template <int W, int H> inline size_t count_lit_region(const StubEffect &fx) {
  size_t n = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++n;
  return n;
}

} // namespace hs_test
