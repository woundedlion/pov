/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_types.h
 * @brief Linear surface colors and optional rendering work counters. */

#include <cstdint>

namespace Scan {
/** @brief Area-weighted linear color. */
struct SurfaceColor {
  float r = 0, g = 0, b = 0;
};

/** @brief Work counts for a radial surface area pass. */
struct SurfaceAreaStats {
  uint64_t pixels = 0, intervals = 0, fragments = 0, clips = 0;
  int max_intervals = 0, max_parts = 0;
  uint64_t edge_pixels = 0, point_intervals = 0, outside_intervals = 0;
  int max_point_intervals = 0;
  uint64_t mask_queries = 0, mask_hits = 0, mask_prepares = 0,
           mask_candidates = 0, mask_pieces = 0, mask_row_writes = 0;
};

} // namespace Scan
