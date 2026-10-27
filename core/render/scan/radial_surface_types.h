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
  float r = 0; ///< Red.
  float g = 0; ///< Green.
  float b = 0; ///< Blue.
};

/** @brief Work counts for a radial surface area pass. */
struct SurfaceAreaStats {
  uint64_t pixels = 0;            ///< Shaded pixels.
  uint64_t intervals = 0;         ///< Interval visits by all queries.
  uint64_t fragments = 0;         ///< Shaded pixel parts.
  uint64_t clips = 0;             ///< Pixels blended across two edges.
  int max_intervals = 0;          ///< Most interval visits for one pixel.
  int max_parts = 0;              ///< Most parts in one pixel.
  uint64_t edge_pixels = 0;       ///< Pixels given edge coverage.
  uint64_t point_intervals = 0;   ///< Interval visits by center queries.
  uint64_t outside_intervals = 0; ///< Interval visits by outside-color queries.
  int max_point_intervals = 0;    ///< Most center-query visits for one pixel.
  uint64_t mask_queries = 0;      ///< Candidate-mask lookups.
  uint64_t mask_hits = 0;         ///< Lookups with candidate intervals.
  uint64_t mask_prepares = 0;     ///< Candidate-mask preparations.
  uint64_t mask_candidates = 0;   ///< Canonical candidate segments prepared.
  uint64_t mask_pieces = 0;       ///< Padded candidate pieces.
  uint64_t mask_row_writes = 0;   ///< Candidate-cell updates.
};

} // namespace Scan
