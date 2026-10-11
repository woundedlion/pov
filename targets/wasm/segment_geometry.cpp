/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file segment_geometry.cpp
 * @brief Pure segment geometry exports compiled to standalone browser JavaScript.
 */
#include "hardware/pov_segment_map.h"

extern "C" {
/**
 * @brief Inclusive segment left column for an admitted layout.
 * @param id Segment index.
 * @param total Total segment count.
 * @param width Canvas columns.
 * @param height Canvas rows.
 * @return Partition boundary in canvas coordinates.
 */
double segment_x0(double id, double total, double width, double height) {
  return pov::segment_range(id, total, width, height).x0;
}
/**
 * @brief Exclusive segment right column for an admitted layout.
 * @param id Segment index.
 * @param total Total segment count.
 * @param width Canvas columns.
 * @param height Canvas rows.
 * @return Partition boundary in canvas coordinates.
 */
double segment_x1(double id, double total, double width, double height) {
  return pov::segment_range(id, total, width, height).x1;
}
/**
 * @brief Inclusive segment top row for an admitted layout.
 * @param id Segment index.
 * @param total Total segment count.
 * @param width Canvas columns.
 * @param height Canvas rows.
 * @return Partition boundary in canvas coordinates.
 */
double segment_y0(double id, double total, double width, double height) {
  return pov::segment_range(id, total, width, height).y0;
}
/**
 * @brief Exclusive segment bottom row for an admitted layout.
 * @param id Segment index.
 * @param total Total segment count.
 * @param width Canvas columns.
 * @param height Canvas rows.
 * @return Partition boundary in canvas coordinates.
 */
double segment_y1(double id, double total, double width, double height) {
  return pov::segment_range(id, total, width, height).y1;
}
}
