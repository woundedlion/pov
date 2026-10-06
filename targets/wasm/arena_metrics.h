/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file arena_metrics.h
 * @brief Arena metrics reporting shared by the WASM binding headers.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/memory.h"

/**
 * @brief Adds one arena's {usage, high_water_mark, lifetime_high_water_mark,
 *        capacity} entry to a report.
 * @param metrics Report object to extend.
 * @param name Key the entry is stored under.
 * @param arena Arena to measure.
 * @details high_water_mark covers the window since the last peak
 *          reset/rebind; lifetime_high_water_mark folds every window in and can
 *          legitimately exceed capacity after a re-split, so an overrun gate
 *          reads the windowed mark.
 */
static void add_arena_metrics(emscripten::val &metrics, const char *name,
                              const Arena &arena) {
  emscripten::val m = emscripten::val::object();
  m.set("usage", arena.get_offset());
  m.set("high_water_mark", arena.get_high_water_mark());
  m.set("lifetime_high_water_mark", arena.get_lifetime_high_water_mark());
  m.set("capacity", arena.get_capacity());
  metrics.set(name, m);
}

/**
 * @brief Builds a {usage, high_water_mark, lifetime_high_water_mark, capacity}
 *        report for the three engine arenas.
 * @return JS object mapping each engine arena name to its metrics, in bytes.
 * @details Engine arenas only; every entry costs an embind round-trip.
 */
static emscripten::val collect_engine_arena_metrics() {
  emscripten::val metrics = emscripten::val::object();
  add_arena_metrics(metrics, "scratch_arena_a", scratch_arena_a);
  add_arena_metrics(metrics, "scratch_arena_b", scratch_arena_b);
  add_arena_metrics(metrics, "persistent_arena", persistent_arena);
  return metrics;
}
