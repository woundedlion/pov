/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file arena_metrics.h
 * @brief Arena metrics reporting shared by the WASM binding headers.
 *
 * The engine arenas are read both per-frame by HolosphereEngine's memory HUD
 * and on demand by the MeshOps tooling HUD, which appends the tooling arenas
 * to the same report. Included only by targets/wasm/wasm.cpp.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/engine/memory.h"

using namespace emscripten;

/**
 * @brief Adds one arena's {usage, high_water_mark, capacity} entry to a report.
 * @param metrics Report object to extend.
 * @param name Key the entry is stored under.
 * @param arena Arena to measure.
 */
static void add_arena_metrics(val &metrics, const char *name, Arena &arena) {
  val m = val::object();
  m.set("usage", arena.get_offset());
  m.set("high_water_mark", arena.get_high_water_mark());
  m.set("capacity", arena.get_capacity());
  metrics.set(name, m);
}

/**
 * @brief Builds a {usage, high_water_mark, capacity} report for the three engine
 *        arenas.
 * @return JS object mapping each engine arena name to its metrics, in bytes.
 * @details The per-frame HUD path. Every entry costs an embind round-trip, so
 *          this covers only the arenas an engine instance can move; the tooling
 *          arenas are reported by collect_arena_metrics().
 */
static val collect_engine_arena_metrics() {
  val metrics = val::object();
  add_arena_metrics(metrics, "scratch_arena_a", scratch_arena_a);
  add_arena_metrics(metrics, "scratch_arena_b", scratch_arena_b);
  add_arena_metrics(metrics, "persistent_arena", persistent_arena);
  return metrics;
}
