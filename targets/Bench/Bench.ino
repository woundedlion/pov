/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file Bench.ino
 * @brief Bench — stationary test image for the segmented rig (288×144).
 *
 * Target: the shipping 4× Teensy 4.0 Phantasm rig, flashed to every board.
 * Runs BenchPattern alone: one colour across the whole canvas, holding on red,
 * green, blue and white with slow ramps between, readable with the sphere at
 * rest. A stable wrong segment ID is invisible under this uniform pattern.
 */

#include "targets/Phantasm/phantasm_target.h"
#include "targets/Bench/bench_pattern.h"

FLASHMEM void setup();
void loop();

namespace {
using Pattern = BenchPattern<CANVAS_W, CANVAS_H>;

constexpr bool STANDALONE =
#ifdef HS_BENCH_STANDALONE
    true;
#else
    false;
#endif

constexpr int SEGMENT_ID =
#ifdef HS_BENCH_SEGMENT_ID
    HS_BENCH_SEGMENT_ID;
#else
    -1;
#endif
static_assert(SEGMENT_ID >= -1 && SEGMENT_ID < NUM_SEGMENTS,
              "HS_BENCH_SEGMENT_ID must be -1 or a valid segment ID");

// One display window opens per arm half-sweep.
constexpr uint32_t WINDOWS_PER_REVOLUTION = 2;

const POV::EffectFactory EFFECT_FACTORIES[] = {&construct_effect<Pattern>};

// One epoch per colour cycle.
static_assert(Pattern::CYCLE_FRAMES % WINDOWS_PER_REVOLUTION == 0);
constexpr uint32_t BENCH_REVOLUTIONS[] = {Pattern::CYCLE_FRAMES /
                                          WINDOWS_PER_REVOLUTION};

static_assert(std::size(BENCH_REVOLUTIONS) == std::size(EFFECT_FACTORIES));

constexpr pov::sync::Config bench_config() {
  auto cfg = pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W, 1);
  cfg.set_effect_revolutions(BENCH_REVOLUTIONS);
  return cfg;
}

static_assert(bench_config().valid() == nullptr,
              "Bench pov::sync::Config invariants violated");
} // namespace

/** @brief Parks the sync output, brings up serial, then starts the driver. */
FLASHMEM void setup() {
  POV::park_sync_out();
  boot_serial();
  log_reset_cause();
  POV::begin(SEGMENT_ID);
  Serial.printf("[Bench] %s | ID %d (%s)\n",
                STANDALONE ? "Standalone; sync transmitter disabled"
                           : "Synchronized",
                POV::segment_index(), SEGMENT_ID >= 0 ? "override" : "straps");
}

/** @brief Runs the show; never returns. */
void loop() {
  // Never returns: runs the single-entry playlist forever.
  POV::run_show<std::size(EFFECT_FACTORIES)>(
      EFFECT_FACTORIES, &BENCH_REVOLUTIONS, nullptr, STANDALONE);
}
