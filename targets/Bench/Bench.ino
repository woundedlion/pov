/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Bench — stationary test image for the segmented rig (288×144)
 *
 * Target: the shipping 4× Teensy 4.0 Phantasm rig, flashed to every board.
 * Runs BenchPattern alone: one colour across the whole canvas, holding on red,
 * green, blue and white with slow ramps between. It reads with the sphere at
 * rest: every LED on every segment shows the same colour at the same instant.
 * A stable wrong segment ID is indistinguishable under this uniform pattern.
 */

#include "../Phantasm/phantasm_target.h"
#include "bench_pattern.h"

FLASHMEM void setup();
void loop();

namespace {
using Pattern = BenchPattern<CANVAS_W, CANVAS_H>;

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

FLASHMEM void setup() {
  POV::park_sync_out();
  boot_serial();
  log_reset_cause();
  POV::begin();
}

void loop() {
  // Never returns: runs the single-entry playlist forever.
  POV::run_show(EFFECT_FACTORIES, &BENCH_REVOLUTIONS);
}
