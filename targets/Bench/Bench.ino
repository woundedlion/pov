/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Bench — stationary test image for the segmented rig (288×144)
 *
 * Target: the shipping 4× Teensy 4.0 Phantasm rig, flashed to every board.
 * Runs BenchPattern alone: one colour across the whole canvas, holding on red,
 * green, blue and white with slow ramps between. Nothing in the image depends
 * on the rotor's angle, so it reads with the sphere at rest — every LED on
 * every segment shows the same colour at the same instant, and a board whose
 * ID straps, LED transport or sync has failed shows as a dark or off-colour
 * arm.
 *
 * ID straps, sync wire and LED transport are the shipping Phantasm ones
 * (targets/Phantasm/phantasm_target.h). Flash the `phantasm` env to return to
 * the show.
 */

#include "../Phantasm/phantasm_target.h"
#include "bench_pattern.h"

FLASHMEM void setup();
void loop();

namespace {
using Pattern = BenchPattern<CANVAS_W, CANVAS_H>;

// One display window opens per arm half-sweep (pov_sync.h), so the rotor's
// revolutions and the pattern's frames are related by this factor alone.
constexpr uint32_t WINDOWS_PER_REVOLUTION = 2;

static_assert(Pattern::FRAMES_PER_SECOND ==
                  RPM / 60 * WINDOWS_PER_REVOLUTION,
              "BenchPattern's hold and ramp lengths are in frames; this rotor "
              "delivers a different number of them per second");

const POV::EffectFactory EFFECT_FACTORIES[] = {&construct_effect<Pattern>};

// One epoch per colour cycle: the rebuild lands on the cycle's own wrap, and
// the commit window's blackout marks it on every board at once.
constexpr uint32_t BENCH_REVOLUTIONS[] = {Pattern::CYCLE_FRAMES /
                                          WINDOWS_PER_REVOLUTION};

static_assert(std::size(BENCH_REVOLUTIONS) == std::size(EFFECT_FACTORIES));

constexpr pov::sync::Config bench_config() {
  auto cfg = pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W, 1);
  cfg.effect_revolutions = BENCH_REVOLUTIONS;
  cfg.effect_revolutions_count = std::size(BENCH_REVOLUTIONS);
  return cfg;
}

static_assert(bench_config().valid() == nullptr,
              "Bench pov::sync::Config invariants violated");
} // namespace

FLASHMEM void setup() {
  POV::park_sync_out();
  boot_serial();
  log_reset_cause();
  create_pov();
}

void loop() {
  // Never returns: runs the single-entry playlist forever.
  g_pov->run_show(EFFECT_FACTORIES, &BENCH_REVOLUTIONS);
}
