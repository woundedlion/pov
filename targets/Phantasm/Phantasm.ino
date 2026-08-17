/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Phantasm — Multi-Teensy segmented POV display (288×144)
 *
 * Target: 4× Teensy 4.0 by default; 8× with PHANTASM_NUM_SEGMENTS=8
 * Total physical LEDs: 288 (72 per default segment, 144 per arm)
 * Virtual canvas: 288×144
 *
 * Each Teensy reads its hardware ID at boot (pins 21–23 as required) to determine
 * which segment of the LED strip it owns.  Segment 0 is the sync master:
 * it emits count-coded symbol bursts on the single sync wire, and every
 * board (master included) generates its own columns from a local flywheel
 * timebase disciplined by those symbols.  The playlist is epoch-counted —
 * the master broadcasts an EPOCH mark when an effect's revolutions elapse
 * and all boards switch in lockstep (docs/specs/phantasm_frame_sync_spec.md).
 *
 * Hardware ID assignment is active-low (ground to set). IDs [0, N/2) map
 * arm A; IDs [N/2, N) map arm B. ID 0 has all straps open and is the master.
 */

#include "../common/phantasm_target.h"

namespace {
// Generated from the Phantasm playlist roster; table order is identical on
// every board (spec §6.1).
#define HS_FACTORY_ONE(name, duration_seconds)                                 \
  &construct_effect<name<CANVAS_W, CANVAS_H>>,
const POV::EffectFactory EFFECT_FACTORIES[] = {
    HS_PHANTASM_EFFECT_LIST(HS_FACTORY_ONE)};
#undef HS_FACTORY_ONE

#define HS_DURATION_ONE(name, duration_seconds)                                \
  static_cast<uint32_t>(duration_seconds) * RPM / 60,
constexpr uint32_t EFFECT_REVOLUTIONS[] = {
    HS_PHANTASM_EFFECT_LIST(HS_DURATION_ONE)};
#undef HS_DURATION_ONE

#define HS_SEED_ONE(name, duration_seconds) hs::stable_effect_seed(#name),
constexpr uint64_t EFFECT_SEEDS[] = {HS_PHANTASM_EFFECT_LIST(HS_SEED_ONE)};
#undef HS_SEED_ONE

static_assert(RPM % 60 == 0,
              "Phantasm show durations require whole revolutions per second");
static_assert(std::size(EFFECT_FACTORIES) == std::size(EFFECT_REVOLUTIONS));
static_assert(std::size(EFFECT_FACTORIES) == std::size(EFFECT_SEEDS));

constexpr pov::sync::Config show_config() {
  auto cfg = pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W,
                                        HS_PHANTASM_EFFECT_COUNT);
  cfg.effect_revolutions = EFFECT_REVOLUTIONS;
  return cfg;
}

// Every input to the sync protocol config is a compile-time constant on this
// board, so reject an inconsistent protocol at the build instead of the boot.
// run_show()'s runtime HS_CHECK still guards any non-constexpr instantiation.
static_assert(show_config().valid() == nullptr,
              "Phantasm pov::sync::Config invariants violated");
} // namespace

FLASHMEM void setup() {
  POV::park_sync_out();
  boot_serial();
  log_reset_cause();
  create_pov();
}

void loop() {
  // Never returns: the driver runs the epoch-synchronized show forever
  // using the per-entry durations from HS_PHANTASM_EFFECT_LIST.
  g_pov->run_show(EFFECT_FACTORIES, &EFFECT_REVOLUTIONS, &EFFECT_SEEDS);
}
