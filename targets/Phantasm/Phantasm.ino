/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Phantasm — Multi-Teensy segmented POV display (288×144)
 *
 * Target: 4× Teensy 4.0
 * Total physical LEDs: 288 (72 per segment, 144 per arm)
 * Virtual canvas: 288×144
 *
 * Each Teensy reads its hardware ID at boot (pins 21, 22) to determine
 * which segment of the LED strip it owns.  Segment 0 is the sync master:
 * it emits count-coded symbol bursts on the single sync wire, and every
 * board (master included) generates its own columns from a local flywheel
 * timebase disciplined by those symbols.  The playlist is epoch-counted —
 * the master broadcasts an EPOCH mark when an effect's revolutions elapse
 * and all boards switch in lockstep (docs/phantasm_frame_sync_spec.md).
 *
 * Hardware ID assignment (active-low, ground to set):
 *   ID 0: both floating  — arm A outer (sync master)
 *   ID 1: pin 21 grounded — arm A inner
 *   ID 2: pin 22 grounded — arm B outer
 *   ID 3: both grounded   — arm B inner
 */

// Select the DMA HD107S output path. Must precede any include that pulls in
// led.h (via pov_segmented.h). Guarded because platformio.ini's phantasm env
// also passes -D USE_DMA_LEDS: an Arduino-IDE/VMicro build sees only this
// #define, a PlatformIO build sees only the flag, and the #ifndef keeps the two
// from colliding into a redefinition warning (the flag's value 1 vs this empty
// define).
#ifndef USE_DMA_LEDS
#define USE_DMA_LEDS
#endif

#include <FastLED.h>
#include <SPI.h>
#include <new> // std::nothrow — fail-fast OOM check on the POV allocation below

#include "pov_segmented.h"
#include "engine/effects.h"

static constexpr int TOTAL_PIXELS = 288;
static constexpr int NUM_SEGMENTS = 4;
static constexpr unsigned int RPM = 480;

using POV = POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>;

#if defined(USE_DMA_LEDS)
// Out-of-line definition for this target's controller, emitted as the required
// DMAMEM explicit specialization (see pov_segmented.h for why a generic template
// definition would silently land in DTCM and break DMA cache coherency).
HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(TOTAL_PIXELS, NUM_SEGMENTS, RPM);
#endif

namespace {
POV *g_pov;  // g_-prefixed: a bare `pov` collides with the hardware `namespace pov`

static constexpr size_t MAX_EFFECT_HEAP_BYTES = 3584;

// Foreground effect constructor for one roster entry: LUTs, arenas, init.
// Called from the driver's show loop during the epoch construction window.
template <typename E> Effect *construct_effect() {
  static_assert(sizeof(E) <= MAX_EFFECT_HEAP_BYTES,
                "Phantasm effect exceeds the heap-object budget");
  // Eager-fill the scanline LUTs before the first frame so the flywheel ISR
  // never observes a half-filled table.
  GeometryResolution<E>::init();
  // Effect constructors must not allocate from the engine arenas: here they run
  // before configure_arenas_default() (WASM's setEffect configures first), so an
  // arena allocation in a ctor would corrupt on Teensy yet pass on WASM.
  E *e = new (std::nothrow) E();
  HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
  configure_arenas_default(); // Reset before init so effects can override
  e->init();
  return e;
}

// Generated from the single-source effect roster (HS_EFFECT_LIST); the table
// order IS the playlist order, identical on all four boards (spec §6.1).
#define HS_FACTORY_ONE(name) &construct_effect<name<288, 144>>,
const POV::EffectFactory EFFECT_FACTORIES[] = {
  HS_EFFECT_LIST(HS_FACTORY_ONE)
};
#undef HS_FACTORY_ONE

// Every input to the sync protocol config is a compile-time constant on this
// board, so reject an inconsistent protocol at the build instead of the boot.
// run_show()'s runtime HS_CHECK still guards any non-constexpr instantiation.
static_assert(
    pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W, HS_EFFECT_COUNT).valid(),
    "Phantasm pov::sync::Config invariants violated");
} // namespace

void setup() {
  Serial.begin(9600); // baud inert on Teensy USB-CDC; initializes Serial only
  delay(1000);        // USB-CDC enumeration settle so early output isn't lost
  // nothrow new + HS_CHECK: a thrown bad_alloc has no handler on Teensy, so
  // fail-fast at the allocation site rather than null-deref in run_show().
  g_pov = new (std::nothrow) POV();
  HS_CHECK(g_pov != nullptr, "POV allocation failed (OOM)");
}

void loop() {
  // Never returns: the driver runs the epoch-synchronized show forever
  // (every effect plays for the same 960 revolutions = 120 s).
  g_pov->run_show(EFFECT_FACTORIES, HS_EFFECT_COUNT);
}
