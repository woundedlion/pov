/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Phantasm — Multi-Teensy segmented POV display (288×144)
 *
 * Target: 4× Teensy 4.1
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

// Select the DMA HD107S output path for this target. Must be defined before any
// include that pulls in led.h (via pov_segmented.h) so the driver compiles the
// non-blocking DMA controller (24 MHz, see pov_segmented.h) rather than the
// engine-default FastLED/WS2801 fallback — the wrong chipset, and its blocking
// show() cannot meet the column ISR budget. The led.h global toggle stays
// commented so WASM/sim and single-board builds keep the FastLED path.
#define USE_DMA_LEDS

#include <FastLED.h>
#include <SPI.h>
#include <new> // std::nothrow — fail-fast OOM check on the POV allocation below

#include "pov_segmented.h"
#include "effects.h"

static constexpr int TOTAL_PIXELS = 288;
static constexpr int NUM_SEGMENTS = 4;
static constexpr unsigned int RPM = 480;

using POV = POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>;

namespace {
POV *pov;

// Foreground effect constructor for one roster entry: LUTs, arenas, init —
// the same sequence the old per-effect show<E>() performed. Called from the
// driver's show loop during the epoch construction window.
template <typename E> Effect *construct_effect() {
  // Eager-fill the scanline LUTs for this effect's resolution before the
  // first frame, so the flywheel ISR never observes a half-filled table.
  GeometryResolution<E>::init();
  E *e = new E();
  configure_arenas_default(); // Reset before init so effects can override
  e->init();
  return e;
}

// The factory table is GENERATED from the single-source effect roster
// (HS_EFFECT_LIST in core/effects.h) rather than hand-maintained, so the
// firmware roster cannot silently drift from the shipped effect set: an
// effect added to the roster joins the show automatically, and one renamed
// or removed from it becomes a compile error here instead of a stale entry.
// The table order IS the playlist order, identical on all four boards —
// which is what lets the epoch symbol say just "advance" (spec §6.1).
#define HS_FACTORY_ONE(name) &construct_effect<name<288, 144>>,
const POV::EffectFactory kEffectFactories[] = {
  HS_EFFECT_LIST(HS_FACTORY_ONE)
};
#undef HS_FACTORY_ONE
} // namespace

void setup() {
  Serial.begin(9600);
  delay(1000);
  // Fail fast at the allocation site if the heap can't hold the driver, matching
  // the project's crash-on-violation rule (see HS_CHECK on the wasm tooling
  // malloc and the arena OOM traps). nothrow new is used so the guard fires
  // regardless of whether the build compiles exceptions: a thrown bad_alloc has
  // no handler on Teensy, and a silently-null pov would be dereferenced by the
  // first run_show() call below — a null deref away from the violation site.
  pov = new (std::nothrow) POV();
  HS_CHECK(pov != nullptr, "POV allocation failed (OOM)");
}

void loop() {
  // Never returns: the driver runs the epoch-synchronized show forever
  // (every effect plays for the same 960 revolutions = 120 s).
  pov->run_show(kEffectFactories, HS_EFFECT_COUNT);
}
