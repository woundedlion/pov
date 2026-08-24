/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Target boilerplate shared by the Phantasm-class sketches (targets/Phantasm,
 * targets/Profile): output-path selection, rotor/canvas geometry, the
 * POVSegmented alias and its LED-controller definition, the boot sequence, and
 * the effect construction sequence.
 *
 * Include it FIRST from the sketch — it selects the LED transport before
 * pulling in led.h — and from exactly ONE translation unit per image, since
 * HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER emits a strong definition.
 */
#pragma once

// Select the DMA HD107S output path. Guarded because the PlatformIO envs also
// pass -D USE_DMA_LEDS: an Arduino-IDE/VMicro build sees only this #define, a
// PlatformIO build sees only the flag, and the #ifndef keeps the two from
// colliding into a redefinition warning (the flag's value 1 vs this empty
// define).
#ifndef USE_DMA_LEDS
#define USE_DMA_LEDS
#endif

#ifndef PHANTASM_NUM_SEGMENTS
#define PHANTASM_NUM_SEGMENTS 4
#endif

#include <FastLED.h>
#include <SPI.h>
#include <new> // std::nothrow — fail-fast OOM check at the allocation sites

#include "pov_segmented.h"
#include "engine/effects.h"

inline constexpr int TOTAL_PIXELS = 288;
inline constexpr int NUM_SEGMENTS = PHANTASM_NUM_SEGMENTS;
inline constexpr unsigned int RPM = 480;

/** Per-effect heap-object budget for the Phantasm playlist, in bytes. */
inline constexpr size_t HS_PHANTASM_EFFECT_HEAP_BYTES = 3584;

using POV = POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>;

// Out-of-line definition for this target's controller, emitted as the required
// DMAMEM explicit specialization (see pov_segmented.h for why a generic template
// definition would silently land in DTCM and break DMA cache coherency).
HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(TOTAL_PIXELS, NUM_SEGMENTS, RPM);

namespace {
POV *g_pov; // g_ prefix: a bare `pov` collides with hardware `namespace pov`

/**
 * @brief Brings up USB serial, first step of setup().
 * @details The baud rate is inert on Teensy USB-CDC and only initializes
 * Serial; the delay lets enumeration settle so early output isn't lost.
 */
FLASHMEM void boot_serial() {
  Serial.begin(9600);
  delay(1000);
  hs::configure_debug_telemetry();
}

/**
 * @brief Logs the SoC reset cause latched since the last boot, then clears it.
 * @details Answers whether a board that stopped streaming reset itself or was
 * power-cycled by hand. A normal upload reboot reads back `por`, so that is the
 * uninformative baseline; `wdog` or `lockup-or-swreset` is the signal. SRC_SRSR
 * is write-1-to-clear and accumulates across resets, so leaving it set would
 * report every earlier boot's cause alongside this one. Bit 1 does not separate
 * a CPU lockup from a software SYSRESETREQ.
 */
FLASHMEM void log_reset_cause() {
  const uint32_t srsr = SRC_SRSR;
  SRC_SRSR = srsr;
  hs::log("reset cause: 0x%03x%s%s%s%s%s%s%s%s%s", (unsigned)srsr,
          (srsr & SRC_SRSR_IPP_RESET_B) ? " por" : "",
          (srsr & SRC_SRSR_LOCKUP_SYSRESETREQ) ? " lockup-or-swreset" : "",
          (srsr & SRC_SRSR_CSU_RESET_B) ? " csu" : "",
          (srsr & SRC_SRSR_IPP_USER_RESET_B) ? " user-reset" : "",
          (srsr & SRC_SRSR_WDOG_RST_B) ? " wdog" : "",
          (srsr & SRC_SRSR_JTAG_RST_B) ? " jtag" : "",
          (srsr & SRC_SRSR_JTAG_SW_RST) ? " jtag-sw" : "",
          (srsr & SRC_SRSR_WDOG3_RST_B) ? " wdog3" : "",
          (srsr & SRC_SRSR_TEMPSENSE_RST_B) ? " tempsense" : "");
}

/**
 * @brief Allocates the segmented POV driver into g_pov.
 * @details nothrow new + HS_CHECK: a thrown bad_alloc has no handler on Teensy,
 * so fail-fast at the allocation site rather than null-deref in run_show().
 */
FLASHMEM void create_pov() {
  g_pov = new (std::nothrow) POV();
  HS_CHECK(g_pov != nullptr, "POV allocation failed (OOM)");
}

/**
 * @brief Builds one playlist entry's effect: LUTs, arenas, construction, init.
 * @tparam E Effect type to instantiate.
 * @tparam MAX_BYTES Heap-object budget the instance must fit within.
 * @return The constructed effect, owned by the caller.
 * @details Called from the driver's show loop during the epoch construction
 * window.
 */
template <typename E, size_t MAX_BYTES = HS_PHANTASM_EFFECT_HEAP_BYTES>
Effect *construct_effect() {
  static_assert(sizeof(E) <= MAX_BYTES,
                "effect exceeds the heap-object budget");
  // Eager-fill the scanline LUTs before the first frame so the flywheel ISR
  // never observes a half-filled table.
  GeometryResolution<E>::init();
  configure_arenas_default(); // Reset before init so effects can override
  E *e = new (std::nothrow) E();
  HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
  e->init();
  return e;
}
} // namespace
