/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Target boilerplate shared by the Phantasm-class sketches.
 *
 * Include it FIRST from the sketch (it selects the LED transport) and from
 * exactly ONE translation unit per image (it emits a strong definition).
 */
#pragma once

// Select the DMA HD107S output path; PlatformIO builds pass -D USE_DMA_LEDS.
#ifndef USE_DMA_LEDS
#define USE_DMA_LEDS
#endif

#ifndef PHANTASM_NUM_SEGMENTS
#define PHANTASM_NUM_SEGMENTS 4
#endif

#include <FastLED.h>
#include <SPI.h>
#include <new>

#include "core/math/geometry.h"
#include "core/memory.h"
#include "pov_segmented.h"
#include "targets/effects.h"

inline constexpr int TOTAL_PIXELS = 288;
inline constexpr int NUM_SEGMENTS = PHANTASM_NUM_SEGMENTS;
inline constexpr unsigned int RPM = 480;
static_assert(HS_SHOW_FRAMES_PER_SECOND == RPM / 60 * 2,
              "one frame per half-revolution");

/** Per-effect heap-object budget for the Phantasm playlist, in bytes. */
inline constexpr size_t HS_PHANTASM_EFFECT_HEAP_BYTES = 3584;

using POV = POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>;

// DMAMEM keeps the TX buffers out of RAM1/DTCM. Cached OCRAM requires a cache
// flush before each DMA transfer.
HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(TOTAL_PIXELS, NUM_SEGMENTS, RPM);

namespace {

/**
 * @brief Brings up USB serial after the sync output is parked.
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
 * @details A normal upload reboot reads back `por`. SRC_SRSR is
 * write-1-to-clear and accumulates across resets. Bit 1 does not separate a
 * CPU lockup from a software SYSRESETREQ.
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
 * @brief Builds one playlist entry's effect: LUTs, arenas, construction, init.
 * @tparam E Effect type to instantiate.
 * @tparam MAX_BYTES Heap-object budget the instance must fit within.
 * @return The constructed effect, owned by the caller.
 */
template <typename E, size_t MAX_BYTES = HS_PHANTASM_EFFECT_HEAP_BYTES>
Effect *construct_effect() {
  static_assert(sizeof(E) <= MAX_BYTES,
                "effect exceeds the heap-object budget");
  // The per-pixel lazy-init guards are non-atomic and rely on this eager fill.
  math::GeometryResolution<E>::init();
  configure_arenas_default(); // Reset before init so effects can override
  E *e = new (std::nothrow) E();
  HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
  e->init();
  return e;
}
} // namespace
