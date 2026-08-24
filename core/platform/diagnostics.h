/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file diagnostics.h
 * @brief The platform surface the profiling layer stands on: the hs::log /
 *        hs::flush_log sink and the HS_OS_CYCLES cycle-counter read.
 */

#include "platform/attributes.h"

#ifdef ARDUINO
#ifndef NDEBUG
#define NDEBUG // Strip assert() to avoid linking newlib's __assert_func → fprintf
#endif
#include <Arduino.h>
#include <cstdarg>
#include <cstdio>

namespace hs {
/**
 * @brief Logs one formatted line to Serial on the device.
 * @param msg printf-style format string; trailing args supply the values.
 * @details Formats into a fixed 256-byte stack buffer (no heap) via the
 *          integer-only vsniprintf, which keeps newlib's float formatter out of
 *          ITCM — the device never logs a float. A longer line is truncated to
 *          255 characters; the host sink below prints it whole.
 */
HS_FLASH_INLINE inline void log(const char *msg, ...)
    __attribute__((format(printf, 1, 2)));
HS_FLASH_INLINE inline void log(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  char buf[256];
  vsniprintf(buf, sizeof(buf), msg, args);
  va_end(args);
  Serial.println(buf);
}

/**
 * @brief Logs a formatted fragment to Serial without terminating the line.
 * @param msg printf-style format string; trailing args supply the values.
 * @details Uses the same bounded integer-only formatting as log().
 */
HS_FLASH_INLINE inline void log_fragment(const char *msg, ...)
    __attribute__((format(printf, 1, 2)));
HS_FLASH_INLINE inline void log_fragment(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  char buf[256];
  vsniprintf(buf, sizeof(buf), msg, args);
  va_end(args);
  Serial.print(buf);
}

/** @brief Blocks until pending Serial output has drained (used before trap). */
inline void flush_log() { Serial.flush(); }
} // namespace hs

#ifdef CORE_TEENSY
#define HS_OS_CYCLES() ARM_DWT_CYCCNT
#else
#define HS_OS_CYCLES() 0
#endif

#else

#include <cstdarg>
#include <cstdio>

namespace hs {
/**
 * @brief Logs one formatted line to stdout on the host.
 * @param fmt printf-style format string; trailing args supply the values.
 * @details Formats straight to stdout, so the line length is unbounded and
 *          float conversions work — both diverge from the device sink above.
 */
inline void log(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
inline void log(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  printf("\n");
}

/**
 * @brief Logs a formatted fragment to stdout without terminating the line.
 * @param fmt printf-style format string; trailing args supply the values.
 */
inline void log_fragment(const char *fmt, ...)
    __attribute__((format(printf, 1, 2)));
inline void log_fragment(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
}

/** @brief Flushes stdout (used before trap so the breadcrumb is not lost). */
inline void flush_log() { fflush(stdout); }
} // namespace hs

#define HS_OS_CYCLES() 0

#endif

namespace hs {
/** @brief Build-time default for foreground telemetry logging. */
#ifdef HS_DEBUG_TELEMETRY
inline constexpr bool DEBUG_TELEMETRY_DEFAULT = true;
#else
inline constexpr bool DEBUG_TELEMETRY_DEFAULT = false;
#endif

/** @brief Global debug-logging toggle, initialized off until target setup. */
inline bool debug = false;

/** @brief Enables or disables debug telemetry from foreground code. */
inline void configure_debug_telemetry(bool enabled = DEBUG_TELEMETRY_DEFAULT) {
  debug = enabled;
}
} // namespace hs
