/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file platform.h
 * @brief Platform layer: invariant checks, timing, clamp and the device/host
 *        split. Pulls in the attribute, diagnostics and PRNG headers so a
 *        single include still covers the whole platform surface.
 */

#include "platform/build_features.h"

#if HS_ENABLE_TEST_HOOKS
#define HS_SOURCE_FILE __FILE__
#elif defined(__FILE_NAME__)
#define HS_SOURCE_FILE __FILE_NAME__
#else
#define HS_SOURCE_FILE __FILE__
#endif

// ---------------------------------------------------------------------------
/**
 * @brief Always-on invariant trap that survives NDEBUG.
 * @param cond Condition that must hold; the macro traps when it is false.
 * @param ... Optional printf-style format string and arguments for the message.
 * @details On failure it logs a located breadcrumb and flushes before
 *          trapping.
 */
#define HS_CHECK(cond, ...)                                                    \
  do {                                                                         \
    if (!(cond))                                                               \
      ::hs::check_fail(HS_CHECK_SITE(#cond) __VA_OPT__(, ) __VA_ARGS__);       \
  } while (0)

#define HS_CHECK_STRINGIZE_IMPL(x) #x
#define HS_CHECK_STRINGIZE(x) HS_CHECK_STRINGIZE_IMPL(x)
/**
 * @brief Folds a check site into one "file:line: (cond)" string literal.
 * @param cond_str String literal naming the failed condition.
 */
#define HS_CHECK_SITE(cond_str)                                                \
  HS_SOURCE_FILE ":" HS_CHECK_STRINGIZE(__LINE__) ": (" cond_str ")"

/**
 * @brief Optional structural audit that traps when enabled.
 * @param cond Condition that must hold; the macro traps when it is false.
 * @param ... Optional printf-style format string and arguments for the message.
 * @details For host-reachable O(n) seam checks; device-reachable capacity and
 *          bounds guards use HS_CHECK.
 */
#if HS_ENABLE_STRUCTURAL_AUDITS
#define HS_AUDIT_CHECK(cond, ...) HS_CHECK(cond __VA_OPT__(, ) __VA_ARGS__)
#else
#define HS_AUDIT_CHECK(cond, ...) ((void)0)
#endif

#include <type_traits>
#include <cstdint>
#include <cstring>

#include "platform/attributes.h"
#include "platform/diagnostics.h"
#include "platform/rng.h"

#ifdef ARDUINO
#include <FastLED.h>
#include <cstdarg>
#include <cstdio>

namespace hs {
/**
 * @brief Wrapped millis() for namespace consistency.
 * @return Milliseconds since boot from the Arduino runtime.
 */
inline unsigned long millis() { return ::millis(); }
/**
 * @brief Wrapped micros() for namespace consistency.
 * @return Microseconds since boot from the Arduino runtime.
 */
inline unsigned long micros() { return ::micros(); }
/** @brief Disables interrupts (Arduino). */
inline void disable_interrupts() { noInterrupts(); }
/** @brief Enables interrupts (Arduino). */
inline void enable_interrupts() { interrupts(); }
/**
 * @brief Disables interrupts and returns the mask state that preceded it.
 * @return PRIMASK as it was on entry; hand it to restore_interrupts().
 * @details Nestable, unlike disable_interrupts()/enable_interrupts(): safe
 *          from an ISR or inside another IRQ-off region.
 */
inline uint32_t save_disable_interrupts() {
  uint32_t primask;
  __asm__ volatile("mrs %0, primask" : "=r"(primask)::"memory");
  noInterrupts();
  return primask;
}
/**
 * @brief Restores the mask state save_disable_interrupts() returned.
 * @param primask The value returned by the matching save_disable_interrupts().
 */
inline void restore_interrupts(uint32_t primask) {
  __asm__ volatile("msr primask, %0" ::"r"(primask) : "memory");
}

} // namespace hs

#else

// Non-Arduino / PC Simulation Platform
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#include <cstdint>
#include <cstdarg>
#include <cmath>
#include <algorithm>

#include <chrono>
#include <cstdio>

// ---------------------------------------------------------------------------
// Test-only injectable clock (host builds only); off by default.
namespace hs {
inline bool use_mock_time =
    false; /**< When true, millis/micros return mock values. */
inline unsigned long mock_millis_value =
    0; /**< Pinned millisecond time when mocking. */
inline unsigned long mock_micros_value =
    0; /**< Pinned microsecond time when mocking. */
/**
 * @brief Pins time to a fixed value for deterministic tests.
 * @param ms Millisecond value millis() should return.
 * @param us Microsecond value micros() should return.
 */
inline void set_mock_time(unsigned long ms, unsigned long us) {
  use_mock_time = true;
  mock_millis_value = ms;
  mock_micros_value = us;
}
/** @brief Restores the real wall clock after set_mock_time(). */
inline void clear_mock_time() { use_mock_time = false; }
/**
 * @brief Returns milliseconds since an arbitrary epoch.
 * @return Monotonic millisecond count.
 */
inline unsigned long millis();
} // namespace hs

#include "platform/arduino_mocks.h"

// Mock EVERY_N_MILLIS using a simple static checker.
// Two-level macro so __COUNTER__ expands before pasting.
/**
 * @brief Pastes two tokens after expanding them (inner stage).
 * @param a Left token.
 * @param b Right token.
 */
#define HS_CONCAT_INNER(a, b) a##b
/**
 * @brief Pastes two tokens, expanding macros such as __COUNTER__ first.
 * @param a Left token.
 * @param b Right token.
 */
#define HS_CONCAT(a, b) HS_CONCAT_INNER(a, b)

#define EVERY_N_MILLIS_I(NAME, N)                                              \
  static hs::EveryNMillis NAME((N));                                           \
  if (NAME)
/**
 * @brief Executes the guarded block at most once every N milliseconds.
 * @param N Interval in milliseconds.
 * @details Expands to a static throttle object plus one `if`, like FastLED's
 * macro, so it cannot serve as the unbraced body of an outer control statement.
 * See hs::EveryNMillis for the timing semantics.
 */
#define EVERY_N_MILLIS(N) EVERY_N_MILLIS_I(HS_CONCAT(hs_every_, __COUNTER__), N)

#define EVERY_N_SECONDS_I(NAME, N)                                             \
  static hs::EveryNSeconds NAME((N));                                          \
  if (NAME)
/**
 * @brief Executes the guarded block at most once every N seconds.
 * @param N Interval in seconds.
 * @details Same shape as EVERY_N_MILLIS; whole-second quantized (see
 * hs::EveryNSeconds).
 */
#define EVERY_N_SECONDS(N)                                                     \
  EVERY_N_SECONDS_I(HS_CONCAT(hs_every_, __COUNTER__), N)
/**
 * @brief Executes the guarded block at most once every N milliseconds (alias).
 * @param N Interval in milliseconds.
 */
#define EVERY_N_MILLISECONDS(N) EVERY_N_MILLIS(N)

namespace hs {
/**
 * @brief Returns milliseconds since an arbitrary epoch (host millis()).
 * @return Monotonic millisecond count, or the injected mock time when enabled.
 * @details steady_clock, narrowed through uint32_t so it wraps at 2^32 ms
 *          (~49 days) like the device.
 */
inline unsigned long millis() {
  if (use_mock_time)
    return mock_millis_value;
  using namespace std::chrono;
  return static_cast<uint32_t>(
      duration_cast<milliseconds>(steady_clock::now().time_since_epoch())
          .count());
}

/**
 * @brief Host throttle backing EVERY_N_MILLIS, mirroring FastLED's CEveryNMillis.
 * @details The first evaluation waits a full period; the stamp is not reset
 * across effect switches.
 */
class EveryNMillis {
public:
  explicit EveryNMillis(unsigned long interval_ms)
      : last(millis()), period(static_cast<uint32_t>(interval_ms)) {}

  /** @brief True at most once per `period` ms; stamps the trigger when it fires. */
  bool ready() {
    unsigned long now = millis();
    // 32-bit modular elapsed: matches the device's uint32 millis() wrap on LP64 hosts.
    if (static_cast<uint32_t>(now - last) >= period) {
      last = now;
      return true;
    }
    return false;
  }

  /** @brief Contextual-bool form so `if (obj)` reads as the throttle gate. */
  explicit operator bool() { return ready(); }

private:
  unsigned long last;
  uint32_t
      period; // 32-bit: matches the device wrap; caps the interval at ~49.7 days.
};

/**
 * @brief Host throttle backing EVERY_N_SECONDS, mirroring FastLED's
 *        CEveryNSeconds.
 * @details Compares whole seconds (FastLED's `seconds16()`), so the first fire
 * can come up to a second early, as on the device. `last` is not reset across
 * effect switches. The host keeps 32-bit stamps; the device truncates to 16.
 */
class EveryNSeconds {
public:
  explicit EveryNSeconds(unsigned long interval_s)
      : last(now_seconds()), period(static_cast<uint32_t>(interval_s)) {}

  /** @brief True at most once per `period` s; stamps the trigger when it fires. */
  bool ready() {
    const uint32_t now = now_seconds();
    if (now - last >= period) {
      last = now;
      return true;
    }
    return false;
  }

  /** @brief Contextual-bool form so `if (obj)` reads as the throttle gate. */
  explicit operator bool() { return ready(); }

private:
  static uint32_t now_seconds() {
    return static_cast<uint32_t>(millis()) / 1000U;
  }

  uint32_t last;
  uint32_t period;
};
/**
 * @brief Returns microseconds since an arbitrary epoch (host micros()).
 * @return Monotonic microsecond count, or the injected mock time when enabled.
 * @details Narrowed through uint32_t so it wraps at 2^32 us (~71 min), matching
 *          the device's 32-bit return on every host.
 */
inline unsigned long micros() {
  if (use_mock_time)
    return mock_micros_value;
  using namespace std::chrono;
  return static_cast<uint32_t>(
      duration_cast<microseconds>(steady_clock::now().time_since_epoch())
          .count());
}
/** @brief Disables interrupts (no-op on host). */
inline void disable_interrupts() {}
/** @brief Enables interrupts (no-op on host). */
inline void enable_interrupts() {}
/**
 * @brief Nestable interrupt disable (no-op on host).
 * @return Zero; the host has no mask to save.
 */
inline uint32_t save_disable_interrupts() { return 0; }
/** @brief Restores a saved interrupt mask (no-op on host). */
inline void restore_interrupts(uint32_t) {}

/** @brief Test-only pixel-mapping height offset (HS_TEST_H_OFFSET). */
#if defined(HS_TEST_H_OFFSET)
inline constexpr int H_OFFSET = HS_TEST_H_OFFSET;
#else
inline constexpr int H_OFFSET = 0;
#endif
} // namespace hs

/**
 * @brief Global millis() alias forwarding to hs::millis().
 * @return Monotonic millisecond count.
 */
inline unsigned long millis() { return hs::millis(); }
/**
 * @brief Global micros() alias forwarding to hs::micros().
 * @return Monotonic microsecond count.
 */
inline unsigned long micros() { return hs::micros(); }

#endif

// ---------------------------------------------------------------------------
// HS_CHECK's trap routine, defined once for both platform branches.
// ---------------------------------------------------------------------------
namespace hs {

/**
 * @brief Backing routine for HS_CHECK: logs a located breadcrumb then traps.
 * @param site Failed site as "file:line: (cond)", built by HS_CHECK_SITE.
 * @param fmt printf-style message format; trailing args supply the values.
 * @details Formats into a fixed stack buffer (no heap), so it is safe in an
 *          OOM context; flushes the log before trapping.
 */
[[noreturn]] HS_FLASH_INLINE __attribute__((format(printf, 2, 3))) inline void
check_fail(const char *site, const char *fmt, ...) {
  char msg[256];
  va_list args;
  va_start(args, fmt);
#ifdef ARDUINO
  // Integer-only formatter, as in hs::log.
  vsniprintf(msg, sizeof(msg), fmt, args);
#else
  vsnprintf(msg, sizeof(msg), fmt, args);
#endif
  va_end(args);
#ifdef __EMSCRIPTEN__
  // fd 2 routes to Module.printErr; the buffer drains on the newline.
  fprintf(stderr, "HS_CHECK failed: %s %s\n", site, msg);
  fflush(stderr);
  // wasm `unreachable` does not unwind the shadow stack, so the module is dead;
  // flag it for a JS caller that catches the RuntimeError.
  EM_ASM({ Module['HS_MODULE_DEAD'] = true; });
#elif defined(ARDUINO)
  hs::log_fragment("HS_CHECK failed: ");
  Serial.print(site);
  Serial.print(" ");
  Serial.println(msg);
#elif HS_ENABLE_TEST_HOOKS
  fprintf(stderr, "HS_CHECK failed: %s %s\n", site, msg);
  fflush(stderr);
#else
  hs::log("HS_CHECK failed: %s %s", site, msg);
#endif
  hs::flush_log();
  __builtin_trap();
}

// HS_CHECK(cond) with no message; "%s" avoids -Wformat-zero-length.
[[noreturn]] HS_FLASH_INLINE inline void check_fail(const char *site) {
  check_fail(site, "%s", "");
}

} // namespace hs

#include "platform/inplace_function.h"
/**
 * @brief Heap-free callable wrapper; invoking an empty callable traps.
 * @tparam Sig Call signature, e.g. void(int).
 * @tparam Cap Inline storage capacity in bytes for the captured state.
 */
template <typename Sig, size_t Cap = 16>
using Fn = hs::inplace_function<Sig, Cap>;

// Detect x86 / x64 architecture (Desktop/Simulator)
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h> // Required for SSE intrinsics
#define HS_ARCH_X86
#endif

namespace hs {

// hs::clamp maps NaN to hi and serves as the saturating guard before float->int
// casts; -ffinite-math-only would fold that guard away.
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ != 0
#error                                                                         \
    "hs::clamp NaN->hi contract requires -fno-finite-math-only: a bare -ffast-math (or -ffinite-math-only) makes the compiler assume no NaN and folds the saturating clamp guard away, reintroducing float->int cast UB engine-wide."
#endif

#ifdef HS_ARCH_X86
// --- x86 / x64 EXPLICIT HARDWARE CLAMP ---
/**
 * @brief Clamps a float to [lo, hi] (x86 SSE backend).
 * @param v Value to clamp; a NaN maps to hi (load-bearing contract).
 * @param lo Lower bound; must not be NaN.
 * @param hi Upper bound; must not be NaN.
 * @return v clamped to [lo, hi]; hi when v is NaN.
 * @details minss returns its second operand on NaN, so v must stay the first
 *          operand of min(v, hi).
 */
inline constexpr __attribute__((always_inline)) float clamp(float v, float lo,
                                                            float hi) {
  // The SSE intrinsics are not constant-evaluable; the NaN-suppressing builtins
  // yield the same NaN -> hi result.
  if (__builtin_is_constant_evaluated())
    return __builtin_fmaxf(lo, __builtin_fminf(v, hi));
  __m128 mv = _mm_set_ss(v);
  __m128 mlo = _mm_set_ss(lo);
  __m128 mhi = _mm_set_ss(hi);
  __m128 res = _mm_max_ss(mlo, _mm_min_ss(mv, mhi));
  return _mm_cvtss_f32(res);
}

#else
// --- NON-X86 CLAMP (Teensy Cortex-M7, WASM) ---
/**
 * @brief Clamps a float to [lo, hi] (Cortex-M7 / WASM backend).
 * @param v Value to clamp; a NaN maps to hi (same contract as the x86 backend).
 * @param lo Lower bound; must not be NaN.
 * @param hi Upper bound; must not be NaN.
 * @return v clamped to [lo, hi]; hi when v is NaN.
 * @details __builtin_fminf/fmaxf return the non-NaN operand, given
 *          -fno-finite-math-only.
 */
inline constexpr __attribute__((always_inline)) float clamp(float v, float lo,
                                                            float hi) {
  return __builtin_fmaxf(lo, __builtin_fminf(v, hi));
}
#endif

/**
 * @brief Clamps an integer to [lo, hi].
 * @param v Value to clamp.
 * @param lo Lower bound.
 * @param hi Upper bound.
 * @return v clamped to [lo, hi].
 */
inline constexpr __attribute__((always_inline)) int clamp(int v, int lo,
                                                          int hi) {
  const int upper = v > hi ? hi : v;
  return upper < lo ? lo : upper;
}

/**
 * @brief Branch-free scalar linear interpolation.
 * @param a Value at t == 0.
 * @param b Value at t == 1.
 * @param t Interpolation parameter (typically in [0, 1]).
 * @return a + (b - a) * t.
 * @details Call as hs::lerp; an unqualified `lerp` may resolve to a leaked
 *          std::lerp.
 */
inline constexpr __attribute__((always_inline)) float lerp(float a, float b,
                                                           float t) {
  return a + (b - a) * t;
}

} // namespace hs

#include "platform/profiling.h"
