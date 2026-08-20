/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file platform.h
 * @brief Platform layer: canvas dimensions, invariant checks, timing, clamp and
 *        the device/host split. Pulls in the attribute, diagnostics and PRNG
 *        headers so a single include still covers the whole platform surface.
 */

#include "platform/build_features.h"

/** @brief Canvas width in pixels (override via -DCANVAS_W). */
#ifndef CANVAS_W
#define CANVAS_W 288
#endif
/** @brief Canvas height in pixels (override via -DCANVAS_H). */
#ifndef CANVAS_H
#define CANVAS_H 144
#endif

// ---------------------------------------------------------------------------
/**
 * @brief Always-on invariant trap that survives NDEBUG and pulls in no stdio.
 * @param cond Condition that must hold; the macro traps when it is false.
 * @param ... Optional printf-style format string and arguments for the message.
 * @details Not stripped by NDEBUG, so it still fires in the optimized device
 *          build. Cold paths only (container growth, arena OOM, capacity guards):
 *          compiles to a single predicted-not-taken branch, never the per-pixel
 *          hot loop. On failure it logs a located breadcrumb and flushes before
 *          trapping.
 */
#define HS_CHECK(cond, ...)                                                    \
  do {                                                                         \
    if (!(cond))                                                               \
      ::hs::check_fail(__FILE__, __LINE__, #cond __VA_OPT__(, ) __VA_ARGS__);  \
  } while (0)

/**
 * @brief Optional structural audit that traps when enabled.
 * @param cond Condition that must hold; the macro traps when it is false.
 * @param ... Optional printf-style format string and arguments for the message.
 * @details Use for host-reachable O(n) seam checks. Device-reachable capacity
 *          and bounds guards remain HS_CHECK.
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
 * @details The nestable form of disable_interrupts(): a bracket closed with
 *          enable_interrupts() unmasks unconditionally, so one entered from an
 *          ISR or from inside another IRQ-off region would open it early.
 */
inline uint32_t save_disable_interrupts() {
  // The Teensy 4 core exposes __disable_irq()/__enable_irq() but not CMSIS's
  // PRIMASK accessors, so read and write the register directly.
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

/**
 * @brief Virtual rows appended below the physical LED ring (device value 3).
 * @details The latitude mapping phi = y * PI / (H + H_OFFSET - 1) lands the
 *          bottom physical row short of PI, clipping (not stretching) the image
 *          where the LEDs stop short of the south pole. The host/sim build sets
 *          H_OFFSET = 0, an intentional device/host divergence; regression tests
 *          inject the hardware value explicitly (see tests/test_geometry.h).
 *          Callers pass H (not H + H_OFFSET) to y_to_phi<H>(), which adds the
 *          offset internally.
 */
inline constexpr int H_OFFSET = 3;
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
// Test-only injectable clock (host builds only). Host timing reads the wall
// clock, so tests pin time to a fixed per-frame schedule to make the cross-run
// determinism check in tests/test_effects.h possible. OFF by default (real wall
// clock, so sim/WASM is bit-for-bit unchanged); absent from the device (ARDUINO)
// branch. The predicted-not-taken branch lives only in millis/micros (per-frame),
// never the per-pixel hot loop.
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
 * @brief Returns milliseconds since an arbitrary epoch (defined below).
 * @return Monotonic millisecond count; the beat/beatsin helpers route through it.
 */
inline unsigned long millis(); // defined below
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

/**
 * @brief Executes the guarded block at most once every N milliseconds.
 * @param N Interval in milliseconds.
 * @details Expands to a static throttle object plus one `if`, the same two-token
 * shape as the device's class-based FastLED macro (`static CEveryNMillis o(N);
 * if (o)`), so the body lines up statement-for-statement with the device's. Like
 * FastLED's, it cannot serve as the *unbraced* body of an outer control statement
 * (a leading `static` decl is not a valid lone substatement). The throttle object
 * is named from `__COUNTER__` so two uses on one source line do not collide. See
 * hs::EveryNMillis for the timing semantics.
 */
#define EVERY_N_MILLIS_I(NAME, N)                                              \
  static hs::EveryNMillis NAME((N));                                           \
  if (NAME)
#define EVERY_N_MILLIS(N) EVERY_N_MILLIS_I(HS_CONCAT(hs_every_, __COUNTER__), N)

/**
 * @brief Executes the guarded block at most once every N seconds.
 * @param N Interval in seconds.
 * @details Same two-token shape and naming scheme as EVERY_N_MILLIS. See
 * hs::EveryNSeconds for the timing semantics, which are whole-second quantized
 * rather than a millisecond throttle.
 */
#define EVERY_N_SECONDS_I(NAME, N)                                             \
  static hs::EveryNSeconds NAME((N));                                          \
  if (NAME)
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
 * @details Uses steady_clock (monotonic) so an NTP step or clock change cannot
 *          make millis() jump backward and wrap the unsigned `now - last` in
 *          EVERY_N_MILLIS. Narrowed through uint32_t so it wraps at 2^32 ms
 *          (~49 days), matching the device's 32-bit return on every host.
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
 * @details Class-based like the device's FastLED macro so EVERY_N_MILLIS expands
 * to a single guarded statement. `last` is seeded to `millis()` at construction
 * so the first evaluation waits a full period, matching the device; the stamp is
 * never reset across effect switches (function-local `static`).
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
 * @details Stamps and compares in whole seconds (FastLED's `seconds()`, i.e.
 * `millis() / 1000`), not milliseconds: a throttle constructed part-way through
 * a second first fires that fraction of a second short of a full period, which
 * is what the device does. `last` is never reset across effect switches
 * (function-local `static`).
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

/**
 * @brief Virtual rows appended below the physical LED ring; 0 on host/sim.
 * @details The simulator has no physical LED ring to clip against, so it maps
 *          the full sphere — an intentional divergence from the device's
 *          H_OFFSET = 3 (see the ARDUINO definition above).
 *
 *          HS_TEST_H_OFFSET override: a dedicated host test executable defines it
 *          to 3 so the whole pipeline compiles with the hardware offset and the
 *          south-pole renorm path runs against an energy-conservation oracle (see
 *          tests/test_h_offset_renorm.h). It MUST live in its own translation
 *          unit: offset-3 and offset-0 instantiations of PhiLUT<H>/TrigLUT<W,H>
 *          have different static-array sizes and would clash under ODR.
 */
#if defined(HS_TEST_H_OFFSET)
inline constexpr int H_OFFSET = HS_TEST_H_OFFSET;
#else
inline constexpr int H_OFFSET = 0;
#endif
} // namespace hs

// Global millis/micros if needed, though prefer namespaced
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
 * @param file Source file of the failed check (typically __FILE__).
 * @param line Source line of the failed check (typically __LINE__).
 * @param cond Stringified failed condition.
 * @param fmt printf-style message format; trailing args supply the values.
 * @details Flushes the log before trapping, so a release/device build records
 *          which invariant fired and where. Formats msg into a fixed stack
 *          buffer (no heap) so it is safe to call from a corrupted-arena / OOM
 *          context. Never returns.
 */
[[noreturn]] HS_FLASH_INLINE __attribute__((format(printf, 4, 5))) inline void
check_fail(const char *file, int line, const char *cond, const char *fmt, ...) {
  char msg[256];
  va_list args;
  va_start(args, fmt);
#ifdef ARDUINO
  // Integer-only formatter keeps newlib's float path out of ITCM (matching hs::log).
  vsniprintf(msg, sizeof(msg), fmt, args);
#else
  vsnprintf(msg, sizeof(msg), fmt, args);
#endif
  va_end(args);
  // Strip the directory so the basename does not crowd out the message in the bounded log buffer.
  const char *base = file;
  for (const char *p = file; *p; ++p) {
    if (*p == '/' || *p == '\\')
      base = p + 1;
  }
#ifdef __EMSCRIPTEN__
  char buf[256];
  snprintf(buf, sizeof(buf), "HS_CHECK failed: %s:%d: (%s) %s", base, line,
           cond, msg);
  EM_ASM({ console.error(UTF8ToString($0)); }, buf);
#else
  hs::log("HS_CHECK failed: %s:%d: (%s) %s", base, line, cond, msg);
#endif
  hs::flush_log();
  __builtin_trap();
}

// HS_CHECK(cond) with no message. Delegates with an empty formatted message
// ("%s", "") rather than passing a literal "" as the format, so no zero-length
// format string ever reaches the printf-format check (gcc -Wformat-zero-length).
[[noreturn]] HS_FLASH_INLINE inline void check_fail(const char *file, int line,
                                                    const char *cond) {
  check_fail(file, line, cond, "%s", "");
}

} // namespace hs

// ---------------------------------------------------------------------------
// Fn<Sig, Cap> — platform-aware callable wrapper. Both backends use heap-free
// inline storage, so a captured closure is never heap-allocated (which, stored in
// an ArenaVector that never destroys its elements, would leak under LSan).
//   Teensy:     teensy::inplace_function
//   Host/WASM:  hs::inplace_function
//
// Cap is a hard inline byte budget: a capture that overflows it is a compile
// error, not a heap allocation. A pointer capture is wider on the 64-bit host, so
// a pointer-capturing callsite picks a fixed byte Cap with headroom for the wider
// host closure (e.g. SpriteFn's 16 B holds two host pointers) rather than
// inflating every Fn here. See SpriteFn in concepts.h.
//
// Invoking an unbound Fn diverges: hs::inplace_function traps via check_fail,
// while the vendored teensy:: one returns a zero-initialized R. Row 9 of
// docs/ledgers/device_host_divergence_ledger.md. That zero return is a
// `static_cast<R>(0)` inside the empty vtable every teensy::inplace_function
// instantiates, so on device R must be constructible from 0: an Fn returning a
// type without that conversion compiles on host and fails only on device.
// ---------------------------------------------------------------------------
#ifdef ARDUINO
#include <inplace_function.h>
/**
 * @brief Platform-aware callable wrapper (Teensy: heap-free inplace_function).
 * @tparam Sig Call signature, e.g. void(int).
 * @tparam Cap Inline storage capacity in bytes for the captured state.
 */
template <typename Sig, size_t Cap = 16>
using Fn = teensy::inplace_function<Sig, Cap>;
#else
#include "engine/inplace_function.h" // hs::inplace_function (declared after check_fail)
/**
 * @brief Platform-aware callable wrapper (host/WASM: heap-free inplace_function).
 * @tparam Sig Call signature, e.g. void(int).
 * @tparam Cap Inline storage capacity in bytes for the captured state.
 */
template <typename Sig, size_t Cap = 16>
using Fn = hs::inplace_function<Sig, Cap>;
#endif

// Detect x86 / x64 architecture (Desktop/Simulator)
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h> // Required for SSE intrinsics
#define HS_ARCH_X86
#endif

namespace hs {

// The clamp NaN->hi contract below is load-bearing for engine-wide float->int
// domain safety: Spherical, vector_to_pixel, blend_alpha, Gradient::get and every
// palette lookup feed a possibly-NaN value through clamp as a saturating guard
// before a float->int cast. In Spherical and vector_to_pixel the guard covers the
// latitude channel only; fast_atan2 carries a NaN through the azimuth unclamped,
// so those two still require a finite input vector from the caller.
// -ffinite-math-only (implied by a bare -ffast-math)
// lets the compiler assume no NaN/Inf and fold the guard away, reintroducing the
// cast UB engine-wide; the WASM build keeps the contract by re-applying
// -fno-finite-math-only after -ffast-math (see CMakeLists.txt). The #error below
// traps at compile time on every target if that protection is ever lost.
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ != 0
#error                                                                         \
    "hs::clamp NaN->hi contract requires -fno-finite-math-only: a bare -ffast-math (or -ffinite-math-only) makes the compiler assume no NaN and folds the saturating clamp guard away, reintroducing float->int cast UB engine-wide."
#endif

#ifdef HS_ARCH_X86
// --- x86 / x64 EXPLICIT HARDWARE CLAMP ---
/**
 * @brief Clamps a float to [lo, hi] (x86 SSE backend).
 * @param v Value to clamp; a NaN maps to hi (load-bearing contract).
 * @param lo Lower bound.
 * @param hi Upper bound.
 * @return v clamped to [lo, hi]; hi when v is NaN.
 * @details CONTRACT (load-bearing): computes max(lo, min(v, hi)) with v as the
 *          FIRST operand to the inner min. The x86 minss instruction returns its
 *          SECOND source operand on any NaN, so min(NaN, hi) == hi. This backend
 *          is REORDER-SENSITIVE: swapping to min(hi, v) would yield NaN and break
 *          the contract. Do not reorder the min operands.
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
 * @param lo Lower bound.
 * @param hi Upper bound.
 * @return v clamped to [lo, hi]; hi when v is NaN.
 * @details On Cortex-M7 compiles directly to VMIN.F32 / VMAX.F32. IEEE
 *          __builtin_fminf/fmaxf are NaN-SUPPRESSING (return the non-NaN operand
 *          regardless of position), so min(NaN, hi) == hi; this backend is
 *          REORDER-INSENSITIVE, operand order kept identical to the x86 overload
 *          only for parity. NaN-suppression relies on -fno-finite-math-only
 *          surviving after -ffast-math (the __FINITE_MATH_ONLY__ preprocessor
 *          guard enforces it).
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
  return (v < lo) ? lo : ((v > hi) ? hi : v);
}

/**
 * @brief Branch-free scalar linear interpolation.
 * @param a Value at t == 0.
 * @param b Value at t == 1.
 * @param t Interpolation parameter (typically in [0, 1]).
 * @return a + (b - a) * t.
 * @details Qualify calls as hs::lerp instead of an unqualified `lerp`: the
 *          latter resolves to a std::lerp global leak, which is not guaranteed
 *          on every build.
 */
inline constexpr __attribute__((always_inline)) float lerp(float a, float b,
                                                           float t) {
  return a + (b - a) * t;
}

} // namespace hs

#include "engine/profiling.h"
