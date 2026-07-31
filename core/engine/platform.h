/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

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

#include <type_traits>
#include <cstdint>
#include <cstring>

namespace hs {
/**
 * @brief Small deterministic PRNG (PCG XSH-RR 64/32) — the process-wide RNG.
 * @details Models a UniformRandomBitGenerator, so hs::rand_* consume it
 *          unchanged. DETERMINISM CONTRACT: device and host both seed this
 *          identical type with 1337, so the draw stream stays bit-identical
 *          across the two builds (the sim/device parity invariant); nothing may
 *          depend on the specific values, only on reproducibility. Consume it
 *          only through hs:: helpers — a \<random\> algorithm or distribution
 *          draws an implementation-defined number of times and breaks the
 *          contract; use hs::shuffle, not std::shuffle. Reference: pcg32 by
 *          Melissa O'Neill.
 */
class Pcg32 {
public:
  using result_type = uint32_t;
  static constexpr result_type min() { return 0u; }
  static constexpr result_type max() { return 0xFFFFFFFFu; }

  explicit Pcg32(uint64_t seed = 1337u) { this->seed(seed); }

  /**
   * @brief Re-initializes the generator to the deterministic state for `s`.
   * @param s Seed value (mirrors std::mt19937::seed).
   */
  void seed(uint64_t s) {
    state = 0u;
    inc = (STREAM_SEQ << 1u) | 1u; // stream id must be odd
    (*this)();
    state += s;
    (*this)();
  }

  result_type operator()() {
    uint64_t old = state;
    state = old * 6364136223846793005ULL + inc;
    uint32_t xorshifted = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = static_cast<uint32_t>(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((0u - rot) & 31u));
  }

private:
  uint64_t state = 0u;
  uint64_t inc = 0u;
  static constexpr uint64_t STREAM_SEQ = 0x14057b7ef767814fULL;
};

/**
 * @brief Derives the shared-RNG seed for effect epoch @p epoch.
 * @param epoch Absolute effect/epoch index (beacon-synchronized on the device).
 * @return 1337 when epoch == 0 — the identity seed the determinism contract
 *         pins; otherwise a splitmix64-mixed value of (1337, epoch).
 * @details Used at effect handoff so every replica derives the same fresh
 *          per-visit draw stream locally from already-shared state (nothing is
 *          distributed). Integer-only, so device and host agree bit-for-bit.
 */
constexpr uint64_t epoch_seed(uint32_t epoch) {
  if (epoch == 0)
    return 1337u;
  uint64_t z = 1337u + uint64_t{epoch} * 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

/**
 * @brief Returns the global deterministic random number generator.
 * @return Reference to the process-wide Pcg32 seeded with 1337.
 * @details DETERMINISM CONTRACT: this `Pcg32(1337)` is the only RNG that is
 *          bit-identical device-vs-simulator; parity-sensitive effects must draw
 *          through it via `hs::random()`/`hs::rand_f`/`hs::rand_int`, not the
 *          FastLED `random8()`/`random16()`/Arduino `random()` path (that
 *          resolves to FastLED's LCG on device but this Pcg32 on the host mocks,
 *          so the two diverge; legacy effects only).
 *
 *          REENTRANCY CONTRACT: the generator is a function-local `static`, so it
 *          is main-loop-only — never call it from an ISR or any preemptive
 *          context.
 */
inline Pcg32 &random() {
  static Pcg32 gen(1337);
  return gen;
}

/** @brief Global debug-logging toggle. */
inline bool debug = false;
} // namespace hs

#ifdef ARDUINO
#ifndef NDEBUG
#define NDEBUG // Strip assert() to avoid linking newlib's __assert_func → fprintf
#endif
#include <Arduino.h>
#include <FastLED.h>
#include <cstdarg>
#include <cstdio>

namespace hs {
/**
 * @brief Logs one formatted line to Serial on the device.
 * @param msg printf-style format string; trailing args supply the values.
 * @details Formats into a fixed 256-byte stack buffer (no heap) via the
 *          integer-only vsniprintf, which keeps newlib's float formatter out of
 *          ITCM — the device never logs a float.
 */
inline void log(const char *msg, ...) __attribute__((format(printf, 1, 2)));
inline void log(const char *msg, ...) {
  va_list args;
  va_start(args, msg);
  char buf[256];
  vsniprintf(buf, sizeof(buf), msg, args);
  va_end(args);
  Serial.println(buf);
}

/** @brief Blocks until pending Serial output has drained (used before trap). */
inline void flush_log() { Serial.flush(); }

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
// rand_f/rand_int and ScanMetrics are defined once in the shared hs namespace
// after the #endif below.

#ifdef CORE_TEENSY
#define HS_OS_CYCLES() ARM_DWT_CYCCNT
#else
#define HS_OS_CYCLES() 0
#endif

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
#ifndef DMAMEM
#define DMAMEM
#endif

#ifndef PROGMEM
#define PROGMEM
#endif

#ifndef FLASHMEM
#define FLASHMEM
#endif

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
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

#include "engine/platform_arduino_mocks.h"

namespace hs {
/**
 * @brief Logs one formatted line to stdout on the host.
 * @param fmt printf-style format string; trailing args supply the values.
 */
inline void log(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
inline void log(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  printf("\n");
}

/** @brief Flushes stdout (used before trap so the breadcrumb is not lost). */
inline void flush_log() { fflush(stdout); }
} // namespace hs

// Mock EVERY_N_MILLIS using a simple static checker.
// Two-level macro so __COUNTER__ expands before pasting.
/**
 * @brief Pastes two tokens after expanding them (inner stage).
 * @param a Left token.
 * @param b Right token.
 */
#define HS_CONCAT_INNER(a, b) a##b
/**
 * @brief Pastes two tokens, expanding macros such as __LINE__ first.
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
#define EVERY_N_MILLIS(N) EVERY_N_MILLIS_I(HS_CONCAT(__every_, __COUNTER__), N)

/**
 * @brief Executes the guarded block at most once every N seconds.
 * @param N Interval in seconds.
 */
#define EVERY_N_SECONDS(N) EVERY_N_MILLIS((N) * 1000UL)
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
// rand_f/rand_int and ScanMetrics are defined once in the shared hs namespace
// after the #endif below.
#define HS_OS_CYCLES() 0

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
// HS_O3_BEGIN / HS_O3_END: compile the enclosed function definitions at -O3 on
// the -Os device image (selective optimization; docs/selective_o3_spec.md).
// Active only for device GCC building at -Os (__OPTIMIZE_SIZE__): the holosphere
// -O3 image, host clang, and WASM see no-ops, so those builds are byte-identical.
// The fast-math flags are restated because GCC 11's optimize pragma rebuilds
// optimization flags from defaults, dropping command-line -ffast-math /
// -fno-finite-math-only for the region (fixed in GCC 12; harmless to restate).
// HS_O3_FN is the shared single-function attribute and backs HS_COLD_MEMBER.
//
// no-unroll-loops is NOT the lever here; -funswitch-loops is. GCC 15 unswitches
// the region's per-pixel loops on their invariant branches (e.g. the feedback
// compositor's pair_on flag), emitting TWO copies of the body where GCC 11 kept
// one. Only one copy ever runs, so it buys no speed and costs pure size: the
// feedback composite lambda doubles 6,356 -> 11,588 B and Phantasm's ITCM goes
// 188,088 -> 200,520 B, overflowing FlexRAM. Disabling it restores GCC 11's
// shape at full -O3 across the roster; measured on the real image, not guessed.
// ---------------------------------------------------------------------------
#if defined(ARDUINO) && defined(__GNUC__) && !defined(__clang__) &&            \
    defined(__OPTIMIZE_SIZE__)
#define HS_O3_BEGIN                                                            \
  _Pragma("GCC push_options") _Pragma(                                         \
      "GCC optimize(\"O3\", \"fast-math\", \"no-finite-math-only\", \"no-unswitch-loops\")")
#define HS_O3_END _Pragma("GCC pop_options")
#define HS_O3_FN                                                               \
  __attribute__((optimize("O3", "fast-math", "no-finite-math-only",            \
                          "no-unswitch-loops")))
#else
#define HS_O3_BEGIN
#define HS_O3_END
#define HS_O3_FN
#endif

// ---------------------------------------------------------------------------
// HS_COLD: keep a setup-only function off the fast ITCM banks. FLASHMEM routes it
// to FLASH; noinline collapses per-call-site inline copies and noclone blocks the
// .constprop/.isra IPA clones (which drop the section attribute and land in ITCM
// regardless). Apply ONLY to internal-linkage (`static`) free functions on cold
// paths (mesh/solid construction): a section attribute on a COMDAT (inline/template
// member) function is a section-type conflict. Off-device it degrades to a no-op.
// HS_COLD_MEMBER is the COMDAT-safe variant for inline/template member functions:
// `cold` prefixes the per-function section (.text.unlikely.*) instead of naming a
// shared one, and tools/phantasm.ld routes .text.unlikely* to FLASH. On the -Os
// device image it also uses HS_O3_FN; elsewhere HS_O3_FN is empty.
// ---------------------------------------------------------------------------
#if defined(__GNUC__) && !defined(__clang__)
#define HS_COLD FLASHMEM __attribute__((noinline, noclone))
#define HS_COLD_MEMBER HS_O3_FN __attribute__((cold, noinline, noclone))
#define HS_NOINLINE_NOCLONE __attribute__((noinline, noclone))
#else
#define HS_COLD FLASHMEM
#define HS_COLD_MEMBER
#define HS_NOINLINE_NOCLONE __attribute__((noinline))
#endif

// ---------------------------------------------------------------------------
// Platform-agnostic hs:: helpers (defined once; hs::random() is defined above,
// ahead of both platform branches).
// ---------------------------------------------------------------------------
namespace hs {

// Defined later in this header; forward-declared so the helpers below can HS_CHECK.
[[noreturn]] inline void check_fail(const char *file, int line,
                                    const char *cond, const char *fmt, ...)
    __attribute__((format(printf, 4, 5)));
// No-message overload (HS_CHECK(cond) with no varargs); see definition below.
[[noreturn]] inline void check_fail(const char *file, int line,
                                    const char *cond);

/**
 * @brief Maps a raw RNG draw in [0, max] onto the half-open interval [0.0, 1.0).
 * @param value Raw RNG draw, in [0, max].
 * @param max Maximum possible draw value.
 * @pre max == UINT32_MAX: the top-band clamp constant is derived for that exact
 *      divisor. rand_f()'s static_assert enforces it for the global RNG.
 * @return A float in [0.0, 1.0), clamped just below 1.0f at the top band.
 * @details The naive value/max can land on exactly 1.0f for the top band of
 *          draws: both value and the divisor 2^32-1 round UP to 2^32 in float32
 *          (2^32-1 is not representable there), so (int)(u * N) would index N —
 *          one past the end. Clamp only those top draws to the float just below
 *          1.0f; every other draw is byte-for-byte unchanged. Pure so the
 *          boundary is unit-testable without driving the global RNG to its max.
 */
inline float random_to_unit(uint32_t value, uint32_t max) {
  float r = static_cast<float>(value) / static_cast<float>(max);
  constexpr float JUST_BELOW_ONE = 0x1.fffffep-1f; // nextafterf(1.0f, 0.0f)
  return r > JUST_BELOW_ONE ? JUST_BELOW_ONE : r;
}

/**
 * @brief Generates a pseudo-random floating-point number in [0.0, 1.0).
 * @return A random float in the half-open range [0.0, 1.0).
 */
inline float rand_f() {
  using Rng = std::remove_reference_t<decltype(hs::random())>;
  static_assert(
      Rng::max() == 0xFFFFFFFFu,
      "rand_f()/random_to_unit assume a 32-bit-range RNG; update them "
      "if the global generator changes.");
  return random_to_unit(static_cast<uint32_t>(hs::random()()),
                        static_cast<uint32_t>(hs::random().max()));
}

/**
 * @brief Generates a pseudo-random float in [min, max).
 * @param min Lower bound (inclusive).
 * @param max Upper bound (exclusive).
 * @return A random float in the half-open range [min, max).
 */
inline float rand_f(float min, float max) {
  return min + rand_f() * (max - min);
}

/**
 * @brief Generates a pseudo-random integer within a specified range.
 * @param min The minimum value (inclusive).
 * @param max The maximum value (exclusive).
 * @return A random integer in the range [min, max).
 * @note Uses `% (max - min)`, so the result is modulo-biased for ranges that do
 * not divide 2^32 evenly. Acceptable here: callers pass small setup-time ranges
 * and rely on `Pcg32` for determinism, not uniformity.
 */
inline int rand_int(int min, int max) {
  if (max > min) {
    return min + (hs::random()() % (max - min));
  }
  return min;
}

/**
 * @brief Shuffles [first, last) using the process-wide RNG.
 * @tparam It Random-access iterator.
 * @param first Range begin.
 * @param last Range end.
 * @note Replaces std::shuffle, whose permutation AND draw count are
 * implementation-defined: the three standard libraries this project builds
 * against (device libstdc++, WASM libc++, host) each produced a different
 * sequence and desynchronized the shared stream, breaking the determinism
 * contract above. Descending Fisher-Yates, exactly one draw per step.
 */
template <typename It> inline void shuffle(It first, It last) {
  const int n = static_cast<int>(last - first);
  for (int i = n - 1; i > 0; --i) {
    const int j = rand_int(0, i + 1);
    const auto tmp = first[i];
    first[i] = first[j];
    first[j] = tmp;
  }
}

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
[[noreturn]] inline void check_fail(const char *file, int line,
                                    const char *cond, const char *fmt, ...) {
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
  hs::flush_log();
#endif
  __builtin_trap();
}

// HS_CHECK(cond) with no message. Delegates with an empty formatted message
// ("%s", "") rather than passing a literal "" as the format, so no zero-length
// format string ever reaches the printf-format check (gcc -Wformat-zero-length).
[[noreturn]] inline void check_fail(const char *file, int line,
                                    const char *cond) {
  check_fail(file, line, cond, "%s", "");
}

/** @brief Scanline profiling counters (platform-agnostic). */
struct ScanMetrics {
  uint32_t plot = 0;          /**< Cycles spent plotting pixels. */
  uint32_t sdf_dist = 0;      /**< Cycles spent evaluating SDF distances. */
  uint32_t frag_shader = 0;   /**< Cycles spent in the fragment shader. */
  uint32_t bounds = 0;        /**< Cycles spent computing bounds. */
  uint32_t face_setup = 0;    /**< Cycles spent on per-face setup. */
  uint32_t scan_loop = 0;     /**< Cycles spent in the scanline loop. */
  uint32_t pixels_tested = 0; /**< Count of pixels tested. */
  uint32_t pixels_culled = 0; /**< Count of pixels culled before shading. */
  uint32_t exact_hits =
      0; /**< Count of full distance evaluations (umbrella: convex, sector, and exact-walk paths). */
  uint32_t convex_hits =
      0; /**< Count of evaluations on the convex half-plane path. */
  uint32_t sector_hits =
      0; /**< Count of evaluations on the concave sector walk (subset of exact_hits). */
  uint32_t lut_hits = 0; /**< Count of class-LUT bilinear serves. */
  uint32_t plot_backstop_hits =
      0; /**< Count of plot() steps_cache capacity-backstop trips. */
  uint32_t shade_candidates =
      0; /**< Count of pixels passing the scan's d < pixel_width test (shading + alpha-rejected). */
  /** @brief Zeroes every counter. */
  void reset() {
    plot = sdf_dist = frag_shader = bounds = face_setup = scan_loop =
        pixels_tested = pixels_culled = exact_hits = convex_hits = sector_hits =
            lut_hits = plot_backstop_hits = shade_candidates = 0;
  }
};
/** @brief Global scanline profiling counters. Compiled in only when
 *  HS_SCAN_METRICS is defined; otherwise HS_SCAN_METRIC(...) expands to nothing
 *  and this would be dead storage. */
#ifdef HS_SCAN_METRICS
inline ScanMetrics g_scan_metrics;
#endif

/**
 * @brief Per-probe cycle buckets splitting one Face::distance probe into its
 *        stages, plus the event counts each bucket divides by.
 * @details Accumulated from raw cycle-counter deltas rather than CycleScope
 * RAII, which at tens of thousands of probes per frame would cost more than the
 * stages being measured. One counter read per stage boundary; `tick` sums a
 * back-to-back read pair per probe so a capture measures its own read cost and
 * the buckets can be discounted by it. Report ratios from such a build, not
 * absolute times.
 */
struct ProbeBreakdown {
  uint32_t point = 0; /**< Cycles: probe entry through the back-face cull. */
  uint32_t project =
      0; /**< Cycles: gnomonic projection through the radius cull. */
  uint32_t edge_lut = 0;    /**< Cycles: class-LUT bilinear serve. */
  uint32_t edge_convex = 0; /**< Cycles: convex half-plane max. */
  uint32_t edge_sector = 0; /**< Cycles: concave sector walk. */
  uint32_t edge_exact = 0;  /**< Cycles: full per-edge walk. */
  uint32_t pack =
      0; /**< Cycles: plane->angle conversion and result packaging. */
  uint32_t alpha = 0;   /**< Cycles: scan-side AA coverage kernel. */
  uint32_t tick = 0;    /**< Cycles: summed back-to-back counter-read pairs. */
  uint32_t n_probe = 0; /**< Probes entered. */
  uint32_t n_cull_cos = 0; /**< Probes leaving at the back-face cull. */
  uint32_t n_cull_r = 0;   /**< Probes leaving at the radius cull. */
  uint32_t n_lut = 0;      /**< Probes served by the class LUT. */
  uint32_t n_convex = 0;   /**< Probes taking the convex path. */
  uint32_t n_sector = 0;   /**< Probes taking the sector walk. */
  uint32_t n_exact = 0;    /**< Probes taking the full edge walk. */
  uint32_t n_alpha = 0;    /**< Probes reaching the AA coverage kernel. */
  /** @brief Zeroes every bucket and count. */
  void reset() {
    point = project = edge_lut = edge_convex = edge_sector = edge_exact = pack =
        alpha = tick = 0;
    n_probe = n_cull_cos = n_cull_r = n_lut = n_convex = n_sector = n_exact =
        n_alpha = 0;
  }
};
/** @brief Global per-probe cycle buckets. Compiled in only when
 *  HS_PROBE_BREAKDOWN is defined; otherwise HS_PROBE_* expand to nothing and
 *  this would be dead storage. */
#ifdef HS_PROBE_BREAKDOWN
inline ProbeBreakdown g_probe_breakdown;
#endif

} // namespace hs

// Per-pixel scan instrumentation is OFF by default: a g_scan_metrics increment is
// a non-atomic global load-modify-store on a shared cache line for every pixel.
// Define HS_SCAN_METRICS to compile the counters in (the native test build does,
// to assert which Face::distance path each sample took); otherwise
// HS_SCAN_METRIC(...) expands to nothing.
#ifdef HS_SCAN_METRICS
#define HS_SCAN_METRIC(stmt)                                                   \
  do {                                                                         \
    (stmt);                                                                    \
  } while (0)
#else
#define HS_SCAN_METRIC(stmt) ((void)0)
#endif

// Per-probe stage timing, OFF by default: each boundary is a cycle-counter read
// plus a global accumulate, which at the scan's probe rate distorts the very
// stages it splits. Define HS_PROBE_BREAKDOWN to compile the buckets in and read
// RATIOS from the capture, discounted by the self-measured `tick` read cost.
// HS_PROBE_MARK opens a rolling timestamp; HS_PROBE_SPAN closes one stage and
// reopens the next off the same read, so a chain of N stages costs N reads.
#ifdef HS_PROBE_BREAKDOWN
#define HS_PROBE_MARK(var) uint32_t var = HS_OS_CYCLES()
#define HS_PROBE_SPAN(field, var)                                              \
  do {                                                                         \
    uint32_t hs_now = HS_OS_CYCLES();                                          \
    hs::g_probe_breakdown.field += hs_now - (var);                             \
    (var) = hs_now;                                                            \
  } while (0)
#define HS_PROBE_COUNT(field)                                                  \
  do {                                                                         \
    ++hs::g_probe_breakdown.field;                                             \
  } while (0)
#define HS_PROBE_TICK()                                                        \
  do {                                                                         \
    uint32_t hs_a = HS_OS_CYCLES();                                            \
    uint32_t hs_b = HS_OS_CYCLES();                                            \
    hs::g_probe_breakdown.tick += hs_b - hs_a;                                 \
  } while (0)
#else
#define HS_PROBE_MARK(var)
#define HS_PROBE_SPAN(field, var) ((void)0)
#define HS_PROBE_COUNT(field) ((void)0)
#define HS_PROBE_TICK() ((void)0)
#endif

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
// docs/device_host_divergence_ledger.md.
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
// before a float->int cast. -ffinite-math-only (implied by a bare -ffast-math)
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
 *          latter resolves to a platform-specific global overload (only the
 *          non-Arduino branch of this header defines one) or a std::lerp global
 *          leak, neither guaranteed on every build.
 */
inline constexpr __attribute__((always_inline)) float lerp(float a, float b,
                                                           float t) {
  return a + (b - a) * t;
}

} // namespace hs

// ---------------------------------------------------------------------------
// Cycle-counting instrumentation
//   CycleCounter — named cumulative accumulator (self-registers for bulk log)
//   CycleScope   — RAII guard that accumulates into a CycleCounter
//   HS_PROFILE   — one-liner convenience macro
// ---------------------------------------------------------------------------
namespace hs {

/**
 * @brief Formats v as decimal into buf.
 * @param v Value to format.
 * @param buf Buffer of at least 21 bytes (20 digits + NUL).
 * @return Pointer to the first digit inside buf.
 * @details Manual conversion because newlib-nano's integer printf (the -Os
 *          device build) has no long-long support.
 */
inline const char *u64_dec(uint64_t v, char *buf) {
  char *p = buf + 20;
  *p = '\0';
  do {
    *--p = static_cast<char>('0' + v % 10);
    v /= 10;
  } while (v);
  return p;
}

/**
 * @brief Named cumulative cycle accumulator. Each instance self-registers into
 *        a static intrusive list at construction so log_all()/reset_all() can
 *        walk every counter without a central registry. Counters nest: a
 *        CycleScope sets `parent` to whichever counter was active when this one
 *        started, giving log_all() a call tree with per-parent percentages.
 * @warning REENTRANCY: the registry head and the `active` nesting pointer are
 *        non-atomic statics (like hs::random()'s generator), so construction and
 *        CycleScope enter/exit are main-loop-only — driving a CycleScope from an
 *        ISR would race the list/active pointer and corrupt the call tree.
 */
struct CycleCounter {
  static constexpr uint32_t CYCLES_PER_US =
      600; /**< Core clock: Teensy 4 @ 600 MHz. */

  const char *name;    /**< Counter label used in log output. */
  uint64_t cycles = 0; /**< Accumulated cycle count. 64-bit because a
                                        32-bit accumulator overflows after only
                                        ~7 s of summed time at 600 MHz, which a
                                        multi-frame profiling run easily exceeds. */
  uint32_t count = 0;  /**< Number of timed invocations. */
  CycleCounter *parent = nullptr; /**< Enclosing counter for tree nesting. */
  CycleCounter *next =
      nullptr; /**< Next link in the intrusive registry list. */

  /**
   * @brief Constructs a named counter and self-registers it for bulk logging.
   * @param n Counter label (must outlive the counter; typically a literal).
   */
  explicit CycleCounter(const char *n) : name(n), next(head) { head = this; }

  /** @brief Zeroes this counter's accumulated cycles and call count. */
  void reset() {
    cycles = 0;
    count = 0;
  }

  /** @brief Logs every root counter (no parent) and its subtree as a tree. */
  static void log_all() {
    hs::log("--- Cycle Counters ---");
    for (auto *c = head; c; c = c->next)
      if (!c->parent && c->count)
        log_node(c, 0);
  }

  /** @brief Zeroes every registered counter (between profiling runs). */
  static void reset_all() {
    for (auto *c = head; c; c = c->next)
      c->reset();
  }

  /**
   * @brief Finds the first registered counter whose name ends with @p suffix.
   * @param suffix Name suffix to match (e.g. "_buffer_wait").
   * @return The counter, or nullptr if none is registered yet (counters
   *         self-register on first scope entry).
   */
  static CycleCounter *find_suffix(const char *suffix) {
    const size_t sl = strlen(suffix);
    for (auto *c = head; c; c = c->next) {
      const size_t nl = strlen(c->name);
      if (nl >= sl && memcmp(c->name + nl - sl, suffix, sl) == 0)
        return c;
    }
    return nullptr;
  }

private:
  static inline CycleCounter *head =
      nullptr; /**< Head of the intrusive registry list. */
  static inline CycleCounter *active =
      nullptr; /**< Currently active counter (for nesting). */
  friend struct CycleScope;

  /**
   * @brief Recursively logs one counter node and its children as a tree.
   * @param node Counter node to log.
   * @param depth Tree depth; drives indentation.
   * @details The reported percentage is this node's cycles over its parent's
   *          (or 100% for a root), and cycles are converted to microseconds via
   *          CYCLES_PER_US.
   */
  static void log_node(const CycleCounter *node, int depth) {
    if (!node->count)
      return;
    uint64_t ref = node->parent ? node->parent->cycles : node->cycles;
    uint32_t pct = ref ? (uint32_t)(node->cycles * 100 / ref) : 100;
    char cyc_buf[21], us_buf[21];
    const char *cyc = hs::u64_dec(node->cycles, cyc_buf);
    const char *us = hs::u64_dec(node->cycles / CYCLES_PER_US, us_buf);
    int indent = depth * 2;
    int name_w = 22 - indent;
    if (name_w < 1)
      name_w = 1;
    hs::log("%*s%-*s %s us (%lu%%)  %lu calls  %s cyc", indent, "", name_w,
            node->name, us, (unsigned long)pct, (unsigned long)node->count,
            cyc);
    for (auto *c = head; c; c = c->next)
      if (c->parent == node)
        log_node(c, depth + 1);
  }
};

/**
 * @brief RAII guard that times its enclosing scope and accumulates the elapsed
 *        cycles into a CycleCounter. On construction it makes its counter the
 *        active one (recording the previously-active counter as parent on first
 *        use) and snapshots the cycle counter; the destructor adds the delta and
 *        restores the previous active counter, rebuilding the nesting tree.
 */
struct CycleScope {
  CycleCounter &counter; /**< Counter this scope accumulates into. */
  CycleCounter
      *prev_active; /**< Counter to restore as active on destruction. */
  uint32_t start;   /**< Cycle snapshot taken at construction (32-bit,
                                    matching the hardware DWT CYCCNT register). */

  /**
   * @brief Begins timing the enclosing scope into the given counter.
   * @param c Counter that receives the elapsed cycles.
   * @details Makes c the active counter (recording the previously-active
   *          counter as its parent on first use) and snapshots the cycle
   *          counter.
   */
  explicit CycleScope(CycleCounter &c) : counter(c), start(HS_OS_CYCLES()) {
    prev_active = CycleCounter::active;
    if (!counter.parent && prev_active)
      counter.parent = prev_active;
    CycleCounter::active = &counter;
  }
  /**
   * @brief Adds the elapsed cycles to the counter and restores the previous one.
   * @pre The scope must not span a full CYCCNT wrap (~7 s at 600 MHz). The delta
   *      below is a 32-bit subtraction (matching the hardware register width),
   *      correct modulo 2^32, so a single scope longer than one wrap reads short
   *      by a multiple of 2^32. Accumulation across scopes is wrap-safe: the
   *      32-bit delta widens into the 64-bit `cycles` accumulator.
   */
  ~CycleScope() {
    counter.cycles += (uint32_t)(HS_OS_CYCLES() - start);
    counter.count++;
    CycleCounter::active = prev_active;
  }

  /**
   * @brief Deleted copy constructor; a scope guard must not be copied.
   */
  CycleScope(const CycleScope &) = delete;
  /**
   * @brief Deleted copy assignment; a scope guard must not be copied.
   * @return Never returns; deleted.
   */
  CycleScope &operator=(const CycleScope &) = delete;
};

/**
 * @brief ISR-safe cycle accumulator: plain single-writer fields, no registry.
 * @details CycleCounter/CycleScope are main-loop-only (non-atomic registry and
 *          nesting pointer), so ISR paths accumulate into one of these
 *          instead. Contract: the ISR is the sole writer; a foreground reader
 *          copies and reset()s under a brief IRQ-off window.
 */
struct IsrCycleStats {
  uint64_t cycles = 0;       /**< Accumulated cycles across all scopes. */
  uint32_t count = 0;        /**< Number of timed scopes. */
  uint32_t min = UINT32_MAX; /**< Shortest single scope, in cycles. */
  uint32_t max = 0;          /**< Longest single scope, in cycles. */

  /** @brief Folds one scope's elapsed cycles into the accumulator. */
  void add(uint32_t dt) {
    cycles += dt;
    ++count;
    if (dt < min)
      min = dt;
    if (dt > max)
      max = dt;
  }
  /** @brief Zeroes the accumulator (foreground, IRQs off). */
  void reset() {
    cycles = 0;
    count = 0;
    min = UINT32_MAX;
    max = 0;
  }
};

/**
 * @brief RAII guard timing its enclosing scope into an IsrCycleStats.
 */
struct IsrCycleScope {
  IsrCycleStats &stats; /**< Accumulator receiving the elapsed cycles. */
  uint32_t start;       /**< Cycle snapshot taken at construction. */

  explicit IsrCycleScope(IsrCycleStats &s) : stats(s), start(HS_OS_CYCLES()) {}
  ~IsrCycleScope() { stats.add((uint32_t)(HS_OS_CYCLES() - start)); }
  IsrCycleScope(const IsrCycleScope &) = delete;
  IsrCycleScope &operator=(const IsrCycleScope &) = delete;
};

} // namespace hs

/**
 * @brief Times the enclosing scope into an IsrCycleStats instance.
 * @param stats An hs::IsrCycleStats lvalue expression. One use per block (the
 *        guard has a fixed name; open a nested block for a second scope).
 * @details The ISR-context sibling of HS_PROFILE; compiled in only under
 *          HS_PROFILE_ENABLE.
 */
#ifdef HS_PROFILE_ENABLE
#define HS_ISR_PROFILE(stats) hs::IsrCycleScope hs_isr_scope(stats)
#else
#define HS_ISR_PROFILE(stats) ((void)0)
#endif

/**
 * @brief Times the enclosing scope into a named cycle counter.
 * @param label Counter name (used both as the identifier suffix and log label).
 * @details Compiled in only under HS_PROFILE_ENABLE; off by default so regular
 *          builds pay nothing for the per-scope bookkeeping and CYCCNT read on
 *          every hot-path face/pixel scope.
 */
#ifdef HS_PROFILE_ENABLE
#define HS_PROFILE(label)                                                      \
  static hs::CycleCounter hs_ctr_##label(#label);                              \
  hs::CycleScope hs_scope_##label(hs_ctr_##label)
#else
#define HS_PROFILE(label) ((void)0)
#endif

#if defined(HS_PROFILE_DEEP_ENABLE) && !defined(HS_PROFILE_ENABLE)
#error "HS_PROFILE_DEEP_ENABLE needs HS_PROFILE_ENABLE (the counter registry)"
#endif

/**
 * @brief Times the enclosing scope, but only in a deep-profile build.
 * @param label Counter name (used both as the identifier suffix and log label).
 * @details The form shared render code must use for any scope entered more than
 *          once per draw call — per pixel, per cell, per face. Such a scope sits
 *          on every effect's hot path, so leaving it in the ordinary `profile`
 *          image would tax the whole roster's numbers to instrument one effect
 *          (the feedback composite's per-pixel set costs ~8 cyc per entry,
 *          ~2.5% of the flush). Off unless HS_PROFILE_DEEP_ENABLE is defined on
 *          top of HS_PROFILE_ENABLE, so a deep run is opt-in per capture
 *          (HS_PROFILE_DEEP=1 in profile_one.sh, deep=1 in `just profile`).
 *          Per-frame scopes stay on plain HS_PROFILE — they are what the
 *          standard reports are built from.
 */
#ifdef HS_PROFILE_DEEP_ENABLE
#define HS_PROFILE_DEEP(label) HS_PROFILE(label)
#else
#define HS_PROFILE_DEEP(label) ((void)0)
#endif
