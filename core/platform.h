/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Canvas resolution — override per target via -D flags if needed
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
 * @details Unlike assert(), HS_CHECK is NOT stripped by NDEBUG, so it still
 *          fires in the optimized device build — where platform.h defines
 *          NDEBUG to keep newlib's __assert_func -> fprintf out of the image.
 *          Use it on COLD paths only (container growth, arena OOM, capacity
 *          guards), where an invariant violation is a logic/sizing bug with no
 *          valid recovery: trapping at the violation site is strictly better
 *          than silently writing out of bounds, and a corrupted arena that ships
 *          garbage is the worst outcome. It compiles to a single
 *          predicted-not-taken branch — never place it in the per-pixel hot
 *          loop. Reserve bounded/soft handling for genuine transient conditions
 *          (DMA overrun, dropped frame); those are not invariant violations.
 *
 *          On failure it logs "HS_CHECK failed: <file>:<line>: (<cond>) <msg>"
 *          and flushes the log before trapping, so a release/device build leaves
 *          a breadcrumb identifying exactly which check fired. Usage:
 *            HS_CHECK(i < n);
 *            HS_CHECK(i < n, "index out of range");
 *            HS_CHECK(i < n, "i=%d n=%d", i, n);
 *          hs::check_fail is defined later in this header; the macro only
 *          expands in translation units that include the whole file, so the
 *          forward use is fine.
 */
#define HS_CHECK(cond, ...)                                                     \
  do {                                                                          \
    if (!(cond))                                                                \
      ::hs::check_fail(__FILE__, __LINE__, #cond, "" __VA_ARGS__);             \
  } while (0)

#include <random>
#include <type_traits>
#include <cstdint>

namespace hs {
/**
 * @brief Small deterministic PRNG (PCG XSH-RR 64/32) — the process-wide RNG.
 * @details Replaces std::mt19937 as the engine behind hs::random(): it cuts the
 *          global generator state from ~2,500 B to 16 B of DTCM (RAM1) and has a
 *          shorter critical path. Models a UniformRandomBitGenerator (result_type
 *          / min() / max() / operator()), so std::shuffle and hs::rand_* consume
 *          it unchanged.
 *
 *          DETERMINISM CONTRACT: device and host both instantiate this identical
 *          type seeded with 1337, so the draw stream stays bit-identical across
 *          the two builds — the sim/device parity invariant is preserved. Only
 *          the algorithm changed, so the sequence differs from the former
 *          mt19937; nothing may depend on the specific values, only on
 *          reproducibility. Reference implementation by Melissa O'Neill (pcg32).
 */
class Pcg32 {
public:
  using result_type = uint32_t;
  static constexpr result_type min() { return 0u; }
  static constexpr result_type max() { return 0xFFFFFFFFu; }

  explicit Pcg32(uint64_t seed = 1337u) { this->seed(seed); }

  /**
   * @brief Re-initializes the generator to the deterministic state for `s`.
   * @param s Seed value (mirrors std::mt19937::seed so callers pinning a known
   *          stream — the determinism tests — work unchanged).
   */
  void seed(uint64_t s) {
    state_ = 0u;
    inc_ = (kStreamSeq << 1u) | 1u; // stream id must be odd
    (*this)();
    state_ += s;
    (*this)();
  }

  result_type operator()() {
    uint64_t old = state_;
    state_ = old * 6364136223846793005ULL + inc_;
    uint32_t xorshifted = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = static_cast<uint32_t>(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((0u - rot) & 31u));
  }

private:
  uint64_t state_ = 0u;
  uint64_t inc_ = 0u;
  static constexpr uint64_t kStreamSeq = 0x14057b7ef767814fULL;
};
} // namespace hs

#ifdef ARDUINO
#ifndef NDEBUG
#define NDEBUG  // Strip assert() to avoid linking newlib's __assert_func → fprintf
#endif
#include <Arduino.h>
#include <FastLED.h>
#include <cstdarg>
#include <cstdio>

namespace hs {
/**
 * @brief Logs one formatted line to Serial on the device.
 * @param msg printf-style format string; trailing args supply the values.
 * @details Formats into a fixed 256-byte stack buffer with vsniprintf (no heap),
 *          then writes one line. Sized to hold a full check_fail() breadcrumb
 *          ("HS_CHECK failed: file:line: (cond) msg") without truncating the
 *          message tail. The integer-only vsniprintf is deliberate: it keeps
 *          newlib's float formatter (_dtoa_r + the %f/%g bignum helpers, ~5 KB)
 *          out of ITCM. The device never logs a float, so this costs nothing.
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
 * @brief Returns the global deterministic random number generator.
 * @return Reference to the process-wide Pcg32 seeded with 1337.
 * @details DETERMINISM CONTRACT: this `Pcg32(1337)` is the only RNG that is
 *          bit-identical device-vs-simulator — the host build (see the `#else`
 *          branch below) returns the same seeded `Pcg32`. Effects that must
 *          render identically on hardware and in the byte-deterministic
 *          simulator therefore draw through this generator: `hs::random()`,
 *          `hs::rand_f`, `hs::rand_int`.
 *
 *          The bare FastLED `random8()`/`random16()`/Arduino `random()` are NOT
 *          covered by this contract: on device they resolve to FastLED's LCG
 *          (seeded by `randomSeed(1337)` in pov_single.h), but the host mocks
 *          route them through this `Pcg32` — a different algorithm, so the two
 *          diverge. Only legacy effects (`effects_legacy.h`) use that path;
 *          modern effects must not, and the cheap LCG is deliberately not
 *          adopted on the hot path because it would break this contract for no
 *          measurable gain (RNG is spawn/setup-time, never per-pixel). Pcg32 is
 *          itself a fast deterministic PRNG, so a future per-pixel RNG need can
 *          draw from here rather than reaching for the platform LCG.
 *
 *          REENTRANCY CONTRACT: the generator is a function-local `static`, so
 *          advancing it mutates shared state with no lock. It is therefore
 *          main-loop-only — never call `hs::random()`/`rand_f`/`rand_int` from an
 *          ISR or any preemptive context. Interleaving a draw from an interrupt
 *          would both corrupt the generator's internal state (a torn multi-word
 *          update) and desync the deterministic stream from the simulator. Safe
 *          today only because every caller runs on the render/setup path.
 */
inline Pcg32& random() {
  static Pcg32 gen(1337);
  return gen;
}
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

// Global state
/** @brief Global debug-logging toggle. */
inline bool debug = false;
// rand_f/rand_int and ScanMetrics are platform-agnostic — defined once in the
// shared hs namespace after the #endif below.

#ifdef CORE_TEENSY
#define HS_OS_CYCLES() ARM_DWT_CYCCNT
#else
#define HS_OS_CYCLES() 0
#endif

/**
 * @brief Virtual rows appended below the physical LED ring (device value 3).
 * @details The latitude mapping is phi = y * PI / (H + H_OFFSET - 1), so the
 *          bottom physical row y = H-1 lands short of PI: the image is clipped
 *          (not stretched) where the LEDs stop short of the south pole. The
 *          host/sim build deliberately sets H_OFFSET = 0 (see the other
 *          definition in the non-Arduino branch below), so the simulator maps
 *          the full sphere and does not reproduce this bottom clipping — an
 *          intentional device/host divergence. Because the native build cannot
 *          observe a non-zero offset, the regression tests inject the hardware
 *          value explicitly (see tests/test_geometry.h). Callers pass H (not
 *          H + H_OFFSET) to y_to_phi<H>(), which adds the offset internally.
 */
static constexpr int H_OFFSET = 3;
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

#include <cstring>
#include <chrono>
#include <iostream>

// ---------------------------------------------------------------------------
// Test-only injectable clock (host builds only).
//
// Animation timing on host reads the wall clock (hs::millis / hs::micros /
// beatsin*), so the same frame rendered twice sees a different time and cannot
// be compared byte-for-byte. Tests enable this seam to pin time to a fixed
// per-frame schedule, which is what makes the cross-run determinism check in
// tests/test_effects.h possible.
//
// OFF by default: with use_mock_time == false the helpers return the real wall
// clock, so the simulator / WASM build is bit-for-bit unchanged. The seam does
// not exist in the device (ARDUINO) branch at all, so hardware is untouched.
// The single predicted-not-taken branch lives only in millis/micros, which are
// per-frame calls — never the per-pixel hot loop.
namespace hs {
inline bool use_mock_time = false;        /**< When true, millis/micros return mock values. */
inline unsigned long mock_millis_value = 0; /**< Pinned millisecond time when mocking. */
inline unsigned long mock_micros_value = 0; /**< Pinned microsecond time when mocking. */
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

/**
 * @brief HSV color structure for non-Arduino platforms.
 * @details Mirrors FastLED's CHSV so host effects compile unchanged.
 */
struct CHSV {
  /** @brief Hue, saturation and value, each in [0, 255]. */
  uint8_t h, s, v;
  /**
   * @brief Constructs a zero-initialised (black) color.
   */
  constexpr CHSV() : h(0), s(0), v(0) {}
  /**
   * @brief Constructs a color from explicit hue, saturation and value.
   * @param h Hue in [0, 255].
   * @param s Saturation in [0, 255].
   * @param v Value/brightness in [0, 255].
   */
  constexpr CHSV(uint8_t h, uint8_t s, uint8_t v) : h(h), s(s), v(v) {}
};

// --- Mock FastLED Types ---
/**
 * @brief RGB color structure mimicking FastLED's CRGB.
 * @details Reproduces FastLED's constructors, operators and helpers so host
 *          effects compile and behave identically to the device.
 */
struct CRGB {
  /** @brief Red, green and blue channels, each in [0, 255]. */
  uint8_t r, g, b;
  /**
   * @brief Constructs a zero-initialised (black) color.
   */
  constexpr CRGB() : r(0), g(0), b(0) {}
  /**
   * @brief Constructs a color from explicit channel values.
   * @param r Red channel in [0, 255].
   * @param g Green channel in [0, 255].
   * @param b Blue channel in [0, 255].
   */
  constexpr CRGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
  /**
   * @brief Constructs a color from a packed 0xRRGGBB code.
   * @param colorcode Packed color; the low 24 bits decode as red, green, blue.
   *        This is a packed color, not a grayscale fill.
   */
  constexpr CRGB(uint32_t colorcode)
      : r((colorcode >> 16) & 0xFF), g((colorcode >> 8) & 0xFF),
        b(colorcode & 0xFF) {}

  /**
   * @brief Constructs an RGB color by converting from HSV.
   * @param hsv Source color in HSV space.
   * @details Basic integer HSV-to-RGB conversion over six hue sectors.
   */
  constexpr CRGB(const CHSV &hsv) {
    unsigned char region, remainder, p, q, t;

    if (hsv.s == 0) {
      r = hsv.v;
      g = hsv.v;
      b = hsv.v;
      return;
    }

    region = hsv.h / 43;
    remainder = (hsv.h - (region * 43)) * 6;

    p = (hsv.v * (255 - hsv.s)) >> 8;
    q = (hsv.v * (255 - ((hsv.s * remainder) >> 8))) >> 8;
    t = (hsv.v * (255 - ((hsv.s * (255 - remainder)) >> 8))) >> 8;

    switch (region) {
    case 0:
      r = hsv.v;
      g = t;
      b = p;
      break;
    case 1:
      r = q;
      g = hsv.v;
      b = p;
      break;
    case 2:
      r = p;
      g = hsv.v;
      b = t;
      break;
    case 3:
      r = p;
      g = q;
      b = hsv.v;
      break;
    case 4:
      r = t;
      g = p;
      b = hsv.v;
      break;
    default:
      r = hsv.v;
      g = p;
      b = q;
      break;
    }
  }

  /**
   * @brief Tests two colors for exact channel equality.
   * @param rhs Color to compare against.
   * @return true if all three channels match.
   */
  bool operator==(const CRGB &rhs) const {
    return r == rhs.r && g == rhs.g && b == rhs.b;
  }
  /**
   * @brief Tests two colors for channel inequality.
   * @param rhs Color to compare against.
   * @return true if any channel differs.
   */
  bool operator!=(const CRGB &rhs) const { return !(*this == rhs); }

  /**
   * @brief Adds another color into this one with per-channel saturation.
   * @param rhs Color to add; each channel clamps to [0, 255].
   * @return Reference to this color after the saturated add.
   */
  CRGB &operator+=(const CRGB &rhs) {
    r = (r + rhs.r > 255) ? 255 : r + rhs.r;
    g = (g + rhs.g > 255) ? 255 : g + rhs.g;
    b = (b + rhs.b > 255) ? 255 : b + rhs.b;
    return *this;
  }

  /**
   * @brief Linearly interpolates toward another color (FastLED lerp16).
   * @param other Target color at frac == 65535.
   * @param frac Interpolation fraction in [0, 65535], where 0 yields *this.
   * @return The interpolated color.
   */
  CRGB lerp16(const CRGB &other, uint16_t frac) const {
    CRGB ret;
    // frac is 0..65535. The >> 16 truncates (no +0x8000 round bias) to match
    // FastLED's integer fixed-point lerp16by16/scale16, which also truncate — so
    // the host mock and the device stay byte-identical. Add a +0x8000 round only
    // if a future FastLED variant this mocks starts rounding.
    ret.r = static_cast<uint8_t>((static_cast<uint32_t>(r) * (65535 - frac) +
                                  static_cast<uint32_t>(other.r) * frac) >>
                                 16);
    ret.g = static_cast<uint8_t>((static_cast<uint32_t>(g) * (65535 - frac) +
                                  static_cast<uint32_t>(other.g) * frac) >>
                                 16);
    ret.b = static_cast<uint8_t>((static_cast<uint32_t>(b) * (65535 - frac) +
                                  static_cast<uint32_t>(other.b) * frac) >>
                                 16);
    return ret;
  }
};

// --- Mock FastLED Functions ---
/**
 * @brief Saturated 8-bit addition (FastLED qadd8).
 * @param i First addend in [0, 255].
 * @param j Second addend in [0, 255].
 * @return i + j clamped to a maximum of 255.
 */
inline uint8_t qadd8(uint8_t i, uint8_t j) {
  int t = i + j;
  if (t > 255)
    t = 255;
  return t;
}

/**
 * @brief Saturated 8-bit subtraction (FastLED qsub8).
 * @param i Minuend in [0, 255].
 * @param j Subtrahend in [0, 255].
 * @return i - j clamped to a minimum of 0.
 */
inline uint8_t qsub8(uint8_t i, uint8_t j) {
  int t = i - j;
  if (t < 0)
    t = 0;
  return t;
}

#ifndef PI // Arduino/FastLED define PI on-device; guard the host definition.
#define PI 3.1415926535897932384626433832795
#endif

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

/**
 * @brief Returns the global deterministic random number generator.
 * @return Reference to the process-wide Pcg32 seeded with 1337.
 * @details Mirrors the device `Pcg32(1337)` so that `hs::random()`/`hs::rand_*`
 *          reproduce hardware bit-for-bit — see the determinism contract on the
 *          ARDUINO branch's `hs::random()` for which paths are (and are not)
 *          covered.
 */
inline Pcg32& random() {
  static Pcg32 gen(1337);
  return gen;
}
} // namespace hs

// --- Mock Arduino Constants/Types ---
/** @brief Arduino-style alias for an unsigned 8-bit byte. */
using byte = uint8_t;
/** @brief Arduino-style alias for a boolean. */
using boolean = bool;

// --- Mock FastLED Constants ---
/** @brief Mock FastLED color/temperature correction selectors. */
enum FastLEDCheck {
  UncorrectedColor,
  TypicalLEDStrip,
  UncorrectedTemperature,
  Candle
};

/**
 * @brief Mock implementation of the FastLED controller for simulation.
 * @details Every method is a no-op; the simulator renders into its own
 *          framebuffer rather than driving real LEDs.
 */
struct FastLEDMock {
  /**
   * @brief Sets the color-correction profile (no-op on host).
   */
  void setCorrection(int) {}
  /**
   * @brief Sets the color-temperature profile (no-op on host).
   */
  void setTemperature(int) {}
  /**
   * @brief Registers an LED strip (no-op on host).
   * @tparam T LED chipset type.
   * @tparam P1 First chipset template parameter.
   * @tparam P2 Second chipset template parameter.
   * @tparam P3 Third chipset template parameter.
   * @tparam P4 Fourth chipset template parameter.
   */
  template <typename T, int P1, int P2, int P3, int P4>
  void addLeds(CRGB *, int) {}
  /**
   * @brief Pushes the framebuffer to the strip (no-op on host).
   */
  void show() {}
  /**
   * @brief Fills the whole strip with a single color (no-op on host).
   */
  void showColor(const CRGB &) {}
};
inline FastLEDMock FastLED;

/** @brief Mock LED chipset selector for addLeds template arguments. */
enum LEDType { WS2801 };
/** @brief Mock LED color-order selector for addLeds template arguments. */
enum ColorOrder { RGB };
#define DATA_RATE_MHZ(x) x

// --- Mock Arduino Functions ---
// These are deliberately at global scope to mirror Arduino/FastLED, which
// expose random()/map() as free globals; unqualified callers (e.g. the legacy
// effects) must resolve identically on host and device.
/**
 * @brief Returns a pseudo-random integer in [0, max) (Arduino random()).
 * @param max Exclusive upper bound.
 * @return A value in [0, max), or 0 when max <= 0.
 * @details Guards a degenerate range like the device: Arduino's random(0)
 *          returns 0, so the host avoids a modulo-by-zero SIGFPE and matches.
 */
inline int random(int max) {
  if (max <= 0) return 0;
  return hs::random()() % max;
}
/**
 * @brief Returns a pseudo-random integer in [min, max) (Arduino random()).
 * @param min Inclusive lower bound.
 * @param max Exclusive upper bound.
 * @return A value in [min, max), or min when max <= min.
 * @details Mirrors the device: random(min, max) with min >= max returns min,
 *          so the host avoids a modulo-by-zero SIGFPE.
 */
inline int random(int min, int max) {
  if (max <= min) return min;
  return min + (hs::random()() % (max - min));
}
/**
 * @brief Re-maps a value from one integer range to another (Arduino map()).
 * @param x Value to map.
 * @param in_min Lower bound of the input range.
 * @param in_max Upper bound of the input range.
 * @param out_min Lower bound of the output range.
 * @param out_max Upper bound of the output range.
 * @return x scaled from [in_min, in_max] onto [out_min, out_max]; out_min when
 *         the input range is degenerate (in_max == in_min).
 * @details Arduino's map() divides by (in_max - in_min) with no guard. A
 *          degenerate input range raises SIGFPE on the host (x86 integer
 *          div-by-zero) while the Cortex-M7 device returns 0 from the divide
 *          (SDIV-by-zero traps are off), so map() there yields out_min. Match
 *          the device rather than crashing only in the simulator.
 */
inline long map(long x, long in_min, long in_max, long out_min, long out_max) {
  if (in_max == in_min) return out_min;
  return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

// --- System Mock ---
/**
 * @brief Mock of Arduino's Serial that writes to std::cout on the host.
 */
struct SerialMock {
  /**
   * @brief Initialises the port (no-op on host).
   */
  void begin(int) {}
  /**
   * @brief Writes a C string without a trailing newline.
   * @param msg Text to write.
   */
  void print(const char *msg) { std::cout << msg; }
  /**
   * @brief Writes a signed integer without a trailing newline.
   * @param val Value to write.
   */
  void print(int val) { std::cout << val; }
  /**
   * @brief Writes an unsigned long without a trailing newline.
   * @param val Value to write.
   */
  void print(unsigned long val) { std::cout << val; }
  /**
   * @brief Writes a float without a trailing newline.
   * @param val Value to write.
   */
  void print(float val) { std::cout << val; }
  /**
   * @brief Writes a C string followed by a newline.
   * @param msg Text to write.
   */
  void println(const char *msg) { std::cout << msg << std::endl; }
  /**
   * @brief Writes a signed integer followed by a newline.
   * @param val Value to write.
   */
  void println(int val) { std::cout << val << std::endl; }
  /**
   * @brief Writes an unsigned long followed by a newline.
   * @param val Value to write.
   */
  void println(unsigned long val) { std::cout << val << std::endl; }
  /**
   * @brief Formats and writes a printf-style message (Arduino Serial.printf).
   * @param fmt printf-style format string.
   * @details Expands into a fixed 256-byte stack buffer (no heap) and emits.
   * @warning Host/device divergence (mirrors the note on hs::log and check_fail):
   *          this host mock uses full vsnprintf, so `%f`/`%g` format here. The
   *          device hs::log/check_fail use integer-only vsniprintf to keep
   *          newlib's float formatter out of ITCM, so a float conversion that
   *          works in the simulator silently drops on hardware. vsniprintf is a
   *          newlib extension and does not exist on the host, so the host cannot
   *          simply match it — avoid `%f`/`%g` in any message destined for the
   *          device path.
   */
  void printf(const char *fmt, ...) {
    char buf[256];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    std::cout << buf;
  }
};
inline SerialMock Serial;

// --- FastLED Mocks ---
// NOTE: these route through hs::random() (Pcg32) only to give the host a
// reproducible stream — they do NOT match the device, where random8/random16
// are FastLED's LCG. This path is therefore per-platform and used only by
// legacy effects; modern effects use hs::rand_* (see the determinism contract
// on the ARDUINO branch's hs::random()).
/**
 * @brief Returns a pseudo-random 8-bit value (FastLED random8).
 * @return A value in [0, 255].
 */
inline uint8_t random8() { return hs::random()() % 256; }
/**
 * @brief Returns a pseudo-random 8-bit value below a limit (FastLED random8).
 * @param top Exclusive upper bound.
 * @return A value in [0, top), or 0 when top == 0.
 * @details FastLED's random8(lim) is a scaled multiply ((r*lim)>>8), so
 *          random8(0)==0 on the device; the host modulo would SIGFPE, so guard
 *          to match rather than crash.
 */
inline uint8_t random8(uint8_t top) {
  if (top == 0) return 0;
  return hs::random()() % top;
}
/**
 * @brief Returns a pseudo-random 16-bit value (FastLED random16).
 * @return A value in [0, 65535].
 */
inline uint16_t random16() { return hs::random()() % 65536; }
/**
 * @brief Mixes additional entropy into the RNG (no-op on host).
 */
inline void random16_add_entropy(uint16_t) {}

// FastLED fixed-point sine + scaling primitives. The device build pulls these
// from <FastLED.h>; the host mock reproduces their exact integer semantics (LUT
// sine, 8.8 beat sawtooth, scale8 range fit) so the simulator predicts the
// device instead of approximating it with a smooth float sine of a different
// shape, phase convention, and parameter order.

/**
 * @brief Unsigned 8-bit fractional scale, scale8(i, sc) = i * (1 + sc) / 256.
 * @param i Value to scale, in [0, 255].
 * @param sc Scale factor, in [0, 255].
 * @return i scaled by sc/256 in the SCALE8_FIXED sense, in [0, 255].
 * @details The (1 + sc) is FastLED's SCALE8_FIXED form (default in modern
 *          releases): it makes scale8(x, 255) == x, so a full-scale fade is the
 *          identity. The device pulls the FIXED variant from <FastLED.h>;
 *          matching it here keeps the simulator bit-exact rather than 1 LSB low
 *          on every fade.
 */
inline uint8_t scale8(uint8_t i, uint8_t sc) {
  return (static_cast<uint16_t>(i) * (1 + static_cast<uint16_t>(sc))) >> 8;
}
/**
 * @brief 16-bit SCALE8_FIXED counterpart, scale16(i, sc) = i * (1 + sc) / 65536.
 * @param i Value to scale, in [0, 65535].
 * @param sc Scale factor, in [0, 65535].
 * @return i scaled by sc/65536 in the SCALE8_FIXED sense, in [0, 65535].
 */
inline uint16_t scale16(uint16_t i, uint16_t sc) {
  return (static_cast<uint32_t>(i) * (1 + static_cast<uint32_t>(sc))) >> 16;
}

/**
 * @brief FastLED's 8-bit LUT sine (sin8_C), bit-exact with the device.
 * @param theta Phase in [0, 255], one full turn.
 * @return Sine in [0, 255] centred on 128 (sin8(0)==128, sin8(64)==255).
 */
inline uint8_t sin8(uint8_t theta) {
  static const uint8_t b_m16_interleave[] = {0, 49, 49, 41, 90, 27, 117, 10};
  uint8_t offset = theta;
  if (theta & 0x40) offset = 255 - offset;
  offset &= 0x3F;
  uint8_t secoffset = offset & 0x0F;
  if (theta & 0x40) secoffset++;
  uint8_t section = offset >> 4;
  const uint8_t *p = b_m16_interleave + section * 2;
  uint8_t b = *p++;
  uint8_t m16 = *p;
  uint8_t mx = (m16 * secoffset) >> 4;
  int8_t y = mx + b;
  if (theta & 0x80) y = -y;
  return static_cast<uint8_t>(y + 128);
}

/**
 * @brief FastLED's 16-bit LUT sine (sin16_C), bit-exact with the device.
 * @param theta Phase in [0, 65535], one full turn.
 * @return Signed sine in [-32767, 32767].
 */
inline int16_t sin16(uint16_t theta) {
  static const uint16_t base[] = {0,     6393,  12539, 18204,
                                  23170, 27245, 30273, 32137};
  static const uint8_t slope[] = {49, 48, 44, 38, 31, 23, 14, 4};
  uint16_t offset = (theta & 0x3FFF) >> 3; // 0..2047
  if (theta & 0x4000) offset = 2047 - offset;
  uint8_t section = offset / 256; // 0..7
  uint16_t b = base[section];
  uint8_t m = slope[section];
  uint8_t secoffset8 = static_cast<uint8_t>(offset) / 2;
  int16_t y = static_cast<int16_t>(m * secoffset8) + b;
  if (theta & 0x8000) y = -y;
  return y;
}

/**
 * @brief Free-running 16-bit sawtooth phase at the given 8.8 fixed-point BPM.
 * @param bpm88 Tempo as an 8.8 fixed-point beats-per-minute value.
 * @param timebase Millisecond offset for the zero of time.
 * @return The current phase in [0, 65535].
 * @details Sourced from hs::millis() so the test time-injection seam keeps
 *          beats deterministic. The * 280 constant is FastLED's
 *          (≈ 65536 * 1000 / 60000) ms→phase scale.
 *
 *          Sim/device wrap: the intermediate `(millis-timebase)*bpm88*280` is
 *          `unsigned long`, which is 32-bit on the device (and Win32 host) but
 *          64-bit on a LP64 host, so the device wraps it mod 2^32 while a LP64
 *          host does not. This does NOT diverge: the function returns only the
 *          uint16_t formed by `>>16`, i.e. bits 16..31 of the product, and by
 *          the modular-multiply identity those bits are identical whether the
 *          product is taken mod 2^32 (device) or in full (LP64) — the high bits
 *          a 64-bit intermediate keeps are exactly the ones the uint16_t cast
 *          discards. Verified exhaustively over random millis/timebase/bpm88
 *          (incl. the timebase>millis underflow case). Hence no uint32_t cast is
 *          needed; the host beat/beatsin phases match the device bit-for-bit.
 */
inline uint16_t beat88(uint16_t bpm88, uint32_t timebase = 0) {
  return ((hs::millis() - timebase) * bpm88 * 280) >> 16;
}
/**
 * @brief Free-running 16-bit sawtooth phase at the given whole BPM.
 * @param bpm Tempo in beats per minute.
 * @param timebase Millisecond offset for the zero of time.
 * @return The current phase in [0, 65535].
 * @pre bpm <= 255. bpm << 8 is truncated to 16 bits, so bpm >= 256 wraps
 *      (e.g. beat16(300) == beat16(44)). Parity-faithful to FastLED's own
 *      beat16; >255 BPM is musically nonsensical, hence a precondition not a
 *      guard.
 */
inline uint16_t beat16(uint16_t bpm, uint32_t timebase = 0) {
  return beat88(static_cast<uint16_t>(bpm << 8), timebase);
}
/**
 * @brief Free-running 8-bit sawtooth phase at the given whole BPM.
 * @param bpm Tempo in beats per minute.
 * @param timebase Millisecond offset for the zero of time.
 * @return The current phase in [0, 255].
 */
inline uint8_t beat8(uint16_t bpm, uint32_t timebase = 0) {
  return beat16(bpm, timebase) >> 8;
}

/**
 * @brief FastLED-faithful beatsin8: a sin8 oscillation at the given BPM.
 * @param bpm Tempo in beats per minute.
 * @param lowest Lower bound of the output, in [0, 255].
 * @param highest Upper bound of the output, in [0, 255].
 * @param timebase Millisecond offset for the zero of time.
 * @param phase_offset Phase shift added to the wave, in [0, 255].
 * @return An 8-bit value oscillating within [lowest, highest].
 * @pre lowest <= highest. The `highest - lowest` span is an unsigned subtraction,
 *      so passing lowest > highest underflows and escapes the documented range.
 *      This matches device <FastLED.h>, which subtracts identically and is not
 *      guarded either — clamping here would diverge the host from the device.
 * @details The parameter order and LUT match <FastLED.h>, and the scale8 range
 *          fit keeps the result within [lowest, highest].
 */
inline uint8_t beatsin8(uint16_t bpm, uint8_t lowest = 0, uint8_t highest = 255,
                        uint32_t timebase = 0, uint8_t phase_offset = 0) {
  uint8_t s = sin8(beat8(bpm, timebase) + phase_offset);
  return lowest + scale8(s, highest - lowest);
}

/**
 * @brief FastLED-faithful beatsin16: the 16-bit counterpart (sin16 + scale16).
 * @param bpm Tempo in beats per minute.
 * @param lowest Lower bound of the output, in [0, 65535].
 * @param highest Upper bound of the output, in [0, 65535].
 * @param timebase Millisecond offset for the zero of time.
 * @param phase_offset Phase shift added to the wave, in [0, 65535].
 * @return A 16-bit value oscillating within [lowest, highest].
 * @pre lowest <= highest. The `highest - lowest` span is an unsigned subtraction,
 *      so passing lowest > highest underflows and escapes the documented range,
 *      matching the unguarded device <FastLED.h> (see beatsin8).
 */
inline uint16_t beatsin16(uint16_t bpm, uint16_t lowest = 0,
                          uint16_t highest = 65535, uint32_t timebase = 0,
                          uint16_t phase_offset = 0) {
  uint16_t s = static_cast<uint16_t>(sin16(beat16(bpm, timebase) + phase_offset) +
                                     32768);
  return lowest + scale16(s, highest - lowest);
}

/**
 * @brief Modular 8-bit addition (FastLED addmod8).
 * @param a First addend.
 * @param b Second addend.
 * @param m Modulus.
 * @return (a + b) mod m, or (a + b) when m == 0.
 * @details A zero modulus SIGFPEs on the host (x86 integer div-by-zero) while
 *          the Cortex-M7 device returns the unreduced sum: UDIV-by-zero yields a
 *          0 quotient (traps off), so the remainder reduces to (a + b). Match the
 *          device rather than crashing only in the simulator, mirroring the
 *          m == 0 guards in map()/random8().
 */
inline uint8_t addmod8(uint8_t a, uint8_t b, uint8_t m) {
  if (m == 0) return static_cast<uint8_t>(a + b);
  return (a + b) % m;
}

/**
 * @brief Maps the full 0..255 input onto [rangeStart, rangeEnd] (FastLED map8).
 * @param in Input value in [0, 255].
 * @param rangeStart Lower bound of the output range.
 * @param rangeEnd Upper bound of the output range.
 * @return in scaled via scale8 onto [rangeStart, rangeEnd].
 * @pre rangeStart <= rangeEnd. The `rangeEnd - rangeStart` span is an unsigned
 *      subtraction, so passing rangeStart > rangeEnd underflows and escapes the
 *      documented range, matching the unguarded device <FastLED.h> (see beatsin8).
 * @details Maps the fixed 0..255 input, not a remap of an arbitrary input range.
 */
inline uint8_t map8(uint8_t in, uint8_t rangeStart, uint8_t rangeEnd) {
  return rangeStart + scale8(in, rangeEnd - rangeStart);
}

/**
 * @brief Symmetric triangle wave over an 8-bit input (FastLED triwave8).
 * @param in Input value in [0, 255].
 * @return A triangle wave value (0->0, 128->254).
 */
inline uint8_t triwave8(uint8_t in) {
  if (in & 0x80) {
    in = 255 - in;
  }
  return in << 1;
}

#define FASTRUN

// Mock EVERY_N_MILLIS using a simple static checker.
// Two-level macro for proper __LINE__ token pasting.
/**
 * @brief Pastes two tokens after expanding them (inner stage).
 * @param a Left token.
 * @param b Right token.
 */
#define HS_CONCAT_(a, b) a##b
/**
 * @brief Pastes two tokens, expanding macros such as __LINE__ first.
 * @param a Left token.
 * @param b Right token.
 */
#define HS_CONCAT(a, b) HS_CONCAT_(a, b)

/**
 * @brief Executes the guarded block at most once every N milliseconds.
 * @param N Interval in milliseconds.
 * @details Expands to a static throttle object plus one `if`, the SAME two-token
 * shape as the device's class-based FastLED macro (`static CEveryNMillis o(N);
 * if (o)`). The previous host mock expanded to THREE statements — a `static`
 * decl, an extra non-static local `__now`, then the `if` — so its body did not
 * line up statement-for-statement with the device's, the divergence this closes
 * (`EVERY_N_MILLIS host-mock parity`). As on the device, the trailing `if` is
 * what the following braced block attaches to; like FastLED's, this still cannot
 * serve as the *unbraced* body of an outer control statement (a leading `static`
 * decl is not a valid lone substatement) — that constraint now matches the
 * device instead of being a third host-only failure mode. See hs::EveryNMillis
 * for the preserved timing semantics.
 */
#define EVERY_N_MILLIS(N)                                                      \
  static hs::EveryNMillis HS_CONCAT(__every_, __LINE__)((N));                  \
  if (HS_CONCAT(__every_, __LINE__))

/**
 * @brief Executes the guarded block at most once every N seconds.
 * @param N Interval in seconds.
 */
#define EVERY_N_SECONDS(N) EVERY_N_MILLIS((N) * 1000)
/**
 * @brief Executes the guarded block at most once every N milliseconds (alias).
 * @param N Interval in milliseconds.
 */
#define EVERY_N_MILLISECONDS(N) EVERY_N_MILLIS(N)

namespace hs {
/**
 * @brief Returns milliseconds since an arbitrary epoch (host millis()).
 * @return Monotonic millisecond count, or the injected mock time when enabled.
 * @details Uses steady_clock (monotonic), matching micros(): a wall-clock
 *          source would let an NTP step or manual clock change make millis()
 *          jump or go backward, so the unsigned `now - last` in EVERY_N_MILLIS
 *          wraps huge and the beat/beatsin phases jump. The device millis() is
 *          monotonic; the simulator must match.
 */
inline unsigned long millis() {
  if (use_mock_time) return mock_millis_value;
  using namespace std::chrono;
  return duration_cast<milliseconds>(steady_clock::now().time_since_epoch())
      .count();
}

/**
 * @brief Host throttle backing EVERY_N_MILLIS, mirroring FastLED's CEveryNMillis.
 * @details Class-based like the device's FastLED macro so EVERY_N_MILLIS can
 * expand to a single guarded statement (`static EveryNMillis o(N); if (o)`)
 * instead of a multi-statement sequence, keeping the host structurally in step
 * with the device and nesting correctly in control flow.
 *
 * Timing is preserved exactly from the prior mock: `last_` starts at 0 so the
 * first evaluation fires (`now - 0 >= period`), and the trigger stamp is never
 * reset across effect switches (the object is a function-local `static`). The
 * period is captured at construction, matching the device's FastLED object
 * (which likewise fixes its period when the static is built).
 */
class EveryNMillis {
public:
  explicit EveryNMillis(unsigned long period) : last_(0), period_(period) {}

  /** @brief True at most once per `period_` ms; stamps the trigger when it fires. */
  bool ready() {
    unsigned long now = millis();
    if (now - last_ >= period_) {
      last_ = now;
      return true;
    }
    return false;
  }

  /** @brief Contextual-bool form so `if (obj)` reads as the throttle gate. */
  explicit operator bool() { return ready(); }

private:
  unsigned long last_;
  unsigned long period_;
};
/**
 * @brief Returns microseconds since an arbitrary epoch (host micros()).
 * @return Monotonic microsecond count, or the injected mock time when enabled.
 */
inline unsigned long micros() {
  if (use_mock_time) return mock_micros_value;
  using namespace std::chrono;
  return (unsigned long)duration_cast<microseconds>(
      steady_clock::now().time_since_epoch()).count();
}
/** @brief Disables interrupts (no-op on host). */
inline void disable_interrupts() {}
/** @brief Enables interrupts (no-op on host). */
inline void enable_interrupts() {}

// Global state
/** @brief Global debug-logging toggle. */
inline bool debug = false;
// rand_f/rand_int and ScanMetrics are platform-agnostic — defined once in the
// shared hs namespace after the #endif below.
#define HS_OS_CYCLES() 0

/**
 * @brief Virtual rows appended below the physical LED ring; 0 on host/sim.
 * @details The simulator has no physical LED ring to clip against, so it maps
 *          the full sphere. This intentionally diverges from the device's
 *          H_OFFSET = 3 (see the CORE_TEENSY definition above for the full
 *          rationale), so vertical sphere coverage is not bit-identical to
 *          hardware.
 *
 *          HS_TEST_H_OFFSET override: the south-pole Y-clip renormalization in
 *          the rasterizer (the bilinear-tap fold in Screen::AntiAlias::plot, and
 *          the H_VIRT handling threaded through geometry/scan/plot/sdf) only does
 *          real work when H_OFFSET > 0, which never happens on a normal host
 *          build. A dedicated host test executable defines HS_TEST_H_OFFSET=3 so
 *          the WHOLE pipeline compiles with the hardware offset and the renorm
 *          path runs against an energy-conservation oracle (see
 *          tests/test_h_offset_renorm.h, tests/CMakeLists.txt). This is the same
 *          recompile-under-device-config tactic as the fastmath_clamp_check TU.
 *          It MUST live in its own translation unit/executable: an offset-3 and
 *          an offset-0 instantiation of the same PhiLUT<H>/TrigLUT<W,H> would
 *          otherwise have different static-array sizes and clash under ODR.
 */
#if defined(HS_TEST_H_OFFSET)
static constexpr int H_OFFSET = HS_TEST_H_OFFSET;
#else
static constexpr int H_OFFSET = 0;
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
// Keep a setup-only function off the fast ITCM FlexRAM banks (RAM1). It carries
// FLASHMEM (the ldscript routes .flashmem to FLASH) plus noinline + noclone, and
// the latter two do the real work for cold operators that the optimizer would
// otherwise replicate across ITCM: noinline collapses the per-call-site inline
// copies a hot caller would stamp out (every solid factory inlining the same
// Conway operator), and noclone blocks the .constprop / .isra IPA clones — which
// drop the section attribute and so land in ITCM regardless, as
// classify_faces_by_topology did with FLASHMEM + noinline alone. Apply ONLY to
// internal-linkage (`static`) free functions on cold paths (mesh/solid
// construction): a section attribute on a COMDAT (inline/template member)
// function is a section-type conflict, and the per-frame render path must never
// pay flash latency. Off-device FLASHMEM is empty and noclone is GCC-only, so
// the macro degrades to a no-op for the host/simulator build.
// ---------------------------------------------------------------------------
#if defined(__GNUC__) && !defined(__clang__)
#define HS_COLD FLASHMEM __attribute__((noinline, noclone))
#else
#define HS_COLD FLASHMEM
#endif

// ---------------------------------------------------------------------------
// Platform-agnostic hs:: helpers (defined once; both branches above provide
// hs::random()). Hoisted out of the per-platform #if to remove duplication.
// ---------------------------------------------------------------------------
namespace hs {

/**
 * @brief Maps a raw RNG draw in [0, max] onto the half-open interval [0.0, 1.0).
 * @param value Raw RNG draw, in [0, max].
 * @param max Maximum possible draw value.
 * @return A float in [0.0, 1.0), clamped just below 1.0f at the top band.
 * @details The naive value/max can land on exactly 1.0f for the top band of
 *          draws because both operands — value and the divisor max, which is
 *          2^32-1 (UINT32_MAX), not 2^32 — round UP to 2^32 in float32 (2^32-1
 *          is not representable there), so (int)(u * N) would
 *          occasionally index N — one past the end. Clamp only those top draws
 *          to the float just below 1.0f; no float is representable between that
 *          constant and 1.0f, so every other draw is byte-for-byte unchanged and
 *          the stream stays deterministic. Pure so the boundary is unit-testable
 *          without driving the global RNG to its (rare) max.
 */
inline float random_to_unit(uint32_t value, uint32_t max) {
  float r = static_cast<float>(value) / static_cast<float>(max);
  constexpr float kJustBelowOne = 0x1.fffffep-1f; // nextafterf(1.0f, 0.0f)
  return r > kJustBelowOne ? kJustBelowOne : r;
}

/**
 * @brief Generates a pseudo-random floating-point number in [0.0, 1.0).
 * @return A random float in the half-open range [0.0, 1.0).
 */
inline float rand_f() {
  // The uint32_t casts below — and random_to_unit's 2^32-1 divisor — are exact
  // only while the global RNG draws fill the full 32-bit range. Trap a future
  // generator swap (e.g. to a 64-bit engine) that would silently truncate both
  // the draw and max() to uint32_t and break the [0,1) guarantee.
  using Rng = std::remove_reference_t<decltype(hs::random())>;
  static_assert(Rng::max() == 0xFFFFFFFFu,
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
 */
inline int rand_int(int min, int max) {
  if (max > min) {
    return min + (hs::random()() % (max - min));
  }
  return min;
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
                                    const char *cond, const char *fmt, ...)
    __attribute__((format(printf, 4, 5)));
[[noreturn]] inline void check_fail(const char *file, int line,
                                    const char *cond, const char *fmt, ...) {
  // Size the message buffer to the downstream 256-byte sink (the EM_ASM buf and
  // hs::log's internal buffer). A smaller intermediate would truncate a long
  // HS_CHECK message *before* the sink even sees it — a second, hidden choke on
  // top of the documented single bound, defeating the "no truncation" rationale.
  char msg[256];
  va_list args;
  va_start(args, fmt);
#ifdef ARDUINO
  // Integer-only formatter: keeps newlib's float path (_dtoa_r + the bignum
  // helpers behind %f/%g) out of ITCM, matching hs::log. Nothing reaches an
  // HS_CHECK with a float argument, so the device loses no diagnostics.
  vsniprintf(msg, sizeof(msg), fmt, args);
#else
  vsnprintf(msg, sizeof(msg), fmt, args);
#endif
  va_end(args);
  // Strip the directory from __FILE__: the compiler bakes in an absolute build
  // path that can be far longer than the basename, and on the bounded log buffer
  // below it would crowd out the cond/message tail — the part that actually
  // identifies the failure. Keep just the filename.
  const char *base = file;
  for (const char *p = file; *p; ++p) {
    if (*p == '/' || *p == '\\') base = p + 1;
  }
#ifdef __EMSCRIPTEN__
  // Route straight to console.error: synchronous (the trap can't drop it) and
  // visible even if stdout is not wired to the page console.
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
  uint32_t lut_hits = 0;      /**< Count of lookup-table hits. */
  uint32_t exact_hits = 0;    /**< Count of exact (non-LUT) computations. */
  /** @brief Zeroes every counter. */
  void reset() { plot = sdf_dist = frag_shader = bounds = face_setup = scan_loop = pixels_tested = pixels_culled = lut_hits = exact_hits = 0; }
};
/** @brief Global scanline profiling counters. Compiled in only when
 *  HS_SCAN_METRICS is defined — every reader goes through HS_SCAN_METRIC(...),
 *  which expands to nothing in the release build, so the global would otherwise
 *  be dead storage there. */
#ifdef HS_SCAN_METRICS
inline ScanMetrics g_scan_metrics;
#endif

} // namespace hs

// Per-pixel scan instrumentation is OFF by default. The README's engineering
// philosophy requires the per-pixel hot loop to stay lean, and a g_scan_metrics
// increment is a non-atomic global load-modify-store on a shared cache line for
// every pixel. Define HS_SCAN_METRICS to compile the counters back in — the
// native test build does, to assert which Face::distance path each sample took;
// the device/WASM release build leaves it undefined so HS_SCAN_METRIC(...)
// expands to nothing and the hot loop pays nothing.
#ifdef HS_SCAN_METRICS
#define HS_SCAN_METRIC(stmt) do { (stmt); } while (0)
#else
#define HS_SCAN_METRIC(stmt) ((void)0)
#endif

// ---------------------------------------------------------------------------
// Fn<Sig, Cap> — platform-aware callable wrapper. BOTH backends use heap-free
// inline storage, so a captured closure is never heap-allocated (which, stored in
// an ArenaVector that never destroys its elements, would leak under LSan).
//   Teensy:     teensy::inplace_function
//   Host/WASM:  hs::inplace_function
//
// Cap is a hard inline byte budget on both backends: a capture that overflows it
// is a compile error, not a heap allocation. Because Cap counts bytes, a closure
// that captures pointers is wider on the 64-bit host than on the 32-bit
// device/WASM. The rare callsite whose device-tuned Cap is too tight for 64-bit
// pointers sizes its Cap by sizeof(void*) rather than inflating every Fn here
// (which would push Fn-bearing types past TimelineEvent::MAX_ANIM_SIZE). See
// SpriteFn in concepts.h.
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
#include "inplace_function.h" // hs::inplace_function (declared after check_fail)
/**
 * @brief Platform-aware callable wrapper (host/WASM: heap-free inplace_function).
 * @tparam Sig Call signature, e.g. void(int).
 * @tparam Cap Inline storage capacity in bytes for the captured state.
 */
template <typename Sig, size_t Cap = 16>
using Fn = hs::inplace_function<Sig, Cap>;
#endif

// Detect x86 / x64 architecture (Desktop/Simulator)
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) ||             \
    defined(_M_IX86)
#include <xmmintrin.h> // Required for SSE intrinsics
#define HS_ARCH_X86
#endif

namespace hs {

// The clamp NaN->hi contract below is load-bearing for engine-wide float->int
// domain safety: Spherical, vector_to_pixel, blend_alpha, Gradient::get, and
// every palette lookup feed a possibly-NaN value through clamp as a saturating
// guard before a float->int cast. That guard survives only under real IEEE
// (non-finite) semantics. -ffinite-math-only — which a bare -ffast-math implies
// — licenses the compiler to assume no NaN/Inf exists and fold the guard away,
// silently reintroducing the cast UB across the whole renderer. The WASM release
// build keeps the contract solely by re-applying -fno-finite-math-only after
// -ffast-math (see CMakeLists.txt). Trap at compile time, on every target, if
// that protection is ever lost (flag reorder, toolchain default change, or an
// -ffinite-math-only sub-target) so the regression is a build error here rather
// than silent corruption on the sphere. The native fastmath clamp test
// (tests/CMakeLists.txt) exercises the contract under the real flags.
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ != 0
#error "hs::clamp NaN->hi contract requires -fno-finite-math-only: a bare -ffast-math (or -ffinite-math-only) makes the compiler assume no NaN and folds the saturating clamp guard away, reintroducing float->int cast UB engine-wide."
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
 *          SECOND source operand whenever either input is NaN, so min(NaN, hi)
 *          == hi and then max(lo, hi) == hi. This backend is therefore
 *          REORDER-SENSITIVE: swapping to min(hi, v) would make min(hi, NaN) ==
 *          NaN and break the contract — hence the do-not-reorder note. (The
 *          fminf backend below reaches the same hi by a different rule, IEEE
 *          NaN-suppression, and is reorder-insensitive; see its note.) Callers
 *          feed a possibly-NaN value through this as a saturating guard before a
 *          float->int cast (e.g. blend_alpha and Gradient::get; see
 *          test_blend_alpha_clamps_before_cast and
 *          test_gradient_get_clamps_out_of_range). Do not reorder the min
 *          operands.
 */
inline __attribute__((always_inline)) float clamp(float v, float lo, float hi) {
  // Load floats into the 128-bit SSE registers
  __m128 mv = _mm_set_ss(v);
  __m128 mlo = _mm_set_ss(lo);
  __m128 mhi = _mm_set_ss(hi);

  // Hardware execution: max(lo, min(v, hi))
  __m128 res = _mm_max_ss(mlo, _mm_min_ss(mv, mhi));

  // Extract the result back to a standard C++ float
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
 * @details On Cortex-M7 compiles directly to VMIN.F32 / VMAX.F32. The NaN
 *          contract holds for a DIFFERENT reason than the x86 backend: IEEE
 *          __builtin_fminf/fmaxf are NaN-SUPPRESSING — they return the non-NaN
 *          operand regardless of its position — so min(NaN, hi) == hi (and
 *          min(hi, NaN) would too), then max(lo, hi) == hi. This backend is thus
 *          REORDER-INSENSITIVE; the operand order is kept identical to the x86
 *          overload only for parity, not for correctness here. See the x86
 *          overload for the caller contract. (NaN-suppression relies on
 *          -fno-finite-math-only surviving after -ffast-math; the compile-time
 *          __FINITE_MATH_ONLY__ #error above guards exactly that.)
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
inline __attribute__((always_inline)) int clamp(int v, int lo, int hi) {
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
 * @brief Named cumulative cycle accumulator. Each instance self-registers into
 *        a static intrusive list at construction so log_all()/reset_all() can
 *        walk every counter without a central registry. Counters nest: a
 *        CycleScope sets `parent` to whichever counter was active when this one
 *        started, giving log_all() a call tree with per-parent percentages.
 * @warning REENTRANCY: the registry head and the `active_` nesting pointer are
 *        non-atomic statics (like hs::random()'s generator), so construction and
 *        CycleScope enter/exit are main-loop-only — driving a CycleScope from an
 *        ISR would race the list/active pointer and corrupt the call tree.
 */
struct CycleCounter {
  static constexpr uint32_t CYCLES_PER_US = 600; /**< Core clock: Teensy 4 @ 600 MHz. */

  const char* name;                /**< Counter label used in log output. */
  uint64_t cycles = 0;             /**< Accumulated cycle count. 64-bit because a
                                        32-bit accumulator overflows after only
                                        ~7 s of summed time at 600 MHz, which a
                                        multi-frame profiling run easily exceeds. */
  uint32_t count = 0;              /**< Number of timed invocations. */
  CycleCounter* parent = nullptr;  /**< Enclosing counter for tree nesting. */
  CycleCounter* next = nullptr;    /**< Next link in the intrusive registry list. */

  /**
   * @brief Constructs a named counter and self-registers it for bulk logging.
   * @param n Counter label (must outlive the counter; typically a literal).
   */
  explicit CycleCounter(const char* n) : name(n), next(head_) { head_ = this; }

  /** @brief Zeroes this counter's accumulated cycles and call count. */
  void reset() { cycles = 0; count = 0; }

  /** @brief Logs every root counter (no parent) and its subtree as a tree. */
  static void log_all() {
    hs::log("--- Cycle Counters ---");
    for (auto* c = head_; c; c = c->next)
      if (!c->parent && c->count) log_node(c, 0);
  }

  /** @brief Zeroes every registered counter (between profiling runs). */
  static void reset_all() {
    for (auto* c = head_; c; c = c->next)
      c->reset();
  }

private:
  static inline CycleCounter* head_ = nullptr;   /**< Head of the intrusive registry list. */
  static inline CycleCounter* active_ = nullptr; /**< Currently active counter (for nesting). */
  friend struct CycleScope;

  /**
   * @brief Recursively logs one counter node and its children as a tree.
   * @param node Counter node to log.
   * @param depth Tree depth; drives indentation.
   * @details The reported percentage is this node's cycles over its parent's
   *          (or 100% for a root), and cycles are converted to microseconds via
   *          CYCLES_PER_US.
   */
  static void log_node(const CycleCounter* node, int depth) {
    if (!node->count) return;
    uint64_t ref = node->parent ? node->parent->cycles : node->cycles;
    uint32_t pct = ref ? (uint32_t)(node->cycles * 100 / ref) : 100;
    uint64_t us = node->cycles / CYCLES_PER_US;
    int indent = depth * 2;
    int name_w = 22 - indent;
    if (name_w < 1) name_w = 1;
    hs::log("%*s%-*s %lu us (%lu%%)  %lu calls",
            indent, "", name_w, node->name,
            (unsigned long)us, (unsigned long)pct, (unsigned long)node->count);
    for (auto* c = head_; c; c = c->next)
      if (c->parent == node) log_node(c, depth + 1);
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
  CycleCounter& counter;       /**< Counter this scope accumulates into. */
  CycleCounter* prev_active;   /**< Counter to restore as active on destruction. */
  uint32_t start;              /**< Cycle snapshot taken at construction (32-bit,
                                    matching the hardware DWT CYCCNT register). */

  /**
   * @brief Begins timing the enclosing scope into the given counter.
   * @param c Counter that receives the elapsed cycles.
   * @details Makes c the active counter (recording the previously-active
   *          counter as its parent on first use) and snapshots the cycle
   *          counter.
   */
  explicit CycleScope(CycleCounter& c) : counter(c), start(HS_OS_CYCLES()) {
    prev_active = CycleCounter::active_;
    if (!counter.parent && prev_active)
      counter.parent = prev_active;
    CycleCounter::active_ = &counter;
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
    CycleCounter::active_ = prev_active;
  }

  /**
   * @brief Deleted copy constructor; a scope guard must not be copied.
   */
  CycleScope(const CycleScope&) = delete;
  /**
   * @brief Deleted copy assignment; a scope guard must not be copied.
   * @return Never returns; deleted.
   */
  CycleScope& operator=(const CycleScope&) = delete;
};

} // namespace hs

/**
 * @brief Times the enclosing scope into a named cycle counter.
 * @param label Counter name (used both as the identifier suffix and log label).
 */
#define HS_PROFILE(label) \
  static hs::CycleCounter hs_ctr_##label(#label); \
  hs::CycleScope hs_scope_##label(hs_ctr_##label)

