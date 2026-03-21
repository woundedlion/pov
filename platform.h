/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <random>

#ifdef ARDUINO
#include <Arduino.h>
#include <FastLED.h>
#include <cstdarg>

namespace hs {
/**
 * @brief On-device logging to Serial (no vsnprintf to avoid ~4KB stdio in ITCM).
 */
inline void log(const char *msg, ...) {
  Serial.println(msg);
}

/**
 * @brief Returns the global deterministic random number generator.
 */
inline std::mt19937& random() {
  static std::mt19937 gen(1337);
  return gen;
}
/** @brief Wrapped millis() for namespace consistency. */
inline unsigned long millis() { return ::millis(); }
/** @brief Disables interrupts (Arduino). */
inline void disable_interrupts() { noInterrupts(); }
/** @brief Enables interrupts (Arduino). */
inline void enable_interrupts() { interrupts(); }

/**
 * @brief Generates a pseudo-random floating-point number between 0.0 and 1.0.
 * @return A random float in the range [0.0, 1.0].
 */
inline float rand_f() {
  return static_cast<float>(::random(0, std::numeric_limits<int32_t>::max())) /
         std::numeric_limits<int32_t>::max();
}

inline float rand_f(float min, float max) {
  return min + rand_f() * (max - min);
}

/**
 * @brief Generates a pseudo-random integer within a specified range.
 * @param min The minimum value (inclusive).
 * @param max The maximum value (exclusive).
 * @return A random integer in the range [min, max).
 */
inline int rand_int(int min, int max) { return ::random(min, max); }

// Global state
inline bool debug = false;

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

template <typename T, typename U> constexpr T lerp(T a, T b, U t) {
  return (T)(a + (b - a) * t);
}

/**
 * @brief HSV Color structure for non-Arduino platforms.
 */
struct CHSV {
  uint8_t h, s, v;
  constexpr CHSV() : h(0), s(0), v(0) {}
  constexpr CHSV(uint8_t h, uint8_t s, uint8_t v) : h(h), s(s), v(v) {}
};

// --- Mock FastLED Types ---
/**
 * @brief RGB Color structure mimicking FastLED's CRGB.
 */
struct CRGB {
  uint8_t r, g, b;
  constexpr CRGB() : r(0), g(0), b(0) {}
  constexpr CRGB(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
  constexpr CRGB(uint8_t gray) : r(gray), g(gray), b(gray) {}

  // Convert HSV to RGB (Basic implementation)
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

  // Operators
  bool operator==(const CRGB &rhs) const {
    return r == rhs.r && g == rhs.g && b == rhs.b;
  }
  bool operator!=(const CRGB &rhs) const { return !(*this == rhs); }

  CRGB &operator+=(const CRGB &rhs) {
    // Saturated add
    r = (r + rhs.r > 255) ? 255 : r + rhs.r;
    g = (g + rhs.g > 255) ? 255 : g + rhs.g;
    b = (b + rhs.b > 255) ? 255 : b + rhs.b;
    return *this;
  }

  // Methods matching FastLED
  CRGB lerp16(const CRGB &other, uint16_t frac) const {
    CRGB ret;
    // frac is 0..65535
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
// qadd8: saturated addition for 8-bit integers
inline uint8_t qadd8(uint8_t i, uint8_t j) {
  int t = i + j;
  if (t > 255)
    t = 255;
  return t;
}

// qsub8: saturated subtraction
inline uint8_t qsub8(uint8_t i, uint8_t j) {
  int t = i - j;
  if (t < 0)
    t = 0;
  return t;
}

#define PI 3.1415926535897932384626433832795

namespace hs {
inline void log(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  printf("\n");
}

/**
 * @brief Returns the global deterministic random number generator.
 */
inline std::mt19937& random() {
  static std::mt19937 gen(1337);
  return gen;
}
} // namespace hs

// --- Mock Arduino Constants/Types ---
using byte = uint8_t;
using boolean = bool;

// --- Mock FastLED Constants ---
enum FastLEDCheck {
  UncorrectedColor,
  TypicalLEDStrip,
  UncorrectedTemperature,
  Candle
};

/**
 * @brief Mock implementation of FastLED controller for simulation.
 */
struct FastLEDMock {
  void setCorrection(int) {}
  void setTemperature(int) {}
  template <typename T, int P1, int P2, int P3, int P4>
  void addLeds(CRGB *data, int nLeds) {}
  void show() {}
  void showColor(const CRGB &) {}
};
static FastLEDMock FastLED;

// Helper for addLeds template args
enum LEDType { WS2801 };
enum ColorOrder { RGB };
#define DATA_RATE_MHZ(x) x

// --- Mock Arduino Functions ---
inline int random(int max) { return hs::random()() % max; }
inline int random(int min, int max) { return min + (hs::random()() % (max - min)); }
inline long map(long x, long in_min, long in_max, long out_min, long out_max) {
  return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

// --- System Mock ---
struct SerialMock {
  void println(const char *msg) { std::cout << msg << std::endl; }
  void println(int val) { std::cout << val << std::endl; }
  void printf(const char *fmt, ...) {
    // simplified
    std::cout << fmt << std::endl;
  }
};
static SerialMock Serial;

// --- FastLED Mocks ---
inline uint8_t random8() { return hs::random()() % 256; }
inline uint8_t random8(uint8_t top) { return hs::random()() % top; }
inline uint16_t random16() { return hs::random()() % 65536; }
inline void random16_add_entropy(uint16_t) {}

/**
 * @brief Generates a sine wave beat using system clock.
 * @details Mocks FastLED's beatsin8: oscillates between low and high at
 * the given BPM. Phase is a 16-bit offset into the wave cycle.
 */
inline uint8_t beatsin8(uint16_t bpm, uint8_t low, uint8_t high,
                        uint16_t phase = 0, uint16_t offset = 0) {
  using namespace std::chrono;
  float t = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() / 1000.0f;
  float beats = t * bpm / 60.0f;
  float wave = (sinf(beats * 2.0f * PI + phase / 65536.0f * 2.0f * PI) + 1.0f) * 0.5f;
  return low + static_cast<uint8_t>(wave * (high - low));
}

/**
 * @brief 16-bit version of beatsin8.
 */
inline uint16_t beatsin16(uint16_t bpm, uint16_t low, uint16_t high,
                          uint16_t phase = 0, uint16_t offset = 0) {
  using namespace std::chrono;
  float t = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() / 1000.0f;
  float beats = t * bpm / 60.0f;
  float wave = (sinf(beats * 2.0f * PI + phase / 65536.0f * 2.0f * PI) + 1.0f) * 0.5f;
  return low + static_cast<uint16_t>(wave * (high - low));
}

inline uint8_t addmod8(uint8_t a, uint8_t b, uint8_t m) { return (a + b) % m; }

inline uint8_t map8(uint8_t x, uint8_t in_min, uint8_t in_max) {
  return (x - in_min) * 255 / (in_max - in_min); // rough approximation
}

inline uint8_t triwave8(uint8_t in) {
  if (in & 0x80) {
    in = 255 - in;
  }
  return in << 1;
}

#define FASTRUN

// Mock EVERY_N_MILLIS using a simple static checker
/**
 * @brief Macro to execute a block of code periodically.
 * @param N Interval in milliseconds.
 */
#define EVERY_N_MILLIS(N)                                                      \
  static unsigned long __last_##__LINE__ = 0;                                  \
  unsigned long __now_##__LINE__ = hs::millis();                               \
  if (__now_##__LINE__ - __last_##__LINE__ >= (N) &&                           \
      (__last_##__LINE__ = __now_##__LINE__))

#define EVERY_N_SECONDS(N) EVERY_N_MILLIS((N) * 1000)
#define EVERY_N_MILLISECONDS(N) EVERY_N_MILLIS(N)

namespace hs {
inline unsigned long millis() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch())
      .count();
}
inline void disable_interrupts() {}
inline void enable_interrupts() {}

/**
 * @brief Generates a pseudo-random floating-point number between 0.0 and 1.0.
 * @return A random float in the range [0.0, 1.0].
 */
inline float rand_f() {
  return static_cast<float>(hs::random()()) / static_cast<float>(hs::random().max());
}

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

// Global state
inline bool debug = false;

static constexpr int H_OFFSET = 0;
} // namespace hs

// Global millis/micros if needed, though prefer namespaced
inline unsigned long millis() { return hs::millis(); }
inline unsigned long micros() { return hs::millis() * 1000; } // approx

#endif

// ---------------------------------------------------------------------------
// Fn<Sig, Cap> — platform-aware callable wrapper
//   Teensy:  teensy::inplace_function  (no heap, tiny codegen)
//   WASM:    std::function             (Cap ignored)
// ---------------------------------------------------------------------------
#include <functional>
#ifdef ARDUINO
#include <inplace_function.h>
template <typename Sig, size_t Cap = 16>
using Fn = teensy::inplace_function<Sig, Cap>;
#else
template <typename Sig, size_t Cap = 16> using Fn = std::function<Sig>;
#endif

// Detect x86 / x64 architecture (Desktop/Simulator)
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) ||             \
    defined(_M_IX86)
#include <xmmintrin.h> // Required for SSE intrinsics
#define HS_ARCH_X86
#endif

namespace hs {

#ifdef HS_ARCH_X86
// --- x86 / x64 EXPLICIT HARDWARE CLAMP ---
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
// --- ARM CORTEX-M7 (TEENSY) HARDWARE CLAMP ---
inline constexpr __attribute__((always_inline)) float clamp(float v, float lo,
                                                            float hi) {
  // Compiles directly to VMIN.F32 and VMAX.F32
  return __builtin_fmaxf(lo, __builtin_fminf(v, hi));
}
#endif

// Standard integer fallback
inline __attribute__((always_inline)) int clamp(int v, int lo, int hi) {
  return (v < lo) ? lo : ((v > hi) ? hi : v);
}

} // namespace hs