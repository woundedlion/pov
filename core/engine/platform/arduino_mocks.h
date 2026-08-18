/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file arduino_mocks.h
 * @brief Emulates the Arduino/FastLED API surface off-device.
 *
 * Included from the host/sim branch of engine/platform.h, which supplies
 * hs::random/hs::millis.
 */

#include <cstdarg>
#include <cstdint>
#include <cstdio>

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
 *          effects compile and behave identically to the device — the lone
 *          exception is CRGB(const CHSV&), which is not bit-identical (see its
 *          @warning).
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
   * @warning NOT bit-identical to the device: FastLED's runtime path is
   *          hsv2rgb_rainbow, while this is a 6-sector integer spectrum giving
   *          visibly different RGB. A legacy-compat conversion only — do not use
   *          on parity-sensitive paths (no modern effect uses CHSV).
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

/** @brief Mock LED chipset selector for addLeds template arguments. */
enum LEDType { WS2801 };
/** @brief Mock LED color-order selector for addLeds template arguments. */
enum ColorOrder { RGB };
#define DATA_RATE_MHZ(x) (x)

/**
 * @brief Mock implementation of the FastLED controller for simulation.
 * @details The simulator renders into its own framebuffer rather than driving
 *          real LEDs, so every method is a no-op except setCorrection and
 *          setTemperature, which record their selector.
 */
struct FastLEDMock {
  /** @brief Selector last passed to setCorrection; -1 before any call. */
  int last_correction = -1;
  /** @brief Selector last passed to setTemperature; -1 before any call. */
  int last_temperature = -1;
  /**
   * @brief Records the color-correction profile.
   * @param correction Correction selector, stored in last_correction.
   */
  void setCorrection(int correction) { last_correction = correction; }
  /**
   * @brief Records the color-temperature profile.
   * @param temperature Temperature selector, stored in last_temperature.
   */
  void setTemperature(int temperature) { last_temperature = temperature; }
  /**
   * @brief Registers an LED strip (no-op on host).
   * @tparam CHIPSET LED chipset selector.
   * @tparam DATA_PIN Data pin number.
   * @tparam CLOCK_PIN Clock pin number.
   * @tparam RGB_ORDER Channel order the chipset expects.
   * @tparam SPI_DATA_RATE SPI clock rate in Hz.
   * @details The parameter kinds mirror FastLED's SPI overload, whose chipset
   *          and order arguments are enumerators rather than types.
   */
  template <LEDType CHIPSET, uint8_t DATA_PIN, uint8_t CLOCK_PIN,
            ColorOrder RGB_ORDER, uint32_t SPI_DATA_RATE>
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

// --- Mock Arduino Functions ---
// Global scope mirrors Arduino/FastLED, which expose random()/map() as free
// globals; unqualified callers (e.g. the legacy effects) must resolve
// identically on host and device.
/**
 * @brief Returns a pseudo-random integer in [0, max) (Arduino random()).
 * @param max Exclusive upper bound.
 * @return A value in [0, max), or 0 when max <= 0.
 * @details Guards a degenerate range like the device: Arduino's random(0)
 *          returns 0, so the host avoids a modulo-by-zero SIGFPE and matches.
 */
inline int random(int max) {
  if (max <= 0)
    return 0;
  return hs::random()() % max;
}
/**
 * @brief Returns a pseudo-random integer in [min, max) (Arduino random()).
 * @param min Inclusive lower bound.
 * @param max Exclusive upper bound.
 * @return A value in [min, max), or min when max <= min.
 * @details Mirrors the device: random(min, max) with min >= max returns min,
 *          so the host avoids a modulo-by-zero SIGFPE.
 * @note The span is computed unsigned: `max - min` overflows int for a span
 *       wider than INT_MAX.
 */
inline int random(int min, int max) {
  if (max <= min)
    return min;
  const uint32_t span = static_cast<uint32_t>(max) - static_cast<uint32_t>(min);
  return static_cast<int>(static_cast<uint32_t>(min) + (hs::random()() % span));
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
 * @details Degenerate input range (in_max == in_min) returns out_min: Arduino's
 *          map() divides with no guard, which SIGFPEs on the host while the
 *          Cortex-M7 returns 0 from the divide (relies on CCR.DIV_0_TRP staying
 *          at its clear reset default). Match the device rather than crashing
 *          only in the simulator.
 */
inline long map(long x, long in_min, long in_max, long out_min, long out_max) {
  // Device computes in 32-bit `long`. Multiply in uint32_t (defined wrap mod
  // 2^32) and reinterpret to int32_t to reproduce its two's-complement
  // truncation without 64-bit widening (LP64) or signed-overflow UB.
  const int32_t divisor = static_cast<int32_t>(in_max - in_min);
  if (divisor == 0)
    return out_min;
  const int32_t product =
      static_cast<int32_t>(static_cast<uint32_t>(x - in_min) *
                           static_cast<uint32_t>(out_max - out_min));
  // INT32_MIN / -1 overflows (host UB) but ARM SDIV returns INT32_MIN; negate
  // mod 2^32 to reproduce the device result without trapping.
  const int32_t scaled =
      divisor == -1 ? static_cast<int32_t>(0u - static_cast<uint32_t>(product))
                    : product / divisor;
  return static_cast<int32_t>(static_cast<uint32_t>(scaled) +
                              static_cast<uint32_t>(out_min));
}

// --- System Mock ---
/**
 * @brief Mock of Arduino's Serial that writes to stdout on the host.
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
  void print(const char *msg) { fputs(msg, stdout); }
  /**
   * @brief Writes a signed integer without a trailing newline.
   * @param val Value to write.
   */
  void print(int val) { printf("%d", val); }
  /**
   * @brief Writes an unsigned long without a trailing newline.
   * @param val Value to write.
   */
  void print(unsigned long val) { printf("%lu", val); }
  /**
   * @brief Writes a float without a trailing newline.
   * @param val Value to write.
   */
  void print(float val) { printf("%g", static_cast<double>(val)); }
  /**
   * @brief Writes a C string followed by a newline.
   * @param msg Text to write.
   */
  void println(const char *msg) { printf("%s\n", msg); }
  /**
   * @brief Writes a signed integer followed by a newline.
   * @param val Value to write.
   */
  void println(int val) { printf("%d\n", val); }
  /**
   * @brief Writes an unsigned long followed by a newline.
   * @param val Value to write.
   */
  void println(unsigned long val) { printf("%lu\n", val); }
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
    fputs(buf, stdout);
  }
};
inline SerialMock Serial;

// --- FastLED Mocks ---
// These route through hs::random() (Pcg32) and do NOT match the device, where
// random8/random16 are FastLED's LCG; used only by legacy effects (modern
// effects use hs::rand_*).
/**
 * @brief Returns a pseudo-random 8-bit value (FastLED random8).
 * @return A value in [0, 255].
 */
inline uint8_t random8() { return hs::random()() % 256; }
/**
 * @brief Returns a pseudo-random 8-bit value below a limit (FastLED random8).
 * @param top Exclusive upper bound.
 * @return A value in [0, top), or 0 when top == 0.
 * @details Host modulo, not the device's scaled multiply ((r*lim)>>8): like the
 *          other FastLED mocks here this is legacy-only and not bit-faithful (see
 *          the note above). The top==0 guard returns 0 to match the device's
 *          scaled form (which yields 0) and to avoid the host modulo's SIGFPE.
 */
inline uint8_t random8(uint8_t top) {
  if (top == 0)
    return 0;
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

// FastLED fixed-point sine + scaling primitives. The device pulls these from
// <FastLED.h>; the host mock reproduces their exact integer semantics (LUT sine,
// 8.8 beat sawtooth, scale8 range fit) so the simulator predicts the device
// rather than approximating with a float sine.

/**
 * @brief Unsigned 8-bit fractional scale, scale8(i, sc) = i * (1 + sc) / 256.
 * @param i Value to scale, in [0, 255].
 * @param sc Scale factor, in [0, 255].
 * @return i scaled by sc/256 in the SCALE8_FIXED sense, in [0, 255].
 * @details The (1 + sc) is FastLED's SCALE8_FIXED form, so scale8(x, 255) == x
 *          (a full-scale fade is the identity). Matching it keeps the simulator
 *          bit-exact rather than 1 LSB low on every fade.
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
  if (theta & 0x40)
    offset = 255 - offset;
  offset &= 0x3F;
  uint8_t secoffset = offset & 0x0F;
  if (theta & 0x40)
    secoffset++;
  uint8_t section = offset >> 4;
  const uint8_t *p = b_m16_interleave + section * 2;
  uint8_t b = *p++;
  uint8_t m16 = *p;
  uint8_t mx = (m16 * secoffset) >> 4;
  int8_t y = mx + b;
  if (theta & 0x80)
    y = -y;
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
  if (theta & 0x4000)
    offset = 2047 - offset;
  uint8_t section = offset / 256; // 0..7
  uint16_t b = base[section];
  uint8_t m = slope[section];
  uint8_t secoffset8 = static_cast<uint8_t>(offset) / 2;
  int16_t y = static_cast<int16_t>(m * secoffset8) + b;
  if (theta & 0x8000)
    y = -y;
  return y;
}

/**
 * @brief Free-running 16-bit sawtooth phase at the given 8.8 fixed-point BPM.
 * @param bpm88 Tempo as an 8.8 fixed-point beats-per-minute value.
 * @param timebase Millisecond offset for the zero of time.
 * @return The current phase in [0, 65535].
 * @details Sourced from hs::millis() so the test time-injection seam keeps beats
 *          deterministic. The * 280 constant is FastLED's ms->phase scale
 *          (~65536 * 1000 / 60000). The `unsigned long` intermediate wraps mod
 *          2^32 on the device but not on a LP64 host; harmless because only bits
 *          16..31 (the `>>16` result) are returned, and those are identical
 *          either way.
 */
inline uint16_t beat88(uint16_t bpm88, uint32_t timebase = 0) {
  return ((hs::millis() - timebase) * bpm88 * 280) >> 16;
}
/**
 * @brief Free-running 16-bit sawtooth phase at the given whole BPM.
 * @param bpm Tempo in beats per minute.
 * @param timebase Millisecond offset for the zero of time.
 * @return The current phase in [0, 65535].
 * @details Matches FastLED: a value below 256 is promoted to 8.8 fixed point,
 *          256 and above is taken as already 8.8.
 */
inline uint16_t beat16(uint16_t bpm, uint32_t timebase = 0) {
  return beat88(bpm < 256 ? static_cast<uint16_t>(bpm << 8) : bpm, timebase);
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
 *      so passing lowest > highest underflows and escapes the documented range,
 *      matching the unguarded device <FastLED.h>.
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
  uint16_t s = static_cast<uint16_t>(
      sin16(beat16(bpm, timebase) + phase_offset) + 32768);
  return lowest + scale16(s, highest - lowest);
}

/**
 * @brief Modular 8-bit addition (FastLED addmod8).
 * @param a First addend.
 * @param b Second addend.
 * @param m Modulus.
 * @return The 8-bit-wrapped sum reduced mod m, or that sum when m == 0.
 * @details FastLED reduces the uint8_t-wrapped sum by repeated subtraction, so
 *          a sum past 255 wraps before the reduction: addmod8(200, 100, 7) is
 *          44 % 7, not 300 % 7. Its m == 0 loop subtracts zero forever on the
 *          device; the host returns the wrapped sum instead of hanging.
 */
inline uint8_t addmod8(uint8_t a, uint8_t b, uint8_t m) {
  uint8_t s = static_cast<uint8_t>(a + b);
  if (m == 0)
    return s;
  return static_cast<uint8_t>(s % m);
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
