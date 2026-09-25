/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file pixel.h
 * @brief Linear-light pixels, alpha, and integer sRGB conversion.
 */

#include <algorithm>
#include "platform/platform.h"
#include "color/srgb_decode.h"

/** @brief Rounds and saturates a linear-light channel to [0, 65535]. */
__attribute__((always_inline)) inline uint16_t
round_linear_channel(float value) {
  return static_cast<uint16_t>(hs::clamp(value + 0.5f, 0.0f, 65535.0f));
}

#if defined(__ARM_FEATURE_DSP)
// Inline assembly avoids a CMSIS header dependency for the saturating add.
__attribute__((always_inline)) inline uint32_t inline_uqadd16(uint32_t a,
                                                              uint32_t b) {
  uint32_t res;
  __asm__ volatile("uqadd16 %0, %1, %2" : "=r"(res) : "r"(a), "r"(b));
  return res;
}
#else
// Portable software model of ARM `uqadd16`: two independent 16-bit unsigned
// saturating adds, one per halfword lane.
inline uint32_t inline_uqadd16(uint32_t a, uint32_t b) {
  uint32_t lo = (a & 0xFFFFu) + (b & 0xFFFFu);
  uint32_t hi = (a >> 16) + (b >> 16);
  if (lo > 0xFFFFu)
    lo = 0xFFFFu;
  if (hi > 0xFFFFu)
    hi = 0xFFFFu;
  return (hi << 16) | lo;
}
#endif

struct Pixel;
// Saturating per-channel add packed into two uqadd16 lanes (g|b in one 32-bit
// word, r alone in another). Used by Pixel::operator+=.
inline Pixel pixel_blend_add_packed(const Pixel &c1, const Pixel &c2);

/**
 * @brief Maps an 8-bit sRGB channel value to its 16-bit linear equivalent.
 * @param srgb sRGB channel value in [0, 255].
 * @return 16-bit linear channel value.
 */
inline uint16_t srgb_to_linear(uint8_t srgb);

/**
 * @brief Represents a 16-bit per channel RGB color (Linear space).
 * @details Used for high-precision mixing and HDR rendering before
 * downsampling/tone-mapping to 8-bit output.
 */
struct Pixel {
  uint16_t r, g, b;

  /**
   * @brief Constructs a black pixel (all channels zero).
   */
  constexpr Pixel() : r(0), g(0), b(0) {}
  /**
   * @brief Constructs a pixel from explicit 16-bit linear channels.
   * @param r Red channel in [0, 65535].
   * @param g Green channel in [0, 65535].
   * @param b Blue channel in [0, 65535].
   */
  constexpr Pixel(uint16_t r, uint16_t g, uint16_t b) : r(r), g(g), b(b) {}

  /**
   * @brief Constructs a pixel from HSV (converts to sRGB then Linear).
   * @param hsv Source color in HSV space.
   */
  Pixel(const CHSV &hsv) {
    CRGB srgb(hsv);
    r = srgb_to_linear(srgb.r);
    g = srgb_to_linear(srgb.g);
    b = srgb_to_linear(srgb.b);
  }

  /**
   * @brief Constructs a pixel from CRGB (converts to Linear).
   * @param c Source color in 8-bit sRGB space.
   */
  Pixel(const CRGB &c) {
    r = srgb_to_linear(c.r);
    g = srgb_to_linear(c.g);
    b = srgb_to_linear(c.b);
  }

  /**
   * @brief Lossy 16-bit-linear -> 8-bit-sRGB downcast.
   * @return The color quantized to an 8-bit sRGB CRGB.
   * @details Explicit so a stray Pixel in a CRGB context is a compile error,
   * not a silent round-trip through 8-bit gamma.
   */
  explicit operator CRGB() const;

  /**
   * @brief Saturated per-channel addition into this pixel.
   * @param rhs Pixel to add.
   * @return Reference to this pixel after the clamped add.
   */
  Pixel &operator+=(const Pixel &rhs) {
#if defined(__ARM_FEATURE_DSP)
    *this = pixel_blend_add_packed(*this, rhs);
#else
    r = (uint16_t)std::min((uint32_t)65535, (uint32_t)r + rhs.r);
    g = (uint16_t)std::min((uint32_t)65535, (uint32_t)g + rhs.g);
    b = (uint16_t)std::min((uint32_t)65535, (uint32_t)b + rhs.b);
#endif
    return *this;
  }

  /**
   * @brief Saturated per-channel sum of two pixels.
   * @param rhs Pixel to add.
   * @return A new pixel with each channel clamped to the 16-bit max.
   */
  Pixel operator+(const Pixel &rhs) const {
    Pixel out = *this;
    out += rhs;
    return out;
  }

  /**
   * @brief Scales every channel by a float factor (saturated).
   * @param s Scale factor; may be any finite float (NaN maps to the hi bound).
   * @return A new pixel with each channel clamped to [0, 65535].
   * @details Rounds to nearest (+0.5f, inside the clamp so the hi bound stays
   * exactly 65535). Clamps in float before the cast: r*s can exceed INT_MAX and
   * float->int is UB out of range; hs::clamp also maps a NaN scale to the hi
   * bound before it can reach the cast.
   */
  Pixel operator*(float s) const {
    return Pixel(round_linear_channel(r * s), round_linear_channel(g * s),
                 round_linear_channel(b * s));
  }

  /**
   * @brief Linearly interpolates 16-bit between this pixel and another.
   * @param other Target pixel at frac == 65535.
   * @param frac Blend weight in [0, 65535]; 0 yields this pixel, 65535 yields other.
   * @return The interpolated pixel, round-to-nearest per channel.
   * @details Round-to-nearest div-by-65535 via shifts:
   *   (x + (x>>16) + 32768) >> 16, within 1 LSB of round(x/65535) and exact at
   * the endpoints (frac 0/65535 -> a/b). Plain 32-bit MACs, not packed `smlad`:
   * smlad's signed 16x16 dual-MAC reads an operand >= 32768 as negative.
   */
  __attribute__((always_inline)) Pixel lerp16(const Pixel &other,
                                              uint16_t frac) const {
    uint16_t inv = 65535 - frac;
    uint32_t xr = (uint32_t)r * inv + (uint32_t)other.r * frac;
    uint32_t xg = (uint32_t)g * inv + (uint32_t)other.g * frac;
    uint32_t xb = (uint32_t)b * inv + (uint32_t)other.b * frac;
    uint32_t r32 = (xr + (xr >> 16) + 32768) >> 16;
    uint32_t g32 = (xg + (xg >> 16) + 32768) >> 16;
    uint32_t b32 = (xb + (xb >> 16) + 32768) >> 16;
    return Pixel((uint16_t)r32, (uint16_t)g32, (uint16_t)b32);
  }

  /**
   * @brief Tests two pixels for exact channel equality.
   * @param rhs Pixel to compare against.
   * @return True if all three channels match.
   */
  bool operator==(const Pixel &rhs) const {
    return r == rhs.r && g == rhs.g && b == rhs.b;
  }

  /**
   * @brief Tests equality against an HSV color (converted to Pixel).
   * @param rhs Color in HSV space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CHSV &rhs) const { return *this == Pixel(rhs); }

  /**
   * @brief Tests equality against a CRGB color (converted to Pixel).
   * @param rhs Color in 8-bit sRGB space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CRGB &rhs) const { return *this == Pixel(rhs); }
};

/**
 * @brief Quantizes a [0,1] interpolation fraction to a 16-bit lerp16 weight.
 * @param frac Blend fraction; clamped to [0, 1].
 * @return The fraction as a 16-bit weight in [0, 65535], rounded.
 */
__attribute__((always_inline)) inline uint16_t frac_to_q16(float frac) {
  return static_cast<uint16_t>(hs::clamp(frac, 0.0f, 1.0f) * 65535.0f + 0.5f);
}

/**
 * @brief The interpolatable pixel of a lookup-table entry.
 */
__attribute__((always_inline)) inline Pixel lut_entry_pixel(const Pixel &e) {
  return e;
}

/**
 * @brief Master-alpha gate: one 8-bit LSB of the user's opacity slider.
 * @details A whole-effect gate on the slider value, not a per-sample floor —
 * the framebuffer is 16-bit linear and sRGB's toe lifts a fragment at this
 * alpha to roughly 13 encoded levels. Per-sample cuts use
 * MIN_ENCODABLE_ALPHA.
 */
inline constexpr float MIN_VISIBLE_ALPHA = 1.0f / 255.0f;

/**
 * @brief Per-sample alpha floor: below this a full-scale fragment encodes to 0.
 * @details linear_to_srgb8 first leaves zero at linear channel 10 of 65535
 * (sRGB's 12.92x toe puts the 0.5/255 rounding step at 1/(510*12.92) linear),
 * so a fragment whose premultiplied peak `max_channel(color) * alpha` stays
 * under 10/65535 cannot change any output pixel. Pinned against the decode
 * table by unit_color's test_min_encodable_alpha_is_the_encode_floor.
 */
inline constexpr float MIN_ENCODABLE_ALPHA = 10.0f / 65535.0f;

/**
 * @brief Represents a color with a STRAIGHT (non-premultiplied) alpha channel.
 * @details `color` holds the un-premultiplied color, `alpha` its coverage;
 * premultiplication happens once, at the final canvas write (`color * alpha`).
 */
struct Color4 {
  Pixel color;
  float alpha;

  /**
   * @brief Constructs a transparent black color (alpha 0.0).
   */
  Color4() : color(Pixel(0, 0, 0)), alpha(0.0f) {}
  /**
   * @brief Constructs a color from a Pixel and alpha.
   * @param p Linear-space pixel color.
   * @param a Alpha in [0, 1]; defaults to fully opaque.
   */
  Color4(Pixel p, float a = 1.0f) : color(p), alpha(a) {}
  /**
   * @brief Constructs a color from 8-bit sRGB channels and alpha.
   * @param r Red channel in [0, 255].
   * @param g Green channel in [0, 255].
   * @param b Blue channel in [0, 255].
   * @param a Alpha in [0, 1]; defaults to fully opaque.
   * @details `explicit` so the sRGB->linear convention is opt-in, not taken by
   *          a braced `{r,g,b}` from a caller modeling Color4 as already-linear.
   */
  explicit Color4(uint8_t r, uint8_t g, uint8_t b, float a = 1.0f)
      : color(Pixel(srgb_to_linear(r), srgb_to_linear(g), srgb_to_linear(b))),
        alpha(a) {}
  /**
   * @brief Constructs a color reusing another's pixel with a new alpha.
   * @param c Source color whose pixel is copied.
   * @param a Alpha in [0, 1] to apply.
   */
  Color4(const Color4 &c, float a) : color(c.color), alpha(a) {}

  /**
   * @brief Interpolates color (16-bit linear) and alpha by t.
   * @param other Target color at t == 1.
   * @param t Blend weight; clamped to [0, 1].
   * @return The interpolated color.
   * @details Alpha is straight, not premultiplied: color and alpha interpolate
   * independently, so a fully transparent endpoint still contributes its RGB to
   * the blend. Endpoints intended to fade out must carry the color they fade
   * towards. t clamped to [0,1] so out-of-range t saturates at an endpoint
   * rather than letting alpha extrapolate while color stays clamped.
   */
  Color4 lerp(const Color4 &other, float t) const {
    const float ct = hs::clamp(t, 0.0f, 1.0f);
    uint16_t frac = frac_to_q16(ct);
    Pixel blended = color.lerp16(other.color, frac);
    float blended_a = alpha + (other.alpha - alpha) * ct;
    return Color4(blended, blended_a);
  }

  /**
   * @brief Converts to 8-bit sRGB CRGB, discarding alpha.
   * @return The pixel downcast to CRGB.
   * @details Explicit so a Color4 never silently round-trips through 8-bit gamma.
   */
  explicit operator CRGB() const { return static_cast<CRGB>(color); }
};

/**
 * @brief The interpolatable pixel of a lookup-table entry.
 */
__attribute__((always_inline)) inline Pixel lut_entry_pixel(const Color4 &e) {
  return e.color;
}

/**
 * @brief Lower entry of a fractional lookup-table index.
 * @param idx Fractional index; must be non-negative and non-NaN.
 * @return The truncated index; callers pin `>= size - 1` to the last entry.
 */
__attribute__((always_inline)) inline int lut_index_lo(float idx) {
  return static_cast<int>(idx);
}

/**
 * @brief The lerp16 weight from a lookup-table index toward entry `lo + 1`.
 * @param idx Fractional index; must be non-negative and non-NaN.
 * @param lo Its lower entry, from lut_index_lo, with `lo + 1` still in range.
 * @return The fractional part quantized to [0, 65535].
 * @details One spelling of this arithmetic for every sampler: -ffast-math may
 * compile two spellings of the same expression differently. The fractional part
 * of a non-negative index is in [0, 1) by construction, so quantizing it needs
 * no clamp.
 */
__attribute__((always_inline)) inline uint16_t lut_index_weight(float idx,
                                                                int lo) {
  const float frac = idx - static_cast<float>(lo);
  return static_cast<uint16_t>(frac * 65535.0f + 0.5f);
}

/**
 * @brief Samples a color lookup table at a fractional index, interpolating
 * between adjacent entries.
 * @tparam Entry Table element type accepted by lut_entry_pixel.
 * @param table Table of at least @p size entries.
 * @param size Entry count.
 * @param idx Fractional index; must be non-negative and non-NaN, and is pinned
 * to the last entry from above.
 * @return The interpolated pixel.
 */
template <typename Entry>
__attribute__((always_inline)) inline Pixel
lut_sample_pixel(const Entry *table, int size, float idx) {
  const int lo = lut_index_lo(idx);
  if (lo >= size - 1)
    return lut_entry_pixel(table[size - 1]);
  return lut_entry_pixel(table[lo]).lerp16(lut_entry_pixel(table[lo + 1]),
                                           lut_index_weight(idx, lo));
}

#include "color/color_luts.h"

inline uint16_t srgb_to_linear(uint8_t srgb) {
  return srgb_to_linear_lut[srgb];
}

/**
 * @brief sRGB float [0,1] -> 16-bit linear, interpolating the 256-entry LUT.
 * @param s_srgb sRGB value; out-of-range or NaN inputs are clamped to [0, 1]
 * internally (required for float->int cast safety).
 * @return 16-bit linear channel value.
 * @details Lerps between the two bracketing LUT entries by the fractional part
 * of s*255 (no powf). Lerping the convex sRGB transfer in linear space adds a
 * small upward (secant) bias versus exact powf.
 */
inline uint16_t srgb_to_linear_interp(float s_srgb) {
  // Clamp before the int cast: NaN/out-of-range would be float->int UB below.
  s_srgb = hs::clamp(s_srgb, 0.0f, 1.0f);
  float f = s_srgb * 255.0f;
  int i = static_cast<int>(f);
  if (i >= 255)
    return srgb_to_linear_lut[255];
  float frac = f - static_cast<float>(i);
  float lo = static_cast<float>(srgb_to_linear_lut[i]);
  float hi = static_cast<float>(srgb_to_linear_lut[i + 1]);
  return static_cast<uint16_t>(lo + (hi - lo) * frac + 0.5f);
}

/**
 * @brief Lossy 16-bit-linear -> 8-bit-sRGB downcast.
 * @return The color quantized to an 8-bit sRGB CRGB.
 */
inline Pixel::operator CRGB() const {
  return CRGB(linear_to_srgb8(r), linear_to_srgb8(g), linear_to_srgb8(b));
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Blending Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Packs g|b into the low add lane and r into a separate lane (its high halfword
// stays 0, so uqadd16's upper add is a harmless 0+0).
inline Pixel pixel_blend_add_packed(const Pixel &c1, const Pixel &c2) {
  uint32_t bg1 = ((uint32_t)c1.g << 16) | c1.b;
  uint32_t bg2 = ((uint32_t)c2.g << 16) | c2.b;
  uint32_t sum_bg = inline_uqadd16(bg1, bg2);

  uint32_t sum_r = inline_uqadd16((uint32_t)c1.r, (uint32_t)c2.r);

  return Pixel((uint16_t)sum_r, (uint16_t)(sum_bg >> 16), (uint16_t)sum_bg);
}

/**
 * @brief Returns a straight-alpha "over" functor.
 * @param a Source coverage in [0, 1]; NaN maps to the hi bound.
 * @return A functor taking (dst, src) Pixels and returning src * a +
 *         dst * (1 - a).
 */
inline auto blend_alpha(float a) {
  uint16_t ai = frac_to_q16(a);
  return
      [ai](const Pixel &dst, const Pixel &src) { return dst.lerp16(src, ai); };
}

/**
 * @brief A constexpr-compatible RGB pixel structure for Flash storage.
 * Layout compatible with CRGB but without non-constexpr constructors.
 */
struct CPixel {
  uint8_t r, g, b;
  /**
   * @brief Constructs a black CPixel (all channels zero).
   */
  constexpr CPixel() : r(0), g(0), b(0) {}
  /**
   * @brief Constructs a CPixel from explicit 8-bit channels.
   * @param r Red channel in [0, 255].
   * @param g Green channel in [0, 255].
   * @param b Blue channel in [0, 255].
   */
  constexpr CPixel(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
  /**
   * @brief Constructs a CPixel from a packed 0xRRGGBB hex value.
   * @param hex Packed color; bits 16-23 red, 8-15 green, 0-7 blue.
   * @details `explicit` to match the file's explicit-cast policy: a packed hex
   * is a deliberate construction (`CPixel{0xRRGGBB}`), so a stray int can't
   * silently decay into a color through an implicit conversion.
   */
  constexpr explicit CPixel(uint32_t hex)
      : r((hex >> 16) & 0xFF), g((hex >> 8) & 0xFF), b(hex & 0xFF) {}
  /**
   * @brief Constructs a CPixel from a FastLED CRGB.
   * @param c Source color in 8-bit sRGB space.
   */
  CPixel(const CRGB &c) : r(c.r), g(c.g), b(c.b) {}

  /**
   * @brief Converts to a 16-bit linear Pixel via CRGB.
   * @return The color promoted to linear-space Pixel.
   */
  operator Pixel() const { return CRGB(r, g, b); }
};
