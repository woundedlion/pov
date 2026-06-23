/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "platform.h"
#include "3dmath.h"
#include "util.h"

#include "memory.h"

#if defined(__ARM_FEATURE_DSP)
// Inline assembly avoids a CMSIS header dependency for the saturating add.
__attribute__((always_inline)) static inline uint32_t inline_uqadd16(uint32_t a, uint32_t b) {
  uint32_t res;
  __asm__ volatile("uqadd16 %0, %1, %2" : "=r"(res) : "r"(a), "r"(b));
  return res;
}
#else
// Portable software model of ARM `uqadd16`: two independent 16-bit unsigned
// saturating adds, one per halfword lane. The host has no uqadd16 instruction,
// so this lets the device's exact lane-packing (pixel16_blend_add_packed)
// compile and be unit-tested natively, pinning the bit layout the asm path
// can't run to verify. Not on any host hot path: the host add operators use
// the plain per-channel std::min form below; this exists for the device-parity
// test only.
inline uint32_t inline_uqadd16(uint32_t a, uint32_t b) {
  uint32_t lo = (a & 0xFFFFu) + (b & 0xFFFFu);
  uint32_t hi = (a >> 16) + (b >> 16);
  if (lo > 0xFFFFu) lo = 0xFFFFu;
  if (hi > 0xFFFFu) hi = 0xFFFFu;
  return (hi << 16) | lo;
}
#endif

struct Pixel16;
// The device's saturating per-channel add, packed into two uqadd16 lanes
// (g|b in one 32-bit word, r alone in another). Shared by Pixel16::operator+=
// and blend_add on the device so the lane layout has a single definition, and
// compiled on the host (via the software uqadd16 above) so a unit test can pin
// that layout against an independent per-channel reference.
inline Pixel16 pixel16_blend_add_packed(const Pixel16 &c1, const Pixel16 &c2);

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
struct Pixel16 {
  uint16_t r, g, b;

  /**
   * @brief Constructs a black pixel (all channels zero).
   */
  constexpr Pixel16() : r(0), g(0), b(0) {}
  /**
   * @brief Constructs a pixel from explicit 16-bit linear channels.
   * @param _r Red channel in [0, 65535].
   * @param _g Green channel in [0, 65535].
   * @param _b Blue channel in [0, 65535].
   */
  constexpr Pixel16(uint16_t _r, uint16_t _g, uint16_t _b)
      : r(_r), g(_g), b(_b) {}

  /**
   * @brief Constructs a pixel from HSV (converts to sRGB then Linear).
   * @param hsv Source color in HSV space.
   */
  Pixel16(const CHSV &hsv) {
    CRGB srgb(hsv);
    r = srgb_to_linear(srgb.r);
    g = srgb_to_linear(srgb.g);
    b = srgb_to_linear(srgb.b);
  }

  /**
   * @brief Constructs a pixel from CRGB (converts to Linear).
   * @param c Source color in 8-bit sRGB space.
   */
  Pixel16(const CRGB &c) {
    r = srgb_to_linear(c.r);
    g = srgb_to_linear(c.g);
    b = srgb_to_linear(c.b);
  }

  /**
   * @brief Lossy 16-bit-linear -> 8-bit-sRGB downcast.
   * @return The color quantized to an 8-bit sRGB CRGB.
   * @details Forward-declared; defined once CRGB is complete. Explicit so a
   * stray Pixel16 in a CRGB context is a compile error rather than a silent
   * round-trip through 8-bit gamma — sink call sites that genuinely need the
   * cast spell it out.
   */
  explicit operator CRGB() const;

  /**
   * @brief Saturated per-channel addition into this pixel.
   * @param rhs Pixel to add.
   * @return Reference to this pixel after the clamped add.
   */
  Pixel16 &operator+=(const Pixel16 &rhs) {
#if defined(__ARM_FEATURE_DSP)
    *this = pixel16_blend_add_packed(*this, rhs);
#else
    r = (uint16_t)std::min((uint32_t)65535, (uint32_t)r + rhs.r);
    g = (uint16_t)std::min((uint32_t)65535, (uint32_t)g + rhs.g);
    b = (uint16_t)std::min((uint32_t)65535, (uint32_t)b + rhs.b);
#endif
    return *this;
  }

  /**
   * @brief Scales every channel by a float factor (saturated).
   * @param s Scale factor; may be any finite float (NaN maps to the hi bound).
   * @return A new pixel with each channel clamped to [0, 65535].
   * @details Rounds to nearest (+0.5f) rather than truncating, matching the
   * +0.5f / +32768 round-to-nearest parity used by the rest of this file (this
   * is the per-frame fade multiply, so truncation would systematically bias
   * every faded channel down each frame). Clamps in float before the cast
   * because r*s can exceed INT_MAX and float->int conversion is UB out of
   * range; hs::clamp also maps a NaN scale to the hi bound instead of letting
   * NaN reach the cast. The +0.5f sits inside the clamp so the saturating hi
   * bound stays exactly 65535.
   */
  Pixel16 operator*(float s) const {
    // Round (+0.5f) and clamp in float before the cast: r*s can exceed INT_MAX
    // for large s, and float->int conversion is UB out of range, so an int
    // clamp would fire too late. hs::clamp (not std::clamp) also maps a NaN
    // scale to the hi bound instead of letting NaN reach the cast.
    return Pixel16((uint16_t)hs::clamp(r * s + 0.5f, 0.0f, 65535.0f),
                   (uint16_t)hs::clamp(g * s + 0.5f, 0.0f, 65535.0f),
                   (uint16_t)hs::clamp(b * s + 0.5f, 0.0f, 65535.0f));
  }

  /**
   * @brief Linearly interpolates 16-bit between this pixel and another.
   * @param other Target pixel at frac == 65535.
   * @param frac Blend weight in [0, 65535]; 0 yields this pixel, 65535 yields other.
   * @return The interpolated pixel, round-to-nearest per channel.
   * @details Round-to-nearest div-by-65535 via shifts:
   *   round(x/65535) == (x + (x>>16) + 32768) >> 16. The (x + (x>>16)) term
   * reconstructs x*65536/65535 to within <1 LSB; adding half the divisor
   * (32768) before the >>16 rounds instead of truncating, matching the +0.5f
   * round-to-nearest parity used by the rest of this file. Still exact at the
   * endpoints (frac 0/65535 -> a/b) because the bias stays strictly below one
   * output quantum. The two channel products use plain 32-bit MACs on every
   * target rather than a packed `smlad`: `smlad` is a *signed* 16x16 dual-MAC,
   * so any operand >= 32768 reads as negative and the product is wrong across
   * the upper half of the range. On a Cortex-M7 a `mul` + `mla` per channel is
   * single-cycle and dual-issuable, so the portable form costs at most a cycle
   * or two. (The unsigned `uqadd16` in operator+= is fine — only signed
   * multiply is unsound for unsigned operands.)
   */
  Pixel16 lerp16(const Pixel16 &other, uint16_t frac) const {
    uint16_t inv = 65535 - frac;
    uint32_t xr = (uint32_t)r * inv + (uint32_t)other.r * frac;
    uint32_t xg = (uint32_t)g * inv + (uint32_t)other.g * frac;
    uint32_t xb = (uint32_t)b * inv + (uint32_t)other.b * frac;
    uint32_t r32 = (xr + (xr >> 16) + 32768) >> 16;
    uint32_t g32 = (xg + (xg >> 16) + 32768) >> 16;
    uint32_t b32 = (xb + (xb >> 16) + 32768) >> 16;
    return Pixel16((uint16_t)r32, (uint16_t)g32, (uint16_t)b32);
  }

  /**
   * @brief Tests two pixels for exact channel equality.
   * @param rhs Pixel to compare against.
   * @return True if all three channels match.
   */
  bool operator==(const Pixel16 &rhs) const {
    return r == rhs.r && g == rhs.g && b == rhs.b;
  }

  /**
   * @brief Tests two pixels for channel inequality.
   * @param rhs Pixel to compare against.
   * @return True if any channel differs.
   */
  bool operator!=(const Pixel16 &rhs) const { return !(*this == rhs); }

  /**
   * @brief Tests equality against an HSV color (converted to Pixel16).
   * @param rhs Color in HSV space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CHSV &rhs) const { return *this == Pixel16(rhs); }

  /**
   * @brief Tests inequality against an HSV color (converted to Pixel16).
   * @param rhs Color in HSV space.
   * @return True if this pixel differs from the converted color.
   */
  bool operator!=(const CHSV &rhs) const { return !(*this == rhs); }

  /**
   * @brief Tests equality against a CRGB color (converted to Pixel16).
   * @param rhs Color in 8-bit sRGB space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CRGB &rhs) const { return *this == Pixel16(rhs); }

  /**
   * @brief Tests inequality against a CRGB color (converted to Pixel16).
   * @param rhs Color in 8-bit sRGB space.
   * @return True if this pixel differs from the converted color.
   */
  bool operator!=(const CRGB &rhs) const { return !(*this == rhs); }
};

using Pixel = Pixel16;

/**
 * @brief Represents a color with a STRAIGHT (non-premultiplied) alpha channel.
 * @details Alpha is straight: `color` holds the un-premultiplied color and
 * `alpha` its coverage/opacity. `lerp` interpolates the two independently,
 * consistent with that convention. The arithmetic operators (`operator*=`,
 * `operator+=`) are NOT alpha compositing — they treat (color, alpha) as a plain
 * 4-vector for the SSAA "scale each sample by 1/N, then sum" averaging algebra,
 * so `*=` deliberately scales alpha alongside color. Premultiplication happens
 * once, at the final canvas write (`color * alpha`); do not read these operators
 * as over/under compositing.
 */
struct Color4 {
  Pixel color;
  float alpha;

  /**
   * @brief Constructs an opaque black color (alpha 1.0).
   */
  Color4() : color(Pixel(0, 0, 0)), alpha(1.0f) {}
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
   * @details `explicit` so the sRGB→linear convention is opt-in: a caller that
   *          models Color4 as already-linear must spell out `Color4(r,g,b)`
   *          rather than have a braced `{r,g,b}` silently take the sRGB path.
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
   * @details t is clamped to [0,1] so an out-of-range t saturates both channels
   * at an endpoint rather than letting alpha extrapolate while color stays
   * clamped.
   */
  Color4 lerp(const Color4 &other, float t) const {
    const float ct = hs::clamp(t, 0.0f, 1.0f);
    uint16_t frac = static_cast<uint16_t>(ct * 65535.0f + 0.5f);
    Pixel blended = color.lerp16(other.color, frac);
    float blended_a = alpha + (other.alpha - alpha) * ct;
    return Color4(blended, blended_a);
  }

  /**
   * @brief Adds another color's pixel and alpha into this one (both saturating).
   * @param rhs Color to add.
   * @return Reference to this color after the add.
   * @details The pixel add saturates at the channel max; alpha saturates at 1.0
   * to match, so the sum's alpha stays a valid blend weight. The SSAA averaging
   * path (scale each sample by 1/SAMPLES, then sum) keeps the running alpha <= 1
   * regardless, so the clamp only ever trims floating-point spillover or a
   * caller that sums un-normalized weights.
   */
  Color4 &operator+=(const Color4 &rhs) {
    color += rhs.color;
    alpha = hs::clamp(alpha + rhs.alpha, 0.0f, 1.0f);
    return *this;
  }

  /**
   * @brief Scales both pixel and alpha by a float factor.
   * @param s Scale factor.
   * @return Reference to this color after scaling.
   */
  Color4 &operator*=(float s) {
    color = color * s;
    alpha *= s;
    return *this;
  }

  /**
   * @brief Converts to 8-bit sRGB CRGB, discarding alpha.
   * @return The pixel downcast to CRGB.
   * @details Explicit (like Pixel16::operator CRGB()) so a Color4 never silently
   * round-trips through 8-bit gamma; callers must opt in with an explicit cast.
   */
  explicit operator CRGB() const { return static_cast<CRGB>(color); }
};

/**
 * @brief Perceptual (OKLab) hue rotation with a precomputed rotation.
 * @param c Source color.
 * @param ca Cosine of the rotation angle.
 * @param sa Sine of the rotation angle.
 * @return The hue-rotated color.
 * @details Defined after the OKLab conversion helpers below. The (ca, sa)
 * overload takes a precomputed rotation so frame-constant callers can hoist the
 * sin/cos out of the per-pixel loop.
 */
inline Color4 hue_rotate(const Color4 &c, float ca, float sa);
/**
 * @brief Perceptual (OKLab) hue rotation by a turn amount.
 * @param c Source color.
 * @param amount Rotation in turns (0..1 = full turn).
 * @return The hue-rotated color.
 */
inline Color4 hue_rotate(const Color4 &c, float amount);

#include "color_luts.h"

// Definition; canonical docs are on the declaration above.
inline uint16_t srgb_to_linear(uint8_t srgb) {
  return srgb_to_linear_lut[srgb];
}

/**
 * @brief sRGB float [0,1] -> 16-bit linear, interpolating the 256-entry LUT.
 * @param s_srgb sRGB value in [0, 1]; callers must clamp before calling.
 * @return 16-bit linear channel value.
 * @details Lerps between the two bracketing LUT entries by the fractional part
 * of s*255, so the sRGB->linear step keeps sub-8-bit precision at the cost of a
 * LUT lookup + lerp (no powf).
 */
inline uint16_t srgb_to_linear_interp(float s_srgb) {
  float f = s_srgb * 255.0f; // [0, 255]
  // The contract is [0,1], but a negative input would make (int)f negative — an
  // OOB LUT read — and a negative frac that wraps the uint16_t cast. Clamp the
  // low end defensively, mirroring the i>=255 high-end guard.
  if (f <= 0.0f)
    return srgb_to_linear_lut[0];
  int i = static_cast<int>(f);
  if (i >= 255)
    return srgb_to_linear_lut[255];
  float frac = f - static_cast<float>(i);
  float lo = static_cast<float>(srgb_to_linear_lut[i]);
  float hi = static_cast<float>(srgb_to_linear_lut[i + 1]);
  return static_cast<uint16_t>(lo + (hi - lo) * frac + 0.5f);
}

/**
 * @brief Lossy 16-bit-linear -> 8-bit-sRGB downcast via the LUT.
 * @return The color quantized to an 8-bit sRGB CRGB.
 */
inline Pixel16::operator CRGB() const {
  return CRGB(linear_to_srgb_lut[r], linear_to_srgb_lut[g],
              linear_to_srgb_lut[b]);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Blending Functions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * @brief Blends two pixels by taking the maximum of each component.
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return The max blend result.
 */
inline Pixel blend_max(const Pixel &c1, const Pixel &c2) {
  return Pixel(std::max(c1.r, c2.r), std::max(c1.g, c2.g),
               std::max(c1.b, c2.b));
}

/**
 * @brief Blends two pixels by overwriting (Painter's algorithm).
 * @param c2 Source pixel.
 * @return c2 (Source).
 */
inline Pixel blend_over(const Pixel &, const Pixel &c2) { return c2; }

/**
 * @brief Blends two pixels by keeping the destination (Under).
 * @param c1 Destination pixel.
 * @return c1 (Destination).
 */
inline Pixel blend_under(const Pixel &c1, const Pixel &) { return c1; }

/**
 * @brief Blends two pixels by additive mixing (saturated).
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return c1 + c2 (clamped).
 */
inline Pixel blend_add(const Pixel &c1, const Pixel &c2) {
  // Saturated Add
#if defined(__ARM_FEATURE_DSP)
  return pixel16_blend_add_packed(c1, c2);
#else
  uint32_t r = (uint32_t)c1.r + c2.r;
  uint32_t g = (uint32_t)c1.g + c2.g;
  uint32_t b = (uint32_t)c1.b + c2.b;
  return Pixel((r > 65535) ? 65535 : (uint16_t)r,
               (g > 65535) ? 65535 : (uint16_t)g,
               (b > 65535) ? 65535 : (uint16_t)b);
#endif
}

// Definition of the device packed-add helper declared above Pixel16. Packs g|b
// into the low add lane and r into a separate lane (its high halfword stays 0,
// so uqadd16's upper add is a harmless 0+0). Identical layout to the asm the
// device runs; on the host the software inline_uqadd16 drives it so the test
// suite can pin the lane assignment.
inline Pixel16 pixel16_blend_add_packed(const Pixel16 &c1, const Pixel16 &c2) {
  uint32_t bg1 = ((uint32_t)c1.g << 16) | c1.b;
  uint32_t bg2 = ((uint32_t)c2.g << 16) | c2.b;
  uint32_t sum_bg = inline_uqadd16(bg1, bg2);

  uint32_t sum_r = inline_uqadd16((uint32_t)c1.r, (uint32_t)c2.r);

  return Pixel16((uint16_t)sum_r, (uint16_t)(sum_bg >> 16), (uint16_t)sum_bg);
}

/**
 * @brief Returns a blend functor that lerps c1->c2 by alpha.
 * @param a Blend weight in [0, 1]; NaN maps to the hi bound.
 * @return A functor taking (c1, c2) Pixels and returning the lerped Pixel.
 * @details Rounds to nearest (+0.5f) and clamps in float before the cast,
 * mirroring Pixel16::operator* and the +0.5f weight conversion in Color4::lerp
 * / Gradient::get: casting an unclamped a*65535 to int is UB when a is NaN or
 * overflows int, and rounding keeps the [0,1]->[0,65535] map consistent with
 * the rest of the file. The +0.5f sits inside the clamp so a == 1 still maps
 * exactly to 65535 (no saturation past the hi bound).
 */
inline auto blend_alpha(float a) {
  // hs::clamp is a hardware min/max that also maps NaN to the hi bound.
  uint16_t ai = (uint16_t)hs::clamp(a * 65535.0f + 0.5f, 0.0f, 65535.0f);
  return [ai](const Pixel &c1, const Pixel &c2) { return c1.lerp16(c2, ai); };
}

/**
 * @brief Blends two pixels by averaging/mean.
 * @param c1 Destination pixel.
 * @param c2 Source pixel.
 * @return (c1 + c2) / 2, rounded to nearest (matching this file's +0.5f /
 *   +32768 round-to-nearest parity rather than truncating downward).
 */
inline Pixel blend_mean(const Pixel &c1, const Pixel &c2) {
  return Pixel((c1.r + c2.r + 1) / 2, (c1.g + c2.g + 1) / 2,
               (c1.b + c2.b + 1) / 2);
}



///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Defines types of color harmonies for generative palettes.
 */
enum class HarmonyType {
  TRIADIC,
  SPLIT_COMPLEMENTARY,
  COMPLEMENTARY,
  ANALOGOUS
};

/**
 * @brief Defines the visual shape or distribution of colors across the palette
 * domain.
 */
enum class GradientShape { STRAIGHT, CIRCULAR, VIGNETTE, FALLOFF };

/**
 * @brief Defines the overall brightness profile across the palette domain.
 */
enum class BrightnessProfile { ASCENDING, DESCENDING, FLAT, BELL, CUP };

/**
 * @brief Defines the saturation profile.
 */
enum class SaturationProfile { PASTEL, MID, VIBRANT };

///////////////////////////////////////////////////////////////////////////////

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

/**
 * @brief High-precision sRGB float [0,1] -> linear float [0,1].
 * @param s sRGB value in [0, 1].
 * @return Linear value in [0, 1].
 * @details Not constexpr: powf is not a constant expression, so a constexpr
 * marking here could never be evaluated at compile time (ill-formed NDR).
 * inline for ODR.
 */
inline float srgb_to_linear_float(float s) {
  return (s <= 0.04045f) ? s / 12.92f : powf((s + 0.055f) / 1.055f, 2.4f);
}

/**
 * @brief Inverse: linear float [0,1] -> sRGB float [0,1].
 * @param l Linear value in [0, 1].
 * @return sRGB value in [0, 1].
 */
inline float linear_to_srgb_float(float l) {
  return (l <= 0.0031308f) ? l * 12.92f : 1.055f * powf(l, 1.0f / 2.4f) - 0.055f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OKLab / OKLCH Color Space (Björn Ottosson, 2020)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * @brief OKLab perceptual color: lightness L and chroma axes a, b.
 */
struct OKLab { float L, a, b; };
/**
 * @brief OKLCH polar color: lightness L, chroma C, hue h (radians).
 */
struct OKLCH { float L, C, h; };

/** @brief Cone-response (LMS) triple, the OKLab intermediate before the cube-root
 *  nonlinearity. */
struct LMS { float l, m, s; };

/**
 * @brief Linear-RGB -> LMS cone response (the first OKLab matrix).
 * @param r Linear red in [0, 1].
 * @param g Linear green in [0, 1].
 * @param b Linear blue in [0, 1].
 * @return The (l, m, s) cone responses, before the cube-root nonlinearity.
 * @details Shared by linear_rgb_to_oklab and hue_rotate so the matrix lives in
 * one place; both then apply their own cube-root (exact cbrtf vs fast_cbrt) and
 * call lms_to_oklab. Inline — zero runtime cost.
 */
inline LMS linear_rgb_to_lms(float r, float g, float b) {
  return {0.4122214708f * r + 0.5363325363f * g + 0.0514459929f * b,
          0.2119034982f * r + 0.6806995451f * g + 0.1073969566f * b,
          0.0883024619f * r + 0.2817188376f * g + 0.6299787005f * b};
}

/**
 * @brief Cube-rooted LMS -> OKLab (the second OKLab matrix).
 * @param l_ Cube-rooted l cone response.
 * @param m_ Cube-rooted m cone response.
 * @param s_ Cube-rooted s cone response.
 * @return The color in OKLab space.
 * @details Takes the already-nonlinearized (cube-rooted) triple so the caller
 * picks the cube-root flavour. Shared with hue_rotate; inline, zero runtime cost.
 */
inline OKLab lms_to_oklab(float l_, float m_, float s_) {
  return {0.2104542553f * l_ + 0.7936177850f * m_ - 0.0040720468f * s_,
          1.9779984951f * l_ - 2.4285922050f * m_ + 0.4505937099f * s_,
          0.0259040371f * l_ + 0.7827717662f * m_ - 0.8086757660f * s_};
}

/**
 * @brief Converts linear RGB [0,1] to OKLab.
 * @param r Linear red in [0, 1].
 * @param g Linear green in [0, 1].
 * @param b Linear blue in [0, 1].
 * @return The color in OKLab space.
 */
inline OKLab linear_rgb_to_oklab(float r, float g, float b) {
  LMS lms = linear_rgb_to_lms(r, g, b);
  return lms_to_oklab(cbrtf(lms.l), cbrtf(lms.m), cbrtf(lms.s));
}

/**
 * @brief Converts OKLab to linear RGB [0,1].
 * @param lab Source color in OKLab space.
 * @param r Out: linear red in [0, 1] (may exit gamut before clamping).
 * @param g Out: linear green in [0, 1].
 * @param b Out: linear blue in [0, 1].
 */
inline void oklab_to_linear_rgb(OKLab lab, float &r, float &g, float &b) {
  float l_ = lab.L + 0.3963377774f * lab.a + 0.2158037573f * lab.b;
  float m_ = lab.L - 0.1055613458f * lab.a - 0.0638541728f * lab.b;
  float s_ = lab.L - 0.0894841775f * lab.a - 1.2914855480f * lab.b;

  float l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_;

  r = +4.0767416621f * l - 3.3077115913f * m + 0.2309699292f * s;
  g = -1.2684380046f * l + 2.6097574011f * m - 0.3413193965f * s;
  b = -0.0041960863f * l - 0.7034186147f * m + 1.7076147010f * s;
}

/**
 * @brief Tests whether a linear-RGB triple lies inside the [0,1] display cube.
 * @param r Linear red.
 * @param g Linear green.
 * @param b Linear blue.
 * @return True if every channel is within [0,1] (with a small epsilon slack).
 * @details The slack absorbs the float rounding that can leave an in-gamut color
 * a hair past 1.0 after the OKLab inverse, so a color in gamut up to ULP noise
 * takes the cheap clamp path rather than spending a chroma search.
 */
inline bool linear_rgb_in_gamut(float r, float g, float b) {
  constexpr float lo = -1e-4f, hi = 1.0f + 1e-4f;
  return r >= lo && r <= hi && g >= lo && g <= hi && b >= lo && b <= hi;
}

/**
 * @brief Maps an out-of-gamut OKLab color into the sRGB cube by reducing chroma,
 * holding hue and lightness fixed (Ottosson's preserve-chroma projection).
 * @param lab Source color, assumed out of gamut (callers gate on
 * linear_rgb_in_gamut) and with L in [0,1].
 * @return An OKLab color with the same L and hue but chroma scaled down until it
 * is just inside the display cube.
 * @details Binary-searches a uniform scale s in [0,1] on the (a,b) chroma axes.
 * Scaling (a,b) uniformly leaves hue = atan2(b,a) and L untouched and shrinks
 * chroma C = hypot(a,b) by s, so this is a pure chroma reduction. s = 0 is the
 * achromatic color at L — always in gamut for L in [0,1], since the gray axis
 * maps to r=g=b=L^3 — which gives the search a guaranteed in-gamut floor; s = 1
 * is the out-of-gamut input. 16 bisections pin s to < 2^-16, below one 16-bit
 * output quantum. This is the chroma-preserving alternative to a per-channel RGB
 * clip, which twists hue on saturated colors — the drift that accumulates in the
 * hue-rotate feedback loop. Cold path only: the per-pixel hot path reaches it
 * solely for the rare pixel that leaves gamut.
 */
inline OKLab gamut_clip_preserve_chroma(OKLab lab) {
  float lo = 0.0f, hi = 1.0f;
  for (int i = 0; i < 16; ++i) {
    float mid = 0.5f * (lo + hi);
    float r, g, b;
    oklab_to_linear_rgb({lab.L, lab.a * mid, lab.b * mid}, r, g, b);
    if (linear_rgb_in_gamut(r, g, b))
      lo = mid; // still inside: push the scale up
    else
      hi = mid; // outside: pull it back
  }
  // Return the largest known-in-gamut scale (lo), never hi.
  return {lab.L, lab.a * lo, lab.b * lo};
}

/**
 * @brief OKLab -> linear RGB with chroma-reduction gamut mapping off the fast
 * path.
 * @param lab Source color in OKLab space.
 * @param r Out: linear red in [0,1] (may sit a hair past the bound; callers
 * still clamp).
 * @param g Out: linear green.
 * @param b Out: linear blue.
 * @details Converts directly first; only when the result leaves the [0,1] cube
 * does it pay for gamut_clip_preserve_chroma and re-convert. In-gamut colors (the
 * overwhelming majority) cost one matrix mul plus the gate test — no search — so
 * this is a drop-in for oklab_to_linear_rgb on paths that want hue held under
 * clipping without taxing the common case.
 */
inline void oklab_to_linear_rgb_gamut(OKLab lab, float &r, float &g, float &b) {
  oklab_to_linear_rgb(lab, r, g, b);
  if (!linear_rgb_in_gamut(r, g, b))
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r, g, b);
}

/**
 * @brief Quantizes a [0,1] linear channel to a 16-bit Pixel component.
 * @param v Linear channel value; clamped to [0, 1].
 * @return The channel as a 16-bit value in [0, 65535].
 * @details Clamps, then rounds (+0.5f) rather than truncating — truncation
 * would bias every channel down by up to ~1/65535. Shared by the float->Pixel
 * helpers below (oklch_to_pixel, srgb_to_linear_interp).
 */
inline uint16_t float_to_pixel16(float v) {
  return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * 65535.0f + 0.5f);
}

// Definition; canonical docs are on the declaration above. Rotates the (a,b)
// chroma plane in the OKLab perceptual color space, which preserves perceived
// lightness (L) and chroma (sqrt(a^2+b^2)) across the shift. A linear-RGB
// rotation cannot: its equal angular steps are not perceptually uniform and its
// perceived brightness drifts with hue. Hot path: this is the default feedback
// color transform (Feedback::hue_fade, run per pixel every frame) and Flyby's
// per-pixel shader. To stay affordable, the forward linear->OKLab nonlinearity
// uses fast_cbrt (~2e-5 rel error) while the inverse is exact cubes, and the
// shift is a direct 2D rotation of (a,b) — no atan2/sqrt of a full OKLCH polar
// round-trip.
inline Color4 hue_rotate(const Color4 &c, float ca, float sa) {
  constexpr float INV16 = 1.0f / 65535.0f;
  float r = c.color.r * INV16, g = c.color.g * INV16, b = c.color.b * INV16;

  // linear RGB -> OKLab, sharing the two matrices with linear_rgb_to_oklab but
  // substituting fast_cbrt for the forward nonlinearity (the hot-path cost cut).
  LMS lms = linear_rgb_to_lms(r, g, b);
  OKLab lab = lms_to_oklab(fast_cbrt(lms.l), fast_cbrt(lms.m), fast_cbrt(lms.s));
  float L = lab.L, A = lab.a, B = lab.b;

  // Rotate the chroma plane by the precomputed (ca, sa) = (cos, sin) of the
  // turn angle — preserves L and |(A,B)|.
  float A2 = A * ca - B * sa;
  float B2 = A * sa + B * ca;

  // OKLab -> linear RGB (exact inverse) -> 16-bit linear. Off the hot path
  // (only the rare pixel that leaves gamut), reduce chroma to hold hue and L
  // instead of letting the per-channel clamp below twist the hue — this is what
  // keeps the hue-rotate feedback loop from drifting saturated colors frame over
  // frame. In-gamut pixels (the common case) skip the search and fall straight
  // through to the clamp.
  float nr, ng, nb;
  oklab_to_linear_rgb_gamut({L, A2, B2}, nr, ng, nb);

  Color4 result = c;
  // Round, don't truncate: float_to_pixel16 applies the +0.5f rule (and the
  // [0,1] gamut clamp) so this per-pixel hot path doesn't bias every channel
  // down by up to ~1/65535 every frame.
  result.color.r = float_to_pixel16(nr);
  result.color.g = float_to_pixel16(ng);
  result.color.b = float_to_pixel16(nb);
  return result;
}

// Definition; canonical docs are on the declaration above. Computes the
// rotation per call. On hot paths where `amount` is constant across the frame
// (the default Feedback::hue_fade) call the (ca, sa) overload with values
// hoisted out of the per-pixel loop instead (Style::sync_hue populates
// Style::hue_ca / hue_sa). Callers with a genuinely per-fragment amount (Flyby)
// use this overload.
inline Color4 hue_rotate(const Color4 &c, float amount) {
  float angle = amount * (2.0f * PI_F);
  return hue_rotate(c, fast_cosf(angle), fast_sinf(angle));
}

/**
 * @brief Converts OKLab (Cartesian a,b) to OKLCH (polar C, h).
 * @param lab Source color in OKLab space.
 * @return The color in OKLCH space; hue h in radians.
 */
inline OKLCH oklab_to_oklch(OKLab lab) {
  float C = sqrtf(lab.a * lab.a + lab.b * lab.b);
  float h = atan2f(lab.b, lab.a);
  return {lab.L, C, h};
}

/**
 * @brief Converts OKLCH (polar) to OKLab (Cartesian).
 * @param lch Source color in OKLCH space.
 * @return The color in OKLab space.
 */
inline OKLab oklch_to_oklab(OKLCH lch) {
  return {lch.L, lch.C * cosf(lch.h), lch.C * sinf(lch.h)};
}

/**
 * @brief Convenience: sRGB [0-255] channels to OKLCH.
 * @param r Red channel in [0, 255].
 * @param g Green channel in [0, 255].
 * @param b Blue channel in [0, 255].
 * @return The color in OKLCH space.
 */
inline OKLCH srgb_to_oklch(uint8_t r, uint8_t g, uint8_t b) {
  float rf = srgb_to_linear_float(r / 255.0f);
  float gf = srgb_to_linear_float(g / 255.0f);
  float bf = srgb_to_linear_float(b / 255.0f);
  return oklab_to_oklch(linear_rgb_to_oklab(rf, gf, bf));
}

/**
 * @brief Converts OKLCH to a 16-bit linear Pixel (gamut-clamped).
 * @param lch Source color in OKLCH space.
 * @return The color as a linear-space Pixel.
 */
inline Pixel oklch_to_pixel(OKLCH lch) {
  float r, g, b;
  // Chroma-reduce rather than per-channel clip when the OKLCH color is past the
  // sRGB cusp (off the fast path: in-gamut stops skip the search), so a vivid
  // palette color keeps its hue and lightness instead of clipping toward a
  // primary.
  oklab_to_linear_rgb_gamut(oklch_to_oklab(lch), r, g, b);
  return Pixel(float_to_pixel16(r), float_to_pixel16(g), float_to_pixel16(b));
}

/**
 * @brief Interpolates two OKLCH colors along the shortest-arc hue.
 * @param a Start color at t == 0.
 * @param b End color at t == 1.
 * @param t Blend weight; may extrapolate outside [0, 1].
 * @return The interpolated OKLCH color with L and C clamped valid.
 * @details An extrapolated t (reachable via the unbounded
 * GenerativePalette::lerp / ColorWipe paths) can overshoot a valid endpoint
 * into an invalid OKLCH: negative L renders near-black, negative C flips the
 * hue 180°. In-range t is unaffected since both endpoints are already valid.
 * Hue is left free to wrap.
 */
inline OKLCH lerp_oklch(OKLCH a, OKLCH b, float t) {
  // Handle achromatic cases (near-zero chroma)
  float h;
  if (a.C < 1e-4f && b.C < 1e-4f) {
    h = 0.0f;
  } else if (a.C < 1e-4f) {
    h = b.h;
  } else if (b.C < 1e-4f) {
    h = a.h;
  } else {
    // Shortest arc interpolation
    float dh = b.h - a.h;
    if (dh > PI_F) dh -= 2.0f * PI_F;
    if (dh < -PI_F) dh += 2.0f * PI_F;
    h = a.h + dh * t;
  }
  // Clamp the magnitude channels (see @details); hue is left free to wrap.
  // C is floored at 0 but intentionally not capped above: for the interpolation
  // contract (t in [0,1]) the result lies between two in-gamut endpoints, so no
  // upper bound is needed; an out-of-[0,1] t could overshoot C and hue-shift
  // after the downstream RGB clamp, but no caller extrapolates.
  float L = hs::clamp(a.L + (b.L - a.L) * t, 0.0f, 1.0f);
  float C = __builtin_fmaxf(0.0f, a.C + (b.C - a.C) * t);
  return {L, C, h};
}

/**
 * @brief Interpolates two sRGB CPixels in OKLCH space to a 16-bit linear Pixel.
 * @param c1 Start color at t == 0.
 * @param c2 End color at t == 1.
 * @param t Blend weight.
 * @return The interpolated color as a linear-space Pixel.
 */
inline Pixel lerp_oklch(const CPixel &c1, const CPixel &c2, float t) {
  OKLCH a = srgb_to_oklch(c1.r, c1.g, c1.b);
  OKLCH b = srgb_to_oklch(c2.r, c2.g, c2.b);
  return oklch_to_pixel(lerp_oklch(a, b, t));
}

/**
 * @brief Converts an OKLCH color to an 8-bit sRGB CPixel, chroma-reducing
 *        out-of-gamut colors to hold hue and lightness (Ottosson clip).
 * @param lch Source color in OKLCH space.
 * @return The gamut-mapped color as an 8-bit sRGB CPixel.
 * @details Shares the gamut-mapped OKLab->linear->sRGB tail with
 *          lerp_oklch_srgb (which now calls this) and the OKLCH palette-key
 *          authoring path: a vivid color past the sRGB cusp gets its chroma
 *          pulled to the boundary, holding hue and lightness, rather than the
 *          per-channel clip that twists hue toward a primary. In-gamut colors
 *          pay only the gate test.
 */
inline CPixel oklch_to_cpixel(OKLCH lch) {
  float r, g, b;
  oklab_to_linear_rgb_gamut(oklch_to_oklab(lch), r, g, b);
  return CPixel(
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(r, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)),
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(g, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)),
    static_cast<uint8_t>(hs::clamp(linear_to_srgb_float(hs::clamp(b, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f, 255.0f)));
}

/**
 * @brief Interpolates two sRGB CPixels in OKLCH space to an sRGB CPixel.
 * @param c1 Start color at t == 0.
 * @param c2 End color at t == 1.
 * @param t Blend weight.
 * @return The interpolated color as an 8-bit sRGB CPixel.
 */
inline CPixel lerp_oklch_srgb(const CPixel &c1, const CPixel &c2, float t) {
  OKLCH a = srgb_to_oklch(c1.r, c1.g, c1.b);
  OKLCH b = srgb_to_oklch(c2.r, c2.g, c2.b);
  return oklch_to_cpixel(lerp_oklch(a, b, t));
}

/**
 * @brief Abstract base for all palettes.
 * @details Provides a uniform interface for color lookup, replacing
 * std::variant + std::visit with a single vtable pointer.
 */
class Palette {
public:
  /**
   * @brief Samples the palette at a coordinate.
   * @param t Lookup coordinate, conventionally in [0, 1].
   * @return The color at t.
   */
  virtual Color4 get(float t) const = 0;
  /**
   * @brief Virtual destructor for polymorphic deletion.
   */
  virtual ~Palette() = default;
};

/**
 * @brief A palette backed by a precomputed 256-entry linear-RGB lookup table,
 * filled at construction by interpolating between the color stops in linear
 * space.
 */
class Gradient : public Palette {
public:
  Pixel entries[256];

  /**
   * @brief Builds the 256-entry LUT by interpolating between color stops.
   * @param points Sorted-ascending (position in [0,1], color) stops.
   * @details Stop positions index entries[256] via static_cast<int>(pos * 255);
   * a pos outside [0,1] is an out-of-bounds write into the table, and segments
   * are only filled when end > start, so a transposed/unsorted pair would
   * silently degenerate. Gradient literals are authored once at construction
   * (cold), so the bounds and ordering are trapped always-on rather than
   * corrupting adjacent memory or rendering wrong under NDEBUG on-device.
   */
  Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {
    Pixel black(0, 0, 0);
    for (int i = 0; i < 256; i++)
      entries[i] = black;

    if (points.size() == 0)
      return;

    // Trap out-of-range or unsorted stops always-on at construction (see
    // @details) rather than corrupting the table or rendering wrong.
    float prevCheck = -1.0f;
    for (const auto &stop : points) {
      HS_CHECK(stop.first >= 0.0f && stop.first <= 1.0f,
               "Gradient stop position out of [0,1]");
      HS_CHECK(stop.first >= prevCheck, "Gradient stops must be sorted ascending");
      prevCheck = stop.first;
    }

    auto it = points.begin();
    float prevPos = it->first;
    CPixel prevColor = it->second;

    // Fill start. Convert via the same srgb_to_linear_float + float_to_pixel16
    // path the interpolated segment body uses (below), not the 8-bit
    // srgb_to_linear() LUT, so a stop's flat-region entry is bit-identical to
    // the value an adjacent interpolated segment computes at that stop — no
    // ~1-LSB precision seam at flat/interpolated boundaries. Cold path.
    int firstStop = static_cast<int>(prevPos * 255);
    Pixel prevLinear(float_to_pixel16(srgb_to_linear_float(prevColor.r / 255.0f)),
                     float_to_pixel16(srgb_to_linear_float(prevColor.g / 255.0f)),
                     float_to_pixel16(srgb_to_linear_float(prevColor.b / 255.0f)));
    for (int i = 0; i <= firstStop; i++)
      entries[i] = prevLinear;

    it++;
    while (it != points.end()) {
      float nextPos = it->first;
      CPixel nextColor = it->second;

      int start = static_cast<int>(prevPos * 255);
      int end = static_cast<int>(nextPos * 255);

      // Two stops that quantize to the same index (equal or near-equal
      // positions) leave end == start and skip this segment: the later stop's
      // color simply overwrites the index in a subsequent segment or the
      // end-fill. This is the intended "hard stop" — an abrupt color boundary
      // with no interpolation — not a dropped stop.
      if (end > start) {
        // Pre-convert stop colors to linear float
        float pr = srgb_to_linear_float(prevColor.r / 255.0f);
        float pg = srgb_to_linear_float(prevColor.g / 255.0f);
        float pb = srgb_to_linear_float(prevColor.b / 255.0f);
        float nr = srgb_to_linear_float(nextColor.r / 255.0f);
        float ng = srgb_to_linear_float(nextColor.g / 255.0f);
        float nb = srgb_to_linear_float(nextColor.b / 255.0f);
        for (int i = start; i <= end; i++) {
          float t = static_cast<float>(i - start) / (end - start);

          // Interpolate in linear space
          float r_lin = pr * (1.0f - t) + nr * t;
          float g_lin = pg * (1.0f - t) + ng * t;
          float b_lin = pb * (1.0f - t) + nb * t;

          entries[i] = Pixel(float_to_pixel16(r_lin), float_to_pixel16(g_lin),
                             float_to_pixel16(b_lin));
        }
      }
      prevPos = nextPos;
      prevColor = nextColor;
      it++;
    }

    // Fill end. Same float conversion path as the start-fill and segment body
    // (see the start-fill note) for bit-identical flat/interpolated parity.
    int lastStop = static_cast<int>(prevPos * 255);
    Pixel lastLinear(float_to_pixel16(srgb_to_linear_float(prevColor.r / 255.0f)),
                     float_to_pixel16(srgb_to_linear_float(prevColor.g / 255.0f)),
                     float_to_pixel16(srgb_to_linear_float(prevColor.b / 255.0f)));
    for (int i = lastStop; i < 256; i++)
      entries[i] = lastLinear;
  }

  /**
   * @brief LUT lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The interpolated color (alpha 1.0).
   * @details The engine's premise is smooth 16-bit gradients, so a
   * nearest-index lookup would band visibly at low t resolution.
   */
  Color4 get(float t) const override {
    // Clamp before scaling: t > 1 would index past the table and t < 0 is UB for
    // the float->int conversion. After the clamp idx is in [0,255], so lo is a
    // valid index and lo+1 only overruns at the t==1 endpoint, handled below.
    float idx = hs::clamp(t, 0.0f, 1.0f) * 255.0f;
    int lo = static_cast<int>(idx);
    if (lo >= 255) return Color4(entries[255], 1.0f);
    float frac = idx - lo;
    // Round (+0.5f) the blend weight via the shared helper rather than
    // truncating, which would bias every interpolated sample down by up to
    // ~1/65535.
    return Color4(entries[lo].lerp16(entries[lo + 1], float_to_pixel16(frac)),
                  1.0f);
  }

  /**
   * @brief Destroys the gradient (LUT is an inline member; nothing to free).
   */
  ~Gradient() {}
};

/**
 * @brief A palette that generates colors based on defined harmony and
 * brightness profiles.
 */
class GenerativePalette : public Palette {
public:
  /**
   * @brief Default-constructs a straight, analogous, flat-brightness palette.
   */
  GenerativePalette()
      : GenerativePalette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                          BrightnessProfile::FLAT) {}

  /**
   * @brief Builds a palette from three pre-authored sRGB key colors.
   * @param gradient_shape Shape/distribution of colors across the domain.
   * @param ka First key color (stop a).
   * @param kb Second key color (stop b).
   * @param kc Third key color (stop c).
   * @details RNG-free: the keys are supplied directly rather than sampled from
   *          profiles. Used by from_hsv_keys (the tooling bridge) and any caller
   *          that already holds the three keys.
   */
  GenerativePalette(GradientShape gradient_shape, const CPixel &ka,
                    const CPixel &kb, const CPixel &kc)
      : gradient_shape(gradient_shape) {
    a = ka;
    b = kb;
    c = kc;
    update_stops();
  }

  /**
   * @brief Builds a 3-key palette from a base hue, harmony, and profiles.
   * @param gradient_shape Shape/distribution of colors across the domain.
   * @param harmony_type Harmony rule deriving the two companion hues.
   * @param profile Brightness profile across the domain.
   * @param sat_profile Saturation profile; defaults to MID.
   * @param manual_seed Base hue in [0, 255] when >= 0, else -1 to auto-seed.
   *        Supplying it is the opt-out from the shared `g_hue_seed` cursor: the
   *        global is then neither read nor advanced, giving an instance-local,
   *        deterministic base hue.
   * @details The base hue comes from manual_seed when >= 0, else from the global
   * golden-ratio hue cursor (which is then advanced).
   */
  GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type,
                    BrightnessProfile profile,
                    SaturationProfile sat_profile = SaturationProfile::MID,
                    int manual_seed = -1)
      : gradient_shape(gradient_shape) {
    uint8_t palette_hue;
    if (manual_seed != -1) {
      palette_hue = static_cast<uint8_t>(manual_seed);
    } else {
      palette_hue = g_hue_seed;
      // Advance the global seed by a fixed integer step: this is an R1 additive
      // recurrence (Weyl sequence) mod 256, x' = (x + s) mod 256. static_cast
      // truncates INV_PHI * 255 = 157.6 to s = 157. The truncation (not rounding
      // to 158) is load-bearing, not incidental: 157 is coprime with 256 so the
      // cursor visits every residue before repeating, whereas 158 is even and
      // would short-cycle. Calling it "low-discrepancy" overstates it — the true
      // low-discrepancy property of an R1 sequence needs the exact irrational
      // increment, which the truncation to 157 discards; equidistribution here
      // is only approximate. Kept as-is for output stability across generations.
      g_hue_seed = static_cast<uint8_t>((static_cast<uint32_t>(g_hue_seed) +
                                         static_cast<uint32_t>(INV_PHI * 255.0f)) %
                                        256);
    }

    uint8_t h1 = palette_hue;
    uint8_t h2, h3;

    calc_hues(h1, h2, h3, harmony_type);

    uint8_t s1 = 0, s2 = 0, s3 = 0;
    switch (sat_profile) {
    case SaturationProfile::PASTEL:
      s1 = s2 = s3 = 100;
      break;
    case SaturationProfile::MID:
      s1 = hs::rand_int(153, 204);
      s2 = hs::rand_int(153, 204);
      s3 = hs::rand_int(153, 204);
      break;
    case SaturationProfile::VIBRANT:
      s1 = s2 = s3 = 255;
      break;
    }

    uint8_t v1 = 0, v2 = 0, v3 = 0;
    switch (profile) {
    case BrightnessProfile::ASCENDING:
      v1 = hs::rand_int(25, 76);
      v2 = hs::rand_int(127, 178);
      v3 = hs::rand_int(204, 255);
      break;
    case BrightnessProfile::DESCENDING:
      v1 = hs::rand_int(204, 255);
      v2 = hs::rand_int(127, 178);
      v3 = hs::rand_int(25, 76);
      break;
    case BrightnessProfile::FLAT:
      v1 = v2 = v3 = 255;
      break;
    case BrightnessProfile::BELL:
      v1 = hs::rand_int(51, 127);
      v2 = hs::rand_int(178, 255);
      v3 = v1;
      break;
    case BrightnessProfile::CUP:
      v1 = hs::rand_int(178, 255);
      v2 = hs::rand_int(51, 127);
      v3 = v1;
      break;
    }

    // Author the three keys directly in OKLCH (see author_key / key_oklch). The
    // (h,s,v) integer draws above are kept untouched -- the global RNG stream,
    // and every other effect's randomness, are unperturbed -- and only
    // reinterpreted as OKLCH coordinates, then baked to the 8-bit sRGB key
    // (gamut-mapped) to preserve the 9-byte snapshot and the lerp machinery.
    a = author_key(h1, s1, v1);
    b = author_key(h2, s2, v2);
    c = author_key(h3, s3, v3);

    update_stops();
  }

  /**
   * @brief Builds a palette from three explicit, fully-resolved HSV key triples
   *        with NO RNG draws.
   * @param gradient_shape Shape/distribution of colors across the domain.
   * @param h1 First key hue in [0,255]; s1/v1 its saturation/value.
   * @param h2 Second key hue; s2/v2 its saturation/value.
   * @param h3 Third key hue; s3/v3 its saturation/value.
   * @return A palette whose keys are the OKLCH-authored (h,s,v) triples.
   * @details Mirrors the profile constructor's key authoring and stop layout
   *          exactly, but takes the (h,s,v) values as parameters instead of
   *          sampling them from the global RNG. This is the entry point the
   *          daydream palette tool's WASM bridge (PaletteOps) uses: the tool owns
   *          its own deterministic profile randomization and asks the engine only
   *          for the exact perceptual color math, so the two never drift.
   */
  static GenerativePalette from_hsv_keys(GradientShape gradient_shape,
                                         uint8_t h1, uint8_t s1, uint8_t v1,
                                         uint8_t h2, uint8_t s2, uint8_t v2,
                                         uint8_t h3, uint8_t s3, uint8_t v3) {
    return GenerativePalette(gradient_shape, author_key(h1, s1, v1),
                             author_key(h2, s2, v2), author_key(h3, s3, v3));
  }

  // Peak OKLCH chroma at mid-lightness. The saturation profile scales this and a
  // sin(pi*L) envelope (key_oklch / get()) tapers it toward the lightness
  // extremes where the sRGB gamut narrows. Tuned so a key at L~0.67 (FLAT) keeps
  // ~its pre-envelope chroma ((s/255)*0.23*sin(0.67pi) ~= (s/255)*0.20).
  static constexpr float kChromaPeak = 0.23f;

  // Hue torsion: radians of hue rotation per unit lightness, applied in get() as
  // h += kHueTorsion * (L - 0.5). Couples a small hue drift to lightness so
  // shadows and highlights shift along the ramp -- the painterly "lit from
  // within" look -- instead of every stop holding one fixed hue. 0.35 rad/L is
  // ~+/-10deg across the [0,1] lightness range (less within the authored band).
  static constexpr float kHueTorsion = 0.35f;

  /**
   * @brief OKLCH coordinates for one HSV-authored key. Single source of truth
   *        for the perceptual key placement.
   * @param hue Key hue in [0,255]; reinterpreted as an OKLCH hue angle.
   * @param sat Key saturation in [0,255]; scaled into OKLCH chroma.
   * @param val Key value in [0,255]; mapped into a perceptual lightness band.
   * @return The key's OKLCH coordinates (before the gamut-mapped sRGB bake).
   */
  static OKLCH key_oklch(uint8_t hue, uint8_t sat, uint8_t val) {
    // Hue: even perceptual spacing -- the integer harmony offsets from
    // calc_hues become true OKLCH-hue offsets (triadic is a real 120deg).
    float h = (hue / 256.0f) * (2.0f * PI_F);
    // Lightness: compress HSV value into a perceptual L band. The ceiling is
    // deliberately well below L=1: perceptual lightness near white starves the
    // sRGB gamut of chroma, so a high-value key (e.g. FLAT, v=255) would clip to
    // a washed-out pastel. Mapping [0,255] into [0.12, 0.67] lands full value in
    // the chroma-rich mid-high zone where the gamut bulges -- keys stay
    // saturated rather than bright-and-gray -- while the floor keeps dark stops
    // off pure black. FLAT is still a genuinely isoluminant shimmer, now at a
    // lightness the LEDs can actually render with saturation.
    float L = 0.12f + (val / 255.0f) * 0.55f;
    // Chroma co-varies with lightness: the saturation profile sets the ceiling
    // (kChromaPeak), and a sin(pi*L) envelope tapers it toward both ends, where
    // the gamut narrows. So a dark or near-white key can't also be fully
    // saturated -- the garish dark+saturated failure mode -- and chroma peaks at
    // mid-L where the gamut bulges. kChromaPeak is tuned so a typical key (FLAT
    // sits at L=0.67) keeps roughly its previous chroma. get() re-applies the
    // same envelope at the interpolated L so midpoints stay on this curve.
    //
    // VIBRANT is deliberately NOT authored at the exact per-hue cusp: get()'s
    // envelope peaks chroma at mid-L *between* keys (it overshoots a key's own
    // chroma when the segment passes through L~0.5), so cusp-authored keys push
    // that mid-segment peak past the gamut boundary -- the in-gamut/chroma-clip
    // path switch in oklch_to_pixel then bands the gradient. The fixed ceiling
    // leaves the headroom that overshoot needs. (The banding-safe chroma sits
    // about where this ceiling already is, so exact-cusp authoring bought almost
    // no extra vividness for the artifact it introduced.)
    float C = (sat / 255.0f) * kChromaPeak * sinf(PI_F * L);
    return {L, C, h};
  }

  /**
   * @brief Authors one key: OKLCH placement baked to a gamut-mapped 8-bit sRGB
   *        CPixel.
   * @param hue Key hue in [0,255].
   * @param sat Key saturation in [0,255].
   * @param val Key value in [0,255].
   * @return The authored key as an 8-bit sRGB CPixel.
   */
  static CPixel author_key(uint8_t hue, uint8_t sat, uint8_t val) {
    return oklch_to_cpixel(key_oklch(hue, sat, val));
  }

  /**
   * @brief Rebuilds stop positions/colors for the current gradient_shape.
   * @details Derives the stops from the three keys a/b/c, then caches their
   * OKLCH forms for get(). Call after any change to a/b/c or gradient_shape.
   */
  void update_stops() {
    const CPixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = {0, 0.1f, 0.5f, 0.9f, 1.0f};
      colors = {vignette_color, a, b, c, vignette_color};
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = {0, 0.5f, 1.0f};
      colors = {a, b, c};
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = {0, 0.33f, 0.66f, 1.0f};
      colors = {a, b, c, a};
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = {0, 0.33f, 0.66f, 0.9f, 1.0f};
      colors = {a, b, c, vignette_color, vignette_color};
      size = 5;
      break;
    }
    for (int i = 0; i < size; ++i) {
      colors_oklch[i] = srgb_to_oklch(colors[i].r, colors[i].g, colors[i].b);
      // Recover the stop's chroma ceiling from its authored (L, C): the key was
      // built as C = cmax * sin(pi*L), so cmax = C / sin(pi*L). get() re-applies
      // the envelope at the interpolated L. Guard the divisor at the L extremes
      // (black vignette stops sit at L~0, sin~0) -- there C~0, so cmax is ~0
      // anyway and the stop stays achromatic.
      float env = sinf(PI_F * colors_oklch[i].L);
      colors_cmax[i] = (env > 1e-3f) ? colors_oklch[i].C / env : 0.0f;
    }
  }

  /**
   * @brief Lightweight snapshot of just the 3 color keys (9 bytes).
   */
  struct Snapshot {
    CPixel a, b, c;
  };

  /**
   * @brief Captures the current three color keys.
   * @return A Snapshot of keys a, b, c.
   */
  Snapshot snapshot() const { return {a, b, c}; }

  /**
   * @brief Sets this palette's keys to the OKLCH interpolation of two palettes.
   * @param from Source palette at amount == 0.
   * @param to Source palette at amount == 1.
   * @param amount Blend weight.
   */
  void lerp(const GenerativePalette &from, const GenerativePalette &to,
            float amount) {
    lerp(from.snapshot(), to.snapshot(), amount);
  }

  /**
   * @brief Sets keys to the per-key OKLCH interpolation of two snapshots.
   * @param from Source snapshot at amount == 0.
   * @param to Source snapshot at amount == 1.
   * @param amount Blend weight.
   * @details Rebuilds the stops after interpolating.
   */
  void lerp(const Snapshot &from, const Snapshot &to, float amount) {
    // Enforce the [0,1] blend contract here rather than leaning on the
    // downstream RGB clamp: lerp_oklch can overshoot an extrapolated amount
    // into an invalid OKLCH (L>1 or C past gamut). Cold path (per-frame at
    // most); hs::clamp also folds NaN to 1.0 (house contract).
    amount = hs::clamp(amount, 0.0f, 1.0f);
    a = lerp_oklch_srgb(from.a, to.a, amount);
    b = lerp_oklch_srgb(from.b, to.b, amount);
    c = lerp_oklch_srgb(from.c, to.c, amount);
    update_stops();
  }

  /**
   * @brief Samples the generated palette at a coordinate.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to 1.0).
   * @return The color at t (alpha 1.0).
   */
  Color4 get(float t) const override {
    // Clamp to [0,1] like Gradient::get / BakedPalette::get. Without this a
    // t < shape[0] (=0) matches no segment and falls through to the size-2
    // fallback below — returning the LAST segment's start (the middle color for
    // STRAIGHT) instead of the first stop, a discontinuity at t=0. hs::clamp
    // also folds NaN to 1.0 (the house contract), so a NaN t yields the last
    // stop deterministically rather than the same stray fallback.
    t = hs::clamp(t, 0.0f, 1.0f);
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0)
      seg = size - 2;

    float start = shape[seg];
    float end = shape[seg + 1];
    CPixel c1 = colors[seg];

    float dist = end - start;
    if (dist < 0.0001f)
      return Color4(c1, 1.0f);

    float p = std::clamp((t - start) / dist, 0.0f, 1.0f);

    // Interpolate in OKLCH for perceptually uniform gradients. Stops are
    // pre-converted in update_stops(), so this avoids the per-sample
    // sRGB->OKLCH cost (load-bearing on GSReactionDiffusion's 4x-SSAA path).
    OKLCH blended = lerp_oklch(colors_oklch[seg], colors_oklch[seg + 1], p);
    // Chroma rides the lightness envelope at the *interpolated* L: interpolate
    // the two stops' chroma ceilings, then C = cmax * sin(pi*L). This keeps a
    // segment that spans a lightness change on the sin curve rather than the
    // straight chord lerp_oklch produces -- no chroma dip where the gradient
    // brightens through mid-L. fast_sinf: get() is per-sample on the non-baked
    // palette path (the +1.86% near-zero error is negligible -- chroma is ~0
    // there). lerp_oklch still owns the L and shortest-arc hue; only C changes.
    float cmax = colors_cmax[seg] + (colors_cmax[seg + 1] - colors_cmax[seg]) * p;
    blended.C = cmax * fast_sinf(PI_F * blended.L);
    // Hue torsion: drift hue with lightness (centered at L=0.5) so the ramp
    // reads as "lit from within" rather than a constant-hue lerp. One multiply-
    // add; oklch_to_oklab takes cos/sin of the angle, so no explicit wrap needed.
    blended.h += kHueTorsion * (blended.L - 0.5f);
    return Color4(oklch_to_pixel(blended), 1.0f);
  }

private:
  /**
   * @brief Wraps an integer hue into [0,255], correct for negative inputs.
   * @param hue Hue value to wrap; may be negative.
   * @return The hue reduced to [0, 255].
   */
  uint8_t wrap_hue(int hue) const { return (hue % 256 + 256) % 256; }

  /**
   * @brief Derives the two companion hues from base hue h1 per a harmony.
   * @param h1 Base hue in [0, 255].
   * @param h2 Out: first companion hue.
   * @param h3 Out: second companion hue.
   * @param harmony_type Harmony rule selecting the offsets.
   * @details Some harmonies add randomized jitter, so output is not pure.
   */
  void calc_hues(uint8_t h1, uint8_t &h2, uint8_t &h3,
                 HarmonyType harmony_type) const {
    const int h1_int = h1;
    switch (harmony_type) {
    case HarmonyType::TRIADIC:
      h2 = wrap_hue(h1_int + 85);
      h3 = wrap_hue(h1_int + 170);
      break;
    case HarmonyType::SPLIT_COMPLEMENTARY: {
      const int complement = wrap_hue(h1_int + 128);
      const int offset = 21;
      h2 = wrap_hue(complement - offset);
      h3 = wrap_hue(complement + offset);
      break;
    }
    case HarmonyType::COMPLEMENTARY: {
      h2 = wrap_hue(h1_int + 128);
      const int offset = hs::rand_int(-7, 8);
      h3 = wrap_hue(h1_int + offset);
      break;
    }
    case HarmonyType::ANALOGOUS:
    default: {
      const int dir = (hs::rand_int(0, 2) == 0) ? 1 : -1;
      const int offset1 = dir * hs::rand_int(11, 22);
      h2 = wrap_hue(h1_int + offset1);
      const int offset2 = dir * hs::rand_int(11, 22);
      h3 = wrap_hue(h2 + offset2);
      break;
    }
    }
  }

  GradientShape gradient_shape;
  CPixel a, b, c;

  // Mutable static cursor shared across all GenerativePalette instances: each
  // auto-seeded construction reads it and advances it so successive palettes get
  // distinct base hues. Non-atomic by design — safe only because palette
  // construction is single-threaded (engine setup / the render thread); it is
  // NOT a concurrency guard and a second constructing thread would race it.
  // OPT-OUT: the global is engaged only on the auto-seed path. Passing an
  // explicit `manual_seed >= 0` to the constructor neither reads nor advances
  // this cursor, so a caller wanting a deterministic, instance-local base hue
  // (or isolation from the shared cursor) seeds through the constructor instead.
  static inline uint8_t g_hue_seed = 0;
  std::array<float, 5> shape;
  std::array<CPixel, 5> colors;
  /**
   * @brief OKLCH forms of `colors`, cached by update_stops().
   * @details Lets the per-sample get() hot path skip the sRGB->OKLCH conversion
   * (6 powf + 6 cbrtf + 2 atan2f per stop). The stops change only on
   * construction and in lerp(), both of which route through update_stops(), so
   * this stays in sync.
   */
  std::array<OKLCH, 5> colors_oklch;
  /**
   * @brief Per-stop chroma ceiling, cached by update_stops().
   * @details Each stop's chroma is authored as C = cmax * sin(pi*L) (the
   * lightness-coupled envelope, see key_oklch). This recovers cmax = C/sin(pi*L)
   * so get() can re-apply the envelope at the interpolated L, keeping midpoints
   * on the sin curve rather than the straight chroma chord a plain lerp gives.
   */
  std::array<float, 5> colors_cmax{};
  int size = 0;

public:
  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~GenerativePalette() {}

  /**
   * @brief Pins the global generative-hue cursor.
   * @param seed Hue value in [0, 255] to set the cursor to; defaults to 0.
   * @details Auto-seeded palettes draw their base hue from `g_hue_seed` and
   * advance it (golden-ratio step) so successive palettes across the live show
   * keep evolving — deliberately stateful, never reset in production. The test
   * harness calls this to restore identical global state before re-rendering an
   * effect, so a cross-run determinism check is not defeated by the drift. Does
   * not affect production, which never invokes it.
   */
  static void reset_hue_seed(uint8_t seed = 0) { g_hue_seed = seed; }
};

/**
 * @brief A palette defined by a mathematical cosine wave function.
 * C(t) = A + B * cos(2 * PI * (C * t + D))
 */
class ProceduralPalette : public Palette {
public:
  /**
   * @brief Default-constructs a palette with all-zero cosine coefficients.
   */
  constexpr ProceduralPalette()
      : a{0, 0, 0}, b{0, 0, 0}, c{0, 0, 0}, d{0, 0, 0} {}
  /**
   * @brief Constructs from the four cosine coefficient vectors.
   * @param a Bias term per RGB channel.
   * @param b Amplitude per RGB channel.
   * @param c Frequency per RGB channel.
   * @param d Phase per RGB channel.
   */
  constexpr ProceduralPalette(std::array<float, 3> a, std::array<float, 3> b,
                              std::array<float, 3> c, std::array<float, 3> d)
      : a(a), b(b), c(c), d(d) {}

  /**
   * @brief Evaluates the cosine palette at a coordinate.
   * @param t Lookup coordinate.
   * @return The color at t (alpha 1.0).
   * @details Determines color in high-precision float sRGB space first, then
   * converts to 16-bit linear via the interpolated LUT — avoids the 8-bit
   * quantization a plain srgb_to_linear(uint8_t) lookup would impose, without a
   * per-channel powf.
   */
  Color4 get(float t) const override {
    // Determine color in high-precision float sRGB space first, then convert to
    // 16-bit linear via the interpolated LUT — avoids the 8-bit quantization a
    // plain srgb_to_linear(uint8_t) lookup would impose, without a per-channel
    // powf.
    float r_srgb =
        hs::clamp(a[0] + b[0] * cosf(2 * PI_F * (c[0] * t + d[0])), 0.0f, 1.0f);
    float g_srgb =
        hs::clamp(a[1] + b[1] * cosf(2 * PI_F * (c[1] * t + d[1])), 0.0f, 1.0f);
    float b_srgb =
        hs::clamp(a[2] + b[2] * cosf(2 * PI_F * (c[2] * t + d[2])), 0.0f, 1.0f);

    Pixel color(srgb_to_linear_interp(r_srgb), srgb_to_linear_interp(g_srgb),
                srgb_to_linear_interp(b_srgb));
    return Color4(color, 1.0f);
  }

protected:
  std::array<float, 3> a, b, c, d;

public:
  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~ProceduralPalette() {}
};

/**
 * @brief A palette that allows continuous mutation between two procedural
 * palettes.
 */
class MutatingPalette : public ProceduralPalette {
public:
  /**
   * @brief Constructs from two endpoint cosine parameter sets.
   * @param a1 Start bias per channel.
   * @param b1 Start amplitude per channel.
   * @param c1 Start frequency per channel.
   * @param d1 Start phase per channel.
   * @param a2 End bias per channel.
   * @param b2 End amplitude per channel.
   * @param c2 End frequency per channel.
   * @param d2 End phase per channel.
   * @details Initializes the active parameters to the start set (mutate(0)).
   */
  MutatingPalette(std::array<float, 3> a1, std::array<float, 3> b1,
                  std::array<float, 3> c1, std::array<float, 3> d1,
                  std::array<float, 3> a2, std::array<float, 3> b2,
                  std::array<float, 3> c2, std::array<float, 3> d2)
      : ProceduralPalette(a1, b1, c1, d1), // Initialize base with start values
        a1(a1), b1(b1), c1(c1), d1(d1), a2(a2), b2(b2), c2(c2), d2(d2) {
    mutate(0.0f);
  }

  /**
   * @brief Sets the active cosine parameters to the endpoint interpolation.
   * @param t Blend weight in [0, 1] between the start and end parameter sets.
   */
  void mutate(float t) {
    for (int i = 0; i < 3; ++i) {
      a[i] = lerp(a1[i], a2[i], t);
      b[i] = lerp(b1[i], b2[i], t);
      c[i] = lerp(c1[i], c2[i], t);
      d[i] = lerp(d1[i], d2[i], t);
    }
  }

private:
  /**
   * @brief Linearly interpolates between two scalars.
   * @param x Value at t == 0.
   * @param y Value at t == 1.
   * @param t Blend weight.
   * @return The interpolated value.
   */
  float lerp(float x, float y, float t) { return x * (1.0f - t) + y * t; }
  std::array<float, 3> a1, b1, c1, d1;
  std::array<float, 3> a2, b2, c2, d2;

public:
  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~MutatingPalette() {}
};

///////////////////////////////////////////////////////////////////////////////
// Palette Modifiers
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Linearly cycles the palette coordinate.
 *
 * The offset driver is optional by design (defaults to null): a null driver is
 * a deliberate "no cycling, static palette" mode, not an error — so modify()
 * passes t through. Contrast BreatheModifier/RippleModifier, whose driver is
 * mandatory and trapped at construction.
 */
struct CycleModifier {
  const float *offset;

  /**
   * @brief Constructs with an optional offset driver.
   * @param driver_offset Pointer to the per-frame offset, or null for static.
   */
  CycleModifier(const float *driver_offset = nullptr) : offset(driver_offset) {}

  /**
   * @brief Shifts the coordinate by the driver offset (pass-through if null).
   * @param t Input coordinate.
   * @return t plus the offset, or t unchanged when no driver is bound.
   */
  float modify(float t) const { return offset ? t + *offset : t; }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier {
  const float *phase;
  float amplitude;
  /**
   * @brief Per-instance memo of fast_sinf(*phase).
   * @details *phase is a frame-constant driver, so the sine is recomputed only
   * when it advances (once per frame) rather than for every pixel. mutable:
   * modify() is const but the memo is not part of the modifier's observable
   * value.
   */
  mutable float cached_phase_ = 0.0f;
  mutable float cached_sin_ = 0.0f;   /**< Memoized sine of cached_phase_. */
  mutable bool primed_ = false;       /**< Whether the memo has been populated. */

  /**
   * @brief Constructs with a mandatory phase driver and amplitude.
   * @param driver_phase Pointer to the per-frame phase; must not be null.
   * @param amp Oscillation amplitude; defaults to 0.1.
   * @details The phase driver is mandatory (no default): a null one would
   * silently freeze the breathe with no diagnostic. Trap at construction (cold
   * path) so the per-pixel modify() can dereference unconditionally.
   */
  BreatheModifier(const float *driver_phase, float amp = 0.1f)
      : phase(driver_phase), amplitude(amp) {
    HS_CHECK(phase, "BreatheModifier: phase driver must not be null");
  }

  /**
   * @brief Oscillates the coordinate by amplitude * sin(phase).
   * @param t Input coordinate.
   * @return t plus the memoized oscillation term.
   */
  float modify(float t) const {
    if (!primed_ || *phase != cached_phase_) {
      cached_phase_ = *phase;
      cached_sin_ = fast_sinf(*phase);
      primed_ = true;
    }
    return t + cached_sin_ * amplitude;
  }
};

/**
 * @brief Distorts the palette spatially with a sine wave, creating a liquid
 * ripple effect. When applied to a spatial coordinate, colors will compress and
 * expand like waves.
 */
struct RippleModifier {
  const float *phase;
  float frequency;
  float amplitude;

  /**
   * @brief Constructs with a mandatory phase driver, frequency, and amplitude.
   * @param phase Pointer to the per-frame phase; must not be null.
   * @param freq Spatial frequency of the ripple; defaults to 3.0.
   * @param amp Distortion amplitude; defaults to 0.1.
   * @details Mandatory phase driver (no default) — trap a null one at
   * construction rather than silently passing t through on every pixel.
   */
  RippleModifier(const float *phase, float freq = 3.0f, float amp = 0.1f)
      : phase(phase), frequency(freq), amplitude(amp) {
    HS_CHECK(phase, "RippleModifier: phase driver must not be null");
  }

  /**
   * @brief Distorts the coordinate with a sine wave of the given frequency.
   * @param t Input coordinate.
   * @return t plus the local sine distortion.
   */
  float modify(float t) const {
    // Calculate a local distortion based on the current t position. The sine
    // argument is per-pixel (depends on t), so unlike BreatheModifier it cannot
    // be memoized; use fast_sinf for the per-pixel path, matching the rest of
    // the per-pixel trig in this file.
    return t + fast_sinf(t * frequency * PI_F * 2.0f + *phase) * amplitude;
  }
};

/**
 * @brief Folds the palette back and forth like a kaleidoscope.
 * A folds value of 2.0 maps [0...1] to [1 -> 0 -> 1] (one full bounce);
 * each unit of folds adds another half-bounce.
 *
 * The phase driver is optional by design (defaults to null) — a null driver is
 * the deliberate "no phase offset" mode (shift = 0), not an error.
 */
struct FoldModifier {
  const float *phase;
  float folds;

  /**
   * @brief Constructs with a fold count and optional phase driver.
   * @param folds Number of bounces; defaults to 2.0 (one full bounce).
   * @param phase Pointer to an optional phase offset, or null for none.
   */
  FoldModifier(float folds = 2.0f, const float *phase = nullptr)
      : phase(phase), folds(folds) {}

  /**
   * @brief Folds the coordinate back and forth via a triangle wave.
   * @param t Input coordinate.
   * @return The folded coordinate in [0, 1].
   */
  float modify(float t) const {
    float shift = phase ? *phase : 0.0f;
    float scaled = (t * folds) + shift;

    // Triangle wave formula: creates symmetrical, continuous bouncing between 0
    // and 1. fmodf keeps the sign of the dividend, so reduce into [0, 2) first
    // — otherwise negative scaled (e.g. a negative phase shift) folds to values
    // above 1.
    float m = fmodf(scaled, 2.0f);
    if (m < 0.0f) m += 2.0f;
    return fabsf(m - 1.0f);
  }
};

/**
 * @brief Pinches or expands the center of the palette.
 * positive tension pulls colors toward the center, negative pushes them to the
 * edges. Creates an elastic, tension-release visual dynamic.
 *
 * The tension driver is optional by design (defaults to null) — a null driver is
 * the deliberate "no pinch" pass-through mode, not an error (matching
 * CycleModifier/FoldModifier).
 */
struct PinchModifier {
  const float *tension; /**< Pinch tension driver; expects roughly -0.9 to 0.9. */

  /**
   * @brief Constructs with an optional tension driver.
   * @param t Pointer to the tension value, or null for pass-through.
   */
  PinchModifier(const float *t = nullptr) : tension(t) {}

  /**
   * @brief Pinches or expands the coordinate around the domain center.
   * @param t Input coordinate.
   * @return The reshaped coordinate, or t unchanged when no driver is bound.
   */
  float modify(float t) const {
    if (!tension)
      return t;

    // Shift t to -1.0 to 1.0 range based on a 0.0 to 1.0 domain block
    float wrapped_t = wrap_t(t);
    float centered = wrapped_t * 2.0f - 1.0f;
    float sign = centered < 0.0f ? -1.0f : 1.0f;

    // Apply power curve
    float amount = std::clamp(*tension, -0.99f, 0.99f);
    float power = (amount < 0.0f) ? (1.0f / (1.0f + std::abs(amount)))
                                  : (1.0f + amount * 3.0f);

    centered = sign * powf(std::abs(centered), power);

    // Map back to 0.0 to 1.0 and add to integer spatial frame. floorf(t) pairs
    // with wrap_t(t) (== t - floorf(t)), so the reshaped fraction is re-anchored
    // to t's own integer cell — including for negative t, where both agree on
    // the same floor (no truncation-toward-zero mismatch).
    return floorf(t) + ((centered + 1.0f) * 0.5f);
  }
};

/**
 * @brief Snaps smooth gradients into harsh, distinct bands (Posterization).
 * By animating the step count or offsetting before quantizing, you get a
 * retro/glitch aesthetic.
 */
struct QuantizeModifier {
  const float *dynamic_steps;
  float base_steps;

  /**
   * @brief Constructs with a base step count and optional dynamic driver.
   * @param steps Base number of quantization steps.
   * @param d_steps Pointer to an animated step count, or null to use base.
   */
  QuantizeModifier(float steps, const float *d_steps = nullptr)
      : dynamic_steps(d_steps), base_steps(steps) {}

  /**
   * @brief Snaps the coordinate to the nearest of N steps.
   * @param t Input coordinate.
   * @return The quantized coordinate.
   */
  float modify(float t) const {
    float s = dynamic_steps ? *dynamic_steps : base_steps;
    if (s < 1.0f)
      s = 1.0f;

    // Round to nearest step in the infinite domain.
    return roundf(t * s) / s;
  }
};

/**
 * @brief Multiplies the palette coordinate, increasing the frequency
 * so the palette repeats multiple times across the domain.
 */
struct ScaleModifier {
  const float *dynamic_scale;
  float base_scale;

  /**
   * @brief Constructs with a base scale and optional dynamic driver.
   * @param s Base scale factor; defaults to 1.0.
   * @param d_scale Pointer to an animated scale, or null to use base.
   */
  ScaleModifier(float s = 1.0f, const float *d_scale = nullptr)
      : dynamic_scale(d_scale), base_scale(s) {}

  /**
   * @brief Multiplies the coordinate by the active scale.
   * @param t Input coordinate.
   * @return The scaled coordinate.
   */
  float modify(float t) const {
    return t * (dynamic_scale ? *dynamic_scale : base_scale);
  }
};

/**
 * @brief Reverses the palette coordinate (t -> 1 - t).
 */
struct ReverseModifier {
  /**
   * @brief Reverses the coordinate.
   * @param t Input coordinate.
   * @return 1 - t.
   */
  float modify(float t) const { return 1.0f - t; }
};

/**
 * @brief Mirrors the coordinate so [0,1] maps to [0,1,0].
 * @details One symmetric bounce, for a seamless loop.
 */
struct MirrorModifier {
  /**
   * @brief Mirrors the coordinate into a symmetric bounce.
   * @param t Input coordinate.
   * @return The mirrored coordinate in [0, 1].
   */
  float modify(float t) const { return 1.0f - fabsf(2.0f * t - 1.0f); }
};

/**
 * @brief Compresses the source domain into an inset window [lo, hi] -> [0, 1].
 * @details Clamps outside so t below lo samples the first stop and t above hi
 * the last. Pairs with EdgeFadeShade / EdgeAlphaShade to build vignettes.
 */
struct InsetModifier {
  float lo, hi;
  /**
   * @brief Constructs the inset window bounds.
   * @param lo Lower domain bound mapped to 0; defaults to 0.2.
   * @param hi Upper domain bound mapped to 1; defaults to 0.8.
   */
  InsetModifier(float lo = 0.2f, float hi = 0.8f) : lo(lo), hi(hi) {
    HS_CHECK(hi > lo, "InsetModifier: hi must be > lo (modify divides by hi - lo)");
  }
  /**
   * @brief Remaps the coordinate from [lo, hi] into [0, 1], clamping outside.
   * @param t Input coordinate.
   * @return The remapped coordinate in [0, 1].
   */
  float modify(float t) const {
    return hs::clamp((t - lo) / (hi - lo), 0.0f, 1.0f);
  }
};

///////////////////////////////////////////////////////////////////////////////
// Color Modifiers — reshape the sample after the source lookup.
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Scales alpha by a caller-supplied falloff curve over the coordinate.
 */
struct AlphaFalloffShade {
  using FalloffFunction = float (*)(float);
  FalloffFunction fn;
  /**
   * @brief Constructs with the falloff function.
   * @param fn Function mapping a coordinate to an alpha multiplier.
   */
  AlphaFalloffShade(FalloffFunction fn) : fn(fn) {}
  /**
   * @brief Scales the sample's alpha by the falloff curve at t.
   * @param c Sample color to reshape.
   * @param t Coordinate passed to the falloff function.
   * @return The sample with alpha scaled.
   */
  Color4 shade(Color4 c, float t) const {
    c.alpha *= fn(t);
    return c;
  }
};

/**
 * @brief Fades the sample color to black near the coordinate edges.
 * @details Opaque vignette. Pair with InsetModifier so the edge bands resolve
 * to the source's first/last stop before fading.
 */
struct EdgeFadeShade {
  float edge;
  /**
   * @brief Constructs with the edge fade width.
   * @param edge Fraction of the domain over which each edge fades; default 0.2.
   */
  EdgeFadeShade(float edge = 0.2f) : edge(edge) {
    HS_CHECK(edge > 0.0f, "EdgeFadeShade: edge must be > 0 (shade divides by it)");
  }
  /**
   * @brief Fades the sample color toward black within the edge bands.
   * @param c Sample color to reshape.
   * @param t Coordinate in [0, 1].
   * @return The sample with its color faded near the edges.
   */
  Color4 shade(Color4 c, float t) const {
    // Pixel (16-bit linear) black, not CRGB: a CRGB temporary here would bind
    // the low-edge blend to CRGB::lerp16 (quantize -> sRGB lerp -> re-expand),
    // banding the fade and making it asymmetric with the 16-bit high edge.
    Pixel black(0, 0, 0);
    if (t < edge)
      return Color4(
          black.lerp16(c.color, float_to_pixel16(quintic_kernel(t / edge))),
          c.alpha);
    if (t >= 1.0f - edge)
      return Color4(c.color.lerp16(black, float_to_pixel16(quintic_kernel(
                                              (t - (1.0f - edge)) / edge))),
                    c.alpha);
    return c;
  }
};

/**
 * @brief Fades the sample alpha (not color) near the coordinate edges.
 * @details Transparent vignette. Pair with InsetModifier as with EdgeFadeShade.
 */
struct EdgeAlphaShade {
  float edge;
  /**
   * @brief Constructs with the edge fade width.
   * @param edge Fraction of the domain over which each edge fades; default 0.2.
   */
  EdgeAlphaShade(float edge = 0.2f) : edge(edge) {
    HS_CHECK(edge > 0.0f, "EdgeAlphaShade: edge must be > 0 (shade divides by it)");
  }
  /**
   * @brief Fades the sample alpha within the edge bands.
   * @param c Sample color to reshape.
   * @param t Coordinate in [0, 1].
   * @return The sample with its alpha faded near the edges.
   */
  Color4 shade(Color4 c, float t) const {
    if (t < edge)
      c.alpha *= quintic_kernel(t / edge);
    else if (t >= 1.0f - edge)
      c.alpha *= quintic_kernel(1.0f - (t - (1.0f - edge)) / edge);
    return c;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Compile-Time Palette Composition
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Concept for a coordinate modifier.
 * @tparam T Type required to expose a const modify(float)->float method.
 * @details Remaps the lookup coordinate before the source is sampled.
 */
template <typename T>
concept CoordMod = requires(const T m, float t) {
  { m.modify(t) } -> std::convertible_to<float>;
};

/**
 * @brief Concept for a color modifier.
 * @tparam T Type required to expose a const shade(Color4, float)->Color4 method.
 * @details Reshapes the sample after the lookup, with the original coordinate
 * in hand.
 */
template <typename T>
concept ColorMod = requires(const T m, Color4 c, float t) {
  { m.shade(c, t) } -> std::convertible_to<Color4>;
};

/**
 * @brief Type-list tag for the coordinate-modifier axis of a StaticPalette.
 * @tparam M Coordinate modifier types.
 */
template <typename... M> struct Coords {};
/**
 * @brief Type-list tag for the color-modifier axis of a StaticPalette.
 * @tparam M Color modifier types.
 */
template <typename... M> struct Colors {};

/**
 * @brief A compile-time composition of a Source palette, a coordinate-modifier
 * chain, and a color-modifier chain.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CoordList Coords<> type-list of coordinate modifiers.
 * @tparam ColorList Colors<> type-list of color modifiers.
 * @tparam Wrap Whether to wrap the coordinate before sampling the source.
 * @details Default constructible to allow safe wiring in init(). Both chains are
 * inlined by fold expression with zero runtime overhead, following the
 * ArenaVector idiom: default construct, then bind(). get() applies the coord
 * mods to t (in order), samples the source (wrapping the coordinate unless
 * Wrap is false), then applies the color mods to the sample with the *original*
 * coordinate. Wrap=false suits inset/falloff pipelines that must reach the
 * source's exact endpoints (wrap_t(1)==0 would otherwise fold the top edge).
 */
template <typename Source, typename CoordList = Coords<>,
          typename ColorList = Colors<>, bool Wrap = true>
class StaticPalette;

/**
 * @brief Partial specialization splitting the two modifier type-lists.
 * @tparam Source Source palette type exposing Color4 get(float) const.
 * @tparam CMods Coordinate modifier types.
 * @tparam XMods Color modifier types.
 * @tparam Wrap Whether to wrap the coordinate before sampling the source.
 */
template <typename Source, typename... CMods, typename... XMods, bool Wrap>
class StaticPalette<Source, Coords<CMods...>, Colors<XMods...>, Wrap> {
  static_assert((CoordMod<CMods> && ...), "Coords<> entries must be CoordMods");
  static_assert((ColorMod<XMods> && ...), "Colors<> entries must be ColorMods");

public:
  /**
   * @brief Default-constructs an unbound composition (bind() before use).
   */
  StaticPalette() = default;

  /**
   * @brief Binds the source and modifier chains by pointer.
   * @param src Source palette; must not be null.
   * @param cms Coordinate-modifier pointers, one per CMods entry; none null.
   * @param xms Color-modifier pointers, one per XMods entry; none null.
   * @details Cold wiring (init only): get() guards source_ with a debug-only
   * assert (stripped on-device), and a null member read does not fault on
   * Teensy 4.x, so the whole chain is trapped always-on here at the cold seam
   * (empty packs fold to true) for zero hot-path cost.
   */
  void bind(const Source *src, const CMods *...cms, const XMods *...xms) {
    HS_CHECK(src != nullptr, "StaticPalette bound to null source");
    HS_CHECK(((cms != nullptr) && ...), "StaticPalette bound to null coord modifier");
    HS_CHECK(((xms != nullptr) && ...), "StaticPalette bound to null color modifier");
    source_ = src;
    coords_ = std::make_tuple(cms...);
    colors_ = std::make_tuple(xms...);
  }

  /**
   * @brief Applies the coord chain, samples the source, then the color chain.
   * @param t Lookup coordinate.
   * @return The fully modified color.
   * @details The coord mods remap t in order; the source is sampled (wrapping
   * the coordinate unless Wrap is false); then the color mods reshape the
   * sample with the *original* coordinate.
   */
  Color4 get(float t) const {
    // Per-pixel hot path: debug-only guard, deliberately NOT HS_CHECK (an
    // always-on branch here costs on every pixel on-device).
    assert(source_ != nullptr && "StaticPalette used before bind()!");

    float ft = t;
    std::apply([&](const auto *...m) { ((ft = m->modify(ft)), ...); }, coords_);

    float u = ft;
    if constexpr (Wrap)
      u = wrap_t(ft);
    Color4 c = source_->get(u);

    std::apply([&](const auto *...m) { ((c = m->shade(c, t)), ...); }, colors_);
    return c;
  }

private:
  const Source *source_ = nullptr;
  std::tuple<const CMods *...> coords_{};
  std::tuple<const XMods *...> colors_{};
};

/**
 * @brief Runtime Palette facade over a compile-time StaticPalette composition.
 * @tparam SP StaticPalette composition type exposing Color4 get(float) const.
 * @details Bridges a zero-overhead StaticPalette into the polymorphic
 * `const Palette*` world (preset tables, BakedPalette::bake). The virtual call
 * is paid only where a Palette* is required — at bake time (cold) — never on the
 * per-pixel path, which runs off the resulting BakedPalette LUT.
 */
template <typename SP> class PaletteFacade : public Palette {
public:
  /**
   * @brief Default-constructs an unbound facade (bind() before use).
   */
  PaletteFacade() = default;
  /**
   * @brief Constructs a facade bound to a composition.
   * @param sp Composition to forward get() to.
   */
  explicit PaletteFacade(const SP *sp) : sp_(sp) {}
  /**
   * @brief Binds the facade to a composition.
   * @param sp Composition to forward get() to; must not be null.
   */
  void bind(const SP *sp) {
    HS_CHECK(sp != nullptr, "PaletteFacade bound to null composition");
    sp_ = sp;
  }
  /**
   * @brief Forwards the lookup to the bound composition.
   * @param t Lookup coordinate.
   * @return The composition's color at t.
   */
  Color4 get(float t) const override {
    assert(sp_ != nullptr && "PaletteFacade used before bind()!");
    return sp_->get(t);
  }

private:
  const SP *sp_ = nullptr;
};

/**
 * @brief Palette that returns one fixed color for every coordinate.
 */
class SolidColorPalette : public Palette {
public:
  /**
   * @brief Constructs with the fixed color.
   * @param color Color returned for every lookup.
   */
  SolidColorPalette(const Color4 &color) : color(color) {}
  /**
   * @brief Returns the fixed color.
   * @return The stored color, regardless of coordinate.
   */
  Color4 get(float) const override { return color; }
  Color4 color;
};

/**
 * @brief Pre-baked 256-entry Pixel16 LUT allocated in an arena.
 * @details Converts any Palette into a fast table lookup with lerp
 * interpolation. Not a Palette subclass — call get(t) directly for
 * zero-overhead lookups.
 */
class BakedPalette {
public:
  static constexpr int LUT_SIZE = 256;

  /**
   * @brief Default-constructs an unbaked palette (bake() before use).
   */
  BakedPalette() = default;

  /**
   * @brief Bakes any source into a 256-entry LUT in the given arena.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param arena Arena to allocate the LUT from.
   * @param source Source palette or composition to sample.
   * @details Works for a runtime Palette or a compile-time StaticPalette alike.
   */
  template <typename Source> void bake(Arena &arena, const Source &source) {
    lut_ = static_cast<Color4 *>(
        arena.allocate(LUT_SIZE * sizeof(Color4), alignof(Color4)));
    rebake(source);
  }

  /**
   * @brief Refills the existing LUT without allocating. Use for animated palettes.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param source Source palette or composition to sample.
   */
  template <typename Source> void rebake(const Source &source) {
    // Cold (per-frame at most, never per-pixel): a null lut_ here means rebake()
    // was called before bake() — a wiring bug with no valid recovery. Trap
    // always-on (survives NDEBUG on-device) rather than write through null.
    HS_CHECK(lut_ != nullptr, "BakedPalette::rebake before bake()");
    for (int i = 0; i < LUT_SIZE; ++i) {
      float t = static_cast<float>(i) / (LUT_SIZE - 1);
      lut_[i] = source.get(t);
    }
  }

  /**
   * @brief Fast lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to the last entry).
   * @return The interpolated color.
   */
  Color4 get(float t) const {
    // Per-pixel hot path: debug-only guard (parity with StaticPalette::get),
    // deliberately NOT HS_CHECK — a branch here costs on every pixel on-device.
    // The unbaked-palette wiring bug is trapped always-on at the cold rebake().
    assert(lut_ != nullptr && "BakedPalette::get before bake()");
    // Clamp in float space before the int cast (parity with the NaN-safe
    // Gradient::get): static_cast<int>(NaN) is UB and a NaN idx would slip past
    // both integer bound guards below. hs::clamp is a branchless hardware
    // min/max that maps NaN to the hi bound, so NaN resolves to the last entry.
    // Valid t in [0,1] is unchanged (no-op clamp), so the per-pixel hot path
    // keeps its cost. The clamp guarantees idx >= 0, so no `lo < 0` guard is
    // needed.
    float idx = hs::clamp(t * (LUT_SIZE - 1), 0.0f,
                          static_cast<float>(LUT_SIZE - 1));
    int lo = static_cast<int>(idx);
    if (lo >= LUT_SIZE - 1) return lut_[LUT_SIZE - 1];
    float frac = idx - lo;
    const Color4 &a = lut_[lo];
    const Color4 &b = lut_[lo + 1];
    // Round (+0.5f) the blend weight rather than truncating, matching
    // float_to_pixel16 / Gradient::get — bare truncation would bias every
    // interpolated sample down by up to ~1/65535. The helper's clamp is skipped
    // here: frac is already in [0,1) after the bounds checks above, and get() is
    // a per-pixel hot path. The +0.5f round applies only to the color channel:
    // its weight is quantized to a uint16_t, so truncation would bias it. Alpha
    // is interpolated with the raw float `frac` in full float precision (no
    // quantization step to bias), so it deliberately needs no rounding term.
    return Color4(a.color.lerp16(b.color,
                                 static_cast<uint16_t>(frac * 65535.0f + 0.5f)),
                  a.alpha + (b.alpha - a.alpha) * frac);
  }

  /**
   * @brief Deep-copies the LUT from another BakedPalette into the given arena.
   * @param src Source palette to copy; must already be baked.
   * @param arena Arena to allocate the new LUT from.
   * @details Used by Persist for arena compaction.
   */
  void clone_from(const BakedPalette &src, Arena &arena) {
    // Cold path (Persist arena compaction): cloning a never-baked source would
    // memcpy from null (UB). Trap always-on like rebake() rather than corrupt.
    HS_CHECK(src.lut_ != nullptr, "BakedPalette::clone_from before src bake()");
    lut_ = static_cast<Color4 *>(
        arena.allocate(LUT_SIZE * sizeof(Color4), alignof(Color4)));
    memcpy(lut_, src.lut_, LUT_SIZE * sizeof(Color4));
  }

private:
  Color4 *lut_ = nullptr;
};

/**
 * @brief Bank of N baked palettes for bulk Persist/clone operations.
 */
struct BakedPaletteBank {
  static constexpr int N = 5;
  BakedPalette entries[N];

  /**
   * @brief Deep-copies all entries into a target arena.
   * @param src Source bank to copy from.
   * @param dst Destination bank to fill.
   * @param arena Arena to allocate the cloned LUTs from.
   * @details Required by Cloneable.
   */
  static void clone(const BakedPaletteBank &src, BakedPaletteBank &dst,
                    Arena &arena) {
    for (int i = 0; i < N; ++i)
      dst.entries[i].clone_from(src.entries[i], arena);
  }
};
