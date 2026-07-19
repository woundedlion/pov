/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

#include "engine/platform.h"
#include "math/3dmath.h"
#include "engine/util.h"

#include "engine/memory.h"

#include "color/gamut_lut.h"

#if defined(__ARM_FEATURE_DSP)
// Inline assembly avoids a CMSIS header dependency for the saturating add.
__attribute__((always_inline)) static inline uint32_t inline_uqadd16(uint32_t a, uint32_t b) {
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
  if (lo > 0xFFFFu) lo = 0xFFFFu;
  if (hi > 0xFFFFu) hi = 0xFFFFu;
  return (hi << 16) | lo;
}
#endif

struct Pixel16;
// Saturating per-channel add packed into two uqadd16 lanes (g|b in one 32-bit
// word, r alone in another). Used by Pixel16::operator+=.
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
   * @details Explicit so a stray Pixel16 in a CRGB context is a compile error,
   * not a silent round-trip through 8-bit gamma.
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
   * @brief Saturated per-channel sum of two pixels.
   * @param rhs Pixel to add.
   * @return A new pixel with each channel clamped to the 16-bit max.
   */
  Pixel16 operator+(const Pixel16 &rhs) const {
    Pixel16 out = *this;
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
  Pixel16 operator*(float s) const {
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
   *   (x + (x>>16) + 32768) >> 16, within 1 LSB of round(x/65535) and exact at
   * the endpoints (frac 0/65535 -> a/b). Plain 32-bit MACs, not packed `smlad`:
   * smlad's signed 16x16 dual-MAC reads an operand >= 32768 as negative.
   */
  __attribute__((always_inline)) Pixel16 lerp16(const Pixel16 &other, uint16_t frac) const {
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
 * @brief Quantizes a [0,1] interpolation fraction to a 16-bit lerp16 weight.
 * @param frac Blend fraction; clamped to [0, 1].
 * @return The fraction as a 16-bit weight in [0, 65535], rounded.
 */
__attribute__((always_inline)) inline uint16_t frac_to_q16(float frac) {
  return static_cast<uint16_t>(hs::clamp(frac, 0.0f, 1.0f) * 65535.0f + 0.5f);
}

/**
 * @brief Represents a color with a STRAIGHT (non-premultiplied) alpha channel.
 * @details `color` holds the un-premultiplied color, `alpha` its coverage. The
 * arithmetic operators (`operator*=`, `operator+=`) are NOT alpha compositing:
 * they treat (color, alpha) as a plain 4-vector for SSAA averaging (scale by
 * 1/N, then sum), so `*=` scales alpha alongside color. Premultiplication
 * happens once, at the final canvas write (`color * alpha`).
 */
struct Color4 {
  Pixel color;
  float alpha;

  /**
   * @brief Constructs a transparent black color (alpha 0.0), the additive/SSAA
   * accumulator identity for `operator+=`.
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
   * @details t clamped to [0,1] so out-of-range t saturates at an endpoint
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
   * @brief Adds another color's pixel and alpha into this one (both saturating).
   * @param rhs Color to add.
   * @return Reference to this color after the add.
   * @details Pixel add saturates at the channel max; alpha saturates at 1.0 so
   * the sum stays a valid blend weight.
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
   * @details Explicit so a Color4 never silently round-trips through 8-bit gamma.
   */
  explicit operator CRGB() const { return static_cast<CRGB>(color); }
};

/**
 * @brief Sums two colors, pixel and alpha both saturating (see operator+=).
 * @param lhs First color (taken by value as the accumulator).
 * @param rhs Color to add.
 * @return The saturated sum.
 */
inline Color4 operator+(Color4 lhs, const Color4 &rhs) { return lhs += rhs; }

/**
 * @brief Scales a color's pixel and alpha by a float (see operator*=).
 * @param lhs Color to scale (taken by value).
 * @param s Scale factor.
 * @return The scaled color.
 */
inline Color4 operator*(Color4 lhs, float s) { return lhs *= s; }

/**
 * @brief Scales a color's pixel and alpha by a float (see operator*=).
 * @param s Scale factor.
 * @param rhs Color to scale.
 * @return The scaled color.
 */
inline Color4 operator*(float s, const Color4 &rhs) { return rhs * s; }

/**
 * @brief Perceptual (OKLab) hue rotation with a precomputed rotation.
 * @param c Source color.
 * @param ca Cosine of the rotation angle.
 * @param sa Sine of the rotation angle.
 * @return The hue-rotated color.
 * @details Precomputed (ca, sa) lets frame-constant callers hoist sin/cos out
 * of the per-pixel loop.
 */
inline Color4 hue_rotate(const Color4 &c, float ca, float sa);
/**
 * @brief Perceptual (OKLab) hue rotation by a turn amount.
 * @param c Source color.
 * @param amount Rotation in turns (0..1 = full turn).
 * @return The hue-rotated color.
 */
inline Color4 hue_rotate(const Color4 &c, float amount);

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

// Packs g|b into the low add lane and r into a separate lane (its high halfword
// stays 0, so uqadd16's upper add is a harmless 0+0).
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
 * @details Rounds to nearest (+0.5f) and clamps in float before the cast:
 * casting an unclamped a*65535 to int is UB when a is NaN or overflows int. The
 * +0.5f sits inside the clamp so a == 1 still maps exactly to 65535.
 */
inline auto blend_alpha(float a) {
  uint16_t ai = (uint16_t)hs::clamp(a * 65535.0f + 0.5f, 0.0f, 65535.0f);
  return [ai](const Pixel &c1, const Pixel &c2) { return c1.lerp16(c2, ai); };
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
 * @details Not constexpr: powf is not a constant expression. inline for ODR.
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

/** @brief Linear-RGB triple in [0,1] (may sit slightly out of gamut before
 *  clamping). */
struct LinRGB { float r, g, b; };

/**
 * @brief Linear-RGB -> LMS cone response (the first OKLab matrix).
 * @param r Linear red in [0, 1].
 * @param g Linear green in [0, 1].
 * @param b Linear blue in [0, 1].
 * @return The (l, m, s) cone responses, before the cube-root nonlinearity.
 * @details Shared by linear_rgb_to_oklab and hue_rotate; each applies its own
 * cube-root (exact cbrtf vs fast_cbrt) then calls lms_to_oklab.
 */
HS_O3_FN inline LMS linear_rgb_to_lms(float r, float g, float b) {
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
 * @details Takes the already-cube-rooted triple so the caller picks the
 * cube-root flavour.
 */
HS_O3_FN inline OKLab lms_to_oklab(float l_, float m_, float s_) {
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
 * @brief Converts OKLab to cube-rooted LMS (the inverse OKLab matrix).
 * @param lab Source color in OKLab space.
 * @param l_ Out: cube-rooted l cone response.
 * @param m_ Out: cube-rooted m cone response.
 * @param s_ Out: cube-rooted s cone response.
 */
HS_O3_FN inline void oklab_to_lms_cbrt(OKLab lab, float &l_, float &m_, float &s_) {
  l_ = lab.L + 0.3963377774f * lab.a + 0.2158037573f * lab.b;
  m_ = lab.L - 0.1055613458f * lab.a - 0.0638541728f * lab.b;
  s_ = lab.L - 0.0894841775f * lab.a - 1.2914855480f * lab.b;
}

/**
 * @brief Converts cube-rooted LMS to linear RGB [0,1] (cube + RGB matrix).
 * @param l_ Cube-rooted l cone response.
 * @param m_ Cube-rooted m cone response.
 * @param s_ Cube-rooted s cone response.
 * @param r Out: linear red (may exit gamut before clamping).
 * @param g Out: linear green.
 * @param b Out: linear blue.
 */
HS_O3_FN
inline void lms_cbrt_to_linear_rgb(float l_, float m_, float s_, float &r,
                                   float &g, float &b) {
  float l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_;

  r = +4.0767416621f * l - 3.3077115913f * m + 0.2309699292f * s;
  g = -1.2684380046f * l + 2.6097574011f * m - 0.3413193965f * s;
  b = -0.0041960863f * l - 0.7034186147f * m + 1.7076147010f * s;
}

/**
 * @brief Converts OKLab to linear RGB [0,1].
 * @param lab Source color in OKLab space.
 * @param r Out: linear red in [0, 1] (may exit gamut before clamping).
 * @param g Out: linear green in [0, 1].
 * @param b Out: linear blue in [0, 1].
 */
HS_O3_FN inline void oklab_to_linear_rgb(OKLab lab, float &r, float &g, float &b) {
  float l_, m_, s_;
  oklab_to_lms_cbrt(lab, l_, m_, s_);
  lms_cbrt_to_linear_rgb(l_, m_, s_, r, g, b);
}

/**
 * @brief Converts OKLab to linear RGB [0,1].
 * @param lab Source color in OKLab space.
 * @return The linear-RGB triple (may exit gamut before clamping).
 */
HS_O3_FN inline LinRGB oklab_to_linear_rgb(OKLab lab) {
  LinRGB rgb;
  oklab_to_linear_rgb(lab, rgb.r, rgb.g, rgb.b);
  return rgb;
}

/**
 * @brief Tests whether a linear-RGB triple lies inside the [0,1] display cube.
 * @param r Linear red.
 * @param g Linear green.
 * @param b Linear blue.
 * @return True if every channel is within [0,1] (with a small epsilon slack).
 * @details The epsilon slack absorbs float rounding that can leave an in-gamut
 * color a hair past 1.0 after the OKLab inverse.
 */
HS_O3_FN inline bool linear_rgb_in_gamut(float r, float g, float b) {
  constexpr float lo = -1e-4f, hi = 1.0f + 1e-4f;
  const float least = __builtin_fminf(__builtin_fminf(r, g), b);
  const float most = __builtin_fmaxf(__builtin_fmaxf(r, g), b);
  return least >= lo && most <= hi;
}

// The tolerance linear_rgb_in_gamut allows. The boundary cubics are solved
// against it so a solved crossing is the same crossing the gate reports.
inline constexpr float GAMUT_SLACK = 1e-4f;

// Chroma pulled back off the solved crossing. Without it the caller's own
// re-conversion rounds a channel a part in a million past the gate; with it the
// worst residual over a dense sweep of L, hue and chroma is comfortably inside.
inline constexpr float GAMUT_CLIP_MARGIN = 2e-5f;

// Newton steps per crossing. Eight is what the slowest crossings need: at the
// blue cube vertex a face runs nearly tangent to the ray, and a near-double
// root drops Newton from quadratic to one bit per step.
inline constexpr int GAMUT_REFINE_STEPS = 8;

/**
 * @brief Real roots of a2 t^2 + a1 t + a0, ordered.
 * @param a2 Quadratic coefficient.
 * @param a1 Linear coefficient.
 * @param a0 Constant coefficient.
 * @param r0 Out: smaller root.
 * @param r1 Out: larger root.
 * @return False when the roots are complex or the quadratic is degenerate.
 */
HS_O3_FN inline bool quadratic_roots(float a2, float a1, float a0, float &r0,
                                     float &r1) {
  if (!(std::fabs(a2) > 1e-20f)) {
    if (!(std::fabs(a1) > 1e-20f))
      return false;
    r0 = r1 = -a0 / a1;
    return true;
  }
  const float disc = a1 * a1 - 4.0f * a2 * a0;
  if (disc < 0.0f)
    return false;
  const float sd = std::sqrt(disc);
  const float inv = 0.5f / a2;
  const float p = (-a1 - sd) * inv, q = (-a1 + sd) * inv;
  r0 = p < q ? p : q;
  r1 = p < q ? q : p;
  return true;
}

/**
 * @brief Chroma at which a channel cubic crosses zero, inside one segment.
 * @param a3 Cubic coefficient.
 * @param a2 Quadratic coefficient.
 * @param a1 Linear coefficient.
 * @param a0 Constant coefficient, already offset by the face being crossed.
 * @param x Segment start, on the in-range side of the crossing.
 * @param y Segment end, past the crossing.
 * @return The crossing chroma, inside [x, y].
 * @details The segment is monotone and of one convexity, so Newton started from
 * the end the cubic curves toward closes on the root from that side. The
 * bracket narrows on every step and a step leaving it becomes a bisection, so a
 * flat slope near a turning point cannot throw the result.
 */
HS_O3_FN inline float cubic_segment_root(float a3, float a2, float a1, float a0,
                                         float x, float y) {
  const float gx = ((a3 * x + a2) * x + a1) * x + a0;
  const float curv = 3.0f * a3 * (x + y) + 2.0f * a2;
  float t = gx * curv > 0.0f ? x : y;

  float lo = x, hi = y;
  for (int i = 0; i < GAMUT_REFINE_STEPS; ++i) {
    const float g = ((a3 * t + a2) * t + a1) * t + a0;
    if (g * gx > 0.0f)
      lo = t;
    else
      hi = t;
    const float dg = (3.0f * a3 * t + 2.0f * a2) * t + a1;
    const float next = std::fabs(dg) > 1e-20f ? t - g / dg : hi;
    // Inclusive bounds: a converged step lands exactly on the bound it just
    // set, and rejecting it would bisect a solved root back apart.
    t = next >= lo && next <= hi ? next : 0.5f * (lo + hi);
  }
  // The iterates close on the root from the far side, so the last residual is
  // subtracted rather than applied: the answer lands just inside the crossing.
  const float g = ((a3 * t + a2) * t + a1) * t + a0;
  const float dg = (3.0f * a3 * t + 2.0f * a2) * t + a1;
  if (!(std::fabs(dg) > 1e-20f))
    return lo;
  // Twice the residual: at a near-tangent crossing the two roots sit close
  // together and a single Newton distance under-reads the gap by half. Where
  // the crossing is clean the doubled pull-back costs well under a part in a
  // million of chroma.
  const float back = 2.0f * g / dg;
  const float pulled = t - (back < 0.0f ? -back : back);
  return pulled > lo ? pulled : lo;
}

/**
 * @brief Smallest chroma at which one RGB channel leaves the display range.
 * @param wl Channel's l weight in the LMS -> linear-RGB matrix.
 * @param wm Channel's m weight.
 * @param ws Channel's s weight.
 * @param L Lightness held fixed along the ray.
 * @param L3 L cubed, the channel's value at zero chroma.
 * @param k_l Rate of change of cube-rooted l with chroma.
 * @param k_m Rate of change of cube-rooted m with chroma.
 * @param k_s Rate of change of cube-rooted s with chroma.
 * @param best Smallest exit chroma found so far; the search stops there.
 * @return `best`, replaced by a smaller exit chroma when one exists.
 * @details The channel is a cubic in chroma. Its two turning points and its
 * inflection -- which always lies between them -- cut the ray into monotone,
 * singly curved segments, so the first segment whose far end is out of range
 * holds the FIRST crossing, whichever face it goes through. Walking the
 * segments rather than refining a seed is what makes the answer the first exit:
 * near the blue cube vertex a channel dips through a face and back, and a
 * seeded refinement settles on the far side of that dip. Kept out of line:
 * three inlined copies per call site overrun the Teensy ITCM budget.
 */
HS_O3_FN __attribute__((noinline)) inline float
gamut_channel_exit(float wl, float wm, float ws, float L, float L3, float k_l,
                   float k_m, float k_s, float best) {
  const float a3 =
      wl * k_l * k_l * k_l + wm * k_m * k_m * k_m + ws * k_s * k_s * k_s;
  const float a2 =
      3.0f * L * (wl * k_l * k_l + wm * k_m * k_m + ws * k_s * k_s);
  const float a1 = 3.0f * L * L * (wl * k_l + wm * k_m + ws * k_s);

  float c0, c1;
  const bool turns = quadratic_roots(3.0f * a3, 2.0f * a2, a1, c0, c1);
  const float infl = std::fabs(a3) > 1e-20f ? -a2 * (1.0f / 3.0f) / a3 : best;

  // Segment ends, ascending and capped at the incumbent: a crossing past it
  // cannot win. The inflection always sits between the turning points.
  float ends[4];
  int n = 0;
  if (turns && c0 > 0.0f && c0 < best)
    ends[n++] = c0;
  if (infl > 0.0f && infl < best && (n == 0 || infl > ends[n - 1]))
    ends[n++] = infl;
  if (turns && c1 > 0.0f && c1 < best && (n == 0 || c1 > ends[n - 1]))
    ends[n++] = c1;
  ends[n++] = best;

  // L3 sits inside the range for every L in [0,1], so the walk starts in range.
  float x = 0.0f;
  for (int i = 0; i < n; ++i) {
    const float y = ends[i];
    const float fy = ((a3 * y + a2) * y + a1) * y + L3;
    if (fy < -GAMUT_SLACK)
      return cubic_segment_root(a3, a2, a1, L3 + GAMUT_SLACK, x, y);
    if (fy > 1.0f + GAMUT_SLACK)
      return cubic_segment_root(a3, a2, a1, L3 - 1.0f - GAMUT_SLACK, x, y);
    x = y;
  }
  return best;
}

/**
 * @brief Maps an out-of-gamut OKLab color into the sRGB cube by solving the
 * channel cubics, holding hue and lightness fixed.
 * @param lab Source color; L is clamped to [0,1] internally.
 * @return An OKLab color with the same L and hue but chroma scaled down until it
 * is just inside the display cube.
 * @details At fixed L and hue each RGB channel is a cubic in chroma, so the
 * boundary is the smallest chroma at which any channel crosses a face of the
 * cube. Each channel is searched over its whole ray rather than around a seed,
 * which is what makes the answer the FIRST exit: near the blue cube vertex a
 * channel dips through a face and back, and a seeded refinement settles on the
 * far side of that dip. This is the fallback path, taken when no owner has
 * armed the boundary table.
 */
HS_O3_FN inline OKLab gamut_clip_analytic(OKLab lab) {
  // The C=0 achromatic floor is in gamut only for L in [0,1]; an out-of-range L
  // would leave even the floor out of gamut and return a non-gamut color.
  lab.L = hs::clamp(lab.L, 0.0f, 1.0f);
  const float C = std::sqrt(lab.a * lab.a + lab.b * lab.b);
  if (!(C > 1e-6f))
    return {lab.L, 0.0f, 0.0f};
  const float inv_C = 1.0f / C;
  const float a_ = lab.a * inv_C, b_ = lab.b * inv_C;

  const float k_l = 0.3963377774f * a_ + 0.2158037573f * b_;
  const float k_m = -0.1055613458f * a_ - 0.0638541728f * b_;
  const float k_s = -0.0894841775f * a_ - 1.2914855480f * b_;
  const float L3 = lab.L * lab.L * lab.L;

  // Seeded with the input chroma: a ray that never exits below it was already
  // in gamut, and the seed also caps the segment walk at a finite chroma.
  float t = C;
  t = gamut_channel_exit(4.0767416621f, -3.3077115913f, 0.2309699292f, lab.L,
                         L3, k_l, k_m, k_s, t);
  t = gamut_channel_exit(-1.2684380046f, 2.6097574011f, -0.3413193965f, lab.L,
                         L3, k_l, k_m, k_s, t);
  t = gamut_channel_exit(-0.0041960863f, -0.7034186147f, 1.7076147010f, lab.L,
                         L3, k_l, k_m, k_s, t);

  t = hs::clamp(t - GAMUT_CLIP_MARGIN, 0.0f, C);
  return {lab.L, a_ * t, b_ * t};
}

/**
 * @brief Arena-resident gamut boundary bracket grid and the scales indexing it.
 * @details `table` null means no owner has armed the grid and the clip path
 * solves the cubics instead. Grouped in one object so the per-pixel path loads
 * a single base address rather than five unrelated globals.
 */
struct GamutLut {
  const uint16_t *table = nullptr; /**< Arena copy, or null when disarmed. */
  int angle_steps = 0;             /**< Diamond-angle buckets over [0, 4). */
  int l_steps = 0;                 /**< Lightness buckets over [0, 1]. */
  float angle_scale = 0.0f;        /**< angle_steps / 4, the index scale. */
  float l_scale = 0.0f;            /**< l_steps, the index scale. */
};

/**
 * @brief The single live boundary grid.
 * @details Arena-resident so the scattered per-pixel reads land in RAM rather
 * than the QSPI flash the master table lives in. The persistent arena is reset
 * between effects, so the owning effect must disarm this in its destructor;
 * Canvas guarantees the outgoing effect is destroyed before the next one is
 * constructed and initialized.
 */
inline GamutLut g_gamut_lut;

/**
 * @brief Downsamples GAMUT_LUT into @p arena and arms the table-driven clip.
 * @param arena Arena to hold the copy; must outlive every clip call that
 *        follows, and the caller must disarm before it is reset.
 * @param angle_steps Diamond-angle buckets; must divide GAMUT_LUT_ANGLE_STEPS.
 * @param l_steps Lightness buckets; must divide GAMUT_LUT_L_STEPS.
 * @details Call after the arenas are configured, from the owning effect's
 * init(). A coarse cell takes the minimum of the merged minima and the maximum
 * of the merged maxima, so the true boundary of every ray in the cell still
 * lies inside the stored bracket at any resolution. Cost in arena bytes is
 * gamut_lut_bytes(angle_steps, l_steps); resolution only sets how wide the
 * bracket starts, and the per-pixel bisection sets how far it is narrowed.
 */
inline void init_gamut_lut(Arena &arena, int angle_steps, int l_steps) {
  HS_CHECK(angle_steps > 0 && l_steps > 0 &&
               GAMUT_LUT_ANGLE_STEPS % angle_steps == 0 &&
               GAMUT_LUT_L_STEPS % l_steps == 0,
           "init_gamut_lut: %d x %d must divide the %d x %d flash master",
           angle_steps, l_steps, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);

  const int sa = GAMUT_LUT_ANGLE_STEPS / angle_steps;
  const int sl = GAMUT_LUT_L_STEPS / l_steps;
  uint16_t *dst = arena.allocate_n<uint16_t>(angle_steps * l_steps * 2);

  for (int l = 0; l < l_steps; ++l) {
    for (int a = 0; a < angle_steps; ++a) {
      uint16_t c_lo = 0xFFFF, c_hi = 0;
      for (int dl = 0; dl < sl; ++dl) {
        const uint16_t *row =
            &GAMUT_LUT[((l * sl + dl) * GAMUT_LUT_ANGLE_STEPS + a * sa) * 2];
        for (int da = 0; da < sa; ++da) {
          c_lo = std::min(c_lo, row[da * 2]);
          c_hi = std::max(c_hi, row[da * 2 + 1]);
        }
      }
      dst[(l * angle_steps + a) * 2] = c_lo;
      dst[(l * angle_steps + a) * 2 + 1] = c_hi;
    }
  }

  g_gamut_lut = {dst, angle_steps, l_steps, angle_steps * 0.25f,
                 static_cast<float>(l_steps)};
}

/**
 * @brief Disarms the table-driven clip path, reverting to the cubic solve.
 * @details Must run before the arena holding the copy is reset or reused.
 */
inline void release_gamut_lut() { g_gamut_lut = GamutLut{}; }

// Equal steps the stored bracket is walked in, looking for the first one that
// leaves the gamut. A walk rather than a straight bisection because the gate's
// tolerance lets the in-gamut set along a ray break into pieces: bisecting a
// bracket that spans a gap converges on the far side of it, which is in gamut
// but past the first exit and discontinuous in L against the neighbouring cell.
inline constexpr int GAMUT_SCAN_STEPS = 4;

// Bisections inside the walk step that straddles the crossing. Residual is the
// bracket width over GAMUT_SCAN_STEPS, halved once per step, so this is an
// accuracy knob and not a cap: three take the 256 x 128 grid's worst
// mid-lightness bracket to 0.0016 chroma.
inline constexpr int GAMUT_BRACKET_STEPS = 3;

/**
 * @brief Largest in-gamut scale of (a, b) inside a bracketed scale range.
 * @param L Lightness held fixed along the ray.
 * @param a OKLab a of the input, unnormalized; the ray is u * (a, b).
 * @param b OKLab b of the input, unnormalized.
 * @param lo Scale known to be at or below the boundary.
 * @param hi Scale at or above it, already capped at the input's own scale.
 * @return The refined scale, in gamut by construction.
 * @details With l_ = L + A*u and (L + X*u)^3 expanded in u, every linear-RGB
 * channel is a cubic in u whose four coefficients depend only on L and the hue
 * direction and are built once here. A refinement step is then three Horner
 * evaluations and the six bound tests, not an OKLab round trip. The bound tests
 * use linear_rgb_in_gamut's own tolerance, so a solved crossing is the crossing
 * the gate reports. `lo` is probed before it is trusted: a cell minimum that
 * over-reads its region drops the search back to zero chroma rather than
 * returning a color outside the cube.
 */
HS_O3_FN __attribute__((noinline)) inline float
gamut_bracket_refine(float L, float a, float b, float lo, float hi) {
  const float ka = 0.3963377774f * a + 0.2158037573f * b;
  const float km = -0.1055613458f * a - 0.0638541728f * b;
  const float ks = -0.0894841775f * a - 1.2914855480f * b;

  const float ka2 = ka * ka, ka3 = ka2 * ka;
  const float km2 = km * km, km3 = km2 * km;
  const float ks2 = ks * ks, ks3 = ks2 * ks;
  const float l3 = L * L * L, q = 3.0f * L * L, c = 3.0f * L;

  const float r3 = 4.0767416621f * ka3 - 3.3077115913f * km3 +
                   0.2309699292f * ks3;
  const float r2 =
      c * (4.0767416621f * ka2 - 3.3077115913f * km2 + 0.2309699292f * ks2);
  const float r1 =
      q * (4.0767416621f * ka - 3.3077115913f * km + 0.2309699292f * ks);
  const float g3 = -1.2684380046f * ka3 + 2.6097574011f * km3 -
                   0.3413193965f * ks3;
  const float g2 =
      c * (-1.2684380046f * ka2 + 2.6097574011f * km2 - 0.3413193965f * ks2);
  const float g1 =
      q * (-1.2684380046f * ka + 2.6097574011f * km - 0.3413193965f * ks);
  const float b3 = -0.0041960863f * ka3 - 0.7034186147f * km3 +
                   1.7076147010f * ks3;
  const float b2 =
      c * (-0.0041960863f * ka2 - 0.7034186147f * km2 + 1.7076147010f * ks2);
  const float b1 =
      q * (-0.0041960863f * ka - 0.7034186147f * km + 1.7076147010f * ks);

  const auto inside = [&](float u) {
    const float rv = ((r3 * u + r2) * u + r1) * u + l3;
    const float gv = ((g3 * u + g2) * u + g1) * u + l3;
    const float bv = ((b3 * u + b2) * u + b1) * u + l3;
    return linear_rgb_in_gamut(rv, gv, bv);
  };

  float x = lo, y = hi;
  if (inside(lo)) {
    const float step = (hi - lo) * (1.0f / GAMUT_SCAN_STEPS);
    int i = 0;
    for (; i < GAMUT_SCAN_STEPS; ++i) {
      y = (i + 1 == GAMUT_SCAN_STEPS) ? hi : x + step;
      if (!inside(y))
        break;
      x = y;
    }
    // The whole bracket held, so the input's own chroma was already in gamut.
    if (i == GAMUT_SCAN_STEPS)
      return hi;
  } else {
    // The cell minimum over-reads its region; fall back to the whole ray.
    x = 0.0f;
    y = lo;
  }

  for (int i = 0; i < GAMUT_BRACKET_STEPS; ++i) {
    const float mid = 0.5f * (x + y);
    if (inside(mid))
      x = mid;
    else
      y = mid;
  }
  return x;
}

/**
 * @brief Maps an out-of-gamut OKLab color into the sRGB cube by reducing chroma,
 * holding hue and lightness fixed (Ottosson's preserve-chroma projection).
 * @param lab Source color, assumed out of gamut (callers gate on
 * linear_rgb_in_gamut); L is clamped to [0,1] internally.
 * @return An OKLab color with the same L and hue but chroma scaled down until it
 * is just inside the display cube.
 * @details The answer is exactly min(C, C_max(hue, L)), and C_max is a static
 * property of the sRGB gamut. With the grid armed, one cell read brackets C_max
 * and gamut_bracket_refine narrows the bracket; a chroma at or below the cell
 * minimum needs neither. Without it the cubics are solved outright.
 */
HS_O3_FN __attribute__((noinline)) inline OKLab
gamut_clip_preserve_chroma(OKLab lab) {
  const GamutLut &lut = g_gamut_lut;
  if (!lut.table)
    return gamut_clip_analytic(lab);

  // The u=0 achromatic floor is in gamut only for L in [0,1]; an out-of-range L
  // would leave even the floor out of gamut and return a non-gamut color.
  lab.L = hs::clamp(lab.L, 0.0f, 1.0f);
  const float c_sq = lab.a * lab.a + lab.b * lab.b;
  // Below this the color is achromatic past 16-bit output precision, and the
  // reciprocal square root would overflow.
  if (!(c_sq > 1e-12f))
    return {lab.L, 0.0f, 0.0f};

  int ai = static_cast<int>(diamond_angle(lab.b, lab.a) * lut.angle_scale);
  int li = static_cast<int>(lab.L * lut.l_scale);
  ai = hs::clamp(ai, 0, lut.angle_steps - 1);
  li = hs::clamp(li, 0, lut.l_steps - 1);
  const uint16_t *cell = &lut.table[(li * lut.angle_steps + ai) * 2];

  const float c_lo = static_cast<float>(cell[0]) * GAMUT_LUT_INV_SCALE;
  if (c_sq <= c_lo * c_lo)
    return lab;

  const float c_hi = static_cast<float>(cell[1]) * GAMUT_LUT_INV_SCALE;
  const float inv_c = fast_rsqrt(c_sq);
  const float hi = std::min(1.0f, c_hi * inv_c);
  float u = gamut_bracket_refine(lab.L, lab.a, lab.b, c_lo * inv_c, hi);
  // Pulled back off the crossing: without it the caller's own re-conversion
  // rounds a channel a part in a million past the gate.
  u -= GAMUT_CLIP_MARGIN * inv_c;
  if (u < 0.0f)
    u = 0.0f;
  return {lab.L, lab.a * u, lab.b * u};
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
 * does it pay for gamut_clip_preserve_chroma and re-convert. In-gamut colors
 * cost one matrix mul plus the gate test, no search.
 */
inline void oklab_to_linear_rgb_gamut(OKLab lab, float &r, float &g, float &b) {
  oklab_to_linear_rgb(lab, r, g, b);
  if (!linear_rgb_in_gamut(r, g, b))
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r, g, b);
}

/**
 * @brief OKLab -> linear RGB with chroma-reduction gamut mapping off the fast
 * path.
 * @param lab Source color in OKLab space.
 * @return The gamut-mapped linear-RGB triple (may sit a hair past the bound;
 * callers still clamp).
 */
inline LinRGB oklab_to_linear_rgb_gamut(OKLab lab) {
  LinRGB rgb;
  oklab_to_linear_rgb_gamut(lab, rgb.r, rgb.g, rgb.b);
  return rgb;
}

/**
 * @brief Builds the cbrt-LMS-space 3x3 equivalent to an OKLab chroma rotation.
 * @param ca Cosine of the rotation angle.
 * @param sa Sine of the rotation angle.
 * @param k Out: row-major 3x3 acting on cube-rooted LMS.
 * @details The OKLab a/b rotation is linear in cube-rooted LMS, so
 * oklab_to_lms_cbrt . rotate . lms_to_oklab folds into one matrix. Columns are
 * derived by pushing basis vectors through those functions.
 */
inline void hue_rotate_lms_matrix(float ca, float sa, float k[9]) {
  for (int i = 0; i < 3; ++i) {
    OKLab lab = lms_to_oklab(i == 0 ? 1.0f : 0.0f, i == 1 ? 1.0f : 0.0f,
                             i == 2 ? 1.0f : 0.0f);
    float a2 = lab.a * ca - lab.b * sa;
    float b2 = lab.a * sa + lab.b * ca;
    float l_, m_, s_;
    oklab_to_lms_cbrt({lab.L, a2, b2}, l_, m_, s_);
    k[i] = l_;
    k[3 + i] = m_;
    k[6 + i] = s_;
  }
}

/**
 * @brief Applies a cbrt-LMS 3x3 (from hue_rotate_lms_matrix, optionally
 * uniformly scaled) and converts to linear RGB with gamut mapping.
 * @param k Row-major 3x3 acting on cube-rooted LMS.
 * @param l_ Cube-rooted l cone response.
 * @param m_ Cube-rooted m cone response.
 * @param s_ Cube-rooted s cone response.
 * @param r Out: gamut-mapped linear red (may sit a hair past the bound;
 * callers still clamp).
 * @param g Out: linear green.
 * @param b Out: linear blue.
 * @details In-gamut colors never leave cbrt-LMS; the OKLab form is recomputed
 * only on the chroma-clip slow path.
 */
HS_O3_FN
inline void lms_cbrt_transform_rgb(const float k[9], float l_, float m_,
                                   float s_, float &r, float &g, float &b) {
  float ul = k[0] * l_ + k[1] * m_ + k[2] * s_;
  float um = k[3] * l_ + k[4] * m_ + k[5] * s_;
  float us = k[6] * l_ + k[7] * m_ + k[8] * s_;
  lms_cbrt_to_linear_rgb(ul, um, us, r, g, b);
  if (!linear_rgb_in_gamut(r, g, b)) {
    HS_PROFILE_DEEP(gamut_clip);
    OKLab lab = lms_to_oklab(ul, um, us);
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r, g, b);
  }
}

/**
 * @brief Quantizes a [0,1] linear channel to a 16-bit Pixel component.
 * @param v Linear channel value; clamped to [0, 1].
 * @return The channel as a 16-bit value in [0, 65535].
 * @details Clamps, then rounds (+0.5f) rather than truncating; truncation
 * would bias every channel down by up to ~1/65535.
 */
HS_O3_FN inline uint16_t float_to_pixel16(float v) {
  return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * 65535.0f + 0.5f);
}

/**
 * @brief Quantizes a [0,1] linear channel to an 8-bit sRGB component.
 * @param l Linear channel value; clamped to [0, 1].
 * @return The channel as an 8-bit sRGB value in [0, 255].
 */
inline uint8_t linear_float_to_srgb8(float l) {
  return static_cast<uint8_t>(hs::clamp(
      linear_to_srgb_float(hs::clamp(l, 0.0f, 1.0f)) * 255.0f + 0.5f, 0.0f,
      255.0f));
}

/**
 * @brief Rotates the (a,b) chroma plane in OKLab on a linear-RGB float triple.
 * @param r In/out: linear red in [0, 1].
 * @param g In/out: linear green.
 * @param b In/out: linear blue.
 * @param ca Cosine of the rotation angle.
 * @param sa Sine of the rotation angle.
 * @details fast_cbrt forward, exact cubes inverse, direct 2D rotation of (a,b)
 * (no atan2/sqrt OKLCH polar round-trip). Preserves lightness to fast_cbrt
 * accuracy, chroma to fast-trig accuracy.
 */
HS_O3_FN inline void hue_rotate_rgb(float &r, float &g, float &b, float ca, float sa) {
  LMS lms = linear_rgb_to_lms(r, g, b);
  OKLab lab = lms_to_oklab(fast_cbrt(lms.l), fast_cbrt(lms.m), fast_cbrt(lms.s));

  float a2 = lab.a * ca - lab.b * sa;
  float b2 = lab.a * sa + lab.b * ca;

  oklab_to_linear_rgb_gamut({lab.L, a2, b2}, r, g, b);
}

inline Color4 hue_rotate(const Color4 &c, float ca, float sa) {
  constexpr float INV16 = 1.0f / 65535.0f;
  float r = c.color.r * INV16, g = c.color.g * INV16, b = c.color.b * INV16;

  hue_rotate_rgb(r, g, b, ca, sa);

  Color4 result = c;
  result.color.r = float_to_pixel16(r);
  result.color.g = float_to_pixel16(g);
  result.color.b = float_to_pixel16(b);
  return result;
}

inline Color4 hue_rotate(const Color4 &c, float amount) {
  float angle = amount * (2.0f * PI_F);
  float ca = fast_cosf(angle);
  float sa = fast_sinf(angle);
  // renormalize fast trig so the rotation preserves chroma (else the scaling
  // compounds per frame under feedback).
  float inv = 1.0f / sqrtf(ca * ca + sa * sa);
  return hue_rotate(c, ca * inv, sa * inv);
}

/**
 * @brief A Color4 pre-converted to OKLab for repeated hue rotations.
 * @details Caches the forward linear-RGB -> cbrt-LMS -> OKLab transform of a
 * fixed base color so each rotation pays only the (a,b) rotation and the
 * inverse transform; hue_rotate(base, amount) matches hue_rotate(c, amount)
 * exactly.
 */
struct HueRotateBase {
  OKLab lab;   /**< Base color in OKLab. */
  Color4 base; /**< Original color; alpha is carried into each result. */
};

/**
 * @brief Precomputes the OKLab form of a color for repeated hue rotations.
 * @param c Base color.
 * @return The precomputed base for hue_rotate(base, amount).
 */
inline HueRotateBase make_hue_rotate_base(const Color4 &c) {
  constexpr float INV16 = 1.0f / 65535.0f;
  float r = c.color.r * INV16, g = c.color.g * INV16, b = c.color.b * INV16;
  LMS lms = linear_rgb_to_lms(r, g, b);
  return {lms_to_oklab(fast_cbrt(lms.l), fast_cbrt(lms.m), fast_cbrt(lms.s)),
          c};
}

/**
 * @brief Perceptual (OKLab) hue rotation of a precomputed base color.
 * @param hb Precomputed base from make_hue_rotate_base().
 * @param amount Rotation in turns (0..1 = full turn).
 * @return The hue-rotated color.
 */
inline Color4 hue_rotate(const HueRotateBase &hb, float amount) {
  float angle = amount * (2.0f * PI_F);
  float ca = fast_cosf(angle);
  float sa = fast_sinf(angle);
  // renormalize fast trig so the rotation preserves chroma (else the scaling
  // compounds per frame under feedback).
  float inv = 1.0f / sqrtf(ca * ca + sa * sa);
  ca *= inv;
  sa *= inv;

  float a2 = hb.lab.a * ca - hb.lab.b * sa;
  float b2 = hb.lab.a * sa + hb.lab.b * ca;
  float r, g, b;
  oklab_to_linear_rgb_gamut({hb.lab.L, a2, b2}, r, g, b);

  Color4 result = hb.base;
  result.color.r = float_to_pixel16(r);
  result.color.g = float_to_pixel16(g);
  result.color.b = float_to_pixel16(b);
  return result;
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
  constexpr float INV16 = 1.0f / 65535.0f;
  float rf = srgb_to_linear_lut[r] * INV16;
  float gf = srgb_to_linear_lut[g] * INV16;
  float bf = srgb_to_linear_lut[b] * INV16;
  return oklab_to_oklch(linear_rgb_to_oklab(rf, gf, bf));
}

/**
 * @brief Converts OKLCH to a 16-bit linear Pixel (gamut-clamped).
 * @param lch Source color in OKLCH space.
 * @return The color as a linear-space Pixel.
 */
inline Pixel oklch_to_pixel(OKLCH lch) {
  LinRGB rgb = oklab_to_linear_rgb_gamut(oklch_to_oklab(lch));
  return Pixel(float_to_pixel16(rgb.r), float_to_pixel16(rgb.g),
               float_to_pixel16(rgb.b));
}

/**
 * @brief Wraps an angle in radians to (-pi, pi].
 * @param x Angle to wrap; any magnitude.
 * @return The equivalent angle in (-pi, pi].
 */
inline float wrap_angle_pi(float x) {
  while (x > PI_F)
    x -= 2.0f * PI_F;
  while (x < -PI_F)
    x += 2.0f * PI_F;
  return x;
}

/**
 * @brief Interpolates two OKLCH colors along the shortest-arc hue.
 * @param a Start color at t == 0.
 * @param b End color at t == 1.
 * @param t Blend weight; may extrapolate outside [0, 1].
 * @return The interpolated OKLCH color with L and C clamped valid.
 * @details An extrapolated t can overshoot a valid endpoint into an invalid
 * OKLCH: negative L renders near-black, negative C flips the hue 180°. L and C
 * are clamped valid; hue is left free to wrap.
 */
inline OKLCH lerp_oklch(OKLCH a, OKLCH b, float t) {
  // An achromatic (gray) endpoint has no meaningful hue angle: take the
  // chromatic endpoint's hue for the whole segment. If both ends are gray,
  // pin it to 0.
  float h;
  if (a.C < 1e-4f && b.C < 1e-4f) {
    h = 0.0f;
  } else if (a.C < 1e-4f) {
    h = b.h;
  } else if (b.C < 1e-4f) {
    h = a.h;
  } else {
    h = a.h + wrap_angle_pi(b.h - a.h) * t;
  }
  float L = hs::clamp(a.L + (b.L - a.L) * t, 0.0f, 1.0f);
  float C = std::max(0.0f, a.C + (b.C - a.C) * t);
  return {L, C, h};
}

/**
 * @brief Converts an OKLCH color to an 8-bit sRGB CPixel, chroma-reducing
 *        out-of-gamut colors to hold hue and lightness (Ottosson clip).
 * @param lch Source color in OKLCH space.
 * @return The gamut-mapped color as an 8-bit sRGB CPixel.
 * @details A vivid color past the sRGB cusp gets its chroma pulled to the
 *          boundary, holding hue and lightness, rather than the per-channel clip
 *          that twists hue toward a primary. In-gamut colors pay only the gate
 *          test.
 */
inline CPixel oklch_to_cpixel(OKLCH lch) {
  LinRGB rgb = oklab_to_linear_rgb_gamut(oklch_to_oklab(lch));
  return CPixel(linear_float_to_srgb8(rgb.r), linear_float_to_srgb8(rgb.g),
                linear_float_to_srgb8(rgb.b));
}

/**
 * @brief Abstract base for all palettes.
 * @details Uniform color-lookup interface via a single vtable pointer.
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
 * filled at construction by interpolating between the color stops in OKLCH
 * perceptual space.
 */
class Gradient : public Palette {
public:
  Pixel entries[256];

  /**
   * @brief Builds the 256-entry LUT by interpolating between color stops.
   * @param points Sorted-ascending (position in [0,1], color) stops.
   * @details Stop positions index entries[256] via static_cast<int>(pos * 255 + 0.5f);
   * a pos outside [0,1] is an out-of-bounds write, and segments fill only when
   * end > start, so an unsorted pair degenerates silently. Bounds and ordering
   * are trapped always-on (construction is cold).
   */
  Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {

    if (points.size() == 0)
      return;

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

    // Flat fills bake through the same OKLCH path as the segment endpoints, so
    // a flat region matches its adjacent segment at the shared stop.
    int firstStop = static_cast<int>(prevPos * 255.0f + 0.5f);
    Pixel prevSolid =
        oklch_to_pixel(srgb_to_oklch(prevColor.r, prevColor.g, prevColor.b));
    for (int i = 0; i <= firstStop; i++)
      entries[i] = prevSolid;

    it++;
    while (it != points.end()) {
      float nextPos = it->first;
      CPixel nextColor = it->second;

      int start = static_cast<int>(prevPos * 255.0f + 0.5f);
      int end = static_cast<int>(nextPos * 255.0f + 0.5f);

      // end == start (two stops quantizing to the same index) is the intended
      // "hard stop" — an abrupt color boundary, not a dropped stop.
      if (end > start) {
        OKLCH a = srgb_to_oklch(prevColor.r, prevColor.g, prevColor.b);
        OKLCH b = srgb_to_oklch(nextColor.r, nextColor.g, nextColor.b);
        for (int i = start; i < end; i++) {
          float t = static_cast<float>(i - start) / (end - start);
          entries[i] = oklch_to_pixel(lerp_oklch(a, b, t));
        }
      }
      prevPos = nextPos;
      prevColor = nextColor;
      it++;
    }

    int lastStop = static_cast<int>(prevPos * 255.0f + 0.5f);
    Pixel lastSolid =
        oklch_to_pixel(srgb_to_oklch(prevColor.r, prevColor.g, prevColor.b));
    for (int i = lastStop; i < 256; i++)
      entries[i] = lastSolid;
  }

  /**
   * @brief LUT lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The interpolated color (alpha 1.0).
   * @details Interpolated, not nearest-index, to avoid visible banding.
   */
  Color4 get(float t) const override {
    // Clamp before the int cast: t < 0 is float->int UB, t > 1 indexes past the
    // table; lo+1 only overruns at the t==1 endpoint, handled below.
    float idx = hs::clamp(t, 0.0f, 1.0f) * 255.0f;
    int lo = static_cast<int>(idx);
    if (lo >= 255) return Color4(entries[255], 1.0f);
    float frac = idx - lo;
    return Color4(entries[lo].lerp16(entries[lo + 1], frac_to_q16(frac)), 1.0f);
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
  HS_COLD_MEMBER GenerativePalette()
      : GenerativePalette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                          BrightnessProfile::FLAT) {}

  /**
   * @brief Builds a palette from three pre-authored sRGB key colors.
   * @param gradient_shape Shape/distribution of colors across the domain.
   * @param ka First key color (stop a).
   * @param kb Second key color (stop b).
   * @param kc Third key color (stop c).
   * @details RNG-free: the keys are supplied directly rather than sampled from
   *          profiles.
   */
  HS_COLD_MEMBER GenerativePalette(GradientShape gradient_shape, const CPixel &ka,
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
  HS_COLD_MEMBER GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type,
                    BrightnessProfile profile,
                    SaturationProfile sat_profile = SaturationProfile::MID,
                    int manual_seed = -1)
      : gradient_shape(gradient_shape) {
    uint8_t palette_hue;
    if (manual_seed >= 0) {
      palette_hue = static_cast<uint8_t>(manual_seed);
    } else {
      palette_hue = g_hue_seed;
      // Weyl recurrence mod 256, step 157 (= trunc(INV_PHI*255)): 157 is coprime
      // with 256 so the cursor visits every residue; rounding to 158 short-cycles.
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
      v3 = hs::rand_int(204, 256);
      break;
    case BrightnessProfile::DESCENDING:
      v1 = hs::rand_int(204, 256);
      v2 = hs::rand_int(127, 178);
      v3 = hs::rand_int(25, 76);
      break;
    case BrightnessProfile::FLAT:
      v1 = v2 = v3 = 255;
      break;
    case BrightnessProfile::BELL:
      v1 = hs::rand_int(51, 127);
      v2 = hs::rand_int(178, 256);
      v3 = v1;
      break;
    case BrightnessProfile::CUP:
      v1 = hs::rand_int(178, 256);
      v2 = hs::rand_int(51, 127);
      v3 = v1;
      break;
    }

    // The (h,s,v) integer draws above are reinterpreted as OKLCH coordinates,
    // not redrawn, so the global RNG stream stays unperturbed.
    a = author_key(h1, s1, v1);
    b = author_key(h2, s2, v2);
    c = author_key(h3, s3, v3);

    update_stops();
  }

  /**
   * @brief Builds a palette from three explicit, fully-resolved HSV key triples
   *        with NO RNG draws.
   * @param gradient_shape Shape/distribution of colors across the domain.
   * @param h1 First key hue in [0,255].
   * @param s1 First key saturation in [0,255].
   * @param v1 First key value in [0,255].
   * @param h2 Second key hue in [0,255].
   * @param s2 Second key saturation in [0,255].
   * @param v2 Second key value in [0,255].
   * @param h3 Third key hue in [0,255].
   * @param s3 Third key saturation in [0,255].
   * @param v3 Third key value in [0,255].
   * @return A palette whose keys are the OKLCH-authored (h,s,v) triples.
   * @details Mirrors the profile constructor's key authoring and stop layout,
   *          but takes (h,s,v) as parameters instead of sampling the global RNG.
   *          Entry point for the daydream palette tool's WASM bridge (PaletteOps),
   *          which owns its own randomization and asks only for the color math.
   */
  static GenerativePalette from_hsv_keys(GradientShape gradient_shape,
                                         uint8_t h1, uint8_t s1, uint8_t v1,
                                         uint8_t h2, uint8_t s2, uint8_t v2,
                                         uint8_t h3, uint8_t s3, uint8_t v3) {
    return GenerativePalette(gradient_shape, author_key(h1, s1, v1),
                             author_key(h2, s2, v2), author_key(h3, s3, v3));
  }

  // Peak OKLCH chroma at mid-lightness; the saturation profile scales it and a
  // sin(pi*L) envelope (key_oklch / get()) tapers it toward the lightness extremes.
  static constexpr float CHROMA_PEAK = 0.23f;

  // Hue torsion: radians of hue rotation per unit lightness, applied in get() as
  // h += HUE_TORSION * (L - 0.5), so shadows and highlights shift along the ramp.
  static constexpr float HUE_TORSION = 0.35f;

  // 8-bit hue wheel: harmony companion hues are fractions of the 256-value ring.
  static constexpr int HUE_WHEEL = 256;
  static constexpr int HUE_TRIADIC = HUE_WHEEL / 3;    // 1/3 turn (85)
  static constexpr int HUE_COMPLEMENT = HUE_WHEEL / 2; // 1/2 turn (128)

  // Perceptual lightness band a key's HSV value maps into: L = LIGHTNESS_FLOOR +
  // (val/255) * LIGHTNESS_SPAN, i.e. [0.12, 0.67] — below L=1 where the sRGB
  // gamut starves of chroma (see key_oklch).
  static constexpr float LIGHTNESS_FLOOR = 0.12f;
  static constexpr float LIGHTNESS_SPAN = 0.55f;

  /**
   * @brief OKLCH coordinates for one HSV-authored key. Single source of truth
   *        for the perceptual key placement.
   * @param hue Key hue in [0,255]; reinterpreted as an OKLCH hue angle.
   * @param sat Key saturation in [0,255]; scaled into OKLCH chroma.
   * @param val Key value in [0,255]; mapped into a perceptual lightness band.
   * @return The key's OKLCH coordinates (before the gamut-mapped sRGB bake).
   */
  static OKLCH key_oklch(uint8_t hue, uint8_t sat, uint8_t val) {
    // Even perceptual hue spacing: integer harmony offsets become true OKLCH-hue
    // offsets (triadic is a real 120deg).
    float h = (hue / 256.0f) * (2.0f * PI_F);
    float L = LIGHTNESS_FLOOR + (val / 255.0f) * LIGHTNESS_SPAN;
    // Chroma co-varies with lightness via a sin(pi*L) envelope; get() re-applies
    // the same envelope at the interpolated L.
    float C = (sat / 255.0f) * CHROMA_PEAK * fast_sinf(PI_F * L);
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
  HS_COLD_MEMBER static CPixel author_key(uint8_t hue, uint8_t sat, uint8_t val) {
    return oklch_to_cpixel(key_oklch(hue, sat, val));
  }

  /**
   * @brief Rebuilds stop positions/colors for the current gradient_shape.
   * @details Derives the stops from the three keys a/b/c, then caches their
   * OKLCH forms for get(). Call after any change to a/b/c or gradient_shape.
   */
  HS_COLD_MEMBER void update_stops() {
    const CPixel vignette_color(0, 0, 0);
    std::array<CPixel, MAX_STOPS> colors;
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
    default:
      HS_CHECK(false, "GenerativePalette: unknown gradient_shape");
    }
    HS_CHECK(size <= MAX_STOPS,
             "GenerativePalette: stop count exceeds parallel-array capacity");
    for (int i = 0; i < size; ++i) {
      colors_oklch[i] = srgb_to_oklch(colors[i].r, colors[i].g, colors[i].b);
      // Recover the stop's chroma ceiling cmax = C / sin(pi*L) so get() can
      // re-apply the envelope at the interpolated L (fast_sinf matches get()).
      // env -> 0 at L near 0/1 forces cmax = 0. A near-gray stop (C below
      // lerp_oklch's gray threshold) carries no meaningful hue; zero its cmax so
      // get()'s envelope + hue torsion can't bloom a tint into midtone grays.
      float env = fast_sinf(PI_F * colors_oklch[i].L);
      colors_cmax[i] = (colors_oklch[i].C >= 1e-4f && env > 1e-3f)
                           ? colors_oklch[i].C / env
                           : 0.0f;
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
   * @brief Sets keys to the OKLCH interpolation of two snapshots.
   * @param from Source snapshot at amount == 0.
   * @param to Source snapshot at amount == 1.
   * @param amount Blend weight.
   * @details Key hues travel coherent arcs: key a takes its shortest arc and
   * b/c rotate with it while their offsets from a morph along their own
   * shortest arcs. Independent per-key shortest arcs can counter-rotate,
   * sweeping adjacent stops past a half-turn mid-fade, where get()'s
   * shortest-arc segment lerp flips direction — a one-frame hue pop. Falls
   * back to per-key arcs when any key is near-gray (no meaningful hue).
   * Rebuilds the stops after interpolating.
   */
  HS_COLD_MEMBER void lerp(const Snapshot &from, const Snapshot &to, float amount) {
    // Clamp: an extrapolated amount overshoots into an invalid OKLCH (L > 1
    // or C past gamut).
    amount = hs::clamp(amount, 0.0f, 1.0f);
    const OKLCH fk[3] = {srgb_to_oklch(from.a.r, from.a.g, from.a.b),
                         srgb_to_oklch(from.b.r, from.b.g, from.b.b),
                         srgb_to_oklch(from.c.r, from.c.g, from.c.b)};
    const OKLCH tk[3] = {srgb_to_oklch(to.a.r, to.a.g, to.a.b),
                         srgb_to_oklch(to.b.r, to.b.g, to.b.b),
                         srgb_to_oklch(to.c.r, to.c.g, to.c.b)};
    CPixel *keys[3] = {&a, &b, &c};
    bool chromatic = true;
    for (int i = 0; i < 3; ++i)
      chromatic = chromatic && fk[i].C >= 1e-4f && tk[i].C >= 1e-4f;
    if (chromatic) {
      float d0 = wrap_angle_pi(tk[0].h - fk[0].h);
      for (int i = 0; i < 3; ++i) {
        float d = i == 0 ? d0
                         : d0 + wrap_angle_pi((tk[i].h - tk[0].h) -
                                              (fk[i].h - fk[0].h));
        OKLCH k = {
            hs::clamp(fk[i].L + (tk[i].L - fk[i].L) * amount, 0.0f, 1.0f),
            std::max(0.0f, fk[i].C + (tk[i].C - fk[i].C) * amount),
            fk[i].h + d * amount};
        *keys[i] = oklch_to_cpixel(k);
      }
    } else {
      for (int i = 0; i < 3; ++i)
        *keys[i] = oklch_to_cpixel(lerp_oklch(fk[i], tk[i], amount));
    }
    update_stops();
  }

  /**
   * @brief Samples the generated palette at a coordinate.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to 1.0).
   * @return The color at t (alpha 1.0).
   */
  Color4 get(float t) const override {
    // Clamp first: a t < shape[0] matches no segment and falls through to the
    // size-2 fallback below, returning the wrong stop (a discontinuity at t=0).
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

    float dist = end - start;
    // Zero-width segment: pin to the left stop (p=0) and render it through the
    // OKLCH path below, avoiding both a divide by ~0 and the 8-bit sRGB key.
    float p = dist < 0.0001f ? 0.0f : hs::clamp((t - start) / dist, 0.0f, 1.0f);

    // Interpolate in OKLCH; stops are pre-converted in update_stops().
    OKLCH blended = lerp_oklch(colors_oklch[seg], colors_oklch[seg + 1], p);
    // Re-derive chroma on the envelope at the interpolated L (see colors_cmax).
    // fast_sinf matches update_stops().
    float cmax = colors_cmax[seg] + (colors_cmax[seg + 1] - colors_cmax[seg]) * p;
    // Floor at 0: a small negative from fast_sinf would flip hue 180° in oklch_to_oklab.
    blended.C = std::max(0.0f, cmax * fast_sinf(PI_F * blended.L));
    // Hue torsion: drift hue with lightness, centered at L=0.5.
    blended.h += HUE_TORSION * (blended.L - 0.5f);
    return Color4(oklch_to_pixel(blended), 1.0f);
  }

  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~GenerativePalette() {}

  /**
   * @brief Pins the global generative-hue cursor.
   * @param seed Hue value in [0, 255] to set the cursor to; defaults to 0.
   * @details Auto-seeded palettes draw their base hue from `g_hue_seed` and
   * advance it (golden-ratio step); deliberately stateful, never reset in
   * production. The test harness calls this to restore identical global state
   * before re-rendering an effect.
   */
  static void reset_hue_seed(uint8_t seed = 0) { g_hue_seed = seed; }

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
      h2 = wrap_hue(h1_int + HUE_TRIADIC);
      h3 = wrap_hue(h1_int + 2 * HUE_TRIADIC);
      break;
    case HarmonyType::SPLIT_COMPLEMENTARY: {
      const int complement = wrap_hue(h1_int + HUE_COMPLEMENT);
      const int offset = 21;
      h2 = wrap_hue(complement - offset);
      h3 = wrap_hue(complement + offset);
      break;
    }
    case HarmonyType::COMPLEMENTARY: {
      h2 = wrap_hue(h1_int + HUE_COMPLEMENT);
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

  // Static cursor shared across instances: each auto-seeded construction
  // (manual_seed < 0) reads and advances it for a distinct base hue. Non-atomic;
  // palette construction is single-threaded.
  static inline uint8_t g_hue_seed = 0;

  // Single capacity bound for the parallel per-stop arrays below; `size` (set
  // only by update_stops) selects the live prefix shared by all of them.
  static constexpr int MAX_STOPS = 5;
  std::array<float, MAX_STOPS> shape;
  /**
   * @brief OKLCH forms of `colors`, cached by update_stops().
   * @details Lets the per-sample get() hot path skip the sRGB->OKLCH conversion
   * (6 powf + 6 cbrtf + 2 atan2f per stop).
   */
  std::array<OKLCH, MAX_STOPS> colors_oklch;
  /**
   * @brief Per-stop chroma ceiling, cached by update_stops().
   * @details Each stop's chroma is authored as C = cmax * sin(pi*L) (see
   * key_oklch). Caching cmax = C/sin(pi*L) lets get() re-apply the envelope at
   * the interpolated L, keeping midpoints on the sin curve, not a chord.
   */
  std::array<float, MAX_STOPS> colors_cmax{};
  int size = 0;
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
   * @details Computes color in float sRGB space, then converts to 16-bit linear
   * via the interpolated LUT, avoiding 8-bit quantization without a per-channel
   * powf.
   */
  Color4 get(float t) const override {
    float r_srgb = a[0] + b[0] * fast_cosf(2 * PI_F * (c[0] * t + d[0]));
    float g_srgb = a[1] + b[1] * fast_cosf(2 * PI_F * (c[1] * t + d[1]));
    float b_srgb = a[2] + b[2] * fast_cosf(2 * PI_F * (c[2] * t + d[2]));

    Pixel color(srgb_to_linear_interp(r_srgb), srgb_to_linear_interp(g_srgb),
                srgb_to_linear_interp(b_srgb));
    return Color4(color, 1.0f);
  }

  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~ProceduralPalette() {}

protected:
  std::array<float, 3> a, b, c, d;
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
      : ProceduralPalette(a1, b1, c1, d1),
        a1(a1), b1(b1), c1(c1), d1(d1), a2(a2), b2(b2), c2(c2), d2(d2) {
    mutate(0.0f);
  }

  /**
   * @brief Sets the active cosine parameters to the endpoint interpolation.
   * @param t Blend weight in [0, 1] between the start and end parameter sets.
   */
  void mutate(float t) {
    for (int i = 0; i < 3; ++i) {
      a[i] = hs::lerp(a1[i], a2[i], t);
      b[i] = hs::lerp(b1[i], b2[i], t);
      c[i] = hs::lerp(c1[i], c2[i], t);
      d[i] = hs::lerp(d1[i], d2[i], t);
    }
  }

private:
  std::array<float, 3> a1, b1, c1, d1;
  std::array<float, 3> a2, b2, c2, d2;

public:
  /**
   * @brief Trivial constexpr destructor.
   */
  constexpr ~MutatingPalette() {}
};

#include "color/composition.h"
