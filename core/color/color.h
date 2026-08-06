/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file color.h
 * @brief Pixel and Color4, the linear-light color types, with blending, sRGB
 *        conversion and the palette core.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <utility>

#include "engine/platform.h"
#include "math/3dmath.h"
#include "engine/util.h"

#include "engine/memory.h"

#include "color/gamut_lut.h"
#include "color/srgb_decode.h"

#if defined(__ARM_FEATURE_DSP)
// Inline assembly avoids a CMSIS header dependency for the saturating add.
__attribute__((always_inline)) static inline uint32_t
inline_uqadd16(uint32_t a, uint32_t b) {
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
   * @param _r Red channel in [0, 65535].
   * @param _g Green channel in [0, 65535].
   * @param _b Blue channel in [0, 65535].
   */
  constexpr Pixel(uint16_t _r, uint16_t _g, uint16_t _b)
      : r(_r), g(_g), b(_b) {}

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
    return Pixel((uint16_t)hs::clamp(r * s + 0.5f, 0.0f, 65535.0f),
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
   * @brief Tests two pixels for channel inequality.
   * @param rhs Pixel to compare against.
   * @return True if any channel differs.
   */
  bool operator!=(const Pixel &rhs) const { return !(*this == rhs); }

  /**
   * @brief Tests equality against an HSV color (converted to Pixel).
   * @param rhs Color in HSV space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CHSV &rhs) const { return *this == Pixel(rhs); }

  /**
   * @brief Tests inequality against an HSV color (converted to Pixel).
   * @param rhs Color in HSV space.
   * @return True if this pixel differs from the converted color.
   */
  bool operator!=(const CHSV &rhs) const { return !(*this == rhs); }

  /**
   * @brief Tests equality against a CRGB color (converted to Pixel).
   * @param rhs Color in 8-bit sRGB space.
   * @return True if this pixel equals the converted color.
   */
  bool operator==(const CRGB &rhs) const { return *this == Pixel(rhs); }

  /**
   * @brief Tests inequality against a CRGB color (converted to Pixel).
   * @param rhs Color in 8-bit sRGB space.
   * @return True if this pixel differs from the converted color.
   */
  bool operator!=(const CRGB &rhs) const { return !(*this == rhs); }
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
 * @brief Alpha below which a fragment cannot move a pixel.
 * @details One 8-bit LSB at full-scale color. A primitive whose peak alpha falls
 * under this rasterizes to nothing, so callers gate on it instead of drawing.
 */
inline constexpr float MIN_VISIBLE_ALPHA = 1.0f / 255.0f;

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

enum class HueMode : uint8_t { HARMONY, SWEEP, CUSTOM };

enum class PaletteHarmony : uint8_t {
  MONOCHROMATIC,
  ANALOGOUS,
  ACCENTED_ANALOGOUS,
  COMPLEMENTARY,
  SPLIT_COMPLEMENTARY,
  TRIADIC,
};

enum class HueDirection : uint8_t { SHORTEST, CLOCKWISE, COUNTERCLOCKWISE };

enum class AxisCurve : uint8_t {
  CONSTANT,
  ASCENDING,
  DESCENDING,
  BELL,
  CUP,
  CUSTOM,
};

enum class ChromaBasis : uint8_t { LOCAL_GAMUT, PATH_MINIMUM, ABSOLUTE };

enum class ColorPath : uint8_t { OKLCH_ARC, OKLAB_CARTESIAN };

enum class PaletteDomain : uint8_t {
  STRAIGHT,
  MIRROR,
  VIGNETTE,
  FALLOFF,
  LOOP,
};

enum class SegmentEase : uint8_t { LINEAR, COSINE, SMOOTHSTEP };

struct HueControls {
  HueMode mode = HueMode::HARMONY;
  PaletteHarmony harmony = PaletteHarmony::ANALOGOUS;
  HueDirection direction = HueDirection::SHORTEST;
  float base_turns = 0.0f;
  float spread_turns = 0.07f;
  float sweep_turns = 1.0f;
  std::array<float, 3> custom_turns{};
};

struct AxisControls {
  AxisCurve curve = AxisCurve::CONSTANT;
  float center = 0.62f;
  float range = 0.0f;
  std::array<float, 3> custom{};
};

struct ChromaControls {
  AxisCurve curve = AxisCurve::CONSTANT;
  ChromaBasis basis = ChromaBasis::LOCAL_GAMUT;
  float center = 0.62f;
  float range = 0.0f;
  float headroom = 0.94f;
  std::array<float, 3> custom{};
};

struct PaletteRecipe {
  static constexpr uint8_t SCHEMA_VERSION = 2;

  uint8_t schema_version = SCHEMA_VERSION;
  PaletteDomain domain = PaletteDomain::STRAIGHT;
  SegmentEase easing = SegmentEase::COSINE;
  ColorPath color_path = ColorPath::OKLCH_ARC;
  HueControls hue;
  AxisControls lightness;
  ChromaControls chroma;
  float hue_torsion = 0.0f;
  float falloff_start = 0.90f;
};

enum class PaletteCompileCode : uint8_t {
  OK,
  INVALID_SCHEMA,
  NON_FINITE,
  INVALID_ENUM,
  HUE_LIMIT,
  NON_INTEGER_LOOP_SWEEP,
  INVALID_FALLOFF_START,
  INCOMPATIBLE_OPTIONS,
};

enum class PaletteRecipeField : uint8_t {
  NONE = 0,
  PALETTE_DOMAIN = 1,
  EASING = 2,
  COLOR_PATH = 3,
  HUE_MODE = 4,
  HARMONY = 5,
  HUE_DIRECTION = 6,
  BASE_TURNS = 7,
  SPREAD_TURNS = 8,
  SWEEP_TURNS = 9,
  CUSTOM_TURNS_0 = 10,
  CUSTOM_TURNS_1 = 11,
  CUSTOM_TURNS_2 = 12,
  LIGHTNESS_CURVE = 13,
  LIGHTNESS_CENTER = 14,
  LIGHTNESS_RANGE = 15,
  LIGHTNESS_CUSTOM_0 = 16,
  LIGHTNESS_CUSTOM_1 = 17,
  LIGHTNESS_CUSTOM_2 = 18,
  CHROMA_CURVE = 19,
  CHROMA_BASIS = 20,
  CHROMA_CENTER = 21,
  CHROMA_RANGE = 22,
  CHROMA_CUSTOM_0 = 23,
  CHROMA_CUSTOM_1 = 24,
  CHROMA_CUSTOM_2 = 25,
  CHROMA_HEADROOM = 26,
  HUE_TORSION = 27,
  FALLOFF_START = 28,
  SCHEMA_VERSION = 29,
};

struct PaletteAdjustments {
  uint64_t wrapped_fields = 0;
  uint64_t clamped_fields = 0;
  uint64_t canonicalized_fields = 0;
};

struct PaletteCompileStatus {
  PaletteCompileCode code = PaletteCompileCode::OK;
  PaletteRecipeField field = PaletteRecipeField::NONE;
  PaletteAdjustments adjustments;
};

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
  return (l <= 0.0031308f) ? l * 12.92f
                           : 1.055f * powf(l, 1.0f / 2.4f) - 0.055f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// OKLab / OKLCH Color Space (Björn Ottosson, 2020)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * @brief OKLab perceptual color: lightness L and chroma axes a, b.
 */
struct OKLab {
  float L, a, b;
};
/**
 * @brief OKLCH polar color: lightness L, chroma C, hue h (radians).
 */
struct OKLCH {
  float L, C, h;
};

/** @brief Cone-response (LMS) triple, the OKLab intermediate before the cube-root
 *  nonlinearity. */
struct LMS {
  float l, m, s;
};

/** @brief Linear-RGB triple in [0,1] (may sit slightly out of gamut before
 *  clamping). */
struct LinRGB {
  float r, g, b;
};

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
 * @param l_cbrt Cube-rooted l cone response.
 * @param m_cbrt Cube-rooted m cone response.
 * @param s_cbrt Cube-rooted s cone response.
 * @return The color in OKLab space.
 * @details Takes the already-cube-rooted triple so the caller picks the
 * cube-root flavour.
 */
HS_O3_FN inline OKLab lms_to_oklab(float l_cbrt, float m_cbrt, float s_cbrt) {
  return {
      0.2104542553f * l_cbrt + 0.7936177850f * m_cbrt - 0.0040720468f * s_cbrt,
      1.9779984951f * l_cbrt - 2.4285922050f * m_cbrt + 0.4505937099f * s_cbrt,
      0.0259040371f * l_cbrt + 0.7827717662f * m_cbrt - 0.8086757660f * s_cbrt};
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
 * @brief Converts linear RGB [0,1] to OKLab through fast_cbrt.
 * @param r Linear red in [0, 1].
 * @param g Linear green in [0, 1].
 * @param b Linear blue in [0, 1].
 * @return The color in OKLab space, accurate to fast_cbrt.
 */
__attribute__((always_inline)) inline OKLab
linear_rgb_to_oklab_fast(float r, float g, float b) {
  LMS lms = linear_rgb_to_lms(r, g, b);
  return lms_to_oklab(fast_cbrt(lms.l), fast_cbrt(lms.m), fast_cbrt(lms.s));
}

/**
 * @brief Converts OKLab to cube-rooted LMS (the inverse OKLab matrix).
 * @param lab Source color in OKLab space.
 * @param l_cbrt Out: cube-rooted l cone response.
 * @param m_cbrt Out: cube-rooted m cone response.
 * @param s_cbrt Out: cube-rooted s cone response.
 */
HS_O3_FN inline void oklab_to_lms_cbrt(OKLab lab, float &l_cbrt, float &m_cbrt,
                                       float &s_cbrt) {
  l_cbrt = lab.L + 0.3963377774f * lab.a + 0.2158037573f * lab.b;
  m_cbrt = lab.L - 0.1055613458f * lab.a - 0.0638541728f * lab.b;
  s_cbrt = lab.L - 0.0894841775f * lab.a - 1.2914855480f * lab.b;
}

/**
 * @brief Converts cube-rooted LMS to linear RGB [0,1] (cube + RGB matrix).
 * @param l_cbrt Cube-rooted l cone response.
 * @param m_cbrt Cube-rooted m cone response.
 * @param s_cbrt Cube-rooted s cone response.
 * @param r Out: linear red (may exit gamut before clamping).
 * @param g Out: linear green.
 * @param b Out: linear blue.
 */
HS_O3_FN
inline void lms_cbrt_to_linear_rgb(float l_cbrt, float m_cbrt, float s_cbrt,
                                   float &r, float &g, float &b) {
  float l = l_cbrt * l_cbrt * l_cbrt, m = m_cbrt * m_cbrt * m_cbrt,
        s = s_cbrt * s_cbrt * s_cbrt;

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
HS_O3_FN inline void oklab_to_linear_rgb(OKLab lab, float &r, float &g,
                                         float &b) {
  float l_cbrt, m_cbrt, s_cbrt;
  oklab_to_lms_cbrt(lab, l_cbrt, m_cbrt, s_cbrt);
  lms_cbrt_to_linear_rgb(l_cbrt, m_cbrt, s_cbrt, r, g, b);
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

// Chroma pulled back off the refined crossing; without it the caller's own
// re-conversion rounds a channel a part in a million past the gate.
inline constexpr float GAMUT_CLIP_MARGIN = 2e-5f;

// Chroma below which an OKLCH color has no usable hue angle and is handled as
// gray. Every path that classifies a color as achromatic tests against this.
inline constexpr float OKLCH_ACHROMATIC_C = 1e-4f;

/**
 * @brief Gamut boundary bracket grid and the scales indexing it.
 * @details Defaults to the full-resolution flash master, so the clip path is
 * always usable and no effect has to opt in to being correct. Grouped in one
 * object so the per-pixel path loads a single base address rather than five
 * unrelated globals.
 */
struct GamutLut {
  /** @brief Flash master, or an arena copy once one is armed. */
  const uint16_t *table = GAMUT_LUT;
  /** @brief Diamond-angle buckets over [0, 4). */
  int angle_steps = GAMUT_LUT_ANGLE_STEPS;
  /** @brief Lightness buckets over [0, 1]. */
  int l_steps = GAMUT_LUT_L_STEPS;
  /** @brief angle_steps / 4, the index scale. */
  float angle_scale = GAMUT_LUT_ANGLE_STEPS * 0.25f;
  /** @brief l_steps, the index scale. */
  float l_scale = GAMUT_LUT_L_STEPS;
};

/**
 * @brief The single live boundary grid.
 * @details Points at the flash master until an effect arms an arena copy, which
 * only buys read latency: the scattered per-pixel reads land in RAM rather than
 * QSPI flash. Worth it only at per-pixel call rates. configure_arenas() restores
 * the flash default before the storage under a copy is handed out again.
 */
inline GamutLut g_gamut_lut;

/**
 * @brief Downsamples GAMUT_LUT into @p arena and points the clip path at it.
 * @param arena Arena to hold the copy; configure_arenas() drops the pointer
 *        before this storage is handed out again.
 * @param angle_steps Diamond-angle buckets; must divide GAMUT_LUT_ANGLE_STEPS.
 * @param l_steps Lightness buckets; must divide GAMUT_LUT_L_STEPS.
 * @details Optional, and only worth its arena bytes at per-pixel clip rates:
 * the flash master already serves the clip correctly and at higher resolution.
 * Call after the arenas are configured, from the owning effect's init(). A
 * coarse cell takes the minimum of the merged minima and the maximum of the
 * merged maxima, so the true boundary of every ray in the cell still lies
 * inside the stored bracket at any resolution. Cost in arena bytes is
 * gamut_lut_bytes(angle_steps, l_steps); resolution only sets how wide the
 * bracket starts, and the per-pixel bisection sets how far it is narrowed.
 */
HS_COLD_MEMBER inline void init_gamut_lut(Arena &arena, int angle_steps,
                                          int l_steps) {
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
 * @brief Drops any arena copy and points the clip path back at the flash master.
 * @details Runs before persistent storage is handed out again, so no owner can
 * leave a pointer into freed storage behind. The clip stays correct across the
 * swap; only read latency changes.
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
 * @details With l_cbrt = L + A*u and (L + X*u)^3 expanded in u, every linear-RGB
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

  const float r3 =
      4.0767416621f * ka3 - 3.3077115913f * km3 + 0.2309699292f * ks3;
  const float r2 =
      c * (4.0767416621f * ka2 - 3.3077115913f * km2 + 0.2309699292f * ks2);
  const float r1 =
      q * (4.0767416621f * ka - 3.3077115913f * km + 0.2309699292f * ks);
  const float g3 =
      -1.2684380046f * ka3 + 2.6097574011f * km3 - 0.3413193965f * ks3;
  const float g2 =
      c * (-1.2684380046f * ka2 + 2.6097574011f * km2 - 0.3413193965f * ks2);
  const float g1 =
      q * (-1.2684380046f * ka + 2.6097574011f * km - 0.3413193965f * ks);
  const float b3 =
      -0.0041960863f * ka3 - 0.7034186147f * km3 + 1.7076147010f * ks3;
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
    // The whole bracket held: the crossing is at or above hi, and hi is capped
    // at the cell maximum, which bounds it from above.
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
 * @brief Returns the first sRGB gamut-boundary chroma at an OKLCH coordinate.
 * @param L OKLab lightness, clamped to [0,1].
 * @param h Hue in radians.
 * @return The first-exit chroma minus the shared numerical margin.
 */
HS_O3_FN __attribute__((noinline)) inline float gamut_max_chroma(float L,
                                                                 float h) {
  const GamutLut &lut = g_gamut_lut;
  L = hs::clamp(L, 0.0f, 1.0f);
  if (L == 0.0f || L == 1.0f)
    return 0.0f;

  const float a = cosf(h);
  const float b = sinf(h);
  int ai = static_cast<int>(diamond_angle(b, a) * lut.angle_scale);
  int li = static_cast<int>(L * lut.l_scale);
  ai = hs::clamp(ai, 0, lut.angle_steps - 1);
  li = hs::clamp(li, 0, lut.l_steps - 1);
  const uint16_t *cell = &lut.table[(li * lut.angle_steps + ai) * 2];
  const float c_lo = static_cast<float>(cell[0]) * GAMUT_LUT_INV_SCALE;
  const float c_hi = static_cast<float>(cell[1]) * GAMUT_LUT_INV_SCALE;
  const float boundary = gamut_bracket_refine(L, a, b, c_lo, c_hi);
  return std::max(0.0f, boundary - GAMUT_CLIP_MARGIN);
}

/**
 * @brief Reduces OKLab chroma to the first sRGB gamut boundary.
 * @param lab Source color; lightness is clamped to [0,1].
 * @return The source color or its fixed-lightness, fixed-hue projection.
 */
HS_O3_FN __attribute__((noinline)) inline OKLab
gamut_clip_preserve_chroma(OKLab lab) {
  lab.L = hs::clamp(lab.L, 0.0f, 1.0f);
  const float c_sq = lab.a * lab.a + lab.b * lab.b;
  if (!(c_sq > 1e-12f))
    return {lab.L, 0.0f, 0.0f};

  const float C = sqrtf(c_sq);
  const float C_MAX = gamut_max_chroma(lab.L, atan2f(lab.b, lab.a));
  if (C <= C_MAX)
    return lab;
  const float scale = C_MAX / C;
  return {lab.L, lab.a * scale, lab.b * scale};
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
  if (!linear_rgb_in_gamut(r, g, b)) {
    HS_PROFILE_DEEP(gamut_clip);
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r, g, b);
  }
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
    float l_cbrt, m_cbrt, s_cbrt;
    oklab_to_lms_cbrt({lab.L, a2, b2}, l_cbrt, m_cbrt, s_cbrt);
    k[i] = l_cbrt;
    k[3 + i] = m_cbrt;
    k[6 + i] = s_cbrt;
  }
}

/**
 * @brief Applies a cbrt-LMS 3x3 (from hue_rotate_lms_matrix, optionally
 * uniformly scaled) and converts to linear RGB with gamut mapping.
 * @param k Row-major 3x3 acting on cube-rooted LMS.
 * @param l_cbrt Cube-rooted l cone response.
 * @param m_cbrt Cube-rooted m cone response.
 * @param s_cbrt Cube-rooted s cone response.
 * @param r Out: gamut-mapped linear red (may sit a hair past the bound;
 * callers still clamp).
 * @param g Out: linear green.
 * @param b Out: linear blue.
 * @details In-gamut colors never leave cbrt-LMS; the OKLab form is recomputed
 * only on the chroma-clip slow path.
 */
HS_O3_FN
inline void lms_cbrt_transform_rgb(const float k[9], float l_cbrt, float m_cbrt,
                                   float s_cbrt, float &r, float &g, float &b) {
  float ul = k[0] * l_cbrt + k[1] * m_cbrt + k[2] * s_cbrt;
  float um = k[3] * l_cbrt + k[4] * m_cbrt + k[5] * s_cbrt;
  float us = k[6] * l_cbrt + k[7] * m_cbrt + k[8] * s_cbrt;
  lms_cbrt_to_linear_rgb(ul, um, us, r, g, b);
  if (!linear_rgb_in_gamut(r, g, b)) {
    HS_PROFILE_DEEP(gamut_clip);
    OKLab lab = lms_to_oklab(ul, um, us);
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r, g, b);
  }
}

/**
 * @brief Two-pixel lms_cbrt_transform_rgb sharing one code path.
 * @param k Row-major 3x3 acting on cube-rooted LMS.
 * @param l0 Cube-rooted l cone response of the first pixel.
 * @param m0 Cube-rooted m of the first pixel.
 * @param s0 Cube-rooted s of the first pixel.
 * @param l1 Cube-rooted l of the second pixel.
 * @param m1 Cube-rooted m of the second pixel.
 * @param s1 Cube-rooted s of the second pixel.
 * @param r0 Out: gamut-mapped linear red of the first pixel.
 * @param g0 Out: linear green of the first pixel.
 * @param b0 Out: linear blue of the first pixel.
 * @param r1 Out: linear red of the second pixel.
 * @param g1 Out: linear green of the second pixel.
 * @param b1 Out: linear blue of the second pixel.
 * @details Results match two lms_cbrt_transform_rgb calls bit for bit; only the
 * statement order differs, so an in-order FPU can overlap the two independent
 * chains. The uniform work stays interleaved and only the rare chroma-clip
 * fixups run per pixel.
 */
HS_O3_FN
inline void lms_cbrt_transform_rgb2(const float k[9], float l0, float m0,
                                    float s0, float l1, float m1, float s1,
                                    float &r0, float &g0, float &b0, float &r1,
                                    float &g1, float &b1) {
  float ul0 = k[0] * l0 + k[1] * m0 + k[2] * s0;
  float ul1 = k[0] * l1 + k[1] * m1 + k[2] * s1;
  float um0 = k[3] * l0 + k[4] * m0 + k[5] * s0;
  float um1 = k[3] * l1 + k[4] * m1 + k[5] * s1;
  float us0 = k[6] * l0 + k[7] * m0 + k[8] * s0;
  float us1 = k[6] * l1 + k[7] * m1 + k[8] * s1;
  lms_cbrt_to_linear_rgb(ul0, um0, us0, r0, g0, b0);
  lms_cbrt_to_linear_rgb(ul1, um1, us1, r1, g1, b1);
  bool ok0 = linear_rgb_in_gamut(r0, g0, b0);
  bool ok1 = linear_rgb_in_gamut(r1, g1, b1);
  if (!ok0) {
    HS_PROFILE_DEEP(gamut_clip);
    OKLab lab = lms_to_oklab(ul0, um0, us0);
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r0, g0, b0);
  }
  if (!ok1) {
    HS_PROFILE_DEEP(gamut_clip);
    OKLab lab = lms_to_oklab(ul1, um1, us1);
    oklab_to_linear_rgb(gamut_clip_preserve_chroma(lab), r1, g1, b1);
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
 * @brief Normalizes a 16-bit linear Pixel to a linear-RGB [0,1] triple.
 * @param p Source color in 16-bit linear space.
 * @return The three channels scaled by 1/65535.
 */
__attribute__((always_inline)) inline LinRGB pixel_to_linrgb(const Pixel &p) {
  constexpr float INV16 = 1.0f / 65535.0f;
  return {p.r * INV16, p.g * INV16, p.b * INV16};
}

/**
 * @brief Quantizes a [0,1] linear channel to an 8-bit sRGB component.
 * @param l Linear channel value; clamped to [0, 1].
 * @return The channel as an 8-bit sRGB value in [0, 255].
 */
inline uint8_t linear_float_to_srgb8(float l) {
  return static_cast<uint8_t>(
      hs::clamp(linear_to_srgb_float(hs::clamp(l, 0.0f, 1.0f)) * 255.0f + 0.5f,
                0.0f, 255.0f));
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
HS_O3_FN inline void hue_rotate_rgb(float &r, float &g, float &b, float ca,
                                    float sa) {
  OKLab lab = linear_rgb_to_oklab_fast(r, g, b);

  float a2 = lab.a * ca - lab.b * sa;
  float b2 = lab.a * sa + lab.b * ca;

  oklab_to_linear_rgb_gamut({lab.L, a2, b2}, r, g, b);
}

inline Color4 hue_rotate(const Color4 &c, float ca, float sa) {
  LinRGB rgb = pixel_to_linrgb(c.color);

  hue_rotate_rgb(rgb.r, rgb.g, rgb.b, ca, sa);

  Color4 result = c;
  result.color.r = float_to_pixel16(rgb.r);
  result.color.g = float_to_pixel16(rgb.g);
  result.color.b = float_to_pixel16(rgb.b);
  return result;
}

/**
 * @brief Cosine/sine of a turn angle from fast trig, renormalized to unit
 * length.
 * @param turns Angle in turns (0..1 = full turn).
 * @param ca Out: cosine of the angle.
 * @param sa Out: sine of the angle.
 * @details fast trig is non-orthonormal; renormalizing keeps a chroma rotation
 * length-preserving (else the scaling compounds per frame under feedback).
 */
__attribute__((always_inline)) inline void
turn_to_unit_cos_sin(float turns, float &ca, float &sa) {
  float angle = turns * (2.0f * PI_F);
  ca = fast_cosf(angle);
  sa = fast_sinf(angle);
  float inv = 1.0f / sqrtf(ca * ca + sa * sa);
  ca *= inv;
  sa *= inv;
}

inline Color4 hue_rotate(const Color4 &c, float amount) {
  float ca, sa;
  turn_to_unit_cos_sin(amount, ca, sa);
  return hue_rotate(c, ca, sa);
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
  LinRGB rgb = pixel_to_linrgb(c.color);
  OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
  return {lab, c};
}

/**
 * @brief Perceptual (OKLab) hue rotation of a precomputed base color.
 * @param hb Precomputed base from make_hue_rotate_base().
 * @param amount Rotation in turns (0..1 = full turn).
 * @return The hue-rotated color.
 */
inline Color4 hue_rotate(const HueRotateBase &hb, float amount) {
  float ca, sa;
  turn_to_unit_cos_sin(amount, ca, sa);

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
 * @brief Converts a 16-bit linear Pixel to OKLCH.
 * @param p Source color in 16-bit linear space.
 * @return The color in OKLCH space.
 */
inline OKLCH pixel_to_oklch(const Pixel &p) {
  LinRGB rgb = pixel_to_linrgb(p);
  return oklab_to_oklch(linear_rgb_to_oklab(rgb.r, rgb.g, rgb.b));
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
 * @brief Wraps an angle in radians to [-pi, pi].
 * @param x Angle to wrap; any magnitude.
 * @return The equivalent angle in [-pi, pi]. Both endpoints are reachable: an
 * exact half turn keeps the sign it arrived with. Callers use the result as a
 * shortest-arc delta, for which -pi and +pi are equivalent.
 */
inline float wrap_angle_pi(float x) {
  // At large |x| the subtraction below rounds away and the loop never finishes.
  // Negated so NaN takes this branch too.
  if (!(fabsf(x) <= 4.0f * PI_F))
    x = fmodf(x, 2.0f * PI_F);
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
  if (a.C < OKLCH_ACHROMATIC_C && b.C < OKLCH_ACHROMATIC_C) {
    h = 0.0f;
  } else if (a.C < OKLCH_ACHROMATIC_C) {
    h = b.h;
  } else if (b.C < OKLCH_ACHROMATIC_C) {
    h = a.h;
  } else {
    h = a.h + wrap_angle_pi(b.h - a.h) * t;
  }
  float L = hs::clamp(a.L + (b.L - a.L) * t, 0.0f, 1.0f);
  float C = std::max(0.0f, a.C + (b.C - a.C) * t);
  return {L, C, h};
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
  /**
   * @brief Builds the 256-entry LUT by interpolating between color stops.
   * @param points Sorted-ascending (position in [0,1], color) stops.
   * @details Emptiness, stop bounds and ordering are trapped always-on
   * (construction is cold).
   */
  Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {
    HS_CHECK(points.size() > 0, "Gradient requires at least one stop");

    float prev_check = -1.0f;
    for (const auto &stop : points) {
      HS_CHECK(stop.first >= 0.0f && stop.first <= 1.0f,
               "Gradient stop position out of [0,1]");
      HS_CHECK(stop.first >= prev_check,
               "Gradient stops must be sorted ascending");
      prev_check = stop.first;
    }

    auto it = points.begin();
    float prev_pos = it->first;
    CPixel prev_color = it->second;

    // Flat fills bake through the same OKLCH path as the segment endpoints, so
    // a flat region matches its adjacent segment at the shared stop.
    int first_stop = static_cast<int>(prev_pos * 255.0f + 0.5f);
    Pixel prev_solid =
        oklch_to_pixel(srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b));
    for (int i = 0; i <= first_stop; i++)
      entries[i] = prev_solid;

    it++;
    while (it != points.end()) {
      float next_pos = it->first;
      CPixel next_color = it->second;

      int start = static_cast<int>(prev_pos * 255.0f + 0.5f);
      int end = static_cast<int>(next_pos * 255.0f + 0.5f);

      // end == start (two stops quantizing to the same index) is the intended
      // "hard stop" — an abrupt color boundary, not a dropped stop.
      if (end > start) {
        OKLCH a = srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b);
        OKLCH b = srgb_to_oklch(next_color.r, next_color.g, next_color.b);
        for (int i = start; i < end; i++) {
          float t = static_cast<float>(i - start) / (end - start);
          entries[i] = oklch_to_pixel(lerp_oklch(a, b, t));
        }
      }
      prev_pos = next_pos;
      prev_color = next_color;
      it++;
    }

    int last_stop = static_cast<int>(prev_pos * 255.0f + 0.5f);
    Pixel last_solid =
        oklch_to_pixel(srgb_to_oklch(prev_color.r, prev_color.g, prev_color.b));
    for (int i = last_stop; i < 256; i++)
      entries[i] = last_solid;
  }

  /**
   * @brief LUT lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The interpolated color (alpha 1.0).
   * @details Interpolated, not nearest-index, to avoid visible banding.
   */
  Color4 get(float t) const override {
    // Clamp before the int cast: t < 0 is float->int UB and NaN maps to the last
    // entry, both of which lut_sample_pixel requires the caller to have excluded.
    return Color4(
        lut_sample_pixel(entries, 256, hs::clamp(t, 0.0f, 1.0f) * 255.0f),
        1.0f);
  }

private:
  Pixel entries[256];
};

#include "color/generative_palette.h"

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
      : ProceduralPalette(a1, b1, c1, d1), a1(a1), b1(b1), c1(c1), d1(d1),
        a2(a2), b2(b2), c2(c2), d2(d2) {
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

#define HS_COLOR_INTERNAL
#include "color/composition.h"
#undef HS_COLOR_INTERNAL
