/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <variant>
#include <vector>
#include "platform.h"
#include "3dmath.h"
#include "util.h"
#include "static_circular_buffer.h"

#ifndef PROGMEM
#define PROGMEM
#endif

// Forward declarations
struct Pixel16;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// LUT Generation (constexpr)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Helper forward declaration
inline uint16_t srgb_to_linear(uint8_t srgb);

// Gamma 2.2 approximation
static constexpr float GAMMA = 2.2f;
static constexpr float INV_GAMMA = 1.0f / 2.2f;

// Helper to calculate sRGB -> Linear (0-255 -> 0-65535)
constexpr uint16_t srgb_to_linear_calc(uint8_t srgb) {
    // We can't use powf in constexpr in C++17 easily without a constexpr math library.
    // However, for a LUT we might need to rely on initialization time or just generic math if C++20.
    // Since we want this to be PROGMEM and potentially compile-time:
    // Simple iterative power or approximate if std::pow isn't constexpr compliant in the environment.
    // BUT checking standard: std::pow is NOT constexpr.
    // We will assume the compiler can evaluate this if we write a simple constexpr pow or just init at runtime?
    // User requested "Store the LUT in PROGMEM as constexpr". 
    // "constexpr" variable implies compile-time eval.
    // Let's try to use a simple recursive/iterative approach or just hardcode for common values?
    // No, 256/65536 is too big.
    // Actually, for C++14/17 we often use a generator struct. 
    // Re-reading user request: "Store the LUT in PROGMEM as constexpr". 
    // If I cannot make it strictly constexpr due to pow(), I will make it const and initialize it.
    // But PROGMEM variables must be const.
    // Let's use a small constexpr pow helper.
    
    // Approximation or simplistic implementation for constexpr:
    // x^2.2 is roughly x^2.2.
    // Let's use standard initialization for now, keeping it simple.
    // If implicit constexpr fails, we'll drop constexpr and just use const PROGMEM.
    
    float v = srgb / 255.0f;
    // float res = pow(v, 2.2f); // Not constexpr
    // Use x*x as decent approx if strict constexpr needed? 
    // User asked for "accurate linear_to_Srgb() conversion".
    // I will generate the tables using a generator class/lambda trick if possible, 
    // OR just define them and populate. 
    // But populating a const PROGMEM array is tricky without constexpr.
    
    // For now, I'll define the arrays and populate them with a helper if supported, 
    // or assume the environment supports non-constexpr init of globals (runtime before main).
    // BUT PROGMEM usually requires compile time constants for flash placement on some archs.
    
    // Actually, keeping it as 'const' with an initializer list is the standard way.
    // Since writing out 65536 entries is impossible, I will use a class with a constexpr constructor that fills a std::array?
    // That works in C++17/20.
    return 0; // Placeholder, logic below uses std::array generator
}

namespace detail {
    constexpr float my_pow(float base, float exp) {
        // Very basic integer-ish power or taylor series? 
        // Too complex for reliable constexpr.
        // Let's presume we can't do full constexpr pow.
        return base; // Dummy
    }
}

// NOTE: Since generating 64KB constant data via template metaprogramming or constexpr 
// complex math is fragile across compilers, I will implement the lookups as 
// `const` arrays populated by a helper class constructor if possible, 
// OR simpler: Use the user's suggestion of "Accuracy" vs "Constexpr" trade-off.
// With C++20, we can allocate vector in constexpr.
// Given constraints, I will implement the Pixel16 class and use runtime initialization 
// for the tables if constexpr pow is not available, 
// UNLESS I include a large precomputed list, which is bad style here.
//
// RETRACTION: The user asked for "constexpr". I will assume they have a modern compiler/setup 
// that can handle `std::array` initialization.
// Since I cannot call `pow` in constexpr, I will use a simple Gamma 2.0 (square) for constexpr 
// OR just leave them as non-constexpr `const` tables initialized at startup (RAM).
// Wait, `PROGMEM` implies flash. On AVR/Teensy, `PROGMEM` variables must be global constants.
// I will provide the tables and a `init_color()` function, OR use a table generator if one handles `pow`.
//
// ALTERNATIVE: Use `pow(x, 2)` (x*x) for sRGB->Linear (Gamma 2.0). 
// And `sqrt(x)` for Linear->sRGB (Gamma 2.0). 
// This is "accurate enough" for many consistency cases and IS constexpr-able.
// BUT user asked for "Accurate linear_to_Srgb". 
// I will stick to the existing plan's logic: 
// 1. `Pixel16` class. 2. `PROGMEM` LUTs. 
// I'll assume we can init them.

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct Pixel16 {
    uint16_t r, g, b;

    constexpr Pixel16() : r(0), g(0), b(0) {}
    constexpr Pixel16(uint16_t _r, uint16_t _g, uint16_t _b) : r(_r), g(_g), b(_b) {}
    
    // Construct from HSV (converts to sRGB then Linear)
    Pixel16(const CHSV& hsv) {
        CRGB srgb(hsv);
        r = srgb_to_linear(srgb.r);
        g = srgb_to_linear(srgb.g);
        b = srgb_to_linear(srgb.b);
    }
    
    // Construct from CRGB (converts to Linear)
    Pixel16(const CRGB& c) {
        r = srgb_to_linear(c.r);
        g = srgb_to_linear(c.g);
        b = srgb_to_linear(c.b);
    }

    // Explicit forward declaration of conversion
    operator CRGB() const;

    // Basic arithmetic
    Pixel16& operator+=(const Pixel16& rhs) {
        r = (uint16_t)std::min((uint32_t)65535, (uint32_t)r + rhs.r);
        g = (uint16_t)std::min((uint32_t)65535, (uint32_t)g + rhs.g);
        b = (uint16_t)std::min((uint32_t)65535, (uint32_t)b + rhs.b);
        return *this;
    }
    
    Pixel16 operator*(float s) const {
        return Pixel16(
            (uint16_t)std::clamp((int)(r * s), 0, 65535),
            (uint16_t)std::clamp((int)(g * s), 0, 65535),
            (uint16_t)std::clamp((int)(b * s), 0, 65535)
        );
    }
    
    // Lerp (Linear interpolation 16-bit)
    Pixel16 lerp16(const Pixel16& other, uint16_t frac) const {
        // frac 0..65535
        // result = this + (other - this) * frac
        // Exact: (a * (65535 - frac) + b * frac) / 65535
        uint32_t r32 = ((uint32_t)r * (65535 - frac) + (uint32_t)other.r * frac) / 65535;
        uint32_t g32 = ((uint32_t)g * (65535 - frac) + (uint32_t)other.g * frac) / 65535;
        uint32_t b32 = ((uint32_t)b * (65535 - frac) + (uint32_t)other.b * frac) / 65535;
        return Pixel16((uint16_t)r32, (uint16_t)g32, (uint16_t)b32);
    }
};

using Pixel = Pixel16;

 // TODO: 3D Palettes

 ///////////////////////////////////////////////////////////////////////////////

 /**
  * @brief Represents a color with an alpha channel.
  */
struct Color4 {
  Pixel color;
  float alpha;

  Color4() : color(Pixel(0,0,0)), alpha(1.0f) {}
  Color4(Pixel p, float a = 1.0f) : color(p), alpha(a) {}
  Color4(uint8_t r, uint8_t g, uint8_t b, float a = 1.0f) : color(Pixel(srgb_to_linear(r), srgb_to_linear(g), srgb_to_linear(b))), alpha(a) {}
  Color4(Color4 c, float a) : color(c.color), alpha(a) {}

  // Implicit conversion to CRGB (Pixel)
  operator CRGB() const { return color; }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// LUT Storage (Declarations)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Include the generated LUTs
// This file is generated by generate_luts.py and contains:
// - srgb_to_linear_lut[256]
// - linear_to_srgb_lut[65536]
#include "color_luts.h"

// Helper: sRGB (8-bit) -> Linear (16-bit)
inline uint16_t srgb_to_linear(uint8_t srgb) {
    return srgb_to_linear_lut[srgb];
}

// Conversion implementation
inline Pixel16::operator CRGB() const {
    return CRGB(
        linear_to_srgb_lut[r],
        linear_to_srgb_lut[g],
        linear_to_srgb_lut[b]
    );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Blending Functions (Updated for Pixel16)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline Pixel blend_max(const Pixel& c1, const Pixel& c2) {
  return Pixel(std::max(c1.r, c2.r), std::max(c1.g, c2.g), std::max(c1.b, c2.b));
}

inline Pixel blend_over(const Pixel& c1, const Pixel& c2) {
  return c2;
}

inline Pixel blend_under(const Pixel& c1, const Pixel& c2) {
  return c1;
}

inline Pixel blend_add(const Pixel& c1, const Pixel& c2) {
    // Saturated Add
    uint32_t r = (uint32_t)c1.r + c2.r;
    uint32_t g = (uint32_t)c1.g + c2.g;
    uint32_t b = (uint32_t)c1.b + c2.b;
    return Pixel(
        (r > 65535) ? 65535 : (uint16_t)r,
        (g > 65535) ? 65535 : (uint16_t)g,
        (b > 65535) ? 65535 : (uint16_t)b
    );
}

inline auto blend_alpha(float a) {
  uint16_t ai = std::clamp((int)(a * 65535), 0, 65535);
  return [ai](const Pixel& c1, const Pixel& c2) {
      return c1.lerp16(c2, ai);
  };
}

inline auto blend_accumulate(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
      // c1 + c2 * a
      Pixel16 added = c2 * a;
      return blend_add(c1, added);
  };
}

inline Pixel blend_over_max(const Pixel& c1, const Pixel& c2) {
  // Magnitude comparison (approximate luminance) in linear space
  // We can just use sum or max component for speed, or euclidean.
  // Linear RGB euclidean:
  float m1 = (float)c1.r * c1.r + (float)c1.g * c1.g + (float)c1.b * c1.b; 
  float m2 = (float)c2.r * c2.r + (float)c2.g * c2.g + (float)c2.b * c2.b;
  if (m2 < 0.001f) return c1;
  
  // Just return max? The logic in original was mixing them?
  // "s = max(m1, m2) / m2; return c2 * s;" -> Rescales c2 to match max brightness
  if (m1 > m2) {
       float s = sqrtf(m1 / m2);
       return c2 * s;
  }
  return c2;
}

inline Pixel blend_over_min(const Pixel& c1, const Pixel& c2) {
  float m1 = (float)c1.r * c1.r + (float)c1.g * c1.g + (float)c1.b * c1.b; 
  float m2 = (float)c2.r * c2.r + (float)c2.g * c2.g + (float)c2.b * c2.b;
  if (m2 < 0.001f) return Pixel(0,0,0);
  
  if (m1 < m2) {
       float s = sqrtf(m1 / m2);
       return c2 * s;
  }
  return c2;
}

inline Pixel blend_mean(const Pixel& c1, const Pixel& c2) {
  return Pixel((c1.r + c2.r) / 2, (c1.g + c2.g) / 2, (c1.b + c2.b) / 2);
}

enum BlendMode : uint8_t {
    BLEND_OVER = 0,
    BLEND_ADD = 1,
    BLEND_MAX = 2
};

struct TaggedColor {
    Pixel color;
    float alpha;
    uint8_t tag;
};

inline auto blend_add_alpha(float a) {
    return [a](const Pixel& dest, const Pixel& src) {
        return blend_add(dest, src * a);
    };
}

inline auto blend_max_alpha(float a) {
    return [a](const Pixel& dest, const Pixel& src) {
        Pixel16 s = src * a;
        return Pixel(std::max(dest.r, s.r), std::max(dest.g, s.g), std::max(dest.b, s.b));
    };
}

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Base class for all color palette systems.
 */
class Palette {
public:
  virtual Color4 get(float t) const = 0;
  constexpr virtual ~Palette() {}
};

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
 * @brief Defines the visual shape or distribution of colors across the palette domain.
 */
enum class GradientShape {
  STRAIGHT,
  CIRCULAR,
  VIGNETTE,
  FALLOFF
};

/**
 * @brief Defines the overall brightness profile across the palette domain.
 */
enum class BrightnessProfile {
  ASCENDING,
  DESCENDING,
  FLAT,
  BELL,
  CUP
};

/**
 * @brief Defines the saturation profile.
 */
enum class SaturationProfile {
  PASTEL,
  MID,
  VIBRANT
};

/**
 * @brief Converts a float in the range [0.0, 1.0] to a 16-bit integer for FastLED lerping.
 */
uint16_t to_short(float zero_to_one) {
  return std::clamp(static_cast<int>(std::round(zero_to_one * 65535.0f)), 0, 65535);
}

///////////////////////////////////////////////////////////////////////////////

// Helper: constexpr linear interpolation for bytes
constexpr uint8_t lerp8(uint8_t a, uint8_t b, float t) {
    return static_cast<uint8_t>(a + (b - a) * t);
}

/**
 * @brief A constexpr-compatible RGB pixel structure for Flash storage.
 * Layout compatible with CRGB but without non-constexpr constructors.
 */
struct CPixel {
    uint8_t r, g, b;
    constexpr CPixel() : r(0), g(0), b(0) {}
    constexpr CPixel(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
    constexpr CPixel(uint32_t hex) : r((hex >> 16) & 0xFF), g((hex >> 8) & 0xFF), b(hex & 0xFF) {}
    
    // Convert to FastLED CRGB (Pixel)
    operator Pixel() const { return CRGB(r, g, b); }
};

/**
 * @brief A class representing a discrete color gradient/lookup table.
 */
class Gradient : public Palette {
public:
    CPixel entries[256];

    // Constexpr constructor calculates the table at compile-time
    constexpr Gradient(std::initializer_list<std::pair<float, CPixel>> points) : entries() {
        // Initialize with black
        for(int i=0; i<256; i++) entries[i] = CPixel(0,0,0);

        if (points.size() == 0) return;

        auto it = points.begin();
        float prevPos = it->first;
        CPixel prevColor = it->second;
        
        // Fill start
        int startIdx = 0;
        int firstStop = static_cast<int>(prevPos * 255);
        for(int i = 0; i <= firstStop; i++) entries[i] = prevColor;

        it++;
        while(it != points.end()) {
            float nextPos = it->first;
            CPixel nextColor = it->second;
            
            int start = static_cast<int>(prevPos * 255);
            int end = static_cast<int>(nextPos * 255);
            
            if (end > start) {
                for (int i = start; i <= end; i++) {
                    float t = static_cast<float>(i - start) / (end - start);
                    entries[i].r = lerp8(prevColor.r, nextColor.r, t);
                    entries[i].g = lerp8(prevColor.g, nextColor.g, t);
                    entries[i].b = lerp8(prevColor.b, nextColor.b, t);
                }
            }
            prevPos = nextPos;
            prevColor = nextColor;
            it++;
        }
        
        // Fill end
        int lastStop = static_cast<int>(prevPos * 255);
        for(int i = lastStop; i < 256; i++) entries[i] = prevColor;
    }

    Color4 get(float t) const override {
        uint8_t index = static_cast<uint8_t>(t * 255);
        CPixel p = entries[index];
        Pixel color(p.r, p.g, p.b);
        return Color4(color, 1.0f);
    }

    constexpr ~Gradient() override {}
};

/**
 * @brief A palette that generates colors based on defined harmony and brightness profiles.
 */
class GenerativePalette : public Palette {
public:

  GenerativePalette() : 
    GenerativePalette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS, BrightnessProfile::FLAT) {}

  GenerativePalette(
    GradientShape gradient_shape,
    HarmonyType harmony_type,
    BrightnessProfile profile,
    SaturationProfile sat_profile = SaturationProfile::MID,
    int manual_seed = -1) :
    gradient_shape(gradient_shape),
    harmony_type(harmony_type)
  {
    if (manual_seed != -1) {
        palette_hue = static_cast<uint8_t>(manual_seed);
    } else {
        palette_hue = g_hue_seed;
        // Advance global seed for next generation (Golden Ratio distribution)
        g_hue_seed = static_cast<uint8_t>(
            (static_cast<uint32_t>(g_hue_seed) + static_cast<uint32_t>(G * 255.0f)) % 256);
    }

    uint8_t h1 = palette_hue;
    uint8_t h2, h3;

    calc_hues(h1, h2, h3, harmony_type);

    uint8_t s1, s2, s3;
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

    uint8_t v1, v2, v3;
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

    a = CHSV(h1, s1, v1);
    b = CHSV(h2, s2, v2);
    c = CHSV(h3, s3, v3);

    update_luts();
  }

  void update_luts() {
    const Pixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = { 0, 0.1f, 0.5f, 0.9f, 1.0f };
      colors = { vignette_color, a, b, c, vignette_color };
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = { 0, 0.5f, 1.0f };
      colors = { a, b, c };
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = { 0, 0.33f, 0.66f, 1.0f };
      colors = { a, b, c, a };
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = { 0, 0.33f, 0.66f, 0.9f, 1.0f };
      colors = { a, b, c, vignette_color };
      size = 4;
      break;
    }
  }

  void lerp(const GenerativePalette& from, const GenerativePalette& to, float amount) {
    uint16_t fract = to_short(amount);
    a = from.a.lerp16(to.a, fract);
    b = from.b.lerp16(to.b, fract);
    c = from.c.lerp16(to.c, fract);
    update_luts();
  }

  Color4 get(float t) const override {
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0) seg = size - 2;

    float start = shape[seg];
    float end = shape[seg + 1];
    Pixel c1 = colors[seg];
    Pixel c2 = colors[seg + 1];

    // Safe division check
    float dist = end - start;
    if (dist < 0.0001f) return Color4(c1, 1.0f);

    float p = (t - start) / dist;
    Pixel color = c1.lerp16(c2, to_short(std::clamp(p, 0.0f, 1.0f)));
    return Color4(color, 1.0f);
  }

private:
  uint8_t wrap_hue(int hue) const { return (hue % 256 + 256) % 256; }

  void calc_hues(uint8_t h1, uint8_t& h2, uint8_t& h3, HarmonyType harmony_type) const {
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
      const int offset = hs::rand_int(-7, 7);
      h3 = wrap_hue(h1_int + offset);
      break;
    }
    case HarmonyType::ANALOGOUS:
    default: {
      const int dir = (hs::rand_int(0, 1) == 0) ? 1 : -1;
      const int offset1 = dir * hs::rand_int(11, 22);
      h2 = wrap_hue(h1_int + offset1);
      const int offset2 = dir * hs::rand_int(11, 22);
      h3 = wrap_hue(h2 + offset2);
      break;
    }
    }
  }

  GradientShape gradient_shape;
  HarmonyType harmony_type;
  uint8_t palette_hue;
  Pixel a, b, c;
  
  static inline uint8_t g_hue_seed = 0;
  std::array<float, 5> shape;
  std::array<Pixel, 5> colors;
  int size = 0;

public:
  constexpr ~GenerativePalette() override {}
};

/**
 * @brief A palette defined by a mathematical cosine wave function.
 * C(t) = A + B * cos(2 * PI * (C * t + D))
 */
class ProceduralPalette : public Palette {
public:
  constexpr ProceduralPalette(
    std::array<float, 3> a,
    std::array<float, 3> b,
    std::array<float, 3> c,
    std::array<float, 3> d) : a(a), b(b), c(c), d(d) {
  }

  Color4 get(float t) const override {
    Pixel color(
      static_cast<uint8_t>(255 * std::clamp(a[0] + b[0] * cosf(2 * PI_F * (c[0] * t + d[0])), 0.0f, 1.0f)),
      static_cast<uint8_t>(255 * std::clamp(a[1] + b[1] * cosf(2 * PI_F * (c[1] * t + d[1])), 0.0f, 1.0f)),
      static_cast<uint8_t>(255 * std::clamp(a[2] + b[2] * cosf(2 * PI_F * (c[2] * t + d[2])), 0.0f, 1.0f))
    );
    return Color4(color, 1.0f);
  }

protected:
  std::array<float, 3> a, b, c, d;

public:
  constexpr ~ProceduralPalette() override {}
};

/**
 * @brief A palette that allows continuous mutation between two procedural palettes.
 */
class MutatingPalette : public ProceduralPalette {
public:
  MutatingPalette(
    std::array<float, 3> a1, std::array<float, 3> b1, std::array<float, 3> c1, std::array<float, 3> d1,
    std::array<float, 3> a2, std::array<float, 3> b2, std::array<float, 3> c2, std::array<float, 3> d2
  ) : ProceduralPalette(a1, b1, c1, d1), // Initialize base with start values
    a1(a1), b1(b1), c1(c1), d1(d1),
    a2(a2), b2(b2), c2(c2), d2(d2)
  {
    mutate(0.0f);
  }

  void mutate(float t) {
    for (int i = 0; i < 3; ++i) {
      a[i] = lerp(a1[i], a2[i], t);
      b[i] = lerp(b1[i], b2[i], t);
      c[i] = lerp(c1[i], c2[i], t);
      d[i] = lerp(d1[i], d2[i], t);
    }
  }

private:
  float lerp(float x, float y, float t) { return x * (1.0f - t) + y * t; }
  std::array<float, 3> a1, b1, c1, d1;
  std::array<float, 3> a2, b2, c2, d2;

public:
  constexpr ~MutatingPalette() override {}
};

///////////////////////////////////////////////////////////////////////////////
// Palette Wrappers
///////////////////////////////////////////////////////////////////////////////

struct PaletteModifier {
    virtual void step() = 0;
    virtual float transform(float t) const = 0;
    virtual ~PaletteModifier() = default;
};

/**
 * @brief Linearly cycles the palette coordinate.
 */
struct CycleModifier : public PaletteModifier {
    float speed;
    float offset = 0.0f;
    
    CycleModifier(float speed = 0.01f) : speed(speed) {}
    
    void step() override {
        offset += speed;
        // Keep offset within [0, 1] to preserve floating point precision over time
        if (offset > 1.0f) offset -= 1.0f;
        if (offset < 0.0f) offset += 1.0f;
    }
    
    float transform(float t) const override {
        // fmodf handles wrapping for us
        return fmodf(t + offset, 1.0f);
    }
};

/**
 * @brief Oscillates the palette coordinate (Breathing).
 */
struct BreatheModifier : public PaletteModifier {
    float freq;
    float amp;
    float phase = 0.0f;
    float val = 0.0f;
    
    BreatheModifier(float freq = 0.05f, float amp = 0.1f) : freq(freq), amp(amp) {}
    
    void step() override {
        phase += freq;
        if (phase > 2 * PI_F) phase -= 2 * PI_F;
        val = sinf(phase) * amp;
    }
    
    float transform(float t) const override {
        return std::clamp(t + val, 0.0f, 1.0f);
    }
};

/**
 * @brief A wrapper that applies a chain of modifiers to the palette parameter 't'.
 */
class AnimatedPalette : public Palette {
public:
  AnimatedPalette(const Palette& source) : source(source) {}

  /**
   * Explicitly connect a modifier to this palette.
   * Stores a pointer. The modifier must outlive the palette.
   */
  AnimatedPalette& add(PaletteModifier* modifier) {
    if (!modifiers.is_full()) {
        modifiers.push_back(modifier);
    }
    return *this;
  }

  Color4 get(float t) const override {
    float final_t = t;
    // Pipe t through the modifier chain: t -> m1 -> m2 ... -> t'
    for (size_t i = 0; i < modifiers.size(); ++i) {
        if (modifiers[i]) {
            final_t = modifiers[i]->transform(final_t);
        }
    }
    return source.get().get(final_t);
  }

private:
  std::reference_wrapper<const Palette> source;
  StaticCircularBuffer<PaletteModifier*, 4> modifiers;
};



class CircularPalette : public Palette {
public:
  CircularPalette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    return palette.get().get(1.0f - std::abs(2.0f * t - 1.0f));
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class ReversePalette : public Palette {
public:
  ReversePalette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    return palette.get().get(1.0f - t);
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class VignettePalette : public Palette {
public:
  VignettePalette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    CRGB vignette_color(0, 0, 0);
    if (t < 0.2f) {
      return Color4(vignette_color.lerp16(palette.get().get(0.0f).color, to_short(quintic_kernel(t / 0.2f))), 1.0f);
    }
    else if (t >= 0.8f) {
      return Color4(palette.get().get(1.0f).color.lerp16(vignette_color, to_short(quintic_kernel((t - 0.8f) / 0.2f))), 1.0f);
    }
    else {
      return palette.get().get((t - 0.2f) / 0.6f);
    }
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class TransparentVignette : public Palette {
public:
  TransparentVignette(const Palette& palette) : palette(palette) {}
  Color4 get(float t) const override {
    Color4 result;
    float factor = 1.0f;
    if (t < 0.2f) {
      result = palette.get().get(0.0f);
      factor = quintic_kernel(t / 0.2f);
    }
    else if (t >= 0.8f) {
      result = palette.get().get(1.0f);
      factor = quintic_kernel(1.0f - (t - 0.8f) / 0.2f);
    }
    else {
      return palette.get().get((t - 0.2f) / 0.6f);
    }
    result.alpha *= factor;
    return result;
  }
private:
  std::reference_wrapper<const Palette> palette;
};

class AlphaFalloffPalette : public Palette {
public:
    using FalloffFn = std::function<float(float)>;
    AlphaFalloffPalette(FalloffFn fn, const Palette& source) : fn(fn), source(source) {}
    Color4 get(float t) const override {
        Color4 c = source.get().get(t);
        c.alpha *= fn(t);
        return c;
    }
private:
    FalloffFn fn;
    std::reference_wrapper<const Palette> source;
};

class SolidColorPalette : public Palette {
public:
  SolidColorPalette(const Color4& color) : color(color) {}
  Color4 get(float t) const override {
    return color;
  }
  Color4 color;
};

/**
 * @brief Variant holding any supported palette type.
 */
using PaletteVariant = std::variant<
  std::monostate,
  GenerativePalette,
  ProceduralPalette,
  Gradient,
  MutatingPalette,
  ReversePalette,
  VignettePalette,
  TransparentVignette,
  AlphaFalloffPalette,
  SolidColorPalette,
  CircularPalette
>;


inline Color4 get_color(const PaletteVariant& pv, float t) {
  return std::visit([t](auto&& arg) -> Color4 {
    using T = std::decay_t<decltype(arg)>;
    if constexpr (std::is_same_v<T, std::monostate>) {
      return Color4(CRGB(0, 0, 0), 0.0f);
    } else {
      return arg.get(t);
    }
  }, pv);
}


///////////////////////////////////////////////////////////////////////////////
