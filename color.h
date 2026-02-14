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

using Pixel = CRGB;

 // TODO: 3D Palettes

 ///////////////////////////////////////////////////////////////////////////////

 /**
  * @brief Represents a color with an alpha channel.
  */
struct Color4 {
  Pixel color;
  float alpha;

  Color4() : color(CRGB(0,0,0)), alpha(1.0f) {}
  Color4(Pixel p, float a = 1.0f) : color(p), alpha(a) {}
  Color4(uint8_t r, uint8_t g, uint8_t b, float a = 1.0f) : color(r, g, b), alpha(a) {}
  Color4(Color4 c, float a) : color(c.color), alpha(a) {}

  // Implicit conversion to CRGB (Pixel)
  operator Pixel() const { return color; }
};

/**
* @brief Scales a Pixel's RGB components by a floating-point scalar.
*/
Pixel operator*(const Pixel& p, float s) {
  return Pixel(
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.r * s), 0, 255)),
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.g * s), 0, 255)),
    static_cast<uint8_t>(std::clamp(static_cast<int>(p.b * s), 0, 255))
  );
}

/**
 * @brief Scales a Pixel's RGB components by a floating-point scalar in place.
 */
Pixel& operator*=(Pixel& p, float s) {
  p = p * s;
  return p;
}

/**
 * @brief sRGB to Linear RGB lookup table (IEC 61966-2-1).
 * Matches THREE.Color.convertSRGBToLinear() implementation.
 */
const uint8_t srgb_to_linear_lut[256] = {
   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,
   1,   1,   2,   2,   2,   2,   2,   2,   2,   2,   3,   3,   3,   3,   3,   3,
   4,   4,   4,   4,   4,   5,   5,   5,   5,   6,   6,   6,   6,   7,   7,   7,
   8,   8,   8,   8,   9,   9,   9,  10,  10,  10,  11,  11,  12,  12,  12,  13,
  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  17,  18,  18,  19,  19,  20,
  20,  21,  22,  22,  23,  23,  24,  24,  25,  25,  26,  27,  27,  28,  29,  29,
  30,  30,  31,  32,  32,  33,  34,  35,  35,  36,  37,  37,  38,  39,  40,  41,
  41,  42,  43,  44,  45,  45,  46,  47,  48,  49,  50,  51,  51,  52,  53,  54,
  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
  71,  72,  73,  74,  76,  77,  78,  79,  80,  81,  82,  84,  85,  86,  87,  88,
  90,  91,  92,  93,  95,  96,  97,  99, 100, 101, 103, 104, 105, 107, 108, 109,
 111, 112, 114, 115, 116, 118, 119, 121, 122, 124, 125, 127, 128, 130, 131, 133,
 134, 136, 138, 139, 141, 142, 144, 146, 147, 149, 151, 152, 154, 156, 157, 159,
 161, 163, 164, 166, 168, 170, 171, 173, 175, 177, 179, 181, 183, 184, 186, 188,
 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220,
 222, 224, 226, 229, 231, 233, 235, 237, 239, 242, 244, 246, 248, 250, 253, 255
};

/**
 * @brief Converts a color from sRGB to Linear space.
 */
Pixel convert_srgb_to_linear(const Pixel& p) {
  return CRGB(
    srgb_to_linear_lut[p.r],
    srgb_to_linear_lut[p.g],
    srgb_to_linear_lut[p.b]
  );
}

///////////////////////////////////////////////////////////////////////////////
// Blending Functions
///////////////////////////////////////////////////////////////////////////////

Pixel blend_max(const Pixel& c1, const Pixel& c2) {
  return Pixel(std::max(c1.r, c2.r), std::max(c1.g, c2.g), std::max(c1.b, c2.b));
}

Pixel blend_over(const Pixel& c1, const Pixel& c2) {
  return c2;
}

Pixel blend_under(const Pixel& c1, const Pixel& c2) {
  return c1;
}

Pixel blend_add(const Pixel& c1, const Pixel& c2) {
  return Pixel(qadd8(c1.r, c2.r), qadd8(c1.g, c2.g), qadd8(c1.b, c2.b));
}

auto blend_alpha(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      c1.r * (1.0f - a) + c2.r * a,
      c1.g * (1.0f - a) + c2.g * a,
      c1.b * (1.0f - a) + c2.b * a);
    };
}

auto blend_accumulate(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r, c2.r * a),
      qadd8(c1.g, c2.g * a),
      qadd8(c1.b, c2.b * a));
    };
}

Pixel blend_over_max(const Pixel& c1, const Pixel& c2) {
  float m1 = sqrtf(c1.r * c1.r + c1.g * c1.g + c1.b * c1.b);
  float m2 = sqrtf(c2.r * c2.r + c2.g * c2.g + c2.b * c2.b);
  if (m2 == 0) return c1;
  float s = std::max(m1, m2) / m2;
  return c2 * s;
}

Pixel blend_over_min(const Pixel& c1, const Pixel& c2) {
  float m1 = sqrtf(c1.r * c1.r + c1.g * c1.g + c1.b * c1.b);
  float m2 = sqrtf(c2.r * c2.r + c2.g * c2.g + c2.b * c2.b);
  if (m2 == 0) return Pixel(0, 0, 0);
  float s = std::min(m1, m2) / m2;
  return c2 * s;
}

Pixel blend_mean(const Pixel& c1, const Pixel& c2) {
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

auto blend_add_alpha(float a) {
    return [a](const Pixel& dest, const Pixel& src) {
        return Pixel(
            qadd8(dest.r, static_cast<uint8_t>(src.r * a)),
            qadd8(dest.g, static_cast<uint8_t>(src.g * a)),
            qadd8(dest.b, static_cast<uint8_t>(src.b * a))
        );
    };
}

auto blend_max_alpha(float a) {
    return [a](const Pixel& dest, const Pixel& src) {
        return Pixel(
            std::max(dest.r, static_cast<uint8_t>(src.r * a)),
            std::max(dest.g, static_cast<uint8_t>(src.g * a)),
            std::max(dest.b, static_cast<uint8_t>(src.b * a))
        );
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
