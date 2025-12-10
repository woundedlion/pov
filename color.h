#pragma once

// TODO: 3D Palettes

///////////////////////////////////////////////////////////////////////////////

/**
* @brief Scales a Pixel's RGB components by a floating-point scalar.
* @param p The Pixel to scale.
* @param s The scaling factor.
* @return The scaled Pixel.
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
 * @param p The Pixel to scale.
 * @param s The scaling factor.
 * @return Reference to the scaled Pixel.
 */
Pixel& operator*=(Pixel& p, float s) {
  p = p * s;
  return p;
}

/**
 * @brief Gamma brightness lookup table (LUT).
 * @details This table is used to convert linear color values (what the program calculates)
 * to gamma-corrected values (what the human eye perceives) for display.
 * Generated using gamma = 2.20, steps = 256, range = 0-255.
 */
 // Gamma brightness lookup table <https://victornpb.github.io/gamma-table-generator>
 // gamma = 2.20 steps = 256 range = 0-255
const uint8_t gamma_lut[256] = {
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
     1,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,
     3,   3,   3,   3,   3,   4,   4,   4,   4,   5,   5,   5,   5,   6,   6,   6,
     6,   7,   7,   7,   8,   8,   8,   9,   9,   9,  10,  10,  11,  11,  11,  12,
    12,  13,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,  19,
    20,  20,  21,  22,  22,  23,  23,  24,  25,  25,  26,  26,  27,  28,  28,  29,
    30,  30,  31,  32,  33,  33,  34,  35,  35,  36,  37,  38,  39,  39,  40,  41,
    42,  43,  43,  44,  45,  46,  47,  48,  49,  49,  50,  51,  52,  53,  54,  55,
    56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
    73,  74,  75,  76,  77,  78,  79,  81,  82,  83,  84,  85,  87,  88,  89,  90,
    91,  93,  94,  95,  97,  98,  99, 100, 102, 103, 105, 106, 107, 109, 110, 111,
   113, 114, 116, 117, 119, 120, 121, 123, 124, 126, 127, 129, 130, 132, 133, 135,
   137, 138, 140, 141, 143, 145, 146, 148, 149, 151, 153, 154, 156, 158, 159, 161,
   163, 165, 166, 168, 170, 172, 173, 175, 177, 179, 181, 182, 184, 186, 188, 190,
   192, 194, 196, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221,
   223, 225, 227, 229, 231, 234, 236, 238, 240, 242, 244, 246, 248, 251, 253, 255,
};

/**
 * @brief Applies gamma correction to a CRGB pixel color.
 * @param p The input CRGB pixel.
 * @return The gamma-corrected CRGB pixel.
 */
Pixel gamma_correct(const Pixel& p) {
  return CRGB(
    gamma_lut[p.r],
    gamma_lut[p.g],
    gamma_lut[p.b]
  );
}
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Base class for all color palette systems.
 */
class Palette {
public:
  /**
   * @brief Gets the color for a given position/time along the palette.
   * @param t Normalized value (0.0 to 1.0).
   * @return The calculated CRGB pixel color.
   */
  virtual Pixel get(float t) const = 0;
};

/**
 * @brief Defines types of color harmonies for generative palettes.
 */
enum class HarmonyType {
  TRIADIC, /**< Three colors equidistant on the color wheel. */
  SPLIT_COMPLEMENTARY, /**< Base color and two colors adjacent to its complement. */
  COMPLEMENTARY, /**< Base color and its direct opposite. */
  ANALOGOUS /**< Closely spaced colors on the color wheel. */
};

/**
 * @brief Defines the visual shape or distribution of colors across the palette domain.
 */
enum class GradientShape {
  STRAIGHT, /**< Linear blend from start to end colors. */
  CIRCULAR, /**< Wraps the blend back to the start color. */
  VIGNETTE, /**< Fades to black at the edges. */
  FALLOFF /**< Fades only at the end. */
};

/**
 * @brief Defines the overall brightness profile across the palette domain.
 */
enum class BrightnessProfile {
  ASCENDING, /**< Brightness increases across the palette. */
  DESCENDING, /**< Brightness decreases across the palette. */
  FLAT, /**< Constant maximum brightness. */
  BELL, /**< Brightness starts low, peaks in the middle, and ends low. */
  CUP /**< Brightness starts high, dips in the middle, and ends high. */
};

/**
 * @brief Blends two CRGB pixels using the additive blend mode (clamping at 255).
 * @param c1 The first pixel.
 * @param c2 The second pixel.
 * @return The resulting blended pixel.
 */
auto blend_add(const Pixel& c1, const Pixel& c2) {
  return Pixel(
    qadd8(c1.r, c2.r),
    qadd8(c1.g, c2.g),
    qadd8(c1.b, c2.b));
}


/**
 * @brief Creates a functor for blending two CRGB pixels using a constant alpha value.
 * @param a The blending factor (0.0 to 1.0).
 * @return A functor that performs the alpha blend: c_result = c1 * (1-a) + c2 * a.
 */
auto blend_alpha(float a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r * (1 - a), c2.r * a),
      qadd8(c1.g * (1 - a), c2.g * a),
      qadd8(c1.b * (1 - a), c2.b * a));
    };
}

/**
 * @brief Converts a float in the range [0.0, 1.0] to a 16-bit integer for FastLED lerping.
 * @param zero_to_one The input float value.
 * @return The corresponding 16-bit integer (0 to 65535).
 */
uint16_t to_short(float zero_to_one) {
  return std::clamp(std::round(zero_to_one * 65535.0f), 0.0f, 65535.0f);
}

/**
 * @brief A palette that generates colors based on defined harmony and brightness profiles.
 * @details This palette stores three primary HSV colors (a, b, c) and interpolates between them
 * based on a gradient shape and a normalized position 't'.
 */
class GenerativePalette : public Palette {
public:

  /**
   * @brief Constructs a generative palette by selecting colors based on harmony and brightness rules.
   * @param gradient_shape The spatial shape of the blend (e.g., STRAIGHT, VIGNETTE).
   * @param harmony_type The relationship between the primary hues (e.g., TRIADIC, ANALOGOUS).
   * @param profile The overall brightness profile (e.g., ASCENDING, BELL).
   */
  GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type, BrightnessProfile profile) :
    gradient_shape(gradient_shape),
    harmony_type(harmony_type),
    seed_hue(static_cast<uint8_t>(hs::rand_int(0, 256)))
  {
    uint8_t h1 = seed_hue;
    uint8_t h2;
    uint8_t h3;

    seed_hue = static_cast<uint8_t>(
      (static_cast<uint32_t>(seed_hue) + static_cast<uint32_t>(G * 255.0)) % 256);
    calc_hues(h1, h2, h3, harmony_type);

    const uint8_t s1 = hs::rand_int(255 * 0.4, 255 * 0.8);
    const uint8_t s2 = hs::rand_int(255 * 0.4, 255 * 0.8);
    const uint8_t s3 = hs::rand_int(255 * 0.4, 255 * 0.8);

    uint8_t v1, v2, v3;
    switch (profile) {
    case BrightnessProfile::ASCENDING:
      v1 = hs::rand_int(255 * 0.1, 255 * 0.3);
      v2 = hs::rand_int(255 * 0.5, 255 * 0.7);
      v3 = hs::rand_int(255 * 0.8, 255 * 1.0);
      break;
    case BrightnessProfile::DESCENDING:
      v1 = hs::rand_int(255 * 0.8, 255);
      v2 = hs::rand_int(255 * 0.5, 255 * 0.7);
      v3 = hs::rand_int(255 * 0.1, 255 * 0.3);
      break;
    case  BrightnessProfile::FLAT:
      v1 = 255;
      v2 = 255;
      v3 = 255;
      break;
    case  BrightnessProfile::BELL:
      v1 = hs::rand_int(255 * 0.2, 255 * 0.5);
      v2 = hs::rand_int(255 * 0.7, 255);
      v3 = v1;
      break;
    case  BrightnessProfile::CUP:
		v1 = hs::rand_int(255 * 0.7, 255);
		v2 = hs::rand_int(255 * 0.2, 255 * 0.5);
		v3 = v1;
        break;
    }

    a = CHSV(h1, s1, v1);
    b = CHSV(h2, s2, v2);
    c = CHSV(h3, s3, v3);

    update_luts();
  }

  /**
   * @brief Recalculates the internal interpolation points based on the current gradient shape.
   */
  void update_luts() {
    const Pixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = { 0, 0.1, 0.5, 0.9, 1 };
      colors = { vignette_color, a, b, c, vignette_color };
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = { 0, 0.5, 1 };
      colors = { a, b, c };
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = { 0, 0.33, 0.66, 1 };
      colors = { a, b, c, a };
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = { 0, 0.33, 0.66, 0.9, 1 };
      colors = { a, b, c, vignette_color };
      size = 4;
      break;
    }
  }

  /**
   * @brief Interpolates the palette's internal HSV components toward a target palette.
   * @param from The starting palette state.
   * @param to The target palette state.
   * @param amount The interpolation factor (0.0 to 1.0).
   */
  void lerp(const GenerativePalette& from, const GenerativePalette& to, float amount) {
    uint16_t fract = static_cast<uint16_t>(amount * 65535.0f);
    a = from.a.lerp16(to.a, fract);
    b = from.b.lerp16(to.b, fract);
    c = from.c.lerp16(to.c, fract);
    update_luts();
  }

  /**
   * @brief Gets the color for a given position 't' by linearly interpolating between the defined points.
   * @param t Normalized value (0.0 to 1.0).
   * @return The calculated CRGB pixel color.
   */
  Pixel get(float t) const override {
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0) {
      seg = size - 2;
    }

    auto start = shape[seg];
    auto end = shape[seg + 1];
    Pixel c1 = colors[seg];
    Pixel c2 = colors[seg + 1];

    auto r = c1.lerp16(c2, std::clamp((t - start) / (end - start) * 65535.0f, 0.0f, 65535.0f));
    return r;
  }

private:

  /**
   * @brief Wraps a hue value to the valid 0-255 range.
   * @param hue The input hue value.
   * @return The wrapped hue (0-255).
   */
  uint8_t wrap_hue(int hue) const {
    return (hue % 256 + 256) % 256;
  }

  /**
   * @brief Calculates two secondary hues based on the primary hue and a harmony type.
   * @param h1 The primary hue.
   * @param h2 The first secondary hue (output).
   * @param h3 The second secondary hue (output).
   * @param harmony_type The rule to use (e.g., TRIADIC).
   */
  void calc_hues(uint8_t h1, uint8_t& h2, uint8_t& h3, HarmonyType harmony_type) const {
    const int h1_int = h1;

    switch (harmony_type) {
    case HarmonyType::TRIADIC: {
      h2 = wrap_hue(h1_int + 256 / 3);
      h3 = wrap_hue(h1_int + (256 * 2) / 3);
      break;
    }

    case HarmonyType::SPLIT_COMPLEMENTARY: {
      const int complement = wrap_hue(h1_int + 256 / 2);
      const int offset = 256 / 12;
      h2 = wrap_hue(complement - offset);
      h3 = wrap_hue(complement + offset);
      break;
    }

    case HarmonyType::COMPLEMENTARY: {
      h2 = wrap_hue(h1_int + 256 / 2);
      const int offset = hs::rand_int(-256 / 36, 256 / 36);
      h3 = wrap_hue(h1_int + offset);
      break;
    }

    case HarmonyType::ANALOGOUS:
    default: {
      const int dir = (hs::rand_int(0, 1) == 0) ? 1 : -1;
      const int offset1 = dir * hs::rand_int(256 / 8, 256 * 3 / 16);
      h2 = wrap_hue(h1_int + offset1);
      const int offset2 = dir * hs::rand_int(256 / 8, 256 * 3 / 16);
      h3 = wrap_hue(h2 + offset2);
      break;
    }
    }
  }

  GradientShape gradient_shape; /**< The spatial shape of the gradient. */
  HarmonyType harmony_type; /**< The rule used for hue selection. */
  uint8_t seed_hue; /**< The base hue used to generate the palette. */
  Pixel a, b, c; /**< The three primary colors of the palette (CHSV internally). */
  std::array<float, 5> shape; /**< Normalized position along the palette domain where colors change. */
  std::array<Pixel, 5> colors; /**< The corresponding colors at the shape points. */
  int size; /**< The number of active points in shape/colors. */
};

/**
 * @brief A palette defined by a mathematical cosine wave function for each RGB channel.
 * @details The function is: C(t) = A + B * cos(2 * PI * (C * t + D)).
 */
class ProceduralPalette : public Palette {
public:

  /**
   * @brief Constructs the procedural palette from its four defining vectors (A, B, C, D).
   * @param a The A vector (Base Value/Offset).
   * @param b The B vector (Amplitude).
   * @param c The C vector (Frequency).
   * @param d The D vector (Phase Shift).
   */
  ProceduralPalette(
    std::array<float, 3> a,
    std::array<float, 3> b,
    std::array<float, 3> c,
    std::array<float, 3> d) :
    a(a),
    b(b),
    c(c),
    d(d)
  {
  }

  /**
   * @brief Gets the color by calculating the cosine wave for each RGB channel at time 't'.
   * @param t Normalized value (0.0 to 1.0).
   * @return The calculated CRGB pixel color.
   */
  Pixel get(float t) const override {
    return Pixel(
      255 * (a[0] + b[0] * cosf(2 * PI_F * (c[0] * t + d[0]))),
      255 * (a[1] + b[1] * cosf(2 * PI_F * (c[1] * t + d[1]))),
      255 * (a[2] + b[2] * cosf(2 * PI_F * (c[2] * t + d[2])))
    );
  }

private:

  std::array<float, 3> a; /**< Base value / Offset. */
  std::array<float, 3> b; /**< Amplitude. */
  std::array<float, 3> c; /**< Frequency. */
  std::array<float, 3> d; /**< Phase shift. */
};

/**
 * @brief Type alias for the variant holding any supported palette type.
 */
using PaletteVariant =
std::variant<
  GenerativePalette,
  ProceduralPalette
>;

/**
 * @brief Creates a color function that applies a falloff/vignette effect to any palette.
 * @param palette The base palette.
 * @return A functor that returns a color dimmed at the edges (t < 0.2 and t > 0.8).
 */
auto vignette(const Palette& palette) {
  CRGB vignette_color(0, 0, 0);
  return [vignette_color, &palette](float t) {
    if (t < 0.2) {
      return vignette_color.lerp16(palette.get(0), to_short(t / 0.2));
    }
    else if (t >= 0.8) {
      return palette.get(1).lerp16(vignette_color, to_short((t - 0.8) / 0.2));
    }
    else {
      return palette.get((t - 0.2) / 0.6);
    }
    };
}

/**
 * @brief A pre-defined ProceduralPalette instance: Ice Melt.
 */
static const ProceduralPalette iceMelt(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.500, 0.500 }, // B
  { 0.083, 0.147, 0.082 }, // C
  { 0.579, 0.353, 0.244 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Embers.
 */
static const ProceduralPalette embers(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.500, 0.500 }, // B
  { 0.265, 0.285, 0.198 }, // C
  { 0.577, 0.440, 0.358 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Rich Sunset.
 */
static const ProceduralPalette richSunset(
  { 0.309, 0.500, 0.500 }, // A
  { 1.000, 1.000, 0.500 }, // B
  { 0.149, 0.148, 0.149 }, // C
  { 0.132, 0.222, 0.521 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Undersea.
 */
static const ProceduralPalette undersea(
  { 0.000, 0.000, 0.000 }, // A
  { 0.500, 0.276, 0.423 }, // B
  { 0.296, 0.296, 0.296 }, // C
  { 0.374, 0.941, 0.000 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Late Sunset.
 */
static const ProceduralPalette lateSunset(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.261, 0.261, 0.261 }, // C
  { 0.153, 0.483, 0.773 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Mango Peel.
 */
static const ProceduralPalette mangoPeel(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.080, 0.500 }, // B
  { 0.431, 0.431, 0.431 }, // C
  { 0.566, 0.896, 0.236 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Lemon Lime.
 */
static const ProceduralPalette lemonLime(
  { 0.455, 0.455, 0.455 }, // A
  { 0.571, 0.151, 0.571 }, // B
  { 0.320, 0.320, 0.320 }, // C
  { 0.087, 0.979, 0.319 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Algae.
 */
static const ProceduralPalette algae(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.134, 0.134, 0.134 }, // C
  { 0.328, 0.658, 0.948 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Fire Glow.
 */
static const ProceduralPalette fireGlow(
  { 0.000, 0.000, 0.000 }, // A
  { 0.560, 0.560, 0.560 }, // B
  { 0.216, 0.346, 0.174 }, // C
  { 0.756, 0.542, 0.279 }  // D
);

/**
 * @brief A pre-defined ProceduralPalette instance: Dark Primary.
 */
static const ProceduralPalette darkPrimary(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.610, 0.500 }, // B
  { 0.746, 0.347, 0.000 }, // C
  { 0.187, 0.417, 0.670 }  // D);
);

/**
 * @brief A pre-defined ProceduralPalette instance: Mauve Fade.
 */
static const ProceduralPalette mauveFade(
  { 0.583, 0.000, 0.583 }, // A
  { 1.000, 0.000, 1.000 }, // B
  { 0.191, 0.348, 0.191 }, // C
  { 0.175, 0.045, 0.150 }  // D
);