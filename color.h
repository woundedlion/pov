#pragma once

// TODO: 2D/3D Palettes

///////////////////////////////////////////////////////////////////////////////

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

Pixel gamma_correct(const Pixel& p) {
  return CRGB(
    gamma_lut[p.r],
    gamma_lut[p.g],
    gamma_lut[p.b]
  );
}
///////////////////////////////////////////////////////////////////////////////

class Palette {
public:
  virtual Pixel get(float_t t) const = 0;
};

enum class HarmonyType {
  TRIADIC,
  SPLIT_COMPLEMENTARY,
  COMPLEMENTARY,
  ANALOGOUS
};

enum class GradientShape {
  STRAIGHT,
  CIRCULAR,
  VIGNETTE,
  FALLOFF
};

auto blend_add(const Pixel& c1, const Pixel& c2) {
  return Pixel(
    qadd8(c1.r, c2.r),
    qadd8(c1.g, c2.g),
    qadd8(c1.b, c2.b));
}


auto blend_alpha(float_t a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r * (1 - a), c2.r * a),
      qadd8(c1.g * (1 - a), c2.g * a),
      qadd8(c1.b * (1 - a), c2.b * a));
    };
}

uint16_t to_short(float_t zero_to_one) {
  return std::clamp(std::round(zero_to_one * 65535), 0.0f, 65535.0f);
}

class GenerativePalette : public Palette {
public:

  GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type) :
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

    const uint8_t v1 = hs::rand_int(255 * 0.1, 255 * 0.1);
    const uint8_t v2 = hs::rand_int(255 * 0.2, 255 * 0.5);
    const uint8_t v3 = hs::rand_int(255 * 0.6, 255 * 0.8);

    a = CHSV(h1, s1, v1);
    b = CHSV(h2, s2, v2);
    c = CHSV(h3, s3, v3);

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

  Pixel get(float_t t) const override {
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

    auto r = c1.lerp16(c2, std::clamp((t - start) / (end - start) * 65535, 0.0f, 65535.0f));
    return r;
  }

private:

  uint8_t wrap_hue(int hue) {
    return (hue % 256 + 256) % 256;
  }

  void calc_hues(uint8_t h1, uint8_t& h2, uint8_t& h3, HarmonyType harmony_type) {
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
      const int offset1 = dir * hs::rand_int(256 / 12, 256 / 6);
      h2 = wrap_hue(h1_int + offset1);
      const int offset2 = dir * hs::rand_int(256 / 12, 256 / 6);
      h3 = wrap_hue(h2 + offset2);
      break;
    }
    }
  }

  GradientShape gradient_shape;
  HarmonyType harmony_type;
  uint8_t seed_hue;
  Pixel a, b, c;
  std::array<float_t, 5> shape;
  std::array<Pixel, 5> colors;
  int size;
};

class ProceduralPalette : public Palette {
public:

  ProceduralPalette(
    std::array<float_t, 3> a,
    std::array<float_t, 3> b,
    std::array<float_t, 3> c,
    std::array<float_t, 3> d) :
    a(a),
    b(b),
    c(c),
    d(d)
  {
  }

  Pixel get(float_t t) const override {
    return Pixel(
      255 * (a[0] + b[0] * cos(2 * PI * (c[0] * t + d[0]))),
      255 * (a[1] + b[1] * cos(2 * PI * (c[1] * t + d[1]))),
      255 * (a[2] + b[2] * cos(2 * PI * (c[2] * t + d[2])))
    );
  }

private:

  std::array<float_t, 3> a;
  std::array<float_t, 3> b;
  std::array<float_t, 3> c;
  std::array<float_t, 3> d;
};

using PaletteVariant =
std::variant<
  GenerativePalette,
  ProceduralPalette
>;

auto vignette(const Palette& palette) {
  CRGB vignette_color(0, 0, 0);
  return [&](float_t t) {
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


static const ProceduralPalette richSunset(
  { 0.309, 0.500, 0.500 }, // A
  { 1.000, 1.000, 0.500 }, // B
  { 0.149, 0.148, 0.149 }, // C
  { 0.132, 0.222, 0.521 }  // D
);

static const ProceduralPalette undersea(
  { 0.000, 0.000, 0.000 }, // A
  { 0.500, 0.276, 0.423 }, // B
  { 0.296, 0.296, 0.296 }, // C
  { 0.374, 0.941, 0.000 }  // D
);

static const ProceduralPalette lateSunset(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.261, 0.261, 0.261 }, // C
  { 0.153, 0.483, 0.773 }  // D
);

static const ProceduralPalette mangoPeel(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.080, 0.500 }, // B
  { 0.431, 0.431, 0.431 }, // C
  { 0.566, 0.896, 0.236 }  // D
);

static const ProceduralPalette lemonLime(
  { 0.455, 0.455, 0.455 }, // A
  { 0.571, 0.151, 0.571 }, // B
  { 0.320, 0.320, 0.320 }, // C
  { 0.087, 0.979, 0.319 }  // D
);

static const ProceduralPalette algae(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.134, 0.134, 0.134 }, // C
  { 0.328, 0.658, 0.948 }  // D
);

