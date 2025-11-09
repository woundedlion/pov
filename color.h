#pragma once
///////////////////////////////////////////////////////////////////////////////

// gamma = 2.20 steps = 256 range = 0-1023
const uint16_t gamma_lut[256] = {
     0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   2,   2,
     2,   3,   3,   3,   4,   4,   5,   5,   6,   6,   7,   7,   8,   9,   9,  10,
    11,  11,  12,  13,  14,  15,  16,  16,  17,  18,  19,  20,  21,  23,  24,  25,
    26,  27,  28,  30,  31,  32,  34,  35,  36,  38,  39,  41,  42,  44,  46,  47,
    49,  51,  52,  54,  56,  58,  60,  61,  63,  65,  67,  69,  71,  73,  76,  78,
    80,  82,  84,  87,  89,  91,  94,  96,  98, 101, 103, 106, 109, 111, 114, 117,
   119, 122, 125, 128, 130, 133, 136, 139, 142, 145, 148, 151, 155, 158, 161, 164,
   167, 171, 174, 177, 181, 184, 188, 191, 195, 198, 202, 206, 209, 213, 217, 221,
   225, 228, 232, 236, 240, 244, 248, 252, 257, 261, 265, 269, 274, 278, 282, 287,
   291, 295, 300, 304, 309, 314, 318, 323, 328, 333, 337, 342, 347, 352, 357, 362,
   367, 372, 377, 382, 387, 393, 398, 403, 408, 414, 419, 425, 430, 436, 441, 447,
   452, 458, 464, 470, 475, 481, 487, 493, 499, 505, 511, 517, 523, 529, 535, 542,
   548, 554, 561, 567, 573, 580, 586, 593, 599, 606, 613, 619, 626, 633, 640, 647,
   653, 660, 667, 674, 681, 689, 696, 703, 710, 717, 725, 732, 739, 747, 754, 762,
   769, 777, 784, 792, 800, 807, 815, 823, 831, 839, 847, 855, 863, 871, 879, 887,
   895, 903, 912, 920, 928, 937, 945, 954, 962, 971, 979, 988, 997,1005,1014,1023,
};

///////////////////////////////////////////////////////////////////////////////

class Palette {
public:
  virtual Pixel get(double t) const = 0;
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


auto blend_alpha(double a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r * (1 - a), c2.r * a),
      qadd8(c1.g * (1 - a), c2.g * a),
      qadd8(c1.b * (1 - a), c2.b * a));
    };
}

uint16_t to_short(double zero_to_one) {
  return std::clamp(std::round(zero_to_one * 65535), 0.0, 65535.0);
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

  Pixel get(double t) const override {
    Serial.println(t);
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

    auto r = c1.lerp16(c2, std::clamp((t - start) / (end - start) * 65535, 0.0, 65535.0));
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
  std::array<double, 5> shape;
  std::array<Pixel, 5> colors;
  int size;
};

class ProceduralPalette : public Palette {
public:

  ProceduralPalette(
    std::array<double, 3> a,
    std::array<double, 3> b,
    std::array<double, 3> c,
    std::array<double, 3> d) :
    a(a),
    b(b),
    c(c),
    d(d)
  {
  }

  Pixel get(double t) const override {
    return Pixel(
      255 * (a[0] + b[0] * cos(2 * PI * (c[0] * t + d[0]))),
      255 * (a[1] + b[1] * cos(2 * PI * (c[1] * t + d[1]))),
      255 * (a[2] + b[2] * cos(2 * PI * (c[2] * t + d[2])))
    );
  }

private:

  std::array<double, 3> a;
  std::array<double, 3> b;
  std::array<double, 3> c;
  std::array<double, 3> d;
};

using PaletteVariant =
std::variant<
  GenerativePalette,
  ProceduralPalette
>;

auto vignette(const Palette& palette) {
  CRGB vignette_color(0, 0, 0);
  return [&](double t) {
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

