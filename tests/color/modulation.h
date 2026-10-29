/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file modulation.h
 * @brief Generic mechanics, identities and palette composition tests. */

/** @brief Mechanics preserve exact identities and straight-alpha masks. */
inline void test_color_modulation_mechanics() {
  const Color4 INPUT(Pixel(4000, 8000, 12000), .37f);
  for (float scalar : {-1.0f, 0.0f, 1.0f}) {
    for (auto kind : {ColorModulation::Mechanic::HUE_ROTATION,
                      ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                      ColorModulation::Mechanic::OPACITY}) {
      const auto result = ColorModulation::apply(INPUT, scalar, {kind, 0.0f});
      HS_EXPECT_EQ(result.color, INPUT.color);
      HS_EXPECT_EQ(result.alpha, INPUT.alpha);
    }
    const auto masked = ColorModulation::Opacity::apply(INPUT, scalar, {1});
    HS_EXPECT_EQ(masked.color, INPUT.color);
    HS_EXPECT_EQ(masked.alpha, INPUT.alpha * .5f * (scalar + 1));
    HS_EXPECT_EQ(ColorModulation::HueRotation::apply(INPUT, scalar, {1}).color,
                 INPUT.color);
  }
  HS_EXPECT_EQ(
      ColorModulation::PositiveLightnessLift::apply(INPUT, -1, {1}).color,
      INPUT.color);
  HS_EXPECT_EQ(
      ColorModulation::PositiveLightnessLift::apply(INPUT, 0, {1}).color,
      INPUT.color);
  const auto raised =
      ColorModulation::PositiveLightnessLift::apply(INPUT, .5f, {.8f});
  const auto rgb = pixel_to_linrgb(INPUT.color);
  const auto base = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
  const auto output = pixel_to_linrgb(raised.color);
  const auto lab = linear_rgb_to_oklab_fast(output.r, output.g, output.b);
  HS_EXPECT_NEAR(lab.L, base.L + (1 - base.L) * .4f, .005f);
  HS_EXPECT_EQ(raised.alpha, INPUT.alpha);
}
struct PatternTestCoordinate {
  float scalar;
};
/** @brief Generic composition accepts a custom typed field and mechanic. */
inline void test_pattern_modulated_palette() {

  struct Field {
    using Coordinate = ::hs_test::color_tests::PatternTestCoordinate;
    struct Prepared {
      int *calls;
    };
    static float sample(const Coordinate &c, const Prepared &p) {
      ++*p.calls;
      return c.scalar;
    }
  };
  struct Mechanic {
    struct Prepared {};
    static Color4 apply(Color4 c, float scalar, const Prepared &) {
      c.alpha *= scalar;
      return c;
    }
  };
  SolidColorPalette source(Color4(Pixel(1234, 5678, 9012), .4f));
  int calls = 0;
  ColorModulation::PatternModulatedPalette<SolidColorPalette, Field, Mechanic>
      custom{&source, {&calls}, {}};
  auto result = custom.get(.21f, {.25f});
  HS_EXPECT_EQ(result.color, source.get(.21f).color);
  HS_EXPECT_EQ(result.alpha, .1f);
  HS_EXPECT_EQ(calls, 1);
  ColorModulation::PatternModulatedPalette<SolidColorPalette, Field,
                                           ColorModulation::Opacity>
      mask{&source, {&calls}, {0}};
  result = mask.get(.213f, {-1});
  HS_EXPECT_EQ(result.color, source.get(.213f).color);
  HS_EXPECT_EQ(result.alpha, .4f);
  HS_EXPECT_EQ(calls, 1);
}

inline void test_hue_rotation_lut_clamps_out_of_range_value() {
  static std::array<Pixel, 3 * HueRotationLutView::SIZE> storage;
  storage.fill(Pixel(0, 0, 65535));
  Pixel *const lut = storage.data() + HueRotationLutView::SIZE;
  for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
    lut[hue] = Pixel(65535, 0, 0);
    lut[(HueRotationLutView::VALUE_STEPS - 1) * HueRotationLutView::HUE_STEPS +
        hue] = Pixel(0, 65535, 0);
  }

  const HueRotationLutView view{lut, true};
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, -0.5f, 0.25f), Pixel(65535, 0, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, 1.5f, 0.25f), Pixel(0, 65535, 0));
  HS_EXPECT_EQ(sample_hue_rotation_lut(view, NAN, 0.25f), Pixel(0, 65535, 0));
  for (float amount : {-4.75f, -1e-7f, 3.5e8f, -3.5e8f, INFINITY, NAN})
    HS_EXPECT_EQ(sample_hue_rotation_lut(view, 0.0f, amount),
                 Pixel(65535, 0, 0));
}
