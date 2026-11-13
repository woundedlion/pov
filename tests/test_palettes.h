/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/color/palettes.h: named procedural palettes and the mesh palette bank.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "core/color/color.h"
#include "core/color/palettes.h"
#include "core/platform/rng.h"
#include "tests/color_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace palettes_tests {

/** @brief Integer-frequency DARK_RAINBOW has matching endpoint colors. */
inline void test_named_procedural_palette_loop_endpoints() {
  const Color4 a = Palettes::DARK_RAINBOW.get(0.0f);
  const Color4 b = Palettes::DARK_RAINBOW.get(1.0f);
  HS_EXPECT_EQ(a.color.r, b.color.r);
  HS_EXPECT_EQ(a.color.g, b.color.g);
  HS_EXPECT_EQ(a.color.b, b.color.b);
}

/**
 * @brief Bakes each palette slot from its source and keeps the slots distinct.
 */
inline void test_mesh_palette_bank_lookup() {
  alignas(std::max_align_t) static uint8_t
      buf[MeshPaletteBank::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  MeshPaletteBank bank;
  bank.bake_all(arena);

  for (int i = 0; i < MeshPaletteBank::N; ++i) {
    for (float t : {0.0f, 1.0f}) {
      const Color4 actual = bank[i].get(t);
      const Color4 expected = MeshPaletteBank::sources()[i]->get(t);
      HS_EXPECT_EQ(actual.color.r, expected.color.r);
      HS_EXPECT_EQ(actual.color.g, expected.color.g);
      HS_EXPECT_EQ(actual.color.b, expected.color.b);
    }
  }

  // Every slot bakes a distinct LUT (no two share a t=0 color).
  for (int i = 0; i < MeshPaletteBank::N; ++i)
    for (int j = i + 1; j < MeshPaletteBank::N; ++j) {
      Color4 ci = bank[i].get(0.0f), cj = bank[j].get(0.0f);
      HS_EXPECT_TRUE(ci.color.r != cj.color.r || ci.color.g != cj.color.g ||
                     ci.color.b != cj.color.b);
    }
}

/**
 * @brief Verifies shuffle_indices yields a permutation of [0, N).
 * @details The global generator is saved and restored around the shuffle.
 */
inline void test_mesh_palette_bank_shuffle_is_permutation() {
  const hs_test::ScopedRng saved;
  hs::random().seed(1337);
  std::array<int, MeshPaletteBank::N> idx{};
  MeshPaletteBank::shuffle_indices(idx);
  std::array<int, MeshPaletteBank::N> seen{};
  for (int v : idx) {
    HS_EXPECT_GE(v, 0);
    HS_EXPECT_LT(v, MeshPaletteBank::N);
    if (v >= 0 && v < MeshPaletteBank::N)
      seen[v]++;
  }
  bool changed = false;
  for (size_t i = 0; i < idx.size(); ++i)
    changed |= idx[i] != static_cast<int>(i);
  HS_EXPECT_TRUE(changed);
  for (int count : seen)
    HS_EXPECT_EQ(count, 1);
}

/** @brief Expects each baked LUT entry within one code of its source sample. */
inline void expect_bake_matches_source(const Palette &palette) {
  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage baked;
  baked.bake(arena, palette);
  for (int index = 0; index < BakedPalette::LUT_SIZE; ++index) {
    const float t = static_cast<float>(index) / (BakedPalette::LUT_SIZE - 1);
    const Pixel source = palette.get(t).color;
    const Pixel lut = baked.get(t).color;
    HS_EXPECT_NEAR(lut.r, source.r, 1);
    HS_EXPECT_NEAR(lut.g, source.g, 1);
    HS_EXPECT_NEAR(lut.b, source.b, 1);
  }
}

/**
 * @brief Procedural palettes evaluate the clamped cosine formula in sRGB and
 *        bake to within one code of their samples.
 */
inline void test_procedural_palette_matches_cosine_oracle() {
  struct Coefficients {
    std::array<float, 3> a, b, c, d;
  };
  constexpr Coefficients CASES[] = {
      // Clamp: channels overshoot 1 and undershoot 0.
      {{0.9f, 0.1f, 0.5f}, {0.5f, 0.5f, 0.8f}, {1, 1, 1}, {0, 0.25f, 0.6f}},
      // Phase wrap: phases outside [0, 1) turn.
      {{0.5f, 0.5f, 0.5f},
       {0.4f, 0.3f, 0.45f},
       {1, 1, 1},
       {1.25f, -0.75f, 2.4f}},
      // Negative amplitude.
      {{0.5f, 0.4f, 0.6f},
       {-0.5f, -0.3f, -0.45f},
       {1, 0.5f, 2},
       {0, 0.1f, 0.3f}},
      // High frequency.
      {{0.5f, 0.5f, 0.5f},
       {0.5f, 0.4f, 0.3f},
       {7, 11.5f, 16},
       {0.1f, 0.2f, 0.3f}},
  };
  // Bhaskara cosine error, plus argument rounding at high frequency.
  constexpr double COS_ERROR = 1.7e-3 + 1e-4;
  // srgb_to_linear_interp's LUT and secant slack, in codes.
  constexpr double LUT_SLACK = 2.0;
  constexpr double TWO_PI = 6.283185307179586;
  const auto linear = [](double s) {
    return 65535.0 *
           srgb_to_linear_float(static_cast<float>(std::clamp(s, 0.0, 1.0)));
  };
  for (size_t index = 0; index < std::size(CASES); ++index) {
    HS_CONTEXT("case", index);
    const Coefficients &k = CASES[index];
    const ProceduralPalette palette(k.a, k.b, k.c, k.d);
    for (int sample = 0; sample <= 256; ++sample) {
      const float t = sample / 256.0f;
      const Color4 actual = palette.get(t);
      HS_EXPECT_EQ(actual.alpha, 1.0f);
      const uint16_t channels[] = {actual.color.r, actual.color.g,
                                   actual.color.b};
      for (int ch = 0; ch < 3; ++ch) {
        const double exact =
            k.a[ch] +
            k.b[ch] * std::cos(TWO_PI * (double(k.c[ch]) * t + k.d[ch]));
        const double spread = std::fabs(k.b[ch]) * COS_ERROR;
        if (exact - spread >= 1.0)
          HS_EXPECT_EQ(channels[ch], 65535);
        else if (exact + spread <= 0.0)
          HS_EXPECT_EQ(channels[ch], 0);
        HS_EXPECT_GE(channels[ch], linear(exact - spread) - LUT_SLACK);
        HS_EXPECT_LE(channels[ch], linear(exact + spread) + LUT_SLACK);
      }
    }
    expect_bake_matches_source(palette);
  }
}

/** @brief Every shipping procedural palette is opaque, lit, and bakes faithfully. */
inline void test_named_procedural_palettes_bake() {
  const Palette *palettes[] = {
#define HS_PALETTE_ENTRY(name, A, B, C, D) &Palettes::name,
      HS_PROCEDURAL_PALETTE_LIST(HS_PALETTE_ENTRY)
#undef HS_PALETTE_ENTRY
  };
  for (size_t index = 0; index < std::size(palettes); ++index) {
    HS_CONTEXT("palette index", index);
    uint16_t brightest = 0;
    for (int sample = 0; sample <= 16; ++sample) {
      const Color4 color = palettes[index]->get(sample / 16.0f);
      HS_EXPECT_EQ(color.alpha, 1.0f);
      brightest =
          std::max({brightest, color.color.r, color.color.g, color.color.b});
    }
    HS_EXPECT_GT(brightest, 0);
    expect_bake_matches_source(*palettes[index]);
  }
}

/**
 * @brief Runs every palettes-module test and reports the aggregate result.
 * @return 0 on success, non-zero on any failure.
 */
inline int run_palettes_tests() {
  hs_test::ModuleFixture fixture("palettes");

  HS_RUN(test_procedural_palette_matches_cosine_oracle);
  HS_RUN(test_named_procedural_palettes_bake);
  HS_RUN(test_named_procedural_palette_loop_endpoints);
  HS_RUN(test_mesh_palette_bank_lookup);
  HS_RUN(test_mesh_palette_bank_shuffle_is_permutation);

  return fixture.result();
}

} // namespace palettes_tests
} // namespace hs_test
