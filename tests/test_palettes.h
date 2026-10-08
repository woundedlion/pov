/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/color/palettes.h — the named generative/OKLCH palette layer.
 */
#pragma once

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

/**
 * @brief Pins named ProceduralPalette endpoints against golden 16-bit colors.
 */
inline void test_named_procedural_palette_endpoints() {
  // darkRainbow has integer frequencies, so t=0 and t=1 give the same color.
  Color4 dr0 = Palettes::DARK_RAINBOW.get(0.0f);
  Color4 dr1 = Palettes::DARK_RAINBOW.get(1.0f);
  HS_EXPECT_EQ(dr0.color.r, 47426);
  HS_EXPECT_EQ(dr0.color.g, 954);
  HS_EXPECT_EQ(dr0.color.b, 954);
  HS_EXPECT_EQ(dr1.color.r, dr0.color.r);
  HS_EXPECT_EQ(dr1.color.g, dr0.color.g);
  HS_EXPECT_EQ(dr1.color.b, dr0.color.b);
  HS_EXPECT_NEAR(dr0.alpha, 1.0f, 1e-6f);

  // mauveFade: red and blue clamp to 1 at t=0.
  Color4 mf0 = Palettes::MAUVE_FADE.get(0.0f);
  HS_EXPECT_EQ(mf0.color.r, 65535);
  HS_EXPECT_EQ(mf0.color.g, 0);
  HS_EXPECT_EQ(mf0.color.b, 65535);

  // fireGlow / peachPop: fully golden-pinned (fast_cosf channels).
  Color4 fg0 = Palettes::FIRE_GLOW.get(0.0f);
  Color4 fg1 = Palettes::FIRE_GLOW.get(1.0f);
  HS_EXPECT_EQ(fg0.color.r, 108);
  HS_EXPECT_EQ(fg0.color.g, 0);
  HS_EXPECT_EQ(fg0.color.b, 0);
  HS_EXPECT_EQ(fg1.color.r, 17340);
  HS_EXPECT_EQ(fg1.color.g, 9961);
  HS_EXPECT_EQ(fg1.color.b, 0);

  Color4 pp0 = Palettes::PEACH_POP.get(0.0f);
  HS_EXPECT_EQ(pp0.color.r, 65535);
  HS_EXPECT_EQ(pp0.color.g, 28156);
  HS_EXPECT_EQ(pp0.color.b, 0);

  Color4 cb0 = Palettes::CORAL_BLUE.get(0.0f);
  Color4 cb1 = Palettes::CORAL_BLUE.get(1.0f);
  HS_EXPECT_EQ(cb0.color.r, 42854);
  HS_EXPECT_EQ(cb0.color.g, 13737);
  HS_EXPECT_EQ(cb0.color.b, 8572);
  HS_EXPECT_EQ(cb1.color.r, 1343);
  HS_EXPECT_EQ(cb1.color.g, 1158);
  HS_EXPECT_EQ(cb1.color.b, 7717);
  HS_EXPECT_NEAR(cb0.alpha, 1.0f, 1e-6f);

  constexpr auto mesh_sources = MeshPaletteBank::sources();
  HS_EXPECT_TRUE(mesh_sources[5] == &Palettes::TIDAL_JADE);

  constexpr ProceduralPalette FORWARD_ORANGE(
      {0.575f, 0.168f, 0.464f}, {0.406f, 0.697f, 0.357f},
      {0.000f, 0.10051f, 0.042778f}, {0.141f, 0.155f, 0.537f});
  for (float t : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
    Color4 reversed = Palettes::ORANGE_CRUSH.get(t);
    Color4 forward = FORWARD_ORANGE.get(1.0f - t);
    HS_EXPECT_NEAR(reversed.color.r, forward.color.r, 1);
    HS_EXPECT_NEAR(reversed.color.g, forward.color.g, 1);
    HS_EXPECT_NEAR(reversed.color.b, forward.color.b, 1);
  }
}

/**
 * @brief Verifies hue interpolates along the short arc for a seam-straddling
 *        named-palette pair.
 * @details undersea's endpoints straddle the +/-PI seam; the midpoint must
 *          cross the seam, not land near the naive average.
 */
inline void test_named_palette_hue_short_arc() {
  OKLCH a = pixel_to_oklch(Palettes::UNDERSEA.get(0.0f).color);
  OKLCH b = pixel_to_oklch(Palettes::UNDERSEA.get(1.0f).color);
  // Precondition: endpoints really do straddle the seam (numeric gap > PI).
  HS_EXPECT_GT(std::fabs(b.h - a.h), math::PI_F);

  OKLCH mid = lerp_oklch(a, b, 0.5f);
  float short_mid = a.h + 0.5f * wrap_hue_delta(b.h - a.h);
  HS_EXPECT_NEAR(wrap_hue_delta(mid.h - short_mid), 0.0f, 1e-3f);
  float naive_mid = 0.5f * (a.h + b.h);
  HS_EXPECT_GT(std::fabs(wrap_hue_delta(mid.h - naive_mid)), 1.0f);
}

/**
 * @brief Verifies MeshPaletteBank bakes its sources and indexes them distinctly.
 * @details Slot 0's baked endpoints match the first source, and every slot
 *          produces a distinct LUT.
 */
inline void test_mesh_palette_bank_lookup() {
  alignas(std::max_align_t) static uint8_t
      buf[MeshPaletteBank::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  MeshPaletteBank bank;
  bank.bake_all(arena);

  // Slot 0 == embers baked: pinned golden endpoints.
  Color4 s0 = bank[0].get(0.0f);
  Color4 s1 = bank[0].get(1.0f);
  HS_EXPECT_EQ(s0.color.r, 307);
  HS_EXPECT_EQ(s0.color.g, 180);
  HS_EXPECT_EQ(s0.color.b, 1906);
  HS_EXPECT_EQ(s1.color.r, 36642);
  HS_EXPECT_EQ(s1.color.g, 9703);
  HS_EXPECT_EQ(s1.color.b, 157);
  // The baked slot reproduces the embers source it was baked from.
  Color4 e0 = Palettes::EMBERS.get(0.0f);
  HS_EXPECT_EQ(s0.color.r, e0.color.r);
  HS_EXPECT_EQ(s0.color.g, e0.color.g);
  HS_EXPECT_EQ(s0.color.b, e0.color.b);

  for (int i = 0; i < MeshPaletteBank::N; ++i) {
    for (float t : {0.0f, 1.0f}) {
      const Color4 ACTUAL = bank[i].get(t);
      const Color4 EXPECTED = MeshPaletteBank::sources()[i]->get(t);
      HS_EXPECT_EQ(ACTUAL.color.r, EXPECTED.color.r);
      HS_EXPECT_EQ(ACTUAL.color.g, EXPECTED.color.g);
      HS_EXPECT_EQ(ACTUAL.color.b, EXPECTED.color.b);
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
  auto saved = hs::random();
  hs::random().seed(1337);
  std::array<int, MeshPaletteBank::N> idx{};
  MeshPaletteBank::shuffle_indices(idx);
  hs::random() = saved;
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

/**
 * @brief Pins upper-byte color samples across the procedural roster.
 * @details Captured by printing this loop's hashes under the native clang
 * toolchain (cmake/toolchain-native-clang.cmake). Re-derive the roster and
 * per-palette hashes the same way after an intentional palette retune.
 */
inline void test_named_procedural_palette_roster() {
  const Palette *palettes[] = {
#define HS_PALETTE_ENTRY(name, A, B, C, D) &Palettes::name,
      HS_PROCEDURAL_PALETTE_LIST(HS_PALETTE_ENTRY)
#undef HS_PALETTE_ENTRY
  };
  constexpr uint64_t PALETTE_HASHES[] = {
      2886858739051273591ull,  12000220049008470606ull,
      13336441504933032662ull, 7882225716504084908ull,
      18404027224185952013ull, 11759190715424212606ull,
      12081953290500887235ull, 6518849557875476664ull,
      8968209609131571649ull,  15181736134788936411ull,
      10913766975351473374ull, 13329204743012588901ull,
      14302435841029644024ull, 11989053834225940204ull,
      17532158822806458534ull, 13333059190711347450ull,
      1983164980008296279ull,  17191927165933662661ull,
      11014468141546408976ull, 15566500396385584882ull,
      7646521029703227679ull,  10063136414861090927ull,
      11972795300560289134ull, 206858551135936371ull,
      7425003234454588729ull,  51388825653918428ull,
      15037618790413056646ull};
  static_assert(std::size(palettes) == std::size(PALETTE_HASHES));
  uint64_t hash = FNV1A64_BASIS;
  for (size_t index = 0; index < std::size(palettes); ++index) {
    HS_CONTEXT("palette index", index);
    uint64_t palette_hash = FNV1A64_BASIS;
    for (int sample = 0; sample <= 16; ++sample) {
      const Pixel pixel = palettes[index]->get(sample / 16.0f).color;
      for (uint16_t channel : {pixel.r, pixel.g, pixel.b}) {
        const auto upper = static_cast<uint8_t>(channel >> 8);
        hash = fnv1a64_byte(hash, upper);
        palette_hash = fnv1a64_byte(palette_hash, upper);
      }
    }
    HS_EXPECT_EQ(palette_hash, PALETTE_HASHES[index]);
  }
  HS_EXPECT_EQ(hash, uint64_t{1253289522805892800});
}

/**
 * @brief Runs every palettes-module test and reports the aggregate result.
 * @return 0 on success, non-zero on any failure.
 */
inline int run_palettes_tests() {
  hs_test::ModuleFixture fixture("palettes");

  test_named_procedural_palette_roster();
  test_named_procedural_palette_endpoints();
  test_named_palette_hue_short_arc();
  test_mesh_palette_bank_lookup();
  test_mesh_palette_bank_shuffle_is_permutation();

  return fixture.result();
}

} // namespace palettes_tests
} // namespace hs_test
