/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Tolerance pin for core/color/triadic_palette_luts.h, the 1.6 MB palette bank
// tools/mindsplatter_palette_gen.cpp emits. Recompiles every entry from the
// recipe of record (EffectPaletteRecipes::mind_splatter) and compares it to the
// committed one, so a recipe edit, a GenerativePalette change or a hand-edit of
// the table fails here instead of shipping a stale bank.
//
// The sibling generated artifacts are byte-diffed against a regenerated copy;
// this one cannot be. The recipe path runs powf/cbrtf through the OKLab gamut
// search, whose last bits differ between libm builds, so a byte-diff would red
// on whichever host did not bake the committed header. Entries are held to a
// per-channel tolerance instead, and the worst observed delta is printed so a
// drift that stays inside it is still visible.
#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "core/color/effect_palette_recipes.h"
#include "core/color/triadic_palette_luts.h"
#include "core/engine/engine.h"

namespace {

// Per-channel slack in 16-bit linear units. Measured worst delta over the bank
// is 1 at -O0 and at -O2 -ffast-math -fno-finite-math-only, so the band is four
// times the observed libm spread — enough for a cbrtf/powf whose last bits
// differ, and short of the 31 units a 0.03% shift in one recipe field reaches.
constexpr int MAX_CHANNEL_DELTA = 4;

// Mean absolute per-channel delta over the whole bank, in thousandths of a
// linear unit. The per-entry band above is blind to a drift that moves every
// channel by less than itself, which the form gate cannot see either: it
// rewrites Pixel(N, N, N) to Pixel(VALUE) on both sides. Libm noise is sparse:
// the measured mean is 0/1000, fewer than one channel in a thousand differing
// at all, so a shift that moves the whole bank cannot hide under it.
constexpr long long MAX_MEAN_DELTA_MILLI = 30;

/** @brief Per-channel differences between two bank entries. */
struct ChannelDelta {
  int worst; /**< Largest single-channel difference. */
  int total; /**< Sum over the three channels. */
};

ChannelDelta channel_delta(const Pixel &a, const Pixel &b) {
  const int dr = std::abs(static_cast<int>(a.r) - static_cast<int>(b.r));
  const int dg = std::abs(static_cast<int>(a.g) - static_cast<int>(b.g));
  const int db = std::abs(static_cast<int>(a.b) - static_cast<int>(b.b));
  return {std::max({dr, dg, db}), dr + dg + db};
}

} // namespace

int main() {
  int worst = 0;
  int worst_palette = 0;
  int worst_index = 0;
  long long delta_sum = 0;
  long long channels = 0;
  for (int hue = 0; hue < MINDSPLATTER_PALETTE_COUNT; ++hue) {
    const PaletteRecipe recipe = EffectPaletteRecipes::mind_splatter(
        EffectPaletteRecipes::mind_splatter_bank_turns(
            hue, MINDSPLATTER_PALETTE_COUNT));
    GenerativePalette palette;
    PaletteRecipe canonical;
    PaletteCompileStatus status;
    if (!GenerativePalette::try_compile(recipe, palette, canonical, status)) {
      std::printf("mindsplatter palette bank: palette %d no longer compiles "
                  "(recipe error %d at field %d)\n",
                  hue, static_cast<int>(status.code),
                  static_cast<int>(status.field));
      return 1;
    }
    for (int index = 0; index < MINDSPLATTER_PALETTE_LUT_SIZE; ++index) {
      const Pixel derived = palette
                                .get(static_cast<float>(index) /
                                     (MINDSPLATTER_PALETTE_LUT_SIZE - 1))
                                .color;
      const ChannelDelta delta =
          channel_delta(derived, MINDSPLATTER_PALETTES[hue][index]);
      delta_sum += delta.total;
      channels += 3;
      if (delta.worst > worst) {
        worst = delta.worst;
        worst_palette = hue;
        worst_index = index;
      }
    }
  }

  const long long mean_milli = channels ? (delta_sum * 1000) / channels : 0;
  if (worst > MAX_CHANNEL_DELTA) {
    std::printf(
        "mindsplatter palette bank: entry [%d][%d] is off by %d "
        "channel units (tolerance %d) — core/color/triadic_palette_luts.h "
        "no longer matches EffectPaletteRecipes::mind_splatter; "
        "regenerate with: cmake --build --preset tests --target "
        "regenerate_mindsplatter_palette\n",
        worst_palette, worst_index, worst, MAX_CHANNEL_DELTA);
    return 1;
  }
  if (mean_milli > MAX_MEAN_DELTA_MILLI) {
    std::printf(
        "mindsplatter palette bank: mean channel delta %lld/1000 units "
        "exceeds %lld/1000 — the bank has drifted as a whole, which the "
        "per-entry band and the form gate both let through; regenerate with: "
        "cmake --build --preset tests --target "
        "regenerate_mindsplatter_palette\n",
        mean_milli, MAX_MEAN_DELTA_MILLI);
    return 1;
  }
  std::printf("mindsplatter palette bank: %d entries re-derived, worst channel "
              "delta %d of %d, mean %lld/1000 of %lld/1000\n",
              MINDSPLATTER_PALETTE_COUNT * MINDSPLATTER_PALETTE_LUT_SIZE, worst,
              MAX_CHANNEL_DELTA, mean_milli, MAX_MEAN_DELTA_MILLI);
  return 0;
}
