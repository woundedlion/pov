/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>

#include "color/color.h"

namespace Palettes {

// Procedural Palettes
inline constexpr ProceduralPalette DARK_RAINBOW({0.367f, 0.367f, 0.367f},
                                               {0.500f, 0.500f, 0.500f},
                                               {1.000f, 1.000f, 1.000f},
                                               {0.000f, 0.330f, 0.670f});
inline constexpr ProceduralPalette BLOOD_STREAM({0.169f, 0.169f, 0.169f},
                                               {0.313f, 0.313f, 0.313f},
                                               {0.231f, 0.231f, 0.231f},
                                               {0.036f, 0.366f, 0.706f});
inline constexpr ProceduralPalette VINTAGE_SUNSET({0.256f, 0.256f, 0.256f},
                                                 {0.500f, 0.080f, 0.500f},
                                                 {0.277f, 0.277f, 0.277f},
                                                 {0.000f, 0.330f, 0.670f});
inline constexpr ProceduralPalette RICH_SUNSET({0.309f, 0.500f, 0.500f},
                                              {1.000f, 1.000f, 0.500f},
                                              {0.149f, 0.148f, 0.149f},
                                              {0.132f, 0.222f, 0.521f});
inline constexpr ProceduralPalette UNDERSEA({0.000f, 0.000f, 0.000f},
                                            {0.500f, 0.276f, 0.423f},
                                            {0.296f, 0.296f, 0.296f},
                                            {0.374f, 0.941f, 0.000f});
inline constexpr ProceduralPalette LATE_SUNSET({0.337f, 0.500f, 0.096f},
                                              {0.500f, 1.000f, 0.176f},
                                              {0.261f, 0.261f, 0.261f},
                                              {0.153f, 0.483f, 0.773f});
inline constexpr ProceduralPalette MANGO_PEEL({0.500f, 0.500f, 0.500f},
                                             {0.500f, 0.080f, 0.500f},
                                             {0.431f, 0.431f, 0.431f},
                                             {0.566f, 0.896f, 0.236f});
inline constexpr ProceduralPalette ICE_MELT({0.500f, 0.500f, 0.500f},
                                           {0.500f, 0.500f, 0.500f},
                                           {0.083f, 0.147f, 0.082f},
                                           {0.579f, 0.353f, 0.244f});
inline constexpr ProceduralPalette LEMON_LIME({0.455f, 0.455f, 0.455f},
                                             {0.571f, 0.151f, 0.571f},
                                             {0.320f, 0.320f, 0.320f},
                                             {0.087f, 0.979f, 0.319f});
inline constexpr ProceduralPalette ALGAE({0.210f, 0.210f, 0.210f},
                                         {0.500f, 1.000f, 0.021f},
                                         {0.086f, 0.086f, 0.075f},
                                         {0.419f, 0.213f, 0.436f});
inline constexpr ProceduralPalette EMBERS({0.500f, 0.500f, 0.500f},
                                          {0.500f, 0.500f, 0.500f},
                                          {0.265f, 0.285f, 0.198f},
                                          {0.577f, 0.440f, 0.358f});
inline constexpr ProceduralPalette FIRE_GLOW({0.000f, 0.000f, 0.000f},
                                            {0.560f, 0.560f, 0.560f},
                                            {0.216f, 0.346f, 0.174f},
                                            {0.756f, 0.542f, 0.279f});
inline constexpr ProceduralPalette DARK_PRIMARY({0.500f, 0.500f, 0.500f},
                                               {0.500f, 0.610f, 0.500f},
                                               {0.746f, 0.347f, 0.000f},
                                               {0.187f, 0.417f, 0.670f});
inline constexpr ProceduralPalette MAUVE_FADE({0.583f, 0.000f, 0.583f},
                                             {1.000f, 0.000f, 1.000f},
                                             {0.191f, 0.348f, 0.191f},
                                             {0.175f, 0.045f, 0.150f});
inline constexpr ProceduralPalette LAVENDER_LAKE({0.473f, 0.473f, 0.473f},
                                                {0.500f, 0.500f, 0.500f},
                                                {0.364f, 0.124f, 0.528f},
                                                {0.142f, 0.378f, 0.876f});
inline constexpr ProceduralPalette DESERT_ROSE({0.500f, 0.500f, 0.500f},
                                              {0.500f, 0.270f, 0.442f},
                                              {0.303f, 1.012f, 0.585f},
                                              {0.985f, 0.720f, 0.212f});
inline constexpr ProceduralPalette BRUISED_MOSS({0.500f, 0.500f, 0.500f},
                                               {0.500f, 0.500f, 0.500f},
                                               {0.142f, 0.252f, 0.000f},
                                               {0.492f, 0.200f, 0.670f});
inline constexpr ProceduralPalette BRUISED_BANANA({0.620f, 0.620f, 0.620f},
                                                 {0.742f, 0.742f, 0.742f},
                                                 {0.162f, 0.286f, 0.012f},
                                                 {0.235f, 0.205f, 0.688f});
inline constexpr ProceduralPalette BRIGHT_SUNRISE({0.620f, 0.620f, 0.620f},
                                                 {0.742f, 0.742f, 0.742f},
                                                 {0.162f, 0.286f, 0.012f},
                                                 {0.090f, 0.205f, 0.688f});
inline constexpr ProceduralPalette FIRE_AND_ICE({0.500f, 0.500f, 0.500f},
                                              {0.500f, 0.500f, 0.500f},
                                              {0.955f, 1.004f, 0.910f},
                                              {0.167f, 0.018f, 0.930f});
inline constexpr ProceduralPalette PEACH_POP({1.000f, 0.144f, 0.175f},
                                            {0.543f, 0.543f, 0.543f},
                                            {0.507f, 0.409f, 0.507f},
                                            {0.001f, 0.002f, 0.620f});

} // namespace Palettes

/**
 * @brief Shared 5-palette "mesh effect" bank used by HankinSolids / IslamicStars
 *        and any future mesh effect.
 * @details Bundles the standard source-palette set, the bake-all step, and the
 *          per-shape index shuffle these effects share. Zero-overhead: the
 *          source list is constexpr and every accessor is a thin inline wrapper
 *          over BakedPaletteBank, so the per-pixel lookup remains
 *          BakedPalette::get() with no added indirection.
 */
struct MeshPaletteBank {
  static constexpr int N = BakedPaletteBank::N; // 5

  /**
   * @brief Arena bytes bake_all() consumes, including worst-case per-palette
   *        alignment padding.
   */
  static constexpr size_t required_arena_bytes() {
    return N * BakedPalette::required_arena_bytes();
  }

  /**
   * @brief Standard source palettes, in slot order.
   * @return Array of pointers to the constexpr source palettes. The size is
   *         deduced from the list (CTAD) — not fixed to N — so a list that
   *         drifts out of sync with N trips the static_assert below instead of
   *         silently nullptr-padding or reading past the bank.
   */
  static constexpr auto sources() {
    return std::array{&Palettes::EMBERS, &Palettes::RICH_SUNSET,
                      &Palettes::BRIGHT_SUNRISE, &Palettes::BRUISED_MOSS,
                      &Palettes::LAVENDER_LAKE};
  }

  /**
   * @brief (Re)bakes every source palette into the arena.
   * @param arena Destination arena receiving N x 256-entry Color4 LUTs.
   */
  void bake_all(Arena &arena) {
    constexpr auto src = sources();
    for (int i = 0; i < N; ++i)
      bank.entries[i].bake(arena, *src[i]);
  }

  /**
   * @brief Fills out[0..N) with a random permutation of [0, N).
   * @param out Receives the per-shape palette-slot assignment.
   */
  static void shuffle_indices(std::array<int, N> &out) {
    for (int i = 0; i < N; ++i)
      out[i] = i;
    std::shuffle(out.begin(), out.end(), hs::random());
  }

  /**
   * @brief Returns the baked LUT for a slot index.
   * @param i Slot index in [0, N).
   * @return Const reference to the baked palette; hot-path lookup is
   *         bank[slot].get(t).
   */
  const BakedPalette &operator[](int i) const {
    assert(i >= 0 && i < N && "MeshPaletteBank index out of range");
    return bank.entries[i];
  }
  /**
   * @brief Returns the baked LUT for a slot index.
   * @param i Slot index in [0, N).
   * @return Mutable reference to the baked palette.
   */
  BakedPalette &operator[](int i) {
    assert(i >= 0 && i < N && "MeshPaletteBank index out of range");
    return bank.entries[i];
  }

  /**
   * @brief Cloneable hook so effects can Persist<MeshPaletteBank> across
   *        compaction.
   * @param src Source bank to copy from.
   * @param dst Destination bank to populate.
   * @param arena Arena receiving the cloned baked LUTs.
   */
  static void clone(const MeshPaletteBank &src, MeshPaletteBank &dst,
                    Arena &arena) {
    BakedPaletteBank::clone(src.bank, dst.bank, arena);
  }

  BakedPaletteBank bank; /**< Underlying baked-palette bank holding N LUTs. */
};

// Ties the source-palette list length to N so editing one without the other
// fails to compile.
static_assert(MeshPaletteBank::sources().size() == MeshPaletteBank::N,
              "MeshPaletteBank::sources() count must equal N");

