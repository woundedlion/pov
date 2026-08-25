/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file palettes.h
 * @brief The named procedural palettes, plus MeshPaletteBank, the palette
 *        bank the mesh effects share.
 */

#include <array>

#include "color/color.h"

/** @brief Unpacks a parenthesized coefficient triple into its three floats. */
#define HS_PALETTE_VEC3(x, y, z) x, y, z

/**
 * @brief X-macro roster of the named procedural palettes.
 * @param X Macro applied as X(name, A, B, C, D), where each coefficient vec3 is
 *          a parenthesized triple unpacked by HS_PALETTE_VEC3.
 * @details Expanded here to declare the Palettes:: instances, and in
 *          targets/wasm/math_exports.h to export name + coefficients, so the
 *          browser tool's mirror of this table is checked against the same
 *          literals the engine compiles.
 */
#define HS_PROCEDURAL_PALETTE_LIST(X)                                          \
  X(DARK_RAINBOW, (0.367f, 0.367f, 0.367f), (0.500f, 0.500f, 0.500f),          \
    (1.000f, 1.000f, 1.000f), (0.000f, 0.330f, 0.670f))                        \
  X(BLOOD_STREAM, (0.169f, 0.169f, 0.169f), (0.313f, 0.313f, 0.313f),          \
    (0.231f, 0.231f, 0.231f), (0.036f, 0.366f, 0.706f))                        \
  X(VINTAGE_SUNSET, (0.256f, 0.256f, 0.256f), (0.500f, 0.080f, 0.500f),        \
    (0.277f, 0.277f, 0.277f), (0.000f, 0.330f, 0.670f))                        \
  X(RICH_SUNSET, (0.309f, 0.500f, 0.500f), (1.000f, 1.000f, 0.500f),           \
    (0.149f, 0.148f, 0.149f), (0.132f, 0.222f, 0.521f))                        \
  X(UNDERSEA, (0.000f, 0.000f, 0.000f), (0.500f, 0.276f, 0.423f),              \
    (0.296f, 0.296f, 0.296f), (0.374f, 0.941f, 0.000f))                        \
  X(LATE_SUNSET, (0.337f, 0.500f, 0.096f), (0.500f, 1.000f, 0.176f),           \
    (0.261f, 0.261f, 0.261f), (0.153f, 0.483f, 0.773f))                        \
  X(MANGO_PEEL, (0.500f, 0.500f, 0.500f), (0.500f, 0.080f, 0.500f),            \
    (0.431f, 0.431f, 0.431f), (0.566f, 0.896f, 0.236f))                        \
  X(ICE_MELT, (0.500f, 0.500f, 0.500f), (0.500f, 0.500f, 0.500f),              \
    (0.083f, 0.147f, 0.082f), (0.579f, 0.353f, 0.244f))                        \
  X(LEMON_LIME, (0.455f, 0.455f, 0.455f), (0.571f, 0.151f, 0.571f),            \
    (0.320f, 0.320f, 0.320f), (0.087f, 0.979f, 0.319f))                        \
  X(ALGAE, (0.210f, 0.210f, 0.210f), (0.500f, 1.000f, 0.021f),                 \
    (0.086f, 0.086f, 0.075f), (0.419f, 0.213f, 0.436f))                        \
  X(EMBERS, (0.500f, 0.500f, 0.500f), (0.500f, 0.500f, 0.500f),                \
    (0.265f, 0.285f, 0.198f), (0.577f, 0.440f, 0.358f))                        \
  X(FIRE_GLOW, (0.000f, 0.000f, 0.000f), (0.560f, 0.560f, 0.560f),             \
    (0.216f, 0.346f, 0.174f), (0.756f, 0.542f, 0.279f))                        \
  X(DARK_PRIMARY, (0.500f, 0.500f, 0.500f), (0.500f, 0.610f, 0.500f),          \
    (0.746f, 0.347f, 0.000f), (0.187f, 0.417f, 0.670f))                        \
  X(MAUVE_FADE, (0.583f, 0.000f, 0.583f), (1.000f, 0.000f, 1.000f),            \
    (0.191f, 0.348f, 0.191f), (0.175f, 0.045f, 0.150f))                        \
  X(LAVENDER_LAKE, (0.473f, 0.473f, 0.473f), (0.500f, 0.500f, 0.500f),         \
    (0.364f, 0.124f, 0.528f), (0.142f, 0.378f, 0.876f))                        \
  X(DESERT_ROSE, (0.500f, 0.500f, 0.500f), (0.500f, 0.270f, 0.442f),           \
    (0.303f, 1.012f, 0.585f), (0.985f, 0.720f, 0.212f))                        \
  X(BRUISED_MOSS, (0.500f, 0.500f, 0.500f), (0.500f, 0.500f, 0.500f),          \
    (0.142f, 0.252f, 0.000f), (0.492f, 0.200f, 0.670f))                        \
  X(BRUISED_BANANA, (0.175f, 0.470f, 0.171f), (1.000f, 0.622f, 0.000f),        \
    (0.191f, 0.191f, 0.191f), (0.629f, -0.417f, 0.444f))                       \
  X(BRIGHT_SUNRISE, (0.620f, 0.620f, 0.620f), (0.742f, 0.742f, 0.742f),        \
    (0.162f, 0.286f, 0.012f), (0.090f, 0.205f, 0.688f))                        \
  X(FIRE_AND_ICE, (0.500f, 0.500f, 0.500f), (0.500f, 0.500f, 0.500f),          \
    (0.955f, 1.004f, 0.910f), (0.167f, 0.018f, 0.930f))                        \
  X(PEACH_POP, (1.000f, 0.144f, 0.175f), (0.543f, 0.543f, 0.543f),             \
    (0.507f, 0.409f, 0.507f), (0.001f, 0.002f, 0.620f))                        \
  X(POPPED_PEACH, (1.000f, 0.144f, 0.175f), (0.543f, 0.543f, 0.543f),          \
    (-0.507f, -0.409f, -0.507f), (0.508f, 0.411f, 1.127f))                     \
  X(BLUE_LAGOON, (0.253f, 0.500f, 1.000f), (0.500f, 0.844f, 1.000f),           \
    (0.232086f, 0.232086f, 0.232086f), (0.279882f, 0.609882f, 0.949882f))      \
  X(ORANGE_CRUSH, (0.575f, 0.168f, 0.464f), (0.406f, 0.697f, 0.357f),          \
    (0.000f, -0.10051f, -0.042778f), (0.141f, 0.25551f, 0.579778f))            \
  X(PLUM_SUNRISE, (0.407f, 0.000f, 0.296f), (0.332f, 0.592f, 0.029f),          \
    (0.358961f, 0.331145f, 0.274519f), (0.500342f, 0.505109f, 0.278634f))      \
  X(CORAL_BLUE, (0.4f, 0.347f, 0.801f), (0.5f, 0.303f, 0.5f),                  \
    (0.75363518f, -0.20031623f, 0.0110030736f),                                \
    (0.9144297f, -0.16868377f, 0.40006184f))                                   \
  X(BRUISED_MANGO, (0.385f, 0.470f, 0.171f), (1.000f, 0.518f, 0.000f),         \
    (0.191f, 0.191f, 5.000f), (0.619f, -0.427f, 0.887f))

namespace Palettes {

#define HS_DECLARE_PALETTE(name, A, B, C, D)                                   \
  inline constexpr ProceduralPalette name(                                     \
      {HS_PALETTE_VEC3 A}, {HS_PALETTE_VEC3 B}, {HS_PALETTE_VEC3 C},           \
      {HS_PALETTE_VEC3 D});
HS_PROCEDURAL_PALETTE_LIST(HS_DECLARE_PALETTE)
#undef HS_DECLARE_PALETTE

} // namespace Palettes

/**
 * @brief Shared mesh-effect palette bank used by HankinSolids / IslamicStars
 *        and any future mesh effect.
 * @details Bundles the standard source-palette set, the bake-all step, and the
 *          per-shape index shuffle these effects share. Zero-overhead: the
 *          source list is constexpr and every accessor is a thin inline wrapper
 *          over BakedPaletteBank, so the per-pixel lookup remains
 *          BakedPalette::get() with no added indirection.
 */
struct MeshPaletteBank {
  static constexpr int N = BakedPaletteBank::N;

  /**
   * @brief Arena bytes bake_all() consumes, including worst-case per-palette
   *        alignment padding.
   */
  static constexpr size_t required_arena_bytes() {
    return N * BakedPalette::required_arena_bytes();
  }

  /** @brief Shared source palettes, in bank-slot order. */
  static constexpr auto sources() {
    return std::array{&Palettes::EMBERS,         &Palettes::RICH_SUNSET,
                      &Palettes::BRIGHT_SUNRISE, &Palettes::BRUISED_MOSS,
                      &Palettes::LAVENDER_LAKE,  &Palettes::POPPED_PEACH};
  }

  /**
   * @brief (Re)bakes every source palette into the arena.
   * @param arena Destination arena receiving each BakedPalette's parallel
   *        256-entry Pixel and Q16-alpha arrays.
   */
  HS_COLD_MEMBER void bake_all(Arena &arena) {
    constexpr auto sources = MeshPaletteBank::sources();
    for (int i = 0; i < N; ++i)
      bank.entries[i].bake(arena, *sources[i]);
  }

  /**
   * @brief Fills out[0..N) with a random permutation of [0, N).
   * @param out Receives the per-shape palette-slot assignment.
   */
  static HS_COLD_MEMBER void shuffle_indices(std::array<int, N> &out) {
    for (int i = 0; i < N; ++i)
      out[i] = i;
    hs::shuffle(out.begin(), out.end());
  }

  /**
   * @brief Folds a mesh topology class onto a bank slot.
   * @param cls Face topology class; class ids are dense but unbounded, so the
   *        fold aliases every N-th class onto the same slot.
   * @return Slot index in [0, N).
   */
  static int slot_of(int cls) { return wrap(cls, N); }

  /**
   * @brief Assigns every face the palette its topology class maps to.
   * @param topology Per-face topology-class indices.
   * @param faces Number of faces to assign.
   * @param slots Class slot -> palette index (see shuffle_indices).
   * @param out Receives one palette index per face.
   */
  static HS_COLD_MEMBER void assign_by_class(const uint16_t *topology,
                                             size_t faces,
                                             const std::array<int, N> &slots,
                                             uint8_t *out) {
    for (size_t f = 0; f < faces; ++f)
      out[f] = static_cast<uint8_t>(slots[slot_of(topology[f])]);
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

static_assert(MeshPaletteBank::sources().size() == MeshPaletteBank::N,
              "MeshPaletteBank source count must equal N");
