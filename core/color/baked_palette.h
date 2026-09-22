/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file baked_palette.h
 * @brief Arena-backed palette tables and crossfade sampling.
 */

#include <cstring>
#include <utility>
#include "color/palette.h"
#include "engine/memory.h"
#include "math/3dmath.h"

/**
 * @brief Read-only view of a 256-entry arena-backed color/alpha table.
 * @details Copies share table storage and never grant mutation rights.
 */
class BakedPalette {
public:
  static constexpr int LUT_SIZE = 256;

  /**
   * @brief Arena bytes a table consumes, including worst-case alignment padding.
   */
  static constexpr size_t required_arena_bytes() {
    return LUT_SIZE * (sizeof(Pixel) + sizeof(uint16_t)) + alignof(Pixel) +
           alignof(uint16_t);
  }

  /**
   * @brief Default-constructs an empty view; bind it before sampling.
   */
  BakedPalette() = default;

  /**
   * @brief Fast lookup with linear interpolation between adjacent entries.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to the last entry).
   * @return The interpolated color.
   */
  Color4 get(float t) const {
    Color4 out;
    sample_into(t, out);
    return out;
  }

  /**
   * @brief Samples only the interpolated RGB channels.
   * @param t Lookup coordinate; clamped to [0, 1].
   * @return The pixel get() gives for the same index, without interpolating
   * alpha.
   */
  __attribute__((always_inline)) Pixel get_color(float t) const {
    assert(colors != nullptr && "BakedPalette::get_color before bake()");
    float idx =
        hs::clamp(t * (LUT_SIZE - 1), 0.0f, static_cast<float>(LUT_SIZE - 1));
    return sample_color_index(idx);
  }

  /**
   * @brief Samples only the interpolated alpha channel.
   * @param t Lookup coordinate; clamped to [0, 1] (NaN folds to the last entry).
   * @return The alpha get() gives for the same index.
   */
  __attribute__((always_inline)) float get_alpha(float t) const {
    assert(alpha_q16 != nullptr && "BakedPalette::get_alpha before bake()");
    float idx =
        hs::clamp(t * (LUT_SIZE - 1), 0.0f, static_cast<float>(LUT_SIZE - 1));
    const int lo = lut_index_lo(idx);
    if (lo >= LUT_SIZE - 1)
      return alpha_q16[LUT_SIZE - 1] * (1.0f / 65535.0f);
    return lerp_q16(alpha_q16[lo], alpha_q16[lo + 1],
                    lut_index_weight(idx, lo)) *
           (1.0f / 65535.0f);
  }

  /**
   * @brief Samples RGB for a coordinate already clamped to [0, 1].
   * @param t Finite lookup coordinate in [0, 1].
   * @return The pixel get() gives for the same index.
   */
  __attribute__((always_inline)) Pixel get_color_unit(float t) const {
    assert(colors != nullptr && "BakedPalette::get_color_unit before bake()");
    // Also traps NaN: both comparisons fail.
    assert(t >= 0.0f && t <= 1.0f);
    return sample_color_index(t * (LUT_SIZE - 1));
  }

  /** @brief Returns a read-only handle to this arena-backed table. */
  BakedPalette view() const { return *this; }

private:
  friend class BakedPaletteStorage;
  static __attribute__((always_inline)) uint16_t lerp_q16(uint16_t a,
                                                          uint16_t b,
                                                          uint16_t weight) {
    const uint32_t inverse = 65535u - weight;
    const uint32_t x =
        static_cast<uint32_t>(a) * inverse + static_cast<uint32_t>(b) * weight;
    return static_cast<uint16_t>((x + (x >> 16) + 32768u) >> 16);
  }

  __attribute__((always_inline)) Pixel sample_color_index(float idx) const {
    if (idx <= 0.0f)
      return colors[0];
    return lut_sample_pixel(colors, LUT_SIZE, idx);
  }

  __attribute__((always_inline)) void sample_into(float t, Color4 &out) const {
    assert(colors != nullptr && alpha_q16 != nullptr &&
           "BakedPaletteStorage::get before bake()");
    // Clamp before the int cast: static_cast<int>(NaN) is UB. hs::clamp maps NaN
    // to the hi bound (last entry) and guarantees idx >= 0.
    float idx =
        hs::clamp(t * (LUT_SIZE - 1), 0.0f, static_cast<float>(LUT_SIZE - 1));
    // Split and weight through the same helpers lut_sample_pixel uses, so this
    // path and get_color/get_color_unit share one spelling of the arithmetic.
    const int lo = lut_index_lo(idx);
    if (lo >= LUT_SIZE - 1) {
      out = Color4(colors[LUT_SIZE - 1],
                   alpha_q16[LUT_SIZE - 1] * (1.0f / 65535.0f));
      return;
    }
    const uint16_t weight = lut_index_weight(idx, lo);
    out = Color4(colors[lo].lerp16(colors[lo + 1], weight),
                 lerp_q16(alpha_q16[lo], alpha_q16[lo + 1], weight) *
                     (1.0f / 65535.0f));
  }
  Pixel *colors = nullptr;
  uint16_t *alpha_q16 = nullptr;
};
/**
 * @brief Exclusive mutation rights to one arena-backed palette table.
 * @details Moves transfer mutation rights. Views remain valid until the arena
 * is reset and observe subsequent rebakes.
 */
class BakedPaletteStorage {
public:
  BakedPaletteStorage() = default;
  BakedPaletteStorage(const BakedPaletteStorage &) = delete;
  BakedPaletteStorage &operator=(const BakedPaletteStorage &) = delete;
  BakedPaletteStorage(BakedPaletteStorage &&other) noexcept
      : table(other.table) {
    other.table.colors = nullptr;
    other.table.alpha_q16 = nullptr;
  }
  BakedPaletteStorage &operator=(BakedPaletteStorage &&other) noexcept {
    if (this != &other) {
      table.colors = std::exchange(other.table.colors, nullptr);
      table.alpha_q16 = std::exchange(other.table.alpha_q16, nullptr);
    }
    return *this;
  }

  /** @brief Borrows the read-only view while this storage object lives. */
  const BakedPalette &view() const & { return table; }
  const BakedPalette &view() const && = delete;
  operator const BakedPalette &() const & { return table; }
  operator const BakedPalette &() const && = delete;

  Color4 get(float t) const { return table.get(t); }
  __attribute__((always_inline)) Pixel get_color(float t) const {
    return table.get_color(t);
  }
  __attribute__((always_inline)) Pixel get_color_unit(float t) const {
    return table.get_color_unit(t);
  }
  __attribute__((always_inline)) float get_alpha(float t) const {
    return table.get_alpha(t);
  }

  /**
   * @brief Bakes any source into a 256-entry LUT in the given arena.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param arena Arena to allocate the LUT from.
   * @param source Source palette or composition to sample.
   * @details Works for a runtime Palette or a compile-time StaticPalette alike.
   */
  template <typename Source>
  HS_COLD_MEMBER void bake(Arena &arena, const Source &source) {
    table.colors = arena.allocate_n<Pixel>(BakedPalette::LUT_SIZE);
    table.alpha_q16 = arena.allocate_n<uint16_t>(BakedPalette::LUT_SIZE);
    rebake(source);
  }

  /**
   * @brief Refills the existing LUT without allocating. Use for animated palettes.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param source Source palette or composition to sample.
   * @details Entry i samples t = i / (BakedPalette::LUT_SIZE - 1), so the last entry lands on
   * t = 1 exactly. A composition with Wrap=true folds that sample back to 0 and
   * collapses its last entry onto its first — bake such sources with Wrap=false.
   * Mirrored sources copy the first half in reverse. Looping sources copy entry
   * zero to entry 255 so the quantized seam is exact.
   */
  template <typename Source> HS_COLD_MEMBER void rebake(const Source &source) {
    static_assert(!palette_wraps_coordinate<Source>(),
                  "BakedPalette cannot rebake a wrapping source");
    HS_CHECK(table.colors != nullptr && table.alpha_q16 != nullptr,
             "BakedPaletteStorage::rebake before bake()");
    bool mirrors = false;
    if constexpr (requires { source.mirrors_domain(); })
      mirrors = source.mirrors_domain();
    bool loops = false;
    if constexpr (requires { source.loops_domain(); })
      loops = source.loops_domain();
    int sample_count = BakedPalette::LUT_SIZE;
    if (mirrors)
      sample_count = BakedPalette::LUT_SIZE / 2;
    else if (loops)
      sample_count = BakedPalette::LUT_SIZE - 1;
    for (int i = 0; i < sample_count; ++i) {
      float t = static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1);
      const Color4 sample = source.get(t);
      table.colors[i] = sample.color;
      table.alpha_q16[i] = frac_to_q16(sample.alpha);
    }
    if (mirrors) {
      for (int i = 0; i < sample_count; ++i) {
        table.colors[BakedPalette::LUT_SIZE - 1 - i] = table.colors[i];
        table.alpha_q16[BakedPalette::LUT_SIZE - 1 - i] = table.alpha_q16[i];
      }
    } else if (loops) {
      table.colors[BakedPalette::LUT_SIZE - 1] = table.colors[0];
      table.alpha_q16[BakedPalette::LUT_SIZE - 1] = table.alpha_q16[0];
    }
  }

  /**
   * @brief Bakes this LUT as the w-blend of two baked palettes.
   * @param arena Arena the blended LUT is allocated from.
   * @param from The w = 0 endpoint; must be baked.
   * @param to The w = 1 endpoint; must be baked.
   * @param w Blend weight in (0, 1).
   * @details Walks the two source LUTs entry-wise with the fixed-point channel
   * lerp — no per-entry float resampling. Neither endpoint may be this palette;
   * the fresh allocation would retarget it before the blend reads it.
   */
  HS_COLD_MEMBER void bake_blend(Arena &arena, const BakedPalette &from,
                                 const BakedPalette &to, float w) {
    HS_CHECK(from.colors && from.alpha_q16 && to.colors && to.alpha_q16,
             "BakedPaletteStorage::bake_blend before bake()");
    HS_CHECK(&from != &table && &to != &table,
             "BakedPaletteStorage::bake_blend endpoint is the output");
    table.colors = arena.allocate_n<Pixel>(BakedPalette::LUT_SIZE);
    table.alpha_q16 = arena.allocate_n<uint16_t>(BakedPalette::LUT_SIZE);
    // Clamp before the cast: w < 0 or NaN is float->int UB, and a NaN weight
    // would otherwise reach every entry's alpha.
    const float wc = hs::clamp(w, 0.0f, 1.0f);
    fill_blend(from, to, wc);
  }

  /**
   * @brief Refills this LUT with a verbatim copy of another baked LUT.
   * @param src Source palette; must be baked and must not alias this storage.
   */
  HS_COLD_MEMBER void rebake_copy(const BakedPalette &src) {
    HS_CHECK(table.colors != nullptr && table.alpha_q16 != nullptr,
             "BakedPaletteStorage::rebake_copy before bake()");
    HS_CHECK(src.colors != nullptr && src.alpha_q16 != nullptr,
             "BakedPaletteStorage::rebake_copy before src bake()");
    HS_CHECK(src.colors != table.colors,
             "BakedPaletteStorage::rebake_copy from itself");
    memcpy(table.colors, src.colors, BakedPalette::LUT_SIZE * sizeof(Pixel));
    memcpy(table.alpha_q16, src.alpha_q16,
           BakedPalette::LUT_SIZE * sizeof(uint16_t));
  }

  /**
   * @brief Refills this LUT with the w-blend of two baked palettes, without
   * allocating.
   * @param from The w = 0 endpoint; must be baked.
   * @param to The w = 1 endpoint; must be baked.
   * @param w Blend weight; clamped to [0, 1] (NaN folds to 1).
   * @details At w <= 0 or w >= 1 the endpoint LUT is copied verbatim, so
   * crossfade boundaries are bit-exact. Neither endpoint may alias this LUT's
   * storage; the endpoints may alias each other.
   */
  HS_COLD_MEMBER void rebake_crossfade(const BakedPalette &from,
                                       const BakedPalette &to, float w) {
    const float wc = hs::clamp(w, 0.0f, 1.0f);
    if (wc == 0.0f) {
      rebake_copy(from);
      return;
    }
    if (wc == 1.0f) {
      rebake_copy(to);
      return;
    }
    HS_CHECK(table.colors != nullptr && table.alpha_q16 != nullptr,
             "BakedPaletteStorage::rebake_crossfade before bake()");
    HS_CHECK(from.colors && from.alpha_q16 && to.colors && to.alpha_q16,
             "BakedPaletteStorage::rebake_crossfade before endpoint bake()");
    HS_CHECK(
        from.colors != table.colors && to.colors != table.colors,
        "BakedPaletteStorage::rebake_crossfade endpoint aliases the output");
    fill_blend(from, to, wc);
  }

  /**
   * @brief Deep-copies the LUT from another BakedPalette into the given arena.
   * @param src Source palette to copy; must already be baked and must not be
   * this palette.
   * @param arena Arena to allocate the new LUT from.
   * @details Used by Persist for arena compaction. The fresh allocation
   * retargets this storage before the copy reads @p src, so a self-clone would
   * memcpy uninitialized arena onto itself.
   */
  void clone_from(const BakedPalette &src, Arena &arena) {
    HS_CHECK(src.colors != nullptr && src.alpha_q16 != nullptr,
             "BakedPaletteStorage::clone_from before src bake()");
    HS_CHECK(&src != &table, "BakedPaletteStorage::clone_from from itself");
    table.colors = arena.allocate_n<Pixel>(BakedPalette::LUT_SIZE);
    table.alpha_q16 = arena.allocate_n<uint16_t>(BakedPalette::LUT_SIZE);
    memcpy(table.colors, src.colors, BakedPalette::LUT_SIZE * sizeof(Pixel));
    memcpy(table.alpha_q16, src.alpha_q16,
           BakedPalette::LUT_SIZE * sizeof(uint16_t));
  }

private:
  BakedPalette table;
  // wc must already be clamped to [0, 1].
  HS_COLD_MEMBER void fill_blend(const BakedPalette &from,
                                 const BakedPalette &to, float wc) {
    const uint16_t weight = frac_to_q16(wc);
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
      table.colors[i] = from.colors[i].lerp16(to.colors[i], weight);
      table.alpha_q16[i] =
          BakedPalette::lerp_q16(from.alpha_q16[i], to.alpha_q16[i], weight);
    }
  }
};

/**
 * @brief Bake-time adapter mapping a LUT coordinate from the cos domain into a
 *        source palette's angle parameter.
 * @tparam Source Type exposing Color4 get(float) const over t = angle/PI.
 * @details Baking through this folds the d -> acos(d)/PI radial mapping into the
 * bake (256 acos per bake, not one per fragment): the fragment lookup keys the
 * LUT by the raw dot product via dot_key(d). dot_key inverts this mapping:
 * u -> d = 1 - 2u, get(u) returns the source at acos(d)/PI.
 */
template <typename Source> struct DotKeyed {
  static constexpr bool WRAPS_COORDINATE = palette_wraps_coordinate<Source>();

  const Source &source;
  Color4 get(float u) const {
    float d = hs::clamp(1.0f - 2.0f * u, -1.0f, 1.0f);
    return source.get(fast_acos(d) / PI_F);
  }
};

/**
 * @brief Wraps a palette source for a DotKeyed bake.
 * @param source Palette sampled over t = angle/PI; must outlive the bake call.
 */
template <typename Source>
inline DotKeyed<Source> dot_keyed(const Source &source) {
  return DotKeyed<Source>{source};
}

/// LUT coordinate for cos-value d = dot(axis, v); inverse of DotKeyed's mapping.
inline float dot_key(float d) {
  return (1.0f - hs::clamp(d, -1.0f, 1.0f)) * 0.5f;
}

/**
 * @brief Resolves a crossfade into a read-only table view.
 * @param arena Receives a fresh table only for an interior blend weight.
 * @param from The w = 0 endpoint.
 * @param to The w = 1 endpoint.
 * @param w Blend weight; endpoint weights share the endpoint table.
 */
inline BakedPalette bake_palette_blend(Arena &arena, const BakedPalette &from,
                                       const BakedPalette &to, float w) {
  if (w <= 0.0f)
    return from;
  if (w >= 1.0f)
    return to;
  BakedPaletteStorage blended;
  blended.bake_blend(arena, from, to, w);
  return blended.view();
}

/**
 * @brief Bank of N baked palettes for bulk Persist/clone operations.
 */
struct BakedPaletteBank {
  static constexpr int N = 6;
  BakedPaletteStorage entries[N];

  /**
   * @brief Deep-copies all entries into a target arena.
   * @param src Source bank to copy from.
   * @param dst Destination bank to fill.
   * @param arena Arena to allocate the cloned LUTs from.
   * @details Required by Cloneable.
   */
  HS_COLD_MEMBER static void clone(const BakedPaletteBank &src,
                                   BakedPaletteBank &dst, Arena &arena) {
    for (int i = 0; i < N; ++i)
      dst.entries[i].clone_from(src.entries[i], arena);
  }
};
