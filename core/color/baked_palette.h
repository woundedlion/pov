/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file baked_palette.h
 * @brief Arena-backed palette tables and crossfade sampling.
 */

#include <cstring>
#include "color/palette.h"
#include "engine/memory.h"
#include "math/3dmath.h"

/**
 * @brief Pre-baked 256-entry color/alpha LUT allocated in an arena.
 * @details Stores parallel Pixel and Q16-alpha tables with lerp
 * interpolation. Not a Palette subclass — call get(t) directly for
 * zero-overhead lookups. A copy is a non-owning handle onto the source's
 * arena storage: rebake() through any copy writes into every alias. Use
 * clone_from() for an independent LUT. A handle marked by mark_aliased()
 * refuses rebake() outright.
 */
class BakedPalette {
public:
  static constexpr int LUT_SIZE = 256;

  /**
   * @brief Arena bytes bake() consumes, including worst-case alignment padding.
   */
  static constexpr size_t required_arena_bytes() {
    return LUT_SIZE * (sizeof(Pixel) + sizeof(uint16_t)) + alignof(Pixel) +
           alignof(uint16_t);
  }

  /**
   * @brief Default-constructs an unbaked palette (bake() before use).
   */
  BakedPalette() = default;

  /**
   * @brief Bakes any source into a 256-entry LUT in the given arena.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param arena Arena to allocate the LUT from.
   * @param source Source palette or composition to sample.
   * @details Works for a runtime Palette or a compile-time StaticPalette alike.
   */
  template <typename Source>
  HS_COLD_MEMBER void bake(Arena &arena, const Source &source) {
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    aliased = false;
    rebake(source);
  }

  /**
   * @brief Refills the existing LUT without allocating. Use for animated palettes.
   * @tparam Source Type exposing Color4 get(float) const.
   * @param source Source palette or composition to sample.
   * @details Entry i samples t = i / (LUT_SIZE - 1), so the last entry lands on
   * t = 1 exactly. A composition with Wrap=true folds that sample back to 0 and
   * collapses its last entry onto its first — bake such sources with Wrap=false.
   * Mirrored sources copy the first half in reverse. Looping sources copy entry
   * zero to entry 255 so the quantized seam is exact.
   */
  template <typename Source> HS_COLD_MEMBER void rebake(const Source &source) {
    static_assert(!palette_wraps_coordinate<Source>(),
                  "BakedPalette cannot rebake a wrapping source");
    HS_CHECK(colors != nullptr && alpha_q16 != nullptr,
             "BakedPalette::rebake before bake()");
    HS_CHECK(!aliased, "BakedPalette::rebake through an aliasing handle");
    bool mirrors = false;
    if constexpr (requires { source.mirrors_domain(); })
      mirrors = source.mirrors_domain();
    bool loops = false;
    if constexpr (requires { source.loops_domain(); })
      loops = source.loops_domain();
    int sample_count = LUT_SIZE;
    if (mirrors)
      sample_count = LUT_SIZE / 2;
    else if (loops)
      sample_count = LUT_SIZE - 1;
    for (int i = 0; i < sample_count; ++i) {
      float t = static_cast<float>(i) / (LUT_SIZE - 1);
      const Color4 sample = source.get(t);
      colors[i] = sample.color;
      alpha_q16[i] = frac_to_q16(sample.alpha);
    }
    if (mirrors) {
      for (int i = 0; i < sample_count; ++i) {
        colors[LUT_SIZE - 1 - i] = colors[i];
        alpha_q16[LUT_SIZE - 1 - i] = alpha_q16[i];
      }
    } else if (loops) {
      colors[LUT_SIZE - 1] = colors[0];
      alpha_q16[LUT_SIZE - 1] = alpha_q16[0];
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
             "BakedPalette::bake_blend before bake()");
    HS_CHECK(&from != this && &to != this,
             "BakedPalette::bake_blend endpoint is the output");
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    aliased = false;
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
    HS_CHECK(colors != nullptr && alpha_q16 != nullptr,
             "BakedPalette::rebake_copy before bake()");
    HS_CHECK(!aliased, "BakedPalette::rebake_copy through an aliasing handle");
    HS_CHECK(src.colors != nullptr && src.alpha_q16 != nullptr,
             "BakedPalette::rebake_copy before src bake()");
    HS_CHECK(src.colors != colors, "BakedPalette::rebake_copy from itself");
    memcpy(colors, src.colors, LUT_SIZE * sizeof(Pixel));
    memcpy(alpha_q16, src.alpha_q16, LUT_SIZE * sizeof(uint16_t));
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
    HS_CHECK(colors != nullptr && alpha_q16 != nullptr,
             "BakedPalette::rebake_crossfade before bake()");
    HS_CHECK(!aliased,
             "BakedPalette::rebake_crossfade through an aliasing handle");
    HS_CHECK(from.colors && from.alpha_q16 && to.colors && to.alpha_q16,
             "BakedPalette::rebake_crossfade before endpoint bake()");
    HS_CHECK(from.colors != colors && to.colors != colors,
             "BakedPalette::rebake_crossfade endpoint aliases the output");
    fill_blend(from, to, wc);
  }

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

  /**
   * @brief Deep-copies the LUT from another BakedPalette into the given arena.
   * @param src Source palette to copy; must already be baked and must not be
   * this palette.
   * @param arena Arena to allocate the new LUT from.
   * @details Used by Persist for arena compaction. The fresh allocation
   * retargets this handle before the copy reads @p src, so a self-clone would
   * memcpy uninitialized arena onto itself.
   */
  void clone_from(const BakedPalette &src, Arena &arena) {
    HS_CHECK(src.colors != nullptr && src.alpha_q16 != nullptr,
             "BakedPalette::clone_from before src bake()");
    HS_CHECK(&src != this, "BakedPalette::clone_from from itself");
    colors = arena.allocate_n<Pixel>(LUT_SIZE);
    alpha_q16 = arena.allocate_n<uint16_t>(LUT_SIZE);
    aliased = false;
    memcpy(colors, src.colors, LUT_SIZE * sizeof(Pixel));
    memcpy(alpha_q16, src.alpha_q16, LUT_SIZE * sizeof(uint16_t));
  }

  /**
   * @brief Marks this handle as a non-owning view onto another palette's LUT.
   * @details rebake() traps afterwards; bake(), bake_blend() and clone_from()
   * clear the mark by giving the handle storage of its own.
   */
  void mark_aliased() { aliased = true; }

private:
  // wc must already be clamped to [0, 1].
  HS_COLD_MEMBER void fill_blend(const BakedPalette &from,
                                 const BakedPalette &to, float wc) {
    const uint16_t weight = frac_to_q16(wc);
    for (int i = 0; i < LUT_SIZE; ++i) {
      colors[i] = from.colors[i].lerp16(to.colors[i], weight);
      alpha_q16[i] = lerp_q16(from.alpha_q16[i], to.alpha_q16[i], weight);
    }
  }

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
           "BakedPalette::get before bake()");
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
  bool aliased = false;
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
 * @brief Resolves a (from, to) baked-LUT pair at one crossfade weight.
 * @param dst Receives the resolved palette (a default-constructed unbaked
 * instance is fine).
 * @param arena Arena receiving the blended LUT when one is baked.
 * @param from The w = 0 endpoint.
 * @param to The w = 1 endpoint.
 * @param w Blend weight.
 * @details Weights at or beyond an endpoint alias that endpoint's LUT storage
 * (bitwise-exact, no allocation), so crossfade boundaries are exact by
 * construction; only 0 < w < 1 bakes a blended LUT into the arena. @p dst is
 * therefore read-only on the endpoint paths — it is marked aliased there, so a
 * rebake that would overwrite the shared endpoint LUT traps instead.
 */
inline void bake_palette_blend(BakedPalette &dst, Arena &arena,
                               const BakedPalette &from, const BakedPalette &to,
                               float w) {
  if (w <= 0.0f) {
    dst = from;
    dst.mark_aliased();
  } else if (w >= 1.0f) {
    dst = to;
    dst.mark_aliased();
  } else {
    dst.bake_blend(arena, from, to, w);
  }
}

/**
 * @brief Bank of N baked palettes for bulk Persist/clone operations.
 */
struct BakedPaletteBank {
  static constexpr int N = 6;
  BakedPalette entries[N];

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
