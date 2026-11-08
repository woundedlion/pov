/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file hue_rotation_row.h
 * @brief Eager and lazy hue-rotation rows over a bound color. */

#include "core/color/color_space.h"
#include "core/math/periodic.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>

/** @brief Hue-rotated color knots over a bound turn domain.
 * @tparam CELLS Interpolation cell count. */
template <int CELLS> struct HueRotationRow {
  static constexpr int KNOTS = CELLS + 1; /**< Knots per table. */
  /** @brief Persistent bytes of the knot storage, alignment included. */
  static constexpr size_t BYTES = KNOTS * sizeof(Pixel) + alignof(Pixel);
  Pixel *knots = nullptr; /**< KNOTS hue-rotated colors of the bound color. */
  uint64_t valid[(CELLS + 64) / 64] =
      {}; /**< Knots sample_lazy() has baked since the last bind(). */
  const HueRotateBase *base = nullptr; /**< Bound color. */
  float domain = 0.0f;                 /**< Hue-turn interval the knots span. */
  bool cyclic = false;                 /**< Amounts past the domain wrap. */

  /**
   * @brief Binds the table to one color and clears knot validity.
   * @param color_base Precomputed OKLab base; must outlive
   * the row's samples.
   * @param hue_domain Hue-turn interval the knots span.
   * @param wraps Whether amounts past the domain wrap.
   */
  HS_HOT_INLINE void bind(const HueRotateBase &color_base, float hue_domain,
                          bool wraps) {
    base = &color_base;
    domain = hue_domain;
    cyclic = wraps;
    for (uint64_t &word : valid)
      word = 0;
  }

  /**
   * @brief Counts the knots a sample in [0, @p max_amount] can read.
   * @param max_amount Largest hue offset (turns) any sample asks for.
   * @return Knots 0..count - 1 are the only ones sample() reads.
   * @details A wrapping lookup that can reach a full turn, or a domain that
   * is not positive, reads anywhere.
   */
  int knots_through(float max_amount) const {
    if (!(domain > 0.0f) || (cyclic && max_amount >= domain))
      return KNOTS;
    const float x = hs::clamp(max_amount / domain, 0.0f, 1.0f) * CELLS;
    return std::min(static_cast<int>(x) + 2, KNOTS);
  }

  /**
   * @brief Bakes the first @p count knots for the bound color.
   * @param count Knots to bake, at most KNOTS.
   */
  HS_HOT_FLASH_MEMBER void prepare(int count = KNOTS) {
    for (int i = 0; i < count; ++i)
      knots[i] =
          hue_rotate(*base, domain * (static_cast<float>(i) / CELLS)).color;
  }

  /** @brief Samples knots already baked by prepare().
   * @param amount Hue offset in turns.
   * @return Interpolated color. */
  HS_O3_FN Pixel sample(float amount) const {
    return sample_with(amount, [](int) {});
  }

  /** @brief Samples the table, baking each knot it reads on first use.
   * @param amount Hue offset in turns.
   * @return Interpolated color. */
  HS_O3_FN Pixel sample_lazy(float amount) {
    return sample_with(amount, [&](int index) HS_HOT_INLINE {
      const uint64_t bit = uint64_t{1} << (index & 63);
      uint64_t &word = valid[index >> 6];
      if (!(word & bit)) {
        knots[index] =
            hue_rotate(*base, domain * (static_cast<float>(index) / CELLS))
                .color;
        word |= bit;
      }
    });
  }

private:
  /**
   * @brief Interpolates the table at a hue offset.
   * @tparam Ensure Knot preparation callback.
   * @param amount Hue offset in turns.
   * @param ensure Called with every knot index read, before the read.
   * @return Interpolated color.
   */
  template <typename Ensure>
  Pixel sample_with(float amount, Ensure ensure) const {
    HS_AUDIT_CHECK(base, "HueRotationRow sampled before bind()");
    float t = amount / domain;
    t = cyclic ? math::wrap_t(t) : hs::clamp(t, 0.0f, 1.0f);
    float x = t * CELLS;
    if (x >= CELLS) {
      ensure(CELLS);
      return knots[CELLS];
    }
    int i = static_cast<int>(x);
    ensure(i);
    ensure(i + 1);
    return knots[i].lerp16(knots[i + 1], frac_to_q16(x - i));
  }
};
