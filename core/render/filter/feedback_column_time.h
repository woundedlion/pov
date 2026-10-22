/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <cstdint>
#include "core/color/color.h"
#include "core/memory.h"
#include "core/render/filter/feedback_style.h"

/**
 * @file feedback_column_time.h
 * @brief Filter::Pixel::FeedbackColumnTime and FeedbackFadeTable: the elapsed
 * time between a feedback target pixel and the source it reads, and the
 * Style's fade at each elapsed time.
 */

namespace Filter {

namespace Pixel {

/**
 * @brief Elapsed-time model of the feedback warp under column-time display.
 * @tparam W Canvas width in pixels.
 * @details Column x is lit at frame fraction t(x) = (x % (W / 2)) / (W / 2),
 * and each frame's buffer holds every column at its own lit time. A target
 * pixel reading source column s of the previous frame therefore spans
 * 1 + t(x) - t(s) frames, a read between two columns their lit times blended
 * as its taps are. A Style warp displacing D columns per frame reads, within
 * the target's arm, the source lit `HALF / (HALF + D)` frames earlier; past
 * the arm's forward cut it reads the arm lit a frame before that, twice as
 * far; behind its backward cut lies content no stored frame holds, which
 * slides from the cut column lit last to the first column of the arm.
 */
template <int W> struct FeedbackColumnTime {
  static_assert(W % 2 == 0, "FeedbackColumnTime: two arms need an even width");
  /** @brief Columns per arm. */
  static constexpr int HALF = W / 2;
  /** @brief Reciprocal of HALF. */
  static constexpr float INVERSE_HALF = 1.0f / HALF;

  /** @brief A resolved source read. */
  struct Read {
    float x;       ///< Source column.
    float y;       ///< Source row.
    float elapsed; ///< Frames between the source's and the target's lit times.
  };

  /** @brief A source position. */
  struct Source {
    float x; ///< Column.
    float y; ///< Row.
  };

  /**
   * @brief Frames elapsed for a one-frame displacement read within its arm.
   * @param columns The Style's one-frame column displacement of the source.
   * @return The scale that displacement takes: the frames between the target
   *         and the source it lands on, at most 2.
   */
  static float elapsed_scale(float columns) {
    return HALF / (HALF + fmaxf(columns, -0.5f * HALF));
  }

  /**
   * @brief First column of the arm holding @p x.
   * @param x Column in [0, W).
   * @return 0 or HALF.
   */
  static int arm_start(int x) { return x >= HALF ? HALF : 0; }

  /**
   * @brief A source column relative to its target's arm.
   * @param x Source column, within one turn below the arm's first column.
   * @param arm First column of the target's arm.
   * @return The column past the arm's first, in [-0.5, W - 0.5): below
   *         HALF - 0.5 the nearest column lies in the arm.
   */
  static HS_HOT_INLINE float arm_column(float x, int arm) {
    const float column = x - arm;
    return column < -0.5f ? column + W : column;
  }

  /**
   * @brief Lit time of a bilinear read, in frames of the previous sweep.
   * @param column Read column past an arm's first, in [-1, 2 * HALF).
   * @return Frame fraction: linear within an arm, and between an arm's last
   *         column and the next arm's first the two taps' times blended by
   *         their weights.
   */
  static HS_HOT_INLINE float lit_time(float column) {
    const float in_arm = column < 0.0f    ? column + HALF
                         : column >= HALF ? column - HALF
                                          : column;
    return in_arm > HALF - 1 ? (HALF - in_arm) * ((HALF - 1) * INVERSE_HALF)
                             : in_arm * INVERSE_HALF;
  }

  /**
   * @brief A read within the target's arm.
   * @param target Target column, with its half-resolution lane offset.
   * @param arm First column of the target's arm.
   * @param column Source's arm_column(), below HALF - 0.5.
   * @param y Source row.
   * @return The read with its elapsed frames.
   */
  static HS_HOT_INLINE Read read_within(float target, int arm, float column,
                                        float y) {
    return {arm + column, y,
            1.0f + (target - arm) * INVERSE_HALF - lit_time(column)};
  }

  /**
   * @brief A read for a source past one of its arm's cuts.
   * @param target Target column, with its half-resolution lane offset.
   * @param y Target row.
   * @param arm First column of the target's arm.
   * @param within Source displaced by the elapsed-scaled warp, past the arm,
   *        its column unwrapped to within half a turn of @p target.
   * @param beyond Returns the source displaced twice as far, its column
   *        unwrapped likewise.
   * @return Past the forward cut, the twice-displaced source; behind the
   *         backward cut, the slide from the previous arm's last column to
   *         the arm's first.
   */
  template <typename BeyondFnT>
  static HS_FLASH_MEMBER Read read_past_cut(float target, float y, int arm,
                                            const Source &within,
                                            BeyondFnT &&beyond) {
    const float first = static_cast<float>(arm);
    const float lit = 1.0f + (target - first) * INVERSE_HALF;
    if (within.x >= first) {
      const Source far = beyond();
      const float column =
          fminf(fmaxf(far.x - first, HALF - 0.5f), 2.0f * HALF - 1.0f);
      return {first + column, far.y, lit - lit_time(column)};
    }
    const float slide =
        fminf(fmaxf((target - first + 0.5f) / (target - within.x), 0.0f), 1.0f);
    return {first - 1.0f + slide, y + (within.y - y) * slide,
            lit - lit_time(slide - 1.0f)};
  }
};

/**
 * @brief A Style's fade, and its hue rotation, at a lattice of elapsed times.
 * @details Entry `l` applies `l / STEPS_PER_FRAME` frames of the Style's
 * per-frame fade and hue shift; entry UNIT_LEVEL equals the per-frame
 * transform the feedback filter applies with column time off. Storage is
 * arena-owned and rebuilt when the Style's fade or hue shift changes.
 */
struct FeedbackFadeTable {
  /** @brief Entries per frame of elapsed time. */
  static constexpr int STEPS_PER_FRAME = 64;
  /** @brief Longest elapsed time the table holds, frames: no read spans two
   *  full frames. */
  static constexpr int MAX_FRAMES = 2;
  /** @brief Entry count. */
  static constexpr int LEVELS = STEPS_PER_FRAME * MAX_FRAMES + 1;
  /** @brief Entry of one frame. */
  static constexpr int UNIT_LEVEL = STEPS_PER_FRAME;
  /** @brief Floats per entry's cbrt-LMS matrix. */
  static constexpr int MATRIX = 9;
  /** @brief Arena bytes for one table. */
  static constexpr size_t BYTES = LEVELS * (1 + MATRIX) * sizeof(float);

  /** @brief A built table's entries. */
  struct View {
    const float *fades; ///< Brightness multiplier per entry.
    const float *hue_k; ///< cbrt-LMS matrix per entry, MATRIX floats each,
                        ///< scaled for hue_fade_apply().
  };

  /**
   * @brief Entry nearest an elapsed time.
   * @param elapsed Frames in [0, MAX_FRAMES], as every FeedbackColumnTime
   *        read spans.
   * @return Entry index.
   */
  static HS_HOT_INLINE int level(float elapsed) {
    return static_cast<int>(elapsed * STEPS_PER_FRAME + 0.5f);
  }

  /**
   * @brief Allocates persistent storage; call again after any arena reset.
   * @param arena Persistent arena supplying BYTES.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    fades = arena.allocate_n<float>(LEVELS);
    hue_k = arena.allocate_n<float>(LEVELS * MATRIX);
    valid = false;
#ifndef NDEBUG
    stamp.record(arena);
#endif
  }

  /**
   * @brief The table for @p style, rebuilt into storage when its fade or hue
   *        shift changed, or built into @p scratch without storage.
   * @param style Style whose sync_hue() ran this frame.
   * @param scratch Arena for a storage-less table.
   * @return The built entries.
   */
  HS_FLASH_MEMBER View acquire(const ::Feedback::Style &style, Arena &scratch) {
    if (!fades) {
      float *scratch_fades = scratch.allocate_n<float>(LEVELS);
      float *scratch_k = scratch.allocate_n<float>(LEVELS * MATRIX);
      build(style, scratch_fades, scratch_k);
      return {scratch_fades, scratch_k};
    }
#ifndef NDEBUG
    HS_ASSERT_BLOCK_ALIVE(stamp, hue_k, LEVELS * MATRIX * sizeof(float),
                          "Pixel::Feedback fade table");
#endif
    if (!valid || style.fade != built_fade ||
        style.hue_shift != built_hue_shift) {
      build(style, fades, hue_k);
      built_fade = style.fade;
      built_hue_shift = style.hue_shift;
      valid = true;
    }
    return {fades, hue_k};
  }

private:
  HS_COLD_MEMBER static void build(const ::Feedback::Style &style,
                                   float *out_fades, float *out_k) {
    const float fade = style.fade;
    const float frame_shift =
        fade == 0.0f ? 0.0f : style.hue_shift * -logf(fade);
    for (int l = 0; l < LEVELS; ++l) {
      const float elapsed = static_cast<float>(l) / STEPS_PER_FRAME;
      float *k = out_k + l * MATRIX;
      if (l == UNIT_LEVEL) {
        out_fades[l] = fade;
        const float sc = math::fast_cbrt(fade * Q16_TO_UNIT);
        for (int i = 0; i < MATRIX; ++i)
          k[i] = style.hue_k[i] * sc;
        continue;
      }
      const float faded = powf(fade, elapsed);
      out_fades[l] = faded;
      const float shift = frame_shift * elapsed;
      float ca = 1.0f, sa = 0.0f;
      if (shift != 0.0f)
        turn_to_unit_cos_sin(shift, ca, sa);
      hue_rotate_lms_matrix(ca, sa, k);
      const float sc = math::fast_cbrt(faded * Q16_TO_UNIT);
      for (int i = 0; i < MATRIX; ++i)
        k[i] *= sc;
    }
  }

  float *fades = nullptr;       /**< Arena-owned fades, LEVELS entries. */
  float *hue_k = nullptr;       /**< Arena-owned matrices, LEVELS * MATRIX. */
  float built_fade = 0.0f;      /**< Fade the stored entries were built for. */
  float built_hue_shift = 0.0f; /**< Hue shift they were built for. */
  bool valid = false;           /**< Whether the stored entries are built. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the storage was allocated. */
#endif
};

} // namespace Pixel

} // namespace Filter
