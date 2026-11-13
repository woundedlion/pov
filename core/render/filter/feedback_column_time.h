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
 * pixel reading source column s spans 1 + (x - s) / HALF frames, s unwrapped
 * to within half a turn of x. A Style warp displacing D columns per frame
 * reads, within the target's arm, the source lit `HALF / (HALF + D)` frames
 * earlier. A source behind the arm's first column lies in the frame before the
 * previous one; one past the arm's last column lies in the frame being drawn.
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
   * @brief The read of a source for a target.
   * @param target Target column, with its half-resolution lane offset.
   * @param source Source column, unwrapped to within half a turn of @p target.
   * @param y Source row.
   * @return The read, spanning the frames between the source's lit time and
   *         the target's.
   */
  static HS_HOT_INLINE Read read(float target, float source, float y) {
    return {source, y, 1.0f + (target - source) * INVERSE_HALF};
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
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
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
  HS_FLASH_MEMBER static void build(const ::Feedback::Style &style,
                                    float *out_fades, float *out_k) {
    const float fade = style.fade;
    const float frame_shift = style.frame_hue_shift();
    for (int l = 0; l < LEVELS; ++l) {
      const float elapsed = static_cast<float>(l) / STEPS_PER_FRAME;
      float *k = out_k + l * MATRIX;
      if (l == UNIT_LEVEL) {
        out_fades[l] = fade;
        ::Feedback::scale_hue_matrix(style.hue_k, fade, k);
        continue;
      }
      const float faded = powf(fade, elapsed);
      out_fades[l] = faded;
      const float shift = frame_shift * elapsed;
      float ca = 1.0f, sa = 0.0f;
      if (shift != 0.0f)
        turn_to_unit_cos_sin(shift, ca, sa);
      hue_rotate_lms_matrix(ca, sa, k);
      ::Feedback::scale_hue_matrix(k, faded, k);
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

namespace detail {

/** @brief Where a lane's source lies relative to its target's arm. */
enum class Side : uint8_t {
  WITHIN, ///< Within the arm: the previous frame.
  BEHIND, ///< Behind the arm's first column: the stale frame before it.
  AHEAD,  ///< From the arm's last column on: the frame being drawn.
};

/** @brief Lane columns of a row whose sources read the frame as plotted. */
struct AheadSpan {
  int16_t first; ///< First such column past each arm's first; HALF if none.
  int16_t last;  ///< Last such column, -1 if none; finish() extends it to
                 ///< the last lane whose step completes a held column.
};

/** @brief Persistent measure_reach() spans of the cached warp field. */
struct Spans {
  int16_t *reach = nullptr; ///< Per row, the reach behind each arm.
  bool valid = false;       ///< Whether they belong to the cached field.
};

} // namespace detail

} // namespace Pixel

} // namespace Filter
