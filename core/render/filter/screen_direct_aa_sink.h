/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/screen_anti_alias.h"
#include "render/canvas.h"
#include "engine/concepts.h"

/**
 * @file screen_direct_aa_sink.h
 * @brief Filter::Screen::DirectAntiAliasSink: terminal four-tap anti-alias sink
 * that writes the framebuffer directly, standing in for a whole Pipeline.
 */

namespace Filter {
namespace Screen {

/**
 * @brief Terminal four-tap anti-alias sink with direct framebuffer writes.
 * @details Opt-in replacement for `Pipeline<W, H, AntiAlias<W, H>>` when no
 * downstream filter is required. It preserves AntiAlias tap ordering and the
 * base sink's q16 source-over blend while resolving rows, columns and clipping
 * once per sample.
 */
HS_O3_BEGIN
template <int W, int H> class DirectAntiAliasSink : public IsPipelineSink {
public:
  static constexpr bool any_crosses_segments = false;
  static constexpr bool any_reads_outside_band = false;
  static constexpr int segment_margin = 1;
  static constexpr int total_segment_margin = segment_margin;
  static constexpr bool any_2d_history = false;
  static constexpr bool any_3d_history = false;
  static constexpr bool any_2d_trail_history = false;
  static constexpr bool any_terminal_history = false;
  static constexpr bool has_world_cull = false;
  static constexpr bool has_world_stage = false;
  static constexpr bool direct_raster_path = true;

  // Pipeline<> derives its stage-level view from the fold; a sink writes both
  // by hand, so pin them to the same identities.
  static_assert(crosses_segments == any_crosses_segments &&
                    reads_outside_band == any_reads_outside_band &&
                    segment_margin == total_segment_margin &&
                    has_history == (any_2d_history || any_3d_history),
                "DirectAntiAliasSink's stage-level traits contradict its "
                "pipeline folds");

  /** @brief Caches the current frame's framebuffer and clip bounds. */
  void prepare(Canvas &cv) {
    base = cv.data();
    const ClipRegion &cr = cv.clip();
    clip_stamp = cr;
    const ClipRegion::XClip xc = cr.x_clip();
    for (int x = 0; x < W; ++x)
      x_visible[x] = !xc.clipped(x);
    const int y_start = cr.render_y_start();
    const int y_end = cr.render_y_end();
    for (int y = 0; y < H; ++y)
      y_visible[y] = y >= y_start && y < y_end;
  }

  /**
   * @brief True when the cached framebuffer base and clip state still match @p cv.
   * @param cv Canvas the caller is about to draw into.
   * @details The base pointer alone repeats every other frame when the canvas
   * double-buffers, so the clip bounds are stamped alongside it.
   */
  bool prepared_for(Canvas &cv) const {
    return base == cv.data() && clip_stamp == cv.clip();
  }

  /** @brief Splats one screen-space sample directly into the Canvas. */
  void plot(Canvas &cv, float x, float y, const ::Pixel &c, float age,
            float alpha) {
    assert(age >= 0.0f && alpha >= 0.0f);
    assert(prepared_for(cv));
    (void)age;
    (void)cv;

    HS_MSP_STALL_START(aa_start);
    const SplatTaps t = splat_taps<W, H>(x, y);

    const bool x0_ok = x_visible[t.x0];
    const bool x1_ok = x_visible[t.x1];
    const bool y0_ok = t.y0_physical && y_visible[t.y0];
    const bool y1_ok = t.y1_physical && y_visible[t.y1];
    const uint8_t tap_mask = static_cast<uint8_t>(
        (y0_ok && x0_ok && t.v00 > SPLAT_TAP_CUTOFF) |
        ((y0_ok && x1_ok && t.v10 > SPLAT_TAP_CUTOFF) << 1) |
        ((y1_ok && x0_ok && t.v01 > SPLAT_TAP_CUTOFF) << 2) |
        ((y1_ok && x1_ok && t.v11 > SPLAT_TAP_CUTOFF) << 3));
    HS_MSP_STALL_STOP(aa_weights, aa_start);

#ifdef HS_PROFILE_MINDSPLATTER_COUNTS
    const unsigned tap_count = static_cast<unsigned>(tap_mask & 0x01) +
                               static_cast<unsigned>((tap_mask >> 1) & 0x01) +
                               static_cast<unsigned>((tap_mask >> 2) & 0x01) +
                               static_cast<unsigned>((tap_mask >> 3) & 0x01);
    ++hs::g_mindsplatter_counts.aa_tap_masks[tap_count];
#endif

    HS_MSP_STALL_START(blend_start);
    if (tap_mask == 0x0f) {
      HS_MSP_COUNT(interior_splats);
      blend_four(base + t.y0 * W, base + t.y1 * W, t.x0, t.x1, c, alpha, t.v00,
                 t.v10, t.v01, t.v11);
      HS_MSP_STALL_STOP(framebuffer_blend, blend_start);
      return;
    }
    HS_MSP_COUNT(clip_boundary_splats);
    blend_masked(base, t.x0, t.x1, t.y0, t.y1, c, alpha, t.v00, t.v10, t.v01,
                 t.v11, tap_mask);
    HS_MSP_STALL_STOP(framebuffer_blend, blend_start);
  }

  /** @brief Integer-coordinate overload matching a filtered Pipeline. */
  void plot(Canvas &cv, int x, int y, const ::Pixel &c, float age,
            float alpha) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha);
  }

  /** @brief Projects a world point, then applies the direct screen-space splat. */
  void plot(Canvas &cv, const Vector &v, const ::Pixel &c, float age,
            float alpha) {
    HS_MSP_STALL_START(projection_start);
    const PixelCoords p = vector_to_pixel<W, H>(v);
    HS_MSP_STALL_STOP(projection, projection_start);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  /**
   * @brief Trail flush (nothing to flush).
   * @tparam TrailFn ScreenTrailFn or WorldTrailFn.
   * @details Dependent-false guard matching the filterless Pipeline: this sink
   * carries no history, so the call would emit nothing.
   */
  template <typename TrailFn> void flush(Canvas &, const TrailFn &, float) {
    static_assert(
        !sizeof(TrailFn *),
        "Wrong flush() domain: DirectAntiAliasSink has no filter stages, so "
        "it carries no history and this overload emits nothing. Drop the "
        "flush() call, or use a Pipeline carrying the history stage it was "
        "meant for.");
  }
  /**
   * @brief Terminal-stage flush (nothing to flush).
   * @tparam T Unused; carries the dependent-false guard.
   */
  template <typename T = void> void flush(Canvas &, float) {
    static_assert(
        !sizeof(T *),
        "Wrong flush() domain: DirectAntiAliasSink has no filter stages, so "
        "it has no terminal stage and this overload emits nothing. Drop the "
        "flush() call, or use a Pipeline carrying Pixel::Feedback.");
  }

  /** @brief Terminal clip-cull predicate forwarding. */
  template <typename Pred>
  bool could_intersect_clip(const Vector &a, const Vector &b,
                            const Basis *planar_basis, Pred &&pred) const {
    return pred(a, b, planar_basis);
  }

private:
  ::Pixel *base = nullptr;
  ClipRegion clip_stamp{};
  // Zero-init: a plot() before prepare() is a masked no-op, not a read of
  // indeterminate bytes through a null base.
  std::array<uint8_t, W> x_visible{};
  std::array<uint8_t, H> y_visible{};

  static __attribute__((always_inline)) uint16_t tap_alpha_q16(float alpha,
                                                               float weight) {
    return static_cast<uint16_t>(
        hs::clamp(alpha * weight * 65535.0f + 0.5f, 0.0f, 65535.0f));
  }

  static __attribute__((always_inline)) void
  blend_four(::Pixel *row0, ::Pixel *row1, int x0, int x1, const ::Pixel &src,
             float alpha, float v00, float v10, float v01, float v11) {
    const uint16_t a00 = tap_alpha_q16(alpha, v00);
    const uint16_t a10 = tap_alpha_q16(alpha, v10);
    const uint16_t a01 = tap_alpha_q16(alpha, v01);
    const uint16_t a11 = tap_alpha_q16(alpha, v11);
    blend_tap(row0 + x0, src, a00);
    blend_tap(row0 + x1, src, a10);
    blend_tap(row1 + x0, src, a01);
    blend_tap(row1 + x1, src, a11);
  }

  HS_NOINLINE_NOCLONE static void blend_masked(::Pixel *base, int x0, int x1,
                                               int y0, int y1,
                                               const ::Pixel &src, float alpha,
                                               float v00, float v10, float v01,
                                               float v11, uint8_t tap_mask) {
    if (tap_mask & 0x01)
      blend_tap(base + y0 * W + x0, src, tap_alpha_q16(alpha, v00));
    if (tap_mask & 0x02)
      blend_tap(base + y0 * W + x1, src, tap_alpha_q16(alpha, v10));
    if (tap_mask & 0x04)
      blend_tap(base + y1 * W + x0, src, tap_alpha_q16(alpha, v01));
    if (tap_mask & 0x08)
      blend_tap(base + y1 * W + x1, src, tap_alpha_q16(alpha, v11));
  }

  static __attribute__((always_inline)) void
  blend_tap(::Pixel *dst, const ::Pixel &src, uint16_t alpha_q16) {
    *dst = dst->lerp16(src, alpha_q16);
  }
};
HS_O3_END

static_assert(PipelineFoldSurface<::Pipeline<8, 8>>);
static_assert(PipelineFoldSurface<::Pipeline<8, 8, AntiAlias<8, 8>>>);
static_assert(PipelineFoldSurface<DirectAntiAliasSink<8, 8>>,
              "DirectAntiAliasSink stands in for a Pipeline<>, so it must "
              "declare every pipeline fold member itself");

} // namespace Screen
} // namespace Filter
