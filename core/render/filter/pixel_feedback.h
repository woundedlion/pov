/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include "core/math/spherical_field.h"
#include "core/memory.h"
#include "core/render/filter/feedback_cap_plane.h"
#include "core/render/filter/feedback_column_time.h"
#include "core/render/filter/feedback_style.h"
#include "core/render/filter/feedback_warp_cache.h"
#include "core/render/filter/pipeline.h"
#include "core/render/spacetime.h"

/**
 * @file pixel_feedback.h
 * @brief Filter::Pixel::Feedback: the terminal filter that warps and
 * composites the previous frame.
 */

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct FeedbackWhiteBox;
}
#endif

namespace Filter {

namespace Pixel {

/**
 * @brief Style-aware terminal feedback filter that warps the previous frame.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The Style's spatial warp is computed on a spherical latitude-ring
 * field, then interpolated within and between rings. Terminal: it must be the
 * last Pipeline stage, and Pipeline::begin_frame() must run before the frame's
 * plot() calls, or Pipeline::begin_split_frame() before them and
 * Pipeline::finish_frame() after; flushing last blanks the frame at
 * alpha >= 1.
 * Rows under POLAR_TARGET_SINE interpolate offsets in their pole's cap plane,
 * where equirect longitude offsets blow up as 1/sin(phi). A pole row is one
 * point, composited once from its source a frame earlier.
 * With temporal sampling on, each pixel reads its source at the time elapsed
 * between their lit columns (FeedbackColumnTime). A source behind the target
 * arm's first column reads the stale frame the drawing buffer still holds;
 * one past the arm's last column reads the frame being drawn after its plots,
 * so the pixel is held: it keeps the part of its composite the previous
 * frames supply, the plots blend over it while its coverage tracks their
 * transmittance, and finish() adds the rest beneath that coverage. A frame
 * without the stale frame, or an unsplit one, reads such a source from the
 * previous frame at its age there.
 */
template <int W, int H> class Feedback : public Is2DWithHistory {
  using SphereField = hs::SphericalFieldLayout<W, H>;
  using CapPlane = FeedbackCapPlane<W, H>;
  using CapOffset = typename CapPlane::CapOffset;
  using CapPoint = typename CapPlane::CapPoint;
  using CapCell = typename CapPlane::CapCell;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::FeedbackWhiteBox;
#endif
  using WarpCache = FeedbackWarpCache<W, H>;
  using WarpField = typename WarpCache::Buffers;
  using ColumnTime = FeedbackColumnTime<W>;
  using FadeTable = FeedbackFadeTable;
  static constexpr int HALF = ColumnTime::HALF;
  /** @brief A colour transform: one sample and its fade-table level to a
   *  pixel. */
  using PixelTransform = FunctionRef<::Pixel(float, float, float, int)>;

  /** @brief Constructs the shared spherical field layout.
   * @param ds Field downsample factor.
   * @return Matching north/south infill and rectangular column count. */
  static HS_HOT_INLINE constexpr SphereField field_for(int ds) {
    return SphereField(ds, ds, ds, W / ds);
  }

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's). Other values render uncached. */
  static constexpr int CACHE_DOWNSAMPLE = ::Feedback::Style{}.downsample;
  static constexpr int CACHE_COLUMNS = W / CACHE_DOWNSAMPLE;
  static constexpr SphereField CACHE_FIELD = field_for(CACHE_DOWNSAMPLE);
  /** @brief Cell count of the cached spherical warp field. */
  static constexpr int CACHE_CELLS = CACHE_COLUMNS * CACHE_FIELD.ring_count();

public:
  /// Pixel domain: runs after every screen stage.
  static constexpr int domain_rank = IsPixel::domain_rank;
  /** @brief Marks this as terminal: flush() writes the Canvas directly. */
  static constexpr bool is_terminal = true;
  /** @brief Latitude sine under which lattice rings carry cap-plane offsets
   *  and rows between two such rings convert each pixel's own target. */
  static constexpr float POLAR_TARGET_SINE = 0.5f;
  /** @brief Opaque store owns the frame: no history stage may precede it. */
  static constexpr bool terminal_replaces = true;

  // Covers only the default-constructed Style; a runtime-swapped style's
  // downsample is checked per flush.
  static_assert(
      ::Feedback::Style{}.downsample > 0 &&
          W % ::Feedback::Style{}.downsample == 0,
      "Feedback<W,H>: default style downsample must be > 0 and divide "
      "W");

  /**
   * @brief Binds the filter to a live feedback Style.
   * @param style Style supplying the spatial warp and color transforms.
   */
  explicit Feedback(::Feedback::Style &style) : feedback_style(&style) {}

  /**
   * @brief Passes a current-frame tap to the next stage, or blends it into a
   *        held pixel of a split frame with its coverage.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  HS_HOT_INLINE void plot(float x, float y, const ::Pixel &color, float age,
                          float alpha, PassFnT &&pass);

  /**
   * @brief Enables or disables feedback.
   * @param value When false, flush() only clears a stale drawing buffer.
   */
  void set_enabled(bool value) { enabled = value; }

  /**
   * @brief Scratch bytes for a full-width uncached flush at downsample ds.
   * @param ds Downsample factor.
   * @return Scratch byte count.
   */
  static size_t uncached_scratch_bytes(int ds);

  /**
   * @brief Allocates the warp-field cache and fade table from the persistent
   * arena.
   * @param arena Persistent arena supplying STORAGE_BYTES bytes.
   * @details Call from effect init(), not the constructor, and again after
   * any arena reset. Without storage every flush needs
   * uncached_scratch_bytes(downsample) scratch bytes, and with temporal
   * sampling on FeedbackFadeTable::BYTES more.
   */
  HS_FLASH_MEMBER void init_storage(Arena &arena);

  /**
   * @brief Accesses the bound Style.
   * @return Mutable reference to the bound feedback Style.
   */
  ::Feedback::Style &style() { return *feedback_style; }
  /**
   * @brief Accesses the bound Style.
   * @return Const reference to the bound feedback Style.
   */
  const ::Feedback::Style &style() const { return *feedback_style; }

private:
  template <int, int, typename...> friend struct ::Pipeline;

  /**
   * @brief Blends the distorted previous frame into the current frame.
   * @param cv Target canvas (reads cv.prev, writes the back (current-draw)
   * buffer).
   * @param alpha Global blend alpha in [0, 1].
   * @details Computes a coarse warp field via the Style's space_fn, bilinearly
   * upsamples it, then composites the warped previous frame, honoring the
   * segment's cylindrical clip. Clears a stale drawing buffer unless column
   * time reads it as the frame two back; a disabled filter only clears it.
   */
  HS_O3_BEGIN
  void flush(Canvas &cv, float alpha);
  HS_O3_END

  /**
   * @brief Clears the frame a disabled filter leaves to the canvas.
   * @param cv Target canvas.
   * @return Whether the filter is disabled.
   */
  bool disabled(Canvas &cv);

  /**
   * @brief Composites a split frame short of the pixels whose sources lie
   *        past their arm's last column, and holds the frame for finish().
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @details The held frame keeps its scratch_arena_a allocations until
   * finish(); scopes opened meanwhile must close first. A frame no pixel of
   * which reads past a cut is released at once.
   */
  HS_FLASH_MEMBER void flush_split(Canvas &cv, float alpha);

  /**
   * @brief Adds the part of each held pixel's composite that reads the frame
   *        as plotted, beneath the coverage its taps left there.
   * @details No-op without a held frame.
   */
  HS_FLASH_MEMBER void finish(Canvas &);
  HS_O3_BEGIN

private:
  /** @brief The lattice rings carrying 3D targets: rings [0, north_rings) and
   *  [south_ring, field_rows), with the sample index where each run ends or
   *  starts. */
  struct PolarRings {
    int north_rings;
    int north_samples;
    int south_ring;
    int south_sample;
    int samples;
    int ring_count;

    constexpr int rows() const { return north_rings + ring_count - south_ring; }
    /** @brief Compacted row of a polar ring. */
    int row(int field_y) const {
      return field_y < north_rings ? field_y
                                   : north_rings + field_y - south_ring;
    }
    /** @brief Compacted slot of a polar ring's lattice sample. */
    int slot(int sample) const {
      return sample < north_samples ? sample
                                    : north_samples + sample - south_sample;
    }
  };

  struct CoarseGrid {
    int downsample;
    SphereField field;
    int columns;
    int field_rows;
    PolarRings polar;
  };

  struct RenderBand {
    int y_begin;
    int y_end;
    int field_y_begin;
    int field_y_end;
    ClipRegion::XClip x_clip;
    std::bitset<W> coarse_columns_used;
  };

  struct WarpControl {
    int16_t x;
    int16_t y;
  };

  struct PixelAccumulator {
    uint32_t r = 0;
    uint32_t g = 0;
    uint32_t b = 0;

    void add(const ::Pixel &pixel) {
      r += pixel.r;
      g += pixel.g;
      b += pixel.b;
    }

    void remove(const ::Pixel &pixel) {
      r -= pixel.r;
      g -= pixel.g;
      b -= pixel.b;
    }

    ::Pixel average(int width) const {
      const uint32_t round = static_cast<uint32_t>(width / 2);
      return ::Pixel(static_cast<uint16_t>((r + round) / width),
                     static_cast<uint16_t>((g + round) / width),
                     static_cast<uint16_t>((b + round) / width));
    }
  };

  struct ColumnRun {
    int begin;
    int end;
  };

  struct ColumnRuns {
    ColumnRun items[6];
    int count;
  };

  enum class FlushMode : uint8_t { SKIP, HUE, PLAIN, GENERAL };

  using Side = detail::Side;

  /** @brief One bilinear sample on the [0, 65535] scale, unquantized. */
  struct Sample {
    float r;
    float g;
    float b;
  };

  /** @brief A lane's resolved read. */
  struct Lane {
    typename ColumnTime::Read read;
    Side side;
  };

  /** @brief A row's interpolation and reconstruction, resolved per frame. */
  struct RowGeometry {
    int16_t field_y0; ///< Ring at or above the row.
    int16_t field_y1; ///< Ring below it, or the same ring on the last.
    float wy1;        ///< Weight of field_y1.
    bool polar : 1;   ///< Both rings carry cap-plane offsets.
    bool north : 1;   ///< A polar row's pole is the north pole.
    bool pairs : 1;   ///< Composites column pairs.
    bool infill : 1;  ///< Reconstructs at its longitude footprint.
    bool pole : 1;    ///< Every column is the pole.
    bool wide : 1;    ///< Covers every column of each arm (cover_slot()).
    /// Lanes past a column that its reconstruction reads: the next pair, or
    /// the infill filter's radius.
    uint8_t reach;
    /// First of the row's finish samples (finish_window()).
    int16_t window_base;
    /// First of the row's coverage slots.
    uint16_t cover_base;
  };

  using AheadSpan = detail::AheadSpan;

  /** @brief Per-frame composite inputs prepare_flush() resolves. */
  struct FlushContext {
    CoarseGrid grid;
    RenderBand band{};
    /// The warp field, each node at its column's lit time.
    WarpField warp{};
    /// The field at each arm's end, which the right nodes of the arm's last
    /// cells read; @c warp itself when the field holds still.
    WarpField end{};
    ::Pixel *filtered_row = nullptr;
    const ::Pixel *previous = nullptr;
    ::Pixel *current = nullptr;
    ::Pixel poles[SphereField::POLE_STORAGE_COUNT]{};
    /// The drawing buffer's stale frame, two frames back, which holds the
    /// content behind each arm's backward cut; null when it was cleared.
    const ::Pixel *back = nullptr;
    /// Shared pole values of the stale frame.
    ::Pixel back_poles[SphereField::POLE_STORAGE_COUNT]{};
    /// Per row, the even column count from each arm's first whose lanes may
    /// read behind its backward cut. Null without @c back.
    const int16_t *reach = nullptr;
    /// Per row, its interpolation and reconstruction.
    const RowGeometry *rows = nullptr;
    /// Per coverage slot, the transmittance of the taps plotted over the
    /// pixel in its upper 15 bits, zero where the pixel is final, and in its
    /// low bit whether the lane's source reads the frame as plotted. Null
    /// when unsplit.
    uint16_t *cover = nullptr;
    /// Per row, the lanes whose sources read the frame as plotted.
    AheadSpan *ahead_span = nullptr;
    /// Shared pole values of the frame being drawn, set by finish().
    ::Pixel current_poles[SphereField::POLE_STORAGE_COUNT]{};
    ColumnRuns runs{};
    float alpha = 0.0f;
    float fade = 0.0f;
    float pole_half_res = 0.0f;
    bool black_skips_color = false;
    /// Whether pixels read their sources at their elapsed column time.
    bool column_time = false;
    FlushMode mode = FlushMode::SKIP;
    float hue_k[9] = {};
    /// Per-elapsed-time transforms; read only with column_time.
    FadeTable::View fades{};
  };

  __attribute__((always_inline)) CoarseGrid
  make_coarse_grid(const Canvas &cv) const;

  static __attribute__((always_inline)) constexpr PolarRings
  polar_rings(const SphereField &field, int rings) {
    PolarRings polar{0, 0, rings, 0, 0, rings};
    bool leading = true;
    int total = 0;
    auto ring = field.ring(0);
    for (int i = 0; i < rings; ++i, ring = field.next_ring(ring)) {
      const bool is_polar =
          SphereField::latitude_sine(ring.y) < POLAR_TARGET_SINE;
      if (leading && is_polar) {
        polar.north_rings = i + 1;
        polar.north_samples = ring.offset + ring.samples;
      } else {
        leading = false;
        if (!is_polar) {
          polar.south_ring = rings;
        } else if (polar.south_ring == rings) {
          polar.south_ring = i;
          polar.south_sample = ring.offset;
        }
      }
      total = ring.offset + ring.samples;
    }
    if (polar.south_ring == rings)
      polar.south_sample = total;
    polar.samples = polar.north_samples + total - polar.south_sample;
    return polar;
  }

  static __attribute__((always_inline)) RenderBand
  make_render_band(const ClipRegion &clip, const CoarseGrid &grid);

  static __attribute__((always_inline)) ColumnRuns
  make_column_runs(const ClipRegion::XClip &clip, bool split_arms = false);

  /** @brief Builds the cache identity for the current warp and render band.
   * @param band Rows and columns the frame renders.
   * @return Warp cache identity. */
  typename WarpCache::Key warp_cache_key(const RenderBand &band) const;

  /**
   * @brief The frame's populated warp field: the cache's when the Style and
   * band allow it, otherwise a scratch field.
   * @param reused Set when the frame reuses the cached field unchanged.
   * @param prior Receives the field one noise step before the returned one.
   */
  __attribute__((always_inline)) WarpField
  acquire_warp_field(Arena &scratch, const CoarseGrid &grid,
                     const RenderBand &band, bool &reused, WarpField &prior);

  /**
   * @brief Fills @p warp's cells for the band from the Style's warp.
   * @param grid Coarse layout of the frame.
   * @param band Rows and columns the frame renders.
   * @param warp Field to fill.
   * @param controls Scratch for every lattice sample's offset.
   * @param sample_caps Scratch for every polar lattice sample's cap offset.
   * @details Under column time a cell's offsets carry their elapsed scale,
   * which caps a forward column offset at a quarter turn; a pole ring's cap
   * offsets stay unscaled.
   */
  __attribute__((always_inline)) void
  populate_warp_field(const CoarseGrid &grid, const RenderBand &band,
                      const WarpField &warp, WarpControl *controls,
                      CapOffset *sample_caps);

  /**
   * @brief Resolves the frame's warp field and composite inputs.
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @param scratch Arena for the frame's uncached warp field and filter row.
   * @param split Whether the frame's held pixels composite after its plots.
   * @return The composite context; mode SKIP when the previous frame is black.
   */
  HS_HOT_FLASH_MEMBER FlushContext prepare_flush(Canvas &cv, float alpha,
                                                 Arena &scratch,
                                                 bool split = false);

  /**
   * @brief Moves each node of a moving field to its column's lit time: an
   *        arm's first column reads the prior field, its end the current one.
   * @param grid Coarse layout of the frame.
   * @param band Rows and columns the frame renders.
   * @param now The current field.
   * @param prior The field one noise step before it; receives the blend.
   */
  HS_HOT_FLASH_MEMBER static void blend_lit_times(const CoarseGrid &grid,
                                                  const RenderBand &band,
                                                  const WarpField &now,
                                                  const WarpField &prior);

  /**
   * @brief Resolves each row's ring pair, weights and reconstruction.
   * @param ctx Frame context with its grid, band and alpha set.
   * @param rows Receives H entries; rows outside the band stay zero.
   */
  HS_FLASH_MEMBER static void describe_rows(const FlushContext &ctx,
                                            RowGeometry *rows);

  /**
   * @brief Composites the warped previous frame through a color transform.
   * @tparam PAIR_PIXELS Whether to interleave two pixels per step through
   * @p transform_pair.
   * @param ctx Frame context from prepare_flush().
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param transform_pair Maps two sampled RGB triples sharing one level to
   *        two output pixels.
   * @param starts Whether to composite only the rows' lanes short of their
   *        reach, on rows reaching no further than a quarter turn, as
   *        between sweep_deep()'s halves; otherwise every row's lanes from its
   *        reach on.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
  void composite_pixels(const FlushContext &ctx,
                        TransformPixelT &&transform_pixel,
                        TransformPairT &&transform_pair, bool starts = false);

  /**
   * @brief Composites a pole row: one read of the pole's source a frame
   *        earlier, written to every column.
   * @param ctx Frame context from prepare_flush().
   * @param y The pole row.
   * @param g The row's geometry.
   * @param transform_pixel Maps a sample and its fade-table level to a pixel.
   */
  HS_FLASH_MEMBER void composite_pole_row(const FlushContext &ctx, int y,
                                          const RowGeometry &g,
                                          PixelTransform transform_pixel);

  /**
   * @brief Reconstructs a dense infill row at its spherical footprint, under
   *        column time within each arm.
   * @param ctx Frame context from prepare_flush().
   * @param y Row, composited into the current frame or, when @p deferred,
   *        into the filter row.
   * @param deferred Whether the row was composited into the filter row.
   */
  HS_HOT_FLASH_MEMBER void filter_infill_row(const FlushContext &ctx, int y,
                                             bool deferred);

  /**
   * @brief Composites one polar row's lanes from their cap-plane offsets.
   * @details Each lane converts its own target back to field coordinates, one
   * lane per step.
   */
  HS_HOT_FLASH_MEMBER void
  composite_polar_row(const FlushContext &ctx, int xs, int xe, int y,
                      const RowGeometry &g, ::Pixel *output, int stride,
                      float lane_offset, bool defer_filter,
                      PixelTransform transform_pixel);

  /** @brief A polar ring's cap offsets in @p field, indexed by coarse
   *  column. */
  static HS_HOT_INLINE const CapOffset *
  polar_caps(const WarpField &field, const FlushContext &ctx, int field_y) {
    return field.cell_caps + ctx.grid.polar.row(field_y) * ctx.grid.columns;
  }

  /** @brief Whether a coarse column is an arm's first, whose node a cell
   *  ending there reads from the field at the arm's end. */
  static HS_HOT_INLINE bool arm_edge(const FlushContext &ctx, int cx) {
    return cx * ctx.grid.downsample % HALF == 0;
  }

  /** @brief A polar row's cell, its nodes at their lit times. */
  static HS_HOT_INLINE CapCell polar_cell(const FlushContext &ctx,
                                          const RowGeometry &g, int cx0,
                                          int cx1);

  /** @brief A row's angle from its pole, in radians. */
  static HS_HOT_INLINE float cap_angle_of(int y, bool north) {
    const float colatitude =
        SphereField::Geometry::row_to_phi(static_cast<float>(y));
    return north ? colatitude : math::PI_F - colatitude;
  }

  /** @brief A cell's offsets: at its left edge and their change across it,
   *  in pixels. */
  struct CellOffsets {
    float left_x;  ///< Column offset at the left edge.
    float slope_x; ///< Column offset change across the cell.
    float left_y;  ///< Row offset at the left edge.
    float slope_y; ///< Row offset change across the cell.
  };

  /** @brief Resolves a cartesian cell lane from its pixel position.
   * @param x,y Target pixel.
   * @param sub Column within the coarse cell.
   * @param lane_offset Half-resolution column offset.
   * @param inverse_downsample Reciprocal coarse-cell width.
   * @param cell Interpolated cell offsets.
   * @return Resolved source and arm side. */
  static HS_HOT_INLINE Lane cartesian_lane(int x, int y, int sub,
                                           float lane_offset,
                                           float inverse_downsample,
                                           const CellOffsets &cell);

  /**
   * @brief A cell's offsets blended between its two rings, its nodes at
   *        their lit times.
   * @param ctx Frame context.
   * @param row0,row1 First cell index of each ring.
   * @param cx0,cx1 The cell's left and right coarse columns.
   * @param wy0,wy1 The rings' weights.
   * @param reach Receives the largest node column offset, WARP_SCALE units.
   * @return The cell's offsets.
   */
  static HS_HOT_INLINE CellOffsets cell_offsets(const FlushContext &ctx,
                                                int row0, int row1, int cx0,
                                                int cx1, float wy0, float wy1,
                                                int &reach);

  /**
   * @brief Resolves a lane's read from its target and source.
   * @param target Target column, with its half-resolution lane offset.
   * @param arm First column of the target's arm.
   * @param source_x Source column, on any wrap branch.
   * @param source_y Source row.
   * @return The read and the side of the arm its source lies on.
   * @details The source unwraps to the shortest path, under a half turn ahead
   * and no further behind than keeps its taps past the other arm's lane of
   * the same column.
   */
  static HS_HOT_INLINE Lane resolve_lane(float target, int arm, float source_x,
                                         float source_y);

  /**
   * @brief Samples a resolved lane from the buffers holding its taps.
   * @param ctx Frame context.
   * @param arm First column of the target's arm.
   * @param lane The lane; AHEAD only in a frame without coverage. Its elapsed
   *        follows the buffer read.
   * @return The bilinear sample.
   * @details A source across a cut reads the frame holding it when the frame
   * has one, else the previous frame, a frame older past the cut and a frame
   * younger behind it.
   */
  HS_HOT_INLINE Sample sample_lane(const FlushContext &ctx, int arm,
                                   Lane &lane) const;

  /**
   * @brief Marks a lane whose source reads the frame as plotted, and opens
   *        the coverage of every column whose reconstruction reads it.
   * @param ctx Frame context with its coverage allocated.
   * @param y Row.
   * @param arm First column of the lane's arm.
   * @param x The lane's column.
   */
  HS_HOT_FLASH_MEMBER void hold_lane(const FlushContext &ctx, int y, int arm,
                                     int x);

  /** @brief Low coverage bit: the lane's source reads the frame as plotted. */
  static constexpr uint16_t AHEAD_BIT = 1;
  /** @brief The coverage bits holding a pixel's transmittance. */
  static constexpr uint16_t TRANSMITTANCE_BITS = 0xFFFE;

  /**
   * @brief Composites a cell's pixels whose sources may lie past their arm's
   * columns, each resolved by FeedbackColumnTime.
   * @param ctx Frame context from prepare_flush().
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param output Row being written.
   * @param x First column.
   * @param end One past the last column, within the cell.
   * @param y Row.
   * @param sub @p x's column within the cell.
   * @param stride Column step.
   * @param lane_offset Half-resolution lane offset of each column.
   * @param cell The cell's elapsed-scaled offsets.
   * @param plain_store Whether pixels overwrite rather than blend.
   * @return The first column past the cell's pixels.
   * @details The colour transforms stay out of line, in the composite's own
   * placement.
   */
  HS_HOT_FLASH_MEMBER int composite_resolved_cell(
      const FlushContext &ctx, PixelTransform transform_pixel, ::Pixel *output,
      int x, int end, int y, int sub, int stride, float lane_offset,
      const CellOffsets &cell, bool plain_store);

  /** @brief A row's state in a column-major sweep: its pole angle and each
   *  arm's current cell. */
  struct SweepRow {
    float cap_angle;    ///< A polar row's angle from its pole.
    int16_t cell_cx[2]; ///< Each arm's current cell's coarse column, or -1.
    union {
      CellOffsets offsets; ///< A cartesian row's offsets.
      CapCell cap;         ///< A polar row's cap-plane offsets.
    } cell[2];             ///< Each arm's current cell.
  };

  /**
   * @brief Each row's sweep state, no cell decoded yet.
   * @param ctx Frame context.
   * @param arena Arena for the H states.
   * @return The states, indexed by row.
   */
  static SweepRow *sweep_rows(const FlushContext &ctx, Arena &arena);

  /**
   * @brief One lane's read in a column-major sweep.
   * @param ctx Frame context.
   * @param g The lane's row geometry.
   * @param y The lane's row.
   * @param x The lane's column.
   * @param row The row's sweep state; decodes the lane's cell when the lane
   *        enters a new one.
   * @return The resolved lane.
   */
  HS_HOT_INLINE Lane sweep_lane(const FlushContext &ctx, const RowGeometry &g,
                                int y, int x, SweepRow &row) const;

  /**
   * @brief Composites the lanes of the rows reaching past a quarter turn
   *        behind their arm's backward cut, column by column.
   * @param ctx Frame context with its stale frame.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param begin,end Columns past each arm's first to sweep.
   * @details A read behind the cut lands past the other arm's lane of the
   * same column, so lane c of a row may be composited once no row's lanes
   * past c that it may read are. These rows' first quarter goes before and
   * their second after the other rows' lanes short of their reach, which
   * never write past a quarter turn and whose reads land past it.
   */
  HS_HOT_FLASH_MEMBER void sweep_deep(const FlushContext &ctx,
                                      PixelTransform transform_pixel, int begin,
                                      int end);

  /**
   * @brief Composites one column of a row in sweep_deep().
   * @param ctx Frame context with its stale frame.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param y Row.
   * @param c Column past each arm's first.
   * @param row The row's sweep state.
   */
  HS_HOT_FLASH_MEMBER void behind_lane(const FlushContext &ctx,
                                       PixelTransform transform_pixel, int y,
                                       int c, SweepRow &row);

  /**
   * @brief Adds every held lane's composite from the frame as plotted beneath
   *        its coverage, in lit-time order so each read lands on columns
   *        already final.
   * @param ctx Held frame context.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @details A read stops at the last final column of its source rows, at
   * most the reconstruction's reach and one column short of the antipode.
   * Each row keeps its last reconstruction window of samples, so a completed
   * column's addition rounds once.
   */
  HS_HOT_FLASH_MEMBER void sweep_ahead(const FlushContext &ctx,
                                       PixelTransform transform_pixel);

  /**
   * @brief Adds the columns a finish sample completes beneath their coverage.
   * @param ctx Held frame context.
   * @param y Row.
   * @param g The row's geometry.
   * @param arm First column of the lane's arm.
   * @param c The lane's column past its arm's first.
   * @param sample The lane's finish composite; black for an unheld lane.
   * @param window The row and arm's recent samples, finish_window() of them.
   */
  HS_HOT_FLASH_MEMBER static void
  complete_columns(const FlushContext &ctx, int y, const RowGeometry &g,
                   int arm, int c, const ::Pixel &sample, ::Pixel *window);

  /**
   * @brief Adds a column's finish composite beneath its coverage.
   * @param ctx Held frame context.
   * @param y Row.
   * @param arm First column of the column's arm.
   * @param x Column.
   * @param under The part of the column's composite the frame as plotted
   *        supplies.
   */
  static HS_HOT_INLINE void add_under(const FlushContext &ctx, int y, int arm,
                                      int x, const ::Pixel &under);

  /**
   * @brief Composites every pixel from the previous frames: the lanes that may
   *        read behind a cut first, in lit-time order, then each row.
   * @param cv Target canvas.
   * @param ctx Frame context from prepare_flush().
   */
  HS_HOT_INLINE void composite_frame(Canvas &cv, FlushContext &ctx);

  /**
   * @brief Composites the frame, or finishes its held pixels, through the
   *        frame's color transform.
   * @param ctx Frame context from prepare_flush().
   * @param finish Whether to add the held pixels' part from the frame as
   *        plotted instead of compositing from the previous frames.
   */
  void composite_pass(const FlushContext &ctx, bool finish) {
    switch (ctx.mode) {
    case FlushMode::HUE: {
      const float *k = ctx.column_time ? ctx.fades.hue_k : ctx.hue_k;
      composite_phase<true>(
          ctx, finish,
          [&](float r, float g, float b, int level) {
            return ::Feedback::hue_fade_apply(k + level * FadeTable::MATRIX, r,
                                              g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              int level, ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k + level * FadeTable::MATRIX, r0, g0,
                                        b0, r1, g1, b1, p0, p1);
          });
      break;
    }
    case FlushMode::PLAIN:
      composite_plain(ctx, finish);
      break;
    case FlushMode::GENERAL:
      composite_general(ctx, finish);
      break;
    case FlushMode::SKIP:
      break;
    }
  }

  /** @brief composite_pass() through a colour transform. */
  template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
  HS_HOT_INLINE void composite_phase(const FlushContext &ctx, bool finish,
                                     TransformPixelT &&transform_pixel,
                                     TransformPairT &&transform_pair);

  /** @brief Composites through the unrotated hue fade: a per-channel scale. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_plain(const FlushContext &ctx, bool finish);

  /** @brief Composites through the Style's type-erased color transform. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_general(const FlushContext &ctx, bool finish);
  HS_O3_END

  // Round-to-nearest fade stalls dim channels short of black.
  static constexpr float NEAR_BLACK = 64.0f;
  static constexpr float WARP_SCALE = 128.0f;
  /** @brief Column offsets in WARP_SCALE units, one full turn apart. */
  static constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;

  /**
   * @brief Whether a column of a half-resolution row stands for its pair.
   * @param ctx Frame context.
   * @param g The row's geometry.
   * @param x Column.
   * @return True except, under column time, for each arm's last two columns,
   *         which are their own lanes.
   */
  static HS_HOT_INLINE bool pair_lane(const FlushContext &ctx,
                                      const RowGeometry &g, int x) {
    return g.pairs &&
           !(ctx.column_time && x - ColumnTime::arm_start(x) >= HALF - 2);
  }

  /**
   * @brief A row's lanes for the row composite.
   * @param ctx Frame context.
   * @param g The row's geometry.
   * @param y Row.
   * @param starts Whether the lanes short of the row's reach are wanted, on a
   *        row reaching no further than a quarter turn.
   * @return The band's runs, or with a stale frame each arm's columns from
   *         the row's reach to its end; a half-resolution row's runs under
   *         column time end short of each arm's last two columns, which form
   *         runs of their own.
   */
  static ColumnRuns row_lanes(const FlushContext &ctx, const RowGeometry &g,
                              int y, bool starts);

  /**
   * @brief Each row's reach behind the arms' backward cuts.
   * @param grid Coarse layout of the frame.
   * @param warp The frame's populated warp field.
   * @param end The field at each arm's end.
   * @param rows The frame's row geometry.
   * @param reach Receives, per row, the even column count from each arm's
   *        first whose lanes' sources may lie behind its backward cut; zero
   *        on a pole row.
   * @details A cell's pixels read no further than its nodes' extreme column
   * offsets in either field; a polar lane's source turns no further than its
   * cell's longest cap offset subtends, and a cap offset reaching the pole may
   * turn it anywhere.
   */
  HS_FLASH_MEMBER void measure_reach(const CoarseGrid &grid,
                                     const WarpField &warp,
                                     const WarpField &end,
                                     const RowGeometry *rows,
                                     int16_t *reach) const;

  /** @brief Latitude sine under which two columns subtend less than a
   *         number of row pitches.
   *  @param pole_half_res Row pitches (Feedback::Style::pole_half_res).
   *  @return The sine; rows below it composite column pairs. */
  static __attribute__((always_inline)) float
  half_res_limit(float pole_half_res) {
    return pole_half_res * W * SphereField::Geometry::RADIANS_PER_ROW *
           (1.0f / (4.0f * math::PI_F));
  }

  /** @brief Whether an unclipped row is an infill row reconstructed at its
   *         spherical longitude footprint.
   *  @param grid Coarse layout of the frame.
   *  @param y Row.
   *  @return True for a row between a pole and the first lattice ring whose
   *          longitude filter spans more than one column. */
  static __attribute__((always_inline)) bool
  filtered_infill(const CoarseGrid &grid, int y);

  /** @brief Whether every column of a row is its pole. */
  static __attribute__((always_inline)) bool is_pole_row(int y) {
    return (SphereField::HAS_NORTH_POLE && y == 0) ||
           (SphereField::HAS_SOUTH_POLE && y == H - 1);
  }

  /** @brief Three quarters of one pixel and a quarter of another, rounded. */
  static ::Pixel blend_quarter(const ::Pixel &own, const ::Pixel &other) {
    return ::Pixel(static_cast<uint16_t>((3u * own.r + other.r + 2u) >> 2),
                   static_cast<uint16_t>((3u * own.g + other.g + 2u) >> 2),
                   static_cast<uint16_t>((3u * own.b + other.b + 2u) >> 2));
  }

  /** @brief Rounded mean of two pixels. */
  static ::Pixel pair_average(const ::Pixel &a, const ::Pixel &b) {
    return ::Pixel(static_cast<uint16_t>((a.r + b.r + 1u) >> 1),
                   static_cast<uint16_t>((a.g + b.g + 1u) >> 1),
                   static_cast<uint16_t>((a.b + b.b + 1u) >> 1));
  }

  /**
   * @brief Expands a half-resolution run's pair samples back to every column.
   * @param output Row whose even columns in [begin, end) hold their pair's
   *        box average, sampled at the pair's midpoint.
   * @param begin First column of the run.
   * @param end One past the run's last column.
   * @param before The pair sample beside @p begin, or null.
   * @param after The pair sample beside @p end, or null.
   * @details Each column takes three quarters of its own pair's sample and a
   * quarter of the neighbouring pair's on its side. Given neighbours stand in
   * at the run's ends; otherwise a full row wraps across the seam and a
   * clipped run repeats its end pairs.
   */
  HS_HOT_FLASH_MEMBER static void
  reconstruct_half_res_run(::Pixel *output, int begin, int end,
                           const ::Pixel *before = nullptr,
                           const ::Pixel *after = nullptr);

  /**
   * @brief Expands a half-resolution row's pair samples back to every column.
   * @param ctx Frame context.
   * @param output Row of pair samples, as for reconstruct_half_res_run().
   * @param y Row index.
   * @details Under column time a full row expands arm by arm short of each
   * arm's last two columns, which are their own lanes: the pair before an
   * arm's first column is the other arm's last two columns of the previous
   * frame, the pair after its last pair its own last two columns. Clipped rows
   * expand each run short of those lanes.
   */
  HS_HOT_FLASH_MEMBER static void expand_half_res_row(const FlushContext &ctx,
                                                      ::Pixel *output, int y);

  /**
   * @brief Lifts @p v onto the wrap branch nearest @p reference.
   * @param v Column offset to unwrap.
   * @param reference Offset whose branch the result must share.
   * @param period Full turn in @p v's units.
   * @return @p v shifted by at most one period, within half a period of
   *         @p reference.
   */
  static __attribute__((always_inline)) float
  unwrap_near(float v, float reference, float period) {
    const float half = period * 0.5f;
    const float delta = v - reference;
    return delta > half ? v - period : (delta < -half ? v + period : v);
  }

  /**
   * @brief unwrap_near() for stored column offsets, in WARP_SCALE units.
   * @param v Stored column offset to unwrap.
   * @param reference Stored offset whose branch the result must share.
   * @return @p v shifted by at most one period, within half a period of
   *         @p reference; identical to the float form on these integers.
   */
  static __attribute__((always_inline)) int unwrap_near(int v, int reference);
  static_assert(static_cast<float>(static_cast<int>(WRAP_PERIOD)) ==
                        WRAP_PERIOD &&
                    static_cast<int>(WRAP_PERIOD) % 2 == 0,
                "Feedback<W,H>: the integer unwrap needs an even, exact "
                "period");

  /**
   * @brief Projected field coordinates of one lattice sample.
   * @param field Spherical layout the lattice belongs to.
   * @param position Lattice position as populate() hands it out.
   * @param point Exact field coordinates of the same sample.
   * @return @p point on a pole row, where the projection is degenerate,
   *         otherwise the projection of @p position, so the offset subtraction
   *         cancels the projection's approximation.
   */
  static typename SphereField::Coordinates
  lattice_origin(const SphereField &field, const math::Vector &position,
                 const typename SphereField::Coordinates &point) {
    const bool pole_row = (SphereField::HAS_NORTH_POLE && point.y == 0.0f) ||
                          (SphereField::HAS_SOUTH_POLE && point.y == H - 1);
    return pole_row ? point : field.project(position);
  }
  // populate_warp_field() canonicalizes the stored column offset onto
  // [-W*WARP_SCALE/2, W*WARP_SCALE/2] and clamps the elapsed-scaled offset
  // to that range before the int16_t cast.
  static_assert(W * WARP_SCALE * 0.5f <= 32767.0f,
                "Feedback<W,H>: canonical warp offset must fit int16_t");
  // Under runtime geometry, row offsets saturate at +/-32767 WARP_SCALE units.
#if !HS_RUNTIME_DISPLAY_GEOMETRY
  static_assert(
      math::PI_F * math::ROWS_PER_RADIAN<H> * WARP_SCALE <= 32767.0f,
      "Feedback<W,H>: fixed-geometry warp row offset must fit int16_t");
#endif

  /**
   * @brief Tests whether the previous frame has any non-black pixel.
   * @param cv Canvas whose previous-frame buffer is scanned.
   * @return True on the first lit pixel found, false if the frame is all black.
   * @details Scans the clip band, which covers the canvas for this full-frame stage.
   */
  static bool any_pixel_lit(const Canvas &cv);

  /**
   * @brief Chooses one color for a longitude-aliased pole row.
   * @param pole_row Base of the pole row in the previous frame.
   */
  HS_O3_FN static ::Pixel select_pole_sample(const ::Pixel *pole_row);

  /**
   * @brief Chooses each pole row's color in a frame.
   * @param frame Frame buffer.
   * @param poles Receives the shared value of each pole row present.
   */
  static void select_poles(const ::Pixel *frame,
                           ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT]);

  /**
   * @brief Samples a read across one of its target arm's cuts.
   * @param ctx Frame context.
   * @param arm First column of the target's arm.
   * @param x Source column, unwrapped near the target.
   * @param y Source row.
   * @return The bilinear sample, the taps on each side of the cut from the
   *         buffer holding that side: the stale frame behind the arm, the
   *         frame being drawn past it, the previous frame within it.
   */
  HS_HOT_FLASH_MEMBER Sample sample_across(const FlushContext &ctx, int arm,
                                           float x, float y) const;

  /** @brief Whether a near-black sample skips its color transform.
   * @param sample Unquantized linear sample.
   * @param enabled Whether near-black samples skip the color transform.
   * @return True when enabled and all channels are below NEAR_BLACK. */
  static HS_O3_FN HS_HOT_INLINE bool skip_near_black(const Sample &sample,
                                                     bool enabled) {
    return enabled && sample.r < NEAR_BLACK && sample.g < NEAR_BLACK &&
           sample.b < NEAR_BLACK;
  }

  /**
   * @brief Bilinearly samples a frame buffer.
   * @param field Spherical topology and interpolation policy.
   * @param prev Base of the frame buffer sampled, row-major with stride W.
   * @param poles Shared values for every aliased column of the pole rows.
   * @param bx Fractional column in [-W, 2W).
   * @param by Fractional row; north crossings reflect with a half-turn.
   * @return The interpolated channels, returned in registers.
   */
  HS_O3_FN __attribute__((noinline)) Sample
  sample_bilinear_prev(const SphereField &field, const ::Pixel *prev,
                       const ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT],
                       float bx, float by) const;

  /**
   * @brief Calls @p transform out of line.
   * @param transform Single-pixel colour transform.
   * @param r Sampled red.
   * @param g Sampled green.
   * @param b Sampled blue.
   * @param level Fade-table level.
   * @return The transformed pixel.
   */
  template <typename TransformPixelT>
  __attribute__((noinline)) static ::Pixel
  transform_pixel_call(TransformPixelT &transform, float r, float g, float b,
                       int level) {
    return transform(r, g, b, level);
  }

  /**
   * @brief Transforms a sample, out of line when the composite pairs pixels.
   * @tparam PAIR_PIXELS Whether the composite pairs pixels.
   * @param transform Single-pixel colour transform.
   * @param s Sample.
   * @param level Fade-table level.
   * @return The transformed pixel.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT>
  static HS_HOT_INLINE ::Pixel transform_sample(TransformPixelT &transform,
                                                const Sample &s, int level) {
    if constexpr (PAIR_PIXELS)
      return transform_pixel_call(transform, s.r, s.g, s.b, level);
    else
      return transform(s.r, s.g, s.b, level);
  }

  /** @brief Calls @p transform out of line. */
  template <typename TransformPairT>
  __attribute__((noinline)) static void
  transform_pair_call(TransformPairT &transform, float r0, float g0, float b0,
                      float r1, float g1, float b1, int level, ::Pixel &p0,
                      ::Pixel &p1) {
    transform(r0, g0, b0, r1, g1, b1, level, p0, p1);
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a ::Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static HS_HOT_INLINE uint16_t quantize16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 65535.0f) + 0.5f);
  }

  /**
   * @brief The bound generator's own configuration, as a cache-key scalar.
   * @details Only the seed is keyed. Fractal and rotation settings on
   * NoiseParams::noise must remain at defaults, or the persistent arena must be
   * reset before re-running init_storage() after changes.
   */
  static uint32_t noise_config_key(const Animation::NoiseParams *noise) {
    return noise ? static_cast<uint32_t>(noise->seed) : 0u;
  }

  /**
   * @brief Blends a tap near an arm's end into a covered pixel, or passes it
   *        on.
   * @param x Column coordinate as plotted.
   * @param column @p x rounded and wrapped as the canvas writes it.
   * @param y Row coordinate.
   * @param color,age,alpha Shading registers.
   * @param pass Downstream 2D callback, the canvas write.
   * @details The color blends as the canvas write blends it, and the
   * transmittance takes the same blend toward zero, so finish() can add the
   * rest of the pixel's composite beneath both.
   */
  template <typename PassFnT>
  HS_HOT_FLASH_MEMBER void plot_near_arm_end(float x, int column, float y,
                                             const ::Pixel &color, float age,
                                             float alpha, PassFnT &pass);

  /** @brief Rows and lattice rings under POLAR_TARGET_SINE, at most, over
   *  the display geometries math::set_display_geometry() admits: its ends
   *  sampled, a row's angle taken a milliradian toward the cap and one row
   *  added at each end; describe_rows() and init_storage() check the bounds. */
  struct CapExtent {
    int rows;
    int rings;
  };
  static constexpr CapExtent CAP_EXTENT = [] {
#if HS_RUNTIME_DISPLAY_GEOMETRY
    constexpr float CAP = math::PI_F / 6.0f + 1.0e-3f;
    constexpr int STEPS = 24;
    CapExtent most{0, 0};
    for (int i = 0; i <= STEPS; ++i)
      for (int j = 0; j <= STEPS; ++j) {
        const float north = 0.25f * math::PI_F * i / STEPS;
        const float south = math::PI_F * (0.75f + 0.25f * j / STEPS);
        const float pitch = (south - north) / (H - 1);
        int top = 0, bottom = 0;
        for (int y = 0; y < H; ++y) {
          const float phi = north + y * pitch;
          top += phi < CAP;
          bottom += phi > math::PI_F - CAP;
        }
        top = std::min(top + 1, H);
        bottom = std::min(bottom + 1, H);
        const int rings = CACHE_FIELD.rings_within(top, bottom);
        most.rows = std::max(most.rows, std::min(top + bottom, H));
        most.rings = std::max(most.rings, rings);
      }
    return most;
#else
    int rows = 0;
    for (int y = 0; y < H; ++y)
      rows += SphereField::latitude_sine(y) < POLAR_TARGET_SINE;
    return CapExtent{rows,
                     polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows()};
#endif
  }();

  /** @brief Coverage slots per arm of a row: every column of a row under
   *  POLAR_TARGET_SINE, where a cap-plane lane may turn any distance;
   *  otherwise the columns a quarter-turn forward offset, the elapsed
   *  scale's cap, carries past the arm's last column, and the pair column
   *  before them. */
  static constexpr int COVER_CARTESIAN = HALF / 2 + 1;
  /** @brief A row's coverage slots per arm. */
  static HS_HOT_INLINE int cover_columns(const RowGeometry &g) {
    return g.wide ? HALF : COVER_CARTESIAN;
  }
  /** @brief Coverage slots of a frame. */
  static constexpr int COVER_SLOTS =
      2 * (CAP_EXTENT.rows * HALF + (H - CAP_EXTENT.rows) * COVER_CARTESIAN);
  static_assert(COVER_SLOTS <= UINT16_MAX,
                "Feedback<W,H>: RowGeometry::cover_base holds a slot");
  /**
   * @brief Coverage slot of a column.
   * @param g The column's row geometry.
   * @param arm First column of the column's arm.
   * @param x Column.
   * @return Index into FlushContext::cover, or -1 for a column too far from
   *         its arm's end to have one.
   */
  static HS_HOT_INLINE int cover_slot(const RowGeometry &g, int arm, int x) {
    const int width = cover_columns(g);
    const int slot = arm + HALF - 1 - x;
    return slot < width ? g.cover_base + (arm != 0 ? width : 0) + slot : -1;
  }

  /**
   * @brief Bound on an infill row's longitude filter width.
   * @param y Row.
   * @return An odd width no smaller than the layout's.
   */
  static constexpr int filter_window_bound(int y) {
    constexpr int WIDEST = (W & 1) ? W : W - 1;
#if HS_RUNTIME_DISPLAY_GEOMETRY
    // math::set_display_geometry() puts the far end at least 3 pi / 4
    // from each pole, so a row d rows from an end lies at least
    // d * (3 pi / 4) / (H - 1) from its pole; a - a^3 / 6 bounds that angle's
    // sine from below.
    const int rows = std::min(y, H - 1 - y);
    const float angle = rows * (0.75f * math::PI_F / (H - 1));
    const float sine = angle - angle * angle * angle * (1.0f / 6.0f);
#else
    const float sine = SphereField::latitude_sine(y);
#endif
    if (sine < 2.0f / W)
      return WIDEST;
    int width = 1;
    while (width * sine < 1.0f)
      width += 2;
    return std::min(width + 2, WIDEST);
  }

  /** @brief Finish samples a row keeps per arm: the lanes its reconstruction
   *  reaches, at most the arm the infill filter narrows to. */
  static constexpr int finish_window(const RowGeometry &g) {
    return g.infill ? std::min(2 * g.reach + 1, HALF) : g.pairs ? 2 : 0;
  }

  /** @brief Finish window samples a frame may keep: two pair samples per arm on
   *  every row, and the infill rows' filter widths. */
  static constexpr int WINDOW_SLOTS = [] {
    int slots = 0;
    for (int y = 0; y < H; ++y) {
#if HS_RUNTIME_DISPLAY_GEOMETRY
      const bool infill = y < CACHE_DOWNSAMPLE || y >= H - CACHE_DOWNSAMPLE;
#else
      const bool infill =
          (y < CACHE_DOWNSAMPLE && (!SphereField::HAS_NORTH_POLE || y > 0)) ||
          (y >= H - CACHE_DOWNSAMPLE &&
           (!SphereField::HAS_SOUTH_POLE || y < H - 1));
#endif
      slots += 2 * (infill ? std::min(filter_window_bound(y), HALF) : 2);
    }
    return slots;
  }();

  /** @brief Traps while a split frame awaits finish(). */
  HS_FLASH_MEMBER void require_unheld() const {
    HS_CHECK(!held, "feedback: a split frame was never finished");
  }

  /** @brief Releases the held frame and its scratch. */
  HS_FLASH_MEMBER void end_split() {
    held = nullptr;
    held_floor = HALF;
    scratch_arena_a.set_offset(held_mark);
  }

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  /** When false, flush() only clears a stale drawing buffer. */
  bool enabled = true;
  /** Split frame flush_split() holds for finish(), or null. */
  FlushContext *held = nullptr;
  /** scratch_arena_a offset finish() rewinds to. */
  size_t held_mark = 0;
  /** Least held column past an arm's first in the held frame; HALF holds no
   *  column. */
  int held_floor = HALF;
  /** The same per row. */
  uint8_t row_floor[H] = {};
  static_assert(HALF <= UINT8_MAX, "Feedback<W,H>: row_floor holds HALF");
  WarpCache warp_cache; /**< Persistent warp field and lattice origins. */
  using Spans = detail::Spans;
  Spans spans;
  FadeTable fade_table; /**< Persistent per-elapsed-time transforms. */

public:
  /** @brief Cap-cell budget: the most polar rings a display geometry
   *  places (CAP_EXTENT). */
  static constexpr int CACHE_CAP_CELLS = CAP_EXTENT.rings * CACHE_COLUMNS;

  /** @brief Lattice-origin budget: the most samples a display geometry
   *  places. Under runtime geometry, row y lies between y * (3 pi / 4) /
   *  (H - 1) and pi / 4 more from the north pole (math::set_display_geometry()),
   *  and each ring takes the largest sine over that span. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int CACHE_ORIGIN_SAMPLES =
      CACHE_FIELD.sample_count([](int y) {
        constexpr float PITCH = 0.75f * math::PI_F / (H - 1);
        const float low = y * PITCH, high = 0.25f * math::PI_F + y * PITCH;
        if (low <= 0.5f * math::PI_F && high >= 0.5f * math::PI_F)
          return 1.0f;
        return std::min(1.0f, std::max(math::constexpr_sin(low),
                                       math::constexpr_sin(high)) +
                                  1.0e-4f);
      });
#else
  static constexpr int CACHE_ORIGIN_SAMPLES = CACHE_FIELD.sample_count();
#endif

  /** @brief Persistent bytes of one warp field set over CACHE_CELLS: int16
   *  column and row offsets and polar cells' int16 cap-plane offsets. */
  static constexpr size_t FIELD_SET_BYTES =
      CACHE_CELLS * 2 * sizeof(int16_t) + CACHE_CAP_CELLS * 2 * sizeof(int16_t);
  /** @brief Persistent bytes init_storage() reserves over CACHE_CELLS: two
   *  alternating field sets, the lattice's projected origins, the fade table
   *  and the rows' reaches. */
  static constexpr size_t STORAGE_BYTES =
      2 * FIELD_SET_BYTES +
      CACHE_ORIGIN_SAMPLES * sizeof(typename SphereField::Coordinates) +
      FadeTable::BYTES + H * sizeof(int16_t);
  /** @brief Arena allocations init_storage() makes; each pads its share of
   *  STORAGE_BYTES by less than alignof(std::max_align_t). */
  static constexpr size_t STORAGE_ALLOCATIONS = 10;

  /**
   * @brief scratch_arena_a bytes a split frame holds from flush_split() to
   *        finish() with its warp field and fade table in storage.
   * @return The frame context, filter row, row geometry, per-row reaches,
   *         coverage and held-lane spans, each allocation padded.
   */
  static constexpr size_t split_scratch_bytes() {
    constexpr size_t PAD = alignof(std::max_align_t);
    return sizeof(FlushContext) + W * sizeof(::Pixel) +
           H * (sizeof(RowGeometry) + sizeof(int16_t) + sizeof(AheadSpan)) +
           COVER_SLOTS * sizeof(uint16_t) + 6 * PAD;
  }

  /**
   * @brief scratch_arena_a bytes a composite or finish() takes past the
   *        frame it holds, and releases before returning.
   * @pre Style downsample is at most CACHE_DOWNSAMPLE.
   * @return finish()'s sample windows, visited rows and sweep states, each
   *         allocation padded; sweep_deep()'s states are fewer.
   */
  static constexpr size_t transient_scratch_bytes() {
    constexpr size_t PAD = alignof(std::max_align_t);
    return WINDOW_SLOTS * sizeof(::Pixel) +
           H * (sizeof(int16_t) + sizeof(SweepRow)) + 3 * PAD;
  }
};

} // namespace Pixel

} // namespace Filter

#include "core/render/filter/pixel_feedback/frame.h"
#include "core/render/filter/pixel_feedback/rows.h"
#include "core/render/filter/pixel_feedback/sweeps.h"
#include "core/render/filter/pixel_feedback/reconstruction.h"
#include "core/render/filter/pixel_feedback/sampling.h"
