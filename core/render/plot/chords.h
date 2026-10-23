/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file chords.h
 * @brief Strokes polylines that are straight in an azimuthal-equidistant
 *        chart, culled conservatively against the clip band.
 */

#include <array>

#include "core/memory.h"
#include "core/math/display_geometry.h"
#include "core/math/pixel_mapping.h"
#include "core/render/canvas.h"
#include "core/render/clip.h"
#include "core/render/plot/cull.h"
#include "core/render/plot/point_route.h"
#include "core/render/plot/polar_walk.h"
#include "core/render/plot/mesh.h"
#include "core/render/plot/raster.h"

namespace Plot {

/**
 * @brief Screen rows and unwrapped columns a splat must touch to reach a clip
 *        region, and the conservative reach test for chart-straight pieces.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct ClipBand {
  float row_lo = 0.0f;   ///< Lowest row, one above the render region.
  float row_hi = 0.0f;   ///< Exclusive end row of the render region.
  float x_center = 0.0f; ///< Band centre column, unwrapped past W.
  float x_half = 0.0f;   ///< Band half-width in columns.
  bool x_active = false; ///< Whether the clip limits columns.

  /** @brief Rows and columns a splat and the fast-trig pixel map reach past a
   *  bound. */
  static constexpr float PAD = 2.0f;

  /** @brief The band a sample must splat into to touch @p clip's render
   *  region.
   *  @param clip Clip region to cover.
   *  @return The band for @p clip. */
  HS_FLASH_MEMBER static ClipBand of(const ClipRegion &clip) {
    ClipBand band;
    const ClipRegion::XClip x_clip = clip.x_clip();
    band.row_lo = static_cast<float>(clip.render_y_start() - 1);
    band.row_hi = static_cast<float>(clip.render_y_end());
    band.x_active = x_clip.active;
    const float band_start = static_cast<float>(x_clip.rs - 1);
    const float band_end =
        static_cast<float>(x_clip.wrap ? x_clip.re + W : x_clip.re);
    band.x_center = 0.5f * (band_start + band_end);
    band.x_half = 0.5f * (band_end - band_start);
    return band;
  }

  /**
   * @brief Whether a chart-straight piece between two projected points can
   *        splat into the band.
   * @param pa Pixel coordinates of the piece's start.
   * @param pb Pixel coordinates of the piece's end.
   * @param arc The piece's chart length, in radians.
   * @return False only when the piece cannot touch the band.
   * @details The azimuthal-equidistant chart-to-sphere map is 1-Lipschitz, so
   * every point of the piece lies within colatitude (phi_a + phi_b +- arc) / 2
   * and, above the smallest colatitude sine, longitude (lambda_a + lambda_b) /
   * 2 +- arc / (2 sin phi).
   */
  __attribute__((always_inline)) bool may_reach(const math::PixelCoords &pa,
                                                const math::PixelCoords &pb,
                                                float arc) const {
    const float half_rows = 0.5f * arc * math::ROWS_PER_RADIAN<H>;
    const float mid_row = 0.5f * (pa.y + pb.y);
    const float lo = mid_row - half_rows;
    const float hi = mid_row + half_rows;
    if (hi + PAD < row_lo || lo - PAD >= row_hi)
      return false;
    if (!x_active)
      return true;
    const float phi_lo = math::DisplayGeometry<H>::row_to_phi(lo);
    const float phi_hi = math::DisplayGeometry<H>::row_to_phi(hi);
    if (phi_lo <= 0.0f || phi_hi >= math::PI_F)
      return true;
    const float sin_min =
        fminf(math::fast_sinf(phi_lo), math::fast_sinf(phi_hi));
    float span_x = pb.x - pa.x;
    if (span_x > W * 0.5f)
      span_x -= W;
    else if (span_x < -W * 0.5f)
      span_x += W;
    // 0.98 absorbs fast_sinf's error.
    const float half_cols =
        0.5f * arc * (W / (2.0f * math::PI_F)) / (0.98f * sin_min);
    const float centered = pa.x + 0.5f * span_x - x_center + 0.5f * W;
    const float offset =
        fabsf(centered - W * floorf(centered * (1.0f / W)) - 0.5f * W);
    return offset <= half_cols + x_half + PAD;
  }
};

/**
 * @brief Cuts a polyline's chart-straight edges into equal chart pieces
 *        and flags the runs of pieces that cannot reach a clip band, so
 *        Plot::rasterize never simulates them.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Points are added only where a run begins or ends. A run's start
 * depends on the clip, so clipped strokes match the unclipped one to a fraction
 * of a pixel, not bit for bit. Edges touching the chart antipode stay whole.
 */
template <int W, int H> class PlanarBandSplit {
public:
  /** @brief Arena bytes init_storage() takes for up to @p max_points points.
   *  @param max_points Point capacity.
   *  @return Byte count. */
  static constexpr size_t storage_bytes(int max_points) {
    return static_cast<size_t>(max_points) + alignof(uint8_t);
  }

  /** @brief Most points split() emits for @p edges edges of @p pieces each.
   *  @param edges Input edge count.
   *  @param pieces Pieces per edge.
   *  @return Point count. */
  static constexpr int64_t max_points(int edges, int pieces) {
    return static_cast<int64_t>(edges) * pieces + 1;
  }

  /** @brief Binds flag storage for polylines of up to @p max_points points.
   *  @param arena Arena the flags are allocated from.
   *  @param max_points Point capacity; at least 2. */
  HS_COLD_MEMBER void init_storage(Arena &arena, int max_points) {
    HS_CHECK(max_points >= 2, "PlanarBandSplit: max_points %d < 2", max_points);
    capacity = max_points;
    flag_storage = arena.allocate_n<uint8_t>(max_points);
  }

  /**
   * @brief Splits @p ring into @p out and flags each output edge.
   * @param out Empty destination receiving positions only; registers are dropped. Bound for
   * max_points(edges, pieces) points.
   * @param ring Open or closed polyline with @p edges + 1 points.
   * @param edges Edge count.
   * @param pieces Chart pieces per edge.
   * @param planar_basis Azimuthal-equidistant chart the edges are straight in.
   * @param band Clip band the pieces are tested against.
   * @return Rasterize edge flags for @p out, one per edge.
   */
  HS_HOT_FLASH_MEMBER std::span<const uint8_t>
  split(Fragments &out, const Fragments &ring, int edges, int pieces,
        const math::Basis &planar_basis, const ClipBand<W, H> &band) {
    HS_CHECK(out.empty(), "PlanarBandSplit::split: out must be empty");
    HS_CHECK(edges >= 1 && pieces >= 1 &&
                 ring.size() > static_cast<size_t>(edges),
             "PlanarBandSplit: %d edges of %d pieces over %d ring points",
             edges, pieces, static_cast<int>(ring.size()));
    HS_CHECK(max_points(edges, pieces) <= capacity,
             "PlanarBandSplit: %d edges of %d pieces exceed capacity %d", edges,
             pieces, capacity);
    Fragment point;
    point.pos = ring[0].pos;
    out.push_back(point);
    auto push = [&](const math::Vector &position, bool visible)
                    __attribute__((always_inline)) {
                      flag_storage[out.size() - 1] =
                          visible ? RasterOptions::EDGE_VISIBLE : uint8_t{0};
                      point.pos = position;
                      out.push_back(point);
                    };
    const float inv_pieces = 1.0f / static_cast<float>(pieces);
    for (int edge = 0; edge < edges; ++edge) {
      const math::Vector &a = ring[edge].pos;
      const math::Vector &b = ring[edge + 1].pos;
      if (math::dot(a, planar_basis.v) < -COS_PLANAR_ANTIPODE ||
          math::dot(b, planar_basis.v) < -COS_PLANAR_ANTIPODE) {
        push(b, true);
        continue;
      }
      const auto pa = azimuthal_project(a, planar_basis);
      const auto pb = azimuthal_project(b, planar_basis);
      const float dx = pb.first - pa.first;
      const float dy = pb.second - pa.second;
      const float piece_arc = sqrtf(dx * dx + dy * dy) * inv_pieces;
      math::Vector piece_start = a;
      math::PixelCoords piece_start_px = math::vector_to_pixel<W, H>(a);
      bool run_visible = false;
      for (int j = 1; j <= pieces; ++j) {
        math::Vector piece_end = b;
        if (j < pieces) {
          const float t = static_cast<float>(j) * inv_pieces;
          piece_end = newton_unit(azimuthal_unproject(
              pa.first + dx * t, pa.second + dy * t, planar_basis));
        }
        const math::PixelCoords piece_end_px =
            math::vector_to_pixel<W, H>(piece_end);
        const bool visible =
            band.may_reach(piece_start_px, piece_end_px, piece_arc);
        if (j > 1 && visible != run_visible)
          push(piece_start, run_visible);
        run_visible = visible;
        piece_start = piece_end;
        piece_start_px = piece_end_px;
      }
      push(b, run_visible);
    }
    return {flag_storage, out.size() - 1};
  }

private:
  uint8_t *flag_storage = nullptr;
  int capacity = 0;
};

#if HS_ENABLE_TEST_HOOKS
/** @brief Whether PlanarChords band-splits its pole runs; tests clear it to
 *  pin the chord walk's exact clip parity. */
inline bool g_planar_chords_split_pole_runs = true;
#endif

/** @brief Raster configuration PlanarChords hands pole runs to. */
inline constexpr RasterConfig PLANAR_CHORD_RASTER_CONFIG{
    .single_pass = true,
    .derive_planar_arc_registers = false,
    .interpolate_registers = false,
    .sampling_policy = RasterSamplingPolicy::SELECTABLE};

/**
 * @brief Strokes closed polylines that are straight in an azimuthal-equidistant
 *        chart as fixed-step chords between projected anchor points.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each edge is subdivided into at most MAX_ANCHOR_INTERVALS azimuthal
 * anchors and walked at TARGET_STEP screen pixels, with no per-sample tangent or
 * pole scaling. Edges and samples that cannot reach the clip band are skipped
 * without moving any drawn sample. Anchor intervals within POLE_PIECE_ROWS of a
 * pole go to Plot::rasterize under PLANAR_CHORD_RASTER_CONFIG.
 */
template <int W, int H> class PlanarChords {
public:
  /** @brief Maximum rotation lever arm over a chart-straight edge.
   * @param a,b Equidistant chart endpoints in radians, within radius pi.
   * @return Maximum distance from the chart pole's axis on the unit sphere. */
  HS_HOT_INLINE static float rotation_radius(math::PixelCoords a,
                                             math::PixelCoords b) {
    const float dx = b.x - a.x, dy = b.y - a.y;
    const float span2 = dx * dx + dy * dy;
    const float nearest =
        span2 > 0.0f ? hs::clamp(-(a.x * dx + a.y * dy) / span2, 0.0f, 1.0f)
                     : 0.0f;
    const float near_x = a.x + dx * nearest, near_y = a.y + dy * nearest;
    const float radius_lo = sqrtf(near_x * near_x + near_y * near_y);
    const float radius_hi =
        sqrtf(fmaxf(a.x * a.x + a.y * a.y, b.x * b.x + b.y * b.y));
    return sinf(hs::clamp(math::PI_F * 0.5f, radius_lo, radius_hi));
  }

  /** @brief Most chart-straight anchor intervals one edge takes. */
  static constexpr int MAX_ANCHOR_INTERVALS = 6;
  /** @brief Most intervals in a moving radial edge profile. */
  static constexpr int MAX_RADIAL_INTERVALS = 48;
  /** @brief Target chart arc between moving radial anchors, radians. */
  static constexpr float RADIAL_ANCHOR_ARC = math::PI_F / 24.0f;
  /** @brief Target chart arc between anchors, in radians.
   * @details Edges longer than MAX_ANCHOR_INTERVALS * MAX_ANCHOR_ARC use
   * edge_arc / MAX_ANCHOR_INTERVALS spacing. */
  static constexpr float MAX_ANCHOR_ARC = math::PI_F / 36.0f;
  /** @brief Rows from a pole, past the anchor gap, that send an edge to the
   *  pole split. */
  static constexpr float POLE_GUARD_ROWS = 3.0f;
  /** @brief Rows from a pole within which an anchor interval of a pole edge
   *  goes to the rasterizer; the rest of the edge keeps the chord walk. */
  static constexpr float POLE_PIECE_ROWS = 11.0f;
  /** @brief Chart pieces a pole run is split into against the clip band. */
  static constexpr int POLE_RUN_PIECES = 8;
  /** @brief Chord sample spacing in screen pixels. */
  static constexpr float TARGET_STEP = 1.2f;
  /** @brief Physical temporal departure selecting the polar chart, row pitches. */
  static constexpr float POLAR_TIME_ERROR_PX = 0.2f;
  /** @brief Empirical trim holding the chord walk's stroke brightness on the
   *  balanced adaptive walk's at TARGET_STEP. */
  static constexpr float ALPHA_GAIN = 1.028f;

  /**
   * @brief Arena bytes init_storage() takes for polylines of up to
   *        @p max_vertices vertices.
   * @param max_vertices Vertex capacity.
   * @return Byte count.
   */
  static constexpr size_t storage_bytes(int max_vertices) {
    return static_cast<size_t>(max_vertices + 1) *
               (2 * sizeof(float) + sizeof(math::PixelCoords)) +
           (MAX_ANCHOR_INTERVALS + MAX_RADIAL_INTERVALS + 2) *
               sizeof(math::PixelCoords) +
           4 * alignof(math::PixelCoords) +
           PlanarBandSplit<W, H>::storage_bytes(POLE_RUN_POINTS);
  }

  /** @brief Peak transient bytes used by rasterize_run().
   *  @return Byte count. */
  static constexpr size_t scratch_a_bytes() {
    static_assert(PLANAR_CHORD_RASTER_CONFIG.single_pass &&
                  !PLANAR_CHORD_RASTER_CONFIG.derive_planar_arc_registers);
    return (2 + POLE_RUN_POINTS) * sizeof(Fragment) + 2 * alignof(Fragment);
  }

  /** @brief Binds chart and scratch storage for up to @p max_vertices
   *  vertices.
   *  @param arena Arena the storage is allocated from.
   *  @param max_vertices Vertex capacity; at least 1. */
  HS_COLD_MEMBER void init_storage(Arena &arena, int max_vertices) {
    HS_CHECK(max_vertices >= 1, "PlanarChords: max_vertices %d < 1",
             max_vertices);
    capacity = max_vertices;
    chart_x_storage = arena.allocate_n<float>(max_vertices + 1);
    chart_y_storage = arena.allocate_n<float>(max_vertices + 1);
    pixels = arena.allocate_n<math::PixelCoords>(max_vertices + 1);
    anchors = arena.allocate_n<math::PixelCoords>(MAX_ANCHOR_INTERVALS + 1);
    radial_coefficients =
        arena.allocate_n<math::PixelCoords>(MAX_RADIAL_INTERVALS + 1);

    pole_split.init_storage(arena, POLE_RUN_POINTS);
  }

  /** @brief Chart x coordinates of the next polyline, one per vertex.
   *  @return Writable array. */
  float *chart_x() { return chart_x_storage; }
  /** @brief Chart y coordinates of the next polyline, one per vertex.
   *  @return Writable array. */
  float *chart_y() { return chart_y_storage; }

  /** @brief Captures the frame's clip band; call once per frame.
   *  @param clip The frame's clip region. */
  void prepare(const ClipRegion &clip) {
    band = ClipBand<W, H>::of(clip);
    clip_stamp = clip;
    prepared = true;
  }

  /**
   * @brief Strokes one closed chart-straight polyline.
   * @tparam PipelineT Screen-space plot sink (plot(canvas, x, y, ...)).
   * @tparam F Fragment-shader callable type, used by the pole runs.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param points Vertex positions, @p vertices plus the closing repeat of the
   * first.
   * @param vertices Vertex count; chart_x()/chart_y() hold their chart
   * coordinates in @p planar_basis.
   * @param planar_basis Azimuthal-equidistant chart the edges are straight in.
   * @param color Stroke color; its alpha is balanced against the step taken.
   * @param fragment_shader Per-fragment shader for the pole runs.
   */
  template <typename PipelineT, typename F>
  HS_O3_FN HS_NOINLINE_NOCLONE void
  draw_closed(PipelineT &pipeline, Canvas &canvas, const Fragments &points,
              int vertices, const math::Basis &planar_basis,
              const Color4 &color, const F &fragment_shader) {
    using Geometry = math::DisplayGeometry<H>;
    require_closed(pipeline, canvas, points, vertices);
    float *const x_chart = chart_x_storage;
    float *const y_chart = chart_y_storage;
    for (int i = 0; i < vertices; ++i)
      pixels[i] = math::vector_to_pixel<W, H>(points[i].pos);
    x_chart[vertices] = x_chart[0];
    y_chart[vertices] = y_chart[0];
    pixels[vertices] = pixels[0];

    auto walk_segment =
        [&](const math::PixelCoords &from,
            const math::PixelCoords &to) __attribute__((always_inline)) {
          const float x_start = from.x;
          const float y_start = from.y;
          const float segment_dx = to.x - x_start;
          const float segment_dy = to.y - y_start;
          const float length =
              sqrtf(segment_dx * segment_dx + segment_dy * segment_dy);
          const int samples =
              std::max(1, static_cast<int>(ceilf(length / TARGET_STEP)));
          const float sample_count = static_cast<float>(samples);
          const float inv_samples = 1.0f / sample_count;
          // Sample index range whose splat can reach the band; plot() masks the
          // taps, so a range wider than exact only costs time.
          float first = 0.0f;
          float last = sample_count;
          auto restrict_range = [&](float start, float delta, float lo,
                                    float hi) __attribute__((always_inline)) {
            if (fabsf(delta) < 1e-4f) {
              if (start < lo - 1.0f || start > hi + 1.0f)
                last = -1.0f;
              return;
            }
            const float scale = sample_count / delta;
            const float s_lo = (lo - start) * scale;
            const float s_hi = (hi - start) * scale;
            first = fmaxf(first, fminf(s_lo, s_hi) - 1.0f);
            last = fminf(last, fmaxf(s_lo, s_hi) + 1.0f);
          };
          restrict_range(y_start, segment_dy, band.row_lo, band.row_hi);
          if (band.x_active && fabsf(segment_dx) < W * 0.25f &&
              fabsf(segment_dx) < W - 2.0f * band.x_half) {
            const float shift =
                W * rintf((x_start + 0.5f * segment_dx - band.x_center) / W);
            restrict_range(x_start, segment_dx,
                           band.x_center - band.x_half + shift,
                           band.x_center + band.x_half + shift);
          }
          if (last < first)
            return;
          const int sample_begin = static_cast<int>(ceilf(fmaxf(first, 0.0f)));
          const int sample_end =
              static_cast<int>(fminf(floorf(last) + 1.0f, sample_count));
          const float step_ratio =
              fminf(TARGET_STEP, length * inv_samples) / SCREEN_STEP_PX;
          const float sample_alpha =
              fminf(1.0f, ALPHA_GAIN *
                              balanced_sample_alpha(color.alpha, step_ratio));
          for (int sample = sample_begin; sample < sample_end; ++sample) {
            const float t = static_cast<float>(sample) * inv_samples;
            const float x = math::fast_wrap(x_start + segment_dx * t, W);
            const float y = y_start + segment_dy * t;
            pipeline.plot(canvas, x, y, color.color, 0.0f, sample_alpha);
          }
        };

    for (int edge = 0; edge < vertices; ++edge) {
      const float dx = x_chart[edge + 1] - x_chart[edge];
      const float dy = y_chart[edge + 1] - y_chart[edge];
      const float edge_arc = sqrtf(dx * dx + dy * dy);
      const math::PixelCoords &pa = pixels[edge];
      const math::PixelCoords &pb = pixels[edge + 1];
      if (!band.may_reach(pa, pb, edge_arc))
        continue;

      const int anchor_intervals =
          hs::clamp(static_cast<int>(ceilf(edge_arc / MAX_ANCHOR_ARC)), 1,
                    MAX_ANCHOR_INTERVALS);
      const float inv_intervals = 1.0f / static_cast<float>(anchor_intervals);
      const float gap_arc = edge_arc * inv_intervals;
      auto chart_position = [&](float t) __attribute__((always_inline)) {
        return newton_unit(azimuthal_unproject(
            x_chart[edge] + dx * t, y_chart[edge] + dy * t, planar_basis));
      };
      float anchor_row_lo = pa.y;
      float anchor_row_hi = pa.y;
      anchors[0] = pa;
      for (int k = 1; k <= anchor_intervals; ++k) {
        anchors[k] = k == anchor_intervals
                         ? pb
                         : math::vector_to_pixel<W, H>(chart_position(
                               static_cast<float>(k) * inv_intervals));
        const float delta = anchors[k].x - anchors[k - 1].x;
        if (delta > W * 0.5f)
          anchors[k].x -= W;
        else if (delta < -W * 0.5f)
          anchors[k].x += W;
        anchor_row_lo = fminf(anchor_row_lo, anchors[k].y);
        anchor_row_hi = fmaxf(anchor_row_hi, anchors[k].y);
      }

      const float row_margin = gap_arc * math::ROWS_PER_RADIAN<H> + 1.0f;
      if (anchor_row_lo - row_margin >=
              Geometry::NORTH_POLE_ROW + POLE_GUARD_ROWS &&
          anchor_row_hi + row_margin <=
              Geometry::SOUTH_POLE_ROW - POLE_GUARD_ROWS) {
        HS_PROFILE_DEEP(plot_chord_walk);
        for (int k = 0; k < anchor_intervals; ++k)
          walk_segment(anchors[k], anchors[k + 1]);
        continue;
      }

      HS_PROFILE_DEEP(plot_chord_pole);
      auto near_pole = [&](int k) __attribute__((always_inline)) {
        return fminf(anchors[k].y, anchors[k + 1].y) - row_margin <
                   Geometry::NORTH_POLE_ROW + POLE_PIECE_ROWS ||
               fmaxf(anchors[k].y, anchors[k + 1].y) + row_margin >
                   Geometry::SOUTH_POLE_ROW - POLE_PIECE_ROWS;
      };
      auto run_position = [&](float t) __attribute__((always_inline)) {
        if (t <= 0.0f)
          return points[edge].pos;
        if (t >= 1.0f)
          return points[edge + 1].pos;
        return chart_position(t);
      };
      for (int k = 0; k < anchor_intervals;) {
        if (!near_pole(k)) {
          walk_segment(anchors[k], anchors[k + 1]);
          ++k;
          continue;
        }
        int end = k + 1;
        while (end < anchor_intervals && near_pole(end))
          ++end;
        const float t0 = static_cast<float>(k) * inv_intervals;
        const float t1 = end == anchor_intervals
                             ? 1.0f
                             : static_cast<float>(end) * inv_intervals;
        k = end;
        rasterize_run(pipeline, canvas, run_position(t0), run_position(t1),
                      planar_basis, fragment_shader);
      }
    }
  }

  /** @brief Captures the frame's clip band and the gate of polylines moving
   *  over the frame; call once per frame.
   *  @param canvas The frame's canvas.
   *  @param motion Bound on any polyline point's motion over the frame,
   *  radians. */
  void prepare(const Canvas &canvas, float motion) {
    prepare(canvas.clip());
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    reach = motion;
  }

  /**
   * @brief Strokes one closed chart-straight polyline turning rigidly over the
   *        frame, each sample at its column's time.
   * @tparam PipelineT Screen-space plot sink (plot(canvas, x, y, ...)).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param points Frame-end vertex positions, @p vertices plus the closing
   * repeat of the first.
   * @param vertices Vertex count; chart_x()/chart_y() hold their chart
   * coordinates in @p planar_basis.
   * @param planar_basis Azimuthal-equidistant chart the edges are straight in.
   * @param color Stroke color; its alpha is balanced against the step taken.
   * @param start Rotation taking a frame-end point to its frame-start
   * position; no point moves farther than prepare()'s motion.
   * @param max_intervals Most anchor intervals an edge takes; an edge's
   * anchors lie at most MAX_ANCHOR_ARC apart up to that count.
   * @param alternating_radii Adjacent edges reverse the radial anchor profile.
   * @param middle Rotation from frame end to mid-frame.
   * @param reference_chords Reference balanced screen chords, with adaptive
   * balanced pole runs; otherwise reference adaptive arc-table sampling.
   * @param contour_speed Bound on angular speed over the contour, radians per frame.
   * @details Edges share the first edge's radial profile, reversed on alternating
   * edges. Screen anchors carry quadratic motion between frame start, middle
   * and end. Coherent samples solve their column crossing directly; other
   * chords emit each supported column at its own time. Pole intervals use
   * quadratic endpoint captures in a display-pole chart.
   */
  template <typename PipelineT>
  HS_HOT_FLASH_MEMBER void
  draw_closed(PipelineT &pipeline, Canvas &canvas, const Fragments &points,
              int vertices, const math::Basis &planar_basis,
              const Color4 &color, const math::RotationMatrix &start,
              int max_intervals, bool alternating_radii,
              const math::RotationMatrix &middle, bool reference_chords,
              float contour_speed) {
    const float target_step =
        reference_chords ? TARGET_STEP : BALANCED_SCREEN_STEP_PX;
    using Geometry = math::DisplayGeometry<H>;
    require_closed(pipeline, canvas, points, vertices);
    ScratchScope length_scope(scratch_arena_a);
    auto &arc_cumul =
        *scratch_arena_a.make<std::array<float, PLANAR_LEN_SAMPLES + 1>>();

    float *const x_chart = chart_x_storage;
    float *const y_chart = chart_y_storage;
    for (int i = 0; i < vertices; ++i) {
      pixels[i] = math::vector_to_pixel<W, H>(points[i].pos);
    }
    x_chart[vertices] = x_chart[0];
    y_chart[vertices] = y_chart[0];
    pixels[vertices] = pixels[0];

    const float first_dx = x_chart[1] - x_chart[0];
    const float first_dy = y_chart[1] - y_chart[0];
    const float edge_arc = sqrtf(first_dx * first_dx + first_dy * first_dy);
    const int anchor_intervals =
        hs::clamp(static_cast<int>(ceilf(edge_arc / (alternating_radii
                                                         ? MAX_ANCHOR_ARC
                                                         : RADIAL_ANCHOR_ARC))),
                  1, std::min(max_intervals, MAX_RADIAL_INTERVALS));
    const float inv_intervals = 1.0f / static_cast<float>(anchor_intervals);
    planar_arc_cumul({x_chart[0], y_chart[0]}, first_dx, first_dy, planar_basis,
                     arc_cumul);
    auto *reference_density =
        scratch_arena_a.allocate_n<float>(anchor_intervals);
    for (int k = 0; k < anchor_intervals; ++k) {
      const float t = (float(k) + 0.5f) * inv_intervals;
      const float x = x_chart[0] + first_dx * t;
      const float y = y_chart[0] + first_dy * t;
      const float radius = sqrtf(x * x + y * y);
      const float dr = radius > math::EPS_GEOMETRIC
                           ? (x * first_dx + y * first_dy) / radius
                           : 0.0f;
      float sine, cosine;
      math::fast_sincosf_0_pi(radius, sine, cosine);
      const float sinc = radius > math::EPS_GEOMETRIC ? sine / radius : 1.0f;
      const float speed = sqrtf(
          fmaxf(0.0f, dr * dr + sinc * sinc * (edge_arc * edge_arc - dr * dr)));
      const int j =
          std::min(PLANAR_LEN_SAMPLES - 1, int(t * PLANAR_LEN_SAMPLES));
      const float slope =
          PLANAR_LEN_SAMPLES * (arc_cumul[j + 1] - arc_cumul[j]);
      reference_density[k] =
          reference_chords || !(speed > 0.0f) ? 1.0f : slope / speed;
    }
    const float radial_bend =
        alternating_radii ? edge_arc * edge_arc * inv_intervals *
                                inv_intervals * math::ROWS_PER_RADIAN<H> / 8.0f
                          : 0.0f;
    for (int k = 1; k < anchor_intervals; ++k) {
      const float t = static_cast<float>(k) * inv_intervals;
      const float x = x_chart[0] + first_dx * t;
      const float y = y_chart[0] + first_dy * t;
      const float radius = sqrtf(x * x + y * y);
      float sine, cosine;
      math::fast_sincosf_0_pi(radius, sine, cosine);
      radial_coefficients[k] = {
          radius > math::EPS_GEOMETRIC ? sine / radius : 1.0f, cosine};
    }

    const float velocity = contour_speed * math::ROWS_PER_RADIAN<H>;
    auto walk_segment = [&](const math::PixelCoords &from0,
                            const math::PixelCoords &to0,
                            const math::PixelCoords &from1,
                            const math::PixelCoords &to1,
                            const math::PixelCoords &from_mid,
                            const math::PixelCoords &to_mid, float density,
                            float curve_bend) __attribute__((always_inline)) {
      auto sample_alpha = [&](float ratio) HS_HOT_INLINE {
        const float alpha = !reference_chords && ratio < 1.0f
                                ? color.alpha * ratio
                                : balanced_sample_alpha(color.alpha, ratio);
        return fminf(1.0f, density * alpha);
      };
      const math::PixelCoords from_curve{
          2.0f * (from0.x + from1.x - 2.0f * from_mid.x),
          2.0f * (from0.y + from1.y - 2.0f * from_mid.y)};
      const math::PixelCoords to_curve{2.0f * (to0.x + to1.x - 2.0f * to_mid.x),
                                       2.0f *
                                           (to0.y + to1.y - 2.0f * to_mid.y)};
      const float x_pad = 0.25f * fmaxf(fabsf(from_curve.x), fabsf(to_curve.x));
      const float y_pad = 0.25f * fmaxf(fabsf(from_curve.y), fabsf(to_curve.y));
      const float dx1 = to1.x - from1.x, dy1 = to1.y - from1.y;
      const float dx0 = to0.x - from0.x, dy0 = to0.y - from0.y;
      const float length = sqrtf(dx1 * dx1 + dy1 * dy1);
      const int samples =
          std::max(1, static_cast<int>(ceilf(length / target_step)));
      const float sample_count = static_cast<float>(samples);
      const float inv_samples = 1.0f / sample_count;
      float first = sample_count, last = 0.0f;
      auto widen = [&](float x0, float y0, float sdx,
                       float sdy) __attribute__((always_inline)) {
        float lo = 0.0f, hi = sample_count;
        auto restrict_range = [&](float start_at, float delta, float low,
                                  float high) __attribute__((always_inline)) {
          if (fabsf(delta) < 1e-4f) {
            if (start_at < low - 1.0f || start_at > high + 1.0f)
              hi = -1.0f;
            return;
          }
          const float scale = sample_count / delta;
          const float s_lo = (low - start_at) * scale;
          const float s_hi = (high - start_at) * scale;
          lo = fmaxf(lo, fminf(s_lo, s_hi) - 1.0f);
          hi = fminf(hi, fmaxf(s_lo, s_hi) + 1.0f);
        };
        restrict_range(y0, sdy, band.row_lo - y_pad, band.row_hi + y_pad);
        if (band.x_active && fabsf(sdx) < W * 0.25f &&
            fabsf(sdx) < W - 2.0f * band.x_half) {
          const float shift = W * rintf((x0 + 0.5f * sdx - band.x_center) / W);
          restrict_range(x0, sdx, band.x_center - band.x_half + shift - x_pad,
                         band.x_center + band.x_half + shift + x_pad);
        }
        if (hi >= lo) {
          first = fminf(first, lo);
          last = fmaxf(last, hi);
        }
      };
      widen(from1.x, from1.y, dx1, dy1);
      widen(from0.x, from0.y, dx0, dy0);
      if (last < first)
        return;
      const int sample_begin = static_cast<int>(ceilf(fmaxf(first, 0.0f)));
      const int sample_end =
          static_cast<int>(fminf(floorf(last) + 1.0f, sample_count));
      const float step_ratio =
          fminf(target_step, length * inv_samples) / SCREEN_STEP_PX;
      const Color4 sample_color(color.color, sample_alpha(step_ratio));
      constexpr float HALF = W / 2.0f;
      const int arm = static_cast<int>(floorf(from0.x / HALF));
      const float arm_start = arm * HALF;
      const float low = fminf(fminf(from0.x, to0.x), fminf(from1.x, to1.x));
      const float high = fmaxf(fmaxf(from0.x, to0.x), fmaxf(from1.x, to1.x));
      const float mx0 = from1.x - from0.x;
      const float my0 = from1.y - from0.y;
      const float dmx = dx1 - dx0, dmy = dy1 - dy0;
      const float variation2 = dmx * dmx + dmy * dmy;
      const float nearest =
          variation2 > 0.0f
              ? hs::clamp(-(dx0 * dmx + dy0 * dmy) / variation2, 0.0f, 1.0f)
              : 0.0f;
      const float near_x = dx0 + dmx * nearest, near_y = dy0 + dmy * nearest;
      const float dx_mid = to_mid.x - from_mid.x,
                  dy_mid = to_mid.y - from_mid.y;
      const float curve_x = 0.5f * (dx0 + dx1 - 2.0f * dx_mid);
      const float curve_y = 0.5f * (dy0 + dy1 - 2.0f * dy_mid);
      const float span_bend = sqrtf(curve_x * curve_x + curve_y * curve_y);
      const float start_length = sqrtf(dx0 * dx0 + dy0 * dy0);
      const float lower =
          fmaxf(0.0f, sqrtf(near_x * near_x + near_y * near_y) - span_bend);
      const float upper = fmaxf(start_length, length) + span_bend;
      const bool steady_count =
          std::max(1, int(ceilf(lower / target_step))) == samples &&
          std::max(1, int(ceilf(upper / target_step))) == samples;
      const float curve_dx = to_curve.x - from_curve.x;
      const float curve_dy = to_curve.y - from_curve.y;
      const float denominator = HALF - fmaxf(mx0 + fabsf(from_curve.x),
                                             mx0 + dmx + fabsf(to_curve.x));
      const float first_derivative =
          denominator > 0.0f
              ? (fmaxf(fabsf(dx0), fabsf(dx1)) + 0.25f * fabsf(curve_dx)) /
                    denominator
              : std::numeric_limits<float>::infinity();
      const float time_bend =
          denominator > 0.0f
              ? (2.0f * (fabsf(dmx) + fabsf(curve_dx)) * first_derivative +
                 2.0f * fmaxf(fabsf(from_curve.x), fabsf(to_curve.x)) *
                     first_derivative * first_derivative) /
                    (8.0f * denominator)
              : std::numeric_limits<float>::infinity();
      if (low - x_pad > arm_start + 1.0f &&
          high + x_pad < arm_start + HALF - 1.0f &&
          (velocity + 4.0f * curve_bend) * time_bend <= 0.002f &&
          steady_count) {
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::chord_capture_solves>(2);
        const float t0 = Spacetime::quadratic_screen_time<W>(
            from0.x, mx0, from_curve.x, arm_start);
        const float dt = Spacetime::quadratic_screen_time<W>(
                             to0.x, mx0 + dmx, to_curve.x, arm_start) -
                         t0;
        const float k0 = t0 * (t0 - 1.0f);
        const float k1 = dt * (2.0f * t0 - 1.0f);
        const float k2 = dt * dt;
        const float x0 = from0.x + mx0 * t0 + from_curve.x * k0;
        const float x1 =
            dx0 + mx0 * dt + dmx * t0 + from_curve.x * k1 + curve_dx * k0;
        const float x2 = dmx * dt + from_curve.x * k2 + curve_dx * k1;
        const float x3 = curve_dx * k2;
        const float y0 = from0.y + my0 * t0 + from_curve.y * k0;
        const float y1 =
            dy0 + my0 * dt + dmy * t0 + from_curve.y * k1 + curve_dy * k0;
        const float y2 = dmy * dt + from_curve.y * k2 + curve_dy * k1;
        const float y3 = curve_dy * k2;
        const float alpha_start = sample_alpha(
            fminf(target_step, start_length * inv_samples) / SCREEN_STEP_PX);
        const float alpha_middle = sample_alpha(
            fminf(target_step,
                  sqrtf(dx_mid * dx_mid + dy_mid * dy_mid) * inv_samples) /
            SCREEN_STEP_PX);
        const float alpha_curve =
            2.0f * (alpha_start + sample_color.alpha - 2.0f * alpha_middle);
        for (int sample = sample_begin; sample < sample_end; ++sample) {
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::chord_fast_samples>();
          const float u = static_cast<float>(sample) * inv_samples;
          const float x = x0 + u * (x1 + u * (x2 + u * x3));
          const float time = t0 + dt * u;
          const float alpha = alpha_start +
                              (sample_color.alpha - alpha_start) * time +
                              alpha_curve * time * (time - 1.0f);
          const float y = y0 + u * (y1 + u * (y2 + u * y3));
          const float vx = mx0 + dmx * u +
                           (from_curve.x + curve_dx * u) * (2.0f * time - 1.0f);
          const float vy = my0 + dmy * u +
                           (from_curve.y + curve_dy * u) * (2.0f * time - 1.0f);
          PointRoute<W, H>::emit(pipeline, canvas, {{x, y}, time},
                                 {sample_color.color, alpha}, vx, vy);
        }
        return;
      }
      const float guard = 0.25f * fmaxf(fabsf(from_curve.x), fabsf(to_curve.x));
      const int owner_begin = static_cast<int>(floorf(low - guard - 1.0f));
      const int owner_end = static_cast<int>(ceilf(high + guard + 1.0f));
      for (int owner = owner_begin; owner <= owner_end; ++owner) {
        const int column = math::fast_wrap(owner, W);
        if (canvas.render_x_clip().clipped(column))
          continue;
        const float t = Spacetime::SamplingSurface<W>::time(column);
        const float bend = t * (t - 1.0f);
        const math::PixelCoords at_from{
            from0.x + (from1.x - from0.x) * t + from_curve.x * bend,
            from0.y + (from1.y - from0.y) * t + from_curve.y * bend};
        const math::PixelCoords at_to{
            to0.x + (to1.x - to0.x) * t + to_curve.x * bend,
            to0.y + (to1.y - to0.y) * t + to_curve.y * bend};
        const float dx = at_to.x - at_from.x, dy = at_to.y - at_from.y;
        const float span = sqrtf(dx * dx + dy * dy);
        const int count =
            std::max(1, static_cast<int>(ceilf(span / target_step)));
        const float inverse = 1.0f / static_cast<float>(count);
        float lo = 0.0f, hi = static_cast<float>(count - 1);
        auto restrict = [&](float start_at, float delta, float lower,
                            float upper) HS_HOT_INLINE {
          if (fabsf(delta) < 1e-6f) {
            if (start_at < lower || start_at > upper)
              hi = -1.0f;
          } else {
            const float scale = static_cast<float>(count) / delta;
            const float a = (lower - start_at) * scale,
                        b = (upper - start_at) * scale;
            lo = fmaxf(lo, fminf(a, b));
            hi = fminf(hi, fmaxf(a, b));
          }
        };
        restrict(at_from.x, dx, float(owner) - 1.0f, float(owner) + 1.0f);
        restrict(at_from.y, dy, band.row_lo, band.row_hi);
        const int begin = static_cast<int>(ceilf(lo));
        const int end = static_cast<int>(floorf(hi));
        const float alpha =
            sample_alpha(fminf(target_step, span * inverse) / SCREEN_STEP_PX);
        for (int sample = begin; sample <= end; ++sample) {
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::chord_screen_samples>();
          const float u = static_cast<float>(sample) * inverse;
          float x = at_from.x + dx * u;
          x += W * rintf((float(column) - x) / W);
          pipeline.plot_column_at(canvas, {x, at_from.y + dy * u}, column,
                                  color.color, 0.0f, alpha, t);
        }
      }
    };

    for (int edge = 0; edge < vertices; ++edge) {
      const float dx = x_chart[edge + 1] - x_chart[edge];
      const float dy = y_chart[edge + 1] - y_chart[edge];
      if (!band.may_reach(pixels[edge], pixels[edge + 1],
                          edge_arc + 2.0f * reach))
        continue;
      const float row_margin =
          edge_arc * inv_intervals * math::ROWS_PER_RADIAN<H> + 1.0f;
      auto chart_position = [&](int k) __attribute__((always_inline)) {
        if (k == 0)
          return points[edge].pos;
        if (k == anchor_intervals)
          return points[edge + 1].pos;
        const float t = static_cast<float>(k) * inv_intervals;
        const int profile =
            alternating_radii && (edge & 1) ? anchor_intervals - k : k;
        const auto &coefficient = radial_coefficients[profile];
        const float x = x_chart[edge] + dx * t;
        const float y = y_chart[edge] + dy * t;
        return newton_unit(planar_basis.v * coefficient.y +
                           (planar_basis.u * x + planar_basis.w * y) *
                               coefficient.x);
      };
      math::Vector from = points[edge].pos;
      math::PixelCoords end_from = pixels[edge];
      math::PixelCoords start_from{}, mid_from{};
      bool captured_from = false;
      math::Vector start_world_from{}, mid_world_from{};
      PolarAnchor polar_from{};
      bool polar_from_ready = false, polar_from_north = true;
      for (int k = 0; k < anchor_intervals; ++k) {
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::chord_intervals>();
        const bool last_interval = k + 1 == anchor_intervals;
        const math::Vector to = chart_position(k + 1);
        math::PixelCoords end_to =
            last_interval ? pixels[edge + 1] : math::vector_to_pixel<W, H>(to);
        end_to.x += W * rintf((end_from.x - end_to.x) * (1.0f / W));
        if (!band.may_reach(end_from, end_to,
                            edge_arc * inv_intervals + 2.0f * reach)) {
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::chord_intervals_culled>();
          from = to;
          end_from = end_to;
          captured_from = false;
          polar_from_ready = false;
          continue;
        }
        if (!captured_from) {
          start_world_from = start.apply(from);
          mid_world_from = middle.apply(from);
          start_from = math::vector_to_pixel<W, H>(start_world_from);
          mid_from = math::vector_to_pixel<W, H>(mid_world_from);
          start_from.x += W * rintf((end_from.x - start_from.x) / W);
          mid_from.x += W * rintf((end_from.x - mid_from.x) / W);
        }
        const math::Vector start_world_to = start.apply(to);
        const math::Vector mid_world_to = middle.apply(to);
        math::PixelCoords start_to =
            math::vector_to_pixel<W, H>(start_world_to);
        math::PixelCoords mid_to = math::vector_to_pixel<W, H>(mid_world_to);
        start_to.x += W * rintf((end_to.x - start_to.x) / W);
        mid_to.x += W * rintf((end_to.x - mid_to.x) / W);
        const float row_lo =
            fminf(fminf(end_from.y, end_to.y), fminf(start_from.y, start_to.y));
        const float row_hi =
            fmaxf(fmaxf(end_from.y, end_to.y), fmaxf(start_from.y, start_to.y));
        const float pole_distance = fminf(row_lo - Geometry::NORTH_POLE_ROW,
                                          Geometry::SOUTH_POLE_ROW - row_hi);
        const int row = hs::clamp(
            static_cast<int>(floorf(Geometry::NORTH_POLE_ROW +
                                    fmaxf(0.0f, pole_distance - row_margin))),
            0, H - 1);
        const float sine = fmaxf(0.001f, math::TrigLUT<W, H>::sin_phi[row]);
        const float bend = contour_speed * contour_speed *
                           math::ROWS_PER_RADIAN<H> / (8.0f * sine);
        const bool polar = pole_distance < POLE_PIECE_ROWS + row_margin ||
                           bend > POLAR_TIME_ERROR_PX || radial_bend > 0.05f;
        if (!polar) {
          polar_from_ready = false;
          walk_segment(start_from, start_to, end_from, end_to, mid_from, mid_to,
                       reference_density[alternating_radii && (edge & 1)
                                             ? anchor_intervals - 1 - k
                                             : k],
                       bend);
        } else {
          const bool north = row_lo + row_hi < Geometry::NORTH_POLE_ROW +
                                                   Geometry::SOUTH_POLE_ROW;
          auto capture = [&](const math::Vector &a, const math::Vector &m,
                             const math::Vector &b) HS_HOT_INLINE {
            return PolarAnchor::capture(PolarWalk<W, H>::project(a, north),
                                        PolarWalk<W, H>::project(m, north),
                                        PolarWalk<W, H>::project(b, north));
          };
          if (!polar_from_ready || polar_from_north != north)
            polar_from = capture(start_world_from, mid_world_from, from);
          const PolarAnchor polar_to =
              capture(start_world_to, mid_world_to, to);
          PolarWalk<W, H>::draw(
              pipeline, canvas, polar_from, polar_to, north, color, target_step,
              reference_chords,
              reference_density[alternating_radii && (edge & 1)
                                    ? anchor_intervals - 1 - k
                                    : k]);
          polar_from = polar_to;
          polar_from_north = north;
          polar_from_ready = true;
        }
        start_world_from = start_world_to;
        mid_world_from = mid_world_to;
        from = to;
        end_from = end_to;
        start_from = start_to;
        mid_from = mid_to;
        captured_from = true;
      }
    }
  }

private:
  /** @brief Checks that a closed polyline includes its repeated first vertex.
   * @param points Polyline storage.
   * @param vertices Distinct vertex count. */
  template <typename PipelineT>
  HS_HOT_INLINE void require_closed(PipelineT &pipeline, Canvas &canvas,
                                    const Fragments &points,
                                    int vertices) const {
    if constexpr (requires { pipeline.prepared_for(canvas); })
      HS_CHECK(pipeline.prepared_for(canvas),
               "direct raster pipeline not prepared for this canvas");
    HS_CHECK(prepared && clip_stamp == canvas.clip(),
             "PlanarChords: prepare() not called for this canvas clip");
    HS_CHECK(vertices >= 1 && vertices <= capacity,
             "PlanarChords: %d vertices outside capacity %d", vertices,
             capacity);
    HS_CHECK(points.size() > static_cast<size_t>(vertices),
             "PlanarChords: closing vertex is missing");
  }

  /** @brief Plot-rasterizes one chart-straight run with balanced sampling. */
  template <typename PipelineT, typename F>
  HS_FLASH_MEMBER void
  rasterize_run(PipelineT &pipeline, Canvas &canvas, const math::Vector &a,
                const math::Vector &b, const math::Basis &planar_basis,
                const F &fragment_shader) {
    ScratchScope guard(scratch_arena_a);
    Fragments run;
    run.bind(scratch_arena_a, 2);
    Fragment point;
    point.pos = a;
    run.push_back(point);
    point.pos = b;
    run.push_back(point);
    bool split = band.x_active;
#if HS_ENABLE_TEST_HOOKS
    split = split && g_planar_chords_split_pole_runs;
#endif
    if (!split) {
      rasterize<W, H, PLANAR_CHORD_RASTER_CONFIG>(
          pipeline, canvas, run, fragment_shader,
          {.projection = RasterProjection::planar(planar_basis),
           .omit_end = true,
           .balanced_sampling = true});
      return;
    }
    Fragments pieces;
    pieces.bind(scratch_arena_a, static_cast<size_t>(POLE_RUN_POINTS));
    const std::span<const uint8_t> flags =
        pole_split.split(pieces, run, 1, POLE_RUN_PIECES, planar_basis, band);
    rasterize<W, H, PLANAR_CHORD_RASTER_CONFIG>(
        pipeline, canvas, pieces, fragment_shader,
        {.projection = RasterProjection::planar(planar_basis, flags),
         .omit_end = true,
         .balanced_sampling = true});
  }

  static constexpr int POLE_RUN_POINTS =
      PlanarBandSplit<W, H>::max_points(1, POLE_RUN_PIECES);

  ClipBand<W, H> band;
  ClipRegion clip_stamp;
  bool prepared = false;
  PlanarBandSplit<W, H> pole_split;
  float *chart_x_storage = nullptr;
  float *chart_y_storage = nullptr;
  math::PixelCoords *pixels = nullptr;
  math::PixelCoords *anchors = nullptr;
  math::PixelCoords *radial_coefficients =
      nullptr;        ///< Shared radial anchor factors.
  float reach = 0.0f; ///< Bound on a point's motion over the frame.
  int capacity = 0;
};

} // namespace Plot
