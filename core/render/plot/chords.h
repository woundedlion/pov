/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file chords.h
 * @brief Plot::PlanarChords: strokes closed polylines that are straight in an
 *        azimuthal-equidistant chart as fixed-step screen chords, culled to the
 *        clip band, and Plot::ClipBand, the conservative chart-segment cull it
 *        uses.
 */

#include "engine/memory.h"
#include "math/display_geometry.h"
#include "math/pixel_mapping.h"
#include "render/canvas.h"
#include "render/clip.h"
#include "render/plot/cull.h"
#include "render/plot/raster.h"

namespace Plot {

/**
 * @brief Screen rows and unwrapped columns a splat must touch to reach a clip
 *        region, and the conservative reach test for chart-straight pieces.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct ClipBand {
  float row_lo = 0.0f;
  float row_hi = 0.0f;
  float x_center = 0.0f;
  float x_half = 0.0f;
  bool x_active = false;

  /** @brief Rows and columns a splat and the fast-trig pixel map reach past a
   *  bound. */
  static constexpr float PAD = 2.0f;

  /** @brief The band a sample must splat into to touch @p clip's render
   *  region. */
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
 * @brief Cuts a closed polyline's chart-straight edges into equal chart pieces
 *        and flags the runs of pieces that cannot reach a clip band, so
 *        Plot::rasterize never simulates them.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Points are added only where a run begins or ends, so a visible run
 * is walked in one piece. Its start depends on the clip, which moves the walk's
 * sample phase inside it: clipped strokes match the unclipped one to a fraction
 * of a pixel, not bit for bit. Edges touching the chart antipode keep the
 * rasterizer's geodesic fallback whole. Flag storage lives in arena storage
 * bound once by init_storage().
 */
template <int W, int H> class PlanarBandSplit {
public:
  /** @brief Arena bytes init_storage() takes for up to @p max_points points. */
  static constexpr size_t storage_bytes(int max_points) {
    return static_cast<size_t>(max_points) + alignof(uint8_t);
  }

  /** @brief Most points split() emits for @p edges edges of @p pieces each. */
  static constexpr int max_points(int edges, int pieces) {
    return edges * pieces + 1;
  }

  /** @brief Binds flag storage for polylines of up to @p max_points points. */
  HS_COLD_MEMBER void init_storage(Arena &arena, int max_points) {
    HS_CHECK(max_points >= 2, "PlanarBandSplit: max_points %d < 2", max_points);
    capacity = max_points;
    flag_storage = arena.allocate_n<uint8_t>(max_points);
  }

  /**
   * @brief Splits @p ring into @p out and flags each output edge.
   * @param out Receives the split polyline; bound for max_points(edges,
   * pieces) points.
   * @param ring Closed polyline: @p edges vertices plus the closing repeat.
   * @param edges Edge count.
   * @param pieces Chart pieces per edge.
   * @param planar_basis Azimuthal-equidistant chart the edges are straight in.
   * @param band Clip band the pieces are tested against.
   * @return Rasterize edge flags for @p out, one per edge.
   */
  HS_HOT_FLASH_MEMBER std::span<const uint8_t>
  split(Fragments &out, const Fragments &ring, int edges, int pieces,
        const math::Basis &planar_basis, const ClipBand<W, H> &band) {
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
 * @details For dense stacks of short planar edges, where the per-edge setup of
 * Plot::rasterize's adaptive walk outweighs the edges themselves. Each edge is
 * subdivided into at most MAX_ANCHOR_INTERVALS azimuthal anchors and walked at
 * TARGET_STEP screen pixels, with no per-sample tangent or pole scaling. Edges
 * that cannot reach the clip band are skipped whole, and each chord walks only
 * the samples whose splat can reach it; both culls leave every drawn sample
 * where the unclipped walk puts it. Anchor intervals within POLE_PIECE_ROWS of
 * a pole, where the rasterizer's pole scaling is what closes the gaps, are
 * handed to Plot::rasterize under PLANAR_CHORD_RASTER_CONFIG with balanced
 * sampling, band-split by PlanarBandSplit into POLE_RUN_PIECES pieces, so a
 * clipped pole run matches the unclipped one to a fraction of a pixel. Chart coordinates and scratch live in arena storage bound once by
 * init_storage().
 */
template <int W, int H> class PlanarChords {
public:
  /** @brief Most chart-straight anchor intervals one edge takes. */
  static constexpr int MAX_ANCHOR_INTERVALS = 6;
  /** @brief Longest chart arc between consecutive anchors, in radians. */
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
  /** @brief Empirical trim holding the chord walk's stroke brightness on the
   *  balanced adaptive walk's at TARGET_STEP. */
  static constexpr float ALPHA_GAIN = 1.028f;

  /**
   * @brief Arena bytes init_storage() takes for polylines of up to
   *        @p max_vertices vertices.
   */
  static constexpr size_t storage_bytes(int max_vertices) {
    return static_cast<size_t>(max_vertices + 1) *
               (2 * sizeof(float) + sizeof(math::PixelCoords)) +
           (MAX_ANCHOR_INTERVALS + 1) * sizeof(math::PixelCoords) +
           3 * alignof(math::PixelCoords) +
           PlanarBandSplit<W, H>::storage_bytes(POLE_RUN_POINTS);
  }

  /** @brief Binds chart and scratch storage for up to @p max_vertices
   *  vertices. */
  HS_COLD_MEMBER void init_storage(Arena &arena, int max_vertices) {
    HS_CHECK(max_vertices >= 1, "PlanarChords: max_vertices %d < 1",
             max_vertices);
    capacity = max_vertices;
    chart_x_storage = arena.allocate_n<float>(max_vertices + 1);
    chart_y_storage = arena.allocate_n<float>(max_vertices + 1);
    pixels = arena.allocate_n<math::PixelCoords>(max_vertices + 1);
    anchors = arena.allocate_n<math::PixelCoords>(MAX_ANCHOR_INTERVALS + 1);
    pole_split.init_storage(arena, POLE_RUN_POINTS);
  }

  /** @brief Chart x coordinates of the next polyline, one per vertex. */
  float *chart_x() { return chart_x_storage; }
  /** @brief Chart y coordinates of the next polyline, one per vertex. */
  float *chart_y() { return chart_y_storage; }

  /** @brief Captures the frame's clip band; call once per frame. */
  void prepare(const ClipRegion &clip) { band = ClipBand<W, H>::of(clip); }

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
  HS_HOT_FLASH_MEMBER void
  draw_closed(PipelineT &pipeline, Canvas &canvas, const Fragments &points,
              int vertices, const math::Basis &planar_basis,
              const Color4 &color, const F &fragment_shader) {
    using Geometry = math::DisplayGeometry<H>;
    HS_CHECK(vertices >= 1 && vertices <= capacity,
             "PlanarChords: %d vertices outside capacity %d", vertices,
             capacity);
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
        return azimuthal_unproject(x_chart[edge] + dx * t,
                                   y_chart[edge] + dy * t, planar_basis);
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
        return newton_unit(chart_position(t));
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

private:
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
  PlanarBandSplit<W, H> pole_split;
  float *chart_x_storage = nullptr;
  float *chart_y_storage = nullptr;
  math::PixelCoords *pixels = nullptr;
  math::PixelCoords *anchors = nullptr;
  int capacity = 0;
};

} // namespace Plot
