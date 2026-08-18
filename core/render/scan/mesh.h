/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <type_traits>
#include <utility>
#include "render/sdf.h"
#include "render/shading.h"
#include "mesh/mesh_class_types.h"
#include "mesh/mesh_state.h"
#include "color/color.h"
#include "render/filter.h"
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#include "engine/platform.h"
#include "render/scan/raster.h"

/**
 * @file mesh.h
 * @brief rasterize_face, the -O3 face-specialized scan loop, and Scan::Mesh over it.
 */

namespace Scan {

HS_O3_BEGIN
/**
 * @brief Face-specialized scan loop for the mesh path, compiled at -O3 on the
 *        -Os device image.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam PipelineT Plotting pipeline type.
 * @tparam MinimalFragment When true, only v1 (raw distance) is refreshed per
 *         pixel. Nothing else ever writes the fragment: pos, v0, v2, v3, size,
 *         age and color hold their struct defaults until the shader itself
 *         writes them, and from then on whatever the previous pixel left. The
 *         shader must therefore write frag.color on every pixel it is called
 *         for, and must not read a register it did not set — with v2 at 0 and
 *         size at 1, mesh_face_index() reports face 0 and fragment_edge_dist()
 *         returns an unnormalized distance for every face.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @details Self-contained (the shared rasterize/scan_region/process_pixel
 * kernel stays -Os; GCC reuses that existing -Os instantiation rather than
 * re-optimizing it into a region caller). A Face's column intervals are
 * row-independent unless it touches a pole, so for every other face the
 * wrap/sort/coalesce pass — per-row in scan_region — and the clip-arc
 * intersection run once; the per-pixel body mirrors process_pixel's solid path.
 * Takes no debug flag and does not read canvas.debug(): the bounding-box tint
 * would need the shared rasterizer instantiated for SDF::Face, +2.5 KB ITCM.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment = false,
          typename FragmentShaderT>
HS_NOINLINE_NOCLONE inline void
rasterize_face(PipelineT &pipeline, Canvas &canvas, const SDF::Face &shape,
               FragmentShaderT &fragment_shader) {
  const auto &cr = canvas.clip();
  const auto xc = cr.x_clip();

  int y_lo, y_hi;
  static constexpr size_t MAX_RUNS = 8;
  std::pair<int, int> runs[MAX_RUNS];
  size_t num_runs = 0;

  {
    HS_PROFILE_DEEP(raster_setup);
    auto bounds = shape.template get_vertical_bounds<H>();
    if (!clamp_rows_to_clip(bounds, cr, y_lo, y_hi))
      return;

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
  }

  auto build_runs = [&](int y) {
    HS_PROFILE_DEEP(raster_setup);
    num_runs = 0;

    // Each span wraps into <= 2 seam pieces, each clip-split into <= 2 runs.
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Face>::value,
                  "intervals must hold a Face's per-row emission");
    static_assert(MAX_RUNS >= 4 * SDF::sdf_max_spans<SDF::Face>::value,
                  "runs must hold every span seam-split and clip-split");
    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });

    // Hand-rolled rather than the shared clip_run/emit_row_runs: this consumer
    // materializes runs into a cached array instead of walking them, and either
    // shared spelling costs 272 B of ITCM here.
    auto add_run = [&](int x1, int x2) {
      auto push = [&](int a, int b) {
        if (a >= b)
          return;
        HS_CHECK(num_runs < MAX_RUNS, "rasterize_face: run buffer overflow");
        runs[num_runs++] = {a, b};
      };
      if (!xc.active) {
        push(x1, x2);
      } else if (xc.wrap) {
        push(x1, std::min(x2, xc.re));
        push(std::max(x1, xc.rs), x2);
      } else {
        push(std::max(x1, xc.rs), std::min(x2, xc.re));
      }
    };

    if (!handled) {
      add_run(0, W);
    } else if (!intervals.is_empty()) {
      bool full_row = false;
      for (const auto &iv : intervals) {
        if (iv.second - iv.first >= static_cast<float>(W)) {
          full_row = true;
          break;
        }
      }
      if (full_row) {
        add_run(0, W);
      } else {
        StaticCircularBuffer<SDF::Interval, 8> norm;
        static_assert(decltype(norm)::CAPACITY ==
                          2 * decltype(intervals)::CAPACITY,
                      "norm must hold 2 spans per input interval (seam split)");
        coalesce_spans<W>(intervals, norm);
        for (const auto &run : norm)
          add_run(static_cast<int>(run.first), static_cast<int>(run.second));
      }
    }
  };

  // A pole-touching face widens its azimuth wedge toward the pole, so its runs
  // describe one row; every other face reuses one set for the whole band.
  const bool per_row = shape.pole_touch;
  if (!per_row) {
    build_runs(y_lo);
    if (num_runs == 0)
      return;
  }

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4;
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  constexpr float pixel_width = 2.0f * PI_F / W;
  const uint32_t probe_flags = shape.probe_flags();

  SDF::DistanceResult res;
  Fragment frag;

  // 1.0002 margin keeps the sqrt-free cull strictly conservative against the
  // widest threshold a probe is tested on below.
  const float base_reject_rad = pixel_width * 1.0002f;
  const float base_reject_dsq = base_reject_rad * base_reject_rad;

  [[maybe_unused]] const float plane_stretch = report_stretch(shape);

  HS_PROFILE_DEEP(raster_scan);
  for (int y = y_lo; y <= y_hi; ++y) {
    if (per_row) {
      build_runs(y);
      if (num_runs == 0)
        continue;
    }
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    const int stride = pole_lod_run(sp);
    float block_slack = 0.0f;
    float reject_dsq = base_reject_dsq;
    if constexpr (pole_lod_blocks<SDF::Face>) {
      block_slack = pole_lod_slack<W, SDF::Face>(stride, sp) * plane_stretch;
      // the block test below reads d, so the cull must not sentinel a probe
      // inside that threshold
      const float reject_rad = (pixel_width + block_slack) * 1.0002f;
      reject_dsq = reject_rad * reject_rad;
    }

    for (size_t r = 0; r < num_runs; ++r) {
      const int rx2 = runs[r].second;
      for (int x = runs[r].first; x < rx2;) {
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        shape.template distance_with_flags<true>(p, res, reject_dsq,
                                                 probe_flags);
        const float d = res.dist;

        // Columns this one shade covers. Only a canvas-aligned block that fits
        // in the run and that the surface cannot cross qualifies; anything
        // holding an edge stays per-column, so coverage is identical either way
        // and only the shade's source column moves. The clear side additionally
        // needs the probe to carry a distance rather than the cull's sentinel;
        // an inside probe is never culled.
        int span = 1;
        if constexpr (pole_lod_blocks<SDF::Face>) {
          if (stride > 1 && x % stride == 0 && x + stride <= rx2 &&
              (d <= -pixel_width - block_slack ||
               pole_lod_block_settles<SDF::Face>(d, pixel_width, block_slack)))
            span = stride;
        }

        if (d >= pixel_width) {
          x += span;
          continue;
        }
        HS_SCAN_METRIC(hs::g_scan_metrics.shade_candidates++);
        HS_PROBE_COUNT(n_alpha);
        HS_PROBE_MARK(hs_ta);
        const float alpha = solid_coverage(d, pixel_width);
        HS_PROBE_SPAN(alpha, hs_ta);
        if (alpha <= MIN_ALPHA) {
          x += span;
          continue;
        }
        HS_PROFILE_DEEP(raster_shade);
        if constexpr (MinimalFragment) {
          frag.v1 = res.raw_dist;
        } else {
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v0 = res.t;
          frag.v1 = res.raw_dist;
          frag.v2 = 0.0f;
          frag.v3 = res.aux;
          frag.size = res.size;
          frag.age = 0;
        }
        fragment_shader(p, frag);
        if (frag.color.alpha > MIN_ALPHA) {
          const float a = frag.color.alpha * alpha;
          for (int px = x; px < x + span; ++px) {
            if constexpr (requires {
                            pipeline.plot_in_bounds(
                                canvas, px, y, frag.color.color, frag.age, a);
                          }) {
              pipeline.plot_in_bounds(canvas, px, y, frag.color.color, frag.age,
                                      a);
            } else {
              pipeline.plot(canvas, px, y, frag.color.color, frag.age, a);
            }
          }
        }
        x += span;
      }
    }
  }
}
HS_O3_END

/**
 * @brief Arena-hosts one per-mesh Face scratch buffer.
 * @param arena Scratch arena supplying the storage.
 * @return The constructed buffer, valid until the caller's ScratchScope closes.
 * @details Default-init: every field is written by SDF::Face before it is read,
 * and value-init would zero ~7.7 KB per mesh draw.
 */
HS_NOINLINE_NOCLONE inline SDF::FaceScratchBuffer *
new_face_scratch(Arena &arena) {
  return arena.make_default<SDF::FaceScratchBuffer>();
}

HS_O3_BEGIN
/**
 * @brief Rasterizes a polygonal mesh by drawing each face as an SDF::Face,
 *        threading the face index through register v2 so the shader can vary
 *        color per face.
 * @details Every face goes through rasterize_face, so canvas.debug() does not
 * tint a mesh.
 */
struct Mesh {
  /**
   * @brief Traps unless every flat face index addresses a vertex in the pool.
   * @param faces Flat per-face vertex index array.
   * @param num_indices Entries in @p faces.
   * @param num_verts Vertices the mesh carries.
   * @details SDF::Face reads the vertex pool through a std::span, whose
   * operator[] only asserts, so a stale index domain would read arbitrary
   * memory as a Vector on device. Checked once per mesh, not per face: the
   * per-face spans are cut from exactly this array.
   */
  HS_COLD_MEMBER
  static void check_face_index_domain(const uint16_t *faces, size_t num_indices,
                                      size_t num_verts) {
    for (size_t k = 0; k < num_indices; ++k)
      HS_CHECK(static_cast<size_t>(faces[k]) < num_verts,
               "mesh face index exceeds the vertex pool");
  }

  /**
   * @brief Rasterizes every face of a mesh.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; receives the face
   *                        index in register v2. v2 is a float, so the index is
   *                        exact only up to 2^24 faces.
   * @param scratch_arena Arena supplying per-face SDF::Face scratch storage.
   * @param bake Optional congruence-class bake for this mesh (null = exact
   *        path everywhere, today's behavior). When present, each face is
   *        aligned to its canonical class shape after construction and the
   *        class distance LUT is bound for the probe loop.
   * @param face_shader_setup Optional callback receiving the face index and
   *        size. When supplied, it runs once before rasterizing the face and
   *        the shader is invoked directly with only v1 refreshed per pixel —
   *        the callback never sees the fragment, so the index and size reach
   *        the shader only through what it hoists here, not through v2 and
   *        size. The fragment is per-face, so a shader on this path must write
   *        frag.color unconditionally: returning early leaves the previous
   *        pixel's color to be plotted again.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_NOINLINE_NOCLONE static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
            FragmentShaderT &fragment_shader, Arena &scratch_arena,
            const MeshOps::MeshClassBake *bake,
            FaceShaderSetupT &face_shader_setup) {
    // Once per mesh, not per face: rasterize_face indexes the phi LUT by the
    // canvas' own rows and hands SDF::Face the template H.
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    // The per-face wrapper below is itself always non-null, so the erased
    // shader it wraps has to be checked here or not at all.
    if constexpr (std::is_same_v<FragmentShaderT, FragmentShaderFn>)
      check_fragment_shader(fragment_shader);

    ScratchScope scope(scratch_arena);
    auto *scratch = new_face_scratch(scratch_arena);

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();
    size_t fi_size = mesh.get_faces_size();

    check_face_index_domain(fi, fi_size, mesh.vertices.size());

    // An empty bake (build skipped) is equivalent to none; a populated one
    // must cover every face — records are indexed by face order.
    if (bake && !bake->face_recs.is_bound())
      bake = nullptr;
    HS_CHECK(!bake || bake->face_recs.size() == num_f,
             "mesh class bake face count disagrees with the mesh");

    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      // Trap malformed mesh data: an offset/count pair disagreeing with the flat
      // index array yields an out-of-bounds span for SDF::Face. Cold per-face check.
      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");

      std::span<const Vector> verts(mesh.vertices.data(), mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      // IIFE so the HS_PROFILE scope measures Face construction alone (prvalue
      // elided in place), not the rasterize below.
      SDF::Face shape = [&] {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(verts, indices, *scratch, H + hs::H_OFFSET, H,
                         &canvas.clip());
      }();

      // Bind the face's congruence-class LUT: a vertex correlation aligns the
      // current projection to the canonical frame. Culled faces (y_min > y_max)
      // skip it; a missing LUT or degenerate alignment stays on the exact path.
      if (bake && shape.y_min <= shape.y_max) {
        const MeshOps::FaceClassRec &rec = bake->face_recs[i];
        if (rec.class_id != MeshOps::NO_CLASS) {
          const MeshOps::CongruenceClass &cls = bake->classes[rec.class_id];
          if (cls.lut.data && cls.n_verts == shape.count)
            shape.bind_class_lut(&cls.lut, cls.canon_xy, rec.vert_offset,
                                 rec.reflected != 0);
        }
      }

      if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
        auto wrapper = [&](const Vector &p, Fragment &f_in) {
          // v2 carries the face index (decoded by mesh_face_index()).
          // Exact for i < 2^24 (float mantissa); meshes never approach that face count.
          f_in.v2 = static_cast<float>(i);
          fragment_shader(p, f_in);
        };
        FragmentShaderFn raster_shader = wrapper;
        {
          HS_PROFILE(scan_mesh_raster);
          rasterize_face<W, H, PipelineT, false>(pipeline, canvas, shape,
                                                 raster_shader);
        }
      } else {
        face_shader_setup(i, shape.size);
        {
          HS_PROFILE(scan_mesh_raster);
          rasterize_face<W, H, PipelineT, true>(pipeline, canvas, shape,
                                                fragment_shader);
        }
      }
    }
  }

  /** @brief Rasterizes a mesh through the type-erased fragment shader path. */
  template <int W, int H, typename PipelineT = PipelineRef,
            typename FaceShaderSetupT = std::nullptr_t>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake = nullptr,
                   FaceShaderSetupT face_shader_setup = nullptr) {
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, face_shader_setup);
  }

  /** @brief Rasterizes a mesh with an inlinable per-face fragment shader. */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  static void
  draw_specialized(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderT &fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake,
                   FaceShaderSetupT &face_shader_setup) {
    FunctionRef<void(size_t, float)> erased_setup = face_shader_setup;
    // Cold (once per draw): the erasure hides a null setup from draw_impl's
    // nullptr_t branch, which would then call a null thunk per face.
    HS_CHECK(erased_setup, "Scan::Mesh::draw_specialized requires a non-null "
                           "face_shader_setup");
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, erased_setup);
  }

  /**
   * @brief Rasterizes a mesh shaded from a per-face blended-ramp table.
   * @tparam ShadingT Table exposing `const BakedPalette &ramp_for(size_t)`.
   * @param pipeline Filter pipeline applied on plot.
   * @param canvas Target canvas.
   * @param mesh Mesh to rasterize, already in camera space.
   * @param shading Per-face ramp table (an Animation::OpLeg::Shading).
   * @param gain Edge-distance shading gain, divided by the face size.
   * @param alpha Coverage applied to every fragment.
   * @param scratch_arena Arena backing the scan's per-face buffers.
   * @details Face-hoisted: the ramp and the gain/inradius gradient scale
   * resolve once per face, leaving a multiply, a clamp and one LUT fetch per
   * fragment.
   */
  template <int W, int H, typename PipelineT, typename ShadingT>
  static void draw_opleg_shading(PipelineT &pipeline, Canvas &canvas,
                                 const MeshState &mesh, const ShadingT &shading,
                                 float gain, float alpha,
                                 Arena &scratch_arena) {
    FacePaletteShader fragment_shader;
    fragment_shader.alpha = alpha;
    auto select_face = [&](size_t fi, float size) {
      fragment_shader.set_palette(&shading.ramp_for(fi));
      fragment_shader.scale = size > math::TOLERANCE ? gain / size : 0.0f;
    };
    draw_specialized<W, H>(pipeline, canvas, mesh, fragment_shader,
                           scratch_arena, nullptr, select_face);
  }
};
HS_O3_END

} // namespace Scan
