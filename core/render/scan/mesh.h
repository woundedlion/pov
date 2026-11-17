/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstring>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include "core/render/sdf.h"
#include "core/render/shading.h"
#include "core/render/sdf/face_classes.h"
#include "core/mesh/mesh_state.h"
#include "core/color/color.h"
#include "core/render/filter/pipeline.h"
#include "core/containers/static_circular_buffer.h"
#include "core/render/canvas.h"
#include "core/platform/platform.h"
#include "core/render/scan/raster.h"
#include "core/render/scan/face_raster.h"
#include "core/render/scan/swept_face.h"
#include "core/render/scan/mesh_capture.h"
#include "core/render/spacetime/crossing.h"

/**
 * @file mesh.h
 * @brief Mesh instance scan dispatch.
 */

namespace Scan {

HS_O3_BEGIN
/**
 * @brief Rasterizes a polygonal mesh by drawing each face as an SDF::PreparedFace.
 * @details Without a face-shader setup callback the face index reaches the
 * shader in register v2; with one, it reaches the shader only through what
 * the callback hoists.
 * canvas.debug() does not tint a mesh. Under the canvas's column-time basis
 * each face's bounds are swept by its arc and every pixel is evaluated along
 * its column's rotated direction. An open TimeWindowScope limits the scan to its
 * columns; in a deformation scan, Scan::ColumnBasis::mesh.selected_faces
 * selects the faces each pass draws.
 */
struct Mesh {
private:
  /** @brief Binds a face index to the shader's v2 input. */
  template <typename Shader>
  HS_HOT_INLINE static auto face_index_shader(size_t index, Shader &shader) {
    return [index, &shader](const math::Vector &p, Fragment &f) {
      f.v2 = static_cast<float>(index);
      shader(p, f);
    };
  }

  /** @brief Draw-scoped setup for admitting and rasterizing scan faces. */
  template <int W, int H> struct FaceScan {
    const Scan::ColumnBasis *time;
    float sweep;
    float bounds_margin;
    FaceProjectionColumn *projections;
    FacePadTables pads;
    const ClipRegion &clip;
    const Spacetime::TimeWindow &window;
    FaceClipPlanes<W, H> clip_planes;

    /** @brief Reserves the basis-dependent scratch for a scan. */
    template <typename PipelineT>
    HS_HOT_INLINE FaceScan(PipelineT &pipeline, Canvas &canvas, Arena &arena,
                           bool reserve_projections)
        : time(column_basis(pipeline, canvas)),
          sweep(time && time->azimuth ? time->sweep_arc : 0.0f),
          bounds_margin(scan_bounds_margin<W>(sweep)),
          projections(reserve_projections && sweep > 0.0f
                          ? new_face_projection_columns<W>(arena)
                          : nullptr),
          pads(FacePadTables::prepare<W, H>(arena, time)),
          clip(source_clip<W, H>(pipeline, canvas)),
          window(canvas.time_window()),
          clip_planes(clip, bounds_margin + SDF::MOVED_AXIS_SLACK) {}

    /** @brief Prepares a face if its projected bounds meet the scan window.
     * @tparam EARLY_LATITUDE Rejects out-of-band faces before pose setup.
     * @param vertices Vertex pool.
     * @param indices Vertex indices of the face.
     * @param storage Reused face storage for this draw.
     * @return Prepared face handle, or nullopt if culled. */
    template <bool EARLY_LATITUDE>
    HS_HOT_INLINE std::optional<SDF::Face>
    admit(std::span<const math::Vector> vertices,
          std::span<const uint16_t> indices, SDF::FaceStorage &storage) const {
      const math::LatitudeGeometry latitude(H);
      if constexpr (EARLY_LATITUDE)
        if (SDF::PreparedFace::phi_band_rejects(vertices, indices, latitude, H,
                                                bounds_margin, &clip))
          return std::nullopt;
      SDF::FaceBoundsPose pose;
      bool culled;
      const bool seeing =
          face_bounds_pose<W, H>(vertices, indices, time, pose, pads.compact,
                                 pads.wide, clip_planes, culled);
      if (culled)
        return std::nullopt;
      SDF::Face face = [&]() HS_HOT_INLINE {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(vertices, indices, storage, latitude, H, &clip,
                         pads.full, bounds_margin, seeing ? &pose : nullptr,
                         EARLY_LATITUDE);
      }();
      const SDF::PreparedFace &shape = face.prepared();
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.face_setups;
#endif
      if (shape.culled() ||
          (window.active && !face_meets_window(shape, window, clip)))
        return std::nullopt;
      return face;
    }
  };

public:
  /**
   * @brief Shades and rasterizes one constructed face.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shape The face.
   * @param i The face's index in its mesh.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param face_shader_setup As draw_impl() takes it.
   * @param time The scan's column-time state; read when @p sweep is positive.
   * @param sweep Arc the scan's column basis turns over the frame, or 0 for
   *        the fixed meridian.
   * @param projections Draw-scoped column projection storage.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_INLINE static void raster_face(
      PipelineT &pipeline, Canvas &canvas, const SDF::PreparedFace &shape,
      size_t i, FragmentShaderT &fragment_shader,
      FaceShaderSetupT &face_shader_setup, const Scan::ColumnBasis *time,
      float sweep, FaceProjectionColumn *projections) {
#if HS_ENABLE_TEST_HOOKS
    face_audit_index = i;
#endif
    if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
      auto wrapper = face_index_shader(i, fragment_shader);
      FragmentShaderRef raster_shader = wrapper;
      {
        HS_PROFILE(scan_mesh_raster);
        if (sweep > 0.0f)
          rasterize_face_column_time<W, H, PipelineT, false>(
              pipeline, canvas, shape, raster_shader, *time, projections);
        else
          rasterize_face<W, H, PipelineT, false>(pipeline, canvas, shape,
                                                 raster_shader);
      }
    } else {
      face_shader_setup(i, shape.inradius());
      {
        HS_PROFILE(scan_mesh_raster);
        if (sweep > 0.0f)
          rasterize_face_column_time<W, H, PipelineT, true>(
              pipeline, canvas, shape, fragment_shader, *time, projections);
        else
          rasterize_face<W, H, PipelineT, true>(pipeline, canvas, shape,
                                                fragment_shader);
      }
    }
  }

  /**
   * @brief Rasterizes one face a deformation scan draws alone, from its
   *        copy's own vertices.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas, its render-x window open on the copy's
   *        arm.
   * @param copy The face copy.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param arena Arena for the face's scratch.
   * @param face_shader_setup As draw_impl() takes it.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  draw_face_copy(PipelineT &pipeline, Canvas &canvas,
                 const Scan::FaceArmCopy &copy,
                 FragmentShaderT &fragment_shader, Arena &arena,
                 FaceShaderSetupT &face_shader_setup) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ScratchScope scope(arena);
    SDF::FaceStorage &storage = *new_face_storage(arena);
    const FaceScan<W, H> scan(pipeline, canvas, arena, false);
    auto face = scan.template admit<false>({copy.vertices, copy.count},
                                           {copy.indices, copy.count}, storage);
    if (!face)
      return;
    auto *projections =
        scan.sweep > 0.0f ? new_face_projection_columns<W>(arena) : nullptr;
    raster_face<W, H>(pipeline, canvas, face->prepared(), copy.face,
                      fragment_shader, face_shader_setup, scan.time, scan.sweep,
                      projections);
  }

  /**
   * @brief Rasterizes a mesh, or the one face a deformation scan draws
   *        alone or sweeps.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh The mesh.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   * @param face_shader_setup As draw_impl() takes it.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_INLINE static void
  dispatch(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
           FragmentShaderT &fragment_shader, Arena &scratch_arena,
           const MeshOps::MeshClassBake *bake,
           FaceShaderSetupT &face_shader_setup) {
    if constexpr (std::is_same_v<FragmentShaderT, FragmentShaderRef>)
      check_fragment_shader(fragment_shader);
    const Scan::ColumnBasis *own = canvas.column_basis();
    if (own)
      check_mesh_pass(own->mesh, mesh.num_faces(), mesh.vertices.size());
    if (own && own->mesh.face_sweep)
      draw_face_sweep<W, H>(pipeline, canvas, *own->mesh.face_sweep,
                            fragment_shader, scratch_arena, face_shader_setup);
    else if (own && own->mesh.face_copy)
      draw_face_copy<W, H>(pipeline, canvas, *own->mesh.face_copy,
                           fragment_shader, scratch_arena, face_shader_setup);
    else
      draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                      bake, face_shader_setup);
  }

  /**
   * @brief Rasterizes a face a deformation scan sweeps, one column pair at a
   *        time, each pair's face standing as at its sweep time.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param sweep The face and its vertex paths.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param arena Arena for the face's fits and per-column state.
   * @param face_shader_setup As draw_impl() takes it, run per column pair
   *        with the face's size there.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  draw_face_sweep(PipelineT &pipeline, Canvas &canvas,
                  const Scan::PoleFaceSweep &sweep,
                  FragmentShaderT &fragment_shader, Arena &arena,
                  FaceShaderSetupT &face_shader_setup) {
    HS_PROFILE(scan_face_sweep);
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ScratchScope scope(arena);
    SweptFace<W, H> face(sweep, arena);
#if HS_ENABLE_TEST_HOOKS && defined(HS_SCAN_METRICS)
    const size_t near = hs::g_scan_metrics.shade_candidates;
#endif
    if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
      auto wrapper = face_index_shader(sweep.face, fragment_shader);
      FragmentShaderRef raster_shader = wrapper;
      scan_swept_face<W, H>(pipeline, canvas, sweep, face, raster_shader,
                            face_shader_setup);
    } else {
      scan_swept_face<W, H>(pipeline, canvas, sweep, face, fragment_shader,
                            face_shader_setup);
    }
#if HS_ENABLE_TEST_HOOKS && defined(HS_SCAN_METRICS)
    mesh_capture_counters.sweep_pixels +=
        hs::g_scan_metrics.shade_candidates - near;
#endif
  }

  /**
   * @brief Scans a swept face column pair by column pair.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param sweep The face and its vertex paths.
   * @param face The face, fitted over the frame.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param face_shader_setup As draw_impl() takes it, run per column pair
   *        with the face's size there.
   * @details Each column pair the clip admits on an arm the face can reach
   * (SweepWedge) stands the face at its sweep time (SweptFace::column()),
   * which bounds the pair's rows; each admitted arm's rows probe the
   * standing face along the pair's rotated basis and shade as the face
   * rasterizer shades, a time-aware shader at the pair's sweep time.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  scan_swept_face(PipelineT &pipeline, Canvas &canvas,
                  const Scan::PoleFaceSweep &sweep, SweptFace<W, H> &face,
                  FragmentShaderT &fragment_shader,
                  FaceShaderSetupT &face_shader_setup) {
    const Scan::ColumnBasis &time = *column_basis(pipeline, canvas);
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    const ClipPieces clip = clip_pieces<W>(cr.x_clip());
    const auto admits = [&](int column) HS_HOT_INLINE {
      for (int r = 0; r < clip.count; ++r)
        if (column >= clip.lo[r] && column < clip.hi[r])
          return true;
      return false;
    };
    const Scan::SweepWedge wedge =
        Scan::SweepWedge::of(sweep, SWEEP_ROW_REACH<W> + time.sweep_arc);
    constexpr bool TIMED_SHADER =
        Spacetime::reads_time<FragmentShaderT, const math::Vector &,
                              Fragment &>;
    SDF::DistanceResult res;
    Fragment frag;
    for (int x = 0; x < W / 2; ++x) {
      const int arms =
          wedge.arms(time.azimuth[x], admits(x), admits(x + W / 2));
      if (!arms)
        continue;
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.sweep_tests;
#endif
      SDF::Bounds rows[2];
      bool stands;
      {
        HS_PROFILE(scan_sweep_column);
        stands = face.column(time, x, arms, rows);
      }
      if (!stands)
        continue;
      if constexpr (!std::is_same_v<FaceShaderSetupT, std::nullptr_t>)
        face_shader_setup(sweep.face, face.size());
      // A time-aware shader shades at the pair's time while its state moves.
      [[maybe_unused]] const bool per_column =
          shades_per_column<TIMED_SHADER>(fragment_shader);
      const math::Vector &e = time.azimuth[x];
      const math::Vector &v = time.up[x];
      [[maybe_unused]] bool counted = false;
      for (int arm = 0; arm < 2; ++arm) {
        int y_lo, y_hi;
        if (rows[arm].y_min > rows[arm].y_max ||
            !clamp_rows_to_clip(rows[arm], cr, y_lo, y_hi))
          continue;
#if HS_ENABLE_TEST_HOOKS
        mesh_capture_counters.sweep_columns += !counted;
        counted = true;
        mesh_capture_counters.sweep_rows +=
            static_cast<size_t>(y_hi - y_lo + 1);
#endif
        const int column = x + arm * (W / 2);
        const float side = arm == 0 ? 1.0f : -1.0f;
        for (int y = y_lo; y <= y_hi; ++y) {
          const float sp = side * math::TrigLUT<W, H>::sin_phi[y];
          const float cp = math::TrigLUT<W, H>::cos_phi[y];
          const math::Vector p = basis_direction(e, v, sp, cp);
#if HS_ENABLE_TEST_HOOKS
          ++mesh_capture_counters.sweep_probes;
#endif
          face.probe(p, res);
          float alpha;
          if (solid_admit<W>(res.dist, alpha))
            shade_solid_span<W, TIMED_SHADER, false, true>(
                pipeline, canvas, column, x, 1, y, p, res, frag,
                fragment_shader, alpha, per_column);
        }
      }
    }
  }

  /**
   * @brief Rasterizes every face of a mesh.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type (deduced; pass a Scan::ScanRef to erase the sink).
   * @tparam FragmentShaderT Per-pixel shader type.
   * @tparam FaceShaderSetupT Optional per-face shader setup callable.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; without
   *                        @p face_shader_setup it receives the face index in
   *                        register v2. v2 is a float, so the index is exact
   *                        only up to 2^24 faces.
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh (null = exact
   *        distance path for every face). When present, each face is
   *        aligned to its canonical class shape after construction and the
   *        class distance LUT is bound for the probe loop.
   * @param face_shader_setup Optional callback receiving the face index and
   *        size, run before each scan of a face, possibly several times per frame
   *        and per column pair for a swept face. When supplied, the shader is
   *        invoked directly with only v1 refreshed per pixel, so the index and
   *        size reach it only through what the callback hoists. The fragment
   *        persists across pixels, so the shader must write frag.color
   *        unconditionally.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_NOINLINE_NOCLONE static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
            FragmentShaderT &fragment_shader, Arena &scratch_arena,
            const MeshOps::MeshClassBake *bake,
            FaceShaderSetupT &face_shader_setup) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    ScratchScope scope(scratch_arena);
    auto *storage = new_face_storage(scratch_arena);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const FaceScan<W, H> scan(pipeline, canvas, scratch_arena, true);

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();
    size_t fi_size = mesh.get_faces_size();

    detail::check_face_offsets(mesh);

    check_face_index_domain(fi, fi_size, mesh.vertices.size());

    // An unbound bake is treated as none; records are indexed by face order.
    if (bake && !bake->face_recs.is_bound())
      bake = nullptr;
    HS_CHECK(!bake || bake->face_recs.size() == num_f,
             "mesh class bake face count disagrees with the mesh");

    const uint8_t *selected_faces =
        scan.time ? scan.time->mesh.selected_faces : nullptr;
    const uint8_t drawn_face = scan.window.active
                                   ? Scan::detail::FacePass::WINDOWED
                                   : Scan::detail::FacePass::SHARED;
    for (size_t i = 0; i < num_f; ++i) {
      if (scan.time && !mesh_visible(scan.time->mesh.visible_faces, i))
        continue;
      size_t count = fc[i];

      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");
      if (selected_faces && selected_faces[i] != drawn_face)
        continue;

      std::span<const math::Vector> verts(mesh.vertices.data(),
                                          mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      auto face = scan.template admit<true>(verts, indices, *storage);
      if (!face)
        continue;
      const SDF::PreparedFace &shape = face->prepared();

      if (bake) {
        const MeshOps::FaceClassRec &rec = bake->face_recs[i];
        if (rec.class_id != MeshOps::NO_CLASS) {
          HS_CHECK(rec.class_id < bake->classes.size(),
                   "mesh class bake face record names an unknown class");
          const MeshOps::CongruenceClass &cls = bake->classes[rec.class_id];
          if (cls.lut.data && cls.n_verts == shape.vertex_count())
            face->bind_class_lut(&cls.lut, cls.canon_xy, rec.vert_offset,
                                 rec.reflected != 0);
        }
      }

      raster_face<W, H>(pipeline, canvas, shape, i, fragment_shader,
                        face_shader_setup, scan.time, scan.sweep,
                        scan.projections);
    }
  }

  /**
   * @brief Rasterizes a mesh through the type-erased fragment shader path.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; receives the
   *        face index in register v2.
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   */
  template <int W, int H, typename PipelineT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderRef fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake = nullptr) {
    std::nullptr_t face_shader_setup = nullptr;
    dispatch<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena, bake,
                   face_shader_setup);
  }

  /**
   * @brief Rasterizes a mesh with an inlinable per-face fragment shader.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT Callable `void(size_t face, float size)`.
   * @param fragment_shader Shader invoked per covered pixel; must write
   *        frag.color unconditionally.
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   * @param face_shader_setup Non-null callback run before each face scan,
   *        possibly several times per frame and per column pair for swept faces.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  static void
  draw_specialized(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderT &fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake,
                   FaceShaderSetupT &face_shader_setup) {
    FunctionRef<void(size_t, float)> erased_setup = face_shader_setup;
    // the erasure hides a null setup from draw_impl's nullptr_t branch
    HS_CHECK(erased_setup, "Scan::Mesh::draw_specialized requires a non-null "
                           "face_shader_setup");
    dispatch<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena, bake,
                   erased_setup);
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
   * @details The ramp and gradient scale resolve once per face. A table
   * carrying its crossfade endpoints (`bank` and `blend_face`) moves each
   * column-time column along the crossfade from the weight drawn a frame ago.
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
      if constexpr (requires { shading.blend_face(fragment_shader, fi); }) {
        if (shading.bank)
          shading.blend_face(fragment_shader, fi);
      }
      fragment_shader.scale = size > math::TOLERANCE ? gain / size : 0.0f;
    };
    draw_specialized<W, H>(pipeline, canvas, mesh, fragment_shader,
                           scratch_arena, nullptr, select_face);
  }
};
HS_O3_END

} // namespace Scan
