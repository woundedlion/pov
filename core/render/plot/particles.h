/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include "math/geometry.h"
#include "render/shading.h"
#include "color/color.h"
#include "platform/constants.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "mesh/triangular_bitset.h"
#include "render/plot/raster.h"
#include "render/plot/shapes.h"

/**
 * @file particles.h
 * @brief Plot::ParticleSystem: the trail stroke over a particle system.
 */

namespace Plot {

/**
 * @brief Particle System trails.
 * Registers:
 *  v0: Trail Progress (0.0=Tail -> 1.0=Head)
 *  v1: Reserved (always 0)
 *  v2: Particle ID, or the particle_v2 mapper's value
 *  v3: Normalized TTL
 */
struct ParticleSystem {
  template <typename SystemT> static consteval int trail_sample_stride() {
    if constexpr (requires { SystemT::TRAIL_SAMPLE_STRIDE; })
      return SystemT::TRAIL_SAMPLE_STRIDE;
    return 1;
  }

  /**
   * @brief Draws each active particle's history as a rasterized trail.
   * @tparam W,H Rasterization resolution.
   * @tparam HoistableCull True when raw point projections are valid for clip
   *         gating because the source pipeline has no world cull stage.
   * @tparam FuseVertex Apply the typed vertex shader as each point is emitted.
   * @tparam FragmentShaderT Fragment shader type.
   * @tparam VertexShaderFn Vertex shader type.
   * @tparam DeferredShaderT Deferred vertex shader type.
   * @tparam SystemT Particle-system type.
   * @tparam ParticleV2Fn Per-particle v2 mapper type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader (position pass).
   * @param deferred_shader Optional second vertex pass, given each fragment and
   *        its original pre-shader position. Under an active clip it runs only
   *        for trails with at least one cull-surviving edge; a skipped trail
   *        renders nothing, so output is identical to an undeferred shader.
   *        Put per-point work that only affects shading registers here.
   * @param particle_v2 Optional mapper from particle and pool index to v2.
   *        The default stores the pool index. Contract: called exactly once per
   *        materialized particle, before that particle's fragments are emitted,
   *        shaded, or rasterized, and before the next particle's call. A mapper
   *        may therefore bind per-particle state the fragment or deferred shader
   *        reads. Particles the traversal skips (drained history, fewer than two
   *        points) are never mapped and never shaded.
   */
  template <int W, int H, bool HoistableCull, bool FuseVertex,
            bool SinglePassRaster, typename PipelineT, typename SystemT,
            typename FragmentShaderT, typename VertexShaderFn,
            typename DeferredShaderT, typename ParticleV2Fn>
  static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const SystemT &system,
            FragmentShaderT fragment_shader, VertexShaderFn vertex_shader,
            DeferredShaderT deferred_shader, ParticleV2Fn particle_v2) {
    int count = system.active();
    if (count == 0)
      return;

    const float max_life = static_cast<float>(system.max_life);
    HS_CHECK(std::isfinite(max_life) && max_life >= 1.0f &&
                 max_life <= 65535.0f,
             "ParticleSystem render max_life must be finite and in [1, 65535]");
    const float inv_max_life = 1.0f / max_life;
    const bool has_deferred_shader = [&] {
      if constexpr (std::is_same_v<std::decay_t<DeferredShaderT>,
                                   std::nullptr_t>)
        return false;
      else if constexpr (requires { static_cast<bool>(deferred_shader); })
        return static_cast<bool>(deferred_shader);
      else
        return true;
    }();

    // Segment-clip state for the trail-level deferred-shader gate below.
    const auto &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const auto xc = cr.x_clip();
    CartesianQuadrantClip cartesian_clip;
    if constexpr (HoistableCull)
      cartesian_clip = make_cartesian_quadrant_clip<W, H>(cr);

    for (int i = 0; i < count; ++i) {
      const auto &p = system.pool[i];
      const size_t trail_len = p.history.length();
      constexpr int sample_stride = trail_sample_stride<SystemT>();
      const bool append_live_tip = [&] {
        if constexpr (sample_stride == 1) {
          return false;
        } else {
          // A max_life lowered below a live particle's remaining life would
          // underflow this uint16_t subtraction.
          const uint16_t age =
              p.life < system.max_life ? system.max_life - p.life : 0;
          return (age - 1) % sample_stride != 0;
        }
      }();
      const size_t point_count = trail_len + (append_live_tip ? 1 : 0);
      HS_MSP_COUNT(resident_particles);
      if (p.life == 0)
        HS_MSP_COUNT(draining_histories);
      else {
        HS_MSP_COUNT(live_particles);
        if (trail_len == static_cast<size_t>(p.history.CAPACITY))
          HS_MSP_COUNT(full_histories);
        else
          HS_MSP_COUNT(partial_histories);
      }
      if (p.life == 0 || point_count < 2)
        continue;
      const float v2 = [&] {
        if constexpr (std::is_same_v<ParticleV2Fn, std::nullptr_t>)
          return static_cast<float>(i);
        else
          return particle_v2(p, i);
      }();
      const float particle_life =
          std::min(static_cast<float>(p.life) * inv_max_life, 1.0f);
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a, point_count);
      // Original (pre-shader) positions, kept for the deferred pass.
      ArenaVector<Vector> orig;
      if (has_deferred_shader)
        orig.bind(scratch_arena_a, point_count);
      {
        HS_PROFILE(plot_ps_tween);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        hs::DwtStallBatch history_batch(
            hs::g_mindsplatter_stalls.history_vertex);
#endif
        auto emit = [&](const Vector &v, float t) {
          trail.emplace_back(
              Fragment{.pos = v, .v0 = t, .v2 = v2, .v3 = particle_life});
          if constexpr (FuseVertex)
            vertex_shader(trail.back());
          if (has_deferred_shader)
            orig.push_back(v);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
          history_batch.step();
#endif
        };
        if constexpr (sample_stride == 1) {
          p.history.tween(emit);
        } else {
          const float denominator = static_cast<float>(point_count - 1);
          for (size_t j = 0; j < trail_len; ++j)
            emit(p.history.get(j), static_cast<float>(j) / denominator);
          if (append_live_tip)
            emit(p.position, 1.0f);
        }
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        history_batch.finish();
#endif
      }

      if (trail.is_empty())
        continue;
      if constexpr (!FuseVertex) {
        HS_PROFILE(plot_ps_vertex);
        apply_vertex_shader(vertex_shader, trail);
      }

      // Trail-level gate: precompute each edge's cull verdict from the
      // position-shaded points. No visible edge means the trail renders
      // nothing, so the optional deferred pass and the rasterize call are
      // skipped whole; the bits feed rasterize so the cull is evaluated once.
      const uint8_t *vis = nullptr;
      const float *dot_rows = nullptr;
      const float *dot_cols = nullptr;
      if (clip_active && trail.size() >= 2) {
        HS_PROFILE(plot_ps_gate);
        const size_t edges = trail.size() - 1;
        auto *bits = static_cast<uint8_t *>(
            scratch_arena_a.allocate(edges, alignof(uint8_t)));
        bool any = false;
        if constexpr (HoistableCull) {
          CartesianTrailGateResult cartesian_result;
          {
            HS_PROFILE(plot_ps_cartesian_gate);
            cartesian_result =
                cartesian_quadrant_trail_gate(cartesian_clip, trail);
          }
          count_cartesian_trail_gate_result(cartesian_result);
          if (cartesian_result != CartesianTrailGateResult::EXACT_FALLBACK)
            continue;

          // No stage re-emits edges, so the predicate sees the raw points:
          // per-point rows/columns are computed once and shared by every edge,
          // and a conservative whole-trail bound rejects fully-invisible
          // trails before any per-edge work.
          const TrailGatePrologue pro =
              trail_gate_prologue<W, H>(cr, xc, trail);
          if (pro.rejected)
            continue;
          const float *rows = pro.rows;
          const float *cols = pro.cols;
          // rasterize's single-dot shortcut takes both projections or neither.
          if (cols != nullptr) {
            dot_rows = rows;
            dot_cols = cols;
          }

          for (size_t e = 0; e < edges; ++e) {
            HS_MSP_STALL_START(edge_gate_start);
            const Vector &ea = trail[e].pos;
            const Vector &eb = trail[e + 1].pos;
            const bool one_dot = edge_fits_one_dot<W, H>(ea, eb);
            count_particle_edge_class(one_dot);
            if (one_dot) {
              bool v = antialiased_dot_visible_in_clip<W, H>(
                  cr, xc, rows[e], cols != nullptr ? cols[e] : 0.0f);
              if (e + 1 == edges)
                v = v || antialiased_dot_visible_in_clip<W, H>(
                             cr, xc, rows[e + 1],
                             cols != nullptr ? cols[e + 1] : 0.0f);
              bits[e] = RasterOptions::EDGE_CLASSIFIED |
                        RasterOptions::EDGE_ONE_DOT |
                        (v ? RasterOptions::EDGE_VISIBLE : uint8_t{0});
              if (!v)
                HS_MSP_COUNT(edge_rejects);
              any = any || v;
              HS_MSP_STALL_STOP(trail_gate, edge_gate_start);
              continue;
            }
            bool v;
            const RawGeodesicGateResult raw = raw_geodesic_edge_gate<W, H>(
                cr, xc, rows[e], rows[e + 1], cols != nullptr ? cols[e] : 0.0f,
                cols != nullptr ? cols[e + 1] : 0.0f, ea, eb);
            if (raw != RawGeodesicGateResult::EXACT_FALLBACK) {
              v = raw == RawGeodesicGateResult::VISIBLE;
            } else {
              count_particle_exact_gate_fallback();
              const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
              v = exact_geodesic_edge_visible_hoisted<W, H>(cr, xc, rows, cols,
                                                            e, ea, eb, es);
            }
            bits[e] = RasterOptions::EDGE_CLASSIFIED |
                      (v ? RasterOptions::EDGE_VISIBLE : uint8_t{0});
            if (!v)
              HS_MSP_COUNT(edge_rejects);
            any = any || v;
            HS_MSP_STALL_STOP(trail_gate, edge_gate_start);
          }
        } else {
          for (size_t e = 0; e < edges; ++e) {
            bits[e] = edge_visible_in_clip<W, H>(pipeline, cr, xc, trail[e].pos,
                                                 trail[e + 1].pos, nullptr)
                          ? RasterOptions::EDGE_VISIBLE
                          : uint8_t{0};
            any = any || bits[e] != 0;
          }
        }
        if (!any)
          continue;
        HS_MSP_COUNT(visible_trails);
        vis = bits;
      }
      if (!clip_active)
        HS_MSP_COUNT(visible_trails);

      if (has_deferred_shader) {
        HS_PROFILE(plot_ps_deferred);
        for (size_t k = 0; k < trail.size(); ++k)
          deferred_shader(trail[k], orig[k]);
      }
      {
        HS_PROFILE(plot_ps_raster);
        rasterize<W, H,
                  RasterConfig{
                      .single_pass = SinglePassRaster && sample_stride == 1,
                      .open_geodesic = true}>(pipeline, canvas, trail,
                                              fragment_shader,
                                              {.edge_flags = vis,
                                               .point_rows = dot_rows,
                                               .point_cols = dot_cols});
      }
    }
  }

  /**
   * @brief Draws particle trails through the shared raster surface.
   * @details The source pipeline's static cull trait is retained separately so
   *          Cartesian, one-dot and raw-edge gates remain available after plot
   *          dispatch is erased. A pipeline explicitly declaring the direct
   *          raster path retains its compile-time plot calls.
   */
  template <int W, int H, typename PipelineT = PipelineRef,
            typename ParticleV2Fn = std::nullptr_t>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader,
                   DeferredShaderRef deferred_shader = {},
                   ParticleV2Fn particle_v2 = nullptr) {
    if constexpr (pipeline_direct_raster_path<PipelineT>()) {
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), false, false>(
          pipeline, canvas, system, fragment_shader, vertex_shader,
          deferred_shader, particle_v2);
    } else {
      PipelineRef erased(pipeline);
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), false, false>(
          erased, canvas, system, fragment_shader, vertex_shader,
          deferred_shader, particle_v2);
    }
  }

  /**
   * @brief Draws particle trails while applying a typed vertex shader during
   *        trail materialization.
   * @details Fuses point materialization and transformation into one traversal.
   * @tparam SinglePassRaster Emit adaptive geodesic samples without a replay
   *         cache. Applies only to a direct raster pipeline.
   * @tparam DeferredShaderT Deferred vertex shader type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   * @param vertex_shader Typed vertex shader.
   * @param deferred_shader Optional deferred vertex shader.
   * @param particle_v2 Optional particle-to-v2 mapper; see draw_impl for the
   *        per-particle ordering the mapper may rely on.
   */
  template <int W, int H, bool SinglePassRaster = false, typename PipelineT,
            typename FragmentShaderT, typename VertexShaderFn,
            typename DeferredShaderT = DeferredShaderRef,
            typename ParticleV2Fn = std::nullptr_t>
  static void draw_fused_vertex(PipelineT &pipeline, Canvas &canvas,
                                const auto &system,
                                FragmentShaderT fragment_shader,
                                VertexShaderFn vertex_shader,
                                DeferredShaderT deferred_shader = {},
                                ParticleV2Fn particle_v2 = nullptr) {
    if constexpr (pipeline_direct_raster_path<PipelineT>()) {
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), true,
                SinglePassRaster>(pipeline, canvas, system, fragment_shader,
                                  vertex_shader, deferred_shader, particle_v2);
    } else {
      PipelineRef erased(pipeline);
      FragmentShaderFn erased_shader(fragment_shader);
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), true, false>(
          erased, canvas, system, erased_shader, vertex_shader, deferred_shader,
          particle_v2);
    }
  }

  /**
   * @brief Draws particle trails without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, system, fragment_shader, {});
  }
};

} // namespace Plot
