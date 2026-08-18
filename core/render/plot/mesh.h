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
#include "engine/constants.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "mesh/triangular_bitset.h"
#include "render/plot/raster.h"
#include "render/plot/shapes.h"

/**
 * @file mesh.h
 * @brief Plot::Mesh: the edge-deduplicated mesh wireframe stroke.
 */

namespace Plot {

/**
 * @brief Mesh drawing.
 * Registers:
 *  v0: Edge Progress t (0.0 -> 1.0) per edge
 *  v1: Cumulative Arc Length (radians) per edge
 *  v2: Edge index
 */
struct Mesh {
  /**
   * @brief Max distinct vertices the edge-dedup bitset can track.
   * @details A mesh exceeding this traps while its faces are walked: at render
   * time for draw() (per frame for MeshFeedback), or at setup for
   * extract_edges(). Sized for a TriangularBitset of 128*127/2 bits = 1016
   * bytes.
   */
  static constexpr int DEDUP_CAPACITY = 128;

  /**
   * @brief Fragments one wireframe edge can produce: its two endpoints plus a
   *        cut at every clip-band boundary it crosses.
   */
  static constexpr int EDGE_MAX_POINTS = GEODESIC_CLIP_MAX_SPLITS + 2;

  /**
   * @brief Sample, shade, and rasterize one wireframe edge.
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh Mesh supplying vertex positions.
   * @param u,v Endpoint vertex indices (assumed in bounds — the callers run the
   *            cold OOB/capacity traps before delegating here).
   * @param edge_index Value written to each fragment's v2 register.
   * @param cb Clip-band cut boundaries, resolved once by the calling draw().
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader; displaces the edge's two
   *        endpoints, which then bound one geodesic — it does not deform the
   *        edge's interior.
   * @details Shared body for both draw() overloads (face-walk and precomputed
   * edge list); keeping it in one place is why the two paths stay bit-identical.
   * The adaptive rasterizer already walks the true great-circle arc, so the edge
   * is handed over whole; under a clip it is cut at the band's boundaries
   * (geodesic_clip_splits) so the outside pieces cost nothing past their gate
   * verdict.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw_edge(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                        int u, int v, int edge_index, const ClipCutBounds &cb,
                        FragmentShaderFn fragment_shader,
                        VertexShaderRef vertex_shader) {
    Fragment fu;
    fu.pos = mesh.vertices[u];
    Fragment fv;
    fv.pos = mesh.vertices[v];

    const ClipRegion &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const ClipRegion::XClip xc = cr.x_clip();

    // A vertex shader moves the endpoints after this test, so it opts out.
    if constexpr (pipeline_hoistable_cull<PipelineT>()) {
      if (clip_active && !vertex_shader &&
          !edge_visible_in_clip<W, H>(pipeline, cr, xc, fu.pos, fv.pos,
                                      nullptr))
        return;
    }

    ScratchScope edge_guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, EDGE_MAX_POINTS);

    const GeodesicEdgeSpan es = make_geodesic_edge_span(fu.pos, fv.pos);
    bool split = false;
    if constexpr (pipeline_hoistable_cull<PipelineT>())
      split = clip_active && !vertex_shader && es.have_axis;

    if (split) {
      float ts[GEODESIC_CLIP_MAX_SPLITS];
      const int cuts = geodesic_clip_splits(fu.pos, fv.pos, es, cb, ts);
      const Vector perp = cross(es.axis, fu.pos);
      points.push_back(Line::sample_point(fu, fv, es, perp, 0.0f));
      for (int i = 0; i < cuts; ++i)
        points.push_back(Line::sample_point(fu, fv, es, perp, ts[i]));
      points.push_back(Line::sample_point(fu, fv, es, perp, 1.0f));
    } else {
      Line::sample(points, fu, fv);
    }

    for (auto &p : points)
      p.v2 = static_cast<float>(edge_index); // Edge Index
    if (vertex_shader) {
      vertex_shader(points[0]);
      vertex_shader(points.back());
    }

    if constexpr (pipeline_hoistable_cull<PipelineT>()) {
      if (clip_active) {
        uint8_t bits[EDGE_MAX_POINTS - 1];
        HS_CHECK(points.size() >= 2 && points.size() <= EDGE_MAX_POINTS);
        if (points.size() == 2 && !vertex_shader) {
          // Uncut, unshaded: the whole-edge test above already ran on it.
          bits[0] = RasterOptions::EDGE_VISIBLE;
        } else if (!gate_trail_edges<W, H>(pipeline, cr, xc, points, bits)) {
          return;
        }
        rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                        {.edge_flags = bits});
        return;
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader);
  }

  /**
   * @brief Walk a mesh's faces and invoke fn(u, v) once per unique edge.
   * @tparam MeshT Mesh type.
   * @tparam Fn Per-edge callback type.
   * @param mesh Mesh whose faces are walked for edges.
   * @param visited Caller-owned dedup bitset; cleared before walking. Held by
   *                the caller so each path picks its own arena/scope.
   * @param fn Invoked as fn(u, v) for the first occurrence of each edge.
   * @details Shared face-walk/edge-dedup loop behind both draw() and
   * extract_edges().
   */
  template <typename MeshT, typename Fn>
  static void for_each_unique_edge(const MeshT &mesh,
                                   TriangularBitset<DEDUP_CAPACITY> &visited,
                                   Fn &&fn) {
    visited.clear();

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    size_t fi_size = mesh.get_faces_size();
    size_t offset = 0;

    for (size_t i = 0; i < num_f; ++i) {
      int count = fc[i];

      // Trap malformed mesh data: an offset/count pair disagreeing with the flat
      // index array yields out-of-bounds reads. Cold per-face check.
      HS_CHECK(offset + static_cast<size_t>(count) <= fi_size,
               "mesh face span exceeds face index array");

      for (int k = 0; k < count; ++k) {
        int u = fi[offset + k];
        int v = fi[offset + (k + 1) % count];
        int small = std::min(u, v);
        int large = std::max(u, v);

        // A vertex index past the dedup bitset's capacity is a mesh-sizing bug;
        // trap at the face-walk boundary rather than drop the edge.
        HS_CHECK(large < DEDUP_CAPACITY);

        if (!visited.test_and_set(small, large))
          fn(u, v);
      }
      offset += count;
    }
  }

  /**
   * @brief Draws a mesh (wireframe).
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh The mesh to draw.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    int edge_index = 0;

    // O(1) edge dedup in a 1016-byte triangular bit matrix, arena-allocated (deep
    // render chain, tight DTCM stack). Held in scratch_arena_b so the per-edge
    // scratch_arena_a scopes below keep their headroom.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *scratch_arena_b.make<TriangularBitset<DEDUP_CAPACITY>>();

    const ClipRegion &cr = canvas.clip();
    const ClipCutBounds cb = make_clip_cut_bounds<W, H>(cr, cr.x_clip());

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      // mesh.vertices[] only asserts in bounds (stripped on device), so guard the
      // per-edge setup boundary here. u,v come from uint16_t face data (non-
      // negative), so max(u,v) in bounds implies both endpoints are valid.
      HS_CHECK(static_cast<size_t>(std::max(u, v)) < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, u, v, edge_index, cb,
                      fragment_shader, vertex_shader);

      edge_index++;
    });
  }

  /**
   * @brief Draws a mesh (wireframe) without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh The mesh to draw.
   * @param fragment_shader Shader function.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, mesh, fragment_shader, {});
  }

  /**
   * @brief Precomputed edge pair for static-topology meshes.
   */
  struct Edge {
    uint16_t u, v; /**< Endpoint vertex indices into the mesh's vertex array. */
  };

  /**
   * @brief Extract unique edges from a mesh (call once at setup time).
   * @tparam MeshT Mesh type.
   * @param mesh Mesh whose faces are walked for unique edges.
   * @param edges Output edge list; deduplicated edges are appended.
   */
  template <typename MeshT>
  static void extract_edges(const MeshT &mesh, ArenaVector<Edge> &edges) {
    // Dedup bitset (1016 B) in the arena, not the stack (deep setup chain). The
    // output `edges` lives in a separate persistent arena, so scratch_arena_b
    // cannot disturb it.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *scratch_arena_b.make<TriangularBitset<DEDUP_CAPACITY>>();

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      edges.push_back({(uint16_t)u, (uint16_t)v});
    });
  }

  /**
   * @brief Draw using a precomputed edge list (skips face walk + dedup).
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh Mesh supplying vertex positions.
   * @param edges Precomputed unique edge list.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param mask Optional dissolve ownership mask (see DissolveMask), keyed on
   *        the edge's endpoint indices; an unowned edge is skipped before any
   *        of its geometry work, so two complementary masks split one
   *        wireframe's cost across the two meshes of a transition.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void
  draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
       const ArenaVector<Edge> &edges, FragmentShaderFn fragment_shader,
       VertexShaderRef vertex_shader = {}, const DissolveMask *mask = nullptr) {
    const ClipRegion &cr = canvas.clip();
    const ClipCutBounds cb = make_clip_cut_bounds<W, H>(cr, cr.x_clip());

    for (size_t ei = 0; ei < edges.size(); ++ei) {
      if (mask && !mask->owns(edges[ei].u, edges[ei].v))
        continue;
      // Setup-boundary OOB guard (see the face-walk overload above): the raw
      // edge list could outlive or mismatch its mesh, and mesh.vertices[] only
      // asserts (compiled out on device).
      HS_CHECK(edges[ei].u < mesh.vertices.size() &&
               edges[ei].v < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, edges[ei].u, edges[ei].v,
                      static_cast<int>(ei), cb, fragment_shader, vertex_shader);
    }
  }
};

} // namespace Plot
