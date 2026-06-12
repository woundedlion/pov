/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Functional output tests for mesh rasterization — the backbone of the
 * Conway/Hankin/IslamicStars effect family. These draw a known Platonic solid
 * all the way into a live Canvas and assert on the lit pixels:
 *
 *   - Plot::Mesh::draw (wireframe): every unique edge's projected geodesic
 *     midpoint is lit, so a dropped face-walk edge or a broken dedup surfaces
 *     here rather than as a subtle missing wireframe line.
 *   - Scan::Mesh::draw / SDF::Face (solid fill): each face's interior (centroid)
 *     is lit and a closed convex solid's faces tile the whole sphere with no
 *     holes — i.e. the per-face row/column bounding cull is conservative and
 *     never clips a row the face actually covers.
 */
#pragma once

#include <cmath>

#include "core/canvas.h"
#include "core/geometry.h"
#include "core/mesh.h"
#include "core/memory.h"
#include "core/plot.h"
#include "core/scan.h"
#include "core/solids.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_raster_tests {

// Dedicated arenas, kept alive for the whole test so the built mesh's
// ArenaVectors stay valid while it is drawn. The global scratch_arena_a that
// Plot::Mesh::draw uses internally for per-edge line sampling is a *different*
// arena (configure_arenas_default()), so the source mesh is never clobbered.
inline uint8_t mr_seed_a[256 * 1024];
inline uint8_t mr_seed_b[256 * 1024];
inline uint8_t mr_geom[256 * 1024];
inline uint8_t mr_scratch[256 * 1024];

// Minimal Effect that owns a Canvas to draw into; no background, no per-frame
// work, so the lit pixels come solely from the mesh draw under test.
struct MeshFx : public Effect {
  MeshFx(int W, int H) : Effect(W, H) {}
  void draw_frame() override {}
  bool show_bg() const override { return false; }
};

inline bool is_black(const Pixel &p) { return p.r == 0 && p.g == 0 && p.b == 0; }

// Fragment shader that paints a fixed near-white; used as the draw color so any
// covered pixel becomes non-black and thus detectable.
inline void white(const Vector &, Fragment &f) {
  f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
}

// A lit pixel within radius r of (px, py), wrapping in x (the canvas is a
// cylinder in longitude) and clamping in y. Tolerance absorbs sampling
// granularity, AA spread, and projection rounding.
template <int W, int H>
inline bool lit_near(const MeshFx &fx, float px, float py, int r) {
  const int cx = static_cast<int>(std::lround(px));
  const int cy = static_cast<int>(std::lround(py));
  for (int dy = -r; dy <= r; ++dy) {
    const int y = cy + dy;
    if (y < 0 || y >= H)
      continue;
    for (int dx = -r; dx <= r; ++dx) {
      const int x = ((cx + dx) % W + W) % W;
      if (!is_black(fx.get_pixel(x, y)))
        return true;
    }
  }
  return false;
}

// Count of non-black pixels across the whole canvas.
template <int W, int H>
inline size_t count_lit(const MeshFx &fx) {
  size_t lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++lit;
  return lit;
}

// ============================================================================
// Plot::Mesh::draw — wireframe: every unique edge actually drawn
// ============================================================================

// Wireframe draw lights every unique edge: extract the octahedron's 12 edges,
// draw, and assert each edge's projected geodesic midpoint is lit. Also checks
// the wireframe covers some but well under half the canvas.
inline void test_wireframe_draws_every_edge() {
  constexpr int W = 288, H = 144;
  configure_arenas_default(); // Plot::Mesh::draw samples edges via scratch_arena_a

  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  // Octahedron: 6 vertices on the ±axes, 12 edges. Edge midpoints are off the
  // poles (every edge joins two orthogonal axes), so each projects to a
  // well-defined, non-singular pixel.
  PolyMesh mesh = Solids::Platonic::octahedron(seed_a, seed_b);

  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(geom, 64);
  Plot::Mesh::extract_edges(mesh, edges);
  HS_EXPECT_EQ(edges.size(), (size_t)12); // octahedron has 12 edges

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe; // bare sink, no AA
    Plot::Mesh::draw<W, H>(pipe, c, mesh, white);
  }
  fx.advance_display();

  // Every unique edge's geodesic midpoint lands on a lit pixel.
  for (size_t e = 0; e < edges.size(); ++e) {
    Vector a = mesh.vertices[edges[e].u];
    Vector b = mesh.vertices[edges[e].v];
    Vector mid = ((a + b) * 0.5f).normalized();
    PixelCoords p = vector_to_pixel<W, H>(mid);
    HS_EXPECT_TRUE((lit_near<W, H>(fx, p.x, p.y, 3)));
  }

  // A wireframe lights some pixels but nowhere near the whole canvas.
  const size_t lit = count_lit<W, H>(fx);
  HS_EXPECT_GT(lit, (size_t)0);
  HS_EXPECT_LT(lit, (size_t)(W * H) / 2);
}

// ============================================================================
// Scan::Mesh::draw / SDF::Face — solid fill: interior lit, tiling has no holes
// ============================================================================

// Solid fill lights every face interior and tiles the whole sphere: each of the
// octahedron's 8 face centroids is lit, the fill covers far more than the
// wireframe, and a closed convex solid leaves essentially no holes.
inline void test_solid_fill_covers_faces_and_tiles_sphere() {
  constexpr int W = 288, H = 144;
  configure_arenas_default();

  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh poly = Solids::Platonic::octahedron(seed_a, seed_b);

  // Wireframe lit count, for the "fill covers far more than the wireframe"
  // contrast.
  MeshFx wire(W, H);
  {
    Canvas c(wire);
    Pipeline<W, H> pipe;
    Plot::Mesh::draw<W, H>(pipe, c, poly, white);
  }
  wire.advance_display();
  const size_t wire_lit = count_lit<W, H>(wire);

  // Compile to a MeshState (the solid scan path needs face_offsets) and fill.
  MeshState mesh;
  MeshOps::compile(poly, mesh, geom);

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe;
    Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch);
  }
  fx.advance_display();

  // (a) Each face's interior is lit: project the face centroid and assert a lit
  //     pixel there. A face whose bounding cull wrongly dropped its rows would
  //     leave its centroid dark.
  const uint8_t *fc = mesh.get_face_counts_data();
  const uint16_t *fi = mesh.get_faces_data();
  const uint16_t *fo = mesh.get_face_offsets_data();
  const size_t num_f = mesh.get_face_counts_size();
  HS_EXPECT_EQ(num_f, (size_t)8); // octahedron has 8 triangular faces
  for (size_t f = 0; f < num_f; ++f) {
    Vector centroid(0, 0, 0);
    for (int k = 0; k < fc[f]; ++k)
      centroid = centroid + mesh.vertices[fi[fo[f] + k]];
    centroid = centroid.normalized();
    PixelCoords p = vector_to_pixel<W, H>(centroid);
    HS_EXPECT_TRUE((lit_near<W, H>(fx, p.x, p.y, 2)));
  }

  // (b) The fill covers far more than the wireframe...
  const size_t fill_lit = count_lit<W, H>(fx);
  HS_EXPECT_GT(fill_lit, wire_lit * 4);

  // ...and a closed convex solid's spherical faces tile the entire sphere, so
  // the fill leaves essentially no holes. A conservative per-face row/column
  // cull is necessary for this: a cull that clipped any covered row would punch
  // a visible gap. Allow a thin slack for seam pixels exactly on a geodesic
  // boundary between two faces (a pixel center can miss both under the bare
  // non-AA fill rule).
  const size_t total = static_cast<size_t>(W) * H;
  HS_EXPECT_GT(fill_lit, total * 99 / 100);
}

// Runs every mesh-rasterization test in this module and returns the failure
// count from end_module.
inline int run_mesh_raster_tests() {
  auto scope = begin_module("mesh_raster");

  test_wireframe_draws_every_edge();
  test_solid_fill_covers_faces_and_tiles_sphere();

  return end_module(scope);
}

} // namespace mesh_raster_tests
} // namespace hs_test
