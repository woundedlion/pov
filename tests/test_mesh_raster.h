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
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_raster_tests {

/**
 * @brief Dedicated test arenas, kept alive for the whole test so the built
 *        mesh's ArenaVectors stay valid while it is drawn.
 * @details The global scratch_arena_a that Plot::Mesh::draw uses internally for
 *          per-edge line sampling is a different arena
 *          (configure_arenas_default()), so the source mesh is never clobbered.
 */
inline uint8_t mr_seed_a[256 * 1024];
inline uint8_t mr_seed_b[256 * 1024];
inline uint8_t mr_geom[256 * 1024];
inline uint8_t mr_scratch[256 * 1024];

/**
 * @brief Minimal Effect that owns a Canvas to draw into.
 * @details No background and no per-frame work, so the lit pixels come solely
 *          from the mesh draw under test.
 */
struct MeshFx : public Effect {
  /**
   * @brief Constructs the effect with a canvas of the given size.
   * @param W Canvas width in pixels.
   * @param H Canvas height in pixels.
   */
  MeshFx(int W, int H) : Effect(W, H) {}
  /**
   * @brief Per-frame draw hook; intentionally a no-op for this test fixture.
   */
  void draw_frame() override {}
  /**
   * @brief Reports whether the effect paints a background.
   * @return Always false, so only the mesh draw produces lit pixels.
   */
  bool strobe_columns() const override { return false; }
};

/**
 * @brief Tests whether a pixel is fully black.
 * @param p Pixel to inspect.
 * @return True when all of the pixel's RGB channels are zero.
 */
inline bool is_black(const Pixel &p) { return p.r == 0 && p.g == 0 && p.b == 0; }

/**
 * @brief Fragment shader that paints a fixed near-white.
 * @param f Fragment whose color is overwritten; used as the draw color so any
 *          covered pixel becomes non-black and thus detectable.
 * @details The unnamed Vector parameter (surface position) is ignored.
 */
inline void white(const Vector &, Fragment &f) {
  f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
}

/**
 * @brief Reports whether any lit pixel lies within radius r of (px, py).
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param fx Effect whose canvas is sampled.
 * @param px Target column in pixels.
 * @param py Target row in pixels.
 * @param r Search radius in pixels.
 * @return True if a non-black pixel is found in the neighborhood.
 * @details Wraps in x (the canvas is a cylinder in longitude) and clamps in y.
 *          Tolerance absorbs sampling granularity, AA spread, and projection
 *          rounding.
 */
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

/**
 * @brief Counts the non-black pixels across the whole canvas.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param fx Effect whose canvas is scanned.
 * @return Number of lit (non-black) pixels.
 */
template <int W, int H>
inline size_t count_lit(const MeshFx &fx) {
  size_t lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++lit;
  return lit;
}

/**
 * @brief Angular distance (radians) from a unit vector to a great-circle arc.
 * @param v Query unit vector.
 * @param a Arc start (unit vector).
 * @param b Arc end (unit vector).
 * @return The smaller of the off-plane distance (when v projects onto the arc
 *         interior) or the nearer endpoint distance (when it projects outside).
 * @details Used as the analytic oracle for "is this lit pixel ON some edge?":
 *          n = a×b is the arc's plane normal, asin(|v·n|) is the angle off the
 *          great circle, and a+b-span test places the projection inside the arc.
 */
inline float arc_distance(const Vector &v, const Vector &a, const Vector &b) {
  const auto ang = [](const Vector &p, const Vector &q) {
    return acosf(hs::clamp(dot(p, q), -1.0f, 1.0f));
  };
  Vector n = cross(a, b);
  float nlen = n.length();
  const float endpoints = std::min(ang(v, a), ang(v, b));
  if (nlen < 1e-6f)
    return endpoints; // degenerate edge (a ∥ b): nearest endpoint
  n = n * (1.0f / nlen);
  float off = asinf(hs::clamp(std::abs(dot(v, n)), 0.0f, 1.0f));
  Vector p = v - n * dot(v, n); // v projected onto the arc's plane
  float plen = p.length();
  if (plen > 1e-6f) {
    p = p * (1.0f / plen);
    if (ang(a, p) + ang(p, b) <= ang(a, b) + 1e-3f)
      return off; // projection lands on the arc interior
  }
  return endpoints;
}

// ============================================================================
// Plot::Mesh::draw — wireframe: every unique edge actually drawn
// ============================================================================

/**
 * @brief Verifies the wireframe draw lights every unique edge.
 * @details Extracts the octahedron's 12 edges, draws, and asserts each edge's
 *          projected geodesic midpoint is lit. Also checks the wireframe covers
 *          some but well under half the canvas.
 */
inline void test_wireframe_draws_every_edge() {
  constexpr int W = 288, H = 144;
  configure_arenas_default(); // Plot::Mesh::draw samples edges via scratch_arena_a

  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  // Octahedron edge midpoints are off the poles, so each projects to a
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

  for (size_t e = 0; e < edges.size(); ++e) {
    Vector a = mesh.vertices[edges[e].u];
    Vector b = mesh.vertices[edges[e].v];
    Vector mid = ((a + b) * 0.5f).normalized();
    PixelCoords p = vector_to_pixel<W, H>(mid);
    HS_EXPECT_TRUE((lit_near<W, H>(fx, p.x, p.y, 3)));
  }

  const size_t lit = count_lit<W, H>(fx);
  HS_EXPECT_GT(lit, (size_t)0);
  HS_EXPECT_LT(lit, (size_t)(W * H) / 2);
}

/**
 * @brief Geometric oracle: every lit wireframe pixel lies ON an edge arc.
 * @details The midpoint test above proves each edge IS drawn; this proves
 *          nothing is drawn OFF the edges. Each lit pixel is mapped back to its
 *          world direction and must fall within ~a few pixels (angular) of the
 *          nearest edge's analytic great-circle arc. A misprojected sample, a
 *          stray seam-wrap plot, or an edge routed to the wrong vertices would
 *          leave lit pixels off every arc and fail here — exactly the off-by-one
 *          / projection class the "nonempty subset" check (above) passes through.
 */
inline void test_wireframe_pixels_lie_on_edges() {
  constexpr int W = 288, H = 144;
  configure_arenas_default();

  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  PolyMesh mesh = Solids::Platonic::octahedron(seed_a, seed_b);
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(geom, 64);
  Plot::Mesh::extract_edges(mesh, edges);

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe; // bare sink, no AA spread
    Plot::Mesh::draw<W, H>(pipe, c, mesh, white);
  }
  fx.advance_display();

  // Tolerance: a few rows of latitude to absorb single-pixel line width,
  // sampling granularity, and vector_to_pixel rounding.
  const float tol = 4.0f * (PI_F / H);
  size_t lit = 0, off_arc = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      if (is_black(fx.get_pixel(x, y)))
        continue;
      ++lit;
      Vector v = pixel_to_vector<W, H>(x, y);
      float best = PI_F;
      for (size_t e = 0; e < edges.size(); ++e) {
        float d = arc_distance(v, mesh.vertices[edges[e].u],
                               mesh.vertices[edges[e].v]);
        if (d < best)
          best = d;
      }
      if (best > tol)
        ++off_arc;
    }
  HS_EXPECT_GT(lit, (size_t)0);
  HS_EXPECT_EQ(off_arc, (size_t)0);
}

// ============================================================================
// Scan::Mesh::draw / SDF::Face — solid fill: interior lit, tiling has no holes
// ============================================================================

/**
 * @brief Verifies the solid fill lights every face interior and tiles the whole
 *        sphere.
 * @details Each of the octahedron's 8 face centroids is lit, the fill covers
 *          far more than the wireframe, and a closed convex solid leaves
 *          essentially no holes.
 */
inline void test_solid_fill_covers_faces_and_tiles_sphere() {
  constexpr int W = 288, H = 144;
  configure_arenas_default();

  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh poly = Solids::Platonic::octahedron(seed_a, seed_b);

  // Wireframe lit count, for the fill-vs-wireframe contrast. `wire` and `fx`
  // alias the same static double buffer, so `wire` must be torn down (count
  // captured) before `fx` is built — only one Effect may be live at once.
  size_t wire_lit;
  {
    MeshFx wire(W, H);
    {
      Canvas c(wire);
      Pipeline<W, H> pipe;
      Plot::Mesh::draw<W, H>(pipe, c, poly, white);
    }
    wire.advance_display();
    wire_lit = count_lit<W, H>(wire);
  }

  // Compile to a MeshState; the solid scan path needs face_offsets.
  MeshState mesh;
  MeshOps::compile(poly, mesh, geom);

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe;
    Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch);
  }
  fx.advance_display();

  // (a) Each face interior is lit: project the centroid and assert a lit pixel
  //     there. A face whose bounding cull dropped its rows leaves it dark.
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

  // (b) The fill covers far more than the wireframe.
  const size_t fill_lit = count_lit<W, H>(fx);
  HS_EXPECT_GT(fill_lit, wire_lit * 4);

  // A closed convex solid tiles the whole sphere: every pixel center lands in
  // some face, so the fill covers the canvas with no holes. Any dark pixel is an
  // edge hole or a clipping artifact.
  const size_t total = static_cast<size_t>(W) * H;
  HS_EXPECT_EQ(fill_lit, total);
}

/**
 * @brief Runs every mesh-rasterization test in this module.
 * @return Failure count reported by end_module.
 */
inline int run_mesh_raster_tests() {
  hs_test::ModuleFixture fixture("mesh_raster");

  test_wireframe_draws_every_edge();
  test_wireframe_pixels_lie_on_edges();
  test_solid_fill_covers_faces_and_tiles_sphere();

  return fixture.result();
}

} // namespace mesh_raster_tests
} // namespace hs_test
