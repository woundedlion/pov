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
#include <vector>

#include "core/render/canvas.h"
#include "core/math/geometry.h"
#include "core/mesh/mesh.h"
#include "core/engine/memory.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/mesh/solids.h"
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

// ============================================================================
// Generic oracles over an arbitrary closed convex solid — exercise the
// wireframe edge-arc and solid-fill tiling paths beyond the octahedron, on
// quad faces (cube), pentagon faces (dodecahedron), and a large mixed-face
// Goldberg mesh (truncated icosahedron: 60 V, 90 E, 32 mixed pentagon/hexagon
// faces). Mixed-face face-walk, larger edge dedup, and the per-face bounding
// cull at higher face counts are pixel-unverified on the octahedron alone.
// ============================================================================

/**
 * @brief Asserts every lit wireframe pixel lies on some edge arc, for any mesh.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param mesh Closed mesh whose edges are extracted and drawn.
 * @param geom Arena backing the extracted edge list.
 * @details Mirrors test_wireframe_pixels_lie_on_edges' oracle: each lit pixel's
 *          world direction must fall within a few rows of latitude of the
 *          nearest analytic great-circle edge arc.
 */
template <int W, int H>
inline void check_wireframe_pixels_on_edges(PolyMesh &mesh, Arena &geom,
                                            size_t edge_cap) {
  configure_arenas_default();
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(geom, edge_cap);
  Plot::Mesh::extract_edges(mesh, edges);

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe;
    Plot::Mesh::draw<W, H>(pipe, c, mesh, white);
  }
  fx.advance_display();

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

/**
 * @brief Asserts a closed convex solid's solid fill lights every face centroid
 *        and tiles the whole sphere with no holes, for any mesh.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param poly Closed convex mesh to compile and fill.
 * @param geom Arena backing the compiled MeshState.
 * @param scratch Arena for the solid scan path.
 */
template <int W, int H>
inline void check_solid_fill_tiles(PolyMesh &poly, Arena &geom, Arena &scratch) {
  configure_arenas_default();
  MeshState mesh;
  MeshOps::compile(poly, mesh, geom);

  MeshFx fx(W, H);
  {
    Canvas c(fx);
    Pipeline<W, H> pipe;
    Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch);
  }
  fx.advance_display();

  const uint8_t *fc = mesh.get_face_counts_data();
  const uint16_t *fi = mesh.get_faces_data();
  const uint16_t *fo = mesh.get_face_offsets_data();
  const size_t num_f = mesh.get_face_counts_size();
  for (size_t f = 0; f < num_f; ++f) {
    Vector centroid(0, 0, 0);
    for (int k = 0; k < fc[f]; ++k)
      centroid = centroid + mesh.vertices[fi[fo[f] + k]];
    centroid = centroid.normalized();
    PixelCoords p = vector_to_pixel<W, H>(centroid);
    HS_EXPECT_TRUE((lit_near<W, H>(fx, p.x, p.y, 2)));
  }

  const size_t total = static_cast<size_t>(W) * H;
  const size_t fill_lit = count_lit<W, H>(fx);
  HS_EXPECT_EQ(fill_lit, total);
}

/**
 * @brief Wireframe + solid-fill oracles on the cube (quad faces).
 */
inline void test_cube_wireframe_and_fill() {
  constexpr int W = 288, H = 144;
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh cube = Solids::Platonic::cube(seed_a, seed_b);
  check_wireframe_pixels_on_edges<W, H>(cube, geom, 64);
  Arena geom2(mr_geom, sizeof(mr_geom)); // fresh: edge list above is dead now
  check_solid_fill_tiles<W, H>(cube, geom2, scratch);
}

/**
 * @brief Wireframe + solid-fill oracles on the dodecahedron (pentagon faces).
 */
inline void test_dodecahedron_wireframe_and_fill() {
  constexpr int W = 288, H = 144;
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh dodec = Solids::Platonic::dodecahedron(seed_a, seed_b);
  check_wireframe_pixels_on_edges<W, H>(dodec, geom, 64);
  Arena geom2(mr_geom, sizeof(mr_geom));
  check_solid_fill_tiles<W, H>(dodec, geom2, scratch);
}

/**
 * @brief Wireframe + solid-fill oracles on a large mixed-face Goldberg mesh.
 * @details The truncated icosahedron (60 V, 90 E, 32 mixed pentagon/hexagon
 *          faces) is the effect-payload scale the rasterizer-bound hot path
 *          actually runs, with non-triangular mixed faces the octahedron lacks.
 */
inline void test_truncated_icosahedron_wireframe_and_fill() {
  constexpr int W = 288, H = 144;
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh goldberg = Solids::Archimedean::truncatedIcosahedron(seed_a, seed_b);
  check_wireframe_pixels_on_edges<W, H>(goldberg, geom, 128); // 90 edges
  Arena seed_a2(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b2(mr_seed_b, sizeof(mr_seed_b));
  Arena geom2(mr_geom, sizeof(mr_geom));
  // Rebuild: the SolidBuilder above consumed seed_a/seed_b as its op pipeline.
  PolyMesh goldberg2 =
      Solids::Archimedean::truncatedIcosahedron(seed_a2, seed_b2);
  check_solid_fill_tiles<W, H>(goldberg2, geom2, scratch);
}

/**
 * @brief Verifies the per-face clip cull never drops an in-band pixel.
 * @details Renders a sphere-tiling mixed-face solid full-canvas, then under
 * several clip bands, and asserts the clipped output equals the full output at
 * every pixel inside each band. The solid tiles the sphere, so faces straddle
 * every row/column boundary including the poles (rows 0 / H-1) and the x=0 seam
 * — exactly the great-circle-arc cases an unsafe vertex-based cull would drop.
 */
inline void test_clip_band_matches_full() {
  constexpr int W = 288, H = 144;
  configure_arenas_default();
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  PolyMesh poly = Solids::Archimedean::truncatedIcosahedron(seed_a, seed_b);
  MeshState mesh;
  MeshOps::compile(poly, mesh, geom);

  std::vector<Pixel> ref(static_cast<size_t>(W) * H);
  {
    MeshFx full(W, H);
    { Canvas c(full); Pipeline<W, H> pipe;
      Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch); }
    full.advance_display();
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        ref[static_cast<size_t>(y) * W + x] = full.get_pixel(x, y);
  }

  struct Band { int y0, y1, x0, x1; };
  const Band bands[] = {
      {0, 72, 0, 144},     // top-left quadrant (north pole + x=0 seam)
      {72, 144, 144, 288}, // bottom-right quadrant (south pole + x=144)
      {0, 144, 0, 144},    // left half (x=144 boundary, full height)
      {0, 72, 0, 288},     // top band, full width (pole, no x cull)
  };
  for (const Band &b : bands) {
    MeshFx fx(W, H); // one live Effect at a time; ref captured above
    fx.set_clip(b.y0, b.y1, b.x0, b.x1);
    { Canvas c(fx); Pipeline<W, H> pipe;
      Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch); }
    fx.advance_display();
    size_t mismatches = 0;
    for (int y = b.y0; y < b.y1; ++y)
      for (int x = b.x0; x < b.x1; ++x) {
        const Pixel &p = fx.get_pixel(x, y);
        const Pixel &r = ref[static_cast<size_t>(y) * W + x];
        if (p.r != r.r || p.g != r.g || p.b != r.b)
          ++mismatches;
      }
    HS_EXPECT_EQ(mismatches, (size_t)0);
  }
}

// ============================================================================
// Congruence-class bake — census invariants + rendered A/B
//
// build_mesh_class_bake clusters a spawned mesh's faces into congruence
// classes (geometric clustering seeded per topology class) and bakes one
// canonical distance LUT per concave shared class. The census over the whole
// registry established: every islamic mesh's faces land 100% in shared
// classes, with <= 24 classes and worst Procrustes residual < 0.25 px at
// W=288. These tests pin those invariants on registry meshes and verify the
// LUT-served render matches the exact render.
// ============================================================================

/**
 * @brief Generates, compiles, and classifies one islamic registry mesh, then
 *        builds its congruence-class bake.
 * @param islamic_idx Index into Solids::Collections::get_islamic_solids().
 * @param seed_a First generation scratch arena.
 * @param seed_b Second generation scratch arena.
 * @param geom Arena receiving the compiled mesh and the bake.
 * @param mesh Output compiled mesh.
 * @param bake Output congruence-class bake.
 */
inline void build_islamic_bake(size_t islamic_idx, Arena &seed_a, Arena &seed_b,
                               Arena &geom, MeshState &mesh,
                               MeshOps::MeshClassBake &bake) {
  constexpr int W = 288;
  const auto islamic = Solids::Collections::get_islamic_solids();
  HS_EXPECT_TRUE(islamic_idx < islamic.size());
  PolyMesh poly = islamic[islamic_idx].generate(seed_a, seed_b);
  MeshOps::compile(poly, mesh, geom);
  MeshOps::classify_faces_by_topology(mesh, seed_a, seed_b, geom);
  MeshOps::build_mesh_class_bake(mesh, seed_a, geom, 2.0f * PI_F / W, bake);
}

/**
 * @brief Verifies the census invariants of the congruence clustering on one
 *        registry mesh.
 * @param islamic_idx Index into the islamic registry.
 * @details 100% of faces in shared classes, class count within the cap, worst
 * accepted residual under the congruence epsilon, and at least one LUT built.
 */
inline void check_class_bake_census(size_t islamic_idx) {
  configure_arenas_default();
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  MeshState mesh;
  MeshOps::MeshClassBake bake;
  build_islamic_bake(islamic_idx, seed_a, seed_b, geom, mesh, bake);

  const size_t F = mesh.num_faces();
  HS_EXPECT_EQ(bake.face_recs.size(), F);
  HS_EXPECT_GT(bake.classes.size(), (size_t)0);
  HS_EXPECT_LE(bake.classes.size(), (size_t)MeshOps::MAX_CONGRUENCE_CLASSES);

  // Every face lands in a class, and every class is shared (>= 2 members).
  size_t assigned = 0, member_sum = 0;
  for (size_t f = 0; f < F; ++f)
    if (bake.face_recs[f].class_id != MeshOps::NO_CLASS)
      ++assigned;
  for (size_t c = 0; c < bake.classes.size(); ++c) {
    HS_EXPECT_GT(bake.classes[c].members, (uint16_t)1);
    member_sum += bake.classes[c].members;
  }
  HS_EXPECT_EQ(assigned, F);
  HS_EXPECT_EQ(member_sum, F);
  HS_EXPECT_EQ(bake.shared_faces, (uint16_t)F);

  HS_EXPECT_LT(bake.worst_residual_px, MeshOps::CONGRUENCE_EPS_PX);
  HS_EXPECT_GT(bake.luts_built, (uint16_t)0);
}

/**
 * @brief Census invariants on two registry meshes.
 */
inline void test_class_bake_census_invariants() {
  check_class_bake_census(0);
  check_class_bake_census(2);
}

/**
 * @brief Fragment shader encoding the signed distance (v1) as brightness, so
 *        interior-gradient deviations become channel deltas.
 * @param f Fragment whose color is overwritten from its v1 register.
 */
inline void shade_by_distance(const Vector &, Fragment &f) {
  float g = hs::clamp(30000.0f + f.v1 * 200000.0f, 0.0f, 60000.0f);
  uint16_t q = static_cast<uint16_t>(g);
  f.color = Color4(Pixel(q, q, q), 1.0f);
}

/**
 * @brief Renders a mesh twice — exact path vs class bake — and asserts the
 *        outputs agree within the interpolation envelope.
 * @param mesh Mesh to draw (possibly deformed after its bake was built).
 * @param bake Spawn-time congruence bake for the mesh.
 * @param min_lut_hits Floor asserting the LUT path actually served probes.
 * @param label Telemetry tag for the printf lines.
 * @details The shader encodes raw distance, so a wrong reflection/offset
 * convention, a mis-rotated alignment, or a sign flip near a deformed edge
 * (face-separation cracks) blows the mean; the per-pixel cap bounds the
 * legitimate bilinear + congruence + deformation-margin deviation.
 * Frame-budget-tight deltas are the offline visual gate's job (400-frame
 * dump), not this smoke's.
 */
inline void check_class_lut_render_matches_exact(const MeshState &mesh,
                                                 const MeshOps::MeshClassBake &bake,
                                                 uint32_t min_lut_hits,
                                                 const char *label) {
  constexpr int W = 288, H = 144;
  Arena scratch(mr_scratch, sizeof(mr_scratch));

  std::vector<Pixel> ref(static_cast<size_t>(W) * H);
  {
    MeshFx exact(W, H);
    { Canvas c(exact); Pipeline<W, H> pipe;
      Scan::Mesh::draw<W, H>(pipe, c, mesh, shade_by_distance, scratch); }
    exact.advance_display();
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        ref[static_cast<size_t>(y) * W + x] = exact.get_pixel(x, y);
  }

  hs::g_scan_metrics.reset();
  MeshFx lutted(W, H);
  { Canvas c(lutted); Pipeline<W, H> pipe;
    Scan::Mesh::draw<W, H>(pipe, c, mesh, shade_by_distance, scratch,
                           /*debug_bb=*/false, &bake); }
  lutted.advance_display();

  // The LUT path actually served a meaningful share of the probes.
  const uint32_t lut_hits = hs::g_scan_metrics.lut_hits;
  const uint32_t exact_hits = hs::g_scan_metrics.exact_hits;
  HS_EXPECT_GT(lut_hits, min_lut_hits);
  std::printf("  [%s] hits=%u exact=%u share=%.1f%%\n", label, lut_hits,
              exact_hits,
              100.0f * lut_hits / std::max(1u, lut_hits + exact_hits));

  uint64_t delta_sum = 0;
  int delta_max = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = lutted.get_pixel(x, y);
      const Pixel &r = ref[static_cast<size_t>(y) * W + x];
      int d = std::abs((int)p.r - (int)r.r);
      delta_sum += d;
      if (d > delta_max)
        delta_max = d;
    }
  const float mean = static_cast<float>(delta_sum) / (W * H);
  std::printf("  [%s] delta mean=%.1f max=%d (of 60000)\n", label, mean,
              delta_max);
  HS_EXPECT_LT(mean, 300.0f);    // ~0.5% FS: catches convention bugs
  HS_EXPECT_LT(delta_max, 6000); // ~10% FS: interpolation envelope; rippled
                                 // must match it too (bent faces go exact)
}

/**
 * @brief Rendered A/B on the undeformed mesh.
 */
inline void test_class_lut_render_matches_exact() {
  configure_arenas_default();
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  MeshState mesh;
  MeshOps::MeshClassBake bake;
  build_islamic_bake(0, seed_a, seed_b, geom, mesh, bake);
  check_class_lut_render_matches_exact(mesh, bake, 5000, "class lut");
}

/**
 * @brief Rendered A/B on a deformed mesh against its spawn-time bake — the
 *        facility's deformation safety net.
 * @details The class-LUT facility is for meshes that hold still (see
 * mesh_classes.h); this pins what happens when a consumer's mesh deforms
 * anyway: the per-frame alignment must drop bent faces to the exact path (or
 * widen its guard) so the output stays inside the static interpolation
 * envelope. A fixed one-cell guard flips signs near the true edges and opens
 * visible cracks between faces.
 */
inline void test_class_lut_render_matches_exact_rippled() {
  configure_arenas_default();
  Arena seed_a(mr_seed_a, sizeof(mr_seed_a));
  Arena seed_b(mr_seed_b, sizeof(mr_seed_b));
  Arena geom(mr_geom, sizeof(mr_geom));

  MeshState mesh;
  MeshOps::MeshClassBake bake;
  build_islamic_bake(0, seed_a, seed_b, geom, mesh, bake);

  // The real ripple transform at its shipped ceiling (amplitude 0.15,
  // thickness 0.7): a Ricker wavelet that slides vertices tangentially away
  // from the origin — the steep wavelet slope shears faces, which is what
  // breaks a rigid canonical alignment. (Radial displacement would be a
  // false-pass: the gnomonic projection divides it out.)
  Animation::RippleParams rp;
  rp.center = Vector(0.3f, 0.8f, -0.52f).normalized();
  rp.amplitude = 0.15f;
  rp.thickness = 0.7f;
  rp.decay = 0.1f;
  rp.phase = 0.9f; // mid-expansion: wavefront crossing plenty of faces
  rp.prepare_thresholds();
  for (size_t i = 0; i < mesh.vertices.size(); ++i)
    mesh.vertices[i] = ripple_transform(mesh.vertices[i], rp);
  check_class_lut_render_matches_exact(mesh, bake, 1000, "class lut rippled");
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
  test_cube_wireframe_and_fill();
  test_dodecahedron_wireframe_and_fill();
  test_truncated_icosahedron_wireframe_and_fill();
  test_clip_band_matches_full();
  test_class_bake_census_invariants();
  test_class_lut_render_matches_exact();
  test_class_lut_render_matches_exact_rippled();

  return fixture.result();
}

} // namespace mesh_raster_tests
} // namespace hs_test
