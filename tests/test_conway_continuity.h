/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Boundary color-continuity tests for the ConwayMorph transition design
 * (docs/conway_morph_spec.md §7.6): every mesh swap at a leg or cycle boundary
 * must exchange geometrically matching meshes whose positive-area faces keep
 * their colors.
 *
 * Coverage:
 *   - Bookend angle pin: driving HankinSolids to a hankin-cycle end forces the
 *     interlace angle of the final drawn frame to exactly 0 (the flat p_corner
 *     branch), where the sweep's own last sample lands ~0.002 rad off flat.
 *   - Bookend swaps (per node, one solid per symmetry family): update_hankin
 *     at angle 0 emits first-F star faces whose boundary vertices lie on the
 *     base face's boundary; an in-memory framebuffer diff of the two renders
 *     under identity face colors changes no pixels, and the zero-area rosette
 *     faces draw none.
 *   - Forward palette carry (per arrival, real effect): the palettes the leg
 *     landed carry verbatim into the new node's displayed base faces.
 *   - Leg swaps (per edge): base vs op(seed, T_EPS) and the reseed swaps
 *     (ADOPT bridge arrival, DUAL_SWAP ambo crossover) framebuffer-diff within
 *     a budget far below one face's area, so a face landing under the wrong
 *     palette mapping fails loudly (pinned by a deliberate wrong-mapping run).
 *   - Palette mapping/crossfade units: every mapping is total and
 *     deterministic; the crossfade is exact at its endpoints — the first sweep
 *     frame (w = 0) shades every surviving face from its inherited palette,
 *     the last (w = 1) from the leg's landed target assignment.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "core/animation/mesh.h"
#include "core/color/palettes.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "core/render/canvas.h"
#include "core/render/scan.h"
#include "effects/HankinSolids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_conway_morph.h" // run_edge_op
#include "tests/test_conway_soak.h"  // HankinWalkProbe
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_continuity_tests {

using ConwayGraph::T_EPS;

inline uint8_t cc_geom_buf[256 * 1024]; /**< Mesh / compiled-hankin arena. */
inline uint8_t cc_temp_buf[256 * 1024]; /**< Op scratch arena. */
inline uint8_t cc_aux_buf[256 * 1024];  /**< Second-mesh arena. */
inline uint8_t cc_scan_buf[256 * 1024]; /**< Rasterizer face-scratch arena. */
inline uint8_t cc_leg_buf[256 * 1024];  /**< ConwayMorph leg arena. */
inline uint8_t cc_bank_buf[64 * 1024];  /**< Baked palette LUT arena. */

// ---------------------------------------------------------------------------
// Framebuffer harness: render one compiled mesh with flat per-face colors and
// capture the frame in memory.
// ---------------------------------------------------------------------------

/** Render size for the framebuffer diffs: the production resolution, where a
 * T_EPS corner cut spans ~1-2 px and the smallest compared face spans
 * thousands, so the diff budgets discriminate. */
constexpr int FB_W = 288;
constexpr int FB_H = 144;

/**
 * @brief Minimal Effect that owns a Canvas to draw into.
 */
struct ContFx : public Effect {
  /**
   * @brief Constructs the effect with a FB_W x FB_H canvas.
   */
  ContFx() : Effect(FB_W, FB_H) {}
  /**
   * @brief Per-frame draw hook; unused, drawing happens via a local Canvas.
   */
  void draw_frame() override {}
};

/**
 * @brief Well-separated flat face color for an emission index.
 * @param i Face emission index (or any small id).
 * @return Fully opaque color, unique per index for the face counts under test.
 */
inline Color4 face_color(int i) {
  const auto ch = [i](int mult) {
    return static_cast<uint16_t>(4111 + ((i + 1) * mult) % 57000);
  };
  return Color4(Pixel(ch(9973), ch(15731), ch(23459)), 1.0f);
}

/**
 * @brief Renders a compiled mesh with a per-face flat color and captures the
 *        framebuffer.
 * @param out Receives the FB_W x FB_H frame, row-major.
 * @param mesh Compiled mesh (needs face_offsets).
 * @param color_of Maps the face emission index to its flat fill color.
 */
template <typename ColorFn>
inline void render_faces(std::vector<Pixel> &out, const MeshState &mesh,
                         ColorFn &&color_of) {
  ContFx fx;
  {
    Canvas c(fx);
    Pipeline<FB_W, FB_H> pipe;
    Arena scan_scratch(cc_scan_buf, sizeof(cc_scan_buf));
    Scan::Mesh::draw<FB_W, FB_H>(
        pipe, c, mesh,
        [&](const Vector &, Fragment &f) {
          f.color = color_of(static_cast<int>(f.v2));
        },
        scan_scratch);
  }
  fx.advance_display();
  out.resize(static_cast<size_t>(FB_W) * FB_H);
  for (int y = 0; y < FB_H; ++y)
    for (int x = 0; x < FB_W; ++x)
      out[static_cast<size_t>(y) * FB_W + x] = fx.get_pixel(x, y);
}

/**
 * @brief Whether a pixel exactly carries a flat fill color.
 */
inline bool color_eq(const Pixel &p, const Color4 &c) {
  return p.r == c.color.r && p.g == c.color.g && p.b == c.color.b;
}

/**
 * @brief Whether a pixel is "flat": black background or an exact member of the
 *        render's fill palette (i.e. a face interior, not a boundary AA blend).
 */
inline bool is_flat(const Pixel &p, const std::vector<Color4> &palette) {
  if (p.r == 0 && p.g == 0 && p.b == 0)
    return true;
  for (const Color4 &c : palette)
    if (color_eq(p, c))
      return true;
  return false;
}

/**
 * @brief Bucketed framebuffer diff.
 * @details The rasterizer anti-aliases face boundaries, and the two compared
 *          meshes tessellate the same outlines with different vertex lists, so
 *          boundary pixels blend at slightly different ratios; those band
 *          pixels say nothing about face-color continuity and are counted
 *          apart. The continuity signal is `flat`: pixels solidly inside a
 *          face on both sides that changed color.
 */
struct DiffStats {
  size_t newborn = 0; /**< Differing pixels the second render fills with the
                           newborn sentinel (legitimate zero-area births). */
  size_t flat = 0;    /**< Differing pixels flat on both sides — a face
                           interior changed color. */
  size_t blend = 0;   /**< Differing pixels involving a boundary AA blend. */
};

/**
 * @brief Buckets the differing pixels of two captures.
 * @param a First frame and its flat fill palette.
 * @param b Second frame and its flat fill palette.
 * @param sentinel Newborn fill color of the second render, or nullptr.
 * @return Bucketed counts.
 */
inline DiffStats diff_stats(const std::vector<Pixel> &a,
                            const std::vector<Color4> &flat_a,
                            const std::vector<Pixel> &b,
                            const std::vector<Color4> &flat_b,
                            const Color4 *sentinel = nullptr) {
  HS_EXPECT_EQ(a.size(), b.size());
  DiffStats st;
  for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
    if (a[i].r == b[i].r && a[i].g == b[i].g && a[i].b == b[i].b)
      continue;
    if (sentinel && color_eq(b[i], *sentinel))
      ++st.newborn;
    else if (is_flat(a[i], flat_a) && is_flat(b[i], flat_b))
      ++st.flat;
    else
      ++st.blend;
  }
  return st;
}

/**
 * @brief Flat fill palette of a face-indexed render.
 * @param n Face count.
 * @return face_color(0 .. n-1).
 */
inline std::vector<Color4> face_palette(int n) {
  std::vector<Color4> p;
  for (int i = 0; i < n; ++i)
    p.push_back(face_color(i));
  return p;
}

// ---------------------------------------------------------------------------
// §7.6a remainder — bookend angle pin: the cycle's final drawn frame is forced
// to exactly angle 0, not the sweep's ~0.002 rad boundary sample.
// ---------------------------------------------------------------------------

/**
 * @brief Drives HankinSolids across a full hankin cycle and asserts the
 *        closing bookend frame's interlace angle is exactly 0.
 * @details Observes the registered "Angle" param after every frame. The pin
 *          frame is the first exact-0 frame whose predecessor is a small
 *          non-zero sweep sample; the following morph leg leaves the angle
 *          untouched, so exact 0 must then hold for the whole leg — which also
 *          rules out a mid-sweep zero crossing matching the predicate.
 */
inline void test_bookend_angle_pin() {
  reset_globals();
  HankinSolids<96, 20> fx;
  fx.init();

  const auto *angle = fx.getParameters().find("Angle");
  HS_EXPECT_TRUE(angle != nullptr);
  if (!angle)
    return;

  std::vector<float> a;
  for (int f = 0; f < 200; ++f) {
    fx.draw_frame();
    fx.advance_display();
    a.push_back(angle->get());
  }

  int pin = -1;
  for (size_t i = 1; i + 40 < a.size(); ++i) {
    if (a[i] != 0.0f || a[i - 1] == 0.0f)
      continue;
    // Cycle end: the sweep's final sample sits just off flat (not the initial
    // default or a mid-sweep write), and the morph leg holds exact 0 after.
    if (a[i - 1] < 0.05f) {
      pin = static_cast<int>(i);
      break;
    }
  }
  HS_EXPECT_GT(pin, 0);
  if (pin <= 0)
    return;
  HS_EXPECT_GT(a[pin - 1], 0.0f);
  HS_EXPECT_LT(a[pin - 1], 0.05f);
  for (int i = pin; i < pin + 40; ++i)
    HS_EXPECT_EQ(a[i], 0.0f);
}

// ---------------------------------------------------------------------------
// §7.6b — bookend swaps (hankin at angle 0 <-> node base mesh), one node per
// symmetry family.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the angle-0 hankin mesh's first-F faces are geometrically
 *        the base faces: every star-face vertex lies on the base face's
 *        boundary vertex set (normalized corners and edge midpoints), and
 *        every corner is hit.
 * @param base Node base mesh.
 * @param flat update_hankin output at angle 0 for the same compiled pattern.
 */
inline void check_flat_star_faces_match_base(const PolyMesh &base,
                                             const PolyMesh &flat) {
  const size_t F = base.face_counts.size();
  HS_EXPECT_GE(flat.face_counts.size(), F);
  constexpr float TOL = 1e-4f;

  size_t base_off = 0;
  size_t flat_off = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    const int bc = base.face_counts[fi];
    const int oc = flat.face_counts[fi];
    HS_EXPECT_EQ(oc, 2 * bc);

    // Candidate boundary points of base face fi on the unit sphere.
    std::vector<Vector> corners, mids;
    for (int k = 0; k < bc; ++k) {
      const Vector c0 = base.vertices[base.faces[base_off + k]];
      const Vector c1 = base.vertices[base.faces[base_off + (k + 1) % bc]];
      corners.push_back(c0.normalized());
      mids.push_back(((c0 + c1) * 0.5f).normalized());
    }

    int corner_hits = 0;
    for (int k = 0; k < oc; ++k) {
      const Vector v = flat.vertices[flat.faces[flat_off + k]];
      bool on_corner = false, on_mid = false;
      for (int j = 0; j < bc; ++j) {
        if ((v - corners[j]).length() <= TOL)
          on_corner = true;
        if ((v - mids[j]).length() <= TOL)
          on_mid = true;
      }
      HS_EXPECT_TRUE(on_corner || on_mid);
      if (on_corner)
        ++corner_hits;
    }
    HS_EXPECT_EQ(corner_hits, bc);

    base_off += bc;
    flat_off += oc;
  }
}

// The spec contract (sections 2.4/2.5) is exact: the angle-0 hankin mesh and
// the base mesh draw the same shapes in the same colors, and the rosette
// faces are zero-area births that draw nothing — so the bookend swap changes
// no pixels at all.

/**
 * @brief Bookend swap for one solid: mesh-level geometric identity plus the
 *        framebuffer diff under identity face colors.
 * @tparam Solid Seed solid descriptor.
 * @details The hankin render colors star face f like base face f (the §2.5
 *          identity mapping) and paints rosette faces a loud sentinel; the
 *          swap must change no pixels, and the zero-area rosettes must draw
 *          none.
 */
template <typename Solid> inline void check_bookend_swap_one() {
  Arena geom(cc_geom_buf, sizeof(cc_geom_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));
  Arena aux(cc_aux_buf, sizeof(cc_aux_buf));

  PolyMesh base;
  build_solid<Solid>(base, geom);
  const int F = static_cast<int>(base.face_counts.size());

  CompiledHankin compiled;
  MeshOps::compile_hankin(base, compiled, geom, temp);
  PolyMesh flat;
  MeshOps::update_hankin(compiled, flat, geom, 0.0f);
  check_flat_star_faces_match_base(base, flat);

  MeshState base_ms;
  MeshOps::compile(base, base_ms, aux, temp);
  MeshState flat_ms;
  MeshOps::update_hankin(compiled, flat_ms, geom, 0.0f);

  const Color4 sentinel(Pixel(65535, 0, 65535), 1.0f);
  std::vector<Pixel> base_px, flat_px;
  render_faces(base_px, base_ms, [](int f) { return face_color(f); });
  render_faces(flat_px, flat_ms,
               [&](int f) { return f < F ? face_color(f) : sentinel; });

  const DiffStats st =
      diff_stats(base_px, face_palette(F), flat_px, face_palette(F), &sentinel);
  std::printf("  [bookend] F=%d: base vs hankin@0 flat=%zu newborn=%zu "
              "blend=%zu px\n",
              F, st.flat, st.newborn, st.blend);
  HS_EXPECT_EQ(st.flat + st.newborn + st.blend, (size_t)0);

  // Rosette isolation: stars painted background-black, rosettes white — any
  // lit pixel is rosette output.
  std::vector<Pixel> rosette_px;
  render_faces(rosette_px, flat_ms, [&](int f) {
    return f < F ? Color4(Pixel(0, 0, 0), 1.0f)
                 : Color4(Pixel(65535, 65535, 65535), 1.0f);
  });
  size_t rosette_lit = 0;
  for (const Pixel &p : rosette_px)
    if (p.r || p.g || p.b)
      ++rosette_lit;
  std::printf("  [bookend] F=%d: rosette contribution %zu px\n", F,
              rosette_lit);
  HS_EXPECT_EQ(rosette_lit, (size_t)0);
}

/**
 * @brief Bookend swap identity for one node per symmetry family.
 */
inline void test_bookend_swaps_per_family() {
  check_bookend_swap_one<Solids::Cube>();         // octahedral
  check_bookend_swap_one<Solids::Dodecahedron>(); // icosahedral
  check_bookend_swap_one<Solids::Tetrahedron>();  // tetrahedral
}

/**
 * @brief Verifies the §2.5/§2.6 forward palette mapping across real leg
 *        arrivals: every base face of the arrived node displays the palette
 *        its corresponding swept face landed with — carried verbatim.
 * @details Drives HankinSolids through several legs; the in-flight leg's
 *          Landing is snapshotted each frame, and on each arrival the
 *          displayed per-face palettes (node_face_palette_) are compared with
 *          the landed ones over the base-face emission prefix, as per-palette
 *          multisets — invariant under any correct provenance mapping, but
 *          broken by a class merge collapsing distinct landed palettes into
 *          one slot.
 */
inline void test_palette_carry_across_arrivals() {
  reset_globals();
  using Probe = conway_soak_tests::HankinWalkProbe;
  constexpr int PALETTES = Animation::ConwayMorph::PALETTES;

  HankinSolids<96, 20> fx;
  fx.init();

  std::array<uint8_t, PALETTES> to_palette{};
  std::vector<int> topo;
  bool have_landing = false;

  int prev_node = Probe::node(fx);
  int arrivals = 0;
  constexpr int TARGET_ARRIVALS = 10;
  for (int frame = 0; frame < 2200 && arrivals < TARGET_ARRIVALS; ++frame) {
    fx.draw_frame();
    fx.advance_display();

    if (const Animation::ConwayMorph::Landing *landing =
            Probe::pending_landing(fx)) {
      to_palette = landing->to_palette;
      topo.assign(landing->topology, landing->topology + landing->faces);
      have_landing = true;
    }

    const int node = Probe::node(fx);
    if (node == prev_node)
      continue;
    prev_node = node;
    if (!have_landing)
      continue;
    have_landing = false;
    ++arrivals;

    // Base faces are the landing's emission prefix (all faces at a full
    // arrival; the primaries at a t = 0 arrival, whose newborn faces die).
    const size_t nf = Probe::node_faces(fx);
    HS_EXPECT_LE(nf, topo.size());
    if (nf > topo.size())
      continue;
    int landed[PALETTES] = {};
    int shown[PALETTES] = {};
    for (size_t f = 0; f < nf; ++f) {
      ++landed[to_palette[wrap(topo[f], PALETTES)]];
      ++shown[Probe::node_face_palette(fx)[f]];
    }
    const int failed_before = hs_test::stats().failed;
    for (int p = 0; p < PALETTES; ++p)
      HS_EXPECT_EQ(shown[p], landed[p]);
    if (hs_test::stats().failed != failed_before)
      std::printf("    [palette-carry] arrival %d at '%s': displayed palettes "
                  "do not carry the landing\n",
                  arrivals, Solids::simple_registry[node].name);
  }
  HS_EXPECT_EQ(arrivals, TARGET_ARRIVALS);
}

// ---------------------------------------------------------------------------
// §7.6c — leg swaps: swap-in (base <-> op(seed, T_EPS)), the ADOPT bridge
// arrival, and the DUAL_SWAP ambo crossover.
// ---------------------------------------------------------------------------

/** Interior-pixel budget for the leg swaps: a face landing under the wrong
 * mapping recolors a whole interior (thousands of pixels at FB_H = 144). */
constexpr size_t SWAP_FLAT_BUDGET = 200;

/** Newborn-pixel budget: faces born at T_EPS cover at most corner cuts and
 * hairline edge bands. */
constexpr size_t SWAP_NEWBORN_BUDGET = 2500;

/** Backstop on boundary-blend jitter around the swept mesh's moved edges. */
constexpr size_t SWAP_BLEND_BUDGET = 12000;

/** Interior floor a wrong emission mapping must exceed (whole faces flip). */
constexpr size_t WRONG_MAPPING_FLAT_FLOOR = 5000;

/** Interior budget for the DUAL_SWAP crossover: both compared meshes sit
 * T_EPS off the clean cuboctahedron in opposite directions, so every
 * square-triangle boundary carries a ~1 px flat strip of real geometry
 * mismatch (measured 1352 px); a class-mapping error flips whole interiors. */
constexpr size_t DUAL_SWAP_FLAT_BUDGET = 4000;

/**
 * @brief Renders a seed and an operator output under the emission-identity
 *        mapping and buckets the framebuffer diff.
 * @param seed Seed mesh (drawn as compiled).
 * @param swept Operator output at a boundary parameter; primaries take the
 *        seed's colors, newborn faces the sentinel.
 * @param rotate_primaries When true, primaries take face_color((f+1) % F) —
 *        the deliberate wrong mapping used to pin the budgets' discrimination.
 * @return Bucketed diff counts.
 */
inline DiffStats swap_diff(const PolyMesh &seed, const PolyMesh &swept,
                           bool rotate_primaries = false) {
  Arena aux(cc_aux_buf, sizeof(cc_aux_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));

  const int F = static_cast<int>(seed.face_counts.size());
  MeshState seed_ms;
  MeshOps::compile(seed, seed_ms, aux, temp);
  std::vector<Pixel> seed_px;
  render_faces(seed_px, seed_ms, [](int f) { return face_color(f); });

  const Color4 sentinel(Pixel(65535, 0, 65535), 1.0f);
  MeshState swept_ms;
  MeshOps::compile(swept, swept_ms, aux, temp);
  std::vector<Pixel> swept_px;
  render_faces(swept_px, swept_ms, [&](int f) {
    if (f >= F)
      return sentinel;
    return face_color(rotate_primaries ? (f + 1) % F : f);
  });

  return diff_stats(seed_px, face_palette(F), swept_px, face_palette(F),
                    &sentinel);
}

/**
 * @brief Verifies the swap-in exchange (seed base mesh vs op at T_EPS) keeps
 *        every surviving face interior's color for each swept operator, and
 *        that the budget discriminates by failing a rotated (wrong) emission
 *        mapping.
 * @details Newborn vertex/edge faces render as a sentinel and may only cover
 *          their T_EPS-sized birth areas.
 */
inline void test_swap_in_framebuffer() {
  using conway_morph_tests::run_edge_op;
  const ConwayGraph::MorphOp OPS[] = {ConwayGraph::MorphOp::TRUNCATE,
                                      ConwayGraph::MorphOp::EXPAND,
                                      ConwayGraph::MorphOp::SNUB};
  const char *names[] = {"truncate", "expand", "snub"};

  for (int oi = 0; oi < 3; ++oi) {
    Arena geom(cc_geom_buf, sizeof(cc_geom_buf));
    Arena temp(cc_temp_buf, sizeof(cc_temp_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, geom);

    ConwayGraph::EdgeSpec e{};
    e.op = OPS[oi];
    PolyMesh swept = run_edge_op(e, cube, geom, temp, T_EPS, 0.0f);

    const DiffStats st = swap_diff(cube, swept);
    std::printf("  [swap-in] %s: flat=%zu newborn=%zu blend=%zu px\n",
                names[oi], st.flat, st.newborn, st.blend);
    HS_EXPECT_LE(st.flat, SWAP_FLAT_BUDGET);
    HS_EXPECT_LE(st.newborn, SWAP_NEWBORN_BUDGET);
    HS_EXPECT_LE(st.blend, SWAP_BLEND_BUDGET);

    if (oi == 0) {
      // Budget discrimination: rotating the primary mapping by one face must
      // recolor whole face interiors and blow far past the flat budget.
      const DiffStats wrong = swap_diff(cube, swept, /*rotate_primaries*/ true);
      std::printf("  [swap-in] wrong mapping: flat=%zu px (floor %zu)\n",
                  wrong.flat, WRONG_MAPPING_FLAT_FLOOR);
      HS_EXPECT_GT(wrong.flat, WRONG_MAPPING_FLAT_FLOOR);
    }
  }
}

/**
 * @brief Verifies the DERIVE_AMBO (ADOPT-row) departure swap: the displayed
 *        cuboctahedron base vs truncate(ambo(cube), T_EPS) — the first frame
 *        of the cuboctahedron -> truncatedCuboctahedron leg — under the
 *        emission-identity mapping on the derived seed.
 */
inline void test_adopt_departure_swap_framebuffer() {
  Arena geom(cc_geom_buf, sizeof(cc_geom_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh cubocta = MeshOps::ambo(cube, geom, temp);
  PolyMesh swept = MeshOps::truncate(cubocta, geom, temp, T_EPS);

  const DiffStats st = swap_diff(cubocta, swept);
  std::printf("  [adopt] cubocta vs truncate(ambo(cube), eps): flat=%zu "
              "newborn=%zu blend=%zu px\n",
              st.flat, st.newborn, st.blend);
  HS_EXPECT_LE(st.flat, SWAP_FLAT_BUDGET);
  HS_EXPECT_LE(st.newborn, SWAP_NEWBORN_BUDGET);
  HS_EXPECT_LE(st.blend, SWAP_BLEND_BUDGET);
}

/**
 * @brief Verifies the ADOPT bridge arrival swap geometry: the leg's last frame
 *        truncate(tetra, 0.5 - T_EPS) and the adopted octahedron base
 *        (ambo(tetra)) cover the same pixels — under one shared fill, so the
 *        residual is pure T_EPS geometry, not emission order (the ambo
 *        short-circuit reorders emission; palettes carry over per class).
 */
inline void test_adopt_bridge_arrival_geometry() {
  Arena geom(cc_geom_buf, sizeof(cc_geom_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));
  Arena aux(cc_aux_buf, sizeof(cc_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, geom);
  PolyMesh swept = MeshOps::truncate(tetra, geom, temp, 0.5f - T_EPS);
  PolyMesh arrived = MeshOps::ambo(tetra, geom, temp);
  HS_EXPECT_EQ(swept.face_counts.size(), arrived.face_counts.size());

  const std::vector<Color4> fill = {face_color(0)};
  MeshState a_ms, b_ms;
  MeshOps::compile(arrived, a_ms, aux, temp);
  MeshOps::compile(swept, b_ms, aux, temp);
  std::vector<Pixel> a_px, b_px;
  render_faces(a_px, a_ms, [](int) { return face_color(0); });
  render_faces(b_px, b_ms, [](int) { return face_color(0); });

  const DiffStats st = diff_stats(a_px, fill, b_px, fill);
  std::printf("  [adopt-arrival] ambo(tetra) vs truncate(tetra, .5-eps): "
              "flat=%zu blend=%zu px\n",
              st.flat, st.blend);
  HS_EXPECT_LE(st.flat, SWAP_FLAT_BUDGET);
  HS_EXPECT_LE(st.blend, SWAP_BLEND_BUDGET);
}

/**
 * @brief Verifies the DUAL_SWAP ambo crossover: truncate(cube, 0.5 - T_EPS)
 *        and truncate(dual(cube), 0.5 - T_EPS) render the same pixels when
 *        faces are colored by their class signature (clean side count) — the
 *        mapping the crossover uses, since emission order flips at the dual.
 */
inline void test_dual_swap_crossover_framebuffer() {
  Arena geom(cc_geom_buf, sizeof(cc_geom_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));
  Arena aux(cc_aux_buf, sizeof(cc_aux_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  PolyMesh octa = Solids::finalize_solid(MeshOps::dual(cube, aux, temp), geom);

  PolyMesh from_cube = MeshOps::truncate(cube, aux, temp, 0.5f - T_EPS);
  PolyMesh from_octa = MeshOps::truncate(octa, geom, temp, 0.5f - T_EPS);

  // Clean side count at the shared cuboctahedron: primaries keep the seed's
  // count, vertex faces the seed's vertex degree.
  const auto class_color = [](const PolyMesh &seed, int f, bool swap_classes) {
    const int F = static_cast<int>(seed.face_counts.size());
    int sides;
    if (f < F) {
      sides = seed.face_counts[f];
    } else {
      // Both seeds here are regular: vertex degree = total corners / vertices.
      sides = static_cast<int>(seed.faces.size() / seed.vertices.size());
    }
    if (swap_classes)
      sides = sides == 3 ? 4 : 3;
    return face_color(sides);
  };
  const std::vector<Color4> fill = {face_color(3), face_color(4)};

  MeshState a_ms, b_ms;
  Arena ms_arena(cc_leg_buf, sizeof(cc_leg_buf));
  MeshOps::compile(from_cube, a_ms, ms_arena, temp);
  std::vector<Pixel> a_px;
  render_faces(a_px, a_ms, [&](int f) { return class_color(cube, f, false); });

  MeshOps::compile(from_octa, b_ms, ms_arena, temp);
  std::vector<Pixel> b_px, b_wrong_px;
  render_faces(b_px, b_ms, [&](int f) { return class_color(octa, f, false); });
  render_faces(b_wrong_px, b_ms,
               [&](int f) { return class_color(octa, f, true); });

  // Both sides sit T_EPS off the clean cuboctahedron in opposite directions,
  // so every square-triangle boundary carries a ~1 px flat strip of genuine
  // geometry mismatch; a wrong class mapping instead flips whole interiors.
  const DiffStats st = diff_stats(a_px, fill, b_px, fill);
  const DiffStats wrong = diff_stats(a_px, fill, b_wrong_px, fill);
  std::printf("  [dual-swap] truncate(cube) vs truncate(dual(cube)) at ambo: "
              "flat=%zu blend=%zu px; wrong classes flat=%zu (floor %zu)\n",
              st.flat, st.blend, wrong.flat, WRONG_MAPPING_FLAT_FLOOR);
  HS_EXPECT_LE(st.flat, DUAL_SWAP_FLAT_BUDGET);
  HS_EXPECT_LE(st.blend, SWAP_BLEND_BUDGET);
  HS_EXPECT_GT(wrong.flat, WRONG_MAPPING_FLAT_FLOOR);
}

// ---------------------------------------------------------------------------
// §7.6c unit items — palette mappings are total and deterministic; the
// crossfade is exact at its endpoints.
// ---------------------------------------------------------------------------

/** Sample coordinates for exact LUT comparisons. */
constexpr float RAMP_SAMPLES[] = {0.0f, 0.37f, 1.0f};
constexpr int NUM_RAMP_SAMPLES = 3;

/**
 * @brief One frame's shading, copied out of the scratch-backed Shading.
 */
struct ShadingSnapshot {
  std::vector<uint8_t> face_ramp; /**< Face -> ramp pair. */
  std::vector<std::array<Color4, NUM_RAMP_SAMPLES>>
      colors; /**< Per-face sampled colors. */
};

/**
 * @brief Asserts two colors are bitwise equal (channels and alpha).
 */
inline void expect_color_eq(const Color4 &a, const Color4 &b) {
  HS_EXPECT_EQ(a.color.r, b.color.r);
  HS_EXPECT_EQ(a.color.g, b.color.g);
  HS_EXPECT_EQ(a.color.b, b.color.b);
  HS_EXPECT_EQ(a.alpha, b.alpha);
}

/**
 * @brief Steps a ConwayMorph one frame, snapshotting the Shading it hands the
 *        draw callback.
 * @param anim Leg under test.
 * @param fx Canvas provider.
 * @param snap Receives the frame's face_ramp and sampled ramp colors.
 */
inline void step_and_snapshot(Animation::ConwayMorph &anim, ContFx &fx,
                              ShadingSnapshot &snap) {
  snap.face_ramp.clear();
  snap.colors.clear();
  {
    Canvas c(fx);
    anim.step(c);
  }
  fx.advance_display();
}

/**
 * @brief Builds the standard test palette handoff arrays for a seed departed
 *        by emission order.
 * @param seed Departed base mesh.
 * @param pal Receives one palette index per face (deterministic pattern).
 * @param sides Receives the per-face side counts.
 */
inline void fill_emission_handoff(const PolyMesh &seed, uint8_t *pal,
                                  uint8_t *sides) {
  for (size_t f = 0; f < seed.face_counts.size(); ++f) {
    pal[f] = static_cast<uint8_t>(f % Animation::ConwayMorph::PALETTES);
    sides[f] = seed.face_counts[f];
  }
}

/**
 * @brief Verifies the crossfade endpoints on an emission-order leg (cube ->
 *        truncatedCube): frame 1 (w = 0) shades every surviving face from the
 *        handed-off palette and every newborn face from its target; the final
 *        frame (w = 1) shades every face from the landed target assignment.
 */
inline void test_crossfade_exact_at_endpoints_emission() {
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(4242u);

  Arena leg(cc_leg_buf, sizeof(cc_leg_buf));
  Arena bank_arena(cc_bank_buf, sizeof(cc_bank_buf));

  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, leg);
  uint8_t pal[16], sides[16];
  fill_emission_handoff(cube, pal, sides);

  const int edge = [] {
    for (int e = 0; e < ConwayGraph::NUM_EDGES; ++e)
      if (ConwayGraph::EDGES[e].from_node == ConwayGraph::CUBE &&
          ConwayGraph::EDGES[e].to_node == ConwayGraph::TRUNCATED_CUBE)
        return e;
    return -1;
  }();
  HS_EXPECT_GE(edge, 0);

  Animation::ConwayMorph::PaletteHandoff handoff{
      &bank.bank, pal, sides, cube.face_counts.size(), false};

  ShadingSnapshot snap;
  auto cb = [&](Canvas &, const MeshState &m,
                const Animation::ConwayMorph::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    snap.face_ramp.assign(sh.face_ramp, sh.face_ramp + sh.faces);
    snap.colors.resize(sh.faces);
    for (size_t f = 0; f < sh.faces; ++f) {
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                   Animation::ConwayMorph::MAX_BLEND_PAIRS);
      for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
        snap.colors[f][s] = sh.ramps[sh.face_ramp[f]].get(RAMP_SAMPLES[s]);
    }
  };

  constexpr int SWEEP = ConwayGraph::SWEEP_FRAMES;
  Animation::ConwayMorph anim(cube, ConwayGraph::EDGES[edge], false, leg, cb,
                              handoff, SWEEP, 0);
  const Animation::ConwayMorph::Landing &landing = anim.landing();
  HS_EXPECT_EQ(landing.primary_faces, cube.face_counts.size());

  ContFx fx;

  // Frame 1: p = 1/SWEEP < 20%, so w == 0 — the mapped from-state, exactly.
  step_and_snapshot(anim, fx, snap);
  HS_EXPECT_EQ(snap.colors.size(), landing.faces);
  for (size_t f = 0; f < snap.colors.size(); ++f) {
    const uint8_t to = landing.to_palette[wrap(
        landing.topology[f], Animation::ConwayMorph::PALETTES)];
    const uint8_t from = f < landing.primary_faces
                             ? pal[f]
                             : to; // newborn faces skip the crossfade
    for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
      expect_color_eq(snap.colors[f][s],
                      bank.bank.entries[from].get(RAMP_SAMPLES[s]));
  }

  // Final frame: p = 1, so w == 1 — the landed target assignment, exactly.
  for (int f = 1; f < SWEEP; ++f)
    step_and_snapshot(anim, fx, snap);
  for (size_t f = 0; f < snap.colors.size(); ++f) {
    const uint8_t to = landing.to_palette[wrap(
        landing.topology[f], Animation::ConwayMorph::PALETTES)];
    for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
      expect_color_eq(snap.colors[f][s],
                      bank.bank.entries[to].get(RAMP_SAMPLES[s]));
  }
}

/**
 * @brief Verifies the class-signature (DUAL_SWAP) mapping end to end on the
 *        cuboctahedron -> octahedron departure: at frame 1 every primary face
 *        (clean side count 3) inherits the departed triangles' palette and
 *        every vertex face (side count 4) the squares' palette.
 */
inline void test_crossfade_class_signature_mapping() {
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(777u);

  Arena leg(cc_leg_buf, sizeof(cc_leg_buf));
  Arena bank_arena(cc_bank_buf, sizeof(cc_bank_buf));
  Arena temp(cc_temp_buf, sizeof(cc_temp_buf));
  Arena aux(cc_aux_buf, sizeof(cc_aux_buf));

  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  // Departed node: the cuboctahedron at the ambo crossover.
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh cubocta = MeshOps::ambo(cube, aux, temp);
  constexpr uint8_t PAL_TRI = 4, PAL_SQ = 2;
  uint8_t pal[32], sides[32];
  HS_EXPECT_LE(cubocta.face_counts.size(), (size_t)32);
  for (size_t f = 0; f < cubocta.face_counts.size(); ++f) {
    sides[f] = cubocta.face_counts[f];
    pal[f] = sides[f] == 3 ? PAL_TRI : PAL_SQ;
  }

  // Departing leg after the dual swap: reverse traversal of octahedron ->
  // cuboctahedron on the octahedron seed.
  const int edge = [] {
    for (int e = 0; e < ConwayGraph::NUM_EDGES; ++e)
      if (ConwayGraph::EDGES[e].from_node == ConwayGraph::OCTAHEDRON &&
          ConwayGraph::EDGES[e].to_node == ConwayGraph::CUBOCTAHEDRON &&
          ConwayGraph::EDGES[e].reseed == ConwayGraph::Reseed::DUAL_SWAP)
        return e;
    return -1;
  }();
  HS_EXPECT_GE(edge, 0);

  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, leg);
  Animation::ConwayMorph::PaletteHandoff handoff{
      &bank.bank, pal, sides, cubocta.face_counts.size(), true};

  ShadingSnapshot snap;
  auto cb = [&](Canvas &, const MeshState &,
                const Animation::ConwayMorph::Shading &sh) {
    snap.face_ramp.assign(sh.face_ramp, sh.face_ramp + sh.faces);
    snap.colors.resize(sh.faces);
    for (size_t f = 0; f < sh.faces; ++f)
      for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
        snap.colors[f][s] = sh.ramps[sh.face_ramp[f]].get(RAMP_SAMPLES[s]);
  };

  Animation::ConwayMorph anim(octa, ConwayGraph::EDGES[edge], /*reverse*/ true,
                              leg, cb, handoff, ConwayGraph::SWEEP_FRAMES, 0);
  const Animation::ConwayMorph::Landing &landing = anim.landing();
  HS_EXPECT_EQ(landing.primary_faces, octa.face_counts.size());

  ContFx fx;
  step_and_snapshot(anim, fx, snap); // frame 1: w == 0
  HS_EXPECT_EQ(snap.colors.size(), landing.faces);
  for (size_t f = 0; f < snap.colors.size(); ++f) {
    // Primary faces are the octahedron's triangles; the rest are the truncate
    // vertex faces (squares at the octahedron's 4-valent vertices).
    const uint8_t from = f < landing.primary_faces ? PAL_TRI : PAL_SQ;
    for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
      expect_color_eq(snap.colors[f][s],
                      bank.bank.entries[from].get(RAMP_SAMPLES[s]));
  }
}

/**
 * @brief Verifies every edge's palette mapping is total: at frame 1 each
 *        face's ramp resolves to a real bank entry — the handed-off palette
 *        for surviving primaries, the landed target for newborn faces — and
 *        the landed assignment is a permutation of the bank slots.
 */
inline void test_palette_mapping_total_all_edges() {
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);

  Arena bank_arena(cc_bank_buf, sizeof(cc_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  ContFx fx;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    const int failed_before = hs_test::stats().failed;

    Arena leg(cc_leg_buf, sizeof(cc_leg_buf));
    Arena seed_a(cc_geom_buf, sizeof(cc_geom_buf));
    Arena seed_b(cc_temp_buf, sizeof(cc_temp_buf));
    PolyMesh seed =
        Solids::simple_registry[e.seed_solid].generate(seed_a, seed_b);

    uint8_t pal[128], sides[128];
    HS_EXPECT_LE(seed.face_counts.size(), (size_t)128);
    // Class-keyed handoff (faces of one side count share a palette), as the
    // effect hands off: per-face-arbitrary palettes would exceed the leg's
    // MAX_BLEND_PAIRS pair budget on the large seeds.
    for (size_t f = 0; f < seed.face_counts.size(); ++f) {
      sides[f] = seed.face_counts[f];
      pal[f] =
          static_cast<uint8_t>(sides[f] % Animation::ConwayMorph::PALETTES);
    }
    Animation::ConwayMorph::PaletteHandoff handoff{
        &bank.bank, pal, sides, seed.face_counts.size(), false};

    ShadingSnapshot snap;
    auto cb = [&](Canvas &, const MeshState &,
                  const Animation::ConwayMorph::Shading &sh) {
      snap.face_ramp.assign(sh.face_ramp, sh.face_ramp + sh.faces);
      snap.colors.resize(sh.faces);
      for (size_t f = 0; f < sh.faces; ++f) {
        HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                     Animation::ConwayMorph::MAX_BLEND_PAIRS);
        for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
          snap.colors[f][s] = sh.ramps[sh.face_ramp[f]].get(RAMP_SAMPLES[s]);
      }
    };

    hs::random().seed(1000u + static_cast<uint32_t>(ei));
    Animation::ConwayMorph anim(seed, e, false, leg, cb, handoff,
                                ConwayGraph::SWEEP_FRAMES,
                                e.settle ? ConwayGraph::SETTLE_FRAMES : 0);
    const Animation::ConwayMorph::Landing &landing = anim.landing();

    // The landed assignment is a permutation of the bank slots.
    std::array<uint8_t, Animation::ConwayMorph::PALETTES> perm =
        landing.to_palette;
    std::sort(perm.begin(), perm.end());
    for (int i = 0; i < Animation::ConwayMorph::PALETTES; ++i)
      HS_EXPECT_EQ(static_cast<int>(perm[i]), i);

    step_and_snapshot(anim, fx, snap); // frame 1: w == 0
    HS_EXPECT_EQ(snap.colors.size(), landing.faces);
    for (size_t f = 0; f < snap.colors.size(); ++f) {
      const uint8_t to = landing.to_palette[wrap(
          landing.topology[f], Animation::ConwayMorph::PALETTES)];
      const uint8_t from = f < landing.primary_faces ? pal[f] : to;
      for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
        expect_color_eq(snap.colors[f][s],
                        bank.bank.entries[from].get(RAMP_SAMPLES[s]));
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [mapping] edge %d: %s -> %s\n", ei,
                  Solids::simple_registry[e.from_node].name,
                  Solids::simple_registry[e.to_node].name);
  }
}

/**
 * @brief Verifies the palette mapping is deterministic: two legs constructed
 *        from the same RNG seed and inputs land identical target assignments,
 *        classifications, and frame-1 shadings.
 */
inline void test_palette_mapping_deterministic() {
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);

  Arena bank_arena(cc_bank_buf, sizeof(cc_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  const int edge = [] {
    for (int e = 0; e < ConwayGraph::NUM_EDGES; ++e)
      if (ConwayGraph::EDGES[e].from_node == ConwayGraph::CUBE &&
          ConwayGraph::EDGES[e].to_node == ConwayGraph::SNUB_CUBE)
        return e;
    return -1;
  }();
  HS_EXPECT_GE(edge, 0);

  ContFx fx;
  std::array<uint8_t, Animation::ConwayMorph::PALETTES> to_palette[2];
  std::vector<int> topo[2];
  ShadingSnapshot snaps[2];

  for (int run = 0; run < 2; ++run) {
    Arena leg(cc_leg_buf, sizeof(cc_leg_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, leg);
    uint8_t pal[16], sides[16];
    fill_emission_handoff(cube, pal, sides);
    Animation::ConwayMorph::PaletteHandoff handoff{
        &bank.bank, pal, sides, cube.face_counts.size(), false};

    ShadingSnapshot &snap = snaps[run];
    auto cb = [&](Canvas &, const MeshState &,
                  const Animation::ConwayMorph::Shading &sh) {
      snap.face_ramp.assign(sh.face_ramp, sh.face_ramp + sh.faces);
      snap.colors.resize(sh.faces);
      for (size_t f = 0; f < sh.faces; ++f)
        for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
          snap.colors[f][s] = sh.ramps[sh.face_ramp[f]].get(RAMP_SAMPLES[s]);
    };

    hs::random().seed(31337u);
    Animation::ConwayMorph anim(cube, ConwayGraph::EDGES[edge], false, leg, cb,
                                handoff, ConwayGraph::SWEEP_FRAMES,
                                ConwayGraph::SETTLE_FRAMES);
    const Animation::ConwayMorph::Landing &landing = anim.landing();
    to_palette[run] = landing.to_palette;
    topo[run].assign(landing.topology, landing.topology + landing.faces);
    step_and_snapshot(anim, fx, snap);
  }

  HS_EXPECT_TRUE(to_palette[0] == to_palette[1]);
  HS_EXPECT_TRUE(topo[0] == topo[1]);
  HS_EXPECT_TRUE(snaps[0].face_ramp == snaps[1].face_ramp);
  HS_EXPECT_EQ(snaps[0].colors.size(), snaps[1].colors.size());
  for (size_t f = 0; f < snaps[0].colors.size(); ++f)
    for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
      expect_color_eq(snaps[0].colors[f][s], snaps[1].colors[f][s]);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all boundary color-continuity tests.
 * @return The module's failure count.
 */
inline int run_conway_continuity_tests() {
  hs_test::ModuleFixture fixture("conway_continuity");

  test_bookend_angle_pin();
  test_bookend_swaps_per_family();
  test_palette_carry_across_arrivals();

  test_swap_in_framebuffer();
  test_adopt_departure_swap_framebuffer();
  test_adopt_bridge_arrival_geometry();
  test_dual_swap_crossover_framebuffer();

  test_crossfade_exact_at_endpoints_emission();
  test_crossfade_class_signature_mapping();
  test_palette_mapping_total_all_edges();
  test_palette_mapping_deterministic();

  return fixture.result();
}

} // namespace conway_continuity_tests
} // namespace hs_test
