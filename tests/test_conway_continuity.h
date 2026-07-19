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
 *   - Strap-slot crossfade across cycle starts (boot seed plus an eight-seed
 *     epoch sweep): a strap-bearing slot opens on the color it displayed in
 *     the previous cycle (or an on-screen star palette when newborn) and
 *     glides to its target in bounded steps via the strap-face LUT; star
 *     faces stay bitwise on the bank entry, including on star-shared slots.
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

  // Quantization-level channel flips are tolerated: the star 2n-gon's edge
  // planes are fp-distinct from the base n-gon's, so a boundary AA blend can
  // round up to 2 counts apart (a 1-ulp coverage difference through the
  // 16-bit lerp); a recolored face moves channels by thousands.
  constexpr int AA_LSB_TOL = 2;
  size_t diff = 0;
  size_t quant_only = 0;
  for (size_t i = 0; i < base_px.size() && i < flat_px.size(); ++i) {
    const int dr = std::abs((int)base_px[i].r - (int)flat_px[i].r);
    const int dg = std::abs((int)base_px[i].g - (int)flat_px[i].g);
    const int db = std::abs((int)base_px[i].b - (int)flat_px[i].b);
    if (dr == 0 && dg == 0 && db == 0)
      continue;
    if (dr <= AA_LSB_TOL && dg <= AA_LSB_TOL && db <= AA_LSB_TOL)
      ++quant_only;
    else
      ++diff;
  }
  std::printf("  [bookend] F=%d: base vs hankin@0 diff=%zu px "
              "(+%zu quantization-level)\n",
              F, diff, quant_only);
  HS_EXPECT_EQ(diff, (size_t)0);

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
// Leg-start seed-frame continuity: over a scripted walk covering every SeedFix
// path (KEEP, DUAL_SWAP both families, DERIVE_AMBO both directions, all three
// REGEN_TETRA reverse bridges, dual-swap wandering before a REGEN), the mesh
// drawn on the first morph frame must overlie the departed node mesh face for
// face — geometry within tolerance and every face's w = 0 color inherited from
// the departed face it covers. Fails on a seed rebuilt in a frame or order the
// held mesh does not have (the orientation/color snap this pins against).
// ---------------------------------------------------------------------------

/** Scripted edge walk (indices into ConwayGraph::EDGES). Bridges out-and-back
 * (edges 19/20/21), dual-swap wandering in both families before each REGEN,
 * DERIVE_AMBO out-and-back (edge 17), settle legs both directions (15/16/21),
 * expand and snub departures (3/4/16/21). */
constexpr int SEEDFRAME_SCRIPT[] = {19, 8,  1,  0,  2,  8,  19, 21, 14, 10,
                                    15, 15, 16, 16, 10, 17, 17, 14, 21, 18,
                                    20, 20, 20, 8,  1,  3,  3,  4,  4};

/** Two visibly distinct bank entries alternated across the departed faces so a
 * misrouted provenance flips a face between them. */
constexpr uint8_t SEEDFRAME_PAL_A = 1;
constexpr uint8_t SEEDFRAME_PAL_B = 4;

/**
 * @brief Unit-sphere vertex-average centroid of face fi.
 */
inline Vector poly_face_centroid(const PolyMesh &m, size_t fi) {
  size_t off = 0;
  for (size_t i = 0; i < fi; ++i)
    off += m.face_counts[i];
  Vector c(0.0f, 0.0f, 0.0f);
  for (int k = 0; k < m.face_counts[fi]; ++k)
    c = c + m.vertices[m.faces[off + k]];
  return c.normalized();
}

/**
 * @brief The clean node mesh at one end of an edge, built from the held seed —
 * HankinSolids::node_mesh_at, replicated for the walk simulation.
 */
inline PolyMesh seedframe_node_mesh_at(const ConwayGraph::EdgeSpec &e,
                                       bool to_end, const PolyMesh &seed_base,
                                       Arena &a, Arena &b) {
  float t = to_end ? e.t_to : e.t_from;
  PolyMesh seed;
  MeshOps::clone(seed_base, seed, a);
  Solids::SolidBuilder builder(std::move(seed), a, b);
  if (!ConwayGraph::is_platonic(e.seed_solid))
    builder.ambo();
  if (t > 0.0f) {
    switch (e.op) {
    case ConwayGraph::MorphOp::TRUNCATE:
      builder.truncate(t);
      break;
    case ConwayGraph::MorphOp::EXPAND:
      builder.expand(t);
      break;
    case ConwayGraph::MorphOp::SNUB:
      builder.snub(t, to_end ? e.twist_to : e.twist_from);
      break;
    }
    if (e.settle && to_end)
      builder.relax(50);
  }
  return builder.build();
}

/**
 * @brief Drives the seed state machine over SEEDFRAME_SCRIPT and pins, per
 * leg start, the geometric face bijection and the inherited w = 0 shading.
 */
inline void test_leg_start_seed_frame_continuity() {
  using namespace ConwayGraph;
  using conway_morph_tests::run_edge_op;
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(90210u);

  Arena bank_arena(cc_bank_buf, sizeof(cc_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  // The alternating entries must actually differ or misroutes are invisible.
  {
    bool differ = false;
    for (int s = 0; s < NUM_RAMP_SAMPLES; ++s) {
      const Color4 a = bank.bank.entries[SEEDFRAME_PAL_A].get(RAMP_SAMPLES[s]);
      const Color4 b = bank.bank.entries[SEEDFRAME_PAL_B].get(RAMP_SAMPLES[s]);
      if (a.color.r != b.color.r || a.color.g != b.color.g ||
          a.color.b != b.color.b)
        differ = true;
    }
    HS_EXPECT_TRUE(differ);
  }

  constexpr float TOL_SQ = 0.15f * 0.15f;

  Arena seed_arena(cc_geom_buf, sizeof(cc_geom_buf));
  PolyMesh seed_base;
  {
    Arena wa(cc_temp_buf, sizeof(cc_temp_buf));
    Arena wb(cc_aux_buf, sizeof(cc_aux_buf));
    seed_base = Solids::finalize_solid(Solids::Platonic::tetrahedron(wa, wb),
                                       seed_arena);
  }
  int node = TETRAHEDRON;
  int seed_identity = TETRAHEDRON;

  Arena node_arena(cc_scan_buf, sizeof(cc_scan_buf));
  PolyMesh node_mesh = Solids::finalize_solid(seed_base, node_arena);

  ContFx fx;
  for (size_t leg = 0; leg < std::size(SEEDFRAME_SCRIPT); ++leg) {
    const int ei = SEEDFRAME_SCRIPT[leg];
    const EdgeSpec &e = EDGES[ei];
    const bool reverse = (e.to_node == node);
    const int failed_before = hs_test::stats().failed;

    // Seed reconciliation, exactly as start_morph_cycle applies it.
    const SeedFix fix = seed_fix_at_start(ei, seed_identity);
    HS_EXPECT_TRUE(fix != SeedFix::INVALID);
    Arena work(cc_temp_buf, sizeof(cc_temp_buf));
    Arena temp(cc_aux_buf, sizeof(cc_aux_buf));
    if (fix == SeedFix::DUAL_SWAP) {
      PolyMesh d = MeshOps::dual(seed_base, work, temp);
      seed_arena.reset();
      seed_base = Solids::finalize_solid(d, seed_arena);
      seed_identity = dual_platonic(seed_identity);
    } else if (fix == SeedFix::REGEN_TETRA) {
      PolyMesh t = Solids::Platonic::tetrahedron(work, temp);
      seed_arena.reset();
      seed_base = Solids::finalize_solid(t, seed_arena);
      seed_identity = TETRAHEDRON;
    }
    PolyMesh derived;
    if (fix == SeedFix::DERIVE_AMBO)
      derived = MeshOps::ambo(seed_base, work, temp);
    const PolyMesh &leg_seed =
        fix == SeedFix::DERIVE_AMBO ? derived : seed_base;

    // Departed-node handoff: alternating palettes, real sides/centroids.
    const size_t prev_faces = node_mesh.face_counts.size();
    uint8_t pal[128], sides[128];
    Vector cents[128];
    HS_EXPECT_LE(prev_faces, (size_t)128);
    for (size_t f = 0; f < prev_faces; ++f) {
      pal[f] = (f % 2) ? SEEDFRAME_PAL_B : SEEDFRAME_PAL_A;
      sides[f] = node_mesh.face_counts[f];
      cents[f] = poly_face_centroid(node_mesh, f);
    }
    Animation::ConwayMorph::PaletteHandoff handoff{
        &bank.bank, pal, sides, prev_faces, fix == SeedFix::DUAL_SWAP, cents};

    // The mesh the first morph frame draws.
    auto clampp = [&](float t) {
      t = std::max(t, T_EPS);
      if (e.op == MorphOp::TRUNCATE)
        t = std::min(t, 0.5f - T_EPS_AMBO);
      return t;
    };
    const float t_start = clampp(reverse ? e.t_to : e.t_from);
    const float tw_start = reverse ? e.twist_to : e.twist_from;
    PolyMesh start_raw =
        run_edge_op(e, leg_seed, work, temp, t_start, tw_start);
    PolyMesh start_rel;
    if (e.settle && reverse)
      start_rel = MeshOps::relax(start_raw, work, temp, 50);
    const PolyMesh &start = (e.settle && reverse) ? start_rel : start_raw;
    const size_t total = start.face_counts.size();
    const size_t primary = leg_seed.face_counts.size();

    // Geometric pin: every start face overlies a departed face (full legs:
    // a bijection) — the seed-frame continuity contract.
    int match_of[128];
    bool used[128] = {};
    HS_EXPECT_LE(total, (size_t)128);
    for (size_t f = 0; f < total; ++f) {
      const Vector c = poly_face_centroid(start, f);
      size_t best = 0;
      float best_d = 1e9f;
      for (size_t g = 0; g < prev_faces; ++g) {
        const Vector d = c - cents[g];
        const float dsq = dot(d, d);
        if (dsq < best_d) {
          best_d = dsq;
          best = g;
        }
      }
      match_of[f] = static_cast<int>(best);
      if (prev_faces == total) {
        HS_EXPECT_LT(best_d, TOL_SQ);
        HS_EXPECT_TRUE(!used[best]);
        used[best] = true;
      } else if (f < primary) {
        // Seed departures: primaries still overlie their own faces.
        HS_EXPECT_LT(best_d, TOL_SQ);
        HS_EXPECT_EQ(match_of[f], (int)f);
      }
    }

    // Shading pin: the first frame's w = 0 colors are the departed faces'.
    ShadingSnapshot snap;
    auto cb = [&](Canvas &, const MeshState &,
                  const Animation::ConwayMorph::Shading &sh) {
      snap.face_ramp.assign(sh.face_ramp, sh.face_ramp + sh.faces);
      snap.colors.resize(sh.faces);
      for (size_t f = 0; f < sh.faces; ++f)
        for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
          snap.colors[f][s] = sh.ramps[sh.face_ramp[f]].get(RAMP_SAMPLES[s]);
    };
    {
      Arena leg_arena(cc_leg_buf, sizeof(cc_leg_buf));
      Animation::ConwayMorph anim(leg_seed, e, reverse, leg_arena, cb, handoff,
                                  SWEEP_FRAMES, e.settle ? SETTLE_FRAMES : 0);
      step_and_snapshot(anim, fx, snap);
      HS_EXPECT_EQ(snap.colors.size(), total);
      for (size_t f = 0; f < snap.colors.size(); ++f) {
        if (prev_faces == total || f < primary) {
          const uint8_t want = prev_faces == total ? pal[match_of[f]] : pal[f];
          for (int s = 0; s < NUM_RAMP_SAMPLES; ++s)
            expect_color_eq(snap.colors[f][s],
                            bank.bank.entries[want].get(RAMP_SAMPLES[s]));
        } else {
          // Births must open in a color already on screen (a departed
          // palette), never pop in as a fresh target color.
          bool matches_a = true, matches_b = true;
          for (int s = 0; s < NUM_RAMP_SAMPLES; ++s) {
            const Color4 got = snap.colors[f][s];
            const Color4 a =
                bank.bank.entries[SEEDFRAME_PAL_A].get(RAMP_SAMPLES[s]);
            const Color4 b =
                bank.bank.entries[SEEDFRAME_PAL_B].get(RAMP_SAMPLES[s]);
            if (got.color.r != a.color.r || got.color.g != a.color.g ||
                got.color.b != a.color.b)
              matches_a = false;
            if (got.color.r != b.color.r || got.color.g != b.color.g ||
                got.color.b != b.color.b)
              matches_b = false;
          }
          HS_EXPECT_TRUE(matches_a || matches_b);
        }
      }
    }

    // Completion: arrival node mesh becomes the next bookend; ADOPT as tabled.
    const bool arrived_at_to = !reverse;
    const int arrived = arrived_at_to ? e.to_node : e.from_node;
    PolyMesh base =
        seedframe_node_mesh_at(e, arrived_at_to, seed_base, work, temp);
    node_arena.reset();
    node_mesh = Solids::finalize_solid(base, node_arena);
    if (e.reseed == Reseed::ADOPT && is_platonic(arrived) && arrived_at_to) {
      seed_arena.reset();
      seed_base = Solids::finalize_solid(node_mesh, seed_arena);
      seed_identity = arrived;
    }
    node = arrived;

    if (hs_test::stats().failed != failed_before)
      std::printf("    [seed-frame] leg %zu edge %d -> '%s' broke\n", leg, ei,
                  Solids::simple_registry[node].name);
  }
  // The script must end where cycle-accounting expects it: back at the cube
  // pendant chain (guards against a silently mis-scripted walk).
  HS_EXPECT_EQ(node, (int)CUBE);
}

// ---------------------------------------------------------------------------
// In-cycle palette stability: the live class-slot assignment the hankin cycle
// draws with may change only at leg completion (finish_morph_cycle), never at
// leg construction or mid-cycle — a live rewrite would recolor the still-
// visible rosette straps before the morph starts.
// ---------------------------------------------------------------------------

/**
 * @brief Pins palette_idx_ constant from each arrival to the next, including
 * across the leg-construction frame.
 */
inline void test_palette_slots_stable_within_cycle() {
  reset_globals();
  using Probe = conway_soak_tests::HankinWalkProbe;

  HankinSolids<96, 20> fx;
  fx.init();

  auto snapshot = Probe::palette_idx(fx);
  int prev_node = Probe::node(fx);
  int arrivals = 0;
  constexpr int TARGET_ARRIVALS = 10;
  for (int frame = 0; frame < 2600 && arrivals < TARGET_ARRIVALS; ++frame) {
    fx.draw_frame();
    fx.advance_display();

    const int node = Probe::node(fx);
    if (node != prev_node) {
      // Leg completion ran this frame: the assignment may legitimately move.
      prev_node = node;
      snapshot = Probe::palette_idx(fx);
      ++arrivals;
      continue;
    }
    if (Probe::palette_idx(fx) != snapshot) {
      std::printf("    [palette-stability] live palette rewrite mid-cycle at "
                  "frame %d ('%s')\n",
                  frame, Solids::simple_registry[node].name);
      HS_EXPECT_TRUE(false);
      return;
    }
  }
  HS_EXPECT_EQ(arrivals, TARGET_ARRIVALS);
}

// ---------------------------------------------------------------------------
// Strap-slot crossfade across cycle starts: a strap-bearing class slot opens
// each hankin cycle at the color its slot displayed in the previous cycle (or
// at a star palette already on screen when the slot had no predecessor) and
// glides to its target in bounded per-frame steps via the strap-face LUT;
// star faces resolve to the assignment's exact bank entry at every frame,
// keeping the bookend star colors bitwise exact even on slots the mod-5 wrap
// shares between a star and a rosette class. Pre-fix, reborn strap slots
// opened directly on the fresh shuffle: an open-vs-previous-close LUT jump of
// up to ~57000/65535 per channel (measured as a 60210 mean-strap-color pop
// within ~4 frames of the bookend on the 30-leg 288x144 strap-snap harness
// tour); star-shared slots stayed exempt and popped by the full carried
// distance (60098 at truncatedIcosidodecahedron, epoch seed 7, same harness).
// ---------------------------------------------------------------------------

/** Ceiling on one frame's smoothstep advance over the 20-frame window (max
 * slope 1.5/20 = 0.075), applied to the armed pair's own LUT distance. */
constexpr float STRAP_STEP_FRAC = 0.08f;

/** Rounding slack for LUT-quantized blend steps, in 16-bit channel counts. */
constexpr int STRAP_STEP_SLACK = 64;

/**
 * @brief Max-channel distance of two baked LUTs over the standard samples.
 */
inline int lut_sample_dist(const BakedPalette &a, const BakedPalette &b) {
  int d = 0;
  for (int s = 0; s < NUM_RAMP_SAMPLES; ++s) {
    const Color4 ca = a.get(RAMP_SAMPLES[s]);
    const Color4 cb = b.get(RAMP_SAMPLES[s]);
    d = std::max(d, std::abs((int)ca.color.r - (int)cb.color.r));
    d = std::max(d, std::abs((int)ca.color.g - (int)cb.color.g));
    d = std::max(d, std::abs((int)ca.color.b - (int)cb.color.b));
  }
  return d;
}

/** Per-seed result of check_strap_crossfade_arrivals. */
struct StrapSweepStats {
  int arrivals = 0;      /**< Arrivals validated. */
  int far_pairs = 0;     /**< Strap slots with a far (from, to) pair. */
  int shared_far = 0;    /**< Far pairs on slots shared with a star class. */
  int would_be_jump = 0; /**< Worst uncrossfaded open jump (16-bit). */
  int shared_jump = 0;   /**< Worst such jump on a star-shared slot. */
};

/**
 * @brief Drives one HankinSolids walk across arrivals and pins the strap
 *        crossfade: open state continuous with the previous cycle, endpoints
 *        bitwise exact, per-frame steps bounded, star faces on the exact bank
 *        entry at the bookends — including slots a star class shares.
 * @param epoch Epoch seed for hs::random() (the live per-visit reseed).
 * @param target_arrivals Arrivals to validate before returning.
 * @param frame_cap Frame budget for the walk.
 * @return Coverage stats for the caller's discrimination pins.
 */
inline StrapSweepStats check_strap_crossfade_arrivals(uint32_t epoch,
                                                      int target_arrivals,
                                                      int frame_cap) {
  reset_globals();
  hs::random().seed(hs::epoch_seed(epoch));
  using Probe = conway_soak_tests::HankinWalkProbe;
  constexpr int PALETTES = Animation::ConwayMorph::PALETTES;

  HankinSolids<96, 20> fx;
  fx.init();
  const BakedPaletteBank &bank = Probe::palette_bank(fx).bank;
  const int BLEND = Probe::strap_blend_frames(fx);
  Arena resolve_arena(cc_temp_buf, sizeof(cc_temp_buf));

  StrapSweepStats st;
  std::array<int, PALETTES> prev_display{};
  bool prev_used[PALETTES] = {};
  bool have_prev = false;
  int prev_node = Probe::node(fx);

  for (int frame = 0; frame < frame_cap && st.arrivals < target_arrivals;
       ++frame) {
    fx.draw_frame();
    fx.advance_display();
    const int node = Probe::node(fx);
    if (node == prev_node)
      continue;
    prev_node = node;
    ++st.arrivals;
    const int failed_before = hs_test::stats().failed;

    // The arrival frame is the new cycle's opening bookend.
    HS_EXPECT_LE(Probe::hankin_cycle_frame(fx), 1);

    // Slot roles at the new node.
    const MeshState &mesh = Probe::mesh(fx);
    const size_t nf = Probe::node_faces(fx);
    bool star[PALETTES] = {}, strap[PALETTES] = {};
    for (size_t f = 0; f < mesh.topology.size(); ++f) {
      const int slot = wrap(mesh.topology[f], PALETTES);
      (f < nf ? star[slot] : strap[slot]) = true;
    }
    const uint8_t mask = Probe::strap_blend_mask(fx);
    const auto &from = Probe::strap_from(fx);
    const auto &idx = Probe::palette_idx(fx);

    BakedPalette blended[PALETTES];
    const BakedPalette *star_by_slot[PALETTES];
    const BakedPalette *strap_by_slot[PALETTES];
    for (int s = 0; s < PALETTES; ++s) {
      // Only strap-bearing slots arm; star faces resolve to the assignment's
      // bank entry, bitwise, at both bookend endpoints of the window.
      if (!strap[s])
        HS_EXPECT_EQ((mask >> s) & 1, 0);
      if (star[s]) {
        for (int cf : {0, BLEND}) {
          resolve_arena.reset();
          Probe::resolve_slot_luts(fx, cf, blended, star_by_slot, strap_by_slot,
                                   resolve_arena);
          HS_EXPECT_TRUE(star_by_slot[s] == &bank.entries[idx[s]]);
        }
      }
      if (!strap[s])
        continue;

      // Crossfade endpoints, on the production resolver, bitwise: the open
      // frame reproduces the from state, the window end the target.
      resolve_arena.reset();
      Probe::resolve_slot_luts(fx, 0, blended, star_by_slot, strap_by_slot,
                               resolve_arena);
      for (int k = 0; k < NUM_RAMP_SAMPLES; ++k)
        expect_color_eq(strap_by_slot[s]->get(RAMP_SAMPLES[k]),
                        bank.entries[from[s]].get(RAMP_SAMPLES[k]));
      resolve_arena.reset();
      Probe::resolve_slot_luts(fx, BLEND, blended, star_by_slot, strap_by_slot,
                               resolve_arena);
      for (int k = 0; k < NUM_RAMP_SAMPLES; ++k)
        expect_color_eq(strap_by_slot[s]->get(RAMP_SAMPLES[k]),
                        bank.entries[idx[s]].get(RAMP_SAMPLES[k]));

      if (have_prev && prev_used[s]) {
        // Continuity across the cycle start: the slot's straps open on the
        // color the slot displayed when the previous cycle closed. (Pre-fix,
        // rosette-only slots opened on idx[s] and star-shared slots were
        // exempt from arming entirely -- the would-be jump below.)
        HS_EXPECT_EQ(from[s], prev_display[s]);
        const int jump = lut_sample_dist(bank.entries[prev_display[s]],
                                         bank.entries[idx[s]]);
        st.would_be_jump = std::max(st.would_be_jump, jump);
        if (jump > 8192) {
          ++st.far_pairs;
          if (star[s]) {
            ++st.shared_far;
            st.shared_jump = std::max(st.shared_jump, jump);
          }
        }
      } else {
        // No predecessor: births open in a color already on screen (some
        // star face's palette), never a fresh pop.
        bool on_screen = false;
        for (size_t f = 0; f < nf; ++f)
          on_screen |= Probe::node_face_palette(fx)[f] == from[s];
        HS_EXPECT_TRUE(on_screen);
      }

      // Bounded per-frame steps across the whole opening window.
      const int pair_dist =
          lut_sample_dist(bank.entries[from[s]], bank.entries[idx[s]]);
      const int step_bound =
          static_cast<int>(STRAP_STEP_FRAC * pair_dist) + STRAP_STEP_SLACK;
      std::array<Color4, NUM_RAMP_SAMPLES> prev_c{};
      for (int cf = 0; cf <= BLEND; ++cf) {
        resolve_arena.reset();
        Probe::resolve_slot_luts(fx, cf, blended, star_by_slot, strap_by_slot,
                                 resolve_arena);
        for (int k = 0; k < NUM_RAMP_SAMPLES; ++k) {
          const Color4 c = strap_by_slot[s]->get(RAMP_SAMPLES[k]);
          if (cf > 0) {
            HS_EXPECT_LE(std::abs((int)c.color.r - (int)prev_c[k].color.r),
                         step_bound);
            HS_EXPECT_LE(std::abs((int)c.color.g - (int)prev_c[k].color.g),
                         step_bound);
            HS_EXPECT_LE(std::abs((int)c.color.b - (int)prev_c[k].color.b),
                         step_bound);
          }
          prev_c[k] = c;
        }
      }
    }

    // This cycle's display becomes the next arrival's from state.
    prev_display = idx;
    for (int s = 0; s < PALETTES; ++s)
      prev_used[s] = star[s] || strap[s];
    have_prev = true;

    if (hs_test::stats().failed != failed_before)
      std::printf("    [strap-crossfade] epoch %u arrival %d at '%s' broke\n",
                  epoch, st.arrivals, Solids::simple_registry[node].name);
  }
  return st;
}

/**
 * @brief Single-walk strap-crossfade pin at the boot seed (epoch 0).
 * @details Red pre-crossfade via the continuity pin: without it a strap
 *          slot's opening LUT is the fresh target, which jumps from the
 *          previous cycle's display by the full palette distance (printed as
 *          the would-be jump; the run must exercise at least one such far
 *          pair for the pin to discriminate).
 */
inline void test_strap_crossfade_across_cycle_start() {
  constexpr int TARGET_ARRIVALS = 10;
  const StrapSweepStats st =
      check_strap_crossfade_arrivals(0, TARGET_ARRIVALS, 2600);
  HS_EXPECT_EQ(st.arrivals, TARGET_ARRIVALS);
  std::printf("  [strap-crossfade] worst would-be pre-fix open jump %d "
              "(16-bit max channel; crossfaded to per-frame steps)\n",
              st.would_be_jump);
  // The walk must have exercised at least one far (from, to) pair, or the
  // continuity pin proved nothing this run.
  HS_EXPECT_GT(st.far_pairs, 0);
}

/** Epoch seeds for the seed-swept strap pin. Every seed's walk is
 * deterministic; 3 and 8 reach the truncatedIcosidodecahedron (the one
 * roster node whose mod-5 wrap puts star faces in all five slots) by arrival
 * 18, and 13 and 15 by arrival 26, so the sweep always exercises star-shared
 * far pairs. */
constexpr uint32_t STRAP_SWEEP_EPOCHS[] = {0, 1, 2, 3, 5, 8, 13, 15};

/** Arrival budget per swept seed (frame cap scales with it). */
constexpr int STRAP_SWEEP_ARRIVALS[] = {6, 6, 6, 18, 6, 18, 26, 26};

/**
 * @brief Seed-swept strap continuity: across eight epoch seeds, every strap
 *        slot of every arrival opens on its previous displayed color and
 *        glides in bounded steps — including slots shared with a star class,
 *        whose star faces stay bitwise on the bank entry at the bookends.
 * @details Red pre-fix: star-shared slots were exempt from arming, so their
 *          straps reopened bitwise on the carried star palette — a full-
 *          distance pop whenever the epoch's shuffle moved the slot (13
 *          continuity-pin failures across the sweep, worst 65532/65535 at
 *          epochs 3 and 15, arriving icosidodecahedron ->
 *          truncatedIcosidodecahedron).
 */
inline void test_strap_crossfade_seed_swept() {
  int shared_far = 0;
  int shared_jump = 0;
  for (size_t i = 0; i < std::size(STRAP_SWEEP_EPOCHS); ++i) {
    const int want = STRAP_SWEEP_ARRIVALS[i];
    const StrapSweepStats st = check_strap_crossfade_arrivals(
        STRAP_SWEEP_EPOCHS[i], want, 140 * (want + 1));
    HS_EXPECT_EQ(st.arrivals, want);
    shared_far += st.shared_far;
    shared_jump = std::max(shared_jump, st.shared_jump);
  }
  std::printf("  [strap-sweep] %d star-shared far pairs, worst would-be "
              "shared open jump %d (crossfaded to per-frame steps)\n",
              shared_far, shared_jump);
  // The sweep must exercise the star-shared aliasing it exists to pin.
  HS_EXPECT_GT(shared_far, 0);
}

// ---------------------------------------------------------------------------
// Strap opening-fade: newborn straps reveal over the opening window instead of
// popping into former star interiors in one frame.
//
// At the opening bookend (angle 0) the straps are zero-area and the frame is
// the base solid's star faces. One frame later the interlace angle steps to
// ~0.0038 rad and the strap faces are born as thin slivers cutting through
// those star interiors; drawn at full coverage they hard-recolor the cut
// pixels (a strap-vs-star sliver pop, measured 1208 px of 41472 at 288x144 on
// truncatedTetrahedron, worst ~3000 px across the tour). Fading strap coverage
// by strap_blend_weight over the window (~0.007 at the first frame) makes the
// birth a bounded reveal: the first strap frame stays close to the bookend and
// the straps ease to full coverage as the window closes.
// ---------------------------------------------------------------------------

/** Per-channel delta above which a pixel counts as a hard recolor (not an AA
 * boundary shimmer), in 16-bit channel counts. */
constexpr int STRAP_OPEN_HARD = 8000;

/** First interlace angle a hankin cycle draws (sin_wave(0, pi/2, 1, 0) at
 * progress 1/64); see start_hankin_cycle. */
constexpr float STRAP_OPEN_ANGLE = 0.0037819f;

/** Hard-recolor pixel count between two captured frames. */
inline int hard_recolor_count(const std::vector<Pixel> &a,
                              const std::vector<Pixel> &b) {
  int n = 0;
  for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
    const int dr = std::abs((int)a[i].r - b[i].r);
    const int dg = std::abs((int)a[i].g - b[i].g);
    const int db = std::abs((int)a[i].b - b[i].b);
    if (std::max({dr, dg, db}) > STRAP_OPEN_HARD)
      ++n;
  }
  return n;
}

/** Renders the effect's current-node hankin mesh at a fixed camera into a full
 * frame, at `angle` with strap opening-fade `fade`. */
template <int W, int H>
inline void capture_opening(HankinSolids<W, H> &fx, float angle, float fade,
                            std::vector<Pixel> &out, float close_blend = 1.0f,
                            float terminal_fade = 1.0f, int cycle_frame = 1,
                            float star_close = 1.0f) {
  using Probe = conway_soak_tests::HankinWalkProbe;
  Arena scratch(cc_temp_buf, sizeof(cc_temp_buf));
  {
    Canvas c(fx);
    Probe::render_at_angle(fx, c, angle, cycle_frame, fade, close_blend,
                           terminal_fade, star_close, scratch);
  }
  fx.advance_display();
  out.resize(static_cast<size_t>(W) * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      out[static_cast<size_t>(y) * W + x] = fx.get_pixel(x, y);
}

/**
 * @brief Pins the strap opening-fade at a cycle start: the newborn straps,
 *        faded, leave the angle-0 bookend nearly unchanged, whereas drawn at
 *        full coverage they hard-recolor the cut star interiors.
 * @details Drives to an arrival that cuts enough star interior to pin, then
 *          renders three frames at the same (unadvanced) camera: the angle-0
 *          bookend, the first strap frame at full coverage (the pre-fix
 *          artifact), and the same frame faded (the fix). The full-coverage pop
 *          must be large and the faded pop small, so the assertion is red
 *          without the fade and green with it.
 */
inline void test_strap_open_fade() {
  reset_globals();
  using Probe = conway_soak_tests::HankinWalkProbe;
  constexpr int W = 288, H = 144;
  HankinSolids<W, H> fx;
  fx.init();

  // How much star interior the newborn straps cut depends on which node the
  // walk arrives at, and some nodes barely cut any, so search arrivals for one
  // whose full-coverage artifact is big enough to pin.
  constexpr int MAX_ARRIVALS = 12;
  constexpr int MIN_POP_FULL = 400;
  float fade = 0.0f;
  int pop_full = 0, pop_faded = 0;
  std::vector<Pixel> bookend, full, faded;
  int arrival = 0;
  for (; arrival < MAX_ARRIVALS; ++arrival) {
    const int prev_node = Probe::node(fx);
    int guard = 0;
    while (Probe::node(fx) == prev_node && guard++ < 400) {
      fx.draw_frame();
      fx.advance_display();
    }
    HS_EXPECT_LT(guard, 400);

    fade = Probe::strap_open_fade(fx, 1);
    capture_opening(fx, 0.0f, 1.0f, bookend);           // straps zero-area
    capture_opening(fx, STRAP_OPEN_ANGLE, 1.0f, full);  // pre-fix: full coverage
    capture_opening(fx, STRAP_OPEN_ANGLE, fade, faded); // fix: faded coverage
    pop_full = hard_recolor_count(bookend, full);
    pop_faded = hard_recolor_count(bookend, faded);
    if (pop_full > MIN_POP_FULL)
      break;
  }
  HS_EXPECT_LT(arrival, MAX_ARRIVALS);

  HS_EXPECT_GT(fade, 0.0f);
  HS_EXPECT_LT(fade, 0.2f); // the first frame is early in the window
  std::printf("  [strap-open-fade] full-coverage pop=%d, faded pop=%d px\n",
              pop_full, pop_faded);

  // The artifact must be real (else the pin is vacuous), and the fade must cut
  // it to a small residue (the inherent star-face deformation at the step).
  HS_EXPECT_GT(pop_full, MIN_POP_FULL);
  HS_EXPECT_LT(pop_faded, pop_full / 4);
  HS_EXPECT_LT(pop_faded, 200);
}

/** Summed absolute channel change between two captures — the visible weight of
 * a transition, which a bare changed-pixel count misses when the same pixels
 * change by less. */
inline long long frame_energy(const std::vector<Pixel> &a,
                              const std::vector<Pixel> &b) {
  long long e = 0;
  for (size_t i = 0; i < a.size() && i < b.size(); ++i)
    e += std::abs((int)a[i].r - b[i].r) + std::abs((int)a[i].g - b[i].g) +
         std::abs((int)a[i].b - b[i].b);
  return e;
}

/**
 * @brief Pins the cycle-close shaping: the last strap frame dissolves into the
 *        host face's rim instead of winking out as a bright line.
 * @details Renders the last strap frame and the closing bookend at a fixed
 *          camera, unshaped and shaped, and compares the transition's energy.
 *          Unshaped, the strap holds its own ramp interior — a colored line
 *          against the rim it sits on — and the whole line disappears in one
 *          frame. Shaped, it is already the rim color and its terminal sliver
 *          is faded, so the step carries far less.
 */
inline void test_strap_close_dissolve() {
  reset_globals();
  using Probe = conway_soak_tests::HankinWalkProbe;
  constexpr int W = 288, H = 144;
  HankinSolids<W, H> fx;
  fx.init();

  int prev_node = Probe::node(fx);
  int guard = 0;
  while (Probe::node(fx) == prev_node && guard++ < 400) {
    fx.draw_frame();
    fx.advance_display();
  }
  HS_EXPECT_LT(guard, 400);

  // Frame adjacent to the closing bookend: the sweep is a full sine period, so
  // it samples the same angle as the opening's first strap frame.
  const int duration = 64, cf = 63;
  const float close_blend = Probe::strap_open_fade(fx, duration - cf);
  const float term = hs::clamp(static_cast<float>(duration - cf) /
                                   Probe::strap_terminal_frames(fx),
                               0.0f, 1.0f);
  HS_EXPECT_LT(close_blend, 0.2f);
  HS_EXPECT_LT(term, 1.0f);

  std::vector<Pixel> bookend, last_plain, last_shaped;
  capture_opening(fx, 0.0f, 1.0f, bookend, 1.0f, 1.0f, cf);
  capture_opening(fx, STRAP_OPEN_ANGLE, 1.0f, last_plain, 1.0f, 1.0f, cf);
  capture_opening(fx, STRAP_OPEN_ANGLE, 1.0f, last_shaped, close_blend, term,
                  cf);

  const long long e_plain = frame_energy(last_plain, bookend);
  const long long e_shaped = frame_energy(last_shaped, bookend);
  std::printf("  [strap-close] wink energy unshaped=%lld shaped=%lld\n",
              e_plain, e_shaped);

  // The wink must be real, and the shaping must materially cut it.
  HS_EXPECT_GT(e_plain, 1000000);
  HS_EXPECT_LT(e_shaped, e_plain * 3 / 4);
}

/** Interlace angle one frame either side of the sweep's midpoint, where the
 * star face has closed to a last sliver (sin_wave at progress 31/64). */
constexpr float STAR_CLOSE_ANGLE = 1.56701f;

/**
 * @brief Pins the mid-sweep star dissolve: the star's last sliver carries the
 *        rosette rim color it is closing into, not its own interior color.
 * @details The mirror of the strap close. At mid-sweep the star face shuts to
 *          zero area and the rosettes hosted inside it fill its place, so the
 *          transition into the fully-closed midpoint must cost less when the
 *          star has already taken their rim color. Rendered at a fixed camera,
 *          shaped against unshaped.
 */
inline void test_star_midpoint_dissolve() {
  reset_globals();
  using Probe = conway_soak_tests::HankinWalkProbe;
  constexpr int W = 288, H = 144;
  HankinSolids<W, H> fx;
  fx.init();

  int prev_node = Probe::node(fx);
  int guard = 0;
  while (Probe::node(fx) == prev_node && guard++ < 400) {
    fx.draw_frame();
    fx.advance_display();
  }
  HS_EXPECT_LT(guard, 400);

  // Drive on until a node whose rosettes carry a different palette from the
  // star hosting them. On some solids the mod-NUM_PALETTES wrap aliases the two
  // onto one slot, and the star then has nothing to cross-fade onto — a real
  // no-op, but one that would make the pin below vacuous.
  const auto has_distinct_rim = [&]() {
    const uint8_t *own = Probe::node_face_palette(fx);
    const uint8_t *rim = Probe::star_rim_palette(fx);
    for (size_t j = 0; j < Probe::node_faces(fx); ++j)
      if (rim[j] != own[j])
        return true;
    return false;
  };
  int hops = 0;
  while (!has_distinct_rim() && hops++ < 40) {
    const int at = Probe::node(fx);
    int spin = 0;
    while (Probe::node(fx) == at && spin++ < 400) {
      fx.draw_frame();
      fx.advance_display();
    }
  }
  HS_EXPECT_TRUE(has_distinct_rim());

  const int duration = 64, mid = duration / 2, cf = mid - 1;
  const float star_blend = Probe::strap_open_fade(fx, mid - cf);
  HS_EXPECT_GT(star_blend, 0.0f);
  HS_EXPECT_LT(star_blend, 1.0f);

  std::vector<Pixel> closed, sliver_plain, sliver_shaped;
  // The midpoint itself: the star is gone, so its shaping cannot apply.
  capture_opening(fx, PI_F / 2.0f, 1.0f, closed, 1.0f, 1.0f, mid, 1.0f);
  capture_opening(fx, STAR_CLOSE_ANGLE, 1.0f, sliver_plain, 1.0f, 1.0f, cf,
                  1.0f);
  capture_opening(fx, STAR_CLOSE_ANGLE, 1.0f, sliver_shaped, 1.0f, 1.0f, cf,
                  star_blend);

  // Only star fragments are shaped, so the pixels where the two renders differ
  // are exactly the star's. Scoring the whole frame drowns them in the
  // rosettes, which fill nearly all of it by the midpoint.
  long long e_plain = 0, e_shaped = 0;
  size_t star_px = 0;
  for (size_t i = 0; i < sliver_plain.size(); ++i) {
    if (sliver_plain[i].r == sliver_shaped[i].r &&
        sliver_plain[i].g == sliver_shaped[i].g &&
        sliver_plain[i].b == sliver_shaped[i].b)
      continue;
    ++star_px;
    e_plain += std::abs((int)sliver_plain[i].r - closed[i].r) +
               std::abs((int)sliver_plain[i].g - closed[i].g) +
               std::abs((int)sliver_plain[i].b - closed[i].b);
    e_shaped += std::abs((int)sliver_shaped[i].r - closed[i].r) +
                std::abs((int)sliver_shaped[i].g - closed[i].g) +
                std::abs((int)sliver_shaped[i].b - closed[i].b);
  }
  std::printf("  [star-mid] %zu shaped px, close energy unshaped=%lld "
              "shaped=%lld\n",
              star_px, e_plain, e_shaped);

  // The shaping must reach the closing star, and must carry it toward the ramp
  // the midpoint leaves in its place.
  HS_EXPECT_GT(star_px, 100u);
  HS_EXPECT_GT(e_plain, 100000);
  HS_EXPECT_LT(e_shaped, e_plain);

  // Symmetric on the far side: the weight keys off distance to the midpoint, so
  // the star comes back out of the rosette ramp as it reopens rather than
  // snapping to its own color.
  const int reopen_cf = mid + 1;
  HS_EXPECT_NEAR(Probe::strap_open_fade(fx, reopen_cf - mid), star_blend,
                 1e-6f);
  std::vector<Pixel> reopen_shaped;
  capture_opening(fx, STAR_CLOSE_ANGLE, 1.0f, reopen_shaped, 1.0f, 1.0f,
                  reopen_cf, star_blend);
  size_t reopen_px = 0;
  for (size_t i = 0; i < reopen_shaped.size(); ++i)
    if (reopen_shaped[i].r != sliver_plain[i].r ||
        reopen_shaped[i].g != sliver_plain[i].g ||
        reopen_shaped[i].b != sliver_plain[i].b)
      ++reopen_px;
  std::printf("  [star-mid] reopening side shaped px=%zu\n", reopen_px);
  HS_EXPECT_GT(reopen_px, 100u);
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

  test_leg_start_seed_frame_continuity();
  test_palette_slots_stable_within_cycle();
  test_strap_crossfade_across_cycle_start();
  test_strap_crossfade_seed_swept();
  test_strap_open_fade();
  test_strap_close_dissolve();
  test_star_midpoint_dissolve();

  return fixture.result();
}

} // namespace conway_continuity_tests
} // namespace hs_test
