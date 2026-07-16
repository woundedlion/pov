/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"
#include "core/mesh/spatial.h"

#include <cmath>
#include <span>

/**
 * @brief Spherical Voronoi effect: scatters sites on the unit sphere, animates
 *        them, and shades each pixel by its nearest site with edge-sharpening
 *        and optional cell borders.
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 */
template <int W, int H> class Voronoi : public Effect {
public:
  /**
   * @brief One Voronoi seed: position on the unit sphere, the rotation axis it
   *        spins about, and the cell fill color.
   */
  struct Site {
    Vector pos;   /**< Position on the unit sphere. */
    Vector axis;  /**< Unit rotation axis the site spins about. */
    Color4 color; /**< Cell fill color (16-bit linear channels). */
  };

  /**
   * @brief Constructs the effect with the templated framebuffer dimensions.
   */
  HS_COLD_MEMBER Voronoi() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Configures arenas, registers GUI params, allocates the sites buffer,
   *        and seeds the initial sites.
   */
  void init() override {
    // Persistent holds the sites buffer; scratch_arena_a holds the per-frame
    // KD-tree (positions + nodes + build indices).
    configure_arenas(GLOBAL_ARENA_SIZE - 64 * 1024, 64 * 1024, 0);

    register_param("Num Sites", &params.num_sites, 1.0f,
                  static_cast<float>(MAX_SITES));
    register_param("Speed", &params.speed, 0.0f, 100.0f);
    register_param("Sharpness", &params.sharpness, 0.0f, 500.0f);
    register_param("Border Thick", &params.border_thickness, 0.0f, 0.1f);

    sites_buffer.bind(persistent_arena, MAX_SITES);
    seed_sites();
  }

  /**
   * @brief Animates the sites, builds a per-frame KD-tree, and shades each
   *        pixel by its nearest site (with edge sharpening and optional
   *        borders).
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(vo_buffer_wait);
      return Canvas(*this);
    }();

    // Re-seed when the GUI changes the site count (integer change only, so
    // dragging within a bucket doesn't thrash).
    {
      HS_PROFILE(vo_animate);
      if (active_site_count() != current_num_sites)
        seed_sites();

      float s = logf(params.speed + 1.0f) * 0.005f;

      for (size_t i = 0; i < sites_buffer.size(); ++i) {
        auto &site = sites_buffer[i];
        // Renormalize: rotate() drifts |pos| off the unit sphere over a long run,
        // and the unit-site invariants (nearest-by-Euclidean == nearest-by-max-dot,
        // and the border acosf(dot) staying in range) require unit vectors.
        Quaternion q = make_rotation(site.axis, s);
        site.pos = rotate(site.pos, q).normalized();
      }
    }

    // Build a KD-tree over the moving site positions once per frame. On the unit
    // sphere nearest-by-Euclidean == nearest-by-max-dot (|p-s|^2 = 2 - 2*p*s),
    // so the k=2 query is exact.
    ScratchScope scope_guard(scratch_arena_a);
    Vector *positions = static_cast<Vector *>(scratch_arena_a.allocate(
        sites_buffer.size() * sizeof(Vector), alignof(Vector)));
    // IIFE times the per-frame KD build while keeping `tree` at frame scope
    // (guaranteed copy elision on the prvalue return).
    KDTree tree = [&]() -> KDTree {
      HS_PROFILE(vo_kdtree);
      for (size_t i = 0; i < sites_buffer.size(); ++i)
        positions[i] = sites_buffer[i].pos;
      return KDTree(scratch_arena_a,
                    std::span<const Vector>(positions, sites_buffer.size()));
    }();

    // One node for all per-pixel work (corner pre-pass + shading loop); never
    // scope inside the per-pixel loop, filter.h counts blended pixels.
    HS_PROFILE(vo_shade);

    // Resolves the final color from the already-identified nearest pair: the
    // nearest site index i0 and its dot d0 (the larger), and — when a second
    // neighbor exists — the second site index i1 and its dot d1.
    auto shade = [&](uint16_t i0, float d0, bool has_second, uint16_t i1,
                     float d1) -> Color4 {
      const Site &best_site = sites_buffer[i0];
      float max_dot1 = d0;
      float max_dot2 = has_second ? d1 : -2.0f;

      Color4 c = best_site.color;

      // Border sharpening: a larger sharpness saturates `factor` for smaller
      // nearest/second-nearest gaps, shrinking the cross-cell blend band.
      if (has_second && params.sharpness > 0.0f) {
        const Site &sec_site = sites_buffer[i1];
        float diff = max_dot1 - max_dot2;
        float factor = std::min(1.0f, diff * params.sharpness);
        factor = quintic_kernel(factor);
        float t = 0.5f + 0.5f * factor;

        uint16_t frac = static_cast<uint16_t>(t * 65535.0f + 0.5f);
        c.color = sec_site.color.color.lerp16(best_site.color.color, frac);
      }

      // Borders — driven entirely by the "Border Thick" slider: a thickness of
      // 0 disables them (and skips the two acosf calls below), any positive
      // value paints the seam between the nearest two sites. d0 is the nearest,
      // so max_dot1 >= max_dot2 → dist1 <= dist2 and the cell gap is non-negative.
      if (params.border_thickness > 0.0f && has_second) {
        float dist1 = acosf(hs::clamp(max_dot1, -1.0f, 1.0f));
        float dist2 = acosf(hs::clamp(max_dot2, -1.0f, 1.0f));
        if (dist2 - dist1 < params.border_thickness) {
          // Paint the seam black. The Scan sink writes color*alpha, so an
          // alpha-0 fragment collapses to (0,0,0) regardless of its RGB.
          c = Color4(0, 0, 0, 0);
        }
      }

      return c;
    };

    // Canonical (order-independent) nearest-pair identity at a sample point:
    // the two query orders along a cell seam (best/second swap) map to the
    // same {lo, hi} set.
    auto classify = [&](const Vector &p) -> CellId {
      auto knn = tree.nearest(p, 2);
      uint16_t a = knn[0].original_index;
      bool has_second = knn.size() > 1;
      uint16_t b = has_second ? knn[1].original_index : a;
      return {std::min(a, b), std::max(a, b), has_second};
    };

    // Coarse-grid coherence: classify the nearest pair once per coarse-grid
    // corner, then shade every pixel of a block from the deduped union of its
    // four corners' pairs (<= 8 candidate sites) by an exact top-2 dot scan —
    // no per-pixel KD query. A site owning any pixel of a block reaches at
    // least one corner's pair in all but degenerate layouts; a cell missed by
    // all four corners is dropped.
    auto &cr = canvas.clip();
    Scan::Shader::check_lut_domain<W, H>(cr);
    const int x0 = cr.x_start;
    const int x1 = cr.x_end;
    const int y0 = cr.render_y_start();
    const int y1 = cr.render_y_end();
    if (x1 <= x0 || y1 <= y0)
      return;

    // Voronoi cell pixel size falls as ~1/sqrt(num_sites) (smallest near the
    // poles, where rows crowd), so a fixed block eventually straddles more than
    // one cell and misses small cells. Shrink the block toward the cell's
    // vertical pixel extent (rows map uniformly over [0,π]); the floor keeps the
    // corner grid inside the scratch budget the static_assert below pins.
    const float cell_px =
        (2.0f * H / PI_F) / sqrtf(static_cast<float>(sites_buffer.size()));
    const int B = hs::clamp(static_cast<int>(cell_px), COHERENCE_BLOCK_MIN,
                            COHERENCE_BLOCK);
    const int nbx = (x1 - x0 - 1) / B + 2; // corner columns spanning [x0, x1)
    const int nby = (y1 - y0 - 1) / B + 2; // corner rows spanning    [y0, y1)
    CellId *cells = static_cast<CellId *>(
        scratch_arena_a.allocate(nbx * nby * sizeof(CellId), alignof(CellId)));
    // Clamp to the last valid pixel (x1/y1 are exclusive) so every classified
    // point indexes the trig LUT.
    auto corner_x = [&](int j) { return std::min(x0 + j * B, x1 - 1); };
    auto corner_y = [&](int k) { return std::min(y0 + k * B, y1 - 1); };

    // Timer spans the corner-classification pre-pass too: its ~W*H/B^2 KD
    // queries are the frame's only KD queries, so getRenderUs must include
    // them.
    Scan::ScopedRenderTimer timer_guard(canvas);
    for (int k = 0; k < nby; ++k)
      for (int j = 0; j < nbx; ++j)
        cells[k * nbx + j] =
            classify(pixel_to_vector<W, H>(corner_x(j), corner_y(k)));

    // One candidate set per block column, rebuilt on each block-row change.
    // Positions are copied in so the per-pixel scan runs over contiguous data.
    const int nblk = nbx - 1;
    CandSet *cands = static_cast<CandSet *>(scratch_arena_a.allocate(
        nblk * sizeof(CandSet), alignof(CandSet)));
    auto build_candidate_row = [&](int ky) {
      for (int jx = 0; jx < nblk; ++jx) {
        CandSet &cs = cands[jx];
        cs.n = 0;
        auto add = [&](uint16_t s) {
          for (uint8_t i = 0; i < cs.n; ++i)
            if (cs.idx[i] == s)
              return;
          cs.idx[cs.n] = s;
          cs.pos[cs.n] = positions[s];
          ++cs.n;
        };
        for (const CellId *c :
             {&cells[ky * nbx + jx], &cells[ky * nbx + jx + 1],
              &cells[(ky + 1) * nbx + jx], &cells[(ky + 1) * nbx + jx + 1]}) {
          add(c->lo);
          add(c->hi);
        }
      }
    };

    // SAMPLES=1: one sample at pixel center. Mirrors Scan::Shader::draw<W, H,
    // 1>'s clip iteration and telemetry (open-coded only because the coarse
    // grid needs the integer pixel coordinates the generic shader callback
    // does not expose).
    int last_ky = -1;
    for (int y = y0; y < y1; ++y) {
      const int ky = (y - y0) / B;
      if (ky != last_ky) {
        build_candidate_row(ky);
        last_ky = ky;
      }
      for (int x = x0; x < x1; ++x) {
        const CandSet &cs = cands[(x - x0) / B];

        Vector p = pixel_to_vector<W, H>(x, y);
        float d0 = -2.0f, d1 = -2.0f;
        uint8_t b0 = 0, b1 = 0;
        for (uint8_t i = 0; i < cs.n; ++i) {
          float d = dot(p, cs.pos[i]);
          if (d > d0) {
            d1 = d0;
            b1 = b0;
            d0 = d;
            b0 = i;
          } else if (d > d1) {
            d1 = d;
            b1 = i;
          }
        }
        Color4 sample = (cs.n >= 2)
                            ? shade(cs.idx[b0], d0, true, cs.idx[b1], d1)
                            : shade(cs.idx[b0], d0, false, cs.idx[b0], d0);
        canvas(x, y) = sample.color * sample.alpha;
      }
    }
  }

  static constexpr int MAX_SITES = 400; /**< Buffer capacity; the sites buffer
                                             is allocated once at this size. */
  static constexpr int COHERENCE_BLOCK = 8; /**< Coarse-coherence block edge in
      pixels at low site counts: each pixel shades from the union of its block
      corners' nearest pairs. Smaller is safer (fewer missed sub-block cells)
      but classifies more corners; the render path shrinks the block toward
      COHERENCE_BLOCK_MIN as the site count rises. */
  static constexpr int COHERENCE_BLOCK_MIN = 4; /**< Smallest adaptive block edge.
      Floors the per-frame block so the corner grid stays within the scratch
      budget (pinned by the static_assert below); ~matches the cell pixel size at
      MAX_SITES. */

  /** @brief Canonical (order-independent) nearest-pair identity at a sample
   *  point; the coarse-grid corner classifier stores one per corner (see the
   *  render path). At class scope so the scratch-budget static_assert below can
   *  size the corner grid against sizeof(CellId). */
  struct CellId {
    uint16_t lo;    /**< min(nearest, second) site index. */
    uint16_t hi;    /**< max(nearest, second) site index. */
    bool has_second; /**< Whether a second neighbor exists (>= 2 sites). */
  };

  /** @brief Shading candidates for one block: the deduped union of its four
   *  corners' CellId pairs, positions copied in for the per-pixel dot scan. */
  struct CandSet {
    Vector pos[8];   /**< Candidate site positions (parallel to idx). */
    uint16_t idx[8]; /**< Candidate site indices into sites_buffer. */
    uint8_t n;       /**< Number of distinct candidates (1..8). */
  };

  // Compile-time high-water check for the 64 KB scratch_arena_a reserve. Two
  // transient peaks share that arena, with positions + KD nodes live across both:
  //   build:   positions + KD nodes + KD build-index scratch
  //   shading: positions + KD nodes + corner cells + candidate-row sets
  static constexpr size_t SCRATCH_A_BYTES = 64 * 1024;
  static constexpr size_t POSITIONS_BYTES = size_t(MAX_SITES) * sizeof(Vector);
  static constexpr size_t KD_NODES_BYTES = size_t(MAX_SITES) * sizeof(KDNode);
  static constexpr size_t KD_BUILD_SCRATCH_BYTES = size_t(MAX_SITES) * sizeof(int);
  // Corner grid spans the full canvas in block-px steps, +2 for the inclusive
  // [0,W)/[0,H) end corners (mirrors render()'s nbx/nby at full clip). Sized at
  // COHERENCE_BLOCK_MIN — the worst case (most corners) the adaptive block hits.
  static constexpr size_t CORNER_COLS = size_t((W - 1) / COHERENCE_BLOCK_MIN + 2);
  static constexpr size_t CORNER_ROWS = size_t((H - 1) / COHERENCE_BLOCK_MIN + 2);
  static constexpr size_t CELLS_BYTES = CORNER_COLS * CORNER_ROWS * sizeof(CellId);
  static constexpr size_t CAND_ROW_BYTES = (CORNER_COLS - 1) * sizeof(CandSet);
  static constexpr size_t SCRATCH_HIGH_WATER =
      POSITIONS_BYTES + KD_NODES_BYTES +
      (KD_BUILD_SCRATCH_BYTES > CELLS_BYTES + CAND_ROW_BYTES
           ? KD_BUILD_SCRATCH_BYTES
           : CELLS_BYTES + CAND_ROW_BYTES);
  static_assert(SCRATCH_HIGH_WATER <= SCRATCH_A_BYTES,
                "Voronoi scratch_arena_a budget (64 KB) too small for MAX_SITES "
                "positions + KD-tree + coarse-grid cells; raise the reserve in "
                "init() or lower MAX_SITES / coarsen COHERENCE_BLOCK");

  int current_num_sites = 0; /**< Count currently seeded; re-seeds (clear +
                                  refill, no realloc) when the slider changes. */
  ArenaVector<Site> sites_buffer; /**< Active Voronoi sites for the frame. */

  /**
   * @brief Returns the active site count from the "Num Sites" slider.
   * @return Slider value rounded to an integer and clamped to [1, MAX_SITES].
   */
  int active_site_count() const {
    return hs::clamp(static_cast<int>(params.num_sites), 1, MAX_SITES);
  }

  /**
   * @brief (Re)seeds the active sites for the current "Num Sites" slider value.
   * @details Clears and refills up to the active count (no re-allocation). Sites
   *          are placed via the Fibonacci-sphere distribution for an even
   *          spread, each with a random spin axis and a palette color.
   */
  HS_COLD_MEMBER void seed_sites() {
    const int n = active_site_count();
    sites_buffer.clear();

    const float golden_angle = PI_F * (3.0f - sqrtf(5.0f));

    for (int i = 0; i < n; i++) {
      // Guard n == 1: a 0 denominator would give NaN y. A single site sits at
      // the pole (y = 1).
      int span = n > 1 ? n - 1 : 1;
      float y = 1.0f - (i / (float)span) * 2.0f;
      float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
      float theta = golden_angle * i;

      float x = cosf(theta) * radius;
      float z = sinf(theta) * radius;

      Vector pos = Vector(x, y, z);

      Vector axis = random_vector();

      float t = i / (float)(n > 1 ? n - 1 : 1);
      Color4 color = Palettes::RICH_SUNSET.get(t);

      sites_buffer.push_back({pos, axis, color});
    }
    current_num_sites = n;
  }

private:
  /**
   * @brief Live-tunable GUI parameters for the Voronoi effect.
   */
  struct Params {
    float num_sites = 200.0f;     /**< Live-tunable site count (GUI slider). */
    float speed = 20.0f;          /**< Site spin rate (GUI slider). */
    float sharpness = 100.0f;     /**< Edge sharpening; larger narrows the
                                       border blend band. */
    float border_thickness = 0.0f; /**< Cell-seam border width; 0 disables. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Voronoi)
